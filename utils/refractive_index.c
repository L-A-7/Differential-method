/*
 *	refractive_index.c
 *
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <complex.h>

#define ROUND(x) ((int)(x<0 ? x-0.5 : x+0.5))
#define STR_SIZE 5000
#define SIZE_STR_BUFFER 200
#define SIZE_LINE_BUFFER 50000
#define CHAR_COMMENT '#'

void err_message(){
	fprintf(stderr,	"usage : refractive_index filename lambda [Cauchy/Lookup]\n where 'filename' points to a file");
}

char *label_search(char *str,const char *label);
void skip_comment(char *str_in_out);
int lire_ligne(FILE *fp, char *line);
int lire_tab(const char *nom_fichier, const char *label, double *tab, int N);
int count_tab(const char *nom_fichier, const char *label);
complex refractive_index(char* name,double lambda, char* method);

int main(int argc, char *argv[]){
	
	int i;
	complex ind;
	double lambda;
	char filename[STR_SIZE], method[STR_SIZE], *endptr;
	
	/* Copie des arguments de la ligne de commande */
	char **argvcp; 
	argvcp = (char **) malloc(sizeof(char*)*argc);
	argvcp[0] = (char*) malloc(sizeof(char)*STR_SIZE*argc);
	for(i=1;i<=argc-1;i++){
		argvcp[i] = argvcp[i-1] + STR_SIZE;
		strncpy(argvcp[i], argv[i],STR_SIZE);
	}
	
	/* Vérification de la présence du nombre minimal d'options */
	if (argc < 2){
		err_message();
		return 1;
	}
	
	/* Reading the arguments */
	strncpy(filename, argvcp[1],STR_SIZE);
	lambda = strtod(argvcp[2], &endptr);
	if (argvcp[2] == endptr){
		fprintf(stderr,"ERROR, %s can't read argument lambda. Exiting.",__FILE__);
		exit(EXIT_FAILURE);
	}
	strncpy(method, "Lookup",STR_SIZE);

/*	strncpy(method, argvcp[3],STR_SIZE);*/ /*Ajouter verifications, etc..*/
	
	/* index determination */
	ind = refractive_index(filename,lambda, method);
	
	fprintf(stdout,"%f +i%f",creal(ind),cimag(ind));

	
	return 0;
}

complex refractive_index(char* filename,double lambda, char* method)
{
	double lambda_angstrom = 10*lambda;
	double cauchy_n[5],cauchy_k[5],n_power[6],k_power[6],index_n,index_k;
	complex index;
	int i,Nb_cauchy = 5;
	
	/* Cauchy method */
	if (!strcmp(method,"Cauchy")){
		lire_tab(filename, "CAUCHY_N", cauchy_n, Nb_cauchy);
		lire_tab(filename, "CAUCHY_K", cauchy_k, Nb_cauchy);
		lire_tab(filename, "N_POWERS", n_power, Nb_cauchy+1);
		lire_tab(filename, "K_POWERS", k_power, Nb_cauchy+1);
/*SaveDbleTab2file (cauchy_n,  Nb_cauchy,"stdout", " ");printf("\n");
SaveDbleTab2file (n_power, Nb_cauchy,"stdout", " ");printf("\n");
SaveDbleTab2file (cauchy_k,  Nb_cauchy,"stdout", " ");printf("\n");
SaveDbleTab2file (k_power, Nb_cauchy,"stdout", " ");printf("\n");*/

		index_n = 0;
		index_k = 0;
		for (i=0;i<=Nb_cauchy-1;i++){
			index_n += cauchy_n[i]*pow(lambda_angstrom,n_power[i]);
			index_k += cauchy_k[i]*pow(lambda_angstrom,k_power[i]);
		}
/*printf("Milieu: %s, lambda = %f, n = %f, k = %f\n",filename,lambda,index_n,index_k);
*/				
		index = index_n + I*index_k; 
		return index;
		
	}
	/* Lookup-Table method */
	if (!strcmp(method,"Lookup")){
		double *tab_n, *tab_k, *tab_lambda, *table_tmp;
		int npoints;
		npoints = ROUND((double) count_tab(filename, "")/3);
		table_tmp= (double *) malloc(sizeof(double)*3*((int)npoints));
		tab_lambda= (double *) malloc(sizeof(double)*((int)npoints));
		tab_n= (double *) malloc(sizeof(double)*((int)npoints));
		tab_k= (double *) malloc(sizeof(double)*((int)npoints));
		lire_tab(filename, "", table_tmp, 3*(int)npoints);
		/* Réarrangement en plusieurs tableaux */
		for (i=0; i<=(int)npoints -1;i++){
			tab_lambda[i] = table_tmp[3*i];
			tab_n[i] = table_tmp[3*i+1];
			tab_k[i] = table_tmp[3*i+2];
		}
		/* Recherche de l'indice de tableau pour lambda */
		int num_min = 0;
		int num_max = (int) npoints-1;
		int numero = num_max>>1;
		while(!(tab_lambda[numero] <= lambda && lambda < tab_lambda[numero+1])){
			numero=num_min+((num_max-num_min)>>1);
			if (lambda < tab_lambda[numero]){
				num_max=numero;
			}
			if (tab_lambda[numero+1] <= lambda){
				num_min=numero;
			}
			if (num_min == num_max){
				fprintf(stderr,"refractive_index, look-up table incompatible with lambda = %f for %s",lambda,filename);
				exit(EXIT_FAILURE);
			}
		}
		/* Interpolation linéaire */
		double l1 = tab_lambda[numero];
		double l2 = tab_lambda[numero+1];
		double n1 = tab_n[numero];
		double n2 = tab_n[numero+1];
		double k1 = tab_k[numero];
		double k2 = tab_k[numero+1];
		index_n= (lambda-l1)*(n2-n1)/(l2-l1) + n1;
		index_k= (lambda-l1)*(k2-k1)/(l2-l1) + k1;

		free(tab_lambda);
		free(tab_n);
		free(tab_k);
		free(table_tmp);
		
/*printf("lambda= %f, lambda1= %f, lambda2= %f\n",lambda,tab_lambda[numero],tab_lambda[numero+1]);
printf("lambda= %f, n = %f, k= %f\n",lambda,index_n,index_k);
*/		index = index_n + I*index_k; 

		return index;
	}
	return -1;/* ERREUR, ne devrait pas arriver là ... */
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int lire_tab(char *nom_fichier, const char *label, double *tab, int N)
 *
 *	\brief	Lit N valeurs de format double dans un fichier et les stocke dans un tableau   \n 
 *			Les valeurs doivent être séparées par un ou plusieurs espaces, tabulations     \n
 *			ou sauts de lignes et précédées d'un label éventuellement suivi d'un signe '='.\n
 *          ex. : (...) tab1 = 3.4  4.5e-3  +46  -7.6e+2 ...                               \n
 *          Remarque : Pour lire un tableau sans label, donner "" comme label.             \n
 *
 *	\return	0 si succès, 1 si le nombre d'éléments lus diffère de N ou si le label n'a pas été trouvé.
 *
 *	\todo	RENDRE PLUS ROBUSTE : PAS DE BUFFER OVERFLOW AU CAS OU IL Y A PLUS DE N LIGNES
 *
 */
/*---------------------------------------------------------------------------------------------*/
int lire_tab(const char *nom_fichier, const char *label, double *tab, int N)
{
	FILE *fp;
	char *pos, *pos2, *endptr, line[SIZE_LINE_BUFFER];
	int cpt=0, line_cpt=0;
	double tmp;
	
	if (!(fp = fopen(nom_fichier,"r"))){
		fprintf(stderr, "%s line %d: ERROR, can't open %s\n",__FILE__, __LINE__,nom_fichier);
		return 1;
	}

	/* Recherche du label */
	while(!feof(fp)){ 
		if (lire_ligne(fp,line) != 0) goto LECTURE_FINIE;
		line_cpt++;
		skip_comment(line);
		if ((pos = label_search(line,label)) != NULL){ /* on cherche le label */	
			pos += strlen(label); 
			if((pos2 = strchr(pos,'=')) != NULL) *pos2 = ' '; /* on remplace '=' par ' ' */
			goto LABEL_TROUVE;
		}
	}

	fprintf(stderr,"%s ligne %d : ERREUR, le label '%s' n'a pas été trouvé dans %s \n", __FILE__, __LINE__, label, nom_fichier);
	fclose(fp);
	return 1;

	/* Lecture des valeurs */
	while(!feof(fp)){ /* Tant qu'on est pas à la fin du fichier */
		if (lire_ligne(fp,line) != 0) {goto LECTURE_FINIE;}
		line_cpt++;
		skip_comment(line); 
		pos = line;
		
	LABEL_TROUVE :
		while(isspace(*pos)) pos++; /* on élimine les espaces */
		while(pos < line+strlen(line)) { /* Tant qu'on est pas à la fin de la ligne */
			tmp = strtod(pos, &endptr); /* on lit le 'double' */
			if (pos == endptr) {goto LECTURE_FINIE;} /* conversion ratée */
			if (cpt <= N+1) tab[cpt] = tmp;
			cpt ++;
			pos = endptr;
			while(isspace(*pos)) pos++; /* on élimine les espaces */
		}
	}

	LECTURE_FINIE:
	if (cpt != N) {
		fprintf(stderr, "%s ligne %d : ERREUR, %s contient %d valeurs au lieu de %d dans %s (ligne %d)\n", __FILE__, __LINE__, label, cpt, N, nom_fichier, line_cpt);
		fclose(fp);
		return cpt;
	}
	fclose(fp);
	return 0;
}
/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int count_tab(char *nom_fichier, const char *label)
 *
 *	\brief	count the number of elements in a tab, 
 *
 */
/*---------------------------------------------------------------------------------------------*/
int count_tab(const char *nom_fichier, const char *label)
{
	FILE *fp;
	char *pos, *pos2, *endptr, line[SIZE_LINE_BUFFER];
	int cpt=0, line_cpt=0;
	double tmp;
	
	if (!(fp = fopen(nom_fichier,"r"))){
		fprintf(stderr, "%s line %d: ERROR, can't open %s\n",__FILE__, __LINE__,nom_fichier);
		return 1;
	}

	/* Recherche du label */
	while(!feof(fp)){ 
		if (lire_ligne(fp,line) != 0) goto LECTURE_FINIE;
		line_cpt++;
		skip_comment(line);
		if ((pos = label_search(line,label)) != NULL){ /* on cherche le label */	
			pos += strlen(label); 
			if((pos2 = strchr(pos,'=')) != NULL) *pos2 = ' '; /* on remplace '=' par ' ' */
			goto LABEL_TROUVE;
		}
	}

	fprintf(stderr,"%s ligne %d : ERREUR, le label '%s' n'a pas été trouvé dans %s \n", __FILE__, __LINE__, label, nom_fichier);
	fclose(fp);
	return 1;

	/* Lecture des valeurs */
	while(!feof(fp)){ /* Tant qu'on est pas à la fin du fichier */
		if (lire_ligne(fp,line) != 0) {goto LECTURE_FINIE;}
		line_cpt++;
		skip_comment(line); 
		pos = line;
		
	LABEL_TROUVE :
		while(isspace(*pos)) pos++; /* on élimine les espaces */
		while(pos < line+strlen(line)) { /* Tant qu'on est pas à la fin de la ligne */
			tmp = strtod(pos, &endptr); /* on lit le 'double' */
			if (pos == endptr) {goto LECTURE_FINIE;} /* conversion ratée */
			cpt ++;
			pos = endptr;
			while(isspace(*pos)) pos++; /* on élimine les espaces */
		}
	}

	LECTURE_FINIE:
		fclose(fp);
		return cpt;
}

/*---------------------------------------------------------------------------------------------*/
/*! \fn		void lire_ligne(FILE *fp, char *line)
 *
 *  \brief	Lit une ligne dans un fichier et la stocke dans une chaine de charactères 
 */
/*---------------------------------------------------------------------------------------------*/
int lire_ligne(FILE *fp, char *line)
{
	if (fgets(line, SIZE_LINE_BUFFER, fp) == NULL) return 1;
	if (strlen(line) == SIZE_LINE_BUFFER-1) {
		fprintf(stderr, "%s ligne %d : ERREUR, taille de buffer insuffisante,impossible de lire plus de "
						"%d caractères par ligne.\n",__FILE__, __LINE__,SIZE_LINE_BUFFER-1);
		exit(EXIT_FAILURE);
	}
	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		char *label_search(char *str,const char *label)
 *
 *	\brief	Cherche un label dans une chaine de caractères, le label doit être isolé, c.a.d, \n
 *          en début de ligne ou précédé d'un espace au sens de isspace() et suivi d'un espace \n
 *          ou d'un signe '='
 *
 *	\return	La position de la 1ere occurence du label dans la chaine ou NULL si le label n'a pas été trouvé
 */
/*---------------------------------------------------------------------------------------------*/
char *label_search(char *str,const char *label)
{
	char *pos;

	/* Cas du label vide */
	if (strlen(label) == 0) return str;

	/* Recherche du label */
	while((pos=strstr(str,label)) != NULL) {
		/* Vérification que le label est en début de ligne ou précédé par un espace */
		if (pos != str && !isspace(*(pos-1))) {
			str = pos + strlen(label);
			continue;
		}
		/* Vérification que le label est suivi par un espace, saut de ligne ou signe '=' */
		if (strlen(pos) > strlen(label)) {
			if (!isspace(*(pos+strlen(label))) && *(pos+strlen(label)) != '=') {
				str = pos + strlen(label);
				continue;
			}
		}
		break;
	}
	return pos;
}


/*---------------------------------------------------------------------------------------------*/
/*! \fn		void skip_comment(char *str_in_out)
 *
 *  \brief	Elimine tout ce qui se trouve après un commentaire '#' dans str_in_out
 */
/*---------------------------------------------------------------------------------------------*/
void skip_comment(char *str_in_out){

	char *pos;
	/* Cherche CHAR_COMMENT et le remplace par le charactère nul '\0' */
	if((pos = strchr(str_in_out,CHAR_COMMENT)) != NULL) {
		*pos = '\0';
	}
 }


