/*
 *	interp.c
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
	fprintf(stderr,	"interp returns the linearly interpolated value f(x) from an array of x[1st column] f(x)[2nd column] given in a file.\nUsage: interp filename x, where 'filename' points to a file\n");
}

char *label_search(char *str,const char *label);
void skip_comment(char *str_in_out);
int lire_ligne(FILE *fp, char *line);
int lire_tab(const char *nom_fichier, const char *label, double *tab, int N);
int count_tab(const char *nom_fichier, const char *label);
double interp(char* name,double x);

int main(int argc, char *argv[]){
	
	int i;
	double x,f;
	char filename[STR_SIZE], *endptr;
	
	/* Copie des arguments de la ligne de commande */
	char **argvcp; 
	argvcp = (char **) malloc(sizeof(char*)*argc);
	argvcp[0] = (char*) malloc(sizeof(char)*STR_SIZE*argc);
	for(i=1;i<=argc-1;i++){
		argvcp[i] = argvcp[i-1] + STR_SIZE;
		strncpy(argvcp[i], argv[i],STR_SIZE);
	}
	
	/* V�rification de la pr�sence du nombre minimal d'options */
	if (argc < 2){
		err_message();
		return 1;
	}
	
	/* Reading the arguments */
	strncpy(filename, argvcp[1],STR_SIZE);
	x = strtod(argvcp[2], &endptr);
	if (argvcp[2] == endptr){
		fprintf(stderr,"ERROR, %s can't read argument x. Exiting.",__FILE__);
		exit(EXIT_FAILURE);
	}

/*	strncpy(method, argvcp[3],STR_SIZE);*/ /*Ajouter verifications, etc..*/
	
	/* f(x) determination */
	f = interp(filename,x);
	
	fprintf(stdout,"%.16f",f);

	
	return 0;
}

double interp(char* filename,double x)
{
	int i;

	double *tab_x, *tab_f, *table_tmp, f;
	int npoints;
	npoints = ROUND((double) count_tab(filename, "")/2);
	table_tmp= (double *) malloc(sizeof(double)*2*((int)npoints));
	tab_x= (double *) malloc(sizeof(double)*((int)npoints));
	tab_f= (double *) malloc(sizeof(double)*((int)npoints));
	lire_tab(filename, "", table_tmp, 2*(int)npoints);
		/* R�arrangement en plusieurs tableaux */
		for (i=0; i<=(int)npoints -1;i++){
			tab_x[i] = table_tmp[2*i];
			tab_f[i] = table_tmp[2*i+1];
		}
		/* Recherche de l'indice de tableau pour x */
		int num_min = 0;
		int num_max = (int) npoints-1;
		int numero = num_max>>1;
		while(!(tab_x[numero] <= x && x < tab_x[numero+1])){
			numero=num_min+((num_max-num_min)>>1);
			if (x < tab_x[numero]){
				num_max=numero;
			}
			if (tab_x[numero+1] <= x){
				num_min=numero;
			}
			if (num_min == num_max){
				fprintf(stderr,"interp incompatible with x = %f for %s",x,filename);
				exit(EXIT_FAILURE);
			}
		}
		/* Interpolation lin�aire */
		double x1 = tab_x[numero];
		double x2 = tab_x[numero+1];
		double f1 = tab_f[numero];
		double f2 = tab_f[numero+1];
		f= (x-x1)*(f2-f1)/(x2-x1) + f1;

		free(tab_x);
		free(tab_f);
		free(table_tmp);
		
/*printf("x= %f, x1= %f, x2= %f\n",x,tab_x[numero],tab_x[numero+1]);
printf("x= %f, n = %f, k= %f\n",x,index_n,index_k);
*/
		return f;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int lire_tab(char *nom_fichier, const char *label, double *tab, int N)
 *
 *	\brief	Lit N valeurs de format double dans un fichier et les stocke dans un tableau   \n 
 *			Les valeurs doivent �tre s�par�es par un ou plusieurs espaces, tabulations     \n
 *			ou sauts de lignes et pr�c�d�es d'un label �ventuellement suivi d'un signe '='.\n
 *          ex. : (...) tab1 = 3.4  4.5e-3  +46  -7.6e+2 ...                               \n
 *          Remarque : Pour lire un tableau sans label, donner "" comme label.             \n
 *
 *	\return	0 si succ�s, 1 si le nombre d'�l�ments lus diff�re de N ou si le label n'a pas �t� trouv�.
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

	fprintf(stderr,"%s ligne %d : ERREUR, le label '%s' n'a pas �t� trouv� dans %s \n", __FILE__, __LINE__, label, nom_fichier);
	fclose(fp);
	return 1;

	/* Lecture des valeurs */
	while(!feof(fp)){ /* Tant qu'on est pas � la fin du fichier */
		if (lire_ligne(fp,line) != 0) {goto LECTURE_FINIE;}
		line_cpt++;
		skip_comment(line); 
		pos = line;
		
	LABEL_TROUVE :
		while(isspace(*pos)) pos++; /* on �limine les espaces */
		while(pos < line+strlen(line)) { /* Tant qu'on est pas � la fin de la ligne */
			tmp = strtod(pos, &endptr); /* on lit le 'double' */
			if (pos == endptr) {goto LECTURE_FINIE;} /* conversion rat�e */
			if (cpt <= N+1) tab[cpt] = tmp;
			cpt ++;
			pos = endptr;
			while(isspace(*pos)) pos++; /* on �limine les espaces */
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

	fprintf(stderr,"%s ligne %d : ERREUR, le label '%s' n'a pas �t� trouv� dans %s \n", __FILE__, __LINE__, label, nom_fichier);
	fclose(fp);
	return 1;

	/* Lecture des valeurs */
	while(!feof(fp)){ /* Tant qu'on est pas � la fin du fichier */
		if (lire_ligne(fp,line) != 0) {goto LECTURE_FINIE;}
		line_cpt++;
		skip_comment(line); 
		pos = line;
		
	LABEL_TROUVE :
		while(isspace(*pos)) pos++; /* on �limine les espaces */
		while(pos < line+strlen(line)) { /* Tant qu'on est pas � la fin de la ligne */
			tmp = strtod(pos, &endptr); /* on lit le 'double' */
			if (pos == endptr) {goto LECTURE_FINIE;} /* conversion rat�e */
			cpt ++;
			pos = endptr;
			while(isspace(*pos)) pos++; /* on �limine les espaces */
		}
	}

	LECTURE_FINIE:
		fclose(fp);
		return cpt;
}

/*---------------------------------------------------------------------------------------------*/
/*! \fn		void lire_ligne(FILE *fp, char *line)
 *
 *  \brief	Lit une ligne dans un fichier et la stocke dans une chaine de charact�res 
 */
/*---------------------------------------------------------------------------------------------*/
int lire_ligne(FILE *fp, char *line)
{
	if (fgets(line, SIZE_LINE_BUFFER, fp) == NULL) return 1;
	if (strlen(line) == SIZE_LINE_BUFFER-1) {
		fprintf(stderr, "%s ligne %d : ERREUR, taille de buffer insuffisante,impossible de lire plus de "
						"%d caract�res par ligne.\n",__FILE__, __LINE__,SIZE_LINE_BUFFER-1);
		exit(EXIT_FAILURE);
	}
	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		char *label_search(char *str,const char *label)
 *
 *	\brief	Cherche un label dans une chaine de caract�res, le label doit �tre isol�, c.a.d, \n
 *          en d�but de ligne ou pr�c�d� d'un espace au sens de isspace() et suivi d'un espace \n
 *          ou d'un signe '='
 *
 *	\return	La position de la 1ere occurence du label dans la chaine ou NULL si le label n'a pas �t� trouv�
 */
/*---------------------------------------------------------------------------------------------*/
char *label_search(char *str,const char *label)
{
	char *pos;

	/* Cas du label vide */
	if (strlen(label) == 0) return str;

	/* Recherche du label */
	while((pos=strstr(str,label)) != NULL) {
		/* V�rification que le label est en d�but de ligne ou pr�c�d� par un espace */
		if (pos != str && !isspace(*(pos-1))) {
			str = pos + strlen(label);
			continue;
		}
		/* V�rification que le label est suivi par un espace, saut de ligne ou signe '=' */
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
 *  \brief	Elimine tout ce qui se trouve apr�s un commentaire '#' dans str_in_out
 */
/*---------------------------------------------------------------------------------------------*/
void skip_comment(char *str_in_out){

	char *pos;
	/* Cherche CHAR_COMMENT et le remplace par le charact�re nul '\0' */
	if((pos = strchr(str_in_out,CHAR_COMMENT)) != NULL) {
		*pos = '\0';
	}
 }


