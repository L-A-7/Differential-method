/*
 *	sub_tanpsicosdelta.c
 *
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <complex.h>

#define STR_SIZE 5000
#define SIZE_STR_BUFFER 200
#define SIZE_LINE_BUFFER 50000
#define CHAR_COMMENT '#'

void err_message(){
	fprintf(stderr,	"usage : sub_tanpsicosdelta Nb_lambda SIMU_FILE_NAME EXP_FILE_NAME\n");
}
char *label_search(char *str,const char *label);
void skip_comment(char *str_in_out);
int lire_ligne(FILE *fp, char *line);
int lire_tab(const char *nom_fichier, const char *label, double *tab, int N);
int SaveDbleTab2file (double *tab, int N, char *filename, char *separateur);

int main(int argc, char *argv[]){
	
	int i,nl,Nlambda;
	double eps, *lambda_simu, *lambda_exp, *tanPsi_simu, *tanPsi_exp, *cosDelta_simu, *cosDelta_exp, *diff;
	char simu_file[STR_SIZE], exp_file[STR_SIZE], *endptr;
	
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
	Nlambda = (int) strtod(argvcp[1], &endptr);
	if (argvcp[1] == endptr){
		fprintf(stderr,"ERROR, %s can't read argument Nlambda. Exiting.",__FILE__);
		exit(EXIT_FAILURE);
	}
	strncpy(simu_file, argvcp[2],STR_SIZE);
	strncpy(exp_file, argvcp[3],STR_SIZE);
	
	/* Memory allocation */
	lambda_simu = (double *) malloc(sizeof(double)*Nlambda);
	lambda_exp = (double *) malloc(sizeof(double)*Nlambda);
	tanPsi_simu = (double *) malloc(sizeof(double)*Nlambda);
	tanPsi_exp = (double *) malloc(sizeof(double)*Nlambda);
	cosDelta_simu = (double *) malloc(sizeof(double)*Nlambda);
	cosDelta_exp = (double *) malloc(sizeof(double)*Nlambda);
	diff = (double *) malloc(sizeof(double)*2*Nlambda);

	/* Reading the data */
	lire_tab(simu_file, "lambda", lambda_simu, Nlambda);	
	lire_tab( exp_file, "lambda", lambda_exp , Nlambda);	
	lire_tab(simu_file, "tan_Psi", tanPsi_simu, Nlambda);	
	lire_tab( exp_file, "tan_Psi", tanPsi_exp , Nlambda);	
	lire_tab(simu_file, "cos_Delta", cosDelta_simu, Nlambda);	
	lire_tab( exp_file, "cos_Delta", cosDelta_exp , Nlambda);	

	/* Checking that the exp & simu values of lambda are the same */
	eps = 1e-5;
	for (nl=0;nl<=Nlambda-1;nl++){
		if (fabs(lambda_simu[nl]-lambda_exp[nl]) > eps){
			fprintf(stderr,"ERROR, %s, lambda_exp and lambda_simu don't have the same values\n",__FILE__);
			exit(EXIT_FAILURE);
		}
	}
	
	/* Calculating the difference between tanPsi + cosDelta exp and simulated */
	for (nl=0;nl<=Nlambda-1;nl++){
		diff[nl] = tanPsi_exp[nl] - tanPsi_simu[nl];
		diff[nl+Nlambda] = cosDelta_exp[nl] - cosDelta_simu[nl];
	}
	SaveDbleTab2file (diff, 2*Nlambda, "stdout", " ");
	
	free(lambda_simu);
	free(lambda_exp);
	free(tanPsi_simu);
	free(tanPsi_exp);
	free(cosDelta_simu);
	free(cosDelta_exp);
	free(diff);
	
	return 0;
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
		fprintf(stderr, "%s ligne %d : ERREUR, impossible d'ouvrir %s\n",__FILE__, __LINE__,nom_fichier);
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
		return 1;
	}
	fclose(fp);
	return 0;
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


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int SaveDbleTab2file (char *filename, double *tab, int N, char *separateur)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
int SaveDbleTab2file (double *tab, int N, char *filename, char *separateur)
{
  int i;
  FILE *fp;
	
  /* Affichage à l'écran si filename = "stdout" */										
  if(!strcmp(filename,"stdout")) fp = stdout;
  else fp = fopen(filename, "w");
	
  for (i=0;i<=N-1;i++){
    fprintf(fp,"%1.6e%s",tab[i],separateur);
  }

  if(strcmp(filename,"stdout")) fclose(fp);

  return 0;
}
