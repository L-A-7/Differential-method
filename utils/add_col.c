/*	lire_tab
 *
 *	Lit les valeurs de format double dans un fichier et les affiche à l'écran.      
 *	Les valeurs doivent être séparées par un ou plusieurs espaces, tabulations     
 *	ou sauts de lignes et précédées d'un label éventuellement suivi d'un signe '='.
 *  ex. : (...) tab1 = 3.4  4.5e-3  +46  -7.6e+2 ...                               
 *
 *	Utilisation :
 *		cat toto.txt | lire_tab label
 *	ou	lire_tab label < toto.txt
 *
 */




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX(a,b) ((a>b)?a:b)
#define SIZE_LINE_BUFFER 5000
#define STRSIZE 100


static int lire_ligne(FILE *fp, char *line);
static void skip_comment(char *str_in_out);

int main(int argc, char *argv[])
{
		
	/**/
	if (argc<2) {
		fprintf(stderr,"Utilisation : cat file | lire_tab N1 N2\n"
		               "Lit les lignes de N1 à N2\n");
		return -1;
	}
				
	/* Lecture du nom du fichier en argument n°1 */
	char nom_fichier[STRSIZE];
	strncpy(nom_fichier,argv[1],STRSIZE);

	char line[SIZE_LINE_BUFFER];
	char line2[SIZE_LINE_BUFFER];
	
	FILE *fp;
	
	/* Ouverture du fichier */	
	if (!(fp = fopen(nom_fichier,"r"))){
		fprintf(stderr, "%s ligne %d : ERREUR, impossible d'ouvrir %s\n",__FILE__, __LINE__,nom_fichier);
		return 1;
	}
	
	int max = 0;
	while(!feof(fp)){ 
		if (lire_ligne(fp,line) != 0)     sprintf(line," ");
		if (lire_ligne(stdin,line2) != 0) sprintf(line2," ");
		if (!strcmp(line," ") && !strcmp(line," ")) break;
		skip_comment(line);
		skip_comment(line2);
		fprintf(stdout,"%s\t%s\n",line,line2);
	}
	
	/* Fermeture des fichiers */
	fclose(fp);
	
	return 0;
}
/**********************************************/


/*! \fn		static void lire_ligne(FILE *fp, char *line)
 *
 *  \brief	Lit une ligne dans un fichier et la stocke dans une chaine de charactères 
 */
static int lire_ligne(FILE *fp, char *line)
{
	if (fgets(line, SIZE_LINE_BUFFER, fp) == NULL) return 1;
	if (strlen(line) == SIZE_LINE_BUFFER-1) {
		fprintf(stderr, "%s ligne %d : ERREUR, taille de buffer insuffisante,impossible de lire plus de "
						"%d caractères par ligne.\n",__FILE__, __LINE__,SIZE_LINE_BUFFER-1);
		exit(EXIT_FAILURE);
	}
	return 0;
}


static void skip_comment(char *str_in_out)
{

	char *pos;
	/* Cherche CHAR_COMMENT et le remplace par le charactère nul '\0' */
	if((pos = strchr(str_in_out,'\n')) != NULL) {
		*pos = '\0';
	}
 }
