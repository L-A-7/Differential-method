/*	lire_ligne
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

#define CHAR_COMMENT '#'
#define SIZE_LINE_BUFFER 50000

static int lire_ligne(FILE *fp, char *line);

int main(int argc, char *argv[])
{
	char line[SIZE_LINE_BUFFER];
	int i,line_N;
		

	/**/
	if (argc<2) {
		fprintf(stderr,"Utilisation : cat file | lire_ligne ligne_n°\n");
		return -1;
	}
			
	/* Lecture de label en argument n°1 */
	line_N = atoi(argv[1]);

	for (i=1;i<=line_N;i++)	{
		lire_ligne(stdin,line);
	}
	
	fprintf(stdout,"%s",line);
	return 0;
}

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

