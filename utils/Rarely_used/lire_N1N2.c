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

#define CHAR_COMMENT '#'
#define SIZE_LINE_BUFFER 50000
#define STRSIZE 100

static int lire_tab(int N1, int N2);
static int lire_ligne(FILE *fp, char *line);

int main(int argc, char *argv[])
{
	int N1, N2;
	
	/**/
	if (argc<3) {
		fprintf(stderr,"Utilisation : cat file | lire_tab N1 N2\n"
		               "Lit les lignes de N1 à N2\n");
		return -1;
	}
				
	/* Lecture de N1 et N2 en argument n°2 et n°3 */
	N1 = atoi(argv[1]);
	N2 = atoi(argv[2]);
	
	lire_tab(N1, N2);
	
	return 0;
}


static int lire_tab(int N1, int N2)
{
	char line[SIZE_LINE_BUFFER];
	int line_cpt=0;

	/* Lecture des valeurs */
	while(!feof(stdin)){ /* Tant qu'on est pas à la fin du fichier */
		if (lire_ligne(stdin,line) != 0) return 1;
		line_cpt++;
		if (line_cpt >= N1){
			if(line_cpt <= N2){
				fprintf(stdout,"%s",line);
			}else{
				return 0;
			}
		}
	}
	
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


