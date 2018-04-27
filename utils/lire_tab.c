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

static int lire_tab(const char *label);
static int lire_ligne(FILE *fp, char *line);
static void skip_comment(char *str_in_out);
char *label_search(char *str, const char *label);

int main(int argc, char *argv[])
{
	/**/
	if (argc<2) {
		fprintf(stderr,"Utilisation : cat file | lire_tab label\n");
		return -1;
	}
			
	/* Lecture de label en argument n°1 */
	char label[1000];
	strncpy(label, argv[1], 1000);
	
	lire_tab(label);
	
	return 0;
}


static int lire_tab(const char *label)
{
	char *pos, *pos2, *endptr, line[SIZE_LINE_BUFFER];
	int cpt=0, line_cpt=0;
	double tmp;
	
	/* Recherche du label */
	while(!feof(stdin)){ 
		if (lire_ligne(stdin,line) != 0) goto LECTURE_FINIE;
		line_cpt++;
		skip_comment(line);
		if ((pos = label_search(line,label)) != NULL){ /* on cherche le label */	
			pos += strlen(label); 
			if((pos2 = strchr(pos,'=')) != NULL) *pos2 = ' '; /* on remplace '=' par ' ' */
			goto LABEL_TROUVE;
		}
	}

	fprintf(stderr,"%s ligne %d : ERREUR, le label '%s' n'a pas été trouvé\n", __FILE__, __LINE__, label);
	return 1;

	/* Lecture des valeurs */
	while(!feof(stdin)){ /* Tant qu'on est pas à la fin du fichier */
		if (lire_ligne(stdin,line) != 0) goto LECTURE_FINIE;
		line_cpt++;
		skip_comment(line); 
		pos = line;
		
LABEL_TROUVE :
		while(isspace(*pos)) pos++; /* on élimine les espaces */
		while(pos < line+strlen(line)) { /* Tant qu'on est pas à la fin de la ligne */
			tmp = strtod(pos, &endptr); /* on lit le 'double' */
			if (pos == endptr) goto LECTURE_FINIE; /* conversion ratée */
			fprintf(stdout,"% 1.12e\n",tmp);
			cpt ++;
			pos = endptr;
			while(isspace(*pos)) pos++; /* on élimine les espaces */
		}
	}

LECTURE_FINIE:
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


/*! \fn		static void skip_comment(char *str_in_out)
 *
 *  \brief	Elimine tout ce qui se trouve après un commentaire '#' dans str_in_out
 */
static void skip_comment(char *str_in_out){

	char *pos;
	/* Cherche CHAR_COMMENT et le remplace par le charactère nul '\0' */
	if((pos = strchr(str_in_out,CHAR_COMMENT)) != NULL) {
		*pos = '\0';
	}
 }


/*!	\fn		char *label_search(char *str,const char *label)
 *
 *	\brief	Cherche un label dans une chaine de caractères, le label doit être isolé, c.a.d, \n
 *          en début de ligne ou précédé d'un espace au sens de isspace() et suivi d'un espace \n
 *          ou d'un signe '='
 *
 *	\return	La position de la 1ere occurence du label dans la chaine ou NULL si le label n'a pas été trouvé
 */
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
