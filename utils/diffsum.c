/*	diffsum
 *
 *	Calcule	la somme des différences absolues entre les éléments d'un ensemble de paires
 *	[A1 A2 B1 B2 C1 C2 ...] --> |A1-A2|+|B1-B2|+|C1-C2|+ ...
 *
 *	Utilisation :
 *		cat tab.txt | diffsum
 *	ou	diffsum < tab.txt
 *
 */




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#define CHAR_COMMENT '#'
#define SIZE_LINE_BUFFER 50000

static int lire_ligne(FILE *fp, char *line);
static void skip_comment(char *str_in_out);
char *label_search(char *str, const char *label);
static int lire_tab(const char *label, double *som_diff, double *som_A1, double *som_A2);


int main(int argc, char *argv[])
{
	double som_diff, som_A1, som_A2;
	
	lire_tab("", &som_diff, &som_A1, &som_A2);
	printf("%f\n", som_diff*100.0);
/*	printf("som A1 = %f\n", som_A1);
	printf("som A2 = %f\n", som_A2);
*/	
	
	return 0;
}



static int lire_tab(const char *label, double *som_diff, double *som_A1, double *som_A2)
{
	char *pos, *pos2, *endptr, line[SIZE_LINE_BUFFER];
	int cpt=0, line_cpt=0, num = 0;
	double A[2], somA[2], tmp, somdiff = 0;
	somA[0] = 0;
	somA[1] = 0;
	
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
			
			A[num] = tmp;
			somA[num] += A[num];
			if (num == 1){
				somdiff += fabs(A[0]-A[1]);
			}
			num = 1 - num; /* num alterne entre 0 et 1 */
		
			pos = endptr;
			while(isspace(*pos)) pos++; /* on élimine les espaces */
		}
	}

LECTURE_FINIE:
	*som_diff = somdiff;
	*som_A1 = somA[0];
	*som_A2 = somA[1];
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
