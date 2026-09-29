#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define CHAR_COMMENT '#'
#define SIZE_LINE_BUFFER 5000

int lire_string(FILE *fp, const char *label, char *value);
static int lire_ligne(FILE *fp, char *line);
static void skip_comment(char *str_in_out);
char *label_search(char *str, const char *label);

int main(int argc, char *argv[])
{
	FILE *fp;
	char result[100];
	
	/* Lecture des options de la ligne de commande */
	if (argc < 3) {
		fprintf(stdout,"Utilisation : lire_string filename label\n");
		return 1;
	}
	char *filename = argv[1];
	char *label    = argv[2];

	/* Ouverture du fichier */
	if (!(fp = fopen(filename,"r"))){
		fprintf(stderr, "%s ligne %d : ERREUR, impossible d'ouvrir %s\n",__FILE__, __LINE__,filename);
		return 1;
	}
	
	/* lecture du label */
	if (lire_string(fp, label, result) != 0) return 1;
	
	/* Ecriture du résultat */
	fprintf(stdout,"%s",result);
	
	fclose(fp);
	return 0;
}








/*! \fn    int lire_string(FILE *fp, char *label, char *value)
 *
 *  \brief	Lit dans le fichier pointé par *fp la valeur entiere 'value' indiquée par 'label' \n
 *			sous la forme label = value (ex.: fichier = toto.dat)
 *  \return	0 si lecture réussie 1 sinon
 *	\todo	REMPLACER sscanf("%s") par qqchose de plus sur	
 */
int lire_string(FILE *fp, const char *label, char *value){

	char *pos, *pos2, stmp[SIZE_LINE_BUFFER];

	rewind(fp);
	while(!feof(fp)){
		if (lire_ligne(fp,stmp) != 0) break;
		skip_comment(stmp); /* Vire les commentaires */
		if ((pos = label_search(stmp,label)) != NULL){ /* Recherche le label */	
			if((pos2 = strchr(pos+strlen(label),'=')) != NULL) *pos2 = ' '; /* remplace '=' par un espace */
			if (sscanf(pos+strlen(label)," %s", value) > 0){ /* Lit la valeur */
				return 0;
			}
		}
	}
	return 1;
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



/*! \fn		void lire_ligne(FILE *fp, char *line)
 *
 *  \brief	Lit une ligne dans un fichier et la stocke dans une chaine de charactères 
 */
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


/*! \fn		void skip_comment(char *str_in_out)
 *
 *  \brief	Elimine tout ce qui se trouve après un commentaire '#' dans str_in_out
 */
void skip_comment(char *str_in_out){

	char *pos;
	/* Cherche CHAR_COMMENT et le remplace par le charactère nul '\0' */
	if((pos = strchr(str_in_out,CHAR_COMMENT)) != NULL) {
		*pos = '\0';
	}
 }
