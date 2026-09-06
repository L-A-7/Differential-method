/*! \file		md_io_utils.c
 *
 *	\brief		Elementary input-output functions
 *
 *  \date		../../2006
 *  \authors	Laurent ARNAUD
 */

#include "md_io_utils.h"


/*---------------------------------------------------------------------------------------------*/
/*! \fn    int lire_int(FILE *fp, char *label, int *value)
 *
 *  \brief	Lit dans le fichier pointé par *fp la valeur entiere 'value' indiquée par 'label' \n
 *			sous la forme label = value (ex.: N2 = 10)
 *  \return	0 si lecture réussie 1 sinon
 */
/*---------------------------------------------------------------------------------------------*/
int lire_int(FILE *fp, const char *label, int *value){

	char stmp[SIZE_LINE_BUFFER];
	char *pos, *pos2;

	rewind(fp);
	while(!feof(fp)){
		if (lire_ligne(fp,stmp) != 0) break;
		skip_comment(stmp); /* Vire les commentaires */
		if ((pos = label_search(stmp,label)) != NULL){ /* Recherche le label */	
			if((pos2 = strchr(pos+strlen(label),'=')) != NULL) *pos2 = ' '; /* remplace '=' par un espace */
			if (sscanf(pos+strlen(label)," %d", value) == 1){ /* Lit la valeur */
				return 0;
			}
		}
	}
	return 1;
}

/*---------------------------------------------------------------------------------------------*/
/*! \fn    int lire_double(FILE *fp, char *label, double *value)
 *
 *  \brief	Lit dans le fichier pointé par *fp la valeur entiere 'value' indiquée par 'label' \n
 *			sous la forme label = value (ex.:  x = 0.12310)
 *  \return	0 si lecture réussie 1 sinon
 */
/*---------------------------------------------------------------------------------------------*/
int lire_double(FILE *fp, const char *label, double *value){

	char *pos, *pos2, stmp[SIZE_LINE_BUFFER], stmp2[SIZE_LINE_BUFFER], *endptr;
	double tmp;

	rewind(fp);
	while(!feof(fp)){
		if (lire_ligne(fp,stmp) != 0) break;
		skip_comment(stmp); /* Vire les commentaires */
		if ((pos = label_search(stmp,label)) != NULL){ /* Recherche le label */	
			if((pos2 = strchr(pos+strlen(label),'=')) != NULL) *pos2 = ' '; /* remplace '=' par un espace */
			if (sscanf(pos+strlen(label)," %s", stmp2) == 1){ /* Lit la valeur */
				tmp = strtod(stmp2, &endptr);
				if (stmp2 != endptr){
					*value = tmp;
					return 0; 
				}
			}
		}
	}
	return 1;
}


/*---------------------------------------------------------------------------------------------*/
/*! \fn    int lire_string(FILE *fp, char *label, char *value)
 *
 *  \brief	Lit dans le fichier pointé par *fp la valeur entiere 'value' indiquée par 'label' \n
 *			sous la forme label = value (ex.: fichier = toto.dat)
 *  \return	0 si lecture réussie 1 sinon
 *	\todo	REMPLACER sscanf("%s") par qqchose de plus sur	
 */
/*---------------------------------------------------------------------------------------------*/
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

/*---------------------------------------------------------------------------------------------*/
/*! \fn    int lire_complex(FILE *fp,  char *label, COMPLEX *value)
 *
 *  \brief	Lit dans le fichier pointé par *fp la valeur entiere 'value' indiquée par 'label' \n
 *			sous la forme label = value (ex.: Z1 = 1.0 + i0.5 )
 *  \return	0 si lecture réussie 1 sinon
 */
/*---------------------------------------------------------------------------------------------*/
int lire_complex(FILE *fp, const char *label, COMPLEX *value){

	char *pos, *pos2, stmp[SIZE_LINE_BUFFER];
	double tmp1, tmp2;

	rewind(fp);
	while(!feof(fp)){
		if (lire_ligne(fp,stmp) != 0) break;
		skip_comment(stmp); /* Vire les commentaires */
		if ((pos = label_search(stmp,label)) != NULL){ /* Recherche le label */	
			if((pos2 = strchr(pos+strlen(label),'=')) != NULL) *pos2 = ' '; /* remplace '=' par un espace */
			if (sscanf(pos+strlen(label)," %lf + i%lf", &tmp1, &tmp2) == 2){
				*value = tmp1 + I*tmp2;
				return 0;
			}else if (sscanf(pos+strlen(label)," %lf+i%lf", &tmp1, &tmp2) == 2){
				*value = tmp1 + I*tmp2;
				return 0;
			}
		}
	}
	return 1;
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

	fprintf(stderr,"%s line %d: ERROR, label '%s' was not found in %s \n", __FILE__, __LINE__, label, nom_fichier);
	fclose(fp);
	return 1;

	/* Lecture des valeurs */
	while(!feof(fp)){
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
		fprintf(stderr, "%s line %d: ERROR, %s contains %d elements instead of %d in %s (line %d)\n", __FILE__, __LINE__, label, cpt, N, nom_fichier, line_cpt);
		fclose(fp);
		return 2;
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
		fprintf(stderr, "%s line %d: ERROR, insufficient buffer size. Can't read more than %d characters per line.\n",__FILE__, __LINE__,SIZE_LINE_BUFFER-1);
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


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int	ecrire_dble_tab(FILE *fp, double *tab, int N, char *separateur1, int Nmax1, char *separateur2)
 *
 *	\brief	Ecrit les valeurs d'un tableau séparées par les 'séparateurs1' (par ex " "), plus par les \n
 *          'séparateurs2' (par ex "\n") une fois tous les Nmax1 éléments.
 */
/*---------------------------------------------------------------------------------------------*/
int ecrire_dble_tab(FILE *fp, double *tab, int N, char *separateur1, int Nmax1, char *separateur2)
{
	int i,k=0;

	while((k+1)*Nmax1 < N) {
		for (i=k*Nmax1; i<=MIN(N-1,(k+1)*Nmax1-1); i++){
			fprintf(fp,"% 1.12e%s",tab[i],separateur1);
		}
		k++;
		fprintf(fp,"%s",separateur2);
	}
	/* Derniere ligne, traitée à part car pas de séparateur2 à la fin ! */
	for (i=k*Nmax1; i<=MIN(N-1,(k+1)*Nmax1-1); i++){
		fprintf(fp,"% 1.12e%s",tab[i],separateur1);
	}
	
	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int	ecrire_cplx_tab(FILE *fp, COMPLEX *tab, int N, int mode, char *separateur1, int Nmax1, char *separateur2)
 *
 *	\brief	Ecrit les valeurs d'un tableau séparées par les 'séparateurs1' (par ex " "), plus par les \n
 *          'séparateurs2' (par ex "\n") une fois tous les Nmax1 éléments, mode RE/IM
 */
/*---------------------------------------------------------------------------------------------*/
int ecrire_cplx_tab(FILE *fp, COMPLEX *tab, int N, int mode, char *separateur1, int Nmax1, char *separateur2)
{
	int i,k=0;
	if (mode == RE){
		while((k+1)*Nmax1 < N) {
			for (i=k*Nmax1; i<=MIN(N-1,(k+1)*Nmax1-1); i++){
				fprintf(fp,"% 1.12e%s",creal(tab[i]),separateur1);
			}
			k++;
			fprintf(fp,"%s",separateur2);
		}
		/* Derniere ligne, traitée à part car pas de séparateur2 à la fin ! */
		for (i=k*Nmax1; i<=MIN(N-1,(k+1)*Nmax1-1); i++){
			fprintf(fp,"% 1.12e%s",creal(tab[i]),separateur1);
		}
	}else	if (mode == IM){
		while((k+1)*Nmax1 < N) {
			for (i=k*Nmax1; i<=MIN(N-1,(k+1)*Nmax1-1); i++){
				fprintf(fp,"% 1.12e%s",cimag(tab[i]),separateur1);
			}
			k++;
			fprintf(fp,"%s",separateur2);
		}
		/* Derniere ligne, traitée à part car pas de séparateur2 à la fin ! */
		for (i=k*Nmax1; i<=MIN(N-1,(k+1)*Nmax1-1); i++){
			fprintf(fp,"% 1.12e%s",cimag(tab[i]),separateur1);
		}
	}else{
		fprintf(stderr,"ERREUR, ecrire_cplx_tab, mode non reconnu");	
	}
	
	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int lire_str_arg(char *dest, char *label, int argc, char **argvcp)
 *
 *	\brief	Lit la valeur de l'argument de la ligne de commande indiqué sous la forme "-label valeur"
 */
/*---------------------------------------------------------------------------------------------*/
int lire_str_arg(char *dest, char *label, int argc, char **argvcp)
{
	int i;

	for(i=1;i<=argc-2;i++){
		if (!strcmp(argvcp[i],label)){
			strncpy(dest, argvcp[i+1], SIZE_STR_BUFFER);
			return 0; 
		}
	}
	return 1;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int lire_dble_arg(double *res, char *label, int argc, char **argvcp)
 *
 *	\brief	Lit la valeur de l'argument de la ligne de commande indiqué sous la forme "-label valeur"
 */
/*---------------------------------------------------------------------------------------------*/
int lire_dble_arg(double *res, char *label, int argc, char **argvcp)
{
	int i;
	char strtmp[SIZE_STR_BUFFER], *endptr;
	double tmp;

	for(i=1;i<=argc-2;i++){
		if (!strcmp(argvcp[i],label)){
			strncpy(strtmp, argvcp[i+1], SIZE_STR_BUFFER);
			tmp = strtod(strtmp, &endptr);
			if (strtmp != endptr){
				*res = tmp;
				return 0; 
			}
		}
	}
	return 1;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn	int lire_int_arg(int *res, char *label, int argc, char **argvcp)
 *
 *	\brief	Lit la valeur de l'argument de la ligne de commande indiqué sous la forme "-label valeur"
 */
/*---------------------------------------------------------------------------------------------*/
int lire_int_arg(int *res, char *label, int argc, char **argvcp)
{
	int i, tmp;
	char strtmp[SIZE_STR_BUFFER], *endptr;

	for(i=1;i<=argc-2;i++){
		if (!strcmp(argvcp[i],label)){
			strncpy(strtmp, argvcp[i+1], SIZE_STR_BUFFER);
			tmp = (int) strtod(strtmp, &endptr);
			if (strtmp != endptr){
				*res = tmp;
				return 0;
			}
		}
	}
	return 1;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int lire_complex_arg(COMPLEX *res, char *label, int argc, char **argvcp)
 *
 *	\brief	Read a complex argument in the command line. The format must be re+iIm, eg. 1.5+i0.1 with no space
 */
/*---------------------------------------------------------------------------------------------*/
int lire_complex_arg(COMPLEX *res, char *label, int argc, char **argvcp)
{
	int i;
	char strtmp[SIZE_STR_BUFFER];
	double Re, Im;

	for(i=1;i<=argc-2;i++){
		if (!strcmp(argvcp[i],label)){
			strncpy(strtmp, argvcp[i+1], SIZE_STR_BUFFER);
			if (sscanf(strtmp, "%lf+i%lf",&Re,&Im) == 2){
				*res= Re +I*Im;
				return 0;
			}
		}
	}
	return 1;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int arg_read(int argc, char **argvcp, FILE *fp, char *type, void *var, char *flag, int exit_or_not)
 *	
 *	\brief	Lecture des paramètres dans un fichier
 *
 *---------------------------------------------------------------------------------------------*/
int arg_read(int argc, char **argvcp, FILE *fp, char *type, void *var, char *flag, int exit_or_not)
{
	int read_error=0;
	char cmd_line_flag[SIZE_STR_BUFFER];
	snprintf(cmd_line_flag, SIZE_STR_BUFFER*sizeof(char), "-%s",flag);

	if (!strcmp(type,"double")){
		if (lire_dble_arg((double *) var, cmd_line_flag, argc, argvcp)) {
			if (lire_double (fp, flag, (double *) var)) read_error = 1;}
	}else if (!strcmp(type,"int")){
		if (lire_int_arg((int *) var, cmd_line_flag, argc, argvcp)) {
			if (lire_int (fp, flag, (int *) var)) read_error = 1;}
	}else if (!strcmp(type,"complex")){
		if (lire_complex_arg((COMPLEX *) var, cmd_line_flag, argc, argvcp)) {
			if (lire_complex (fp, flag, (COMPLEX *) var)) read_error = 1;}
	}else if (!strcmp(type,"string")){
		if (lire_str_arg((char *) var, cmd_line_flag, argc, argvcp)) {
			if (lire_string (fp, flag, (char *) var)) read_error = 1;}
	}else{
		fprintf(stderr, "%s line %d: ERROR, no type \"%s\"\n",__FILE__,__LINE__,type);
		if (exit_or_not==EXIT_ON_ERROR) exit(EXIT_FAILURE);
		return -2;
	}

	if (read_error != 0){
		if (exit_or_not == EXIT_ON_ERROR){
			fprintf(stderr, "%s line %d: ERROR, can't read \"%s\"\n",__FILE__,__LINE__,flag);
			exit(EXIT_FAILURE);
			return -1;
		}else{
			return 1;
		}
	}
	
	return 0;
}


#if 0
/*PAS FINIE, pas utile pour l'insant*/
int ecrire_col(double *tab, char *nomtab, char *nom_fichier) 
{
	char line[SIZE_LINE_BUFFER];
	FILE *fp, *fp_tmp;
	
	/* Ouverture du fichier */	
	if (!(fp = fopen(nom_fichier,"w+"))){
		fprintf(stderr, "%s ligne %d : ERREUR, impossible d'ouvrir %s\n",__FILE__, __LINE__,nom_fichier);
		return 1;
	}
	
	/* Copie dans un fichier temporaire */
	char nom_fichier_tmp[] = "md3D_fichier_tmp_68gIg78GUgkd.tmp";
	if (!(fp_tmp = fopen(nom_fichier_tmp,"w+"))){
		fprintf(stderr, "%s ligne %d : ERREUR, impossible d'ouvrir %s\n",__FILE__, __LINE__,nom_fichier_tmp);
		return 1;
	}

	
	/* Comptage du nombre de caractères de la plus longue ligne */
	int max = 0;
	while(!feof(fp)){ 
		if (lire_ligne(fp,line) != 0) break;
		max = MAX(max,strlen(line));
	}

	/* Fermeture des fichiers */

	return 0;
}
#endif
