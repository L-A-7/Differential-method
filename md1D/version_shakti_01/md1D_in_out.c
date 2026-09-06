/*! \file		md1D_in_out.c
 *
 *	\brief		Routines de gestion des entrées-sorties pour md1D
 *
 *	\version	0.1
 *  \date		../../2004
 *  \authors	Laurent ARNAUD
 */

#include "md1D_in_out.h"

/*!	\fn		int md1D_lire_config(const char *nom_fichier, char *fichier_param)
 *	
 *	\brief	Lecture de la configuration
 */
int md1D_lire_config(const char *nom_fichier, char *fichier_param){

	FILE *fp;
	
	/* Ouverture du fichier */
	if (!(fp = fopen(nom_fichier,"r"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nom_fichier);
		exit(EXIT_FAILURE);
	}

	/* Lecture des paramètres */
	if (lire_string (fp, "fichier_param", fichier_param)){
		fprintf(stderr, "%s ligne %d : Erreur, impossible de lire \"fichier_param\" dans %s\n",__FILE__, __LINE__,nom_fichier);
		exit(EXIT_FAILURE);
	}

	/* fermeture du fichier */
	fclose(fp);

	return 0;
}




/*!	\fn		int md1D_lire_param(const char *nom_fichier, char *fichier_profil, struct Param_struct *par){
 *	
 *	\brief	Lecture des paramètres dans un fichier
 */
int md1D_lire_param(const char *nom_fichier, char *fichier_profil, struct Param_struct *par){

	FILE *fp;
	int ret=0;
	char *erreur="NO_ERROR                     ";
	char str_profil[50], str_pola[10], str_type_calcul[50];

	

	/* Lecture des grandeurs dans fichier_param */
	if (par->verbose) fprintf(stdout,"Lecture des paramètres dans %s\n",nom_fichier);
	if (!(fp = fopen(nom_fichier,"r"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nom_fichier);
		exit(EXIT_FAILURE);
	}
	if (lire_string (fp, "polarisation", str_pola)) erreur="polarisation";
		if      (!strcmp(str_pola,"TE")){par->pola = TE;}
		else if (!strcmp(str_pola,"TM")){par->pola = TM;}
		else                            {erreur = "polarisation";}	
	if (lire_string (fp, "type_calcul", str_type_calcul)) erreur="type_calcul";
		if      (!strcmp(str_type_calcul,"STD"))      {par->type_calcul = STD;}
		else if (!strcmp(str_type_calcul,"VAR_I"))    {par->type_calcul =VAR_I;}
		else if (!strcmp(str_type_calcul,"VAR_I_BIS")){par->type_calcul =VAR_I_BIS;}
		else                                   {erreur = "type_calcul";}
	if (lire_complex(fp, "n_super", &(par->n_super))) erreur="n_super";
	if (lire_complex(fp, "n_sub", &(par->n_sub) )) erreur="n_sub";
	if (lire_double (fp, "L", &(par->L) )) erreur="L";
	if (lire_double (fp, "h", &(par->h) )) erreur="h";
	if (lire_double (fp, "lambda", &(par->lambda) )) erreur="lambda";
	if (lire_double (fp, "angle_i", &(par->angle_i) )) erreur="angle_i";
		par->angle_i *= PI/180.0;
	if (lire_int    (fp, "N", &(par->N) )) erreur="N";
	if (lire_int    (fp, "NS", &(par->NS) )) erreur="NS";
	if (lire_double (fp, "delta_h_approx", &(par->delta_h_approx) )) erreur="delta_h_approx";
	if (lire_string (fp, "nom_profil", par->nom_profil)) erreur="nom_profil";
	if (!strcmp(fichier_profil, NON_LU)){
		if (lire_string (fp, "fichier_profil", fichier_profil)) erreur="fichier_profil";}
	fclose(fp);
	/* Vérification de l'absence d'erreurs de lecture */
	if (strcmp(erreur,"NO_ERROR                     ")){
		fprintf(stderr, "%s : Erreur, probleme de lecture de \"%s\"\n",__FILE__,erreur);
		exit(EXIT_FAILURE);
	}

	/* Lecture des paramètres dans fichier profil */
	if (par->verbose) fprintf(stdout,"Lecture des paramètres dans %s\n",fichier_profil);
	if (!(fp = fopen(fichier_profil,"r"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,fichier_profil);
		exit(EXIT_FAILURE);
	}
	erreur="NO_ERROR                     ";
	if (lire_string (fp, "type_profil", str_profil)) {
		erreur="type_profil";
	}else{
		if       (!strcmp(str_profil,"H_X"))         {par->type_profil = H_X;
			par->N_layers = 0;
		}else if (!strcmp(str_profil,"MULTICOUCHES")){par->type_profil = MULTICOUCHES;
			if (lire_int (fp, "N_layers", &(par->N_layers) )) erreur="N_layers";
		}else if (!strcmp(str_profil,"N_XYZ"))       {par->type_profil = N_XYZ;
		}else                                        {erreur = "type_profil_bis";}
	}
	
	if (lire_int (fp, "N_profil", &(par->N_profil) ))        erreur="N_profil";
	
	fclose(fp);
	/* Vérification de l'absence d'erreurs de lecture */
	if (strcmp(erreur,"NO_ERROR                     ")){
		fprintf(stderr, "%s : Erreur, probleme de lecture de \"%s\"\n",__FILE__,erreur);
		exit(EXIT_FAILURE);
	}

	return ret;
}

int md1D_affiche_valeurs_lues(struct Param_struct *par, struct Noms_fichiers *nomfichier)
{
printf("Eh Oh ! Je suis là !!!\npar->verbose=%d",par->verbose);

	/* Affichage des valeurs lues */
	if (par->verbose){
		fprintf(stdout,"n_super = %f + i%f\n", creal((par->n_super)), cimag((par->n_super)));
		fprintf(stdout,"n_sub   = %f + i%f\n", creal((par->n_sub)), cimag((par->n_sub)));
		fprintf(stdout,"L       = %f\n",par->L);
		fprintf(stdout,"h       = %f\n",par->h);
		fprintf(stdout,"lambda  = %f\n",par->lambda);
		fprintf(stdout,"angle_i = %f rad (%f°)\n",par->angle_i,par->angle_i*180.0/PI);
		fprintf(stdout,"N       = %d\n",par->N);
		fprintf(stdout,"NS      = %d\n",par->NS);
		fprintf(stdout,"Polarisation   : %s\n",(par->pola==TE ? "TE" : "TM"));
		fprintf(stdout,"type_calcul    = %s\n",(par->type_calcul==STD ? "STD" :
		                                        (par->type_calcul==VAR_I ? "VAR_I":"VAR_I_BIS")));

		fprintf(stdout,"delta_h_approx = %f\n",par->delta_h_approx);
		fprintf(stdout,"nom_profil     = %s\n",par->nom_profil);
		fprintf(stdout,"fichier_profil = %s\n",nomfichier->fichier_profil);
		fprintf(stdout,"N_x            = %d\n",par->N_profil);
		fprintf(stdout,"type_profil    = %s\n",(par->type_profil==H_X ? "H_X" :
		                                        (par->type_profil==N_XYZ ? "N_XYZ":"MULTICOUCHES")));
		fprintf(stdout,"N_couches      = %d\n",par->N_layers);
		fflush(stdout);
	}
	return 0;
}



/*!	\fn	int md1D_lire_profil_H_X(const char *nom_fichier, struct Param_struct *par)
 *
 *	\brief	Fonction lisant les valeurs décrivant un profil h(x) dans un fichier. \n
 *		Les valeurs stockées dans le fichier doivent varier entre 0 et 1,     \n
 *		les valeurs lues sont multipliées par h, pour avoir un profil variant \n
 *		entre 0 et h.
 */
int md1D_lire_profil_H_X(const char *nom_fichier, struct Param_struct *par)
{

	int i;
	int N_profil = par->N_profil;
	double h = par->h;
	double *profil = par->profil[0];

	/* Entrée des k2 et invk2 du sub et du super dans les tableaux (inv)k2_layer, pour compatibilité avec multicouches */
	par->k2_layer[0] = (par->n_super*par->k0)*(par->n_super*par->k0);
	par->k2_layer[par->N_layers+1] = (par->n_sub*par->k0)*(par->n_sub*par->k0);
	par->invk2_layer[0] = 1/((par->n_super*par->k0)*(par->n_super*par->k0));
	par->invk2_layer[par->N_layers+1] = 1/((par->n_sub*par->k0)*(par->n_sub*par->k0));


	/* Lecture du profil */
	if (par->verbose) fprintf(stdout,"Lecture du profil %s : ",nom_fichier);
	if (lire_tab(nom_fichier, "profil", profil, N_profil) == 0) {
		if (par->verbose) fprintf(stdout,"OK\n");
	}else{
		fprintf(stderr,"ERREUR de lecture du profil\n");
		exit(EXIT_FAILURE);
	}
	
/***************//***************/
/* Calcul de h */
/*printf("*************************\nATTENTION RECALCUL DE h\n*************************\n");
double mmin = 0, mmax=0;
for (i=0;i<=par->N_profil-1;i++){
	mmin = MIN(mmin,profil[i]);
	mmax = MAX(mmax,profil[i]);
}
par->h = mmax - mmin;
printf("*** h = %f ***\n\n",par->h);
*/
/***************//***************/		

	/* Vérification que le profil varie dans l'intervale [0 1] */
	double max = profil[0];
	double min = profil[0];
	for (i=1; i<=N_profil-1; i++) {
		max = MAX(max, profil[i]);
		min = MIN(min, profil[i]);
	}
	double eps = 1.0e-10; 
	/* Si profil compris dans [0 1], avertissement seulement */
	if ((max<1-eps && min>=-eps) || (max<=1+eps && min>eps)) {
		fprintf(stderr,"ATTENTION, le profil %s varie dans l'intervale [%f %f] et non [0 1] ! \n",nom_fichier,min,max);
	}
	/* Si profil non compris dans [0 1], avertissement et renormalisation de 0 à 1 */
	if (max > 1+eps || min < -eps) {
		fprintf(stderr,"ATTENTION, le profil %s varie dans l'intervale [%f %f] et non [0 1] ! \n",nom_fichier,min,max);
		fprintf(stderr,"=> Normalisation entre 0 et 1\n");
		for (i=0; i<=N_profil-1; i++) {
			profil[i] = (profil[i]-min)/(max-min);
		}
	}
		
	/* Normalisation entre 0 et h */
	for (i=0; i<=N_profil-1; i++) {
		profil[i] *= h;
	}
	if (par->verbose) fprintf(stdout,"Normalisation du profil entre 0 et h : OK\n");


	return 0;
}


/*!	\fn		int md1D_lire_profil_MULTI(const char *nom_fichier, struct Param_struct *par)

 *
 *	\brief	Fonction lisant les valeurs décrivant un profil h(x) dans un fichier. \n
 *		Les valeurs stockées dans le fichier doivent varier entre 0 et 1,     \n
 *		les valeurs lues sont multipliées par h, pour avoir un profil variant \n
 *		entre 0 et h.
 */
int md1D_lire_profil_MULTI(const char *nom_fichier, struct Param_struct *par)
{

	int i, nx, n_layer;
	char nom_indice[10];
	char *erreur="NO_ERROR                     ";
	FILE *fp;
	double *profil_tmp, **profil = par->profil;
	int N_profil = par->N_profil;
	int N_layers = par->N_layers;
	double h = par->h;
	complex indice;
	
	profil_tmp = (double *) malloc(sizeof(double)*N_profil*(N_layers+1));

	if (!(fp = fopen(nom_fichier,"r"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nom_fichier);
		exit(EXIT_FAILURE);
	}
	/* Lecture des indices des couches et calculs des k2 et 1/k2 */
	for (i=1; i<=N_layers; i++){
		sprintf(nom_indice,"n%d",i);
		if (lire_complex(fp, nom_indice, &indice)) erreur=nom_indice;
		par->k2_layer[i]    = (indice*par->k0)*(indice*par->k0); 
		par->invk2_layer[i] = 1/par->k2_layer[i]; 
	}
	par->k2_layer[0]             = (par->n_super*par->k0)*(par->n_super*par->k0);
	par->k2_layer[N_layers+1]    = (par->n_sub*par->k0)*(par->n_sub*par->k0);
	par->invk2_layer[0]          = 1/par->k2_layer[0];
	par->invk2_layer[N_layers+1] = 1/par->k2_layer[N_layers+1];

	fclose(fp);	
	/* Vérification de l'absence d'erreurs de lecture */
	if (strcmp(erreur,"NO_ERROR                     ")){
		fprintf(stderr, "%s : Erreur, probleme de lecture de \"%s\"\n",__FILE__,erreur);
		exit(EXIT_FAILURE);
	}

	
	/* Lecture des profils */
	/* On lit comme un seul tableau, en lisant les lignes les unes à la suite des autres,  */
	/* en considérant qu'une colonne représente les coordonnées d'une interface. On sépare */
	/* ensuite les données en autant de tableaux qu'il y a d'interfaces                    */
	if (par->verbose) fprintf(stdout,"Lecture du profil %s : ",nom_fichier);
	if (lire_tab(nom_fichier, "profil", profil_tmp, N_profil*(N_layers+1)) == 0) {
		if (par->verbose) fprintf(stdout,"OK\n");
	}else{
		fprintf(stderr,"ERREUR de lecture du profil\n");
		exit(EXIT_FAILURE);
	}

	/* Vérification que le profil varie dans l'intervale [0 1] */
	double max = profil_tmp[0];
	double min = profil_tmp[0];
	for (i=1; i<=N_profil*(N_layers+1)-1; i++) {
		max = MAX(max, profil_tmp[i]);
		min = MIN(min, profil_tmp[i]);
	}
	double eps = 1.0e-10; 
	/* Si profil compris dans [0 1], avertissement seulement */
	if ((max<1-eps && min>=-eps) || (max<=1+eps && min>eps)) {
		fprintf(stderr,"ATTENTION, le profil %s varie dans l'intervale [%f %f] et non [0 1] ! \n",nom_fichier,min,max);
	}
	/* Si profil non compris dans [0 1], avertissement et renormalisation de 0 à 1 */
	if (max > 1+eps || min < -eps) {
		fprintf(stderr,"ATTENTION, le profil %s varie dans l'intervale [%f %f] et non [0 1] ! \n",nom_fichier,min,max);
		fprintf(stderr,"=> Normalisation entre 0 et 1\n");
		for (i=0; i<=N_profil-1; i++) {
			profil_tmp[i] = (profil_tmp[i]-min)/(max-min);
		}
	}

	/* Vérification que les profils ne se chevauchent pas */


	/* Normalisation entre 0 et h */
	for (i=0; i<=N_profil*(N_layers+1)-1; i++) {
		profil_tmp[i] *= h;
	}
	if (par->verbose) fprintf(stdout,"Normalisation du profil entre 0 et h : OK\n");

	/* Réarrangement en plusieurs tableaux */
	double eps2 = h*1e-10;
	for (nx=0; nx<=N_profil-1;nx++){
		/* "haut du superstrat", z=0 */
		profil[0][nx] = 0-eps2;
		/* Couches */
		for (n_layer=1; n_layer<=N_layers+1; n_layer++){
			profil[n_layer][nx] = profil_tmp[nx*(N_layers+1)+n_layer-1];
		}
		/* "Bas du substrat", z=h */
		profil[N_layers+2][nx] = h+eps2;
	}
			
	free(profil_tmp);

	/* On inverse tout ! (convention MAP2) */

	return 0;
}


int md1D_ecrire_results(char *nom_fichier, struct Param_struct *par, struct Efficacites_struct *eff)
{
	FILE *fp;
	int LMAX = 1000;



	/* Ouverture du fichier */
	if (!(fp = fopen(nom_fichier,"w"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nom_fichier);
		exit(EXIT_FAILURE);
	}

	fprintf(fp,"#    \n");
	fprintf(fp,"#    \n\n");

	fprintf(fp,    "#------------ Paramètres ------------\n");
	fprintf(fp,"Pola    : %s\n",(par->pola==TE ? "TE" : "TM"));	
	fprintf(fp,"n_super = %f + i%f\n", creal((par->n_super)), cimag((par->n_super)));
	fprintf(fp,"n_sub   = %f + i%f\n", creal((par->n_sub)), cimag((par->n_sub)));
	fprintf(fp,"L       = %f\n",par->L);
	fprintf(fp,"h       = %f\n",par->h);
	fprintf(fp,"lambda  = %f\n",par->lambda);
	fprintf(fp,"angle_i = %f\n",par->angle_i);
	fprintf(fp,"N       = %d\n",par->N);
	fprintf(fp,"NS      = %d\n",par->NS);
	fprintf(fp,"delta_h = %f\n",par->delta_h);
	fprintf(fp,"Nstep   = %d\n",par->Nstep);
/*	fprintf(fp,"fichier_profil = %s\n",fichier_profil);*/
	fprintf(fp,"N_profil = %d\n",par->N_profil);

/*	fprintf(fp,  "\n#--------------- Champs ---------------\n");
	fprintf(fp,  "Ai = "); ecrire_dble_tab(fp, par->Ai, 2*par->N+1, " ", LMAX,"\n");
	fprintf(fp,"\nA0 = "); ecrire_dble_tab(fp, par->A0, 2*par->N+1, " ", LMAX,"\n");
	fprintf(fp,"\nAh = "); ecrire_dble_tab(fp, par->Ah, 2*par->N+1, " ", LMAX,"\n");
*/

	fprintf(fp,  "\n#------------ Efficacités -------------\n");
	fprintf(fp,"somm_eff          = % 1.6le\n",eff->somm_eff);
	fprintf(fp,"un_moins_somm_eff = % 1.6e\n",1.0-eff->somm_eff);
	fprintf(fp,"somm_eff_R        = % 1.6e\n",eff->somm_eff_R);
	fprintf(fp,"somm_eff_T        = % 1.6e\n",eff->somm_eff_T);
	fprintf(fp,  "eff_R       = "); ecrire_dble_tab(fp, eff->eff_R, eff->Nmax_super-eff->Nmin_super+1, " ", LMAX,"\n");
	fprintf(fp,"\nN_eff_R     = "); ecrire_dble_tab(fp, eff->N_eff_R, eff->Nmax_super-eff->Nmin_super+1, " ", LMAX,"\n");
	fprintf(fp,"\ntheta_eff_R = "); ecrire_dble_tab(fp, eff->theta_eff_R, eff->Nmax_super-eff->Nmin_super+1, " ", LMAX,"\n");
	fprintf(fp,"\neff_T       = "); ecrire_dble_tab(fp, eff->eff_T, eff->Nmax_sub-eff->Nmin_sub+1, " ", LMAX,"\n");
	fprintf(fp,"\nN_eff_T     = "); ecrire_dble_tab(fp, eff->N_eff_T, eff->Nmax_sub-eff->Nmin_sub+1, " ", LMAX,"\n");
	fprintf(fp,"\ntheta_eff_T = "); ecrire_dble_tab(fp, eff->theta_eff_T, eff->Nmax_sub-eff->Nmin_sub+1, " ", LMAX,"\n");

	fprintf(fp,"\n\n#--------------- Profil ---------------\n");
	fprintf(fp,"profil = "); ecrire_dble_tab(fp, par->profil[0], par->N_profil, " ", LMAX,"\n");

	fprintf(fp,    "\n\n");

	fclose(fp);
	return 0;
}

/*!	\fn		int md1D_genere_nom_fichier_results(char *nomfichier_results, struct Param_struct *par)
 *
 *	\brief
 */
int md1D_genere_nom_fichier_results(char *nomfichier_results, struct Param_struct *par)
{
	time_t ptime;
	time(&ptime);
	struct tm  temps;
	localtime_r(&ptime, &temps);
	
	/* Génération d'un nom de la forme  nom_2004_10_12_16h34.mdi */
/*	sprintf(nomfichier_results,"results/%s_%d_%02d_%02d_%02d%s%02d.txt",\
		par->nom_profil, temps.tm_year+1900, temps.tm_mon+1, temps.tm_mday, temps.tm_hour, "h", temps.tm_min);
*/
	sprintf(nomfichier_results,"%s.txt",par->nom_profil);

	return 1;

}

/*!	\fn		int	md1D_ecrire_config(struct Param_struct *par, struct Noms_fichiers *nomfichier)
 *
 *	\brief
 */
int	md1D_ecrire_config(struct Param_struct *par, struct Noms_fichiers *nomfichier)
{
	FILE *fp;

	if (!(fp = fopen(nomfichier->fichier_config,"w"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nomfichier->fichier_config);
		return 1;
	}

	fprintf(fp,"#\n#   Fichier de configuration pour le programme md1D\n#\n\n");
	fprintf(fp,"# Fichier contenant les paramètres du système :\n");
	fprintf(fp,"fichier_param = %s\n", nomfichier->fichier_param);
	fprintf(fp,"# Fichier contenant les résultats du dernier calcul :\n");
	fprintf(fp,"fichier_results = %s\n", nomfichier->fichier_results);


	fclose(fp);
	return 0;
}

/*! \fn    int lire_int(FILE *fp, char *label, int *value)
 *
 *  \brief	Lit dans le fichier pointé par *fp la valeur entiere 'value' indiquée par 'label' \n
 *			sous la forme label = value (ex.: N2 = 10)
 *  \return	0 si lecture réussie 1 sinon
 */
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

/*! \fn    int lire_double(FILE *fp, char *label, double *value)
 *
 *  \brief	Lit dans le fichier pointé par *fp la valeur entiere 'value' indiquée par 'label' \n
 *			sous la forme label = value (ex.:  x = 0.12310)
 *  \return	0 si lecture réussie 1 sinon
 */
int lire_double(FILE *fp, const char *label, double *value){

	char *pos, *pos2, stmp[SIZE_LINE_BUFFER];
	float tmp;

	rewind(fp);
	while(!feof(fp)){
		if (lire_ligne(fp,stmp) != 0) break;
		skip_comment(stmp); /* Vire les commentaires */
		if ((pos = label_search(stmp,label)) != NULL){ /* Recherche le label */	
			if((pos2 = strchr(pos+strlen(label),'=')) != NULL) *pos2 = ' '; /* remplace '=' par un espace */
			if (sscanf(pos+strlen(label)," %f", &tmp) == 1){ /* Lit la valeur */
				*value = (double) tmp;
				return 0;
			}
		}
	}
	return 1;
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

/*! \fn    int lire_complex(FILE *fp,  char *label, complex *value)
 *
 *  \brief	Lit dans le fichier pointé par *fp la valeur entiere 'value' indiquée par 'label' \n
 *			sous la forme label = value (ex.: Z1 = 1.0 + i0.5 )
 *  \return	0 si lecture réussie 1 sinon
 */
int lire_complex(FILE *fp, const char *label, complex *value){

	char *pos, *pos2, stmp[SIZE_LINE_BUFFER];
	float tmp1, tmp2;

	rewind(fp);
	while(!feof(fp)){
		if (lire_ligne(fp,stmp) != 0) break;
		skip_comment(stmp); /* Vire les commentaires */
		if ((pos = label_search(stmp,label)) != NULL){ /* Recherche le label */	
			if((pos2 = strchr(pos+strlen(label),'=')) != NULL) *pos2 = ' '; /* remplace '=' par un espace */
			if (sscanf(pos+strlen(label)," %f + i%f", &tmp1, &tmp2) == 2){ /* Lit les valeurs */
				*value = c_omplex(tmp1, tmp2);
				return 0;
			}
		}
	}
	return 1;
}


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


/*!	\fn		int	ecrire_dble_tab(FILE *fp, double *tab, int N, char *separateur1, int Nmax1, char *separateur2)
 *
 *	\brief	Ecrit les valeurs d'un tableau séparées par les 'séparateurs1' (par ex " "), plus par les \n
 *          'séparateurs2' (par ex "\n") une fois tous les Nmax1 éléments.
 */
int	ecrire_dble_tab(FILE *fp, double *tab, int N, char *separateur1, int Nmax1, char *separateur2)
{
	int i,k=0;

	while((k+1)*Nmax1 < N) {
		for (i=k*Nmax1; i<=MIN(N-1,(k+1)*Nmax1-1); i++){
			fprintf(fp,"% 1.6e%s",tab[i],separateur1);
		}
		k++;
		fprintf(fp,"%s",separateur2);
	}
	/* Derniere ligne, traitée à part car pas de séparateur2 à la fin ! */
	for (i=k*Nmax1; i<=MIN(N-1,(k+1)*Nmax1-1); i++){
		fprintf(fp,"% 1.6e%s",tab[i],separateur1);
	}
	
	return 0;
}

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
	char nom_fichier_tmp[] = "md1D_fichier_tmp_68gIg78GUgkd.tmp";
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

