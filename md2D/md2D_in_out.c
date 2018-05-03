/*! \file		md2D_in_out.c
 *
 *	\brief		Routines de gestion des entrées-sorties pour md2D
 *
 *	\version	0.1
 *  \date		../../2004
 *  \authors	Laurent ARNAUD
 */
#include "md2D_in_out.h"

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_lire_param(struct Param_struct *par){
 *	
 *	\brief	Lecture des paramètres dans un fichier
 *
 *---------------------------------------------------------------------------------------------*/
int md2D_lire_param(struct Noms_fichiers *nomfichier, struct Param_struct *par){

	FILE *fp;
	int ret=0, argc=par->argc;
	double n_super_re, n_super_im, n_sub_re, n_sub_im, nu_layer_re, nu_layer_im, tmp_dble;
	char *erreur="NO_ERROR                     ";
	char **argvcp=par->argvcp, str_profil[SIZE_STR_BUFFER], str_tmp[SIZE_STR_BUFFER];

	/* Nom de fichier param : en ligne de commande ou par défaut */
	if (lire_str_arg(nomfichier->fichier_param, "-param", argc, argvcp)) {
		sprintf(nomfichier->fichier_param,"md2D_param.txt");}

	/* Ouverture de fichioer_param */
	if (!(fp = fopen(nomfichier->fichier_param,"r"))){
		fprintf(stderr, "%s line %d: ERROR, can't open %s\n",__FILE__, __LINE__,nomfichier->fichier_param);
		exit(EXIT_FAILURE);
	}
	
	/* Verbosity level reading */
	if (lire_int_arg(&par->verbosity, "-verbosity", argc, argvcp)) {
		if (lire_int (fp, "verbosity", &(par->verbosity) )) par->verbosity = 2;}
	if (par->verbosity >=2 ) fprintf(stdout,"Reading parameters in %s\n",nomfichier->fichier_param);
	
	/* Lecture des paramètres, d'abord en ligne de commande, si rien en ligne de commande,
	   lecture dans fichier_param, sinon erreur et arret du programme */
	arg_read(argc, argvcp, fp, "string", (void *) nomfichier->profile_file, "profile_file", EXIT_ON_ERROR);
	arg_read(argc, argvcp, fp, "string", (void *) par->calcul_type, "calcul_type", EXIT_ON_ERROR);
	arg_read(argc, argvcp, fp, "string", (void *) par->calcul_method, "calcul_method", EXIT_ON_ERROR);
	arg_read(argc, argvcp, fp, "string", (void *) par->profile_name, "profile_name", EXIT_ON_ERROR);
	arg_read(argc, argvcp, fp, "string", (void *) par->i_field_mode, "i_field_mode", EXIT_ON_ERROR);
	arg_read(argc, argvcp, fp, "string", (void *) str_tmp, "pola", EXIT_ON_ERROR);
	if      (!strcmp(str_tmp,"TE")){
		par->pola = TE;
	}else if (!strcmp(str_tmp,"TM")){
		par->pola = TM;
	}else{
		fprintf(stderr, "%s line %d: ERROR, pola must be either 'TE' or 'TM' (instead of %s)\n",__FILE__, __LINE__,str_tmp);
		exit(EXIT_FAILURE);
	}
	arg_read(argc, argvcp, fp, "double", (void *) &par->L, "L", EXIT_ON_ERROR);
	arg_read(argc, argvcp, fp, "double", (void *) &par->lambda, "lambda", EXIT_ON_ERROR);
	arg_read(argc, argvcp, fp, "double", (void *) &par->theta_i, "theta_i", EXIT_ON_ERROR);
	par->theta_i *= PI/180.0;
	if (arg_read(argc, argvcp, fp, "complex", (void *) &par->sigma0_normed, "sigma0_normed", CONTINUE_ON_ERROR) != 0) par->sigma0_normed = SIGMA0_NORMED_NOT_DEFINED;
/*printf("sigma0_normed = %f + i%f\n",creal(par->sigma0_normed),cimag(par->sigma0_normed));*/


/*	arg_read(argc, argvcp, fp, "complex", (void *) &par->n_super, "n_super", EXIT_ON_ERROR);
*/
	if (!strcmp(par->calcul_type,"SWIFTS") || !strcmp(par->calcul_type,"I_SWIFTS")){
		arg_read(argc, argvcp, fp, "double", (void *) &par->h_layer, "h_layer", EXIT_ON_ERROR);
		if (lire_dble_arg(&nu_layer_re, "-nu_layer_re", argc, argvcp)) {
			if (lire_complex(fp, "nu_layer", &(par->nu_layer))) erreur="nu_layer";
		}else{
			if (lire_dble_arg(&nu_layer_im, "-nu_layer_im", argc, argvcp)) {
				nu_layer_im =0;}
			par->nu_layer = c_omplex(nu_layer_re, nu_layer_im);}
	}

	if (lire_dble_arg(&n_super_re, "-n_super_re", argc, argvcp)) {
		if (lire_complex(fp, "n_super", &(par->n_super))) erreur="n_super";
	}else{
		if (lire_dble_arg(&n_super_im, "-n_super_im", argc, argvcp)) {
			n_super_im =0;}
		par->n_super = c_omplex(n_super_re, n_super_im);}

	if (lire_dble_arg(&n_sub_re, "-n_sub_re", argc, argvcp)) {
		if (lire_complex(fp, "n_sub", &(par->n_sub))) erreur="n_sub";
	}else{
		if (lire_dble_arg(&n_sub_im, "-n_sub_im", argc, argvcp)) {
			n_sub_im =0;}
		par->n_sub = c_omplex(n_sub_re, n_sub_im);}

	/* imposed_Ssteps */
	if(arg_read(argc, argvcp, fp, "string", (void *) par->imposed_S_steps, "imposed_S_steps", CONTINUE_ON_ERROR) != 0){
		strcpy(par->imposed_S_steps,"NONE");
	}else{
		if(arg_read(argc, argvcp, fp, "double", (void *) &par->hmin_Sstep, "hmin_Sstep", CONTINUE_ON_ERROR) != 0) par->hmin_Sstep = SUPER_BIG_HMIN_SSTEP;
printf("hmin_Sstep = %f\n",par->hmin_Sstep);
		if (!strcmp(par->imposed_S_steps,"FROM_FILE")){
			arg_read(argc, argvcp, fp, "string", (void *) par->imposed_S_steps_filename, "imposed_S_steps_filename", EXIT_ON_ERROR);
			if (lire_tab(par->imposed_S_steps_filename, "N_imposed_S_steps", &tmp_dble, 1) != 0) {
				fprintf(stderr, "%s line %d: ERROR, can't read N_imposed_S_steps in %s\n",__FILE__, __LINE__,par->imposed_S_steps_filename);
				exit(EXIT_FAILURE);
			}else{
				par->N_imposed_S_steps = ROUND(tmp_dble);
			}
		}
	}

/*	if (lire_int_arg(&par->imposed_S_steps, "-imposed_S_steps", argc, argvcp)) {
		if (lire_int (fp, "imposed_S_steps", &(par->imposed_S_steps) )) par->imposed_S_steps = 0;}
	if (par->imposed_S_steps){
		if (lire_str_arg(par->imposed_S_steps_filename, "-imposed_S_steps_filename", argc, argvcp)) {
			if (lire_string (fp, "imposed_S_steps_filename", par->imposed_S_steps_filename)) erreur="imposed_S_steps_filename";}
		if (lire_int_arg(&par->N_imposed_S_steps, "-N_imposed_S_steps", argc, argvcp)) {
			if (lire_int (fp, "N_imposed_S_steps", &(par->N_imposed_S_steps) )) erreur="N_imposed_S_steps";}
	}*/

/********************/
/* NE PAS EFFACER pour l'instant(peut servir) */
/********************/
/*	if (lire_int_arg(&par->mode_extract_S, "-mode_extract_S", argc, argvcp)) {
		if (lire_int (fp, "mode_extract_S", &(par->mode_extract_S) )) par->mode_extract_S = 0;}
	if (lire_int_arg(&par->READ_MAT_S, "-READ_MAT_S", argc, argvcp)) {
		if (lire_int (fp, "READ_MAT_S", &(par->READ_MAT_S) )) par->READ_MAT_S = 0;}
	if (par->mode_extract_S){ 	
		if (lire_dble_arg(&par->h_extract_S, "-h_extract_S", argc, argvcp)) {
			if (lire_double (fp, "h_extract_S", &(par->h_extract_S) )) erreur="h_extract_S";}}	
	if (par->READ_MAT_S){ 	
		if (lire_str_arg(par->mat_S_file, "-mat_S_file", argc, argvcp)) {
			if (lire_string (fp, "mat_S_file", par->mat_S_file)) erreur="mat_S_file";}}*/


/* TabNS */
/*	if (lire_int_arg(&par->READ_tab_NS, "-READ_tab_NS", argc, argvcp)) {
		if (lire_int (fp, "READ_tab_NS", &(par->READ_tab_NS) )) par->READ_tab_NS = 0;}
	if (par->READ_tab_NS){
		if (lire_str_arg(par->tab_NS_filename, "-tab_NS_filename", argc, argvcp)) {
			if (lire_string (fp, "tab_NS_filename", par->tab_NS_filename)) erreur="tab_NS_filename";}
	}*/


	
/*	if (lire_dble_arg(&par->phi_i, "-phi_i", argc, argvcp)) {
		if (lire_double (fp, "phi_i", &(par->phi_i) )) erreur="phi_i";}
		par->phi_i *= PI/180.0;
	if (lire_dble_arg(&par->psi, "-psi", argc, argvcp)) {
		if (lire_double (fp, "psi", &(par->psi) )) erreur="psi";}
		par->psi *= PI/180.0;*/


	/* h */
	if (lire_dble_arg(&par->h, "-h", argc, argvcp)) {
		if (lire_str_arg(str_tmp, "-h", argc, argvcp)) {
			if (lire_double (fp, "h", &(par->h) )) {
				if (lire_string (fp, "h", str_tmp)) {
					erreur = "h";
				}else if (!strcmp(str_tmp,"AUTO")) {
					par->h = AUTO;
				}else {
					erreur = "h";
				}
			}	
		}else if (!strcmp(str_tmp,"AUTO")) {
			par->h = AUTO;
		}else{
			erreur = "h";
		}
	}

	/* delta_h */
	if (lire_dble_arg(&par->delta_h, "-delta_h", argc, argvcp)) {
		if (lire_str_arg(str_tmp, "-delta_h", argc, argvcp)) {
			if (lire_double (fp, "delta_h", &(par->delta_h) )) {
				if (lire_string (fp, "delta_h", str_tmp)) {
					erreur = "delta_h";
				}else if (!strcmp(str_tmp,"AUTO")) {
					par->delta_h = AUTO;
				}else {
					erreur = "delta_h";
				}
			}	
		}else if (!strcmp(str_tmp,"AUTO")) {
			par->delta_h = AUTO;
		}else{
			erreur = "delta_h";
		}
	}	
	/* N */
	if (lire_int_arg(&par->N, "-N", argc, argvcp)) {
		if (lire_str_arg(str_tmp, "-N", argc, argvcp)) {
			if (lire_int (fp, "N", &(par->N) )) {
				if (lire_string (fp, "N", str_tmp)) {
					erreur = "N";
				}else if (!strcmp(str_tmp,"AUTO")) {
					par->N = AUTO;
				}else {
					erreur = "N";
				}
			}	
		}else if (!strcmp(str_tmp,"AUTO")) {
			par->N = AUTO;
		}else{
			erreur = "N";
		}
	}	
	/* NS */
	if (lire_int_arg(&par->NS, "-NS", argc, argvcp)) {
		if (lire_str_arg(str_tmp, "-NS", argc, argvcp)) {
			if (lire_int (fp, "NS", &(par->NS) )) {
				if (lire_string (fp, "NS", str_tmp)) {
					erreur = "NS";
				}else if (!strcmp(str_tmp,"AUTO")) {
					par->NS = AUTO;
				}else {
					erreur = "NS";
				}
			}	
		}else if (!strcmp(str_tmp,"AUTO")) {
			par->NS = AUTO;
		}else{
			erreur = "NS";
		}
	}	

	fclose(fp);

	/* Vérification de l'absence d'erreurs de lecture */
	if (strcmp(erreur,"NO_ERROR                     ")){
		fprintf(stderr, "%s: ERROR, can't read \"%s\"\n",__FILE__,erreur);
		exit(EXIT_FAILURE);
	}

	/* Lecture des paramètres dans fichier profil */
	/* (Les paramètres liés à la nature ou indisociables du profil sont contenus dans profile_file) */
	if (par->verbosity >= 2) fprintf(stdout,"Reading parameters in %s\n",nomfichier->profile_file);
	if (!(fp = fopen(nomfichier->profile_file,"r"))){
		fprintf(stderr, "%s line %d: ERROR, can't open %s\n",__FILE__, __LINE__,nomfichier->profile_file);
		exit(EXIT_FAILURE);
	}
	erreur="NO_ERROR                     ";
	if (lire_string (fp, "type_profil", str_profil)) {
		erreur="type_profil";
	}else{
		if       (!strcmp(str_profil,"H_X"))         {par->type_profil = H_X;
			if (lire_int (fp, "N_x", &(par->N_x) )) erreur="N_x";
			par->N_layers = 0;
		}else if (!strcmp(str_profil,"MULTICOUCHES")){par->type_profil = MULTICOUCHES;
			if (lire_int (fp, "N_x", &(par->N_x) )) erreur="N_x";
			if (lire_int (fp, "N_layers", &(par->N_layers) )) erreur="N_layers";
		}else if (!strcmp(str_profil,"N_XYZ"))       {par->type_profil = N_XYZ;
			if (lire_int (fp, "N_x", &(par->N_x) )) erreur="N_x";
			if (lire_int (fp, "N_z", &(par->N_z) )) erreur="N_z";
		}else                                        {erreur = "type_profil_bis";}
	}

	fclose(fp);
	
	/* Vérification de l'absence d'erreurs de lecture */
	if (strcmp(erreur,"NO_ERROR                     ")){
		fprintf(stderr, "%s: ERROR, can't read \"%s\"\n",__FILE__,erreur);
		exit(EXIT_FAILURE);
	}

	return ret;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_affiche_valeurs_param(struct Param_struct *par, struct Noms_fichiers *nomfichier)
 *	
 *	\brief	Show parameters values on the screen
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_affiche_valeurs_param(struct Param_struct *par, struct Noms_fichiers *nomfichier)
{

	/* Affichage des valeurs lues */
	if (par->verbosity >= 2){
		fprintf(stdout,"n_super = %f + i%f\n", creal((par->n_super)), cimag((par->n_super)));
		fprintf(stdout,"n_sub   = %f + i%f\n", creal((par->n_sub)), cimag((par->n_sub)));
		fprintf(stdout,"lambda  = %f\n",par->lambda);
		fprintf(stdout,"theta_i = %f rad (%f deg)\n",par->theta_i,par->theta_i*180.0/PI);
		fprintf(stdout,"L       = %f\n",par->L);
		fprintf(stdout,"h       = %f\n",par->h);
		fprintf(stdout,"delta_h = %f\n",par->delta_h);
		fprintf(stdout,"N       = %d\n",par->N);
		fprintf(stdout,"NS      = %d\n",par->NS);
		fprintf(stdout,"N_x     = %d\n",par->N_x);
		fprintf(stdout,"pola    : %s\n",(par->pola==TE ? "TE" : "TM"));
		fprintf(stdout,"calcul_type    = %s\n",par->calcul_type);
		fprintf(stdout,"calcul_method  = %s\n",par->calcul_method);
		fprintf(stdout,"i_field_mode   = %s\n",par->i_field_mode);
		fprintf(stdout,"profile_name     = %s\n",par->profile_name);
		fprintf(stdout,"profile_file = %s\n",nomfichier->profile_file);
		fprintf(stdout,"type_profil    = %s\n",(par->type_profil==H_X ? "H_X" :
		                                        (par->type_profil==N_XYZ ? "N_XYZ":"MULTICOUCHES")));
		fprintf(stdout,"N_couches      = %d\n",par->N_layers);
		fflush(stdout);
	}
	return 0;
}



/*---------------------------------------------------------------------------------------------*/
/*!	\fn	int md2D_lire_profil_H_X(const char *nom_fichier, struct Param_struct *par)
 *
 *	\brief	Reads the h(x) profile, determine the height h, gives an alert in case it is
 * 			different from the indicated value, the indicated values is prioritary.
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_lire_profil_H_X(const char *nom_fichier, struct Param_struct *par)
{

	int i;
	int N_x = par->N_x;
	double h_tmp, *profil = par->profil[0];

	/* Entry of k2 & invk2 for substrate and superstrate in arrays (inv)k2_layer, for compatibility with multicouches */
	par->k2_layer[0] = (par->k_super)*(par->k_super);
	par->k2_layer[par->N_layers+1] = (par->k_sub)*(par->k_sub);
	par->invk2_layer[0] = 1/((par->k_super)*(par->k_super));
	par->invk2_layer[par->N_layers+1] = 1/((par->k_sub)*(par->k_sub));

	/* Reading the profile */
	if (par->verbosity >= 2) fprintf(stdout,"Reading the profile %s : ",nom_fichier);
	if (lire_tab(nom_fichier, "profil", profil, N_x) == 0) {
		if (par->verbosity >= 2) fprintf(stdout,"OK\n");
	}else{
		fprintf(stderr,"Profile reading ERROR\n");
		exit(EXIT_FAILURE);
	}

	/* Determination of the profile height from read h(x) points */
	double h_min = profil[0];
	double h_max = profil[0];
	double eps = 1.0e-10; 
	for (i=0;i<=par->N_x-1;i++){
		h_min = MIN(h_min,profil[i]);
		h_max = MAX(h_max,profil[i]);
	}
/*	if par->SPECIAL_H_MIN_H_MAX { *//* For some special purposes, eg near field map, it is usefull to set h_min & h_max outside of the limit of the modulated zone */
	lire_tab(nom_fichier, "h_min", &h_min, 1);
	lire_tab(nom_fichier, "h_max", &h_max, 1);
	if (par->verbosity > 0) {
		fprintf(stdout,"SPECIAL_H_MIN_H_MAX, h_min=%f, h_max=%f\n",h_min, h_max);
	}
/*	}*/
	h_tmp = h_max - h_min;

	/* Checking that h(determined) = h(indicated) and setting par->h to determined value if h = AUTO */
	if (par->h == AUTO) {
		par->h = h_tmp;
	}else if (fabs(par->h - h_tmp) > eps*par->h && par->verbosity >= 2) {
		fprintf(stdout,  "+---------------------------------------------------------------------\
				\n|                        CAUTION !\
				\n| h determined for %s is %f and h indicated %f !\
				\n+---------------------------------------------------------------------\
				\n",nom_fichier,h_tmp, par->h);
	}
	/* Normalising beetween 0 and h (instead of [h0,h0+h]*/
	for (i=0; i<=N_x-1; i++) {
		profil[i] -= h_min;
	}

	/* Multiplication by h_wanted / h_determined */
	if (h_tmp != 0.0) {
		for (i=0; i<=N_x-1; i++) {
			profil[i] *= par->h / h_tmp;
		}
	}

	if (par->verbosity >= 2) fprintf(stdout,"Profile normalisation between 0 and h : OK\n");


	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_lire_profil_MULTI(const char *nom_fichier, struct Param_struct *par)

 *
 *	\brief	Fonction lisant les valeurs décrivant un profil multi
 *
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_lire_profil_MULTI(const char *nom_fichier, struct Param_struct *par)
{

	int i, nx, n_layer;
	char nom_indice[SIZE_STR_BUFFER];
	char *erreur="NO_ERROR                     ";
	FILE *fp;
	double *profil_tmp, **profil = par->profil;
	int N_x = par->N_x;
	int N_layers = par->N_layers;
	double h_tmp;
	COMPLEX indice;
	
	profil_tmp = (double *) malloc(sizeof(double)*N_x*(N_layers+1));

	if (!(fp = fopen(nom_fichier,"r"))){
		fprintf(stderr, "%s line %d: ERROR, can't open %s\n",__FILE__, __LINE__,nom_fichier);
		exit(EXIT_FAILURE);
	}
	/* Lecture des indices des couches et calculs des k2 et 1/k2 */
	for (i=1; i<=N_layers; i++){
		sprintf(nom_indice,"n%d",i);
		if (lire_complex(fp, nom_indice, &indice)) erreur=nom_indice;
		par->k2_layer[i]    = (indice*par->k_super/par->n_super)*(indice*par->k_super/par->n_super); 
		par->invk2_layer[i] = 1/par->k2_layer[i]; 
	}
	par->k2_layer[0]             = (par->k_super)*(par->k_super);
	par->k2_layer[N_layers+1]    = (par->k_sub)*(par->k_sub);
	par->invk2_layer[0]          = 1/par->k2_layer[0];
	par->invk2_layer[N_layers+1] = 1/par->k2_layer[N_layers+1];

	fclose(fp);	
	/* Vérification de l'absence d'erreurs de lecture */
	if (strcmp(erreur,"NO_ERROR                     ")){
		fprintf(stderr, "%s: ERROR, can't read \"%s\"\n",__FILE__,erreur);
		exit(EXIT_FAILURE);
	}

	
	/* Lecture des profils */
	/* On lit comme un seul tableau, en lisant les lignes les unes à la suite des autres,  */
	/* en considérant qu'une colonne représente les coordonnées d'une interface. On sépare */
	/* ensuite les données en autant de tableaux qu'il y a d'interfaces                    */
	if (par->verbosity >= 2) fprintf(stdout,"Reading profile %s: ",nom_fichier);
	if (lire_tab(nom_fichier, "profil", profil_tmp, N_x*(N_layers+1)) == 0) {
		if (par->verbosity >= 2) fprintf(stdout,"OK\n");
	}else{
		fprintf(stderr,"ERROR of profile reading\n");
		exit(EXIT_FAILURE);
	}
/**************************************************/
	/* Détermination de la hauteur du profil à partir des points h(x) */
	double h_min = profil_tmp[0];
	double h_max = profil_tmp[0];
	double eps = 1.0e-10; 
	for (i=0;i<=par->N_x*(N_layers+1)-1;i++){
		h_min = MIN(h_min,profil_tmp[i]);
		h_max = MAX(h_max,profil_tmp[i]);
	}
/*	if par->SPECIAL_H_MIN_H_MAX { *//* For some special purposes, eg near field map, it is usefull to set h_min & h_max outside of the limit of the modulated zone */
	lire_tab(nom_fichier, "h_min", &h_min, 1);
	lire_tab(nom_fichier, "h_max", &h_max, 1);

	h_tmp = h_max - h_min;

	/* Vérification que h(deteminé) = h(indiqué) et attribution si h = AUTO */
	if (par->h == AUTO) {
		par->h = h_tmp;
	}else if (fabs(par->h - h_tmp) > eps*par->h) {
		fprintf(stderr,  "+---------------------------------------------------------------------\
				\n|                        WARNING !\
				\n| h determined for %s is %f while h indicated is %f !\
				\n+---------------------------------------------------------------------\
				\n",nom_fichier,h_tmp, par->h);
	}
	/* Normalisation entre 0 et h (au lieu de [h0,h0+h]*/
	for (i=0; i<=N_x*(N_layers+1)-1; i++) {
		profil_tmp[i] -= h_min;
	}
	/* Vérification que les profils ne se chevauchent pas */

			
	/* Multiplication par h_voulu / h_mesuré */
	if (h_tmp != 0.0) {
		for (i=0; i<=N_x*(N_layers+1)-1; i++) {
			profil_tmp[i] *= par->h / h_tmp;
		}
	}
	
	if (par->verbosity >= 2) fprintf(stdout,"Normalisation of profile between 0 and h: OK\n");

/**************************************************/

/*	double max = profil_tmp[0];
	double min = profil_tmp[0];
	for (i=1; i<=N_x*(N_layers+1)-1; i++) {
		max = MAX(max, profil_tmp[i]);
		min = MIN(min, profil_tmp[i]);
	}
	double eps = 1.0e-10; 
	if ((max<1-eps && min>=-eps) || (max<=1+eps && min>eps)) {
		fprintf(stderr,"ATTENTION, le profil %s varie dans l'intervale [%f %f] et non [0 1] ! \n",nom_fichier,min,max);
	}
	if (max > 1+eps || min < -eps) {
		fprintf(stderr,"ATTENTION, le profil %s varie dans l'intervale [%f %f] et non [0 1] ! \n",nom_fichier,min,max);
		fprintf(stderr,"=> Normalisation entre 0 et 1\n");
		for (i=0; i<=N_x*(N_layers+1)-1; i++) {
			profil_tmp[i] = (profil_tmp[i]-min)/(max-min);
		}
	}
	for (i=0; i<=N_x*(N_layers+1)-1; i++) {
		profil_tmp[i] *= h;
	}
	if (par->verbosity) fprintf(stdout,"Normalisation du profil entre 0 et h : OK\n");
*/
/**************************************************/
	
	/* Réarrangement en plusieurs tableaux */
	double eps2 = 0*par->h*1e-10;
	for (nx=0; nx<=N_x-1;nx++){
		/* superstrat, z=h */
		profil[0][nx] = par->h+eps2;
		/* Couches */
		for (n_layer=1; n_layer<=N_layers+1; n_layer++){
			profil[n_layer][nx] = profil_tmp[nx*(N_layers+1)+n_layer-1];
		}
		/* substrat, z=0 */
		profil[N_layers+2][nx] = 0-eps2;
	}
			
	free(profil_tmp);

	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn	int md2D_lire_profil_N_XYZ(const char *nom_fichier, struct Param_struct *par)
 *
 *	\brief	Fonction lisant les valeurs de indice(x,z) décrivant un profil volumique
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_lire_profil_N_XYZ(const char *nom_fichier, struct Param_struct *par)
{
	int nx, nz, N_x=par->N_x, N_z=par->N_z;
	COMPLEX **n_xyz = par->n_xyz;
	

	/* Allocations de mémoire pour variables temporaires */
	double *Re_n_xyz = malloc(sizeof(double)*N_x*N_z);
	double *Im_n_xyz = malloc(sizeof(double)*N_x*N_z);

	/* Lecture de la partie réelle */
	if (par->verbosity >= 2) fprintf(stdout,"Reading profile %s:\nreal part: ",nom_fichier);
	if (lire_tab(nom_fichier, "Re_n_xyz", Re_n_xyz, N_x*N_z) == 0) {
		if (par->verbosity >= 2) fprintf(stdout,"OK\n");
	}else{
		fprintf(stderr,"Reading profile ERROR\n");
		exit(EXIT_FAILURE);
	}
	/* Lecture de la partie imaginaire */
	if (par->verbosity >= 2) fprintf(stdout,"imaginary part: ");fflush(stdout);
	if (lire_tab(nom_fichier, "Im_n_xyz", Im_n_xyz, N_x*N_z) == 0) {
		if (par->verbosity >= 2) fprintf(stdout,"OK\n");
	}else{
		fprintf(stderr,"NO IMAGINARY PART (dielectric profile)\n");
		for (nx=0; nx<=N_x*N_z-1; nx++){
			Im_n_xyz[nx] = 0;
		}
	}
	
	/* Création de la matrice COMPLEXe n_xyz */
	for (nz=0; nz<=N_z-1; nz++){
		for (nx=0; nx<=N_x-1; nx++){
			n_xyz[nz][nx] = Re_n_xyz[nx+N_x*nz] + I*Im_n_xyz[nx+N_x*nz];
		}
	}

	free(Re_n_xyz);
	free(Im_n_xyz);

	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn	
 *
 *	\brief	Writing results
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_ecrire_results(char *nom_fichier, struct Param_struct *par, struct Efficacites_struct *eff)
{
	FILE *fp;
	int LMAX = 1000;



	/* Ouverture du fichier */
	if (!(fp = fopen(nom_fichier,"w"))){
		fprintf(stderr, "%s line %d: ERROR, can't open %s\n",__FILE__, __LINE__,nom_fichier);
		exit(EXIT_FAILURE);
	}

	fprintf(fp,"#    \n");
	fprintf(fp,"#    \n\n");

	fprintf(fp,    "#------------ Parameters ------------\n");
	fprintf(fp,"calcul_type   : %s\n", par->calcul_type);	
	fprintf(fp,"calcul_method : %s\n", par->calcul_method);	
	fprintf(fp,"pola    : %s\n",(par->pola==TE ? "TE" : "TM"));	
	fprintf(fp,"theta_i = %f deg\n",par->theta_i*180.0/PI);

	fprintf(fp,"n_super = %f + i%f\n", creal((par->n_super)), cimag((par->n_super)));
	fprintf(fp,"n_sub   = %f + i%f\n", creal((par->n_sub)), cimag((par->n_sub)));
	fprintf(fp,"L       = %f\n",par->L);
	fprintf(fp,"h       = %f\n",par->h);
	fprintf(fp,"lambda  = %f\n",par->lambda);
	fprintf(fp,"N       = %d\n",par->N);
	fprintf(fp,"NS      = %d\n",par->NS);
	fprintf(fp,"delta_h = %f\n",par->delta_h);
	fprintf(fp,"N_steps   = %d\n",par->N_steps);
/*	fprintf(fp,"profile_file = %s\n",profile_file);*/
	fprintf(fp,"N_x = %d\n",par->N_x);
	fprintf(fp,"Calcul_duration = %f s \n", md_chrono(par));


	fprintf(fp,  "\n#------------ Results -------------\n");
	fprintf(fp,"sum_eff           = % 1.6le\n",eff->sum_eff);
	fprintf(fp,"one_minus_sum_eff = % 1.6e\n",1.0-eff->sum_eff);
	fprintf(fp,"somm_eff_r        = % 1.6e\n",eff->sum_eff_r);
	fprintf(fp,"somm_eff_t        = % 1.6e\n",eff->sum_eff_t);
	fprintf(fp,"\ntheta_r = "); ecrire_dble_tab(fp, eff->theta_r, eff->Nmax_super-eff->Nmin_super+1, " ", LMAX,"\n");
	fprintf(fp,"\nN_eff_r = "); ecrire_dble_tab(fp, eff->N_eff_r, eff->Nmax_super-eff->Nmin_super+1, " ", LMAX,"\n");
	fprintf(fp,"\neff_r = "); ecrire_dble_tab(fp, eff->eff_r, eff->Nmax_super-eff->Nmin_super+1, " ", LMAX,"\n");
	fprintf(fp,"\ntheta_t = "); ecrire_dble_tab(fp, eff->theta_t, eff->Nmax_sub-eff->Nmin_sub+1, " ", LMAX,"\n");
	fprintf(fp,"\nN_eff_t = "); ecrire_dble_tab(fp, eff->N_eff_t, eff->Nmax_sub-eff->Nmin_sub+1, " ", LMAX,"\n");
	fprintf(fp,"\neff_t = "); ecrire_dble_tab(fp, eff->eff_t, eff->Nmax_sub-eff->Nmin_sub+1, " ", LMAX,"\n");


	fclose(fp);
	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_genere_nom_fichier_results(char *nomfichier_results, struct Param_struct *par)
 *
 *	\brief
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_genere_nom_fichier_results(char *nomfichier_results, struct Param_struct *par)
{
	int localtime_r(time_t *, struct tm *);

	time_t ptime;
	time(&ptime);
	struct tm  temps;
	localtime_r(&ptime, &temps);
	
	/* Génération d'un nom de la forme  nom_2004_10_12_16h34.mdi */
/*	sprintf(nomfichier_results,"results/%s_%d_%02d_%02d_%02d%s%02d.txt",\
		par->profile_name, temps.tm_year+1900, temps.tm_mon+1, temps.tm_mday, temps.tm_hour, "h", temps.tm_min);
*/
	sprintf(nomfichier_results,"%s.txt",par->profile_name);

	return 1;

}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		COMPLEX md2D_index(char* name,double lambda, char* method)
 *
 *	\brief
 */
/*---------------------------------------------------------------------------------------------*/
COMPLEX md2D_index(char* name,double lambda, char* method)
{
	char filename[SIZE_STR_BUFFER];
	double lambda_angstrom = 10*lambda;
	double cauchy_n[5],cauchy_k[5],n_power[6],k_power[6],index_n,index_k;
	COMPLEX index;
	int i,Nb_cauchy = 5;
	
	/* File name */
	snprintf(filename, SIZE_STR_BUFFER*sizeof(char), "%s_index.txt", name);
		
	/* Cauchy method */
	if (!strcmp(method,"Cauchy")){
		lire_tab(filename, "CAUCHY_N", cauchy_n, Nb_cauchy);
		lire_tab(filename, "CAUCHY_K", cauchy_k, Nb_cauchy);
		lire_tab(filename, "N_POWERS", n_power, Nb_cauchy+1);
		lire_tab(filename, "K_POWERS", k_power, Nb_cauchy+1);
SaveDbleTab2file (cauchy_n,  Nb_cauchy,"stdout", " ",10000000,"");printf("\n");
SaveDbleTab2file (n_power, Nb_cauchy,"stdout", " ",10000000,"");printf("\n");
SaveDbleTab2file (cauchy_k,  Nb_cauchy,"stdout", " ",10000000,"");printf("\n");
SaveDbleTab2file (k_power, Nb_cauchy,"stdout", " ",10000000,"");printf("\n");

		index_n = 0;
		index_k = 0;
		for (i=0;i<=Nb_cauchy-1;i++){
			index_n += cauchy_n[i]*pow(lambda_angstrom,n_power[i]);
			index_k += cauchy_k[i]*pow(lambda_angstrom,k_power[i]);
		}
printf("Milieu: %s, lambda = %f, n = %f, k = %f\n",name,lambda,index_n,index_k);
				
		index = index_n + I*index_k; 
		return index;
		
	}
	/* Lookup-Table method */
	if (!strcmp(method,"Lookup")){
		double *tab_n, *tab_k, *tab_lambda, npoints, *table_tmp;
		lire_tab(filename, "NPOINTS", &npoints, 1);
		table_tmp= (double *) malloc(sizeof(double)*3*((int)npoints));
		tab_lambda= (double *) malloc(sizeof(double)*((int)npoints));
		tab_n= (double *) malloc(sizeof(double)*((int)npoints));
		tab_k= (double *) malloc(sizeof(double)*((int)npoints));
		lire_tab(filename, "TABLE", table_tmp, 3*(int)npoints);
		/* Réarrangement en plusieurs tableaux */
		for (i=0; i<=(int)npoints -1;i++){
			tab_lambda[i] = table_tmp[3*i];
			tab_n[i] = table_tmp[3*i+1];
			tab_k[i] = table_tmp[3*i+2];
		}
		/* Recherche de l'indice de tableau pour lambda */
		int num_min = 0;
		int num_max = (int) npoints-1;
		int numero = num_max>>1;
		while(!(tab_lambda[numero] <= lambda && lambda < tab_lambda[numero+1])){
			numero=num_min+((num_max-num_min)>>1);
			if (lambda < tab_lambda[numero]){
				num_max=numero;
			}
			if (tab_lambda[numero+1] <= lambda){
				num_min=numero;
			}
			if (num_min == num_max){
				fprintf(stderr,"md2D_indice, look-up table incompatible avec lambda = %f pour %s",lambda,name);
				exit(EXIT_FAILURE);
			}
		}
		/* Interpolation linéaire */
		double l1 = tab_lambda[numero];
		double l2 = tab_lambda[numero+1];
		double n1 = tab_n[numero];
		double n2 = tab_n[numero+1];
		double k1 = tab_k[numero];
		double k2 = tab_k[numero+1];
		index_n= (lambda-l1)*(n2-n1)/(l2-l1) + n1;
		index_k= (lambda-l1)*(k2-k1)/(l2-l1) + k1;

		free(tab_lambda);
		free(tab_n);
		free(tab_k);
		free(table_tmp);
		
/*printf("lambda= %f, lambda1= %f, lambda2= %f\n",lambda,tab_lambda[numero],tab_lambda[numero+1]);
printf("lambda= %f, n = %f, k= %f\n",lambda,index_n,index_k);
*/		index = index_n + I*index_k; 

		return index;
	}
	return -1;/* ERREUR, ne devrait pas arriver là ... */
}
