/*! \file		md3D_in_out.c
 *
 *	\brief		Routines de gestion des entrées-sorties pour md3D
 *
 *  \date		../../2007
 *  \authors	Laurent ARNAUD
 */
#include "md3D_in_out.h"

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md3D_lire_param(struct Param_struct *par){
 *	
 *	\brief	Lecture des paramètres dans un fichier
 */
/*---------------------------------------------------------------------------------------------*/
int md3D_lire_param(struct Noms_fichiers *nomfichier, struct Param_struct *par){

	FILE *fp;
	int ret=0, argc=par->argc;
	double nu_super_re, nu_super_im, nu_sub_re, nu_sub_im;
	char *erreur="NO_ERROR                     ";
	char **argvcp=par->argvcp, str_profil[SIZE_STR_BUFFER], str_tmp[SIZE_STR_BUFFER];

	/* Nom de fichier param : en ligne de commande ou par défaut */
	if (lire_str_arg(nomfichier->param_file, "-param", argc, argvcp)) {
		sprintf(nomfichier->param_file,"md3D_param.txt");}

	/* Ouverture de param_file */
	if (!(fp = fopen(nomfichier->param_file,"r"))){
		fprintf(stderr, "%s ligne %d : Error, can't open %s\n",__FILE__, __LINE__,nomfichier->param_file);
		exit(EXIT_FAILURE);
	}
	
	/* Verbosity level (default=2) */
	if (lire_int_arg(&par->verbosity, "-verbosity", argc, argvcp)) {
		if (lire_int (fp, "verbosity", &(par->verbosity) )) par->verbosity = 2;}
	if (par->verbosity >=2 ) fprintf(stdout,"Reading parameters in %s\n",nomfichier->param_file);
	

	/* Lecture des paramètres, d'abord en ligne de commande, si rien en ligne de commande,
	   lecture dans param_file, sinon erreur et arret du programme */
	if (lire_str_arg(nomfichier->profile_file, "-profile_file", argc, argvcp)){
		if (lire_string (fp, "profile_file", nomfichier->profile_file)) erreur="profile_file";}
	if (lire_str_arg(par->profile_file, "-profile_file", argc, argvcp)){
		if (lire_string (fp, "profile_file", par->profile_file)) erreur="profile_file";}
	if (lire_int_arg(&par->Nx, "-Nx", argc, argvcp)) {
		if (lire_int (fp, "Nx", &(par->Nx) )) erreur="Nx";}
	if (lire_int_arg(&par->Ny, "-Ny", argc, argvcp)) {
		if (lire_int (fp, "Ny", &(par->Ny) )) erreur="Ny";}
	if (lire_dble_arg(&par->Lx, "-Lx", argc, argvcp)){
		if (lire_double (fp, "Lx", &(par->Lx) )) erreur="Lx";}
	if (lire_dble_arg(&par->Ly, "-Ly", argc, argvcp)){
		if (lire_double (fp, "Ly", &(par->Ly) )) erreur="Ly";}
	if (lire_str_arg(par->calcul_type, "-calcul_type", argc, argvcp)) {
		if (lire_string (fp, "calcul_type", par->calcul_type)) erreur="calcul_type";}
	if (lire_str_arg(par->calcul_method, "-calcul_method", argc, argvcp)) {
		if (lire_string (fp, "calcul_method", par->calcul_method)) erreur="calcul_method";}
	if (lire_dble_arg(&par->lambda, "-lambda", argc, argvcp)) {
		if (lire_double (fp, "lambda", &(par->lambda) )) erreur="lambda";}
	if (lire_str_arg(par->profile_name, "-profile_name", argc, argvcp)) {
		if (lire_string (fp, "profile_name", par->profile_name)) erreur="profile_name";}
	if (lire_dble_arg(&par->coef_h, "-coef_h", argc, argvcp)) {
		if (lire_double (fp, "coef_h", &(par->coef_h) )) par->coef_h = 1;}
	if (lire_dble_arg(&par->theta_i, "-theta_i", argc, argvcp)) {
		if (lire_double (fp, "theta_i", &(par->theta_i) )) erreur="theta_i";}
		par->theta_i *= PI/180.0;
	if (lire_dble_arg(&par->phi_i, "-phi_i", argc, argvcp)) {
		if (lire_double (fp, "phi_i", &(par->phi_i) )) erreur="phi_i";}
		par->phi_i *= PI/180.0;
	if (lire_dble_arg(&par->psi, "-psi", argc, argvcp)) {
		if (lire_double (fp, "psi", &(par->psi) )) erreur="psi";}
		par->psi *= PI/180.0;
	if (lire_str_arg(par->i_field_mode, "-i_field_mode", argc, argvcp)) {
		if (lire_string (fp, "i_field_mode", par->i_field_mode)) erreur="i_field_mode";}
	if (lire_int_arg(&par->READ_tab_NS, "-READ_tab_NS", argc, argvcp)) {
		if (lire_int (fp, "READ_tab_NS", &(par->READ_tab_NS) )) par->READ_tab_NS = 0;}
	if (par->READ_tab_NS){
		if (lire_str_arg(par->tab_NS_filename, "-tab_NS_filename", argc, argvcp)) {
			if (lire_string (fp, "tab_NS_filename", par->tab_NS_filename)) erreur="tab_NS_filename";}
	}
	if (lire_int_arg(&par->imposed_S_steps, "-imposed_S_steps", argc, argvcp)) {
		if (lire_int (fp, "imposed_S_steps", &(par->imposed_S_steps) )) par->imposed_S_steps = 0;}
	if (par->imposed_S_steps){
		if (lire_str_arg(par->imposed_S_steps_filename, "-imposed_S_steps_filename", argc, argvcp)) {
			if (lire_string (fp, "imposed_S_steps_filename", par->imposed_S_steps_filename)) erreur="imposed_S_steps_filename";}
		if (lire_int_arg(&par->N_imposed_S_steps, "-N_imposed_S_steps", argc, argvcp)) {
			if (lire_int (fp, "N_imposed_S_steps", &(par->N_imposed_S_steps) )) erreur="N_imposed_S_steps";}
	}
	if (lire_dble_arg(&nu_super_re, "-nu_super_re", argc, argvcp)) {
		if (lire_complex(fp, "nu_super", &(par->nu_super))) erreur="nu_super";
	}else{
		if (lire_dble_arg(&nu_super_im, "-nu_super_im", argc, argvcp)) {
			nu_super_im =0;}
		par->nu_super = nu_super_re +I*nu_super_im;}
	if (lire_dble_arg(&nu_sub_re, "-nu_sub_re", argc, argvcp)) {
		if (lire_complex(fp, "nu_sub", &(par->nu_sub))) erreur="nu_sub";
	}else{
		if (lire_dble_arg(&nu_sub_im, "-nu_sub_im", argc, argvcp)) {
			nu_sub_im =0;}
		par->nu_sub = nu_sub_re + I*nu_sub_im;}
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
		fprintf(stderr, "%s: Error: can't read \"%s\"\n",__FILE__,erreur);
		exit(EXIT_FAILURE);
	}

	/* Lecture des paramètres dans fichier profil */
	/* (Les paramètres liés à la nature ou indisociables du profil sont contenus dans profile_file) */
	if (par->verbosity >= 2) fprintf(stdout,"Reading parameters in %s\n",nomfichier->profile_file);
	if (!(fp = fopen(nomfichier->profile_file,"r"))){
		fprintf(stderr, "%s ligne %d : Error, can't open %s\n",__FILE__, __LINE__,nomfichier->profile_file);
		exit(EXIT_FAILURE);
	}
	erreur="NO_ERROR                     ";
	if (lire_string (fp, "profile_type", str_profil)) {
		erreur="profile_type";
	}else{
		if       (!strcmp(str_profil,"H_XY"))         {par->profile_type = H_XY;
			if (lire_int (fp, "Nprx", &(par->Nprx) )) erreur="Nprx";
			if (lire_int (fp, "Npry", &(par->Npry) )) erreur="Npry";
			par->N_layers = 0;
		}else if (!strcmp(str_profil,"MULTICOUCHES")){par->profile_type = MULTICOUCHES;
			if (lire_int (fp, "Nprx", &(par->Nprx) )) erreur="Nprx";
			if (lire_int (fp, "Npry", &(par->Npry) )) erreur="Npry";
			if (lire_int (fp, "N_layers", &(par->N_layers) )) erreur="N_layers";
		}else if (!strcmp(str_profil,"N_XYZ"))       {par->profile_type = N_XYZ;
			if (lire_int (fp, "Nprx", &(par->Nprx) )) erreur="Nprx";
			if (lire_int (fp, "Npry", &(par->Npry) )) erreur="Npry";
			if (lire_int (fp, "Nprz", &(par->Nprz) )) erreur="Nprz";
			par->N_layers = 0;
		}else if (!strcmp(str_profil,"N_XY_ZINVAR"))       {par->profile_type = N_XY_ZINVAR;
			if (lire_int (fp, "Nprx", &(par->Nprx) )) erreur="Nprx";
			if (lire_int (fp, "Npry", &(par->Npry) )) erreur="Npry";
			if (lire_int (fp, "Nprz", &(par->Nprz) )) erreur="Nprz";
			par->N_layers = 0;
		}else if (!strcmp(str_profil,"H_XY_plus_STACK")){par->profile_type = H_XY_plus_STACK;
			if (lire_int (fp, "Nprx", &(par->Nprx) )) erreur="Nprx";
			if (lire_int (fp, "Npry", &(par->Npry) )) erreur="Npry";
			par->N_layers = 0;
			if (lire_int(fp, "N_stack", &(par->N_stack))) erreur="N_stack";
			if (lire_int(fp, "n_patterned_layer", &(par->n_patterned_layer))) erreur="n_patterned_layer";
		}else                                        {erreur = "profile_type_bis";}
	}
	fclose(fp);
	
	/* Vérification de l'absence d'erreurs de lecture */
	if (strcmp(erreur,"NO_ERROR                     ")){
		fprintf(stderr, "%s: Error, can't read \"%s\"\n",__FILE__,erreur);
		exit(EXIT_FAILURE);
	}

	return ret;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md3D_affiche_valeurs_param(struct Param_struct *par, struct Noms_fichiers *nomfichier)
 *	
 *	\brief	Show parameters values on the screen
 */
/*---------------------------------------------------------------------------------------------*/
int md3D_affiche_valeurs_param(struct Param_struct *par, struct Noms_fichiers *nomfichier)
{

	/* Affichage des valeurs lues */
	if (par->verbosity >= 2){
		fprintf(stdout,"nu_super = %f + i%f\n", creal((par->nu_super)), cimag((par->nu_super)));
		fprintf(stdout,"nu_sub   = %f + i%f\n", creal((par->nu_sub)), cimag((par->nu_sub)));
		fprintf(stdout,"lambda  = %f\n",par->lambda);
		fprintf(stdout,"theta_i = %f rad (%f deg)\n",par->theta_i,par->theta_i*180.0/PI);
		fprintf(stdout,"phi_i   = %f rad (%f deg)\n",par->phi_i,par->phi_i*180.0/PI);
		fprintf(stdout,"psi     = %f rad (%f deg)\n",par->psi,par->psi*180.0/PI);
		fprintf(stdout,"Lx       = %f\n",par->Lx);
		fprintf(stdout,"Ly       = %f\n",par->Ly);
		fprintf(stdout,"h       = %f\n",par->h);
		fprintf(stdout,"coef_h  = %f\n",par->coef_h);
		fprintf(stdout,"delta_h = %f\n",par->delta_h);
		fprintf(stdout,"Nx       = %d\n",par->Nx);
		fprintf(stdout,"Ny       = %d\n",par->Ny);
		fprintf(stdout,"NS      = %d\n",par->NS);
		fprintf(stdout,"Nprx     = %d\n",par->Nprx);
		fprintf(stdout,"Npry     = %d\n",par->Npry);
		fprintf(stdout,"calcul_type    = %s\n",par->calcul_type);
		fprintf(stdout,"calcul_method  = %s\n",par->calcul_method);
		fprintf(stdout,"i_field_mode   = %s\n",par->i_field_mode);
		fprintf(stdout,"profile_name     = %s\n",par->profile_name);
		fprintf(stdout,"profile_file = %s\n",nomfichier->profile_file);
		if (par->profile_type==H_XY){
			fprintf(stdout,"profile_type    = %s\n","H_XY");
		}else if (par->profile_type==N_XY_ZINVAR){
			fprintf(stdout,"profile_type    = %s\n","N_XY_ZINVAR");
		}else{
			fprintf(stdout,"*** WARNING ***:\nprofile_type    = %s\n***************\n","UNKNOWN");
		}
		fprintf(stdout,"N_couches      = %d\n",par->N_layers);
		fflush(stdout);
	}
	return 0;
}



/*---------------------------------------------------------------------------------------------*/
/*!	\fn	int md3D_lire_profil_H_XY(const char *nom_fichier, struct Param_struct *par)
 *
 *	\brief	Reads the h(x) profile, determine the height h, gives an alert in case it is
 * 			different from the indicated value, the indicated values is prioritary.
 * 			The profile and h are also multiplied by an eventually !=1 coef_h 
 */
/*---------------------------------------------------------------------------------------------*/
int md3D_lire_profil_H_XY(const char *nom_fichier, struct Param_struct *par)
{
	int i;
	char index_name[SIZE_STR_BUFFER], h_name[SIZE_STR_BUFFER];
	double hstack_tmp;
	COMPLEX index_tmp;
	int Nprx = par->Nprx;
	int Npry = par->Npry;
	FILE *fp;
	double h_tmp, *profil = par->profil[0];
	int Npts = Nprx*Npry;

	/* Entry of k2 & invk2 for substrate and superstrate in arrays (inv)k2_layer, for compatibility with multicouches */
	par->k2_layer[0] = (par->k_super)*(par->k_super);
	par->k2_layer[par->N_layers+1] = (par->k_sub)*(par->k_sub);
	par->invk2_layer[0] = 1/((par->k_super)*(par->k_super));
	par->invk2_layer[par->N_layers+1] = 1/((par->k_sub)*(par->k_sub));

	if (par->profile_type == H_XY_plus_STACK) {
		if (!(fp = fopen(nom_fichier,"r"))){
			fprintf(stderr, "%s line %d: Error, can't open %s\n",__FILE__, __LINE__,nom_fichier);
			exit(EXIT_FAILURE);
		}
		for (i=0; i<=par->N_stack-1; i++){
			if (i != par->n_patterned_layer){
				sprintf(index_name,"n%d",i);
				if (lire_complex(fp, index_name, &index_tmp)) fprintf(stderr, "%s line %d: Error, can't read %s\n",__FILE__, __LINE__,index_name);
				par->nu_stack[i] = index_tmp;
				sprintf(h_name,"h%d",i);
				if (lire_double(fp, h_name, &hstack_tmp))  fprintf(stderr, "%s line %d: Error, can't read %s\n",__FILE__, __LINE__,h_name);
				par->h_stack[i] = hstack_tmp;
			}
		}
		fclose(fp);
	}

	/* Reading the profile */
	if (par->verbosity >= 2) fprintf(stdout,"Reading the profile %s : ",nom_fichier);
	if (lire_tab(nom_fichier, "profil", profil, Npts) == 0) {
		if (par->verbosity >= 2) fprintf(stdout,"OK\n");
	}else{
		fprintf(stderr,"Profile reading ERROR\n");
		exit(EXIT_FAILURE);
	}

	/* Determination of the profile height from read h(x) points */
	double h_min = profil[0];
	double h_max = profil[0];
	double eps = 1.0e-10; 
	for (i=0;i<=Npts-1;i++){
		h_min = MIN(h_min,profil[i]);
		h_max = MAX(h_max,profil[i]);
	}
	h_tmp = h_max - h_min;

	/* Checking that h(detemined) = h(indicated) and setting par->h to determined value if h = AUTO */
	if (par->h == AUTO) {
		par->h = h_tmp;
	}else if (fabs(par->h - h_tmp) > eps*par->h) {
		if (par->verbosity>-1) {
			fprintf(stderr,  "+---------------------------------------------------------------------\
				\n|                        CAUTION !\
				\n| h determined for %s is %f and h indicated %f !\
				\n+---------------------------------------------------------------------\
				\n",nom_fichier,h_tmp, par->h);
		}
	}
	/* Normalising beetween 0 and h (instead of [h0,h0+h]*/
	for (i=0; i<=Npts-1; i++) {
		profil[i] -= h_min;
	}
		
	/* Multiplication by coef_h x h_wanted / h_determined */
	if (h_tmp != 0.0) {
		for (i=0; i<=Npts-1; i++) {
			profil[i] *= par->coef_h * par->h / h_tmp;
		}
	}
	/* Multiplication of par->h by coef_h */
	par->h *= par->coef_h;
	
	if (par->verbosity >= 2) fprintf(stdout,"Profile normalisation between 0 and (h x coef_h) : OK\n");


	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md3D_lire_profil_MULTI(const char *nom_fichier, struct Param_struct *par)

 *
 *	\brief	Fonction lisant les valeurs décrivant un profil multi
 *
 */
/*---------------------------------------------------------------------------------------------*/
int md3D_lire_profil_MULTI(const char *nom_fichier, struct Param_struct *par)
{

	int i, nx, n_layer;
	char nom_indice[SIZE_STR_BUFFER];
	char *erreur="NO_ERROR                     ";
	FILE *fp;
	double *profil_tmp, **profil = par->profil;
	int Nprx = par->Nprx;
	int Npry = par->Npry;
	int Npts = Nprx*Npry;
	int N_layers = par->N_layers;
	double h_tmp;
	COMPLEX indice;
	
	profil_tmp = (double *) malloc(sizeof(double)*Npts*(N_layers+1));

	if (!(fp = fopen(nom_fichier,"r"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nom_fichier);
		exit(EXIT_FAILURE);
	}
	/* Lecture des indices des couches et calculs des k2 et 1/k2 */
	for (i=1; i<=N_layers; i++){
		sprintf(nom_indice,"n%d",i);
		if (lire_complex(fp, nom_indice, &indice)) erreur=nom_indice;
		par->k2_layer[i]    = (indice*par->k_super/par->nu_super)*(indice*par->k_super/par->nu_super); 
		par->invk2_layer[i] = 1/par->k2_layer[i]; 
	}
	par->k2_layer[0]             = (par->k_super)*(par->k_super);
	par->k2_layer[N_layers+1]    = (par->k_sub)*(par->k_sub);
	par->invk2_layer[0]          = 1/par->k2_layer[0];
	par->invk2_layer[N_layers+1] = 1/par->k2_layer[N_layers+1];

	fclose(fp);	
	/* Vérification de l'absence d'erreurs de lecture */
	if (strcmp(erreur,"NO_ERROR                     ")){
		fprintf(stderr, "%s: Error, can't read \"%s\"\n",__FILE__,erreur);
		exit(EXIT_FAILURE);
	}

	
	/* Lecture des profils */
	/* On lit comme un seul tableau, en lisant les lignes les unes à la suite des autres,  */
	/* en considérant qu'une colonne représente les coordonnées d'une interface. On sépare */
	/* ensuite les données en autant de tableaux qu'il y a d'interfaces                    */
	if (par->verbosity >= 2) fprintf(stdout,"Lecture du profil %s : ",nom_fichier);
	if (lire_tab(nom_fichier, "profil", profil_tmp, Npts*(N_layers+1)) == 0) {
		if (par->verbosity >= 2) fprintf(stdout,"OK\n");
	}else{
		fprintf(stderr,"ERREUR de lecture du profil\n");
		exit(EXIT_FAILURE);
	}
/**************************************************/
	/* Détermination de la hauteur du profil à partir des points h(x) */
	double h_min = profil_tmp[0];
	double h_max = profil_tmp[0];
	double eps = 1.0e-10; 
	for (i=0;i<=Npts*(N_layers+1)-1;i++){
		h_min = MIN(h_min,profil_tmp[i]);
		h_max = MAX(h_max,profil_tmp[i]);
	}
	h_tmp = h_max - h_min;

	/* Vérification que h(deteminé) = h(indiqué) et attribution si h = AUTO */
	if (par->h == AUTO) {
		par->h = h_tmp;
	}else if (fabs(par->h - h_tmp) > eps*par->h) {
		fprintf(stderr,  "+---------------------------------------------------------------------\
				\n|                        CAUTION !\
				\n| h determined for %s is %f and h indicated %f !\
				\n+---------------------------------------------------------------------\
				\n",nom_fichier,h_tmp, par->h);
	}
	/* Normalisation entre 0 et h (au lieu de [h0,h0+h]*/
	for (i=0; i<=Npts*(N_layers+1)-1; i++) {
		profil_tmp[i] -= h_min;
	}
	/* Vérification que les profils ne se chevauchent pas */

			
	/* Multiplication par coef_h x h_voulu / h_mesuré */
	if (h_tmp != 0.0) {
		for (i=0; i<=Npts*(N_layers+1)-1; i++) {
			profil_tmp[i] *= par->coef_h * par->h / h_tmp;
		}
	}
	/* Multiplication de par->h par coef_h */
	par->h *= par->coef_h;
	
	if (par->verbosity >= 2) fprintf(stdout,"Normalisation du profil entre 0 et (h x coef_h) : OK\n");

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
	double eps2 = par->h*1e-10;
	for (nx=0; nx<=Npts-1;nx++){
		/* "haut du superstrat", z=h */
		profil[0][nx] = par->h+eps2;
		/* Couches */
		for (n_layer=1; n_layer<=N_layers+1; n_layer++){
			profil[n_layer][nx] = profil_tmp[nx*(N_layers+1)+n_layer-1];
		}
		/* "Bas du substrat", z=h */
		profil[N_layers+2][nx] = 0-eps2;
	}
			
	free(profil_tmp);

	return 0;
}

/*!---------------------------------------------------------------------------------------------
 * \fn int md3D_lire_profil_N_XYZ(const char *nom_fichier, struct Param_struct *par)
 *
 * \brief Read profil in case of n(x,y,z) index distribution
 *
 * Two ways of entering the refractive index in the stack for this function:
 * - first way: enter an array of values (1,2,3,...) corresponding to COMPLEX refractive indices
 *   n1,n2,n3,... which values are defined in the profile file.
 *   Ex.:   n1 = 1.5 +i0.5
 *          n2 = 1.0 +i0.0
 *          n_xyz = 1 1 2 1 2 2 ...
 * - second way: enter two arrays, one for the real part and one for the imaginary part (facultative)
 *   of the refractive index.
 *   Ex.:   Re_n_xyz = 1.5 1.5 1.0 1.5 1.5 ...
 *          Im_n_xyz = 0.5 0.5 0.0 0.5 0.5 ...
 *---------------------------------------------------------------------------------------------*/
int md3D_lire_profil_N_XYZ(const char *nom_fichier, struct Param_struct *par)
{
	int i, i_max, nx, nz, value_is_attributed, Nprx=par->Nprx, Npry=par->Npry, Nprz=par->Nprz, Nptxy = Nprx*Npry;
	FILE *fp;	
	COMPLEX **n_xyz = par->n_xyz;
	char nu_name[SIZE_STR_BUFFER];
	
	/* Check Nprz for ZINVAR case */
	if (par->profile_type == N_XY_ZINVAR && Nprz != 1){
		fprintf(stderr, "%s line %d: ERROR, Nprz must = 1 for n_xy_zinvar. For not z-invariant profile, use profile_type = n_xyz (with uppercase).\n",__FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	/* Memory allocations */
	double *Re_n_xyz = malloc(sizeof(double)*Nptxy*Nprz);
	double *Im_n_xyz = malloc(sizeof(double)*Nptxy*Nprz);
	double *int_n_xyz = malloc(sizeof(double)*Nptxy*Nprz);
	COMPLEX *nu = malloc(sizeof(COMPLEX)*SIZE_INT_BUFFER);
	

	/* Reading index distribution as n_xyz, plus values n1, n2,... */
	if (lire_tab(nom_fichier, "n_xyz", int_n_xyz, Nptxy*Nprz) == 0) {
		if (par->verbosity >= 2) fprintf(stdout,"Reading the N_XYZ index distribution of the type n_xyz in %s... OK\n",nom_fichier);
		/* Reading n1, n2, n3,... indices COMPLEX values */
		if (!(fp = fopen(nom_fichier,"r"))){
			fprintf(stderr, "%s line %d: Error, can't open %s\n",__FILE__, __LINE__,nom_fichier);
			exit(EXIT_FAILURE);
		}
		for (i=1;i>0;i++){
 			snprintf(nu_name, SIZE_STR_BUFFER*sizeof(char), "n%d",i);
			if (lire_complex(fp, nu_name, &nu[i-1])){
				i_max = i;
				break;
			}
		}
		if (i_max > SIZE_INT_BUFFER){
			fprintf(stderr, "%s line %d: ERROR, can't handle more than %d indices (Edit the source code if you need more).\n",__FILE__, __LINE__,SIZE_INT_BUFFER);
			exit(EXIT_FAILURE);
		}
		fclose(fp);
		/* Attributing the COMPLEX values to the index array */
		for (nz=0; nz<=Nprz-1; nz++){
			for (nx=0; nx<=Nptxy-1; nx++){
				value_is_attributed = 0;
				for (i=1;i<=i_max;i++){
					if ((int)int_n_xyz[nx+Nptxy*nz] == i){
						n_xyz[nz][nx] = nu[i-1];
						value_is_attributed = 1;
						break;
					}
				}
				if (!value_is_attributed){
					fprintf(stderr, "%s line %d: Error, can't attribute value to n_xyz for element number %d\n",__FILE__, __LINE__,nx);
					exit(EXIT_FAILURE);
				}
			}
		}
	}else if (lire_tab(nom_fichier, "Re_n_xyz", Re_n_xyz, Nptxy*Nprz) == 0) { /* Reading real part of n_xyz */
		if (par->verbosity >= 2) fprintf(stdout,"Reading the N_XYZ index distribution %s:\nreal part, Re_n_xyz... OK\n",nom_fichier);
		/* Reading imaginary part */
		if (lire_tab(nom_fichier, "Im_n_xyz", Im_n_xyz, Nptxy*Nprz) == 0) {
			if (par->verbosity >= 2) fprintf(stdout,"Reading the N_XYZ index distribution %s:\nImaginary part, Im_n_xyz... OK\n",nom_fichier);
		}else{
			fprintf(stderr,"NO IMAGINARY PART: Dielectric profile\n");
			for (nx=0; nx<=Nptxy*Nprz-1; nx++){
				Im_n_xyz[nx] = 0;
			}
		}
		/* Creating COMPLEX matrix n_xyz */
		for (nz=0; nz<=Nprz-1; nz++){
			for (nx=0; nx<=Nptxy-1; nx++){
				n_xyz[nz][nx] = Re_n_xyz[nx+Nptxy*nz] + I*Im_n_xyz[nx+Nptxy*nz];
			}
		}
	}else{
		fprintf(stderr, "%s line %d: ERROR, can't read either n_xyz or Re_n_xyz in %s\n",__FILE__, __LINE__,nom_fichier);
		exit(EXIT_FAILURE);
	}
	

	free(Re_n_xyz);
	free(Im_n_xyz);
	free(int_n_xyz);
	free(nu);

	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn	
 *
 *	\brief	Writing results
 */
/*---------------------------------------------------------------------------------------------*/
int md3D_ecrire_results(char *nom_fichier, struct Param_struct *par, struct Efficacites_struct *eff)
{
	FILE *fp;



	/* Ouverture du fichier */
	if (!(fp = fopen(nom_fichier,"w"))){
		fprintf(stderr, "%s ligne %d : Error, can't open %s\n",__FILE__, __LINE__,nom_fichier);
		exit(EXIT_FAILURE);
	}

	fprintf(fp,"#    \n");
	fprintf(fp,"#    \n\n");

	fprintf(fp,    "#------------ Parameters ------------\n");
	fprintf(fp,"calcul_type   : %s\n", par->calcul_type);	
	fprintf(fp,"calcul_method : %s\n", par->calcul_method);	
	fprintf(fp,"theta_i = %f deg\n",par->theta_i*180.0/PI);
	fprintf(fp,"phi_i   = %f deg\n",par->phi_i*180.0/PI);
	fprintf(fp,"psi     = %f deg\n",par->psi*180.0/PI);

	fprintf(fp,"nu_super = %f + i%f\n", creal((par->nu_super)), cimag((par->nu_super)));
	fprintf(fp,"nu_sub   = %f + i%f\n", creal((par->nu_sub)), cimag((par->nu_sub)));
	fprintf(fp,"Lx       = %f\n",par->Lx);
	fprintf(fp,"Ly       = %f\n",par->Ly);
	fprintf(fp,"h       = %f\n",par->h);
	fprintf(fp,"coef_h  = %f\n",par->coef_h);
	fprintf(fp,"lambda  = %f\n",par->lambda);
	fprintf(fp,"Nx       = %d\n",par->Nx);
	fprintf(fp,"Ny       = %d\n",par->Ny);
	fprintf(fp,"NS      = %d\n",par->NS);
	fprintf(fp,"delta_h = %f\n",par->delta_h);
	fprintf(fp,"N_steps   = %d\n",par->N_steps);
/*	fprintf(fp,"profile_file = %s\n",profile_file);*/
	fprintf(fp,"Nprx = %d\n",par->Nprx);
	fprintf(fp,"Npry = %d\n",par->Npry);
	fprintf(fp,"Calcul_duration = %f s \n", md_chrono(par));

	fprintf(fp,  "\n#------------ Results -------------\n");
	fprintf(fp,"sum_eff           = % 1.6le\n",eff->sum_eff);
	fprintf(fp,"one_minus_sum_eff = % 1.6e\n",1.0-eff->sum_eff);
	fprintf(fp,"somm_eff_r        = % 1.6e\n",eff->sum_eff_r);
	fprintf(fp,"somm_eff_t        = % 1.6e\n",eff->sum_eff_t);
	fprintf(fp,"\neff_r = "); ecrire_dble_tab(fp, eff->eff_r[0], (eff->Nxmax_super-eff->Nxmin_super+1)*(eff->Nymax_super-eff->Nymin_super+1), " ", (eff->Nxmax_super-eff->Nxmin_super+1),"\n");
	fprintf(fp,"\nnx_eff_r = "); ecrire_dble_tab(fp, eff->nx_eff_r[0], (eff->Nxmax_super-eff->Nxmin_super+1)*(eff->Nymax_super-eff->Nymin_super+1), " ", (eff->Nxmax_super-eff->Nxmin_super+1),"\n");
	fprintf(fp,"\nny_eff_r = "); ecrire_dble_tab(fp, eff->ny_eff_r[0], (eff->Nxmax_super-eff->Nxmin_super+1)*(eff->Nymax_super-eff->Nymin_super+1), " ", (eff->Nxmax_super-eff->Nxmin_super+1),"\n");
	fprintf(fp,"\neff_t = "); ecrire_dble_tab(fp, eff->eff_t[0], (eff->Nxmax_sub-eff->Nxmin_sub+1)*(eff->Nymax_sub-eff->Nymin_sub+1), " ", (eff->Nxmax_sub-eff->Nxmin_sub+1),"\n");
	fprintf(fp,"\nnx_eff_t = "); ecrire_dble_tab(fp, eff->nx_eff_t[0], (eff->Nxmax_sub-eff->Nxmin_sub+1)*(eff->Nymax_sub-eff->Nymin_sub+1), " ", (eff->Nxmax_sub-eff->Nxmin_sub+1),"\n");
	fprintf(fp,"\nny_eff_t = "); ecrire_dble_tab(fp, eff->ny_eff_t[0], (eff->Nxmax_sub-eff->Nxmin_sub+1)*(eff->Nymax_sub-eff->Nymin_sub+1), " ", (eff->Nxmax_sub-eff->Nxmin_sub+1),"\n");

	fclose(fp);
	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md3D_genere_nom_fichier_results(char *nomfichier_results, struct Param_struct *par)
 *
 *	\brief
 */
/*---------------------------------------------------------------------------------------------*/
int md3D_genere_nom_fichier_results(char *nomfichier_results, struct Param_struct *par)
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
/*!	\fn		COMPLEX md3D_index(char* name,double lambda, char* method)
 *
 *	\brief
 */
/*---------------------------------------------------------------------------------------------*/
COMPLEX md3D_index(char* name,double lambda, char* method)
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
/*SaveDbleTab2file (cauchy_n,  Nb_cauchy,"stdout", " ");printf("\n");
SaveDbleTab2file (n_power, Nb_cauchy,"stdout", " ");printf("\n");
SaveDbleTab2file (cauchy_k,  Nb_cauchy,"stdout", " ");printf("\n");
SaveDbleTab2file (k_power, Nb_cauchy,"stdout", " ");printf("\n");*/

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
				fprintf(stderr,"md3D_indice, look-up table incompatible avec lambda = %f pour %s",lambda,name);
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
