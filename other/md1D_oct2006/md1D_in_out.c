/*! \file		md1D_in_out.c
 *
 *	\brief		Routines de gestion des entrées-sorties pour md1D
 *
 *	\version	0.1
 *  \date		../../2004
 *  \authors	Laurent ARNAUD
 */

#include "md1D_in_out.h"

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_lire_param(struct Param_struct *par){
 *	
 *	\brief	Lecture des paramètres dans un fichier
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_lire_param(struct Noms_fichiers *nomfichier, struct Param_struct *par){

	FILE *fp;
	int ret=0, argc=par->argc;
	double n_super_re, n_super_im, n_sub_re, n_sub_im;
	char *erreur="NO_ERROR                     ";
	char **argvcp=par->argvcp, str_profil[SIZE_STR_BUFFER], str_tmp[SIZE_STR_BUFFER];

	/* Nom de fichier param : en ligne de commande ou par défaut */
	if (lire_str_arg(nomfichier->fichier_param, "-param", argc, argvcp)) {
		sprintf(nomfichier->fichier_param,"md1D_param.txt");}

	if (par->verbose) fprintf(stdout,"Lecture des paramètres dans %s\n",nomfichier->fichier_param);
	
	/* Ouverture de fichioer_param */
	if (!(fp = fopen(nomfichier->fichier_param,"r"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nomfichier->fichier_param);
		exit(EXIT_FAILURE);
	}

	/* Lecture des paramètres, d'abord en ligne de commande, si rien en ligne de commande,
	   lecture dans fichier_param, sinon erreur et arret du programme */
	if (lire_str_arg(nomfichier->fichier_profil, "-fichier_profil", argc, argvcp)){
		if (lire_string (fp, "fichier_profil", nomfichier->fichier_profil)) erreur="fichier_profil";}
	if (lire_dble_arg(&par->L, "-L", argc, argvcp)){
		if (lire_double (fp, "L", &(par->L) )) erreur="L";}
	if (lire_str_arg(par->calcul_type, "-calcul_type", argc, argvcp)) {
		if (lire_string (fp, "calcul_type", par->calcul_type)) erreur="calcul_type";}
	if (lire_dble_arg(&par->lambda, "-lambda", argc, argvcp)) {
		if (lire_double (fp, "lambda", &(par->lambda) )) erreur="lambda";}
	if (lire_str_arg(par->nom_profil, "-nom_profil", argc, argvcp)) {
		if (lire_string (fp, "nom_profil", par->nom_profil)) erreur="nom_profil";}
	if (lire_dble_arg(&par->coef_h, "-coef_h", argc, argvcp)) {
		if (lire_double (fp, "coef_h", &(par->coef_h) )) par->coef_h = 1;}
	if (lire_int_arg(&par->mode_extract_S, "-mode_extract_S", argc, argvcp)) {
		if (lire_int (fp, "mode_extract_S", &(par->mode_extract_S) )) par->mode_extract_S = 0;}
	if (par->mode_extract_S){ 	
	if (lire_dble_arg(&par->h_extract_S, "-h_extract_S", argc, argvcp)) {
		if (lire_double (fp, "h_extract_S", &(par->h_extract_S) )) erreur="h_extract_S";}}	
	if (lire_int_arg(&par->READ_MAT_S, "-READ_MAT_S", argc, argvcp)) {
		if (lire_int (fp, "READ_MAT_S", &(par->READ_MAT_S) )) par->READ_MAT_S = 0;}
	if (par->READ_MAT_S){ 	
	if (lire_str_arg(par->mat_S_file, "-mat_S_file", argc, argvcp)) {
		if (lire_string (fp, "mat_S_file", par->mat_S_file)) erreur="mat_S_file";}}	
	if (lire_dble_arg(&par->theta_i, "-theta_i", argc, argvcp)) {
		if (lire_double (fp, "theta_i", &(par->theta_i) )) erreur="theta_i";}
		par->theta_i *= PI/180.0;
	if (lire_str_arg(par->i_field_mode, "-i_field_mode", argc, argvcp)) {
		if (lire_string (fp, "i_field_mode", par->i_field_mode)) erreur="i_field_mode";}
	if (!strcmp(par->calcul_type,"ALEAT_T_ELLIPSO")){
		if (lire_dble_arg(&par->L_segment, "-coef_h", argc, argvcp)) {
			if (lire_double (fp, "L_segment", &(par->L_segment) )) erreur="L_segment";}
		if (lire_dble_arg(&par->ecart_type_segment, "-ecart_type_segment", argc, argvcp)) {
			if (lire_double (fp, "ecart_type_segment", &(par->ecart_type_segment) )) erreur="ecart_type_segment";}
		if (lire_dble_arg(&par->h_total_aleat_T, "-h_total_aleat_T", argc, argvcp)) {
			if (lire_double (fp, "h_total_aleat_T", &(par->h_total_aleat_T) )) erreur="h_total_aleat_T";}
	}
	if (lire_str_arg(str_tmp, "-pola", argc, argvcp)) {
		if (lire_string (fp, "pola", str_tmp)) erreur="pola";}
		if      (!strcmp(str_tmp,"TE")){par->pola = TE;}
		else if (!strcmp(str_tmp,"TM")){par->pola = TM;}
		else                            {erreur = "pola";}	
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
		fprintf(stderr, "%s : Erreur, probleme de lecture de \"%s\"\n",__FILE__,erreur);
		exit(EXIT_FAILURE);
	}

	/* Lecture des paramètres dans fichier profil */
	/* (Les paramètres liés à la nature ou indisociables du profil sont contenus dans fichier_profil) */
	if (par->verbose) fprintf(stdout,"Lecture des paramètres dans %s\n",nomfichier->fichier_profil);
	if (!(fp = fopen(nomfichier->fichier_profil,"r"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nomfichier->fichier_profil);
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
		fprintf(stderr, "%s : Erreur, probleme de lecture de \"%s\"\n",__FILE__,erreur);
		exit(EXIT_FAILURE);
	}

	return ret;
}

/*---------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------*/
int md1D_affiche_valeurs_param(struct Param_struct *par, struct Noms_fichiers *nomfichier)
{

	/* Affichage des valeurs lues */
	if (par->verbose){
		fprintf(stdout,"n_super = %f + i%f\n", creal((par->n_super)), cimag((par->n_super)));
		fprintf(stdout,"n_sub   = %f + i%f\n", creal((par->n_sub)), cimag((par->n_sub)));
		fprintf(stdout,"lambda  = %f\n",par->lambda);
		fprintf(stdout,"theta_i = %f rad (%f°)\n",par->theta_i,par->theta_i*180.0/PI);
		fprintf(stdout,"L       = %f\n",par->L);
		fprintf(stdout,"h       = %f\n",par->h);
		fprintf(stdout,"coef_h  = %f\n",par->coef_h);
		fprintf(stdout,"delta_h = %f\n",par->delta_h);
		fprintf(stdout,"N       = %d\n",par->N);
		fprintf(stdout,"NS      = %d\n",par->NS);
		fprintf(stdout,"N_x     = %d\n",par->N_x);
		fprintf(stdout,"pola    : %s\n",(par->pola==TE ? "TE" : "TM"));
		fprintf(stdout,"calcul_type    = %s\n",par->calcul_type);
		fprintf(stdout,"i_field_mode   = %s\n",par->i_field_mode);
		fprintf(stdout,"nom_profil     = %s\n",par->nom_profil);
		fprintf(stdout,"fichier_profil = %s\n",nomfichier->fichier_profil);
		fprintf(stdout,"type_profil    = %s\n",(par->type_profil==H_X ? "H_X" :
		                                        (par->type_profil==N_XYZ ? "N_XYZ":"MULTICOUCHES")));
		fprintf(stdout,"N_couches      = %d\n",par->N_layers);
		fflush(stdout);
	}
	return 0;
}



/*---------------------------------------------------------------------------------------------*/
/*!	\fn	int md1D_lire_profil_H_X(const char *nom_fichier, struct Param_struct *par)
 *
 *	\brief	Fonction lisant les valeurs décrivant un profil h(x) dans un fichier. \n
 *		Les valeurs stockées dans le fichier doivent varier entre 0 et 1,     \n
 *		les valeurs lues sont multipliées par h, pour avoir un profil variant \n
 *		entre 0 et h.
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_lire_profil_H_X(const char *nom_fichier, struct Param_struct *par)
{

	int i;
	int N_x = par->N_x;
	double h_tmp, *profil = par->profil[0];

	/* Entrée des k2 et invk2 du sub et du super dans les tableaux (inv)k2_layer, pour compatibilité avec multicouches */
	par->k2_layer[0] = (par->n_super*par->k0)*(par->n_super*par->k0);
	par->k2_layer[par->N_layers+1] = (par->n_sub*par->k0)*(par->n_sub*par->k0);
	par->invk2_layer[0] = 1/((par->n_super*par->k0)*(par->n_super*par->k0));
	par->invk2_layer[par->N_layers+1] = 1/((par->n_sub*par->k0)*(par->n_sub*par->k0));


	/* Lecture du profil */
	if (par->verbose) fprintf(stdout,"Lecture du profil %s : ",nom_fichier);
	if (lire_tab(nom_fichier, "profil", profil, N_x) == 0) {
		if (par->verbose) fprintf(stdout,"OK\n");
	}else{
		fprintf(stderr,"ERREUR de lecture du profil\n");
		exit(EXIT_FAILURE);
	}

	/* Détermination de la hauteur du profil à partir des points h(x) */
	double h_min = profil[0];
	double h_max = profil[0];
	double eps = 1.0e-10; 
	for (i=0;i<=par->N_x-1;i++){
		h_min = MIN(h_min,profil[i]);
		h_max = MAX(h_max,profil[i]);
	}
	h_tmp = h_max - h_min;

	/* Vérification que h(deteminé) = h(indiqué) et attribution si h = AUTO */
	if (par->h == AUTO) {
		par->h = h_tmp;
	}else if (fabs(par->h - h_tmp) > eps*par->h) {
		fprintf(stderr,  "+---------------------------------------------------------------------\
				\n|                        ATTENTION !\
				\n| h determiné pour %s vaut %f et h indiqué %f !\
				\n+---------------------------------------------------------------------\
				\n",nom_fichier,h_tmp, par->h);
	}
	/* Normalisation entre 0 et h (au lieu de [h0,h0+h]*/
	for (i=0; i<=N_x-1; i++) {
		profil[i] -= h_min;
	}
		
	/* Multiplication par coef_h x h_voulu / h_mesuré */
	if (h_tmp != 0.0) {
		for (i=0; i<=N_x-1; i++) {
			profil[i] *= par->coef_h * par->h / h_tmp;
		}
	}
	/* Multiplication de par->h par coef_h */
	par->h *= par->coef_h;
	
	if (par->verbose) fprintf(stdout,"Normalisation du profil entre 0 et (h x coef_h) : OK\n");


	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_lire_profil_MULTI(const char *nom_fichier, struct Param_struct *par)

 *
 *	\brief	Fonction lisant les valeurs décrivant un profil h(x) dans un fichier. \n
 *		Les valeurs stockées dans le fichier doivent varier entre 0 et 1,     \n
 *		les valeurs lues sont multipliées par h, pour avoir un profil variant \n
 *		entre 0 et h.
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_lire_profil_MULTI(const char *nom_fichier, struct Param_struct *par)
{

	int i, nx, n_layer;
	char nom_indice[SIZE_STR_BUFFER];
	char *erreur="NO_ERROR                     ";
	FILE *fp;
	double *profil_tmp, **profil = par->profil;
	int N_x = par->N_x;
	int N_layers = par->N_layers;
	double h_tmp;
	complex indice;
	
	profil_tmp = (double *) malloc(sizeof(double)*N_x*(N_layers+1));

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
	if (lire_tab(nom_fichier, "profil", profil_tmp, N_x*(N_layers+1)) == 0) {
		if (par->verbose) fprintf(stdout,"OK\n");
	}else{
		fprintf(stderr,"ERREUR de lecture du profil\n");
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
	h_tmp = h_max - h_min;

	/* Vérification que h(deteminé) = h(indiqué) et attribution si h = AUTO */
	if (par->h == AUTO) {
		par->h = h_tmp;
	}else if (fabs(par->h - h_tmp) > eps*par->h) {
		fprintf(stderr,  "+---------------------------------------------------------------------\
				\n|                        ATTENTION !\
				\n| h determiné pour %s vaut %f et h indiqué %f !\
				\n+---------------------------------------------------------------------\
				\n",nom_fichier,h_tmp, par->h);
	}
	/* Normalisation entre 0 et h (au lieu de [h0,h0+h]*/
	for (i=0; i<=N_x*(N_layers+1)-1; i++) {
		profil_tmp[i] -= h_min;
	}
	/* Vérification que les profils ne se chevauchent pas */

			
	/* Multiplication par coef_h x h_voulu / h_mesuré */
	if (h_tmp != 0.0) {
		for (i=0; i<=N_x*(N_layers+1)-1; i++) {
			profil_tmp[i] *= par->coef_h * par->h / h_tmp;
		}
	}
	/* Multiplication de par->h par coef_h */
	par->h *= par->coef_h;
	
	if (par->verbose) fprintf(stdout,"Normalisation du profil entre 0 et (h x coef_h) : OK\n");

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
	if (par->verbose) fprintf(stdout,"Normalisation du profil entre 0 et h : OK\n");
*/
/**************************************************/
	
	/* Réarrangement en plusieurs tableaux */
	double eps2 = par->h*1e-10;
	for (nx=0; nx<=N_x-1;nx++){
		/* "haut du superstrat", z=0 */
		profil[0][nx] = 0-eps2;
		/* Couches */
		for (n_layer=1; n_layer<=N_layers+1; n_layer++){
			profil[n_layer][nx] = profil_tmp[nx*(N_layers+1)+n_layer-1];
		}
		/* "Bas du substrat", z=h */
		profil[N_layers+2][nx] = par->h+eps2;
	}
			
	free(profil_tmp);

	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn	int md1D_lire_profil_N_XYZ(const char *nom_fichier, struct Param_struct *par)
 *
 *	\brief	Fonction lisant les valeurs de indice(x,z) décrivant un profil volumique
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_lire_profil_N_XYZ(const char *nom_fichier, struct Param_struct *par)
{
	int nx, nz, N_x=par->N_x, N_z=par->N_z;
	complex **n_xyz = par->n_xyz;
	

	/* Allocations de mémoire pour variables temporaires */
	double *Re_n_xyz = malloc(sizeof(double)*N_x*N_z);
	double *Im_n_xyz = malloc(sizeof(double)*N_x*N_z);

	/* Lecture de la partie réelle */
	if (par->verbose) fprintf(stdout,"Lecture du profil %s :\npartie réelle : ",nom_fichier);
	if (lire_tab(nom_fichier, "Re_n_xyz", Re_n_xyz, N_x*N_z) == 0) {
		if (par->verbose) fprintf(stdout,"OK\n");
	}else{
		fprintf(stderr,"ERREUR de lecture du profil\n");
		exit(EXIT_FAILURE);
	}
	/* Lecture de la partie imaginaire */
	if (par->verbose) fprintf(stdout,"partie imaginaire : ");fflush(stdout);
	if (lire_tab(nom_fichier, "Im_n_xyz", Im_n_xyz, N_x*N_z) == 0) {
		if (par->verbose) fprintf(stdout,"OK\n");
	}else{
		fprintf(stderr,"PAS DE PARTIE IMAGINAIRE (profil diélectrique)\n");
		for (nx=0; nx<=N_x*N_z-1; nx++){
			Im_n_xyz[nx] = 0;
		}

	}
	
	/* Création de la matrice complexe n_xyz */
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
/*---------------------------------------------------------------------------------------------*/
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
	fprintf(fp,"calcul_type : %s\n", par->calcul_type);	
	fprintf(fp,"Pola    : %s\n",(par->pola==TE ? "TE" : "TM"));	
	fprintf(fp,"n_super = %f + i%f\n", creal((par->n_super)), cimag((par->n_super)));
	fprintf(fp,"n_sub   = %f + i%f\n", creal((par->n_sub)), cimag((par->n_sub)));
	fprintf(fp,"L       = %f\n",par->L);
	fprintf(fp,"h       = %f\n",par->h);
	fprintf(fp,"coef_h  = %f\n",par->coef_h);
	fprintf(fp,"lambda  = %f\n",par->lambda);
	fprintf(fp,"theta_i = %f\n",par->theta_i);
	fprintf(fp,"N       = %d\n",par->N);
	fprintf(fp,"NS      = %d\n",par->NS);
	fprintf(fp,"delta_h = %f\n",par->delta_h);
	fprintf(fp,"Nstep   = %d\n",par->Nstep);
/*	fprintf(fp,"fichier_profil = %s\n",fichier_profil);*/
	fprintf(fp,"N_x = %d\n",par->N_x);

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

	fprintf(fp,"\nRe_Ai = "); ecrire_cplx_tab(fp, par->Ai, 2*par->N+1, RE, " ", LMAX,"\n");
	fprintf(fp,"\nIm_Ai = "); ecrire_cplx_tab(fp, par->Ai, 2*par->N+1, IM, " ", LMAX,"\n");
	fprintf(fp,"\nRe_A0 = "); ecrire_cplx_tab(fp, par->A0, 2*par->N+1, RE, " ", LMAX,"\n");
	fprintf(fp,"\nIm_A0 = "); ecrire_cplx_tab(fp, par->A0, 2*par->N+1, IM, " ", LMAX,"\n");
	fprintf(fp,"\nRe_Ah = "); ecrire_cplx_tab(fp, par->Ah, 2*par->N+1, RE, " ", LMAX,"\n");
	fprintf(fp,"\nIm_Ah = "); ecrire_cplx_tab(fp, par->Ah, 2*par->N+1, IM, " ", LMAX,"\n");

	
/*	fprintf(fp,"\n\n#--------------- Profil ---------------\n");
	fprintf(fp,"profil = "); ecrire_dble_tab(fp, par->profil[0], par->N_x, " ", LMAX,"\n");

	fprintf(fp,    "\n\n");
*/
	fclose(fp);
	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_genere_nom_fichier_results(char *nomfichier_results, struct Param_struct *par)
 *
 *	\brief
 */
/*---------------------------------------------------------------------------------------------*/
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

complex md1D_indice(char* name,double lambda, char* method)
{
	char filename[SIZE_STR_BUFFER];
	double lambda_angstrom = 10*lambda;
	double cauchy_n[5],cauchy_k[5],n_power[6],k_power[6],index_n,index_k;
	complex index;
	int i,Nb_cauchy = 5;
	
	/* File name */
	snprintf(filename, SIZE_STR_BUFFER*sizeof(char), "%s_index.txt", name);
		
	/* Cauchy method */
	if (!strcmp(method,"Cauchy")){
		lire_tab(filename, "CAUCHY_N", cauchy_n, Nb_cauchy);
		lire_tab(filename, "CAUCHY_K", cauchy_k, Nb_cauchy);
		lire_tab(filename, "N_POWERS", n_power, Nb_cauchy+1);
		lire_tab(filename, "K_POWERS", k_power, Nb_cauchy+1);
SaveDbleTab2file ("stdout", cauchy_n,  Nb_cauchy," ");printf("\n");
SaveDbleTab2file ("stdout", n_power, Nb_cauchy," ");printf("\n");
SaveDbleTab2file ("stdout", cauchy_k,  Nb_cauchy," ");printf("\n");
SaveDbleTab2file ("stdout", k_power, Nb_cauchy," ");printf("\n");
		
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
				fprintf(stderr,"md1D_indice, look-up table incompatible avec lambda = %f pour %s",lambda,name);
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
