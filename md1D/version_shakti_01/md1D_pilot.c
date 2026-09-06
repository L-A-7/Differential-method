/*---------------------------------------------------------------------------------------------*/
/*!	\file		md1D_pilot.c
 *
 * 	\brief		Pilotage de md1D 
 */
/*---------------------------------------------------------------------------------------------*/


#include "md1D_pilot.h"
#define STRSIZE 100
/*---------------------------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------------------------*/
int main(int argc, char *argv[]){

	struct Noms_fichiers nomfichier;
	struct Param_struct param;
	struct Efficacites_struct effic;

	param.verbose =1;

	/* Copie des arguments de la ligne de commande */
	int i;
	char **argvcp; 
	argvcp = (char **) malloc(sizeof(char*)*argc);
	argvcp[0] = (char*) malloc(sizeof(char)*STRSIZE*argc);
	for(i=1;i<=argc-1;i++){
		argvcp[i] = argvcp[i-1] + STRSIZE;
		strncpy(argvcp[i], argv[i],STRSIZE);
	}


	/* Initialisation du programme : lecture des données, allocation de mémoire, etc. */
	md1D_init(&param, &effic, &nomfichier, argc, argvcp);
/*	
(*param.k_2)(param.k2, param.profil, param.N_profil, param.k2_layer, 50.0);
SaveCplxTab2file (param.k2, param.N_profil, "Re", "stdout"); printf("\n");

printf("\nprofil[0] : \n"); SaveDbleTab2file ("stdout",param.profil[0], param.N_profil, "\n"); printf("\n");
printf("\nprofil[1] : \n"); SaveDbleTab2file ("stdout",param.profil[1], param.N_profil, "\n"); printf("\n");
printf("\nprofil[2] : \n"); SaveDbleTab2file ("stdout",param.profil[2], param.N_profil, "\n"); printf("\n");
printf("\nprofil[3] : \n"); SaveDbleTab2file ("stdout",param.profil[3], param.N_profil, "\n"); printf("\n");

printf("%f\n",creal(param.k2_layer[0]));
printf("%f\n",creal(param.k2_layer[1]));
printf("%f\n",creal(param.k2_layer[2]));
return 0;*/

	/* Choix du type de calcul */
	switch (param.type_calcul){
		case STD :
			md1D_standard (&param, &effic, &nomfichier);
			break;
		case VAR_I :
			md1D_variation_incidence_B (&param, &effic, &nomfichier);
			break;
		case VAR_I_BIS :
			md1D_variation_incidence_A (&param, &effic, &nomfichier);
			break;	
/*		case ELLIPSO :
			md1D_ellipso (&param, &effic, &nomfichier);
			break;	
		case VAR_I_ELLIPSO :
			md1D_var_i_ellipso (&param, &effic, &nomfichier);
			break;	
*/	}

printf("Somme des efficacités   : %1.10f\n",effic.somm_eff);
printf("1-Somme des efficacités : %e\n",1-effic.somm_eff);


/*************************************** TESTS *********************************************************/
/*complex *FFT_invk2_directe(double z, complex *TF_invk2, struct Param_struct *par);
complex *FFT_k2_directe(double z, complex *TF_k2, struct Param_struct *par);
complex *FFT_k2_stockee(double z, complex *TF_k2, struct Param_struct *par);

int i,j, Nb=300;
double z;
complex *TF_invk2  = (complex *) malloc(sizeof(complex)*(4*param.N+1));
complex *TF_k2     = (complex *) malloc(sizeof(complex)*(4*param.N+1));
TF_invk2 = FFT_invk2_directe(0.5*param.h, TF_invk2, &param);

SaveCplxTab2file (TF_invk2, 4*param.N+1, "Re", "stdout");*/
/*int N = param.N;
SaveMatrix2file (param.S12, 2*N+1, 2*N+1, "Re", "stdout");
SaveMatrix2file (param.S22, 2*N+1, 2*N+1, "Re", "stdout");
SaveMatrix2file (param.S12, 2*N+1, 2*N+1, "Im", "stdout");
SaveMatrix2file (param.S22, 2*N+1, 2*N+1, "Im", "stdout");*/

/************************************* FIN TESTS *******************************************************/


	/* Libération de la mémoire */
	md1D_free(&param, &effic);
	
	/* Affichage du temps écoulé */
	fprintf(stdout,"Temps écoulé : %f s \n", md1D_chrono(&param));
	
	return 0;		
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_standard (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *	\brief
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_standard (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
{
	/* Calcul de la matrice S de la surface */
	matrice_S(par);

	/* Calcul des amplitudes */
	int j, N=par->N;
	for (j=-N; j<=N; j++) { /* Champ incident */
		par->Ai[j+N] = 0;
	}
	par->Ai[N] = 1;
	md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);

	/* Calcul des efficacités */
	md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);

	/* Ecriture des résultats dans fichier_results */
	md1D_genere_nom_fichier_results(nomfichier->fichier_results, par);
	md1D_ecrire_results(nomfichier->fichier_results, par, eff);

	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *	\brief
 */
/*---------------------------------------------------------------------------------------------*/
#if 0
int md1D_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
{
	complex **S12_TE, **S12_TM, **S22_TE, **S22_TM;

	/* Calculs cas TE */
	par->pola = TE;
	matrice_S(par);
	S12_TE = allocate_CplxMatrix(2*par->N+1, 2*par->N+1);
	S22_TE = allocate_CplxMatrix(2*par->N+1, 2*par->N+1);
	copy_M(S12_TE, S12, 2*par->N+1, 2*par->N+1);
	copy_M(S12_TE, S12, 2*par->N+1, 2*par->N+1);

	/* Calculs cas TM */
	par->pola = TM;
	matrice_S(par);
	S12_TM = S12;
	S22_TM = S22;	
	
	/* Calcul des amplitudes */
	int j, N=par->N;
	for (j=-N; j<=N; j++) { /* Champ incident */
		par->Ai[j+N] = 0;
	}
	par->Ai[N] = 1;
	md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);

	/* Calcul des efficacités */
	md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);

	/* Ecriture des résultats dans fichier_results */
	md1D_genere_nom_fichier_results(nomfichier->fichier_results, par);
	md1D_ecrire_results(nomfichier->fichier_results, par, eff);

	free(S12_TE); free(S22_TE);

	return 0;
}
#endif

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_variation_incidence_A (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier)
 *
 *	\brief	Fait varier l'incidence en utilisant une seule matrice S et plusieurs vecteurs chp incident \n
 *			Ai. Adapté aux grands N (ex : calculs de diffusion)
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_variation_incidence_A (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier) 
{
	int n,j;
	FILE *fp;
	int N = par->N;
	
	/* Calcul de la matrice S de la surface */
	matrice_S(par);


	/* Appel de md1D_efficacités pour calcul de Nmax_super et Nmin_super */
	md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);

	for (n=eff->Nmin_super; n<=eff->Nmax_super; n++) {
		for (j=-N; j<=N; j++) {
			par->Ai[j+N] = 0;
		}
		par->Ai[n+N] = 1;
		
		md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);
		md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);
		
		/* Récupération de la grandeur */
		par->var_i[n-eff->Nmin_super] = eff->eff_R[n-eff->Nmin_super]; /* Efficacité faisceau réfléchi */
		
	}

	/* Ecriture des résultats dans fichier_results */
	md1D_genere_nom_fichier_results(nomfichier->fichier_results, par);
	md1D_ecrire_results(nomfichier->fichier_results, par, eff);

	/* Ajout d'une ligne contenant var_i */
	if (!(fp = fopen(nomfichier->fichier_results,"a"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nomfichier->fichier_results);
		exit(EXIT_FAILURE);
	}
	int LMAX = 100; /* NORMALEMENT UNE MACRO */
	fprintf(fp,"\nvar_i = "); ecrire_dble_tab(fp, par->var_i, eff->Nmax_super-eff->Nmin_super+1, " ", LMAX,"\n");
	fclose(fp);


	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_variation_incidence_B (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier) 
 *
 *	\brief	Fait varier l'incidence en calculant une matrice S ppour chaque angle_i. \n
 *			Adapté aux petits N (ex : réseaux diélectriques de faibles pas)
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_variation_incidence_B (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier) 
{
	int j,k;
	FILE *fp;
	int N = par->N;
	double i_min = 0;
	double i_max = 89;
	double delta_i = 1;
	int N_i = ROUND((i_max - i_min)/delta_i);
	double *var_i_angle, *var_i_effR, *var_i_effR_p1, *var_i_effR_m1, *var_i_effR_m2, *var_i_modR, *var_i_argR;

	var_i_angle = (double *) malloc(sizeof(double)*(N_i+1));
	var_i_effR  = (double *) malloc(sizeof(double)*(N_i+1));
	var_i_effR_p1  = (double *) malloc(sizeof(double)*(N_i+1));
	var_i_effR_m1  = (double *) malloc(sizeof(double)*(N_i+1));
	var_i_effR_m2  = (double *) malloc(sizeof(double)*(N_i+1));
	var_i_modR  = (double *) malloc(sizeof(double)*(N_i+1));
	var_i_argR  = (double *) malloc(sizeof(double)*(N_i+1));


	/* Boucle sur l'angle d'incidence */
	for (j=0; j<=N_i; j++){
	
		par->angle_i = (i_min + delta_i*(double)j)*PI/180.0;
		par->sigma0 = par->k0*sin(par->angle_i);
		var_i_angle[j] = (i_min + delta_i*(double)j);
			
		fprintf(stdout,"\r i = %3.0f   ",par->angle_i*180.0/PI);fflush(stdout);
		par->verbose = 0;
	
		/* Calcul de la matrice S de la surface */
		matrice_S(par);

		/* Calcul des amplitudes */
		for (k=-N; k<=N; k++) {
			par->Ai[k+N] = 0;
		}
		par->Ai[N] = 1;
		md1D_amplitudes(par->Ai, par->A0, par->Ah, par->S12, par->S22, par);

		/* Calcul des efficacités */
		md1D_efficacites(par->Ai, par->A0, par->Ah, par, eff);
			
		/* Récupération de la grandeur */
		var_i_effR[j] = eff->eff_R[-eff->Nmin_super];       /* Efficacité faisceau réfléchi */
		var_i_effR_p1[j] = eff->eff_R[-eff->Nmin_super+1];       /* Efficacité ordre 1 */
		var_i_effR_m1[j] = eff->eff_R[-eff->Nmin_super-1];       /* Efficacité ordre -1 */
		var_i_effR_m2[j] = eff->eff_R[-eff->Nmin_super-2];       /* Efficacité ordre -1 */
		var_i_modR[j] = cabs(par->S12[par->N][par->N]);     /* module du facteur de réflexion du spéculaire */
		var_i_argR[j] = carg(par->S12[par->N][par->N]); /* argument du facteur de réflexion complexe du spéculaire */

	}
	
	/* Ecriture des résultats dans fichier_results */
	md1D_genere_nom_fichier_results(nomfichier->fichier_results, par);
	md1D_ecrire_results(nomfichier->fichier_results, par, eff);

	/* Ajout d'une ligne contenant var_i */
	if (!(fp = fopen(nomfichier->fichier_results,"a"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nomfichier->fichier_results);
		exit(EXIT_FAILURE);
	}
	int LMAX = 1000; /* NORMALEMENT UNE MACRO */
	
	fprintf(fp,"\nvar_i_angle = "); ecrire_dble_tab(fp, var_i_angle, N_i+1, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_effR  = "); ecrire_dble_tab(fp, var_i_effR, N_i+1, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_effR_p1  = "); ecrire_dble_tab(fp, var_i_effR_p1, N_i+1, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_effR_m1  = "); ecrire_dble_tab(fp, var_i_effR_m1, N_i+1, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_effR_m2  = "); ecrire_dble_tab(fp, var_i_effR_m2, N_i+1, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_modR  = "); ecrire_dble_tab(fp, var_i_modR, N_i+1, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_argR  = "); ecrire_dble_tab(fp, var_i_argR, N_i+1, " ", LMAX,"\n");
		
	fclose(fp);


	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_init (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier) 
 *
 *	\brief
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_init (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier, 
				int argc, char **argvcp) 
{

	
	/* Nom du fichier param par defaut */
	sprintf(nomfichier->fichier_param,"md1D_param.txt");

	/* Lecture de fichier param  et fichier profil en argument de la ligne de commande */
	lire_str_arg(nomfichier->fichier_param, "-param", argc, argvcp);
	if (lire_str_arg(nomfichier->fichier_profil, "-fichier_profil", argc, argvcp) != 0){
		sprintf(nomfichier->fichier_profil, NON_LU);}

	/* Lecture des paramètres par defaut dans fichier_param */
	md1D_lire_param(nomfichier->fichier_param, nomfichier->fichier_profil, par);

	/* Lecture des paramètres en arguments de la ligne de commande */
	md1D_lire_args(nomfichier->fichier_profil, par, argc, argvcp);

	/* Affichage des valeurs lues */
	md1D_affiche_valeurs_lues(par, nomfichier);

	/* Initialisations de certaines variables */
	md1D_variables_init(par, eff);
	
	/* Allocation de mémoire pour les tableaux */
	md1D_alloc(par, eff);

	/* Lecture du profil h(x) décrivant la surface */
	(*par->md1D_lire_profil)(nomfichier->fichier_profil, par);

	return 0;
}



/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_variables_init(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *	\brief
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_variables_init(struct Param_struct *par, struct Efficacites_struct *eff)
{
	par->clock0 = clock(); /* Initialisation des chronomètres                  */
	time(&(par->time0));   /* clock0 (courte durées) et time0 (longues durées) */
	par->last_clock = clock();
	time(&(par->last_time)); 


	par->k0 = 2*PI*par->n_super/par->lambda;
	par->kh = 2*PI*par->n_sub/par->lambda;
	par->delta_sigma = 2*PI/par->L;
	par->sigma0 = par->k0*sin(par->angle_i);
	
	par->Nstep_S = ROUND(ceil((par->h/par->NS)/par->delta_h_approx));
	par->delta_h = (par->h/par->NS)/par->Nstep_S;
	par->Nstep = ROUND(par->h/par->delta_h);

	/* Alignement des pointeurs de fonction */
	par->matrice_T = (par->pola == TE ? matrice_T_TE : matrice_T_TM);
	switch (par->type_profil) {
		case H_X          : 
			par->md1D_lire_profil = md1D_lire_profil_H_X;
			par->k_2 = k2_H_X;
			par->invk_2 = invk2_H_X;
			break;
		case MULTICOUCHES : 
			par->md1D_lire_profil = md1D_lire_profil_MULTI;
			par->k_2 = k2_MULTI;
			par->invk_2 = invk2_MULTI;
			break;
	/*	case N_XYZ        : par->md1D_lire_profil = md1D_lire_profil_N_XYZ; */
	}

	par->verbose = 1; /* Si VERBOSE = 1 : affiche plus d'infos sur le terminal */
	
	par->STOCKER_TF = 1; /* Optimisation : On stocke ou non les TFs, + rapide, mais demande + de mémoire */
	par->BLOC_TAILLE_z2n = (int) 10*par->Nstep; /* PAS OPTIMISÉ, à REVOIR */
	par->TAILLE_z2n = par->BLOC_TAILLE_z2n;
	par->N_z2n = 0;

	return 0;
}



/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_alloc(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *	\brief
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_alloc(struct Param_struct *par, struct Efficacites_struct *eff)
{

	par->profil = allocate_DbleMatrix(par->N_layers+3,par->N_profil);
/*	par->indice      = (complex *) malloc(sizeof(complex)*(par->N_layers+2));*/
	par->k2_layer    = (complex *) malloc(sizeof(complex)*(par->N_layers+2));
	par->invk2_layer = (complex *) malloc(sizeof(complex)*(par->N_layers+2));
	par->k2        = (complex *) malloc(sizeof(complex)*par->N_profil);
	par->invk2     = (complex *) malloc(sizeof(complex)*par->N_profil);
	par->TF_k2     = (complex *) malloc(sizeof(complex)*(4*par->N+1));
	par->TF_invk2  = (complex *) malloc(sizeof(complex)*(4*par->N+1));

	par->S12 = allocate_CplxMatrix(2*par->N+1,2*par->N+1);
	par->S22 = allocate_CplxMatrix(2*par->N+1,2*par->N+1);

	par->Ai = (complex *) malloc(sizeof(complex)*(2*par->N+1));
	par->A0 = (complex *) malloc(sizeof(complex)*(2*par->N+1));
	par->Ah = (complex *) malloc(sizeof(complex)*(2*par->N+1));	

	eff->eff_R =       (double *) malloc(sizeof(double)*(2*par->N+1));
	eff->eff_T =       (double *) malloc(sizeof(double)*(2*par->N+1));
	eff->N_eff_R =     (double *) malloc(sizeof(double)*(2*par->N+1));
	eff->N_eff_T =     (double *) malloc(sizeof(double)*(2*par->N+1));
	eff->theta_eff_R = (double *) malloc(sizeof(double)*(2*par->N+1));
	eff->theta_eff_T = (double *) malloc(sizeof(double)*(2*par->N+1));

	par->var_i =       (double *) malloc(sizeof(double)*1000); /* A REVOIR */
	par->var_i2 =       (double *) malloc(sizeof(double)*1000); /* A REVOIR */

	if (par->STOCKER_TF) {
		par->z2n[0] = (long int *) malloc(sizeof(long int)*(par->BLOC_TAILLE_z2n));
		par->z2n[1] = (long int *) malloc(sizeof(long int)*(par->BLOC_TAILLE_z2n));
		par->tab_TF_k2 = allocate_CplxMatrix(par->BLOC_TAILLE_z2n, 4*par->N+1);
		if (par->pola == TM || par->type_calcul == ELLIPSO){
			par->tab_TF_invk2 = allocate_CplxMatrix(par->BLOC_TAILLE_z2n, 4*par->N+1);
		}
	}
	
	return 0;
}



/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_free(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *	\brief
 */
/*---------------------------------------------------------------------------------------------*/
int md1D_free(struct Param_struct *par, struct Efficacites_struct *eff)
{

	free(par->profil);
	free(par->TF_k2);
	free(par->S12[0]);
	free(par->S12);
	free(par->S22[0]);
	free(par->S22);
	free(par->A0);
	free(par->Ai);
	free(par->Ah);
	free(eff->eff_R);
	free(eff->eff_T);
	free(eff->N_eff_R);
	free(eff->N_eff_T);
	free(eff->theta_eff_R);
	free(eff->theta_eff_T);
	
	return 0;
}


int md1D_lire_args(char *fichier_profil, struct Param_struct *par, int argc, char **argvcp)
{

	lire_dble_arg(&par->L, "-L", argc, argvcp);
	lire_dble_arg(&par->h, "-h", argc, argvcp);
/*	lire_int_arg(&par->pola, "-pola", argc, argvcp);*/
	lire_dble_arg(&par->lambda, "-lambda", argc, argvcp);
	lire_dble_arg(&par->angle_i, "-angle_i", argc, argvcp);
	lire_dble_arg(&par->delta_h_approx, "-delta_h_approx", argc, argvcp);
	lire_int_arg(&par->N, "-N", argc, argvcp);
	lire_int_arg(&par->NS, "-NS", argc, argvcp);
	lire_str_arg(par->nom_profil, "-nom_profil", argc, argvcp);

	return 0;	
}


/* METTRE CE QUI SUIT DANS md1D_in_out ou utils */

/* Lit la valeur de l'argument de la ligne de commande indiqué sous la forme "-label valeur" */
int lire_str_arg(char *dest, char *label, int argc, char **argvcp)
{
	int i;

	for(i=1;i<=argc-2;i++){
		if (!strcmp(argvcp[i],label)){
			strncpy(dest, argvcp[i+1], STRSIZE);
			return 0; 
		}
	}
	return 1;
}

int lire_dble_arg(double *res, char *label, int argc, char **argvcp)
{
	int i;
	char strtmp[STRSIZE], *endptr;
	double tmp;

	for(i=1;i<=argc-2;i++){
		if (!strcmp(argvcp[i],label)){
			strncpy(strtmp, argvcp[i+1], STRSIZE);
			tmp = strtod(strtmp, &endptr);
			if (strtmp != endptr){
				*res = tmp;
				return 0; 
			}
		}
	}
	return 1;
}

int lire_int_arg(int *res, char *label, int argc, char **argvcp)
{
	int i, tmp;
	char strtmp[STRSIZE], *endptr;

	for(i=1;i<=argc-2;i++){
		if (!strcmp(argvcp[i],label)){
			strncpy(strtmp, argvcp[i+1], STRSIZE);
			tmp = (int) strtod(strtmp, &endptr);
			if (strtmp != endptr){
				*res = tmp;
				return 0; 
			}
		}
	}
	return 1;
}
