/*---------------------------------------------------------------------------------------------*/
/*!	\file		md2D_pilot.c
 *
 * 	\brief		Pilotage de md2D 
 */
/*---------------------------------------------------------------------------------------------*/
#include<stdio.h>
#include<sched.h>
#include<stdlib.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>

#include "md2D_pilot.h"
/*---------------------------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------------------------*/
int main(int argc, char *argv[]){

	struct Noms_fichiers nomfichier;
	struct Param_struct param;
	struct Efficacites_struct effic;

#if SCHED_FIFO
        /* Set scheduler policy */
        struct sched_param sp;

        /* Get max priority for this processus */
        sp.sched_priority = sched_get_priority_max(SCHED_FIFO);
        sched_setscheduler(getpid(), SCHED_FIFO, &sp);
#endif

	/*param.verbosity =1;*/
	
	/* Copie des arguments de la ligne de commande */
	int i;
	char **argvcp; 
	argvcp = (char **) malloc(sizeof(char*)*argc);
	argvcp[0] = (char*) malloc(sizeof(char)*SIZE_STR_BUFFER*argc);
	for(i=1;i<=argc-1;i++){
		argvcp[i] = argvcp[i-1] + SIZE_STR_BUFFER;
		strncpy(argvcp[i], argv[i],SIZE_STR_BUFFER);
	}
	param.argc = argc;
	param.argvcp = argvcp;

	/* Initialisation du programme : lecture des données, allocation de mémoire, etc. */
	md2D_init(&param, &effic, &nomfichier);

/*********************************** DEBUG *************************************************************/

/*******************************************************************************************************/

	/* Choix du type de calcul */
	if (!strcmp(param.calcul_type,"STD")){
		md2D_conical_FFF (&param, &effic, &nomfichier);
	}else if (!strcmp(param.calcul_type,"ELLIPSO")){
		md2D_conical_FFF_ellipso (&param, &effic, &nomfichier);
	}else if (!strcmp(param.calcul_type,"VAR_LAMBDA_ELLIPSO")){
		md2D_var_lambda_ellipso (&param, &effic, &nomfichier);
	}else{
		fprintf(stderr, "%s, line %d : ERROR, unknown calculation type (\"%s\")\n",__FILE__,__LINE__,param.calcul_type);
		exit(EXIT_FAILURE);
	}


/*************************************** TESTS *********************************************************/


/************************************* FIN TESTS *******************************************************/


	/* Freeing memory */
	md2D_free(&param, &effic);
	
	/* Showing used time & efficiencies sum */
	if (param.verbosity >= 1){
		printf("Efficiencies summ   : %1.10f\n",effic.sum_eff);
		printf("1-Efficiencies summ : %e\n",1-effic.sum_eff);
		fprintf(stdout,"Temps écoulé : %f s \n", md2D_chrono(&param));
	}
	
	return 0;		
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		md2D_conical_FFF (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *		\brief	Complex field calculation, in the case of a 1D structure with conical incidence.
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_conical_FFF (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
{

	/* Calcul des limites des modes propagatifs */
	md2D_propagativ_limits(par, eff);
		
	/* Amplitude du champ incident */
	md2D_incident_field(par,eff);

	S_matrix(par);

	/* Calcul des amplitudes */
	md2D_amplitudes(par->Ai, par->Ar, par->At, par->S12, par->S22, par);

	/* Calcul des efficacités */
	md2D_efficiencies(par->Ai, par->Ar, par->At, par, eff);

	/* Ecriture des résultats dans fichier_results */
	md2D_genere_nom_fichier_results(nomfichier->fichier_results, par);

	md2D_ecrire_results(nomfichier->fichier_results, par, eff);


	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		md2D_conical_FFF_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *		\brief	Calcule efficacité TE, TM et déphasage (la valeur initiale de psi n'a pas d'influence).
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_conical_FFF_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
{
	FILE *fp;
	int n;
	int vec_size = par->vec_size;
	int mid = par->vec_middle;
	REAL complex *Ai_TE, *Ar_TE, *At_TE;
	double *eff_r_TE, *eff_t_TE, *delta_r, *delta_t;
	REAL complex rp, rs, tp, ts;
	double spec_tan_Psi, spec_cos_delta;

	Ai_TE = (REAL complex *) malloc(sizeof(REAL complex)*2*vec_size);
	Ar_TE = (REAL complex *) malloc(sizeof(REAL complex)*2*vec_size);
	At_TE = (REAL complex *) malloc(sizeof(REAL complex)*2*vec_size);
	eff_r_TE = (double *) malloc(sizeof(double)*vec_size);
	eff_t_TE = (double *) malloc(sizeof(double)*vec_size);
	delta_r = (double *) malloc(sizeof(double)*vec_size);
	delta_t = (double *) malloc(sizeof(double)*vec_size);
		
	/* Calcul des limites des modes propagatifs */
	md2D_propagativ_limits(par, eff);
	int Nmin_super = eff->Nmin_super;
	int Nmax_super = eff->Nmax_super;
	int Nmin_sub = eff->Nmin_sub;
	int Nmax_sub = eff->Nmax_sub;
		
	/**/	
	S_matrix(par);
		
	/* cas TE, psi = 0° */
	par->psi = 0.0;
	md2D_incident_field(par,eff);
	md2D_amplitudes(par->Ai, par->Ar, par->At, par->S12, par->S22, par);
	md2D_efficiencies(par->Ai, par->Ar, par->At, par, eff);
	CopyCplxTab(Ai_TE, par->Ai, 2*vec_size);
	CopyCplxTab(Ar_TE, par->Ar, 2*vec_size);
	CopyCplxTab(At_TE, par->At, 2*vec_size);
	CopyDbleTab(eff_r_TE, eff->eff_r, Nmax_super-Nmin_super+1);
	CopyDbleTab(eff_t_TE, eff->eff_t, Nmax_sub-Nmin_sub+1);
	
	/* cas TM, psi = 90° = pi/2 */
	par->psi = PI/2.0;
	md2D_incident_field(par,eff);
	md2D_amplitudes(par->Ai, par->Ar, par->At, par->S12, par->S22, par);
	md2D_efficiencies(par->Ai, par->Ar, par->At, par, eff);

	/* Dephasage */
	for (n=0;n<=vec_size-1;n++){
/*		delta_r_s[n] = carg(par->Ar_TE[n])*180.0/PI; */			/* Déphasage du champ E ...*/
/*		delta_r_p[n] = carg(par->Ar[n+vec_size])*180.0/PI; */ /* !!! Déphasage du champ H !!!*/
		rp = par->Ar[n+vec_size]/par->Ai[n+vec_size];
		rs = Ar_TE[n]/Ai_TE[n];
		tp = par->At[n+vec_size]/par->Ai[n+vec_size];
		ts = At_TE[n]/Ai_TE[n];
		delta_r[n] = carg(rp*conj(rs))*180.0/PI;
		delta_t[n] = carg(tp*conj(ts))*180.0/PI;
	}

	rp = par->Ar[mid+vec_size]/par->Ai[mid+vec_size];
	rs = Ar_TE[mid]/Ai_TE[mid];
	spec_tan_Psi =	  cabs(rp/rs);
	spec_cos_delta = -creal(rp/rs)/cabs(rp/rs);
	
	/* Ecriture des résultats dans fichier_results */
	md2D_genere_nom_fichier_results(nomfichier->fichier_results, par);
	md2D_ecrire_results(nomfichier->fichier_results, par, eff);

	/* Writing values TE, TM, delta to the results file */
	if (!(fp = fopen(nomfichier->fichier_results,"a"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nomfichier->fichier_results);
		exit(EXIT_FAILURE);
	}
	int LMAX = 1000; /* NORMALEMENT UNE MACRO */
	fprintf(fp,"\neff_r_TE   = "); ecrire_dble_tab(fp, eff_r_TE, Nmax_super-Nmin_super+1, " ", LMAX,"\n");
	fprintf(fp,"\neff_r_TM   = "); ecrire_dble_tab(fp, eff->eff_r, Nmax_super-Nmin_super+1, " ", LMAX,"\n");
	fprintf(fp,"\ndelta_r    = "); ecrire_dble_tab(fp, &delta_r[mid+Nmin_super], Nmax_super-Nmin_super+1, " ", LMAX,"\n");
	fprintf(fp,"\neff_t_TE   = "); ecrire_dble_tab(fp, eff_t_TE, Nmax_sub-Nmin_sub+1, " ", LMAX,"\n");
	fprintf(fp,"\neff_t_TM   = "); ecrire_dble_tab(fp, eff->eff_t, Nmax_sub-Nmin_sub+1, " ", LMAX,"\n");
	fprintf(fp,"\ndelta_t    = "); ecrire_dble_tab(fp, &delta_t[mid+Nmin_sub], Nmax_sub-Nmin_sub+1, " ", LMAX,"\n");
	fprintf(fp,"\nspec_tan_Psi = % 1.12e",spec_tan_Psi);
	fprintf(fp,"\nspec_cos_delta = % 1.12e",spec_cos_delta);
	fclose(fp);

	free(Ar_TE);
	free(Ai_TE);
	free(At_TE);
	free(eff_r_TE);
	free(eff_t_TE);
	free(delta_r);
	free(delta_t);

	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_var_lambda_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *	\brief	Calcul des amplitudes, efficacités, déphasages en TE et TM, et du dephasage polarimetrique
 *
 *	\todo	A améliorer, faire plus clean !
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_var_lambda_ellipso (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier)
{
	int i;
	FILE *fp;
	int N_layers = par->N_layers;
	char material[10][SIZE_STR_BUFFER];
	double lambda_min = 210;
	double lambda_max = 770;
	double delta_lambda = 10;
	par->Ni = ROUND((lambda_max - lambda_min)/delta_lambda)+1;
	int Ni = par->Ni;
	double *var_lambda, *var_tan_Psi, *var_cos_delta, *var_i_effR_s, *var_i_effR_p, *var_i_delta_s, *var_i_delta_p, *var_i_delta;
	REAL complex **S12 = par->S12;
	REAL complex index;
	int vec_size = par->vec_size;

	if (par->type_profil != H_X && par->type_profil != MULTICOUCHES) {
		fprintf(stderr,"type de profil incompatible avec fonction var_lambda");
		return -1;}
		
	sprintf(material[0],"AIR");
	sprintf(material[1],"POLY03");
	sprintf(material[2],"OXIDE_THERM");
	sprintf(material[3],"SI_CRISTAL");
	
	var_lambda    = (double *) malloc(sizeof(double)*Ni);
	var_tan_Psi   = (double *) malloc(sizeof(double)*Ni);
	var_cos_delta = (double *) malloc(sizeof(double)*Ni);
	var_i_effR_s  = (double *) malloc(sizeof(double)*Ni);
	var_i_effR_p  = (double *) malloc(sizeof(double)*Ni);
	var_i_delta_s = (double *) malloc(sizeof(double)*Ni);
	var_i_delta_p = (double *) malloc(sizeof(double)*Ni);
	var_i_delta   = (double *) malloc(sizeof(double)*Ni);

		
	/* Amplitude du champ incident */
	md2D_incident_field(par,eff);

	/* Boucle sur lambda */
	for (par->ni=0; par->ni<=Ni-1; (par->ni)++){
	
		/* Reinitializing all values dependent on lambda */
		par->lambda = lambda_min + delta_lambda*(double)par->ni;
		par->n_sub = md2D_index("SI_CRISTAL",par->lambda,"Lookup");
		var_lambda[par->ni] = par->lambda;
		par->k_super = 2*PI*par->n_super/par->lambda;
		par->k_sub = 2*PI*par->n_sub/par->lambda;	
		par->Delta_sigma = 2*PI/par->L;
		par->sigma0 = par->k_super*sin(par->theta_i)*cos(par->phi_i);
		par->ky_0   = par->k_super*sin(par->theta_i)*sin(par->phi_i);
		md2D_arrays_init(par, eff);
			
		/* Détermination des indices pour le lambda considéré */
		for (i=0; i<=N_layers+1; i++){
			index = md2D_index(material[i],par->lambda,"Lookup");
			par->k2_layer[i]    = (index*par->k_super)*(index*par->k_super); /* Asserts n_super = 1 */ 
			par->invk2_layer[i] = 1/par->k2_layer[i]; 
		}
		
/*		fprintf(stdout,"lambda = %3.0f\n",par->lambda);fflush(stdout);*/
		par->verbosity = 0;

		/* Calcul des limites des modes propagatifs */
		md2D_propagativ_limits(par, eff);
		
		/**/	
		S_matrix(par);
		
		/* cas TE, psi = 0° */
		par->psi = 0.0;
		md2D_incident_field(par,eff);
		md2D_amplitudes(par->Ai, par->Ar, par->At, par->S12, par->S22, par);
		md2D_efficiencies(par->Ai, par->Ar, par->At, par, eff);
		/*var_i_Ar_s[par->ni] = par->Ar[par->N];*/ /* Oth order reflected field */
		/*var_i_Ai_s[par->ni] = par->Ai[par->N];*/ /* incident field */
		var_i_effR_s[par->ni] = eff->eff_r[-eff->Nmin_super];       /* Efficacité faisceau réfléchi */
		var_i_delta_s[par->ni] = carg(par->S12[par->vec_middle][par->vec_middle])*180.0/PI;       

		/* cas TM, psi = 90° = pi/2 */
		par->psi = PI/2.0;
		md2D_incident_field(par,eff);
		md2D_amplitudes(par->Ai, par->Ar, par->At, par->S12, par->S22, par);
		md2D_efficiencies(par->Ai, par->Ar, par->At, par, eff);
		/*var_i_Ar_p[par->ni] = par->Ar[par->N+vec_size];*/ /* Oth order reflected field */
		/*var_i_Ai_p[par->ni] = par->Ai[par->N+vec_size];*/ /* incident field */
		var_i_effR_p[par->ni] = eff->eff_r[-eff->Nmin_super];       /* Efficacité faisceau réfléchi */
		var_i_delta_p[par->ni] = carg(par->S12[par->vec_middle+vec_size][par->vec_middle+vec_size])*180.0/PI;       

		/* Récupération des grandeurs */
		var_i_delta[par->ni] = carg(-S12[par->N][par->N]*conj(S12[par->N+vec_size][par->N+vec_size]))*180.0/PI;       
		var_tan_Psi[par->ni] = cabs( S12[par->N+vec_size][par->N+vec_size] / S12[par->N][par->N] );
		var_cos_delta[par->ni] = -creal( S12[par->N+vec_size][par->N+vec_size] / S12[par->N][par->N] ) / var_tan_Psi[par->ni];
	}

	/* Ecriture des résultats dans fichier_results */
	md2D_genere_nom_fichier_results(nomfichier->fichier_results, par);
	md2D_ecrire_results(nomfichier->fichier_results, par, eff);

	/* Ajout d'une ligne contenant var_i */
	if (!(fp = fopen(nomfichier->fichier_results,"a"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,nomfichier->fichier_results);
		exit(EXIT_FAILURE);
	}
	int LMAX = 1000; /* NORMALEMENT UNE MACRO */
	
	fprintf(fp,"\nvar_lambda    = "); ecrire_dble_tab(fp, var_lambda,    Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_tan_Psi   = "); ecrire_dble_tab(fp, var_tan_Psi,   Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_cos_delta = "); ecrire_dble_tab(fp, var_cos_delta, Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_effR_s  = "); ecrire_dble_tab(fp, var_i_effR_s,  Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_effR_p  = "); ecrire_dble_tab(fp, var_i_effR_p,  Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_delta_s = "); ecrire_dble_tab(fp, var_i_delta_s, Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_delta_p = "); ecrire_dble_tab(fp, var_i_delta_p, Ni, " ", LMAX,"\n");
	fprintf(fp,"\nvar_i_delta   = "); ecrire_dble_tab(fp, var_i_delta,   Ni, " ", LMAX,"\n");
	fprintf(fp,"\nNb_lambda = %d",par->Ni);
		
	fclose(fp);

	free(var_lambda);
	free(var_tan_Psi);
	free(var_cos_delta);
	free(var_i_effR_s);
	free(var_i_effR_p);
	free(var_i_delta_s);
	free(var_i_delta_p);
	free(var_i_delta);

	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_init (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier) 
 *
 *	\brief	Initialisation du programme
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_init (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier) 
{

	/* Lecture des paramètres par defaut dans fichier_param */
	md2D_lire_param(nomfichier, par);

	/* Allocation de mémoire pour le profil */
	md2D_alloc_init_profil(par);

	/* Lecture du profil h(x) décrivant la surface */
	(*par->md2D_lire_profil)(nomfichier->fichier_profil, par);

	/* Initialisations de certaines variables */
	md2D_variables_init(par, eff);

	/* Allocation de mémoire pour les tableaux */
	md2D_alloc(par, eff);

	/* Initialising some arrays */
	md2D_arrays_init(par, eff);
	
	/* Affichage des paramètres lus et calculés */
	md2D_affiche_valeurs_param(par, nomfichier);

	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn	int md2D_alloc_init_profil(struct Param_struct *par)
 *
 *	\brief	Memory allocation for profil
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_alloc_init_profil(struct Param_struct *par)
{
	/* Initilisation de variables */
	par->k_super = 2*PI*par->n_super/par->lambda;
	par->k_sub = 2*PI*par->n_sub/par->lambda;
	par->Delta_sigma = 2*PI/par->L;
	par->sigma0 = par->k_super*sin(par->theta_i)*cos(par->phi_i);
	par->ky_0   = par->k_super*sin(par->theta_i)*sin(par->phi_i);
			
	/* Allocations */	
	if (par->type_profil == N_XYZ) {
		par->n_xyz = allocate_CplxMatrix(par->N_z,par->N_x);
	}else{
		par->profil = allocate_DbleMatrix(par->N_layers+3,par->N_x);
		par->k2_layer    = (REAL complex *) malloc(sizeof(REAL complex)*(par->N_layers+2));
		par->invk2_layer = (REAL complex *) malloc(sizeof(REAL complex)*(par->N_layers+2));
	}


	/* Alignement des pointeurs de fonction */
/*	par->matrice_T = (par->pola == TE ? matrice_T_TE : matrice_T_TM);
*/	switch (par->type_profil) {
		case H_X          : 
			par->md2D_lire_profil = md2D_lire_profil_H_X;
			par->k_2 = k2_H_X;
			par->invk_2 = invk2_H_X;
			par->Normal_function = Normal_H_X;
			break;
		case MULTICOUCHES : 
			par->md2D_lire_profil = md2D_lire_profil_MULTI;
			par->k_2 = k2_MULTI;
			par->invk_2 = invk2_MULTI;
			break;
		case N_XYZ        : 
			par->md2D_lire_profil = md2D_lire_profil_N_XYZ;
			par->k_2 = k2_N_XYZ;
			par->invk_2 = invk2_N_XYZ;
			break;
	}
	
	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_variables_init(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *	\brief	Initialisations des variables
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_variables_init(struct Param_struct *par, struct Efficacites_struct *eff)
{
	double K,sig0;
	par->clock0 = clock(); /* Initialisation des chronomètres                  */
	time(&(par->time0));   /* clock0 (courte durées) et time0 (longues durées) */
	par->last_clock = clock();
	time(&(par->last_time)); 

	/* Initialisation des variables en mode "AUTO" */
	/* delta_h = lambda x (n_re + n_im) / 1000 */
	if (par->delta_h == AUTO) {
		par->delta_h = par->lambda/(cabs(par->n_sub)*1000);
	}
	/* NS = (h/lambda) x (n_re + n_im) x 5 */
	if (par->NS == AUTO) {
		par->NS = CEIL( (par->h / par->lambda)*cabs(par->n_sub)*5 );
	}
	/* N : on ajoute 10% de modes evanescents (3 au minimum) */
	if (par->N == AUTO) {
		K = 2.0*PI/par->L;
		fprintf(stderr, "WARNING, Automatic N determination");
		if (par->phi_i != 0){fprintf(stderr, "WARNING phi_i = %f deg, not taken into account for automatic N determination\n",par->phi_i*180/PI);}
		if (!strcmp(par->calcul_type,"STD") || !strcmp(par->calcul_type,"ELLIPSO") ){
			sig0 = par->sigma0;
		}else if (!strcmp(par->calcul_type,"VAR_I") || !strcmp(par->calcul_type,"VAR_I_BIS") || \
			!strcmp(par->calcul_type,"VAR_I_ELLIPSO")){
			sig0 = par->k_super; /* correspond à sigma0 pour theta_i = 90° => N constant et suffisant de 0 à 90°*/
		}else{
			fprintf(stderr, "%s ligne %d, Calcul AUTO de N : %s, type calcul inconnu\n",__FILE__, __LINE__,par->calcul_type);
			exit(EXIT_FAILURE);
		}
		int N_limit =  FLOOR(( MAX(creal(par->k_super),creal(par->k_sub)) + fabs(sig0))/par->Delta_sigma);
		/* On ajoute 10% de modes evanescents (3 au minimum) */
		int N_evanesc = ROUND(MAX(3,0.1*N_limit));
		par->N = N_limit + N_evanesc; 
	}

	/* vec_size : size of most matrices and vetcors */	
	par->vec_size = 2*par->N+1; /* For classical 2D case (not for 3D) */

	/* vec_middle, ex.: sigma[vec_middle] = sigma0 ...*/
	par->vec_middle = par->N;
	/* for H_X profile the normal components need to be calculated only once */
	par->HX_Normal_CALCULATED = 0;
		
	/* ni et Ni, pour compteurs en angle_i */
	par->ni = 0;
	par->Ni = 1; /* Pour compatibilité (une autre valeur sera affectée par les fonctions var_i_... )*/
	
	/* Calcul de la valeur exacte de delta_h de sorte qu'il y en ait un nb entier à chaque étape Matrice-S*/
	par->Nstep_S = ROUND(ceil((par->h/par->NS)/par->delta_h));
	par->delta_h = (par->h/par->NS)/par->Nstep_S;
	par->Nstep = ROUND(par->h/par->delta_h);

	/* tab_NS_ENABLED */
	if (par->READ_tab_NS || par->imposed_S_steps){
		par->tab_NS_ENABLED = 1;
	}else{
		par->tab_NS_ENABLED = 0;
	}
	
	/* imposed_S_steps, NS estimation for malloc */
	if (par->imposed_S_steps){
		par->NS = CEIL(10*par->h/par->lambda)+par->N_imposed_S_steps;
	}
	
	/*par->verbosity = 2;*/ /* Si verbosity > 0 : affiche plus d'infos sur le terminal */
	
	/* functions pointers alignements */
	if (!strcmp(par->calcul_method,"DM")){
			par->T_Matrix = dm_T_Matrix;
	}else if (!strcmp(par->calcul_method,"RCWA")){
			par->T_Matrix = rcwa_T_Matrix;
	}else{
		fprintf(stderr, "%s, line %d : ERROR, unknown calculation method (\"%s\")\n",__FILE__,__LINE__,par->calcul_method);
		exit(EXIT_FAILURE);
	}
	
	par->STOCKER_TF = 0; /* Optimisation : On stocke ou non les TFs, + rapide, mais demande + de mémoire */
	par->BLOC_TAILLE_z2n = (int) 10*par->Nstep; /* PAS OPTIMISÉ, à REVOIR */
	par->TAILLE_z2n = par->BLOC_TAILLE_z2n;
	par->N_z2n = 0;

	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_arrays_init(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *		\brief	Some arrays initialisations
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_arrays_init(struct Param_struct *par, struct Efficacites_struct *eff)
{
	int j;
	
	REAL complex k_super2 = par->k_super*par->k_super;
	REAL complex k_sub2 = par->k_sub*par->k_sub;
	REAL complex ky_02 = par->ky_0*par->ky_0;
	
	for(j=0;j<=par->vec_size-1;j++){
		par->sigma[j] = (j-par->vec_middle)*par->Delta_sigma + par->sigma0;
		par->kz_super[j] = csqrt(k_super2 - par->sigma[j]*par->sigma[j] - ky_02);
		par->kz_sub[j]   = csqrt(k_sub2   - par->sigma[j]*par->sigma[j] - ky_02);
	}
	md2D_PsiMatrix(par->Psi_super, par->k_super, par->kz_super, par);
	md2D_PsiMatrix(par->Psi_sub, par->k_sub, par->kz_sub, par);
	invM(par->invPsi_super, par->Psi_super, 4*par->vec_size);

	if (par->READ_tab_NS){
		if (lire_tab(par->tab_NS_filename, "tab_NS", par->tab_NS, par->NS+1) != 0) {
			fprintf(stderr,"tab_NS reading ERROR\n");
			exit(EXIT_FAILURE);
		}
	}
	if (par->imposed_S_steps){
		if (lire_tab(par->imposed_S_steps_filename, "", par->tab_imposed_S_steps, par->N_imposed_S_steps) != 0) {
			fprintf(stderr,"tab_imposed_S_steps reading ERROR\n");
			exit(EXIT_FAILURE);
		}
		md2D_make_tab_S_steps(par);
	}
/*SaveDbleTab2file(par->tab_NS, par->NS+1, "stdout", " ");*/
/*printf("\nRe(par->Psi_super) :\n");
SaveMatrix2file (par->Psi_super, 4*par->vec_size, 4*par->vec_size, "Re", "stdout");
printf("\nRe(par->Psi_sub) :\n");
SaveMatrix2file (par->Psi_sub, 4*par->vec_size, 4*par->vec_size, "Re", "stdout");
printf("\nRe(par->invPsi_super) :\n");
SaveMatrix2file (par->invPsi_super, 4*par->vec_size, 4*par->vec_size, "Re", "stdout");
*/		
	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_alloc(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *	\brief	Memory allocation
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_alloc(struct Param_struct *par, struct Efficacites_struct *eff)
{
	int i;
	
	par->k2        = (REAL complex *) malloc(sizeof(REAL complex)*par->N_x);
	par->invk2     = (REAL complex *) malloc(sizeof(REAL complex)*par->N_x);
	par->Nx2       = (REAL complex *) malloc(sizeof(REAL complex)*par->N_x);
	par->Nz2       = (REAL complex *) malloc(sizeof(REAL complex)*par->N_x);
	par->NxNz      = (REAL complex *) malloc(sizeof(REAL complex)*par->N_x);

	par->sigma    = (REAL complex *) malloc(sizeof(REAL complex)*par->vec_size);
	par->kz_super = (REAL complex *) malloc(sizeof(REAL complex)*par->vec_size);
	par->kz_sub   = (REAL complex *) malloc(sizeof(REAL complex)*par->vec_size);

	par->TF_k2     = (REAL complex *) malloc(sizeof(REAL complex)*(4*par->N+1));
	par->TF_invk2  = (REAL complex *) malloc(sizeof(REAL complex)*(4*par->N+1));
	par->TF_Nx2    = (REAL complex *) malloc(sizeof(REAL complex)*(4*par->N+1));
	par->TF_Nz2    = (REAL complex *) malloc(sizeof(REAL complex)*(4*par->N+1));
	par->TF_NxNz   = (REAL complex *) malloc(sizeof(REAL complex)*(4*par->N+1));

	par->tmp_tf_k2    = (REAL complex *) malloc(sizeof(REAL complex)*par->N_x);   
	par->tmp_tf_invk2 = (REAL complex *) malloc(sizeof(REAL complex)*par->N_x);
	par->tmp_tf_Nx2   = (REAL complex *) malloc(sizeof(REAL complex)*par->N_x);
	par->tmp_tf_Nz2   = (REAL complex *) malloc(sizeof(REAL complex)*par->N_x);
	par->tmp_tf_NxNz  = (REAL complex *) malloc(sizeof(REAL complex)*par->N_x);
	par->M_tmp1        = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->M_tmp2        = allocate_CplxMatrix(par->vec_size,par->vec_size);

	par->Toep_k2       = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Toep_invk2    = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Toep_Nx2      = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Toep_NxNz     = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Toep_Nz2      = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->invToep_k2    = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->invToep_invk2 = allocate_CplxMatrix(par->vec_size,par->vec_size);

	par->S12 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	par->S22 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	par->S21 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	par->S11 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
/*** DM *****/

	par->Qxx = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Qyy = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Qzz = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Qxz = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Qzz_1 = allocate_CplxMatrix(par->vec_size,par->vec_size);

	par->QxzEx         = (REAL complex *) malloc(sizeof(REAL complex)*par->vec_size);
	par->Qzz_1QxzEx    = (REAL complex *) malloc(sizeof(REAL complex)*par->vec_size);
	par->Qzz_1Hpx      = (REAL complex *) malloc(sizeof(REAL complex)*par->vec_size);
	par->sigmaHpy      = (REAL complex *) malloc(sizeof(REAL complex)*par->vec_size); 
	par->ky0Qzz_1Hpx   = (REAL complex *) malloc(sizeof(REAL complex)*par->vec_size);
	par->Qzz_1sigmaHpy = (REAL complex *) malloc(sizeof(REAL complex)*par->vec_size);
	par->QxxEx         = (REAL complex *) malloc(sizeof(REAL complex)*par->vec_size);
	par->QyyEy         = (REAL complex *) malloc(sizeof(REAL complex)*par->vec_size);

	par->QxzVtmp1      = (REAL complex *) malloc(sizeof(REAL complex)*par->vec_size);
	par->V_tmp1        = (REAL complex *) malloc(sizeof(REAL complex)*par->vec_size);


/**RCWA********/

/**********/

	par->T         = allocate_CplxMatrix(4*par->vec_size,4*par->vec_size);

	par->T11 = (REAL complex **) malloc(sizeof(REAL complex *)*2*par->vec_size);
	par->T12 = (REAL complex **) malloc(sizeof(REAL complex *)*2*par->vec_size);
	par->T21 = (REAL complex **) malloc(sizeof(REAL complex *)*2*par->vec_size);
	par->T22 = (REAL complex **) malloc(sizeof(REAL complex *)*2*par->vec_size);
/*	par->T11 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	par->T12 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	par->T21 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	par->T22 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);*/
	for (i=0;i<=2*par->vec_size-1;i++){
			par->T11[i] = &par->T[i][0];
			par->T12[i] = &par->T[i][2*par->vec_size];
			par->T21[i] = &par->T[i+2*par->vec_size][0];
			par->T22[i] = &par->T[i+2*par->vec_size][2*par->vec_size];
	}

	par->M         = allocate_CplxMatrix(4*par->vec_size,4*par->vec_size);

	par->eig_values = (REAL complex *) malloc(sizeof(REAL complex)*4*par->vec_size);
	par->EigVectors = allocate_CplxMatrix(4*par->vec_size,4*par->vec_size);
	par->eig_buffer = (REAL complex *) malloc(sizeof(REAL complex)*50*4*par->vec_size);
	
	par->invEigVec         = allocate_CplxMatrix(4*par->vec_size,4*par->vec_size);
	par->invVec_Psi        = allocate_CplxMatrix(4*par->vec_size,4*par->vec_size);
	par->M_invVec_Psi      = allocate_CplxMatrix(4*par->vec_size,4*par->vec_size);
	par->M_buffer_4vecsize = allocate_CplxMatrix(4*par->vec_size,4*par->vec_size);
	par->M_sol = (REAL complex *) malloc(sizeof(REAL complex)*4*par->vec_size);


	par->Psi_sub      = allocate_CplxMatrix(4*par->vec_size,4*par->vec_size);
	par->Psi_super    = allocate_CplxMatrix(4*par->vec_size,4*par->vec_size);
	par->invPsi_super = allocate_CplxMatrix(4*par->vec_size,4*par->vec_size);

	par->Ai = (REAL complex *) malloc(sizeof(REAL complex)*(2*par->vec_size));
	par->Ar = (REAL complex *) malloc(sizeof(REAL complex)*(2*par->vec_size));
	par->At = (REAL complex *) malloc(sizeof(REAL complex)*(2*par->vec_size));	
	par->Vi = (REAL complex *) malloc(sizeof(REAL complex)*(2*par->vec_size));
	par->Vr = (REAL complex *) malloc(sizeof(REAL complex)*(2*par->vec_size));
	par->Vt = (REAL complex *) malloc(sizeof(REAL complex)*(2*par->vec_size));	

	par->Exi  = (REAL complex *) malloc(sizeof(REAL complex)*(par->vec_size));	
	par->Exr  = (REAL complex *) malloc(sizeof(REAL complex)*(par->vec_size));	
	par->Ext  = (REAL complex *) malloc(sizeof(REAL complex)*(par->vec_size));	
	par->Hpxi = (REAL complex *) malloc(sizeof(REAL complex)*(par->vec_size));	
	par->Hpxr = (REAL complex *) malloc(sizeof(REAL complex)*(par->vec_size));	
	par->Hpxt = (REAL complex *) malloc(sizeof(REAL complex)*(par->vec_size));	

	eff->eff_r   = (double *) malloc(sizeof(double)*(par->vec_size));
	eff->eff_t   = (double *) malloc(sizeof(double)*(par->vec_size));
	eff->N_eff_r = (double *) malloc(sizeof(double)*(par->vec_size));
	eff->N_eff_t = (double *) malloc(sizeof(double)*(par->vec_size));
	eff->theta_r = (double *) malloc(sizeof(double)*(par->vec_size));
	eff->theta_t = (double *) malloc(sizeof(double)*(par->vec_size));
	eff->phi_r   = (double *) malloc(sizeof(double)*(par->vec_size));
	eff->phi_t   = (double *) malloc(sizeof(double)*(par->vec_size));

	par->var_i  = (double *) malloc(sizeof(double)*1000); /* A REVOIR */
	par->var_i2 = (double *) malloc(sizeof(double)*1000); /* A REVOIR */

	if (par->STOCKER_TF) {
		par->z2n[0] = (long int *) malloc(sizeof(long int)*(par->BLOC_TAILLE_z2n));
		par->z2n[1] = (long int *) malloc(sizeof(long int)*(par->BLOC_TAILLE_z2n));
		par->tab_TF_k2 = allocate_CplxMatrix(par->BLOC_TAILLE_z2n, 4*par->N+1);
		par->tab_TF_invk2 = allocate_CplxMatrix(par->BLOC_TAILLE_z2n, 4*par->N+1);
	}
	if (par->tab_NS_ENABLED){
		par->tab_NS = (double *) malloc(sizeof(double)*(par->NS+1));
	}
	if (par->imposed_S_steps){
		par->tab_imposed_S_steps = (double *) malloc(sizeof(double)*(par->N_imposed_S_steps));
	}
	
	return 0;
}



/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_free(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *	\brief	Libération de la mémoire
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_free(struct Param_struct *par, struct Efficacites_struct *eff)
{
	if (par->verbosity>=3) fprintf(stdout,"Libération de la mémoire : "); fflush(stdout);

	
	if (par->type_profil == N_XYZ) {
		free(par->n_xyz[0]);
		free(par->n_xyz);
	}else{
		free(par->profil[0]);
		free(par->profil);
		free(par->k2_layer);
		free(par->invk2_layer);
	}
	if (par->tab_NS_ENABLED){
		free(par->tab_NS);
	}
	if (par->imposed_S_steps){
		free(par->tab_imposed_S_steps);
	}

	if (par->STOCKER_TF) {
		free(par->z2n[0]);
		free(par->z2n[1]);
		free(par->tab_TF_k2[0]);
		free(par->tab_TF_k2);
		free(par->tab_TF_invk2[0]);
		free(par->tab_TF_invk2);
	}

	free(par->k2);
	free(par->invk2);
	free(par->Nx2);
	free(par->Nz2);
	free(par->NxNz);

	free(par->sigma);
	free(par->kz_super);
	free(par->kz_sub);
		
	free(par->TF_k2);
	free(par->TF_invk2);
	free(par->TF_Nx2);
	free(par->TF_Nz2);
	free(par->TF_NxNz);
	
	free(par->tmp_tf_k2);
	free(par->tmp_tf_invk2);
	free(par->tmp_tf_Nx2);
	free(par->tmp_tf_Nz2);
	free(par->tmp_tf_NxNz);

	free(par->Toep_k2[0]);
	free(par->Toep_k2);
	free(par->Toep_invk2[0]);
	free(par->Toep_invk2);
	free(par->invToep_k2[0]);
	free(par->invToep_k2);
	free(par->invToep_invk2[0]);
	free(par->invToep_invk2);
	free(par->Toep_Nx2[0]);
	free(par->Toep_Nx2);
	free(par->Toep_NxNz[0]);
	free(par->Toep_NxNz);
	free(par->Toep_Nz2[0]);
	free(par->Toep_Nz2);
	free(par->M_tmp1[0]);
	free(par->M_tmp1);
	free(par->M_tmp2[0]);
	free(par->M_tmp2);

	free(par->Qxx[0]);
	free(par->Qyy[0]);
	free(par->Qzz[0]);
	free(par->Qxz[0]);
	free(par->Qzz_1[0]);
	free(par->Qxx);
	free(par->Qyy);
	free(par->Qzz);
	free(par->Qxz);
	free(par->Qzz_1);

	free(par->QxzEx);
	free(par->Qzz_1QxzEx);
	free(par->Qzz_1Hpx);
	free(par->V_tmp1);
	free(par->sigmaHpy); 
	free(par->ky0Qzz_1Hpx);
	free(par->Qzz_1sigmaHpy);
	free(par->QxxEx);
	free(par->QyyEy);
	free(par->QxzVtmp1);

	free(par->S12[0]);
	free(par->S12);
	free(par->S22[0]);
	free(par->S22);
	free(par->S11[0]);
	free(par->S11);
	free(par->S21[0]);
	free(par->S21);

	free(par->T11);
	free(par->T12);
	free(par->T21);
	free(par->T22);
	
	free(par->M[0]);
	free(par->M);

	free(par->T[0]);
	free(par->T);
	free(par->EigVectors[0]);
	free(par->EigVectors);
	free(par->invEigVec[0]);
	free(par->invEigVec);
	free(par->invVec_Psi[0]);
	free(par->invVec_Psi);
	free(par->M_invVec_Psi[0]);
	free(par->M_invVec_Psi);
	free(par->M_buffer_4vecsize[0]);
	free(par->M_buffer_4vecsize);
	free(par->eig_values);
	free(par->eig_buffer);
	free(par->M_sol);

	free(par->Psi_sub[0]);
	free(par->Psi_sub);
	free(par->Psi_super[0]);
	free(par->Psi_super);
	free(par->invPsi_super[0]);
	free(par->invPsi_super);
	
	free(par->Ai);
	free(par->Ar);
	free(par->At);
	free(par->Vi);
	free(par->Vr);
	free(par->Vt);
	free(par->Exi);
	free(par->Exr);
	free(par->Ext);
	free(par->Hpxi);
	free(par->Hpxr);
	free(par->Hpxt);

	free(eff->eff_r);
	free(eff->eff_t);
	free(eff->N_eff_r);
	free(eff->N_eff_t);
	free(eff->theta_r);
	free(eff->theta_t);
	free(eff->phi_r);
	free(eff->phi_t);
	
	if (par->verbosity>=3) fprintf(stdout,"OK\n"); fflush(stdout);

	
	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn	int md2D_read_mat_S(struct Param_struct *par)
 *
 *	\brief	Lecture de matrices S
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_read_mat_S(struct Param_struct *par)
{
	double **re_tmp, **im_tmp;
	int i,j;
	char *nom_fichier = par->mat_S_file;
	char *erreur="NO_ERROR                     ";
	FILE *fp;
	
	/* Ouverture du fichier */
	if (!(fp = fopen(nom_fichier,"r"))){
		fprintf(stderr, "%s ligne %d : ERREUR, impossible d'ouvrir %s\n",__FILE__, __LINE__,nom_fichier);
		return 1;
	}
	/* Lecture des paramètres */
	if (lire_string (fp, "mat_S_name", par->mat_S_name)) erreur="mat_S_name";
	if (lire_double (fp, "h_partial", &(par->h) )) erreur="h_partial";
/*	if (lire_double (fp, "L", &(par->L) )) erreur="L";
	if (lire_double (fp, "Re_k0", &Re_k0 )) erreur="Re_k0";
	if (lire_double (fp, "Re_kh", &Re_kh )) erreur="Re_kh";
	if (lire_double (fp, "Im_k0", &Re_kh )) erreur="Im_k0";
	if (lire_double (fp, "Im_kh", &Re_kh )) erreur="Im_kh";
	if (lire_double (fp, "delta_sigma", &(par->delta_sigma) )) erreur="delta_sigma";
	if (lire_double (fp, "sigma0", &(par->sigma0) )) erreur="sigma0";
	if (lire_int (fp, "N", &(par->N) )) erreur="N";
*/	fclose(fp);
	/* Vérification de l'absence d'erreurs de lecture */
	if (strcmp(erreur,"NO_ERROR                     ")){
		fprintf(stderr, "%s : Erreur, probleme de lecture de \"%s\"\n",__FILE__,erreur);
		exit(EXIT_FAILURE);
	}
	
	/* Lecture des éléments de matrice S */	
	re_tmp = allocate_DbleMatrix(2*par->N+1,2*par->N+1);
	im_tmp = allocate_DbleMatrix(2*par->N+1,2*par->N+1);
	/* S12 */
	lire_tab(nom_fichier, "Re_S12", re_tmp[0], (2*par->N+1)*(2*par->N+1));
	lire_tab(nom_fichier, "Im_S12", im_tmp[0], (2*par->N+1)*(2*par->N+1));
	for (j=0;j<=2*par->N;j++){
		for (i=0;i<=2*par->N;i++){
			par->S12[i][j] = c_omplex(re_tmp[i][j],im_tmp[i][j]);
		}
	}
	/* S22 */
	lire_tab(nom_fichier, "Re_S22", re_tmp[0], (2*par->N+1)*(2*par->N+1));
	lire_tab(nom_fichier, "Im_S22", im_tmp[0], (2*par->N+1)*(2*par->N+1));
	for (j=0;j<=2*par->N;j++){
		for (i=0;i<=2*par->N;i++){
			par->S22[i][j] = c_omplex(re_tmp[i][j],im_tmp[i][j]);
		}
	}

	free(re_tmp[0]);free(re_tmp);
	free(im_tmp[0]);free(im_tmp);
	
	return 0;		
}




