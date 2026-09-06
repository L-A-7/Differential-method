/*---------------------------------------------------------------------------------------------*/
/*!	\file		md2D_pilot.c
 *
 * 	\brief		Piloting md2D 
 */
/*---------------------------------------------------------------------------------------------*/
#include "md2D_pilot.h"
#include "md2D.h"

/*---------------------------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------------------------*/
int main(int argc, char *argv[]){

	struct Noms_fichiers nomfichier;
	struct Param_struct param;
	struct Efficacites_struct effic;

	/*param.verbosity =1;*/
	
	/* Command line arguments copy */
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

	/* Program initialisation: inputs reading, memory allocations, etc. */
	md2D_init(&param, &effic, &nomfichier);

/*********************************** DEBUG *************************************************************/
/*check_complex(1);
matrix_operations_check();
SaveMatrix2file (Psi, 2*vec_size, 2*vec_size, "Re", "stdout");
SaveCplxTab2file (Psi22, vec_size, "Re", "stdout", " ", 10000, "\n");*/
/*
read_S_steps_from_profile(&param);
md2D_make_tab_S_steps(&param);

printf("N_imposed_S_steps = %d\n",param.N_imposed_S_steps);
printf("NS = %d\n",param.NS);
printf("NS_estimated = %d\n",param.NS_estimated);
printf("tab_NS_ENABLED = %d\n",param.tab_NS_ENABLED);
printf("hmin_Sstep = %f\n",param.hmin_Sstep);
printf("imposed_S_steps = %s\n",param.imposed_S_steps);
printf("imposed_S_steps_filename = %s\n",param.imposed_S_steps_filename);
fflush(stdout);

printf("\ntab_imposed_S_steps\n");SaveDbleTab2file (param.tab_imposed_S_steps, param.N_imposed_S_steps, "stdout", " ", 100, "\n");
printf("\ntab_NS\n");SaveDbleTab2file (param.tab_NS, param.NS, "stdout", " ", 100, "\n");
read_S_steps_from_profile(&param);
md2D_make_tab_S_steps(&param);
printf("\ntab_imposed_S_steps\n");SaveDbleTab2file (param.tab_imposed_S_steps, param.N_imposed_S_steps, "stdout", " ", 100, "\n");
printf("\ntab_NS\n");SaveDbleTab2file (param.tab_NS, param.NS, "stdout", " ", 100, "\n");

int vec_size=param.vec_size;

/**** sigma checking ****/
/*printf("\nRe(sigma)\n");SaveCplxTab2file (param.sigma, vec_size, "Re", "stdout", " ", 10000, "\n");
printf("\nIm(sigma)\n");SaveCplxTab2file (param.sigma, vec_size, "Im", "stdout", " ", 10000, "\n");*/

/**** M-matrix checking ****/
/*M_matrix_TE(param.M, param.h/2.0, &param);
M_matrix_TM(param.M, param.h/2.0, &param);
zinvar_M_matrix_TM(param.M, param.h/2.0, &param);
printf("\nRe(M)\n");SaveMatrix2file (param.M, 2*vec_size, 2*vec_size, "Re", "stdout");
printf("\nIm(M)\n");SaveMatrix2file (param.M, 2*vec_size, 2*vec_size, "Im", "stdout");*/

/**** Psi-matrix checking ****/
/*COMPLEX *Psi11, *Psi12, *Psi21, *Psi22;
Psi11 = malloc(sizeof(COMPLEX)*vec_size); Psi12 = malloc(sizeof(COMPLEX)*vec_size);
Psi21 = malloc(sizeof(COMPLEX)*vec_size); Psi22 = malloc(sizeof(COMPLEX)*vec_size);
PsiMatrix(Psi11, Psi12, Psi21, Psi22, param.kz_super, param.k_super, vec_size, param.pola);
printf("\nRe(Psi11)\n");SaveCplxTab2file (Psi11, vec_size, "Re", "stdout", " ", 10000, "\n");
printf("\nIm(Psi11)\n");SaveCplxTab2file (Psi11, vec_size, "Im", "stdout", " ", 10000, "\n");
printf("\nRe(Psi12)\n");SaveCplxTab2file (Psi12, vec_size, "Re", "stdout", " ", 10000, "\n");
printf("\nIm(Psi12)\n");SaveCplxTab2file (Psi12, vec_size, "Im", "stdout", " ", 10000, "\n");
printf("\nRe(Psi21)\n");SaveCplxTab2file (Psi21, vec_size, "Re", "stdout", " ", 10000, "\n");
printf("\nIm(Psi21)\n");SaveCplxTab2file (Psi21, vec_size, "Im", "stdout", " ", 10000, "\n");
printf("\nRe(Psi22)\n");SaveCplxTab2file (Psi22, vec_size, "Re", "stdout", " ", 10000, "\n");
printf("\nIm(Psi22)\n");SaveCplxTab2file (Psi22, vec_size, "Im", "stdout", " ", 10000, "\n");*/

/**** S-matrix checking ****/
/*S_matrix(&param);
printf("\nRe(S12)\n");SaveMatrix2file (param.S12, vec_size, vec_size, "Re", "stdout");
printf("\nIm(S12)\n");SaveMatrix2file (param.S12, vec_size, vec_size, "Im", "stdout");
printf("\nRe(S22)\n");SaveMatrix2file (param.S22, vec_size, vec_size, "Re", "stdout");
printf("\nIm(S22)\n");SaveMatrix2file (param.S22, vec_size, vec_size, "Im", "stdout");
*/
/*	(param.k_2)(&param, param.k2, param.h/2);
	(param.invk_2)(&param, param.invk2, param.h/2);
	printf("\nRe(k2)\n");SaveCplxTab2file (param.k2, param.N_x, "Re", "stdout", "\n", 10000, "\n");*/
/*	return 0;*/
/*******************************************************************************************************/

	/* Calcul type. Use STD or NEAR_FIELD unless you know what you're doing */
	if (!strcmp(param.calcul_type,"STD")){
		md2D_classical_FFF (&param, &effic, &nomfichier);
	}else if (!strcmp(param.calcul_type,"GUIDED")){
		md2D_guided (&param, &effic, &nomfichier);
	}else if (!strcmp(param.calcul_type,"SWIFTS") || !strcmp(param.calcul_type,"I_SWIFTS")){
		md2D_swifts (&param, &effic, &nomfichier);
	}else if (!strcmp(param.calcul_type,"NEAR_FIELD")){
		md2D_near_field (&param, &effic, &nomfichier);
	}else{
		fprintf(stderr, "%s, line %d : ERROR, unknown calculation type (\"%s\")\n",__FILE__,__LINE__,param.calcul_type);
		exit(EXIT_FAILURE);
	}


/*************************************** TESTS *********************************************************/
/*fprintf(stdout,"\nReS12 :\n");SaveMatrix2file (param.S12, param.vec_size, param.vec_size, "Re", "stdout");
fprintf(stdout,"\nReS22 :\n");SaveMatrix2file (param.S22, param.vec_size, param.vec_size, "Re", "stdout");
fprintf(stdout,"\nRePsisuper_TE :\n");SaveMatrix2file (param.Psi_super_TE, 2*param.vec_size, 2*param.vec_size, "Re", "stdout");
fprintf(stdout,"\nRePsisuper_TM :\n");SaveMatrix2file (param.Psi_super_TM, 2*param.vec_size, 2*param.vec_size, "Re", "stdout");*/


/************************************* FIN TESTS *******************************************************/

	/* Freeing memory */
	md2D_free(&param, &effic);
	
	/* Showing used time & efficiencies sum */
	if (param.verbosity >= 1){
		printf("Efficiencies summ  : %1.10f\n",effic.sum_eff);
		printf("1-Efficiencies summ: %e\n",1-effic.sum_eff);
		fprintf(stdout,"Ellapsed time: %f s \n", md_chrono(&param));
	}
	
	return 0;		
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		md2D_classical_FFF (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *		\brief	Complex field calculation, in the case of a 1D structure with conical incidence.
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_classical_FFF (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
{
	/* Propagative orders limits calculation */
	md2D_propagativ_limits(par, eff);
		
	/* Incident field amplitude */
	md2D_incident_field(par,eff);

	S_matrix(par);

	/* Amplitudes calculation */
	md2D_amplitudes(par->Ai, par->Ar, par->At, par->S12, par->S22, par);

	/* Efficiencies calculation */
	md2D_efficiencies(par->Ai, par->Ar, par->At, par, eff);

	/* Writting results in fichier_results */
	md2D_genere_nom_fichier_results(nomfichier->fichier_results, par);

	md2D_ecrire_results(nomfichier->fichier_results, par, eff);

	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		md2D_guided (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *		\brief	Specific case for guided structures
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_guided (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
{
	/* Propagative orders limits calculation */
	md2D_propagativ_limits(par, eff);
		
	S_matrix(par);

	/* Energy flux calculation (which are stored as efficiencies to minimize code modifications) */
	md2D_energy_flux(par, eff);

	/* Writting results in fichier_results */
	md2D_genere_nom_fichier_results(nomfichier->fichier_results, par);

	md2D_ecrire_results(nomfichier->fichier_results, par, eff);

	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_swifts (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *		\brief	Specific case for guided structures
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_swifts (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
{
/*	COMPLEX n_guide=1.53;
	double h_guide=1000.0;*/

 
	/* Propagative orders limits calculation */
	md2D_propagativ_limits(par, eff);
		
	S_matrix(par);

	/* Energy flux calculation (which are stored as efficiencies to minimize code modifications) */
	swifts_md2D_energy_flux(par->h_layer,par, eff);

	/* Writting results in fichier_results */
	md2D_genere_nom_fichier_results(nomfichier->fichier_results, par);

	md2D_ecrire_results(nomfichier->fichier_results, par, eff);

	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		md2D_near_field (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *		\brief	Complex field calculation, in the case of a 1D structure with conical incidence.
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_near_field (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
{
	par->tab_Z   = allocate_CplxMatrix_3(par->vec_size, par->vec_size, par->NS);
	par->tab_S12 = allocate_CplxMatrix_3(par->vec_size, par->vec_size, par->NS);
	COMPLEX *V0m;

	/* Propagative modes limits calculation */
	md2D_propagativ_limits(par, eff);

	/* Incident field amplitudes */
	md2D_incident_field(par,eff);

	S_matrix(par);
	md2D_amplitudes(par->Ai, par->Ar, par->At, par->S12, par->S22, par);
	md2D_efficiencies(par->Ai, par->Ar, par->At, par, eff);
	md2D_genere_nom_fichier_results(nomfichier->fichier_results, par);
	md2D_ecrire_results(nomfichier->fichier_results, par, eff);

	V0m=par->At; /* because h0=0 */
	md2D_local_field_map(par->tab_S12, par->tab_Z, par); /* Normal configuration (far field incident) */
/*	md2D_local_field_map_by_T_products(V0m, par);*/ /* for configuration with waveguide */
/*	md2D_near_field_map_evanescent(par);*/


	free(par->tab_Z[0][0]);
	free(par->tab_Z[0]);
	free(par->tab_Z);
	free(par->tab_S12[0][0]);
	free(par->tab_S12[0]);
	free(par->tab_S12);
	

	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_init (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier) 
 *
 *	\brief	Program nitialisation
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_init (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier) 
{

	/* Reading default parameters in fichier_param */
	md2D_lire_param(nomfichier, par);

	/* Memory allocation for the profile */
	md2D_alloc_init_profil(par);

	/* Reading the profile h(x) of the surface */
	(*par->md2D_lire_profil)(nomfichier->profile_file, par);

	/* Initialisating some variables */
	md2D_variables_init(par, eff);

	/* Allocating memory for arrays */
	md2D_alloc(par, eff);

	/* Initialising some arrays */
	md2D_arrays_init(par, eff);
	
	/* Showing on screen some calculated and initial parameters */
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
	/* Variable initilisation */
	par->k_super = 2*PI*par->n_super/par->lambda;
	par->k_sub = 2*PI*par->n_sub/par->lambda;
	par->Delta_sigma = 2*PI/par->L;
	par->sigma0 = par->k_super*sin(par->theta_i)*cos(par->phi_i);
	if (cabs(par->sigma0_normed/SIGMA0_NORMED_NOT_DEFINED-1) > EPS){
		fprintf(stdout,"*** sigma0_normed = %1.10f + i%f ***\n",creal(par->sigma0_normed),cimag(par->sigma0_normed));
		par->sigma0 = 2*PI*par->sigma0_normed/par->lambda;
	}
/*	par->ky_0   = par->k_super*sin(par->theta_i)*sin(par->phi_i);*/
	par->ky_0   = par->sigma0*tan(par->phi_i); /* TODO: A verifier */
			
	/* Allocations */	
	if (par->type_profil == N_XYZ) {
		par->n_xyz = allocate_CplxMatrix(par->N_z,par->N_x);
	}else{
		par->profil = allocate_DbleMatrix(par->N_layers+3,par->N_x);
		par->k2_layer    = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->N_layers+2));
		par->invk2_layer = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->N_layers+2));
	}


	/* Function's pointers alignement */
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
			par->Normal_function = Normal_H_X_Multi;
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
 *	\brief	Variable initialisations
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_variables_init(struct Param_struct *par, struct Efficacites_struct *eff)
{
	double K,sig0;
	par->clock0 = clock(); /* Initialisation of chronometers                  */
	time(&(par->time0));   /* clock0 (short duration) and time0 (long durations) */
	par->last_clock = clock();
	time(&(par->last_time)); 

	/* Initialisating variables in mode "AUTO" */
	/* delta_h = lambda x (n_re + n_im) / 1000 */
	if (par->delta_h == AUTO) {
		par->delta_h = par->lambda/(cabs(par->n_sub)*1000);
	}
	/* NS = (h/lambda) x (n_re + n_im) x 5 */
	if (par->NS == AUTO) {
		par->NS = CEIL( (par->h / par->lambda)*cabs(par->n_sub)*5 );
	}
	/* N : we add 10% of evanescent modes (and 3 at minimum) */
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
		/* we add 10% of evanescent modes (and 3 at minimum) */
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
	par->Ni = 1; /* Pour compatibilite (une autre valeur sera affectee par les fonctions var_i_... )*/
	
	/* Calculation of the value of delta_h so that there is an integer value of S-Matrix steps */
	par->N_steps = 0; /* Counter initialisation */
/*	par->Nstep_S = ROUND(ceil((par->h/par->NS)/par->delta_h));
	par->delta_h = (par->h/par->NS)/par->Nstep_S;
	par->Nstep = ROUND(par->h/par->delta_h);*/

	/* tab_NS_ENABLED */
	if (!strcmp(par->imposed_S_steps,"AUTO") || !strcmp(par->imposed_S_steps,"FROM_FILE")){
		par->tab_NS_ENABLED = 1;
	}else{
		par->tab_NS_ENABLED = 0;
	}
	if (!strcmp(par->imposed_S_steps,"FROM_FILE")){
		if (lire_tab(par->imposed_S_steps_filename, "tab_imposed_S_steps", par->tab_imposed_S_steps, par->N_imposed_S_steps) != 0){
			fprintf(stderr, "%s line %d: ERROR, can't read tab_imposed_S_steps in %s\n",__FILE__, __LINE__,par->imposed_S_steps_filename);
			exit(EXIT_FAILURE);
		}
		par->NS_estimated = par->N_imposed_S_steps + ROUND(par->h/par->hmin_Sstep); 
	}

	if (!strcmp(par->imposed_S_steps,"AUTO")){
		if (fabs(par->hmin_Sstep/SUPER_BIG_HMIN_SSTEP -1) > EPS){ /* if hmin_Sstep has been defined */
			par->NS_estimated = par->N_x + CEIL(par->h/par->hmin_Sstep); 
		}else{
			par->NS_estimated = par->N_x; 
		}
	}
	
	/*par->verbosity = 2;*/ /* Si verbosity > 0 : affiche plus d'infos sur le terminal */

	if (!strcmp(par->calcul_type,"SWIFTS") || !strcmp(par->calcul_type,"I_SWIFTS")){
		strncpy(par->calculation_on_layer, "ON_1_LAYER", SIZE_STR_BUFFER);
	}
	if (!strcmp(par->calcul_type,"I_SWIFTS")){
		strncpy(par->S_matrix_blocks_calculation, "ALL_4_S_MATRIX_BLOCKS", SIZE_STR_BUFFER);
	}

	/* functions pointers alignements */
	if (par->pola == TE){
			par->M_matrix = M_matrix_TE;
	}else if (!strcmp(par->calcul_method,"Z_INVAR") && par->pola == TM){
			par->M_matrix = zinvar_M_matrix_TM;
	}else{
			par->M_matrix = M_matrix_TM;
	}

	if (!strcmp(par->calcul_method,"RK4")){
			par->P_matrix = rk4_P_matrix;
/*	}else if (!strcmp(par->calcul_method,"SHOOTING_METHOD")){
			par->P_matrix = shooting_P_matrix;*/
	}else if (!strcmp(par->calcul_method,"IMPLICIT_RK")){
			par->P_matrix = implicit_rk_P_matrix;
	}else if (!strcmp(par->calcul_method,"Z_INVAR")){
			par->P_matrix = zinvar_P_matrix;
	}else if (!strcmp(par->calcul_method,"IMPROVED_RCWA")){
			par->P_matrix = zinvar_P_matrix;
	}else if (!strcmp(par->calcul_method,"EULER")){
			par->P_matrix = euler_P_matrix;
	}else if (!strcmp(par->calcul_method,"EULER2")){
			par->P_matrix = euler2_P_matrix;
	}else{
		fprintf(stderr, "%s, line %d: ERROR, unknown calculation method (\"%s\")\n",__FILE__,__LINE__,par->calcul_method);
		exit(EXIT_FAILURE);
	}

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
	
	COMPLEX k_super2 = par->k_super*par->k_super;
	COMPLEX k_sub2 = par->k_sub*par->k_sub;
	COMPLEX ky_02 = par->ky_0*par->ky_0;
	
	for(j=0;j<=par->vec_size-1;j++){
		par->sigma[j] = (j-par->vec_middle)*par->Delta_sigma + par->sigma0;
		par->kz_super[j] = csqrt(k_super2 - par->sigma[j]*par->sigma[j] - ky_02);
		par->kz_sub[j]   = csqrt(k_sub2   - par->sigma[j]*par->sigma[j] - ky_02);
	}

	if (!strcmp(par->imposed_S_steps,"FROM_FILE")){
		md2D_make_tab_S_steps(par);
	}

	if (!strcmp(par->imposed_S_steps,"AUTO")){
		read_S_steps_from_profile(par);
		md2D_make_tab_S_steps(par);
	}

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
	
	par->k2        = (COMPLEX *) malloc(sizeof(COMPLEX)*par->N_x);
	par->invk2     = (COMPLEX *) malloc(sizeof(COMPLEX)*par->N_x);
	par->Nx2       = (COMPLEX *) malloc(sizeof(COMPLEX)*par->N_x);
	par->Nz2       = (COMPLEX *) malloc(sizeof(COMPLEX)*par->N_x);
	par->NxNz      = (COMPLEX *) malloc(sizeof(COMPLEX)*par->N_x);
	if (par->smoothing){
		par->k2_tmp = (COMPLEX *) malloc(sizeof(COMPLEX)*par->N_x);}

	par->sigma    = (COMPLEX *) malloc(sizeof(COMPLEX)*par->vec_size);
	par->kz_super = (COMPLEX *) malloc(sizeof(COMPLEX)*par->vec_size);
	par->kz_sub   = (COMPLEX *) malloc(sizeof(COMPLEX)*par->vec_size);

	par->TF_k2     = (COMPLEX *) malloc(sizeof(COMPLEX)*(4*par->N+1));
	par->TF_invk2  = (COMPLEX *) malloc(sizeof(COMPLEX)*(4*par->N+1));
	par->TF_Nx2    = (COMPLEX *) malloc(sizeof(COMPLEX)*(4*par->N+1));
	par->TF_Nz2    = (COMPLEX *) malloc(sizeof(COMPLEX)*(4*par->N+1));
	par->TF_NxNz   = (COMPLEX *) malloc(sizeof(COMPLEX)*(4*par->N+1));

	par->tmp_tf_k2    = (COMPLEX *) malloc(sizeof(COMPLEX)*par->N_x);   
	par->tmp_tf_invk2 = (COMPLEX *) malloc(sizeof(COMPLEX)*par->N_x);
	par->tmp_tf_Nx2   = (COMPLEX *) malloc(sizeof(COMPLEX)*par->N_x);
	par->tmp_tf_Nz2   = (COMPLEX *) malloc(sizeof(COMPLEX)*par->N_x);
	par->tmp_tf_NxNz  = (COMPLEX *) malloc(sizeof(COMPLEX)*par->N_x);
	par->M_tmp1        = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->M_tmp2        = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->M_tmp3        = allocate_CplxMatrix(par->vec_size,par->vec_size);

	par->Toep_k2       = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Toep_invk2    = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Toep_Nx2      = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Toep_NxNz     = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Toep_Nz2      = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->invToep_k2    = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->invToep_invk2 = allocate_CplxMatrix(par->vec_size,par->vec_size);

	par->S12 = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->S22 = allocate_CplxMatrix(par->vec_size,par->vec_size);
	if (!strcmp(par->S_matrix_blocks_calculation,"ALL_4_S_MATRIX_BLOCKS")){
		par->S21 = allocate_CplxMatrix(par->vec_size,par->vec_size);
		par->S11 = allocate_CplxMatrix(par->vec_size,par->vec_size);
	}

	par->Qxx = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Qyy = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Qzz = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Qxz = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Qzz_1 = allocate_CplxMatrix(par->vec_size,par->vec_size);

	par->QxzEx         = (COMPLEX *) malloc(sizeof(COMPLEX)*par->vec_size);
	par->Qzz_1QxzEx    = (COMPLEX *) malloc(sizeof(COMPLEX)*par->vec_size);
	par->Qzz_1Hpx      = (COMPLEX *) malloc(sizeof(COMPLEX)*par->vec_size);
	par->sigmaHpy      = (COMPLEX *) malloc(sizeof(COMPLEX)*par->vec_size); 
	par->ky0Qzz_1Hpx   = (COMPLEX *) malloc(sizeof(COMPLEX)*par->vec_size);
	par->Qzz_1sigmaHpy = (COMPLEX *) malloc(sizeof(COMPLEX)*par->vec_size);
	par->QxxEx         = (COMPLEX *) malloc(sizeof(COMPLEX)*par->vec_size);
	par->QyyEy         = (COMPLEX *) malloc(sizeof(COMPLEX)*par->vec_size);

	par->QxzVtmp1      = (COMPLEX *) malloc(sizeof(COMPLEX)*par->vec_size);
	par->V_tmp1        = (COMPLEX *) malloc(sizeof(COMPLEX)*par->vec_size);

	par->T         = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);

	par->T11 = (COMPLEX **) malloc(sizeof(COMPLEX *)*par->vec_size);
	par->T12 = (COMPLEX **) malloc(sizeof(COMPLEX *)*par->vec_size);
	par->T21 = (COMPLEX **) malloc(sizeof(COMPLEX *)*par->vec_size);
	par->T22 = (COMPLEX **) malloc(sizeof(COMPLEX *)*par->vec_size);
	for (i=0;i<=par->vec_size-1;i++){
			par->T11[i] = &par->T[i][0];
			par->T12[i] = &par->T[i][par->vec_size];
			par->T21[i] = &par->T[i+par->vec_size][0];
			par->T22[i] = &par->T[i+par->vec_size][par->vec_size];
	}

	par->P         = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	par->M         = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);

	par->M_tmp11 = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->M_tmp12 = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->M_tmp21 = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->M_tmp22 = allocate_CplxMatrix(par->vec_size,par->vec_size);
	
	par->Ai = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));
	par->Ar = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));
	par->At = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));	
	par->Vi = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));
	par->Vr = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));
	par->Vt = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));	

	par->Exi  = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));	
	par->Exr  = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));	
	par->Ext  = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));	
	par->Hpxi = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));	
	par->Hpxr = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));	
	par->Hpxt = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));	

	eff->eff_r   = (double *) malloc(sizeof(double)*(par->vec_size));
	eff->eff_t   = (double *) malloc(sizeof(double)*(par->vec_size));
	eff->arg_Ar   = (double *) malloc(sizeof(double)*(par->vec_size));
	eff->arg_At   = (double *) malloc(sizeof(double)*(par->vec_size));
	eff->N_eff_r = (double *) malloc(sizeof(double)*(par->vec_size));
	eff->N_eff_t = (double *) malloc(sizeof(double)*(par->vec_size));
	eff->theta_r = (double *) malloc(sizeof(double)*(par->vec_size));
	eff->theta_t = (double *) malloc(sizeof(double)*(par->vec_size));
	eff->phi_r   = (double *) malloc(sizeof(double)*(par->vec_size));
	eff->phi_t   = (double *) malloc(sizeof(double)*(par->vec_size));

	par->var_i  = (double *) malloc(sizeof(double)*1000); /* A REVOIR */
	par->var_i2 = (double *) malloc(sizeof(double)*1000); /* A REVOIR */

	if (par->tab_NS_ENABLED){
		par->tab_NS = (double *) malloc(sizeof(double)*(par->NS_estimated+2));
	}
	if (!strcmp(par->imposed_S_steps,"AUTO")){
		par->tab_imposed_S_steps = (double *) malloc(sizeof(double)*(par->N_x+2));
	}else if (!strcmp(par->imposed_S_steps,"FROM_FILE")){
		par->tab_imposed_S_steps = (double *) malloc(sizeof(double)*(par->N_imposed_S_steps+2));
	}
printf("pointer par->tab_imposed_S_steps = %p\n",par->tab_imposed_S_steps);
	
	return 0;
}



/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_free(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *	\brief	Freeing memory
 */
/*---------------------------------------------------------------------------------------------*/
int md2D_free(struct Param_struct *par, struct Efficacites_struct *eff)
{
	if (par->verbosity>=3) fprintf(stdout,"Freeing memory: "); fflush(stdout);

	
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
	if (!strcmp(par->imposed_S_steps,"AUTO") || !strcmp(par->imposed_S_steps,"FROM_FILE")){
		free(par->tab_imposed_S_steps);
	}

	free(par->k2);
	free(par->invk2);
	free(par->Nx2);
	free(par->Nz2);
	free(par->NxNz);
	if (par->smoothing){
		free(par->k2_tmp);}

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
	free(par->M_tmp3[0]);
	free(par->M_tmp3);

	free(par->S12[0]);
	free(par->S12);
	free(par->S22[0]);
	free(par->S22);
	if (!strcmp(par->S_matrix_blocks_calculation,"ALL_4_S_MATRIX_BLOCKS")){
		free(par->S11[0]);
		free(par->S11);
		free(par->S21[0]);
		free(par->S21);
	}

	free(par->T11);
	free(par->T12);
	free(par->T21);
	free(par->T22);
	
	free(par->P[0]);
	free(par->P);
	free(par->M[0]);
	free(par->M);

	free(par->T[0]);
	free(par->T);

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

	free(par->M_tmp11[0]);
	free(par->M_tmp12[0]);
	free(par->M_tmp21[0]);
	free(par->M_tmp22[0]);
	free(par->M_tmp11);
	free(par->M_tmp12);
	free(par->M_tmp21);
	free(par->M_tmp22);
	
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
	free(eff->arg_Ar);
	free(eff->arg_At);
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
 *	\brief	Reading of S-Matrices
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
	/* Verifying absence of reading errors */
	if (strcmp(erreur,"NO_ERROR                     ")){
		fprintf(stderr, "%s : Erreur, probleme de lecture de \"%s\"\n",__FILE__,erreur);
		exit(EXIT_FAILURE);
	}
	
	/* Reading elements of the S-Matrix */	
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




