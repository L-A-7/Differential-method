/*---------------------------------------------------------------------------------------------*/
/*!	\file		md3D_pilot.c
 *
 * 	\brief		Pilotage de md3D 
 */
/*---------------------------------------------------------------------------------------------*/
#include "md3D.h"
/*---------------------------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------------------------*/
int main(int argc, char *argv[]){

	struct Noms_fichiers nomfichier;
	struct Param_struct param;
	struct Efficacites_struct effic;

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
	md3D_init(&param, &effic, &nomfichier);

/*********************************** DEBUG *************************************************************/
/*int M_matrix(COMPLEX **M, double z, struct Param_struct *par);
double z = 0.5*param.h;
int md3D_QMatrix(double z, COMPLEX **Qxx, COMPLEX **Qxy, COMPLEX **Qxz, COMPLEX **Qyy, COMPLEX **Qyz, COMPLEX **Qzz,COMPLEX **Qzz_1, struct Param_struct *par);
int Normal_H_XY(COMPLEX **norm_x, COMPLEX **norm_y, COMPLEX **norm_z, double z, struct Param_struct *par);*/

/*md3D_QMatrix(z, param.Qxx, param.Qxy, param.Qxz, param.Qyy, param.Qyz, param.Qzz,param.Qzz_1, &param);*/

/*M_matrix(param.M, z, &param);
COMPLEX *buffer_NprxNpry = malloc(sizeof(COMPLEX)*param.Npry*param.Nprx);*/

/*(*param.k_2)(&param, buffer_NprxNpry, z);
SaveCplxTab2file (buffer_NprxNpry, param.Nprx*param.Npry, "Re", "stdout", " ", param.Nprx, "\n");*/


/*SaveMatrix2file (param.M, 4*param.vec_size, 4*param.vec_size, "Im", "stdout");*/
/*SaveDbleTab2file (param.profil[0], param.Nprx*param.Npry, "stdout", " ", param.Nprx, "\n");*/
/*printf("\nIm(M)\n");SaveMatrix2file (param.M, 4*param.vec_size, 4*param.vec_size, "Im", "stdout");*/

/*return -12;*/
/*******************************************************************************************************/

	/* Choix du type de calcul */
	if (!strcmp(param.calcul_type,"STD")){
		md3D_std (&param, &effic, &nomfichier);
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
	md3D_free(&param, &effic);
	
	/* Showing used time & efficiencies sum */
	if (param.verbosity >= 1){
		printf("Efficiencies summ   : %1.10f\n",effic.sum_eff);
		printf("1-Efficiencies summ : %e\n",1-effic.sum_eff);
		fprintf(stdout,"Ellapsed time : %f s \n", md_chrono(&param));
	}
	
	return 0;		
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		md3D_classical_FFF (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
 *
 *		\brief	Complex field calculation, in the case of a 3D structure.
 */
/*---------------------------------------------------------------------------------------------*/
int md3D_std (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier)
{

	/* Amplitude du champ incident */
	md3D_incident_field(par,eff);

#if 0
	if (par->profile_type == H_XY_plus_STACK){
		S_matrix_stack(par); /* Case of multilayer structure with one structured stack situated at n_patterned_layer */
	}else{
		S_matrix(par);
	}
#endif
	
	S_matrix(par);

	/* Calcul des amplitudes */
	md3D_amplitudes(par->Ai, par->Ar, par->At, par->S12, par->S22, par);

	/* Calcul des efficacités */
	md3D_efficiencies(par->Ai, par->Ar, par->At, par, eff);

	/* Ecriture des résultats dans fichier_results */
	md3D_genere_nom_fichier_results(nomfichier->fichier_results, par);

	md3D_ecrire_results(nomfichier->fichier_results, par, eff);


	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md3D_init (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier) 
 *
 *	\brief	Initialisation du programme
 */
/*---------------------------------------------------------------------------------------------*/
int md3D_init (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier) 
{

	/* Lecture des paramètres par defaut dans param_file */
	md3D_lire_param(nomfichier, par);

	/* Allocation de mémoire pour le profil */
	md3D_alloc_init_profil(par);

	/* Lecture du profil h(x) décrivant la surface */
	(*par->md3D_lire_profil)(nomfichier->profile_file, par);

	/* Initialisations de certaines variables */
	md3D_variables_init(par, eff);

	/* Calcul des limites des modes propagatifs */
	md3D_propagativ_limits(par, eff);

	/* Allocation de mémoire pour les tableaux */
	md3D_alloc(par, eff);

	/* Initialising some arrays */
	md3D_arrays_init(par, eff);
	
	/* Affichage des paramètres lus et calculés */
	md3D_affiche_valeurs_param(par, nomfichier);

	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn	int md3D_alloc_init_profil(struct Param_struct *par)
 *
 *	\brief	Memory allocation for profil
 */
/*---------------------------------------------------------------------------------------------*/
int md3D_alloc_init_profil(struct Param_struct *par)
{
	/* Initilisation de variables */
	par->k_super = 2*PI*par->nu_super/par->lambda;
	par->k_sub = 2*PI*par->nu_sub/par->lambda;
	par->Delta_sigma_x = 2*PI/par->Lx;
	par->Delta_sigma_y = 2*PI/par->Ly;
	par->sigma_x0 = par->k_super*sin(par->theta_i)*cos(par->phi_i);
	par->sigma_y0 = par->k_super*sin(par->theta_i)*sin(par->phi_i);
			
	/* Allocations */	
	if (par->profile_type == N_XYZ || par->profile_type == N_XY_ZINVAR) {
		par->n_xyz = allocate_CplxMatrix(par->Nprz,par->Nprx*par->Npry);
	}else if(par->profile_type == H_XY || par->profile_type == MULTICOUCHES){
		par->profil = allocate_DbleMatrix(par->N_layers+3,par->Nprx*par->Npry);
		par->k2_layer    = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->N_layers+2));
		par->invk2_layer = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->N_layers+2));
	}else{
		fprintf(stderr, "%s, line %d : ERROR, unknown profile_type\n",__FILE__,__LINE__);
		exit(EXIT_FAILURE);
	}


	/* Alignement des pointeurs de fonction */
/*	par->matrice_T = (par->pola == TE ? matrice_T_TE : matrice_T_TM);
*/	switch (par->profile_type) {
		case H_XY : 
			par->md3D_lire_profil = md3D_lire_profil_H_XY;
			par->k_2 = k2_H_XY;
			par->invk_2 = invk2_H_XY;
			par->Normal_function = Normal_H_XY;
			break;
/*		case H_XY_plus_STACK : 
			par->md3D_lire_profil = md3D_lire_profil_H_XY;
			par->k_2 = k2_H_XY;
			par->invk_2 = invk2_H_XY;
			par->Normal_function = Normal_H_XY;
			break;*/
		case MULTICOUCHES : 
			par->md3D_lire_profil = md3D_lire_profil_MULTI;
			par->k_2 = k2_MULTI;
			par->invk_2 = invk2_MULTI;
			par->Normal_function = Normal_MULTI;
			break;
		case N_XYZ        : 
			par->md3D_lire_profil = md3D_lire_profil_N_XYZ;
			par->k_2 = k2_N_XYZ;
			par->invk_2 = invk2_N_XYZ;
/*			par->Normal_function = Normal_N_XYZ;*/
printf("\n***************\nWARNING: Normal_N_XYZ not defined yet.\n******************\n\n");
			break;
		case N_XY_ZINVAR : 
			par->md3D_lire_profil = md3D_lire_profil_N_XYZ;
			par->k_2 = k2_N_XYZ;
			par->invk_2 = invk2_N_XYZ;
			par->Normal_function = Normal_N_XY_ZINVAR;
			break;
		default:
			fprintf(stderr, "%s, line %d : ERROR, (%d) is an unknown profile_type\n",__FILE__,__LINE__,par->profile_type);
			exit(EXIT_FAILURE);

	}
	if (par->profile_type == H_XY_plus_STACK){
		printf("Alloc nu_stack et h_stack, N_stack = %d",par->N_stack);fflush(stdout);
		par->nu_stack =  (COMPLEX *) malloc(sizeof(COMPLEX)*(par->N_stack));
		par->h_stack =  (double *) malloc(sizeof(double)*(par->N_stack));
	}	

	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md3D_variables_init(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *	\brief	Initialisations des variables
 */
/*---------------------------------------------------------------------------------------------*/
int md3D_variables_init(struct Param_struct *par, struct Efficacites_struct *eff)
{
	par->clock0 = clock(); /* Initialisation des chronomètres                  */
	time(&(par->time0));   /* clock0 (courte durées) et time0 (longues durées) */
	par->last_clock = clock();
	time(&(par->last_time)); 

	/* Initialisation des variables en mode "AUTO" */
	/* delta_h = lambda x (n_re + n_im) / 1000 */
	if (par->delta_h == AUTO) {
		par->delta_h = par->lambda/(cabs(par->nu_sub)*1000);
	}
	/* NS = (h/lambda) x (n_re + n_im) x 5 */
	if (par->NS == AUTO) {
		par->NS = CEIL( (par->h / par->lambda)*cabs(par->nu_sub)*5 );
	}

	/* vec_size : size of most matrices and vetcors */	
	par->vec_size = (2*par->Nx+1)*(2*par->Ny+1);
	par->vec_middle = 2*par->Nx*par->Ny + par->Nx + par->Ny;

	/* for H_XY profile the normal components need to be calculated only once */
	par->HXY_Normal_CALCULATED = 0;
	par->toepNorm_CALCULATED = 0;
	
	/* Calcul de la valeur exacte de delta_h de sorte qu'il y en ait un nb entier à chaque étape Matrice-S*/
	par->N_steps = 0; /* Counter initialisation */
	int Nstep_S = ROUND(ceil((par->h/par->NS)/par->delta_h));
	par->delta_h = (par->h/par->NS)/Nstep_S;
fprintf(stderr,"\
**************************************\n\
* CAUTION: parameter delta_h is not used.\n\
* The z discretization is only done through NS so far.\n\
* TO DO : implementing P matrix multiplication at each delta_h\n\
* %s, line %d\n\
**************************************\n",__FILE__,__LINE__);
	/*par->Nstep = ROUND(par->h/par->delta_h);*/

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
	if (!strcmp(par->calcul_method,"Z_INVAR")){
			par->M_matrix = zinvar_M_matrix;
	}else{
			par->M_matrix = M_matrix;
	}

	if (!strcmp(par->calcul_method,"RK4")){
			par->P_matrix = rk4_P_matrix;
/*			par->P_matrix = zinvar_P_matrix;printf("CAUTION ! P_matrix = zinvar_P_matrix DEBUGGING...");*/
	}else if (!strcmp(par->calcul_method,"Z_INVAR")){
			par->P_matrix = zinvar_P_matrix;
	}else if (!strcmp(par->calcul_method,"IMPROVED_RCWA")){
			par->P_matrix = zinvar_P_matrix;
	}else{
		fprintf(stderr, "%s, line %d : ERROR, unknown calculation method (\"%s\")\n",__FILE__,__LINE__,par->calcul_method);
		exit(EXIT_FAILURE);
	}


	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md3D_arrays_init(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *		\brief	Some arrays initialisations
 */
/*---------------------------------------------------------------------------------------------*/
int md3D_arrays_init(struct Param_struct *par, struct Efficacites_struct *eff)
{
	int i,j, Nx=par->Nx, Ny=par->Ny;
	
	COMPLEX k_super2 = par->k_super*par->k_super;
	COMPLEX k_sub2 = par->k_sub*par->k_sub;

	for (i=-Ny;i<=Ny;i++){
		for (j=-Nx;j<=Nx;j++){
			par->nx[(i+Ny)*(2*Nx+1)+(j+Nx)] = j;
			par->ny[(i+Ny)*(2*Nx+1)+(j+Nx)] = i;
		}
	}
/*printf("\n nx : \n");
for (i=0;i<par->vec_size;i++){
printf("%d ",par->nx[i]);
}*/
	for(j=0;j<=par->vec_size-1;j++){
		par->sigma_x[j] = par->nx[j]*par->Delta_sigma_x + par->sigma_x0;
		par->sigma_y[j] = par->ny[j]*par->Delta_sigma_y + par->sigma_y0;

		par->kz_super[j] = csqrt(k_super2 - par->sigma_x[j]*par->sigma_x[j] - par->sigma_y[j]*par->sigma_y[j]);
		par->kz_sub[j]   = csqrt(k_sub2   - par->sigma_x[j]*par->sigma_x[j] - par->sigma_y[j]*par->sigma_y[j]);
	}
	PsiMatrix(par->Psi_super, par->k_super, par->kz_super, par->sigma_x, par->sigma_y, par->vec_size);
	PsiMatrix(par->Psi_sub, par->k_sub, par->kz_sub, par->sigma_x, par->sigma_y, par->vec_size);
	invM(par->invPsi_super, par->Psi_super, 4*par->vec_size);
/*printf("\n sigma_x : \n");SaveCplxTab2file (par->sigma_x, par->vec_size, "Re", "stdout", " ",2*par->Nx+1,"\n");
printf("\n sigma_y : \n");SaveCplxTab2file (par->sigma_y, par->vec_size, "Re", "stdout", " ",2*par->Nx+1,"\n");
*/
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
		md3D_make_tab_S_steps(par);
	}
/*SaveDbleTab2file(par->tab_NS, par->NS+1, "stdout", " ");*/
/*printf("\nRe(par->Psi_super) :\n");
SaveMatrix2file (par->Psi_super, 2*par->vec_size, 2*par->vec_size, "Re", "stdout");
printf("\nRe(par->Psi_sub) :\n");
SaveMatrix2file (par->Psi_sub, 2*par->vec_size, 2*par->vec_size, "Re", "stdout");
printf("\nRe(par->invPsi_super) :\n");
SaveMatrix2file (par->invPsi_super, 2*par->vec_size, 2*par->vec_size, "Re", "stdout");
*/		
	return 0;
}


/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md3D_alloc(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *	\brief	Memory allocation
 */
/*---------------------------------------------------------------------------------------------*/
int md3D_alloc(struct Param_struct *par, struct Efficacites_struct *eff)
{
	
	par->nx = (int *) malloc(sizeof(int)*par->vec_size);
	par->ny = (int *) malloc(sizeof(int)*par->vec_size);
	par->sigma_x  = (COMPLEX *) malloc(sizeof(COMPLEX)*par->vec_size);
	par->sigma_y  = (COMPLEX *) malloc(sizeof(COMPLEX)*par->vec_size);
	par->kz_super = (COMPLEX *) malloc(sizeof(COMPLEX)*par->vec_size);
	par->kz_sub   = (COMPLEX *) malloc(sizeof(COMPLEX)*par->vec_size);

	par->S12 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	par->S22 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	par->S21 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	par->S11 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);

	par->Nxx = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Nxy = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Nxz = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Nyy = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Nyz = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Nzz = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Qxx = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Qxy = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Qxz = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Qyy = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Qyz = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Qzz = allocate_CplxMatrix(par->vec_size,par->vec_size);
	par->Qzz_1 = allocate_CplxMatrix(par->vec_size,par->vec_size);

	par->T11 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	par->T12 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	par->T21 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	par->T22 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);

	par->P         = allocate_CplxMatrix(4*par->vec_size,4*par->vec_size);
	par->M         = allocate_CplxMatrix(4*par->vec_size,4*par->vec_size);

	par->Psi_sub      = allocate_CplxMatrix(4*par->vec_size,4*par->vec_size);
	par->Psi_super    = allocate_CplxMatrix(4*par->vec_size,4*par->vec_size);
	par->invPsi_super = allocate_CplxMatrix(4*par->vec_size,4*par->vec_size);

	par->Ai = (COMPLEX *) malloc(sizeof(COMPLEX)*(2*par->vec_size));
	par->Ar = (COMPLEX *) malloc(sizeof(COMPLEX)*(2*par->vec_size));
	par->At = (COMPLEX *) malloc(sizeof(COMPLEX)*(2*par->vec_size));	
	par->Vi = (COMPLEX *) malloc(sizeof(COMPLEX)*(2*par->vec_size));
	par->Vr = (COMPLEX *) malloc(sizeof(COMPLEX)*(2*par->vec_size));
	par->Vt = (COMPLEX *) malloc(sizeof(COMPLEX)*(2*par->vec_size));	

	par->Exr  = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));	
	par->Ext  = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));	
	par->Hpxi = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));	
	par->Hpxr = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));	
	par->Hpxt = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));	
	par->Exi  = (COMPLEX *) malloc(sizeof(COMPLEX)*(par->vec_size));	

	eff->eff_r    = allocate_DbleMatrix(eff->Nymax_super-eff->Nymin_super+1,eff->Nxmax_super-eff->Nxmin_super+1);
	eff->nx_eff_r = allocate_DbleMatrix(eff->Nymax_super-eff->Nymin_super+1,eff->Nxmax_super-eff->Nxmin_super+1);
	eff->ny_eff_r = allocate_DbleMatrix(eff->Nymax_super-eff->Nymin_super+1,eff->Nxmax_super-eff->Nxmin_super+1);
	eff->eff_t    = allocate_DbleMatrix(eff->Nymax_sub-eff->Nymin_sub+1,eff->Nxmax_sub-eff->Nxmin_sub+1);
	eff->nx_eff_t = allocate_DbleMatrix(eff->Nymax_sub-eff->Nymin_sub+1,eff->Nxmax_sub-eff->Nxmin_sub+1);
	eff->ny_eff_t = allocate_DbleMatrix(eff->Nymax_sub-eff->Nymin_sub+1,eff->Nxmax_sub-eff->Nxmin_sub+1);

	if (par->tab_NS_ENABLED){
		par->tab_NS = (double *) malloc(sizeof(double)*(par->NS+1));
	}
	if (par->imposed_S_steps){
		par->tab_imposed_S_steps = (double *) malloc(sizeof(double)*(par->N_imposed_S_steps));
	}

	return 0;
}



/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int md3D_free(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *	\brief	Libération de la mémoire
 */
/*---------------------------------------------------------------------------------------------*/
int md3D_free(struct Param_struct *par, struct Efficacites_struct *eff)
{
	if (par->verbosity>=3) fprintf(stdout,"Freeing memory : "); fflush(stdout);

	if (par->profile_type == H_XY_plus_STACK){
		par->nu_stack =  (COMPLEX *) malloc(sizeof(COMPLEX)*(par->N_stack));
		par->h_stack =  (double *) malloc(sizeof(double)*(par->N_stack));
	}	

	if (par->profile_type == N_XYZ || par->profile_type == N_XY_ZINVAR) {
		free(par->n_xyz[0]);
		free(par->n_xyz);
	}else if(par->profile_type == H_XY || par->profile_type == MULTICOUCHES){
		free(par->profil[0]);
		free(par->profil);
		free(par->k2_layer);
		free(par->invk2_layer);
	}else{
		fprintf(stderr, "%s, line %d : ERROR, unknown profile_type\n",__FILE__,__LINE__);
		exit(EXIT_FAILURE);
	}
	if (par->tab_NS_ENABLED){
		free(par->tab_NS);
	}
	if (par->imposed_S_steps){
		free(par->tab_imposed_S_steps);
	}

	free(par->sigma_x);
	free(par->sigma_y);
	free(par->kz_super);
	free(par->kz_sub);
		
	free(par->S12[0]);
	free(par->S12);
	free(par->S22[0]);
	free(par->S22);
	free(par->S11[0]);
	free(par->S11);
	free(par->S21[0]);
	free(par->S21);

	free(par->T11[0]);
	free(par->T12[0]);
	free(par->T21[0]);
	free(par->T22[0]);
	free(par->T11);
	free(par->T12);
	free(par->T21);
	free(par->T22);
	
	free(par->P[0]);
	free(par->P);
	free(par->M[0]);
	free(par->M);

	free(par->Nxx[0]);
	free(par->Nxy[0]);
	free(par->Nxz[0]);
	free(par->Nyy[0]);
	free(par->Nyz[0]);
	free(par->Nzz[0]);
	free(par->Nxx);
	free(par->Nxy);
	free(par->Nxz);
	free(par->Nyy);
	free(par->Nyz);
	free(par->Nzz);

	free(par->Qxx[0]);
	free(par->Qxy[0]);
	free(par->Qxz[0]);
	free(par->Qyy[0]);
	free(par->Qyz[0]);
	free(par->Qzz[0]);
	free(par->Qzz_1[0]);
	free(par->Qxx);
	free(par->Qxy);
	free(par->Qxz);
	free(par->Qyy);
	free(par->Qyz);
	free(par->Qzz);
	free(par->Qzz_1);
	
	free(par->invPsi_super[0]);
	free(par->invPsi_super);
	free(par->Psi_sub[0]);
	free(par->Psi_sub);
	free(par->Psi_super[0]);
	free(par->Psi_super);
	
	free(par->Ai);
	free(par->Ar);
	free(par->At);
	free(par->Vi);
	free(par->Vr);
	free(par->Vt);
	free(par->Exi);
	free(par->Exr);
	free(par->Hpxi);
	free(par->Hpxr);
	free(par->Hpxt);
	free(par->Ext);

	if (par->verbosity>=3) fprintf(stdout,"OK\n"); fflush(stdout);

	
	return 0;
}


