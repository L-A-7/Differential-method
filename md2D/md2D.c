/*!	\file		md2D.c
 *
 * 	\brief		Differential method \n
 *				2-Dimensions (as opposed to 3-D) \n
 *				S-Matrices algorithm \n
 * 				FFF algorithm \n
 * 
 *
 *	\date		nov 2007
 *	\authors	Laurent ARNAUD
 */


#include "md2D.h"


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_incident_field(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *		\brief	Incident field determination
 */
/*-------------------------------------------------------------------------------------*/
int md2D_incident_field(struct Param_struct *par, struct Efficacites_struct *eff)
{
	int n;

	/* Plane wave */
	if (!strcmp(par->i_field_mode,"PLANE_WAVE")){

		for (n=0; n<=par->vec_size-1; n++) {
			par->Ai[n] = 0;
		}
		par->Ai[par->vec_middle] = 1;

/*		par->Ai[par->vec_middle] = -par->k_super;*/
		 
	}else{
		fprintf(stderr, "%s, ligne %d : ERROR, \"%s\" : unsupported i_field_mode type\n",__FILE__,__LINE__,par->i_field_mode);
		exit(EXIT_FAILURE);
	}
	
	return 0;	
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_propagativ_limits(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *	\brief		Calculation of Propagatives modes limits
 */
/*-------------------------------------------------------------------------------------*/
int md2D_propagativ_limits(struct Param_struct *par, struct Efficacites_struct *eff)
{
	int N = par->N;
	double Delta_sigma = par->Delta_sigma;
	COMPLEX k_super = par->k_super;
	COMPLEX k_sub = par->k_sub;
	double sigma0 = par->sigma0;
	int Nmin_super, Nmax_super, Nmin_sub, Nmax_sub;

	/* Propagativ modes limits */
	Nmax_super =  FLOOR( ( creal(k_super) - sigma0 )/Delta_sigma );
	Nmin_super = -FLOOR( ( creal(k_super) + sigma0 )/Delta_sigma );
	Nmax_sub   =  FLOOR( ( creal(k_sub) - sigma0 )/Delta_sigma );
	Nmin_sub   = -FLOOR( ( creal(k_sub) + sigma0 )/Delta_sigma );
	int Nlimit = MAX(MAX(Nmax_super,Nmax_sub),MAX(-Nmin_super,-Nmin_sub));
	if (N < Nlimit) {
		if(par->verbosity >= 1){
			fprintf(stderr,"WARNING, N too small to represent all propagatives orders\n");
			fprintf(stderr,"minimum value: N = %d \n",Nlimit);
		}
		if (N < Nmax_super)  Nmax_super =  N;
		if (N < Nmax_sub)    Nmax_sub   =  N;
		if (Nmin_super < -N) Nmin_super = -N;
		if (Nmin_sub   < -N) Nmin_sub   = -N;
	}
	eff->Nmin_super = Nmin_super;
	eff->Nmax_super = Nmax_super;
	eff->Nmin_sub = Nmin_sub;
	eff->Nmax_sub = Nmax_sub;

	return 0;	
}

	
/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_efficiencies(COMPLEX *Ai, COMPLEX *A0, COMPLEX *Ah, struct Param_struct *par,  struct Efficacites_struct *eff)
 *
 *	\brief		Efficiencies calculation
 */
/*-------------------------------------------------------------------------------------*/
int md2D_efficiencies(COMPLEX *Ai, COMPLEX *Ar, COMPLEX *At, struct Param_struct *par, struct Efficacites_struct *eff)
{
	int n;
	int vec_size = par->vec_size;
	int mid = par->vec_middle;
	
	COMPLEX k_super = par->k_super;
	COMPLEX k_sub = par->k_sub;
	COMPLEX k_super2 = k_super*k_super;
	COMPLEX k_sub2 = k_sub*k_sub;
/*	COMPLEX ky_0  = par->ky_0;
	COMPLEX ky_02 = ky_0*ky_0;*/
	
	COMPLEX *kz_super, *kz_sub, *sigma;
	kz_super = par->kz_super;
	kz_sub = par->kz_sub;
	sigma = par->sigma;

	double Pz_i, Pz_r, Pz_t, sumPzi;
	COMPLEX *Eyi, *Hpyi, *Eyr, *Hpyr, *Eyt, *Hpyt;
	COMPLEX *Exi, *Hpxi, *Exr, *Hpxr, *Ext, *Hpxt;
	Exi = par->Exi;
	Exr = par->Exr;
	Ext = par->Ext;
	Hpxi = par->Hpxi;
	Hpxr = par->Hpxr;
	Hpxt = par->Hpxt;
	
	int Nmin_super = eff->Nmin_super;
	int Nmax_super = eff->Nmax_super;
	int Nmin_sub   = eff->Nmin_sub;
	int Nmax_sub   = eff->Nmax_sub;

	/* y field components */
	Eyi  = par->Ai; /* TE */
	Eyr  = par->Ar;
	Eyt  = par->At;
	Hpyi = par->Ai; /* TM */
	Hpyr = par->Ar;
	Hpyt = par->At;

	/* x field components calculations */
	if (par->pola == TE){
		for(n=0;n<=vec_size-1;n++){
			Exi[n]  = 0;
			Exr[n]  = 0;
			Ext[n]  = 0;
			Hpxi[n] = kz_super[n]*Eyi[n];
			Hpxr[n] = -kz_super[n]*Eyr[n];
			Hpxt[n] = kz_sub[n]*Eyt[n];
		}
	}else if (par->pola == TM){
		for(n=0;n<=vec_size-1;n++){
			Exi[n]  = -kz_super[n]*Hpyi[n]/k_super2;
			Exr[n]  =  kz_super[n]*Hpyr[n]/k_super2;
			Ext[n]  = -kz_sub[n]*Hpyt[n]/k_sub2;
			Hpxi[n] = 0;
			Hpxr[n] = 0;
			Hpxt[n] = 0;
		}
	}else{
		fprintf(stderr, "%s, line %d: ERROR, pola = %d, must be \"TE\" or \"TM\".Exiting\n",__FILE__,__LINE__,par->pola);
		exit(EXIT_FAILURE);
	}

	/* Total incident energie */
	sumPzi = 0;
	for (n=Nmin_super; n<=Nmax_super; n++) {
		Pz_i = fabs(creal(Exi[n+mid]*CONJ(Hpyi[n+mid]) - Eyi[n+mid]*CONJ(Hpxi[n+mid])));
		sumPzi += Pz_i;
	}
/*sumPzi=1;
fprintf(stderr,"WARNING, sumPzi=%f\n",sumPzi);*/

	/* Reflexion : efficiencies and directions */
	for (n=Nmin_super; n<=Nmax_super; n++) {
		Pz_r = fabs(creal(Exr[n+mid]*CONJ(Hpyr[n+mid]) - Eyr[n+mid]*CONJ(Hpxr[n+mid])));
		eff->eff_r[n-Nmin_super] = Pz_r/sumPzi;
		eff->N_eff_r[n-Nmin_super] = (double) n;
		eff->theta_r[n-Nmin_super] =  SIGN(sigma[n+mid])*acos(kz_super[n+mid]/k_super) * 180.0/PI;
	}

	/* Transmission : efficiencies and directions */
	for (n=Nmin_sub; n<=Nmax_sub; n++) {
		Pz_t = fabs(creal(Ext[n+mid]*CONJ(Hpyt[n+mid]) - Eyt[n+mid]*CONJ(Hpxt[n+mid])));
		eff->eff_t[n-Nmin_sub] = Pz_t/sumPzi;
		eff->N_eff_t[n-Nmin_sub] = (double) n;
		eff->theta_t[n-Nmin_sub] = SIGN(sigma[n+mid])*acos(kz_sub[n+mid]/k_sub) * 180.0/PI;
	}

	/* Sum of efficiencies */
	eff->sum_eff_r = 0;
	eff->sum_eff_t = 0;
	for (n=Nmin_super; n<=Nmax_super; n++) {
		eff->sum_eff_r += eff->eff_r[n-Nmin_super];
	}
	for (n=Nmin_sub; n<=Nmax_sub; n++) {
		eff->sum_eff_t += eff->eff_t[n-Nmin_sub];
	}
	eff->sum_eff = eff->sum_eff_r + eff->sum_eff_t;

	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_efficiencies(COMPLEX *Ai, COMPLEX *A0, COMPLEX *Ah, struct Param_struct *par,  struct Efficacites_struct *eff)
 *
 *	\brief		Special function used for calculation of scatterd flux by scatterer on waveguide (energy flux are stored in eff->eff_r or eff->eff_t for minimizing the code modifications, but are not efficiencies)
 */
/*-------------------------------------------------------------------------------------*/
int md2D_energy_flux(struct Param_struct *par, struct Efficacites_struct *eff)
{
	int n;
	int vec_size = par->vec_size;
	int mid = par->vec_middle;
	
	COMPLEX k_super = par->k_super;
	COMPLEX k_sub = par->k_sub;
	COMPLEX k_super2 = k_super*k_super;
	COMPLEX k_sub2 = k_sub*k_sub;
	
	COMPLEX *At, *Ar, *kz_super, *kz_sub, *sigma, Coeff_r, Coeff_t;
	kz_super = par->kz_super;
	kz_sub = par->kz_sub;
	sigma = par->sigma;
	At = par->At;
	Ar = par->Ar;
	
	int Nmin_super = eff->Nmin_super;
	int Nmax_super = eff->Nmax_super;
	int Nmin_sub   = eff->Nmin_sub;
	int Nmax_sub   = eff->Nmax_sub;

	/* Vi = [0 0 ... 0 1 0 ... 0 0] */
	for (n=0;n<=vec_size-1;n++){
		par->Vi[n] = 0;
	}
	par->Vi[par->vec_middle] = 1;
	/* Vt = S22*Vi */
	M_x_V (par->Vt, par->S22, par->Vi, vec_size, vec_size);

	/* Vr = S12*Vi */	
	M_x_V (par->Vr, par->S12, par->Vi, vec_size, vec_size);

	/* Ar = Vr*cexp(-I*kz_super*h) */
	for (n=0;n<=vec_size-1;n++){
		Ar[n] = par->Vr[n]*cexp(-I*kz_super[n]*par->h);
	}
	/* At = Vt */
	for (n=0;n<=vec_size-1;n++){
		At[n] = par->Vt[n];
	}
printf("\nVt2_mid= %f\n", cabs(par->Vt[mid]*conj(par->Vt[mid])));
printf("Vr2_mid= %f\n", cabs(par->Vr[mid]*conj(par->Vr[mid])));
printf("At2_mid= %f\n", cabs(At[mid]*conj(At[mid])));
printf("Ar2_mid= %f\n", cabs(Ar[mid]*conj(Ar[mid])));
/*printf("\nRe(Ar)\n");SaveCplxTab2file (Ar, vec_size, "Re", "stdout", " ", 10000, "\n");
printf("\nIm(Ar)\n");SaveCplxTab2file (Ar, vec_size, "Im", "stdout", " ", 10000, "\n");
printf("\nRe(At)\n");SaveCplxTab2file (At, vec_size, "Re", "stdout", " ", 10000, "\n");
printf("\nIm(At)\n");SaveCplxTab2file (At, vec_size, "Im", "stdout", " ", 10000, "\n");*/

	if (par->pola == TE){
		Coeff_r=1;
		Coeff_t=1;
	}else if (par->pola == TM){
		Coeff_r=1/k_super2;
		Coeff_t=1/k_sub2;
	}else{
		fprintf(stderr, "%s, line %d: ERROR, pola = %d, must be \"TE\" or \"TM\".Exiting\n",__FILE__,__LINE__,par->pola);
		exit(EXIT_FAILURE);
	}

	/* Reflexion */
	for (n=Nmin_super; n<=Nmax_super; n++) {
		eff->eff_r[n-Nmin_super] = Coeff_r*0.5*Ar[n+mid]*conj(Ar[n+mid])*creal(kz_super[n+mid])*par->L;
		eff->N_eff_r[n-Nmin_super] = (double) n;
		eff->theta_r[n-Nmin_super] =  SIGN(sigma[n+mid])*acos(kz_super[n+mid]/k_super) * 180.0/PI;
	}

	/* Transmission */
	for (n=Nmin_sub; n<=Nmax_sub; n++) {
		eff->eff_t[n-Nmin_sub] = Coeff_t*0.5*At[n+mid]*conj(At[n+mid])*creal(kz_sub[n+mid])*par->L;
		eff->N_eff_t[n-Nmin_sub] = (double) n;
		eff->theta_t[n-Nmin_sub] = SIGN(sigma[n+mid])*acos(kz_sub[n+mid]/k_sub) * 180.0/PI;
	}

	/* Sum of efficiencies */
	eff->sum_eff_r = 0;
	eff->sum_eff_t = 0;
	for (n=Nmin_super; n<=Nmax_super; n++) {
		eff->sum_eff_r += eff->eff_r[n-Nmin_super];
	}
	for (n=Nmin_sub; n<=Nmax_sub; n++) {
		eff->sum_eff_t += eff->eff_t[n-Nmin_sub];
	}
	eff->sum_eff = eff->sum_eff_r + eff->sum_eff_t;
printf("totalFlux_r=%f\n", eff->sum_eff_r);
printf("totalFlux_t=%f\n", eff->sum_eff_t);

	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int swifts_md2D_energy_flux(double h_guide, struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *	\brief		Special function used for calculation of scatterd flux by scatterer on waveguide (energy flux are stored in eff->eff_r or eff->eff_t for minimizing the code modifications, but are not efficiencies)
 */
/*-------------------------------------------------------------------------------------*/
int swifts_md2D_energy_flux(double h_guide, struct Param_struct *par, struct Efficacites_struct *eff)
{
	int n;
	int vec_size = par->vec_size;
	int mid = par->vec_middle;
	
	COMPLEX k_super = par->k_super;
	COMPLEX k_sub = par->k_sub;
	COMPLEX k_super2 = k_super*k_super;
	COMPLEX k_sub2 = k_sub*k_sub;
	
	COMPLEX *At, *Ar, *kz_super, *kz_sub, *sigma, Coeff_r, Coeff_t;
	kz_super = par->kz_super;
	kz_sub = par->kz_sub;
	sigma = par->sigma;
	At = par->At;
	Ar = par->Ar;
	
	int Nmin_super = eff->Nmin_super;
	int Nmax_super = eff->Nmax_super;
	int Nmin_sub   = eff->Nmin_sub;
	int Nmax_sub   = eff->Nmax_sub;

	if (!strcmp(par->S_matrix_blocks_calculation,"ALL_4_S_MATRIX_BLOCKS")){ /* Lightening from below */
		/* Vi = [0 0 ... 0 1 0 ... 0 0] */
		for (n=0;n<=vec_size-1;n++){
			par->Vi[n] = 0;
		}
		par->Vi[par->vec_middle] = 1;
		/* Vt = S21*V0p */
		M_x_V (par->Vt, par->S21, par->Vi, vec_size, vec_size);
		/* Vr = S11*V0p */	
		M_x_V (par->Vr, par->S11, par->Vi, vec_size, vec_size);
	}else{ /* Lightening from top */
		/* Vi = [0 0 ... 0 1 0 ... 0 0] */
		for (n=0;n<=vec_size-1;n++){
			par->Vi[n] = 0;
		}
		par->Vi[par->vec_middle] = 1;
		/* Vt = S22*Vi */
		M_x_V (par->Vt, par->S22, par->Vi, vec_size, vec_size);
		/* Vr = S12*Vi */	
		M_x_V (par->Vr, par->S12, par->Vi, vec_size, vec_size);
	}
	/* Ar = Vr*cexp(-I*kz_super*h) */
	for (n=0;n<=vec_size-1;n++){
		Ar[n] = par->Vr[n]*cexp(-I*kz_super[n]*(par->h+h_guide));
	}
	/* At = Vt */
	for (n=0;n<=vec_size-1;n++){
/*		At[n] = par->Vt[n]*cexp(I*kz_super[n]*h_guide);
*/		At[n] = par->Vt[n];
	}

printf("\nVt2_mid= %f\n", cabs(par->Vt[mid]*conj(par->Vt[mid])));
printf("Vr2_mid= %f\n", cabs(par->Vr[mid]*conj(par->Vr[mid])));
printf("At2_mid= %f\n", cabs(At[mid]*conj(At[mid])));
printf("Ar2_mid= %f\n", cabs(Ar[mid]*conj(Ar[mid])));
/*printf("\nRe(Ar)\n");SaveCplxTab2file (Ar, vec_size, "Re", "stdout", " ", 10000, "\n");
printf("\nIm(Ar)\n");SaveCplxTab2file (Ar, vec_size, "Im", "stdout", " ", 10000, "\n");
printf("\nRe(At)\n");SaveCplxTab2file (At, vec_size, "Re", "stdout", " ", 10000, "\n");
printf("\nIm(At)\n");SaveCplxTab2file (At, vec_size, "Im", "stdout", " ", 10000, "\n");*/

	if (par->pola == TE){
		Coeff_r=1;
		Coeff_t=1;
	}else if (par->pola == TM){
		Coeff_r=1/k_super2;
		Coeff_t=1/k_sub2;
	}else{
		fprintf(stderr, "%s, line %d: ERROR, pola = %d, must be \"TE\" or \"TM\".Exiting\n",__FILE__,__LINE__,par->pola);
		exit(EXIT_FAILURE);
	}

	/* Reflexion */
	for (n=Nmin_super; n<=Nmax_super; n++) {
		eff->eff_r[n-Nmin_super] = Coeff_r*0.5*Ar[n+mid]*conj(Ar[n+mid])*creal(kz_super[n+mid])*par->L;
		eff->N_eff_r[n-Nmin_super] = (double) n;
		eff->theta_r[n-Nmin_super] =  SIGN(sigma[n+mid])*acos(kz_super[n+mid]/k_super) * 180.0/PI;
	}

	/* Transmission */
	for (n=Nmin_sub; n<=Nmax_sub; n++) {
		eff->eff_t[n-Nmin_sub] = Coeff_t*0.5*At[n+mid]*conj(At[n+mid])*creal(kz_sub[n+mid])*par->L;
		eff->N_eff_t[n-Nmin_sub] = (double) n;
		eff->theta_t[n-Nmin_sub] = SIGN(sigma[n+mid])*acos(kz_sub[n+mid]/k_sub) * 180.0/PI;
	}

	/* Sum of efficiencies */
	eff->sum_eff_r = 0;
	eff->sum_eff_t = 0;
	for (n=Nmin_super; n<=Nmax_super; n++) {
		eff->sum_eff_r += eff->eff_r[n-Nmin_super];
	}
	for (n=Nmin_sub; n<=Nmax_sub; n++) {
		eff->sum_eff_t += eff->eff_t[n-Nmin_sub];
	}
	eff->sum_eff = eff->sum_eff_r + eff->sum_eff_t;
printf("totalFlux_r=%f\n", eff->sum_eff_r);
printf("totalFlux_t=%f\n", eff->sum_eff_t);

	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_amplitudes(COMPLEX *Ai, COMPLEX *Ar, COMPLEX *At, COMPLEX **S12, 
 *                               COMPLEX **S22, struct Param_struct *par)
 *
 *	\brief	Field amplitudes calculations
 */
/*-------------------------------------------------------------------------------------*/
int md2D_amplitudes(COMPLEX *Ai, COMPLEX *Ar, COMPLEX *At, COMPLEX **S12, COMPLEX **S22, struct Param_struct *par)
{
	int n;
	int vec_size = par->vec_size;
	COMPLEX *kz_super;
	kz_super = par->kz_super;
		
	/* Vi = Ai*cexp(-I*kz_super*h) */
	for (n=0;n<=vec_size-1;n++){
		par->Vi[n] = Ai[n]*cexp(-I*kz_super[n]*par->h);
/*		par->Vi[n+vec_size] = Ai[n+vec_size]*cexp(-I*kz_super[n]*par->h);*/
	}
	/* Vt = S22*Vi */
	M_x_V (par->Vt, S22, par->Vi, vec_size, vec_size);

	/* Vr = S12*Vi */	
	M_x_V (par->Vr, S12, par->Vi, vec_size, vec_size);

	/* Ar = Vr*cexp(-I*kz_super*h) */
	for (n=0;n<=vec_size-1;n++){
		Ar[n] = par->Vr[n]*cexp(-I*kz_super[n]*par->h);
/*		Ar[n+vec_size] = par->Vr[n+vec_size]*cexp(-I*kz_super[n]*par->h);*/
	}
	
	/* At = Vt */
	/*for (n=0;n<=2*vec_size-1;n++){*/
	for (n=0;n<=vec_size-1;n++){
		At[n] = par->Vt[n];
	}
/*printf("Ai\n");SaveCplxTab2file (par->Ai, par->vec_size, "Re", "stdout", " ", 100, "\n");
printf("Ar\n");SaveCplxTab2file (par->Ar, par->vec_size, "Re", "stdout", " ", 100, "\n");
printf("At\n");SaveCplxTab2file (par->At, par->vec_size, "Re", "stdout", " ", 100, "\n");*/

	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn	S_matrix(struct Param_struct *par)
 *
 *		\brief	S-matrix
 */
/*-------------------------------------------------------------------------------------*/
int S_matrix(struct Param_struct *par)
{

	if (par->verbosity) fprintf(stdout,"S-Matrix calculation\n");
	int nS;
	int vec_size = par->vec_size;
	int NS = par->NS;
	COMPLEX **S12, **S22, **S21, **S11, **T11, **T12, **T21, **T22, **Z, **T_tmp, **T_tmp2, **T_tmp3;
	S12 = par->S12;
	S22 = par->S22;
	S21 = par->S21;
	S11 = par->S11;
	T11 = par->T11;
	T12 = par->T12;
	T21 = par->T21;
	T22 = par->T22;

	Z      = allocate_CplxMatrix(par->vec_size,par->vec_size);
	T_tmp  = allocate_CplxMatrix(par->vec_size,par->vec_size);
	T_tmp2 = allocate_CplxMatrix(par->vec_size,par->vec_size);
	T_tmp3 = allocate_CplxMatrix(par->vec_size,par->vec_size);
	
	/* Initialisations */
	S_matrix_init(par);
	
	/* Iterations */
	for (nS=0; nS<=NS-1; nS++) {

		/* T-Matrix calculation */
		T_Matrix(T11, T12, T21, T22, nS, par);

		/* Z = inv(T11 + T12*S12) */
		invM(Z, add_M(T_tmp2, 
			T11, M_x_M(T_tmp,
				T12,S12,	vec_size, vec_size), vec_size, vec_size), vec_size);
		/* S12 = (T21 +T22*S12)*Z */
		M_x_M(S12,
			add_M(T_tmp2, T21, M_x_M(T_tmp,
					T22,S12, vec_size, vec_size), vec_size, vec_size),
			Z, vec_size, vec_size);
		/* S22 = S22*Z */
		M_equals(T_tmp,S22, vec_size, vec_size);
		M_x_M(S22,T_tmp,Z, vec_size, vec_size);

		/* Following blocks only usefull when some light is coming from below [the order (S22, S12),S21,S11 must not be changed] */
		if (!strcmp(par->S_matrix_blocks_calculation,"ALL_4_S_MATRIX_BLOCKS")){
			/* S21 = S21 - S22*T12*S11 */
			sub_M(S21, S21, M_x_M(T_tmp, 
				S22, 
				M_x_M(T_tmp2, T12, S11, vec_size, vec_size), vec_size, vec_size), vec_size, vec_size);
	
			/* S11 = (T22 - S12*T12)*S11 */
			M_equals(T_tmp, S11, vec_size, vec_size);
			M_x_M(S11,
				sub_M(T_tmp2, T22, M_x_M(T_tmp3,
						S12,T12, vec_size, vec_size), vec_size, vec_size),
				T_tmp, vec_size, vec_size);
		}

		/* if NEAR_FIELD, field components are saved in the Near_field_matrix */
		if (!strcmp(par->calcul_type,"NEAR_FIELD")){
			md2D_save_near_field(S12, Z, vec_size, nS, par);
/*			md2D_save_T(T11, T12, T21, T22, vec_size, nS, par);*/ /* Optional but seem to work better for guided modes than more sophisticated local field retriaval algorithm */
		}
		
	}
	
	free(Z[0]);
	free(Z);
	free(T_tmp[0]);
	free(T_tmp);
	free(T_tmp2[0]);
	free(T_tmp2);
	free(T_tmp3[0]);
	free(T_tmp3);
	
	if(par->verbosity >0) {fprintf(stdout,"\n");}
	
	return 0;
}

#if 0
/*-------------------------------------------------------------------------------------*/
/*!	\fn	int S_matrix(struct Param_struct *par)
 *
 *	\brief	S matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int S_matrix(struct Param_struct *par)
{

	if (par->verbosity) fprintf(stdout,"S-Matrix calculation\n");
	int i, j, nS;
	int vec_size = par->vec_size;
	int NS = par->NS;
	COMPLEX **S12, **S22, **T11, **T12, **T21, **T22, **Z, **T_tmp, **T_tmp2;
	S12 = par->S12;
	S22 = par->S22;
	T11 = par->T11;
	T12 = par->T12;
	T21 = par->T21;
	T22 = par->T22;

	Z      = allocate_CplxMatrix(par->vec_size,par->vec_size);
	T_tmp  = allocate_CplxMatrix(par->vec_size,par->vec_size);
	T_tmp2 = allocate_CplxMatrix(par->vec_size,par->vec_size);
	
	/* Initialisations */
	for (i=0; i<=vec_size-1; i++) {
		for (j=0; j<=vec_size-1; j++) {
			S12[i][j] = 0;
			S22[i][j] = 0;
		}
		S22[i][i] = 1;
	}

	/* Iterations */
	for (nS=0; nS<=NS-1; nS++) {

		/* T-Matrix calculation */
		T_Matrix(T11, T12, T21, T22, nS, par);

		/* Z = inv(T11 + T12*S12) */
		invM(Z, add_M(T_tmp2, 
			T11, M_x_M(T_tmp,
				T12,S12,	vec_size, vec_size), vec_size, vec_size), vec_size);
		/* S12 = (T21 +T22*S12)*Z */
		M_x_M(S12,
			add_M(T_tmp2, T21, M_x_M(T_tmp,
					T22,S12, vec_size, vec_size), vec_size, vec_size),
			Z, vec_size, vec_size);
		/* S22 = S22*Z */
		M_equals(T_tmp,S22, vec_size, vec_size);
		M_x_M(S22,T_tmp,Z, vec_size, vec_size);

		/* if NEAR_FIELD, field components are saved in the Near_field_matrix */
		if (!strcmp(par->calcul_type,"NEAR_FIELD")){
			md2D_save_near_field(S12, Z, vec_size, nS, par);
/*			md2D_save_T(T11, T12, T21, T22, vec_size, nS, par);*/ /* Optional but seem to work better for guided modes than more sophisticated local field retriaval algorithm */
		}
		
	}
	
	free(Z[0]);
	free(Z);
	free(T_tmp[0]);
	free(T_tmp);
	free(T_tmp2[0]);
	free(T_tmp2);
	
	if(par->verbosity >0) {fprintf(stdout,"\n");}
	
	return 0;
}
#endif

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int T_Matrix(COMPLEX **T11, COMPLEX **T12, COMPLEX **T21, COMPLEX **T22, int nS, struct Param_struct *par)
 *
 *	\brief T matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int T_Matrix(COMPLEX **T11, COMPLEX **T12, COMPLEX **T21, COMPLEX **T22, int nS, struct Param_struct *par)
{
	int vec_size = par->vec_size;
	int i, j, v = vec_size;
	int pola = par->pola;
	COMPLEX *kz_sub, *kz_super, k_sub, k_super;
	kz_sub = par->kz_sub;
	kz_super = par->kz_super;
	k_sub = par->k_sub;
	k_super = par->k_super;

	double hmin, hmax, z, Delta_z;
	COMPLEX **P, **PPsi11, **PPsi12, **PPsi21, **PPsi22, *Psi22, *Psi11, *Psi12, *Psi21, *iPsi11, *iPsi12, *iPsi21, *iPsi22;

	Psi11 = malloc(sizeof(COMPLEX)*vec_size); Psi12 = malloc(sizeof(COMPLEX)*vec_size);
	Psi21 = malloc(sizeof(COMPLEX)*vec_size); Psi22 = malloc(sizeof(COMPLEX)*vec_size);
	iPsi11 = malloc(sizeof(COMPLEX)*vec_size); iPsi12 = malloc(sizeof(COMPLEX)*vec_size);
	iPsi21 = malloc(sizeof(COMPLEX)*vec_size); iPsi22 = malloc(sizeof(COMPLEX)*vec_size);
	PPsi11 = allocate_CplxMatrix(vec_size,vec_size); PPsi12 = allocate_CplxMatrix(vec_size,vec_size);
	PPsi21 = allocate_CplxMatrix(vec_size,vec_size); PPsi22 = allocate_CplxMatrix(vec_size,vec_size);
	P = allocate_CplxMatrix(2*vec_size,2*vec_size);

	/* Psi matrix */
	if (nS==0){ /* 1st S-Matrix iteration : we are in the substrate */
		PsiMatrix(Psi11, Psi12, Psi21, Psi22, kz_sub, k_sub, vec_size, pola);
	}else{     /* Following iterations : we are in the superstrate */
		PsiMatrix(Psi11, Psi12, Psi21, Psi22, kz_super, k_super, vec_size, pola);
	}
	invPsiMatrix(iPsi11, iPsi12, iPsi21, iPsi22, kz_super, k_super, vec_size, pola);

	/* z of the considered T matrix slice */
	hmin = par->h*nS/par->NS;
	hmax = par->h*(nS+1)/par->NS;
	if (par->tab_NS_ENABLED){
		hmin = par->tab_NS[nS];
		hmax = par->tab_NS[nS+1];
	}
	Delta_z = hmax - hmin;
	z = (hmax + hmin)/2;
/*	z = hmin;*/

	/* P_matrix calculation */
	(*par->P_matrix)(P, z, Delta_z, par);

	/* T matrix : T = inv(Psi_super) * P_matrix * Psi */
/*	M_x_M(par->T,
			invPsi_super, M_x_M(par->M_buffer_2vecsize, par->P, Psi, 2*vec_size, 2*vec_size), 2*vec_size, 2*vec_size);
*/

	/* PPsi = P * Psi */
	for (i=0;i<=vec_size-1;i++){
		for (j=0;j<=vec_size-1;j++){
			/* PPsi11 = P11 * Psi11 + P12 * Psi21 */
			PPsi11[i][j] = P[i][j]*Psi11[j] + P[i][j+v]*Psi21[j];
			/* PPsi12 = P11 * Psi12 + P12 * Psi22 */
			PPsi12[i][j] = P[i][j]*Psi12[j] + P[i][j+v]*Psi22[j];
			/* PPsi21 = P21 * Psi11 + P22 * Psi21 */
			PPsi21[i][j] = P[i+v][j]*Psi11[j] + P[i+v][j+v]*Psi21[j];
			/* PPsi22 = P21 * Psi12 + P22 * Psi22 */
			PPsi22[i][j] = P[i+v][j]*Psi12[j] + P[i+v][j+v]*Psi22[j];
		}
	}

	/* T = inv(Psi_super) * P * Psi */
	for (i=0;i<=vec_size-1;i++){
		for (j=0;j<=vec_size-1;j++){
			/* T11 = iPsi11 * PPsi11 + iPsi12 * PPsi21 */
			T11[i][j] = iPsi11[i]*PPsi11[i][j] + iPsi12[i]*PPsi21[i][j];
			/* PPsi12 = iPsi11 * PPsi12 + iPsi12 * PPsi22 */
			T12[i][j] = iPsi11[i]*PPsi12[i][j] + iPsi12[i]*PPsi22[i][j];
			/* PPsi21 = iPsi21 * PPsi11 + iPsi22 * PPsi21 */
			T21[i][j] = iPsi21[i]*PPsi11[i][j] + iPsi22[i]*PPsi21[i][j];
			/* PPsi22 = iPsi21 * PPsi12 + iPsi22 * PPsi22 */
			T22[i][j] = iPsi21[i]*PPsi12[i][j] + iPsi22[i]*PPsi22[i][j];
		}
	}
/*printf("\nRe(par->T) :\n");
SaveMatrix2file (par->T, 2*par->vec_size, 2*par->vec_size, "Re", "stdout");
printf("\nIm(par->T) :\n");
SaveMatrix2file (par->T, 2*par->vec_size, 2*par->vec_size, "Im", "stdout");*/

	/*Affichage du temps restant à l'écran */
	md2D_affichTemps(par->N,par->N,nS,par->NS,par->ni,par->Ni,par);
	
	free(Psi11);free(Psi12);free(Psi21);free(Psi22);
	free(iPsi11);free(iPsi12);free(iPsi21);free(iPsi22);
	free(PPsi11[0]);free(PPsi11);free(PPsi12[0]);free(PPsi12);
	free(PPsi21[0]);free(PPsi21);free(PPsi22[0]);free(PPsi22);
	free(P[0]);free(P);

	return 0;
}

#if 0
/*-------------------------------------------------------------------------------------*/
/*!	\fn	int swifts_S_matrix(struct Param_struct *par)
 *
 *	\brief	S matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int swifts_S_matrix(COMPLEX n_guide, double h_guide, struct Param_struct *par)
{

	if (par->verbosity) fprintf(stdout,"swifts S-Matrix calculation\n");
	int nS;
	int vec_size = par->vec_size;
	int NS = par->NS;
	COMPLEX **S12, **S22, **T11, **T12, **T21, **T22, **Z, **T_tmp, **T_tmp2;
	S12 = par->S12;
	S22 = par->S22;
	T11 = par->T11;
	T12 = par->T12;
	T21 = par->T21;
	T22 = par->T22;

	Z      = allocate_CplxMatrix(par->vec_size,par->vec_size);
	T_tmp  = allocate_CplxMatrix(par->vec_size,par->vec_size);
	T_tmp2 = allocate_CplxMatrix(par->vec_size,par->vec_size);
	
	/* Initialisations */
	swifts_monolayer_S_matrix(n_guide,h_guide,par);

	/* Iterations */
	for (nS=0; nS<=NS-1; nS++) {

		/* T-Matrix calculation */
		T_Matrix(T11, T12, T21, T22, nS, par);

		/* Z = inv(T11 + T12*S12) */
		invM(Z, add_M(T_tmp2, 
			T11, M_x_M(T_tmp,
				T12,S12,	vec_size, vec_size), vec_size, vec_size), vec_size);
		/* S12 = (T21 +T22*S12)*Z */
		M_x_M(S12,
			add_M(T_tmp2, T21, M_x_M(T_tmp,
					T22,S12, vec_size, vec_size), vec_size, vec_size),
			Z, vec_size, vec_size);
		/* S22 = S22*Z */
		M_equals(T_tmp,S22, vec_size, vec_size);
		M_x_M(S22,T_tmp,Z, vec_size, vec_size);

		/* if NEAR_FIELD, field components are saved in the Near_field_matrix */
		if (!strcmp(par->calcul_type,"NEAR_FIELD")){
			md2D_save_near_field(S12, Z, vec_size, nS, par);
/*			md2D_save_T(T11, T12, T21, T22, vec_size, nS, par);*/ /* Optional but seem to work better for guided modes than more sophisticated local field retriaval algorithm */
		}
	}
	
	free(Z[0]);
	free(Z);
	free(T_tmp[0]);
	free(T_tmp);
	free(T_tmp2[0]);
	free(T_tmp2);
	
	if(par->verbosity >0) {fprintf(stdout,"\n");}
	
	return 0;
}
#endif


/*-------------------------------------------------------------------------------------*/
/*!	\fn	int S_matrix_init(struct Param_struct *par)
 *
 *		\brief	Calculation of S matrix initial value
 */
/*-------------------------------------------------------------------------------------*/
int S_matrix_init(struct Param_struct *par)
{
	int i,j;
	int vec_size=par->vec_size;

	if (!strcmp(par->calculation_on_layer,"ON_1_LAYER")){
		monolayer_S_matrix(par->nu_layer, par->h_layer, par);
	}else{
		/* S12 and S22 */
		for (i=0; i<=vec_size-1; i++) {
			for (j=0; j<=vec_size-1; j++) {
				par->S12[i][j] = 0;
				par->S22[i][j] = 0;
			}
			par->S22[i][i] = 1;
		}
		/* S21 and S11 */
		if (!strcmp(par->S_matrix_blocks_calculation,"ALL_4_S_MATRIX_BLOCKS")){
			for (i=0; i<=vec_size-1; i++) {
				for (j=0; j<=vec_size-1; j++) {
					par->S11[i][j] = 0;
					par->S21[i][j] = 0;
				}
				par->S11[i][i] = 1;
			}
		}
	}

	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn	int monolayer_S_matrix(struct Param_struct *par)
 *
 *		\brief	S matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int monolayer_S_matrix(COMPLEX n_guide, double h_guide, struct Param_struct *par)
{
	int i, j;
	int vec_size = par->vec_size;
	COMPLEX **S12, **S22, **S21, **S11;
	COMPLEX g0,g1,g2,ep,em,A,B,C,D,E,F,G,H,k1_2;
	S12 = par->S12;
	S22 = par->S22;
	S21 = par->S21;
	S11 = par->S11;
	COMPLEX *kz_sub = par->kz_sub;
	COMPLEX *sigma = par->sigma;
	k1_2=((2*PI/par->lambda)*n_guide)*((2*PI/par->lambda)*n_guide);
	/* S12 and S22 */
	for (i=0; i<=vec_size-1; i++) {
		for (j=0; j<=vec_size-1; j++) {
			S12[i][j] = 0;
			S22[i][j] = 0;
		}
		g0=kz_sub[i];
		g1=csqrt(k1_2-sigma[i]*sigma[i]);
		g2=g0; /* Because T algorithm consider the substrate to be the medium for the first iteration */
/*printf("\nn_guide: %f+i%f\n",creal(n_guide),cimag(n_guide));
printf("\nRe{g0, g1, g2}: %f,%f,%f\n",creal(g0),creal(g1),creal(g2));
printf("\nIm{g0, g1, g2}: %f,%f,%f\n",cimag(g0),cimag(g1),cimag(g2));*/
		ep=cexp(I*g1*h_guide);
		em=cexp(-I*g1*h_guide);
		/* S12 */
		A=(g1-g0)/(g1+g0)*ep + em;
		B=g1*(g1-g0)/(g1+g0)*ep -g1*em;
		S12[i][i] = (g2*A+B)/(g2*A-B);
		/* S22 */
		C=(g1+g2)/(2*g1);
		D=(g1-g2)/(2*g1);
		E=D*em+C*ep;
		F=C*em+D*ep;
		G=D*em-C*ep;
		H=C*em-D*ep;
		S22[i][i] = (H*E-F*G)/(H+g0*F/g1);
	}
	/* S21 and S11 */
	if (!strcmp(par->S_matrix_blocks_calculation,"ALL_4_S_MATRIX_BLOCKS")){
		for (i=0; i<=vec_size-1; i++) {
			for (j=0; j<=vec_size-1; j++) {
				S11[i][j] = 0;
				S21[i][j] = 0;
			}
			g0=kz_sub[i];
			g1=csqrt(k1_2-sigma[i]*sigma[i]);
			g2=g0; /* Because T algorithm consider the substrate to be the medium for the first iteration */
			em=cexp(I*g1*h_guide); /* Only difference with S12 & S22 is here */
			ep=cexp(-I*g1*h_guide);
			/* S11 */
			A=(g1-g0)/(g1+g0)*ep + em;
			B=g1*(g1-g0)/(g1+g0)*ep -g1*em;
			S11[i][i] = (g2*A+B)/(g2*A-B);
			/* S21 */
			C=(g1+g2)/(2*g1);
			D=(g1-g2)/(2*g1);
			E=D*em+C*ep;
			F=C*em+D*ep;
			G=D*em-C*ep;
			H=C*em-D*ep;
			S21[i][i] = (H*E-F*G)/(H+g0*F/g1);
		}
	}


/*printf("\nswifts_S_guide:\n");
printf("\nRe(S12)\n");SaveMatrix2file (par->S12, vec_size, vec_size, "Re", "stdout");
printf("\nIm(S12)\n");SaveMatrix2file (par->S12, vec_size, vec_size, "Im", "stdout");
printf("\nRe(S22)\n");SaveMatrix2file (par->S22, vec_size, vec_size, "Re", "stdout");
printf("\nIm(S22)\n");SaveMatrix2file (par->S22, vec_size, vec_size, "Im", "stdout");*/


	
	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int PsiMatrix(COMPLEX *Psi11, COMPLEX *Psi12, COMPLEX *Psi21, COMPLEX *Psi22, COMPLEX *kz, COMPLEX k, int vec_size, int pola)
 *
 *		\brief	Psi block diagonal matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int PsiMatrix(COMPLEX *Psi11, COMPLEX *Psi12, COMPLEX *Psi21, COMPLEX *Psi22, COMPLEX *kz, COMPLEX k, int vec_size, int pola)
{
	int i;

	if (pola == TE){
		for(i=0;i<=vec_size-1;i++){
			Psi11[i] =  1;
			Psi12[i] =  1;
			Psi21[i] =  kz[i];
			Psi22[i] = -kz[i];
		}
	}else if (pola == TM){
		for(i=0;i<=vec_size-1;i++){
			Psi11[i] =  -kz[i]/(k*k);
			Psi12[i] =  kz[i]/(k*k);
			Psi21[i] = 1;
			Psi22[i] = 1;
		}
	}else{
		fprintf(stderr, "%s, line %d : ERROR, unknown polarization... (\"%d\")\n",__FILE__,__LINE__,pola);
		exit(EXIT_FAILURE);
	}

	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int invPsiMatrix(COMPLEX *iPsi11, COMPLEX *iPsi12, COMPLEX *iPsi21, COMPLEX *iPsi22, COMPLEX *kz, COMPLEX k, int vec_size, int pola)
 *
 *		\brief	Psi block diagonal matrix inverse
 */
/*-------------------------------------------------------------------------------------*/
int invPsiMatrix(COMPLEX *iPsi11, COMPLEX *iPsi12, COMPLEX *iPsi21, COMPLEX *iPsi22, COMPLEX *kz, COMPLEX k, int vec_size, int pola)
{
	int i;

	if (pola == TE){
		for(i=0;i<=vec_size-1;i++){
			iPsi11[i] =  0.5;
			iPsi12[i] =  0.5/kz[i];
			iPsi21[i] =  0.5;
			iPsi22[i] = -0.5/kz[i];
		}
	}else if (pola == TM){
		for(i=0;i<=vec_size-1;i++){
			iPsi11[i] = -0.5*k*k/kz[i];
			iPsi12[i] =  0.5;
			iPsi21[i] =  0.5*k*k/kz[i];
			iPsi22[i] =  0.5;
		}
	}else{
		fprintf(stderr, "%s, line %d : ERROR, unknown polarization... (\"%d\")\n",__FILE__,__LINE__,pola);
		exit(EXIT_FAILURE);
	}

	return 0;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn		int PsiMatrixTE(COMPLEX **Psi, COMPLEX *kz, struct Param_struct *par)
 *
 *		\brief	Psi matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int PsiMatrixTE(COMPLEX **Psi, COMPLEX k, COMPLEX *kz, struct Param_struct *par)
{
	int i,j;
	int vec_size = par->vec_size;

	for(i=0;i<=2*vec_size-1;i++){
		for(j=0;j<=2*vec_size-1;j++){
			Psi[i][j] = 0;
		}
	}
	for(i=0;i<=vec_size-1;i++){
		Psi[i][i]            = 1;
		Psi[i][i+  vec_size] = 1;
		Psi[i+vec_size][i]          = kz[i];
		Psi[i+vec_size][i+vec_size] = -kz[i];
	}

	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int PsiMatrixTE(COMPLEX **Psi, COMPLEX *kz, struct Param_struct *par)
 *
 *		\brief	Psi matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int PsiMatrixTM(COMPLEX **Psi, COMPLEX k, COMPLEX *kz, struct Param_struct *par)
{
	int i,j;
	int vec_size = par->vec_size;

	for(i=0;i<=2*vec_size-1;i++){
		for(j=0;j<=2*vec_size-1;j++){
			Psi[i][j] = 0;
		}
	}
	for(i=0;i<=vec_size-1;i++){
		Psi[i][i]            =  -kz[i]/(k*k);
		Psi[i][i+  vec_size] =  kz[i]/(k*k);
		Psi[i+vec_size][i]          = 1;
		Psi[i+vec_size][i+vec_size] = 1;
	}

	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int M_matrix(COMPLEX **M, double z, struct Param_struct *par)
 *
 *		\brief	M matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int M_matrix_TE(COMPLEX **M, double z, struct Param_struct *par)
{
	int i,j;
	int vec_size = par->vec_size;
	
	COMPLEX *sigma = par->sigma;
	COMPLEX *TF_k2    = par->TF_k2;

	/* Calculations of the Fourier transforms of the necessary 'grandeurs' */
	TF_k2 = FFT_k2_directe(z, TF_k2, par);

	/* Matrix initialisation */
	for(i=0;i<=2*vec_size-1;i++){
		for(j=0;j<=2*vec_size-1;j++){
			M[i][j] = 0;
		}
	}
	for(i=0;i<=vec_size-1;i++){
		/* M12 = -I*Id */
		M[i][i+vec_size] = -I;
		for(j=0;j<=vec_size-1;j++){
			/* M21 = I*sigma2*Id - I*k2(n-m) */
			M[i+vec_size][j] =  -I*TF_k2[i-j+2*par->N];
		}
		/* M21 = I*sigma2*Id - I*k2(n-m) */
		M[i+vec_size][i] += I*sigma[i]*sigma[i];
	}

	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int M_matrix(COMPLEX **M, double z, struct Param_struct *par)
 *
 *		\brief	M matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int M_matrix_TM(COMPLEX **M, double z, struct Param_struct *par)
{
	int i,j;
	int vec_size = par->vec_size;

	COMPLEX *sigma = par->sigma;
	COMPLEX **M11,**M12,**M21,**M22;
	COMPLEX **M_tmp1, **M_tmp2, **M_tmp3;
	COMPLEX **Qxx, **Qyy, **Qxz, **Qzz, **invQzz, **invQzzQxz, **invQzzsigma;
	Qxx = par->Qxx;
	Qyy = par->Qyy;
	Qxz = par->Qxz;
	Qzz = par->Qzz;
	invQzz = par->Qzz_1;
	M11 = par->M_tmp11;
	M12 = par->M_tmp12;
	M21 = par->M_tmp21;
	M22 = par->M_tmp22;
	M_tmp1 = par->M_tmp1;
	M_tmp2 = par->M_tmp2;
	M_tmp3 = par->M_tmp3;
	
	/* Toeplitz matrices calculations */
	md2D_QMatrix(z, par->Toep_k2, par->invToep_invk2, Qxx, Qyy, Qxz, Qzz, invQzz, par);

	/* invQzzQxz = invQzz*Qxz */
	invQzzQxz = M_tmp2;
	M_x_M(invQzzQxz, invQzz, Qxz, vec_size, vec_size);

	/* M11 = -sigma*invQzzQxz;*/
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M11[i][j] = -sigma[i]*invQzzQxz[i][j];
		}
	}

	/* M21 = Qxx - Qxz*invQzzQxz */
	sub_M(M21, Qxx, M_x_M(M_tmp1, Qxz, invQzzQxz, vec_size, vec_size), vec_size, vec_size);

	/* invQzzsigma = invQzz*sigma */
	invQzzsigma = M_tmp2;
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			invQzzsigma[i][j] = invQzz[i][j]*sigma[j];
		}
	}

	/* M12 = Oeil - sigma * invQzzsigma */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M12[i][j] = -sigma[i]*invQzzsigma[i][j];
		}
		M12[i][i] += 1;
	}

	/* M22 = -Qxz * invQzzsigma */
	M_x_M(M22, minus_M(M_tmp1,Qxz,vec_size,vec_size), invQzzsigma, vec_size, vec_size);

	/* putting it all together */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M[i][j]                   = I * M11[i][j];
			M[i][j+vec_size]          = I * M12[i][j];
			M[i+vec_size][j]          = I * M21[i][j];
			M[i+vec_size][j+vec_size] = I * M22[i][j];
		}
	}

	return 0;
}
/*-------------------------------------------------------------------------------------*/
/*!	\fn		int zinvar_M_Matrix_TM(COMPLEX **M, double z, struct Param_struct *par)
 *
 *		\brief	M matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int zinvar_M_matrix_TM(COMPLEX **M, double z, struct Param_struct *par)
{
	int i,j;
	int vec_size = par->vec_size;

	COMPLEX *sigma = par->sigma;
	COMPLEX **M12, **M_tmp2;
	COMPLEX **invToep_invk2, **invToep_k2, **Toep_k2;
/*	COMPLEX **Qxx, **invQzz,**M21,**M_tmp1,**M_tmp3;
	Qxx = par->Qxx;
	invQzz = par->Qzz_1;
	M_tmp1 = par->M_tmp1;
	M_tmp3 = par->M_tmp3;
	M21 = par->M_tmp21;*/
	M12 = par->M_tmp12;
	M_tmp2 = par->M_tmp2;
	Toep_k2       = par->Toep_k2;
	invToep_invk2 = par->invToep_invk2;

	/* Toeplitz matrices calculations */
	md2D_zinvarQMatrix(z, Toep_k2, invToep_invk2, par);

	/* M21 = Qxx */
/*	M_equals(M21, invToep_invk2, vec_size, vec_size);
*/

	/* M12 = Oeil - sigma * invToep_k2 * sigma */
	invToep_k2 = M_tmp2;
	invM(invToep_k2, Toep_k2, vec_size);
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M12[i][j] = -sigma[i]*invToep_k2[i][j]*sigma[j];
		}
		M12[i][i] += 1;
	}

	/* putting it all together */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M[i][j]                   = 0;
			M[i][j+vec_size]          = I * M12[i][j];
			M[i+vec_size][j]          = I * invToep_invk2[i][j];
			M[i+vec_size][j+vec_size] = 0;
		}
	}

	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_QMatrix(double z, COMPLEX **Toep_k2, COMPLEX **invToep_invk2, COMPLEX **Qxx, COMPLEX **Qyy, COMPLEX **Qxz, COMPLEX **Qzz, COMPLEX **Qzz_1, struct Param_struct *par)
 *
 *		\brief	Q matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int md2D_QMatrix(double z, COMPLEX **Toep_k2, COMPLEX **invToep_invk2, COMPLEX **Qxx, COMPLEX **Qyy, COMPLEX **Qxz, COMPLEX **Qzz, COMPLEX **Qzz_1, struct Param_struct *par)
{
	int i,j,var_tmp01;
	fftw_plan plan_TFk2, plan_TFinvk2, plan_TFNx2, plan_TFNz2, plan_TFNxNz;
	int N_x = par->N_x;
	int N_tf = 2*par->N;
	int vec_size = par->vec_size;
	double coefnorm;
	COMPLEX *tmp_k2, *tmp_invk2, *tmp_Nx2, *tmp_Nz2, *tmp_NxNz;
	tmp_k2    = par->tmp_tf_k2;
	tmp_invk2 = par->tmp_tf_invk2;
	tmp_Nx2   = par->tmp_tf_Nx2;
	tmp_Nz2   = par->tmp_tf_Nz2;
	tmp_NxNz  = par->tmp_tf_NxNz;

	COMPLEX *TF_k2, *TF_invk2, *TF_Nx2, *TF_Nz2, *TF_NxNz;
	TF_k2    = par->TF_k2;
	TF_invk2 = par->TF_invk2;
	TF_Nx2   = par->TF_Nx2;
	TF_Nz2   = par->TF_Nz2;
	TF_NxNz  = par->TF_NxNz;
	
	COMPLEX **Toep_invk2, **Toep_Nx2, **Toep_NxNz, **Toep_Nz2, **M_tmp1, **M_tmp2;
	Toep_invk2    = par->Toep_invk2;
	Toep_Nx2      = par->Toep_Nx2;
	Toep_NxNz     = par->Toep_NxNz;
	Toep_Nz2      = par->Toep_Nz2;
	M_tmp1 = par->M_tmp1;
	M_tmp2 = par->M_tmp2;
		
	/*--- Calculations of the Fourier transforms of the necessary 'grandeurs' ---*/

	/* Calculating the values in the direct space, using
	function pointors to adapt to surface definition type */
	(*par->k_2)(par, par->k2, z);
	(*par->invk_2)(par, par->invk2, z);
/*printf("\nRe(k2) : \n");SaveCplxTab2file (par->k2, par->N_x, "Re", "stdout"," ");
printf("\nIm(k2) : \n");SaveCplxTab2file (par->k2, par->N_x, "Im", "stdout"," ");
*/	(*par->Normal_function)(par, par->Nx2, par->NxNz, par->Nz2, z);
	
	/* Creating the 'plans' for the FFTW */
	plan_TFk2    = fftw_plan_dft_1d(N_x, (fftw_complex *)par->k2,    (fftw_complex *)tmp_k2,    FFTW_FORWARD, FFTW_ESTIMATE);	
	plan_TFinvk2 = fftw_plan_dft_1d(N_x, (fftw_complex *)par->invk2, (fftw_complex *)tmp_invk2, FFTW_FORWARD, FFTW_ESTIMATE);	
	plan_TFNx2   = fftw_plan_dft_1d(N_x, (fftw_complex *)par->Nx2,   (fftw_complex *)tmp_Nx2,   FFTW_FORWARD, FFTW_ESTIMATE);	
	plan_TFNz2   = fftw_plan_dft_1d(N_x, (fftw_complex *)par->Nz2,   (fftw_complex *)tmp_Nz2,   FFTW_FORWARD, FFTW_ESTIMATE);	
	plan_TFNxNz  = fftw_plan_dft_1d(N_x, (fftw_complex *)par->NxNz,  (fftw_complex *)tmp_NxNz,  FFTW_FORWARD, FFTW_ESTIMATE);	

	/* Calculating the FFT */
	/* (NOTICE : Real DFT could be used for Nx2, Nz2, etc, which would slightly increase speed but also code COMPLEXity) */	
	fftw_execute(plan_TFk2); 
	fftw_execute(plan_TFinvk2); 
	fftw_execute(plan_TFNx2); 
	fftw_execute(plan_TFNz2); 
	fftw_execute(plan_TFNxNz); 

	/* Keeping only the components between -N_tf & +N_tf */ 
	/* and normalizing by 1/N_x */
	coefnorm = 1.0/N_x;
	for (i=0;i<=N_tf-1;i++){
		TF_k2   [i] = tmp_k2   [N_x-N_tf+i] * coefnorm;
		TF_invk2[i] = tmp_invk2[N_x-N_tf+i] * coefnorm;
		TF_Nx2  [i] = tmp_Nx2  [N_x-N_tf+i] * coefnorm;
		TF_Nz2  [i] = tmp_Nz2  [N_x-N_tf+i] * coefnorm;
		TF_NxNz [i] = tmp_NxNz [N_x-N_tf+i] * coefnorm;
		TF_k2   [i+N_tf] = tmp_k2   [i] * coefnorm;
		TF_invk2[i+N_tf] = tmp_invk2[i] * coefnorm;
		TF_Nx2  [i+N_tf] = tmp_Nx2  [i] * coefnorm;
		TF_Nz2  [i+N_tf] = tmp_Nz2  [i] * coefnorm;
		TF_NxNz [i+N_tf] = tmp_NxNz [i] * coefnorm;
	}
	TF_k2   [2*N_tf] = tmp_k2   [N_tf] * coefnorm;
	TF_invk2[2*N_tf] = tmp_invk2[N_tf] * coefnorm;
	TF_Nx2  [2*N_tf] = tmp_Nx2  [N_tf] * coefnorm;
	TF_Nz2  [2*N_tf] = tmp_Nz2  [N_tf] * coefnorm;
	TF_NxNz [2*N_tf] = tmp_NxNz [N_tf] * coefnorm;

	/* Freeing memory */
	fftw_destroy_plan(plan_TFk2);
	fftw_destroy_plan(plan_TFinvk2);
	fftw_destroy_plan(plan_TFNx2);
	fftw_destroy_plan(plan_TFNz2);
	fftw_destroy_plan(plan_TFNxNz);
	
	/* Calculating the Toeplitz matrix for k2, 1/k2, Nx2, Nxz & Nz2 */
	for(i=0;i<=vec_size-1;i++){
		var_tmp01 = vec_size-1+i;
		for(j=0;j<=vec_size-1;j++){
			Toep_k2   [i][j] = TF_k2   [var_tmp01-j];
			Toep_invk2[i][j] = TF_invk2[var_tmp01-j];
			Toep_Nx2  [i][j] = TF_Nx2  [var_tmp01-j];
			Toep_NxNz [i][j] = TF_NxNz [var_tmp01-j];
			Toep_Nz2  [i][j] = TF_Nz2  [var_tmp01-j];
		}
	}
	/* Inverting Toep_invk2 to obtain invToep_invk2*/
	invM(invToep_invk2, Toep_invk2, vec_size);
	
	/*--- Building Q matrices ---*/
	
	/* Building Qyy = Toep_k2 */
	M_equals(Qyy, Toep_k2, vec_size, vec_size);

	/* Building Qxx = Toep_k2*Toep_Nz2 + invToep_invk2*Toep_Nx2 */
	M_x_M(M_tmp1, Toep_k2, Toep_Nz2, vec_size, vec_size);
	M_x_M(M_tmp2, invToep_invk2, Toep_Nx2, vec_size, vec_size);
	add_M(Qxx, M_tmp1, M_tmp2, vec_size, vec_size);

	/* Building Qzz = Toep_k2*Toep_Nx2 + invToep_invk2*Toep_Nz2 */
	M_x_M(M_tmp1, Toep_k2, Toep_Nx2, vec_size, vec_size);
	M_x_M(M_tmp2, invToep_invk2, Toep_Nz2, vec_size, vec_size);
	add_M(Qzz, M_tmp1, M_tmp2, vec_size, vec_size);

	/* Building Qxz = (invToep_invk2 - Toep_k2)*Toep_NxNz */
	sub_M(M_tmp1, invToep_invk2, Toep_k2, vec_size, vec_size);
	M_x_M(Qxz, M_tmp1, Toep_NxNz, vec_size, vec_size);

	/* Calcultating Qzz_1 from Qzz */
	invM(Qzz_1, Qzz, vec_size);
		
	return 0;
}
/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_zinvarQMatrix(double z, struct Param_struct *par)
 *
 *		\brief	Q matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int md2D_zinvarQMatrix(double z, COMPLEX **Toep_k2, COMPLEX **invToep_invk2, struct Param_struct *par)
{
	int i,j,var_tmp01;
	fftw_plan plan_TFk2, plan_TFinvk2;
	int N_x = par->N_x;
	int N_tf = 2*par->N;
	int vec_size = par->vec_size;
	double coefnorm;
	COMPLEX *tmp_k2, *tmp_invk2;
	tmp_k2    = par->tmp_tf_k2;
	tmp_invk2 = par->tmp_tf_invk2;

	COMPLEX *TF_k2, *TF_invk2;
	TF_k2    = par->TF_k2;
	TF_invk2 = par->TF_invk2;
	
	COMPLEX **Toep_invk2;
	Toep_invk2    = par->Toep_invk2;
		
	/*--- Calculations of the Fourier transforms of the necessary 'grandeurs' ---*/

	/* Calculating the values in the direct space, using
	function pointors to adapt to surface definition type */
	(*par->k_2)(par, par->k2, z);
	(*par->invk_2)(par, par->invk2, z);
	
	/* Creating the 'plans' for the FFTW */
	plan_TFk2    = fftw_plan_dft_1d(N_x, (fftw_complex *)par->k2,    (fftw_complex *)tmp_k2,    FFTW_FORWARD, FFTW_ESTIMATE);	
	plan_TFinvk2 = fftw_plan_dft_1d(N_x, (fftw_complex *)par->invk2, (fftw_complex *)tmp_invk2, FFTW_FORWARD, FFTW_ESTIMATE);	

	/* Calculating the FFT */
	/* (NOTICE : Real DFT could be used for Nx2, Nz2, etc, which would slightly increase speed but also code COMPLEXity) */	
	fftw_execute(plan_TFk2); 
	fftw_execute(plan_TFinvk2); 

	/* Keeping only the components between -N_tf & +N_tf */ 
	/* and normalizing by 1/N_x */
	coefnorm = 1.0/N_x;
	for (i=0;i<=N_tf-1;i++){
		TF_k2   [i] = tmp_k2   [N_x-N_tf+i] * coefnorm;
		TF_invk2[i] = tmp_invk2[N_x-N_tf+i] * coefnorm;
		TF_k2   [i+N_tf] = tmp_k2   [i] * coefnorm;
		TF_invk2[i+N_tf] = tmp_invk2[i] * coefnorm;
	}
	TF_k2   [2*N_tf] = tmp_k2   [N_tf] * coefnorm;
	TF_invk2[2*N_tf] = tmp_invk2[N_tf] * coefnorm;

	/* Freeing memory */
	fftw_destroy_plan(plan_TFk2);
	fftw_destroy_plan(plan_TFinvk2);
	
	/* Calculating the Toeplitz matrix for k2, 1/k2, Nx2, Nxz & Nz2 */
	for(i=0;i<=vec_size-1;i++){
		var_tmp01 = vec_size-1+i;
		for(j=0;j<=vec_size-1;j++){
			Toep_k2   [i][j] = TF_k2   [var_tmp01-j];
			Toep_invk2[i][j] = TF_invk2[var_tmp01-j];
		}
	}
	/* Inverting Toep_invk2 to obtain invToep_invk2*/
	invM(invToep_invk2, Toep_invk2, vec_size);

	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int k2_H_X(struct Param_struct *par, COMPLEX *k2_1D, double z)
 *
 *	\brief	Détermine le tableau de COMPLEXes k^2(x) pour un z donné
 *
 */
/*-------------------------------------------------------------------------------------*/
int k2_H_X(struct Param_struct *par, COMPLEX *k2_1D, double z)
{
	int i;
	
	for (i=0;i<=par->N_x-1;i++){
		if (z > par->profil[0][i]) 
			k2_1D[i] = par->k2_layer[0]; /* Superstrat */
		else
			k2_1D[i] = par->k2_layer[1]; /* Substrat */
	}
	
	return 0;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn		int k2_MULTI(struct Param_struct *par, COMPLEX *k2_1D, double z)
 *
 *	\brief	Détermine le tableau de COMPLEXes k^2(x) pour un z donné, pour un multicouches
 */
/*-------------------------------------------------------------------------------------*/
int k2_MULTI(struct Param_struct *par, COMPLEX *k2_1D, double z)
{
	int nx, n_layer=0;
	double eps = EPS*par->h;

	for (nx=0; nx<=par->N_x-1; nx++){
		do{
			if (z <= par->profil[n_layer][nx] + eps){
				if (z >= par->profil[n_layer+1][nx] - eps){
					k2_1D[nx] = par->k2_layer[n_layer];
					break;
				}else{
					n_layer++;
				}
			}else{
				n_layer--;
			}	
		}while (1);
	}
	
	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn			
 *
 *	\brief	Détermine le tableau de COMPLEXes 1/k^2(x) pour un z donné, pour un multicouches
 */
/*-------------------------------------------------------------------------------------*/
int invk2_MULTI(struct Param_struct *par, COMPLEX *invk2_1D, double z)
{
	int nx, n_layer=0;
	double eps = EPS*par->h;

	for (nx=0; nx<=par->N_x-1; nx++){
		do{
			if (z <= par->profil[n_layer][nx] + eps){
				if (z >= par->profil[n_layer+1][nx] - eps){
					invk2_1D[nx] = par->invk2_layer[n_layer];
					break;
				}else{
					n_layer++;
				}
			}else{
				n_layer--;
			}	
		}while (1);
	}
	
	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		k2_N_XYZ(struct Param_struct *par, COMPLEX *k2_1D, double z)	
 *
 *		\brief	Détermine le tableau de COMPLEXes k^2(x) pour un z donné
 */
/*-------------------------------------------------------------------------------------*/
int k2_N_XYZ(struct Param_struct *par, COMPLEX *k2_1D, double z)
{
	int i, nz;
	double DeuxPisurLambda2 = (2*PI/par->lambda)*(2*PI/par->lambda);
	double z_inv = par->h - z;
	
	nz = (int)((z_inv*par->N_z)/par->h);
	if (nz > par->N_z-1) nz = par->N_z-1;
	if (nz < 0)          nz = 0;
	
	for (i=0;i<=par->N_x-1;i++){
		k2_1D[i] = par->n_xyz[nz][i]*par->n_xyz[nz][i]*DeuxPisurLambda2; 
	}
	
	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn	invk2_N_XYZ(struct Param_struct *par, COMPLEX *invk2_1D, double z)	
 *
 *	\brief	Détermine le tableau de COMPLEXes invk^2(x) pour un z donné
 */
/*-------------------------------------------------------------------------------------*/
int invk2_N_XYZ(struct Param_struct *par, COMPLEX *invk2_1D, double z)
{
	int i, nz;
	double invDeuxPisurLambda2 = 1/(2*PI/par->lambda*2*PI/par->lambda);
	double z_inv = par->h - z;
	
	nz = (int)((z_inv*par->N_z)/par->h);
	if (nz > par->N_z-1) nz = par->N_z-1;
	if (nz < 0)          nz = 0;
	
	for (i=0;i<=par->N_x-1;i++){
		invk2_1D[i] = invDeuxPisurLambda2/(par->n_xyz[nz][i]*par->n_xyz[nz][i]); 
	}
	
	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int invk_2(COMPLEX *invk2_1D, double *profil, struct Param_struct *par, double z)
 *
 *	\brief	Détermine le tableau de COMPLEXes 1/k^2(x) pour un z donné
 *
 *	\todo	Prend pour l'instant en compte seulement un profil de type h(x)\n
 *			Doit être plus polyvalent : accepter aussi les profils de type n(x,z)
 */
/*-------------------------------------------------------------------------------------*/
int invk2_H_X(struct Param_struct *par, COMPLEX *invk2_1D, double z)
{
	
	int i;
	
	for (i=0;i<=par->N_x-1;i++){
		if (z > par->profil[0][i]) 
			invk2_1D[i] = par->invk2_layer[0]; /* Superstrat */
		else
			invk2_1D[i] = par->invk2_layer[1]; /* Substrat */
	}
	
	return 0;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn		int Normal_H_X(struct Param_struct *par, COMPLEX *k2_1D, double z)
 *
 *		\brief	Determine the Nx2, Nz2 and NxNz arrays, defined as \n
 * 				Nx2  = norm_x^2,  											\n
 * 				Nz2  = norm_z^2,  											\n
 * 				NxNz = norm_x norm_z 										\n
 * 				norm_x = -dhdx/(csqrt(1+dgdx^2)),      	   		\n
 * 				norm_z = 1/(csqrt(1+dgdx^2)),         					\n
 * 				dhdx being the h(x) derivative calculated as 		\n
 * 				dhdx = [h(x+dx)-h(x-dx)]/2dx
 *
 * 	\todo 	The PRECISION can be IMPROVED with a HIGHER ORDER calculation
 */
/*-------------------------------------------------------------------------------------*/
int Normal_H_X(struct Param_struct *par, COMPLEX *Nx2, COMPLEX *NxNz, COMPLEX *Nz2, double z)
{
	int i;
	int Nx = par->N_x; 	/* Caution : this 'Nx', corresponds to the number of points in x              */
								/* while 'Nx2', corresponds to the square of the surface normal x component */
	double dhdx, norm_x, norm_z;
	double two_dx = 2*par->L/Nx;

	double *profil = par->profil[0];

	if (par->HX_Normal_CALCULATED == 1){
/*		Nx2  = par->Nx2;
		NxNz = par->NxNz;
		Nz2  = par->Nz2;*/
		return 0;
	}else{		
		dhdx = (profil[1] - profil[Nx-1])/two_dx;
		norm_x = -dhdx/(csqrt(1+dhdx*dhdx));
		norm_z = 1.0/(csqrt(1+dhdx*dhdx));
		Nx2[0] = norm_x*norm_x;
		NxNz[0]= norm_x*norm_z;
		Nz2[0] = norm_z*norm_z;
		for (i=1;i<=Nx-2;i++){
			dhdx = (profil[i+1] - profil[i-1])/two_dx;
			norm_x = -dhdx/(csqrt(1+dhdx*dhdx));
			norm_z = 1.0/(csqrt(1+dhdx*dhdx));
			Nx2[i] = norm_x*norm_x;
			NxNz[i]= norm_x*norm_z;
			Nz2[i] = norm_z*norm_z;
		}
		dhdx = (profil[0] - profil[Nx-2])/two_dx;
		norm_x = -dhdx/(csqrt(1+dhdx*dhdx));
		norm_z = 1.0/(csqrt(1+dhdx*dhdx));
		Nx2[Nx-1] = norm_x*norm_x;
		NxNz[Nx-1]= norm_x*norm_z;
		Nz2[Nx-1] = norm_z*norm_z;

		par->HX_Normal_CALCULATED = 1;
	}
	return 0;
}
/*-------------------------------------------------------------------------------------*/
/*!	\fn		int Normal_H_X_Multi(struct Param_struct *par, COMPLEX *k2_1D, double z)
 *
 *		\brief	Determine the Nx2, Nz2 and NxNz arrays, defined as \n
 * 				Nx2  = norm_x^2,  											\n
 * 				Nz2  = norm_z^2,  											\n
 * 				NxNz = norm_x norm_z 										\n
 * 				norm_x = -dhdx/(csqrt(1+dgdx^2)),      	   		\n
 * 				norm_z = 1/(csqrt(1+dgdx^2)),         					\n
 * 				dhdx being the h(x) derivative calculated as 		\n
 * 				dhdx = [h(x+dx)-h(x-dx)]/2dx
 *
 * 	\todo 	The PRECISION can be IMPROVED with a HIGHER ORDER calculation
 */
/*-------------------------------------------------------------------------------------*/
int Normal_H_X_Multi(struct Param_struct *par, COMPLEX *Nx2, COMPLEX *NxNz, COMPLEX *Nz2, double z)
{
	int nx, nx_left, nx_right;
	int Nx = par->N_x; 	/* Caution : this 'Nx', corresponds to the number of points in x              */
								/* while 'Nx2', corresponds to the square of the surface normal x component */
	double dhdx, norm_x, norm_z, z1, z2, dhdx1, dhdx2;
	double two_dx = 2*par->L/Nx;

	double **profil = par->profil;

/* Arrays of Normal vector for each layer */

/* interpolation */
	int n_layer=0;
	double eps = EPS*par->h;

	for (nx=0; nx<=Nx-1; nx++){
		do{
			if (z <= profil[n_layer][nx] + eps){
				if (z >= profil[n_layer+1][nx] - eps){

					z1=profil[n_layer+1][nx];
					z2=profil[n_layer][nx];
					/* Slope calculated with (f(x+dx)-f(x-dx))/2dx */	
					/* Because of periodicity: indices Nx+1=0 and 0-1=Nx */	
					if (nx==0){
						nx_left  = Nx;
						nx_right = nx+1;
					}else if (nx==Nx-1){
						nx_left  = nx-1;
						nx_right = 0;
					}else{
						nx_left  = nx-1;
						nx_right = nx+1;
					} 
					/* Slopes of profiles surrounding the point under consideraiton */
					dhdx1 = (profil[n_layer+1][nx_right] - profil[n_layer+1][nx_left])/two_dx;
					dhdx2 = (profil[n_layer][nx_right] - profil[n_layer][nx_left])/two_dx;
					/* Interpolation of both slopes */	
					if (fabs(z2-z1) <= EPS*par->h){ /* in the particular case where z2==z1, use the average of the two slopes */
						dhdx = dhdx2/2 + dhdx1/2;
					}else{
						dhdx = (1/(z2-z1)) * ((z-z1)*dhdx2+(z2-z)*dhdx1); /* Otherwise, Linear interpolation */
					}
					norm_x = -dhdx/(csqrt(1+dhdx*dhdx));
					norm_z = 1.0/(csqrt(1+dhdx*dhdx));
					Nx2[nx] = norm_x*norm_x;
					NxNz[nx]= norm_x*norm_z;
					Nz2[nx] = norm_z*norm_z;
					break;
				}else{
					n_layer++;
				}
			}else{
				n_layer--;
			}	
		}while (1);
	}

	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		COMPLEX *FFT_k2_directe(double z, COMPLEX *TF_k2, struct Param_struct *par)
 *
 *		\brief	
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX *FFT_k2_directe(double z, COMPLEX *TF_k2, struct Param_struct *par)
{

	int i;
	COMPLEX *tmp;
	fftw_plan plan_TFk2;

	int N_x = par->N_x;
	int N_tf = 2*par->N;
	double coefnorm = 1.0/N_x;

	/* Calcul de k^2(x) à z fixé, à partir du profil */
	(*par->k_2)(par, par->k2, z);

	/* Calcul de la TF de k2, avec N_x points */
	
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!   ALLOUER A L'EXTERIEUR   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
	tmp = (COMPLEX *) malloc(sizeof(COMPLEX) * N_x);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

	plan_TFk2 = fftw_plan_dft_1d(N_x, (fftw_complex *)par->k2, (fftw_complex *)tmp, FFTW_FORWARD, FFTW_ESTIMATE);	
	fftw_execute(plan_TFk2); 

	/* On ne garde que les composantes entre -N_tf et +N_tf */ 
	/* et on normalise par 1/N_x */
	for (i=0;i<=N_tf-1;i++){
		TF_k2[i]   = tmp[N_x-N_tf+i] * coefnorm;
		TF_k2[i+N_tf] = tmp[i] * coefnorm;
	}
	TF_k2[2*N_tf] = tmp[N_tf] * coefnorm;

	fftw_destroy_plan(plan_TFk2);
	free(tmp);
/********** DEBUG **********/
/*printf("\nRe(k2)\n");SaveCplxTab2file (par->k2, par->N_x, "Re", "stdout", " ", 10000, "\n");
printf("\nIm(k2)\n");SaveCplxTab2file (par->k2, par->N_x, "Im", "stdout", " ", 10000, "\n");*/
/*printf("\nRe(TF_k2)\n");SaveCplxTab2file (TF_k2, 2*N_tf+1, "Re", "stdout", " ", 10000, "\n");
printf("\nIm(TF_k2)\n");SaveCplxTab2file (TF_k2, 2*N_tf+1, "Im", "stdout", " ", 10000, "\n");*/
/***************************/

	
	return TF_k2;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn		COMPLEX *FFT_invk2_directe(double z, COMPLEX *TF_invk2, struct Param_struct *par)
 *
 *	\brief	Calcule la TF de 1/k^2(x) pour un z donné et la tronque entre -N et +N
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX *FFT_invk2_directe(double z, COMPLEX *TF_invk2, struct Param_struct *par)
{

	int i;
	COMPLEX *tmp;
	fftw_plan plan_TFk2;

	int N_x = par->N_x;
	int N_tf = 2*par->N;
	double coefnorm = 1.0/N_x;

	/* Calcul de 1 / k^2(x) à z fixé, à partir du profil */
	(*par->invk_2)(par, par->invk2, z);

	/* Calcul de la TF de invk2, avec N_x points */
	tmp = (COMPLEX *) malloc(sizeof(COMPLEX) * N_x); /*TODO : ALLOUER A L'EXTERIEUR */
	plan_TFk2 = fftw_plan_dft_1d(N_x, (fftw_complex *)par->invk2, (fftw_complex *)tmp, FFTW_FORWARD, FFTW_ESTIMATE);	
	fftw_execute(plan_TFk2); 

	/* On ne garde que les composantes entre -N_tf et +N_tf */ 
	/* et on normalise par 1/N_x */
	for (i=0;i<=N_tf-1;i++){
		TF_invk2[i]   = tmp[N_x-N_tf+i] * coefnorm;
		TF_invk2[i+N_tf] = tmp[i] * coefnorm;
	}
	TF_invk2[2*N_tf] = tmp[N_tf] * coefnorm;

	fftw_destroy_plan(plan_TFk2);
	free(tmp);
	
	return TF_invk2;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_affichTemps(int n, int N, int nS, int NS, struct Param_struct *par)
 *
 *	\brief	Affichage du temps restant estimé en cours de calculs
 */
/*-------------------------------------------------------------------------------------*/
int md2D_affichTemps(int n, int N, int nS, int NS, int ni, int Ni, struct Param_struct *par)
{

	/* Si moins de 5 secondes depuis le dernier affichage, on ne change rien */
	if (CHRONO(clock(), par->last_clock) < 5){
		return 0;
	/* Sinon, estimation et affichage de la durée restante */
	}else /*if (par->VERBOSE >= 1)*/{
		int i, n_total;
		time(&par->last_time);
		par->last_clock = clock();
		float t_ecoule = difftime(par->last_time,par->time0);
		float t_total; 

		n_total = 4*(2*N+1);
/*		if (!strcmp(par->calcul_method,"DM")){
			t_total = t_ecoule*( Ni*NS*n_total)/(ni*NS*n_total+nS*n_total+n);
		}else if (!strcmp(par->calcul_method,"RCWA") || !strcmp(par->calcul_method,"IMPROVED_RCWA")){
			t_total = t_ecoule*NS/(nS+1);
		}else{
			fprintf(stderr, "%s, line %d : ERROR, unknown calculation method (\"%s\")\n",__FILE__,__LINE__,par->calcul_method);
		}*/
/*		t_total = t_ecoule*NS/(nS+1);*/
		t_total = t_ecoule*Ni*NS/( ni*NS+nS+1 );

		float t_restant = t_total - t_ecoule;
		int pourcent = ROUND(100.0*t_ecoule/t_total);
				
		fprintf(stdout,"\r");
		if (par->Ni > 1){
			 fprintf(stdout,"%3d %%, i = %d deg [%ds ", pourcent,ROUND(par->theta_i*180/PI),ROUND(t_ecoule));
		}else{
			fprintf(stdout,"%3d %% [%ds ", pourcent,ROUND(t_ecoule));
		}
		for (i=0;i<pourcent/5;i++) {fprintf(stdout,">");}
		for (i=pourcent/5;i<20;i++) {fprintf(stdout," ");}
		fprintf(stdout," %ds] %ds     ",ROUND(t_restant), ROUND(t_total));
		fflush(stdout);
	}
	
	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn	int md2D_save_S_matrix(double h_partial, struct Param_struct *par)
 *
 *	\brief	Sauvegarde la matrice S dans un fichier
 */
/*-------------------------------------------------------------------------------------*/
int md2D_save_S_matrix(double h_partial, struct Param_struct *par)
{
	int i,j,Nlign,Ncol;
	FILE *fp;
	char filename[SIZE_STR_BUFFER];
	Nlign = 2*par->N+1;
	Ncol  = 2*par->N+1;
	
	/* Nom de fichier */
	if (par->pola == TM){
		sprintf(filename,"S_TM_%s_h%f.txt",par->profile_name,h_partial);
	}else{
		sprintf(filename,"S_TE_%s_h%f.txt",par->profile_name,h_partial);
	}	
	fp = fopen(filename, "w");

	/* Ecriture de certains parametres */
	fprintf(fp,"mat_S_name = %s_h%f\n",par->profile_name,h_partial);
	fprintf(fp,"N = %d\n",par->N);
	fprintf(fp,"h_partial = %1.6e\n",h_partial);
	fprintf(fp,"Re_k_super = %1.6e\n",creal(par->k_super));
	fprintf(fp,"Im_k_super = %1.6e\n",cimag(par->k_super));
	fprintf(fp,"Re_k_sub = %1.6e\n",creal(par->k_sub));
	fprintf(fp,"Im_k_sub = %1.6e\n",cimag(par->k_sub));
	fprintf(fp,"L = %1.6e\n",par->L);
	fprintf(fp,"Delta_sigma = %1.6e\n",par->Delta_sigma);
	fprintf(fp,"sigma0 = %1.6e\n",creal(par->sigma0));
	/* Re(S12)*/
	fprintf(fp,"Re_S12 = ");
	for (j=0;j<=Nlign-1;j++){/* Remarque : lignes et colonnes sont permutees, permet un lecture plus facile acev lire_tab */
		for (i=0;i<=Ncol-1;i++){
			fprintf(fp,"% 1.12e  ",creal(par->S12[j][i]));
		}
		fprintf(fp,"\n");
	}
	/* Ima(S12) */
	fprintf(fp,"Im_S12 = ");
	for (j=0;j<=Nlign-1;j++){
		for (i=0;i<=Ncol-1;i++){
			fprintf(fp,"% 1.12e  ",cimag(par->S12[j][i]));
		}
		fprintf(fp,"\n");
	}
	/* Re(S22)*/
	fprintf(fp,"Re_S22 =");
	for (j=0;j<=Nlign-1;j++){
		for (i=0;i<=Ncol-1;i++){
			fprintf(fp,"% 1.12e  ",creal(par->S22[j][i]));
		}
		fprintf(fp,"\n");
	}
	/* Im(S22) */
	fprintf(fp,"Im_S22 = ");
	for (j=0;j<=Nlign-1;j++){
		for (i=0;i<=Ncol-1;i++){
			fprintf(fp,"% 1.12e  ",cimag(par->S22[j][i]));
		}
		fprintf(fp,"\n");
	}

	fclose(fp);

	return 0;	
}

/*!-------------------------------------------------------------------------------------
 *   \fn     int md2D_make_tab_S_steps(struct Param_struct* par)
 *
 *    \brief  Complete imposed S-steps by eventually adding S-step(s) so that there is at least one step every hmin_Sstep.
 *    Note 1: S-steps correspond to the discretization used by S-matrix algorithm.
 *    Note 2: Imposed S-steps are usefull for profiles with well defined horizontal interfaces (eg. multilayers) as opposed to profiles like rough surfaces or with sinusoidal shape for instance for which imposed S-steps are not usefull.
 *-------------------------------------------------------------------------------------*/
int md2D_make_tab_S_steps(struct Param_struct *par)
{

	int N_subDiv, num, nsub, ms;
	double eps = 1e-10;
	double *tab_imposed_S_steps = par->tab_imposed_S_steps;
	double *tab_NS = par->tab_NS;
	int N_imposed_S_steps = par->N_imposed_S_steps;
	double h0, h = par->h;
	double hmin_Ssteps=par->hmin_Sstep;

	/* h=0 */
	num = 0;
	if (fabs(tab_imposed_S_steps[0]) > eps){ /* if the first element of tab_imposed_S_steps is not 0. The first element of tabNS is set to 0 (otherwise it will be set, in the next loop, to the first element of tab_imposed_S_steps, which is also 0 */
		tab_NS[num++] = 0;
	}

	/* between h>0 and max(tab_imposed_S_step) */
	h0 = 0;
	for (ms=0;ms<=N_imposed_S_steps-1;ms++){
		N_subDiv=CEIL((tab_imposed_S_steps[ms]-h0)/hmin_Ssteps);
		for (nsub=2;nsub<=N_subDiv;nsub++){
			tab_NS[num++] = h0 + (nsub-1)*(tab_imposed_S_steps[ms]-h0)/N_subDiv;
		}
		tab_NS[num++]=tab_imposed_S_steps[ms];
		h0=tab_imposed_S_steps[ms];
	}

	/* Same thing between max(tab_imposed_S_steps) and h */
	if (fabs(tab_imposed_S_steps[N_imposed_S_steps-1]-h) > eps){ /* if the last element of tab_imposed_S_steps is not h */
		N_subDiv=CEIL((h-tab_imposed_S_steps[N_imposed_S_steps-1])/hmin_Ssteps);
		for (nsub=2;nsub<=N_subDiv;nsub++){
			tab_NS[num++] = tab_imposed_S_steps[N_imposed_S_steps-1] + (nsub-1)*(h-tab_imposed_S_steps[N_imposed_S_steps-1])/N_subDiv;
		}
		tab_NS[num++] = h;

	}
	par->NS=num-1;

/*printf("\ntab_NS[par->NS-1]=%1.10f\n",tab_NS[par->NS-1]);
printf("\npar->h=%1.10f\n",par->h);
printf("\ntab_imposed_S_steps[N_imposed_S_steps-1]=%f\n",tab_imposed_S_steps[N_imposed_S_steps-1]);
printf("\nN_subDiv=%d\n",N_subDiv);
printf("\ntab_NS=");SaveDbleTab2file(tab_NS, num, "stdout", " ", 1000, "\n");*/
	return 0;
}

/*!-------------------------------------------------------------------------------------
 *	\fn   int read_S_steps_from_profile(struct Param_struct *par)
 *
 *	\brief  Automatic determination of imposed S-steps from a profile.
 *
 *-------------------------------------------------------------------------------------*/
int read_S_steps_from_profile(struct Param_struct *par)
{

	int compare_doubles (const void *a, const void *b);
	double *tab_imposed_S_steps = par->tab_imposed_S_steps;
	int Ntab = par->N_imposed_S_steps;
	double **profil = par->profil;
	int N_layers = par->N_layers;
	int nx, n_layer, ntab, Nx = par->N_x;
	int value_in_tab=0;
	double eps=1e-10;

	/* Select all different values in profile to constitute tab_imposed_S_steps */
	Ntab=0;
	for (n_layer=0;n_layer<=N_layers+2;n_layer++){
		for (nx=0;nx<Nx;nx++){
			value_in_tab=0;
			for (ntab=0;ntab<Ntab;ntab++){
				if (fabs(profil[n_layer][nx]-tab_imposed_S_steps[ntab])< eps){
					value_in_tab=1;
				}
			}
			if (value_in_tab==0){
				tab_imposed_S_steps[Ntab]=profil[n_layer][nx];
				Ntab++;
			}
		}
	}

	/* Sort tab_imposed_S_steps */
	qsort (tab_imposed_S_steps, Ntab, sizeof(double), compare_doubles);

	par->N_imposed_S_steps = Ntab;
/*for (n_layer=0;n_layer<=N_layers+2;n_layer++){
printf("\nprofile=");SaveDbleTab2file(profil[n_layer], Nx, "stdout", " ", 1000, "\n");}
printf("\ntab_imposed_S_steps=");SaveDbleTab2file(tab_imposed_S_steps, Ntab, "stdout", " ", 1000, "\n");*/
	return 0;
}


int compare_doubles (const void *a, const void *b)
{
  double* a1 = (double*) a;
  double* b1 = (double*) b;
  double temp = *a1 - *b1;
  if (temp > 0)
    return 1;
  else if (temp < 0)
    return -1;
  else
    return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn	int md2D_save_near_field(COMPLEX **S12, COMPLEX **Z, int vec_size, int nS)
 *
 *	\brief	Save Z and S12 matrices for near field reconstruction
 */
/*-------------------------------------------------------------------------------------*/
int md2D_save_near_field(COMPLEX **S12, COMPLEX **Z, int vec_size, int nS, struct Param_struct *par)
{
	int i, j;
	
	for (i=0;i<vec_size;i++){
		for (j=0;j<vec_size;j++){
			par->tab_Z[nS][i][j]   = Z[i][j];
			par->tab_S12[nS][i][j] = S12[i][j];
		}
	}

	return 0;	
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn	int md2D_save_near_field(COMPLEX **S12, COMPLEX **Z, int vec_size, int nS)
 *
 *	\brief	Save Z and S12 matrices for near field reconstruction
 */
/*-------------------------------------------------------------------------------------*/
int md2D_save_T(COMPLEX **T11, COMPLEX **T12, COMPLEX **T21, COMPLEX **T22, COMPLEX ***saved_T, int vec_size, int nS, struct Param_struct *par)
{
printf("function not implemented. Uncomment following lines\n");
#if 0	
	int i, j;
	for (i=0;i<vec_size;i++){
		for (j=0;j<vec_size;j++){
			par->saved_T[nS][i]         [j]            = T11[i][j]; /* A VERIFIER : Attention à ne pas inverser ligne et col */
			par->saved_T[nS][i]         [j+vec_size]   = T12[i][j];
			par->saved_T[nS][i+vec_size][j]            = T21[i][j];
			par->saved_T[nS][i+vec_size][j+vec_size]   = T22[i][j];
		}
	}
#endif
	return 0;	
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn	int md2D_near_field_map(par)
 *
 *	\brief	Map field reconstruction
 */
/*-------------------------------------------------------------------------------------*/
int md2D_local_field_map(COMPLEX ***tab_S12, COMPLEX ***tab_Z, struct Param_struct *par)
{
	int q,i,j,n,nx, NS, N_x, vec_size;
	COMPLEX *Vi, *Vq, *Fq, *Vq_m, *Vq_p, *Eq, *Hq, *EHqz, *sigma, *kz_super, k_super2, k_super, *kz_sub, k_sub, k_sub2, k2;
	COMPLEX *Psi11, *Psi12, *Psi21, *Psi22, **Z_prod, **Z_prod_tmp, **Hpx, **Hpy, **Hpz, **Ex, **Ey, **Ez;
	double x;

	NS = par->NS;
	N_x = par->N_x;
	vec_size = par->vec_size;
	sigma = par->sigma;
	kz_super = par->kz_super;
	k_super = par->k_super;
	kz_sub = par->kz_sub;
	k_sub = par->k_sub;
	k_super2 = par->k_super*par->k_super;
	k_sub2 = par->k_sub*par->k_sub;
	
	Z_prod     = allocate_CplxMatrix(vec_size, vec_size);
	Z_prod_tmp = allocate_CplxMatrix(vec_size, vec_size);
	Psi11 = malloc(sizeof(COMPLEX)*vec_size); Psi12 = malloc(sizeof(COMPLEX)*vec_size);
	Psi21 = malloc(sizeof(COMPLEX)*vec_size); Psi22 = malloc(sizeof(COMPLEX)*vec_size);
	
	if (par->pola == TE){
		Ex = NULL;
		Ez = NULL;
		Hpy = NULL;
		Hpx = allocate_CplxMatrix(NS,N_x);
		Hpz = allocate_CplxMatrix(NS,N_x);
		Ey  = allocate_CplxMatrix(NS,N_x);
	}else{ /* TM */
		Ex = allocate_CplxMatrix(NS,N_x);
		Ez = allocate_CplxMatrix(NS,N_x);
		Hpy = allocate_CplxMatrix(NS,N_x);
		Hpx = NULL;
		Hpz = NULL;
		Ey = NULL;
	}
	
	Vi = (COMPLEX *) malloc(sizeof(COMPLEX)*vec_size);
	Vq = (COMPLEX *) malloc(sizeof(COMPLEX)*2*vec_size);
	Fq = (COMPLEX *) malloc(sizeof(COMPLEX)*2*vec_size);
	EHqz = (COMPLEX *) malloc(sizeof(COMPLEX)*vec_size);
	Vq_m = Vq;
	Vq_p = Vq + vec_size;
	Eq = Fq;
	Hq = Fq + vec_size;

	/* Incident field */
	for (i=0;i<vec_size;i++){
		Vi[i]=0;
	}
	Vi[par->vec_middle] = 1.0;

	/* Z_prod = Id */
	for (i=0;i<vec_size;i++){
		for (j=0;j<vec_size;j++){
			Z_prod[i][j] = 0.0;
		}
		Z_prod[i][i] = 1.0;
	}
	
	for (q=NS-1;q>=0;q--){
		/* Z_prod_tmp = Zprod */
		M_equals(Z_prod_tmp, Z_prod, vec_size, vec_size);
		
		/* Z_prod = Zq*Z_prod_tmp */
		M_x_M(Z_prod, tab_Z[q], Z_prod_tmp, vec_size, vec_size);
		
		/* Vq_m = Z_prod*Vi */
		M_x_V(Vq_m, Z_prod, Vi, vec_size, vec_size);
		
		/* Vq_p = S12q*Vq_m */
		M_x_V(Vq_p, tab_S12[q], Vq_m, vec_size, vec_size);

		/* Fq = Psi*Vq */
		if (q==0){
			PsiMatrix(Psi11, Psi12, Psi21, Psi22, kz_sub, k_sub, vec_size, par->pola);
			k2=k_sub2;
		}else{
			PsiMatrix(Psi11, Psi12, Psi21, Psi22, kz_super, k_super, vec_size, par->pola);
			k2=k_super2;
		}
		for (i=0;i<=vec_size-1;i++){
			/* Eq = Psi11 * Vq_m + Psi12 * Vq_p */
			Eq[i] = Psi11[i]*Vq_m[i] + Psi12[i]*Vq_p[i];
			/* Hq = Psi21 * Vq_m + Psi22 * Vq_p */
			Hq[i] = Psi21[i]*Vq_m[i] + Psi22[i]*Vq_p[i];
			if (par->pola == TE){
				EHqz[i] = sigma[i]*(Vq_m[i] + Vq_p[i]);
			}else{ /* TM */
				EHqz[i] = sigma[i]*(Vq_m[i] + Vq_p[i])/k2;
			}
		}

		/* Near field in the real space */
		if (par->pola == TE){
			for (nx=0;nx<=N_x-1;nx++){
				x = nx*par->L/N_x;
				Ey[q][nx]  = 0+0*I;
				Hpx[q][nx] = 0+0*I;
				Hpz[q][nx] = 0+0*I;
				for (n=0;n<vec_size;n++){
					Ey[q][nx]  += Eq[n]*cexp(I*sigma[n]*x);
					Hpx[q][nx] += Hq[n]*cexp(I*sigma[n]*x);
					Hpz[q][nx] += EHqz[n]*cexp(I*sigma[n]*x);
				}
			}
		}else{ /* TM */
			for (nx=0;nx<=N_x-1;nx++){
				x = nx*par->L/N_x;
				Ex[q][nx] = 0+0*I;
				Ez[q][nx] = 0+0*I;
				Hpy[q][nx] = 0+0*I;
				for (n=0;n<vec_size;n++){
					Ex[q][nx]  += Eq[n]*cexp(I*sigma[n]*x);
					Ez[q][nx]  += EHqz[n]*cexp(I*sigma[n]*x);
					Hpy[q][nx] += Hq[n]*cexp(I*sigma[n]*x); 
				}
			}
		}
	}

	/* Writing results */
	if (par->pola == TE){
		fprintf(stdout,"ReEy = \n"); SaveMatrix2file(Ey, NS, N_x, "Re", "stdout");
		fprintf(stdout,"ImEy = \n"); SaveMatrix2file(Ey, NS, N_x, "Im", "stdout");
		fprintf(stdout,"ReHpx = \n");SaveMatrix2file(Hpx, NS, N_x, "Re", "stdout");
		fprintf(stdout,"ImHpx = \n");SaveMatrix2file(Hpx, NS, N_x, "Im", "stdout");
		fprintf(stdout,"ReHpz = \n");SaveMatrix2file(Hpz, NS, N_x, "Re", "stdout");
		fprintf(stdout,"ImHpz = \n");SaveMatrix2file(Hpz, NS, N_x, "Im", "stdout");
	}else{ /* TM */
		fprintf(stdout,"ReEx = \n"); SaveMatrix2file(Ex, NS, N_x, "Re", "stdout");
		fprintf(stdout,"ImEx = \n"); SaveMatrix2file(Ex, NS, N_x, "Im", "stdout");
		fprintf(stdout,"ReEz = \n"); SaveMatrix2file(Ez, NS, N_x, "Re", "stdout");
		fprintf(stdout,"ImEz = \n"); SaveMatrix2file(Ez, NS, N_x, "Im", "stdout");
		fprintf(stdout,"ReHpy = \n");SaveMatrix2file(Hpy, NS, N_x, "Re", "stdout");
		fprintf(stdout,"ImHpy = \n");SaveMatrix2file(Hpy, NS, N_x, "Im", "stdout");
	}

	/* freeing memory... */
	if (par->pola == TE){
		free(Hpx[0]);
		free(Hpx);
		free(Hpz[0]);
		free(Hpz);
		free(Ey[0]);
		free(Ey);
	}else{ /* TM */
		free(Ex[0]);
		free(Ex);
		free(Ez[0]);
		free(Ez);
		free(Hpy[0]);
		free(Hpy);
	}

	free(Z_prod[0]);
	free(Z_prod);
	free(Z_prod_tmp[0]);
	free(Z_prod_tmp);
	free(Vi);
	free(Vq);
	free(Fq);
	free(EHqz);
	free(Psi11);free(Psi12);free(Psi21);free(Psi22);

	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn	int md2D_local_field_by_T_products(struct Param_struct *par)
 *
 *	\brief	Field map reconstruction using simple T-matrix integration
 */
/*-------------------------------------------------------------------------------------*/
int md2D_local_field_map_by_T_products(COMPLEX *V0m, struct Param_struct *par)
{
	COMPLEX **T,**T11,**T12,**T21,**T22,**T_tmp,**T_local,**T_tmp2,**T11_tmp,**T12_tmp,**T21_tmp,**T22_tmp;
	int q,i,j,n,nx, NS, N_x, vec_size;
	COMPLEX *V1_m, *Vq, *Fq, *Vq_m, *Vq_p, *Eq, *Hq, *EHqz, *sigma, *kz_super, k_super2, k_super, *kz_sub, k_sub, k_sub2, k2;
	COMPLEX *Psi11, *Psi12, *Psi21, *Psi22, **Hpx, **Hpy, **Hpz, **Ex, **Ey, **Ez;
	double x;

	NS = par->NS;
	N_x = par->N_x;
	vec_size = par->vec_size;
	sigma = par->sigma;
	kz_super = par->kz_super;
	k_super = par->k_super;
	kz_sub = par->kz_sub;
	k_sub = par->k_sub;
	k_super2 = par->k_super*par->k_super;
	k_sub2 = par->k_sub*par->k_sub;
	
	Psi11 = malloc(sizeof(COMPLEX)*vec_size); Psi12 = malloc(sizeof(COMPLEX)*vec_size);
	Psi21 = malloc(sizeof(COMPLEX)*vec_size); Psi22 = malloc(sizeof(COMPLEX)*vec_size);
	
	if (par->pola == TE){
		Ex = NULL;
		Ez = NULL;
		Hpy = NULL;
		Hpx = allocate_CplxMatrix(NS,N_x);
		Hpz = allocate_CplxMatrix(NS,N_x);
		Ey  = allocate_CplxMatrix(NS,N_x);
	}else{ /* TM */
		Ex = allocate_CplxMatrix(NS,N_x);
		Ez = allocate_CplxMatrix(NS,N_x);
		Hpy = allocate_CplxMatrix(NS,N_x);
		Hpx = NULL;
		Hpz = NULL;
		Ey = NULL;
	}
	
	V1_m = (COMPLEX *) malloc(sizeof(COMPLEX)*vec_size);
	Vq = (COMPLEX *) malloc(sizeof(COMPLEX)*2*vec_size);
	Fq = (COMPLEX *) malloc(sizeof(COMPLEX)*2*vec_size);
	EHqz = (COMPLEX *) malloc(sizeof(COMPLEX)*vec_size);
	Vq_m = Vq;
	Vq_p = Vq + vec_size;
	Eq = Fq;
	Hq = Fq + vec_size;

	T = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	T11 = (COMPLEX **) malloc(sizeof(COMPLEX *)*par->vec_size);
	T12 = (COMPLEX **) malloc(sizeof(COMPLEX *)*par->vec_size);
	T21 = (COMPLEX **) malloc(sizeof(COMPLEX *)*par->vec_size);
	T22 = (COMPLEX **) malloc(sizeof(COMPLEX *)*par->vec_size);
	for (i=0;i<=par->vec_size-1;i++){
			T11[i] = &T[i][0];
			T12[i] = &T[i][par->vec_size];
			T21[i] = &T[i+par->vec_size][0];
			T22[i] = &T[i+par->vec_size][par->vec_size];
	}
	T_tmp = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	T11_tmp = (COMPLEX **) malloc(sizeof(COMPLEX *)*par->vec_size);
	T12_tmp = (COMPLEX **) malloc(sizeof(COMPLEX *)*par->vec_size);
	T21_tmp = (COMPLEX **) malloc(sizeof(COMPLEX *)*par->vec_size);
	T22_tmp = (COMPLEX **) malloc(sizeof(COMPLEX *)*par->vec_size);
	for (i=0;i<=par->vec_size-1;i++){
			T11_tmp[i] = &T_tmp[i][0];
			T12_tmp[i] = &T_tmp[i][par->vec_size];
			T21_tmp[i] = &T_tmp[i+par->vec_size][0];
			T22_tmp[i] = &T_tmp[i+par->vec_size][par->vec_size];
	}
	T_tmp2 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	T_local=T_tmp;

	/* Field at bottom of structure */
	for (i=0;i<vec_size;i++){
		V1_m[i]=V0m[i];
	}

	/* T-Matrix initialisation, T= Id */
	for (i=0; i<=2*vec_size-1; i++) {
		for (j=0; j<=2*vec_size-1; j++) {
			T[i][j] = 0;
		}
		T[i][i] = 1;
	}

	/* Iterations */
/*NS--;*/
	for (q=0; q<=NS-1; q++) {

		/* T-Matrix calculation */
		T_Matrix(T11_tmp, T12_tmp, T21_tmp, T22_tmp, q, par);

		/* T-Matrix incrementation */
		M_x_M(T_tmp2,T_local,T,2*vec_size, 2*vec_size);
		M_equals(T,T_tmp2,2*vec_size, 2*vec_size);

		/* Vq_m = T11*V1_m */
		M_x_V(Vq_m, T11, V1_m, vec_size, vec_size);
		
		/* Vq_p = T21*V1_m */
		M_x_V(Vq_p, T21, V1_m, vec_size, vec_size);

		/* Fq = Psi*Vq */
		if (q==0){
			PsiMatrix(Psi11, Psi12, Psi21, Psi22, kz_sub, k_sub, vec_size, par->pola);
			k2=k_sub2;
		}else{
			PsiMatrix(Psi11, Psi12, Psi21, Psi22, kz_super, k_super, vec_size, par->pola);
			k2=k_super2;
		}
		for (i=0;i<=vec_size-1;i++){
			/* Eq = Psi11 * Vq_m + Psi12 * Vq_p */
			Eq[i] = Psi11[i]*Vq_m[i] + Psi12[i]*Vq_p[i];
			/* Hq = Psi21 * Vq_m + Psi22 * Vq_p */
			Hq[i] = Psi21[i]*Vq_m[i] + Psi22[i]*Vq_p[i];
			if (par->pola == TE){
				EHqz[i] = sigma[i]*(Vq_m[i] + Vq_p[i]);
			}else{ /* TM */
				EHqz[i] = sigma[i]*(Vq_m[i] + Vq_p[i])/k2;
			}
		}

		/* Near field in the real space */
		if (par->pola == TE){
			for (nx=0;nx<=N_x-1;nx++){
				x = nx*par->L/N_x;
				Ey[q][nx]  = 0+0*I;
				Hpx[q][nx] = 0+0*I;
				Hpz[q][nx] = 0+0*I;
				for (n=0;n<vec_size;n++){
					Ey[q][nx]  += Eq[n]*cexp(I*sigma[n]*x);
					Hpx[q][nx] += Hq[n]*cexp(I*sigma[n]*x);
					Hpz[q][nx] += EHqz[n]*cexp(I*sigma[n]*x);
				}
			}
		}else{ /* TM */
			for (nx=0;nx<=N_x-1;nx++){
				x = nx*par->L/N_x;
				Ex[q][nx] = 0+0*I;
				Ez[q][nx] = 0+0*I;
				Hpy[q][nx] = 0+0*I;
				for (n=0;n<vec_size;n++){
					Ex[q][nx]  += Eq[n]*cexp(I*sigma[n]*x);
					Ez[q][nx]  += EHqz[n]*cexp(I*sigma[n]*x);
					Hpy[q][nx] += Hq[n]*cexp(I*sigma[n]*x); 
				}
			}
		}
	}

/* ICI Calculer les amplitudes puis les efficacités, puis écrire les résultats... */
	for (i=0;i<=vec_size-1;i++){
		par->Ar[i]=Vq_p[i];
		par->At[i]=V1_m[i];
		par->Ai[i]=1;
	}

	/* Writing the sum of T11 middle column (for root searching purpose) */
	double sumT11=0;
	for (i=0;i<=vec_size-1;i++){
		sumT11 += cabs(Vq_m[i]);
	}
	fprintf(stdout,"sumT11 = %lf\n",sumT11);
	
	/* Writing results */
	if (par->pola == TE){
		fprintf(stdout,"ReEy = \n"); SaveMatrix2file(Ey, NS, N_x, "Re", "stdout");
		fprintf(stdout,"ImEy = \n"); SaveMatrix2file(Ey, NS, N_x, "Im", "stdout");
		fprintf(stdout,"ReHpx = \n");SaveMatrix2file(Hpx, NS, N_x, "Re", "stdout");
		fprintf(stdout,"ImHpx = \n");SaveMatrix2file(Hpx, NS, N_x, "Im", "stdout");
		fprintf(stdout,"ReHpz = \n");SaveMatrix2file(Hpz, NS, N_x, "Re", "stdout");
		fprintf(stdout,"ImHpz = \n");SaveMatrix2file(Hpz, NS, N_x, "Im", "stdout");
	}else{ /* TM */
		fprintf(stdout,"ReEx = \n"); SaveMatrix2file(Ex, NS, N_x, "Re", "stdout");
		fprintf(stdout,"ImEx = \n"); SaveMatrix2file(Ex, NS, N_x, "Im", "stdout");
		fprintf(stdout,"ReEz = \n"); SaveMatrix2file(Ez, NS, N_x, "Re", "stdout");
		fprintf(stdout,"ImEz = \n"); SaveMatrix2file(Ez, NS, N_x, "Im", "stdout");
		fprintf(stdout,"ReHpy = \n");SaveMatrix2file(Hpy, NS, N_x, "Re", "stdout");
		fprintf(stdout,"ImHpy = \n");SaveMatrix2file(Hpy, NS, N_x, "Im", "stdout");
	}

	/* freeing memory... */
	if (par->pola == TE){
		free(Hpx[0]);
		free(Hpx);
		free(Hpz[0]);
		free(Hpz);
		free(Ey[0]);
		free(Ey);
	}else{ /* TM */
		free(Ex[0]);
		free(Ex);
		free(Ez[0]);
		free(Ez);
		free(Hpy[0]);
		free(Hpy);
	}

	free(T[0]);
	free(T);
	free(T_tmp[0]);
	free(T_tmp);
	free(T_tmp2[0]);
	free(T_tmp2);
	free(T11);
	free(T12);
	free(T21);
	free(T22);
	free(T11_tmp);
	free(T12_tmp);
	free(T21_tmp);
	free(T22_tmp);
	free(V1_m);
	free(Vq);
	free(Fq);
	free(EHqz);
	free(Psi11);free(Psi12);free(Psi21);free(Psi22);

	return 0;
}


#if 0
/*-------------------------------------------------------------------------------------*/
/*!	\fn		int oldT_Matrix(COMPLEX **T11, COMPLEX **T12, COMPLEX **T21, COMPLEX **T22,
				int nS, struct Param_struct *par)
 *
 *	\brief Calcul de la matrice T en polarisation TE
 */
/*-------------------------------------------------------------------------------------*/
int oldT_Matrix(COMPLEX **T11, COMPLEX **T12, COMPLEX **T21, COMPLEX **T22, int nS, struct Param_struct *par)
{
	int vec_size = par->vec_size;
	
	double hmin, hmax, z, Delta_z;
	COMPLEX **Psi, **invPsi_super;

	/* [F] vectors are defined by 	 	[V] vectors by
	[F] = |[Ex ]|								[V] = |[VE-]|
			|[Ey ]|								      |[VH-]|
			|[H'x]|								      |[VE+]|
			|[H'y]|  							      |[VH+]|  
	with H'=omega mu H											*/

	/* Psi matrix */
	if (nS==0){ /* 1st S-Matrix iteration : we are in the substrat */
		if (par->pola == TE){
			Psi = par->Psi_sub_TE;
		}else{
			Psi = par->Psi_sub_TM;}
	}else{     /* Following iterations : we are in the superstrat */
		if (par->pola == TE){
			Psi = par->Psi_super_TE;
		}else{
			Psi = par->Psi_super_TM;}
	}

	if (par->pola == TE){
		invPsi_super = par->invPsi_super_TE;
	}else{ /* TM */
		invPsi_super = par->invPsi_super_TM;
	}

	/* z of the considered T matrix slice */
	hmin = par->h*nS/par->NS;
	hmax = par->h*(nS+1)/par->NS;
	if (par->tab_NS_ENABLED){
		hmin = par->tab_NS[nS];
		hmax = par->tab_NS[nS+1];
	}
	Delta_z = hmax - hmin;
	z = (hmax + hmin)/2;


	/* P_matrix calculation */
	(*par->P_matrix)(par->P, z, Delta_z, par);

	/* T matrix : T = inv(Psi_super) * P_matrix * Psi */
	M_x_M(par->T,
			invPsi_super, M_x_M(par->M_buffer_2vecsize, par->P, Psi, 2*vec_size, 2*vec_size), 2*vec_size, 2*vec_size);
/*printf("\nRe(par->T) :\n");
SaveMatrix2file (par->T, 2*par->vec_size, 2*par->vec_size, "Re", "stdout");
printf("\nIm(par->T) :\n");
SaveMatrix2file (par->T, 2*par->vec_size, 2*par->vec_size, "Im", "stdout");*/


	/*Affichage du temps restant à l'écran */
	md2D_affichTemps(par->N,par->N,nS,par->NS,par->ni,par->Ni,par);
	
	return 0;
}
#endif

