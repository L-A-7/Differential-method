/*!	\file		md2D.c
 *
 * 	\brief		Differential method \n
 *						2-Dimensions (as opposed to 3-D) \n
 *						S-Matrices algorithm \n
 * 					FFF algorithm \n
 * 					conical incidence
 *
 *	\date		january 2007
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
	complex Eyi, Hpyi;
	double phi_i = par->phi_i;
	double psi = par->psi;
	double theta_i = par->theta_i;
	complex k_super = par->k_super;
	int vec_size = par->vec_size;
	int vec_mid = par->vec_middle;
	/* Plane wave */
	if (!strcmp(par->i_field_mode,"PLANE_WAVE")){

		Eyi  = cos(phi_i)*cos(psi) + sin(phi_i)*cos(theta_i)*sin(psi);
		Hpyi = (sin(phi_i)*cos(theta_i)*cos(psi) - cos(phi_i)*sin(psi))*k_super;

		for (n=0; n<=2*vec_size-1; n++) {
			par->Ai[n] = 0;
		}
		par->Ai[vec_mid] = Eyi;
		par->Ai[vec_mid+vec_size] = Hpyi;
		 
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
	complex k_super2 = par->k_super*par->k_super;
	complex k_sub2 = par->k_sub*par->k_sub;
	complex ky_02 = par->ky_0*par->ky_0;
	
	double sigma0 = par->sigma0;
	int Nmin_super, Nmax_super, Nmin_sub, Nmax_sub;

	/* Propagativ modes limits */
	Nmax_super =  FLOOR( ( sqrt(cabs(k_super2 - ky_02)) - sigma0 )/Delta_sigma );
	Nmin_super = -FLOOR( ( sqrt(cabs(k_super2 - ky_02)) + sigma0 )/Delta_sigma );
	Nmax_sub   =  FLOOR( ( sqrt(cabs(k_sub2    - ky_02)) - sigma0 )/Delta_sigma );
	Nmin_sub   = -FLOOR( ( sqrt(cabs(k_sub2    - ky_02)) + sigma0 )/Delta_sigma );
	int Nlimit = MAX(MAX(Nmax_super,Nmax_sub),MAX(-Nmin_super,-Nmin_sub));
	if (N < Nlimit) {
		fprintf(stderr,"ATTENTION, N trop petit pour représenter l'ensemble des modes propagatifs\n");
		fprintf(stderr,"valeur minimale : N = %d \n",Nlimit);
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
/*!	\fn		int md2D_efficiencies(complex *Ai, complex *A0, complex *Ah, struct Param_struct *par,  struct Efficacites_struct *eff)
 *
 *	\brief		Efficiencies calculation
 */
/*-------------------------------------------------------------------------------------*/
int md2D_efficiencies(complex *Ai, complex *Ar, complex *At, struct Param_struct *par,  struct Efficacites_struct *eff)
{
	int n;
	int vec_size = par->vec_size;
	int mid = par->vec_middle;
	
	complex k_super = par->k_super;
	complex k_sub = par->k_sub;
	complex k_super2 = k_super*k_super;
	complex k_sub2 = k_sub*k_sub;
	complex ky_0  = par->ky_0;
	complex ky_02 = ky_0*ky_0;
	
	complex *kz_super, *kz_sub, *sigma;
	kz_super = par->kz_super;
	kz_sub = par->kz_sub;
	sigma = par->sigma;

	double Pz_i, Pz_r, Pz_t, sumPzi;
	complex *Eyi, *Hpyi, *Eyr, *Hpyr, *Eyt, *Hpyt;
	complex *Exi, *Hpxi, *Exr, *Hpxr, *Ext, *Hpxt;
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
	Eyi  = par->Ai;
	Hpyi = par->Ai + vec_size;
	Eyr  = par->Ar;
	Hpyr = par->Ar + vec_size;
	Eyt  = par->At;
	Hpyt = par->At + vec_size;

	/* x field components calculations */
	complex C_super = (1/(k_super2 - ky_02));
	complex C_sub   = (1/(k_sub2   - ky_02));
	for(n=0;n<=vec_size-1;n++){
		Exi[n]  = C_super * (-kz_super[n]*Hpyi[n] - ky_0*sigma[n]*Eyi[n]);
		Exr[n]  = C_super * ( kz_super[n]*Hpyr[n] - ky_0*sigma[n]*Eyr[n]);
		Ext[n]  = C_sub   * (-kz_sub[n]*Hpyt[n]   - ky_0*sigma[n]*Eyt[n]);
		Hpxi[n] = C_super * ( k_super2*kz_super[n]*Eyi[n] - ky_0*sigma[n]*Hpyi[n]);
		Hpxr[n] = C_super * (-k_super2*kz_super[n]*Eyr[n] - ky_0*sigma[n]*Hpyr[n]);
		Hpxt[n] = C_sub   * ( k_sub2*kz_sub[n]*Eyt[n]     - ky_0*sigma[n]*Hpyt[n]);
	}

	/* Total incident energie */
	sumPzi = 0;
	for (n=Nmin_super; n<=Nmax_super; n++) {
/*		Pz_i = cabs(Exi[n+mid]*Hpyi[n+mid] - Eyi[n+mid]*Hpxi[n+mid]);*/
		Pz_i = fabs(creal(Exi[n+mid]*CONJ(Hpyi[n+mid]) - Eyi[n+mid]*CONJ(Hpxi[n+mid])));
		sumPzi += Pz_i;
	}

	/* Reflexion : efficiencies and directions */
	for (n=Nmin_super; n<=Nmax_super; n++) {
/*		Pz_r = cabs(Exr[n+mid]*Hpyr[n+mid] - Eyr[n+mid]*Hpxr[n+mid]);*/
		Pz_r = fabs(creal(Exr[n+mid]*CONJ(Hpyr[n+mid]) - Eyr[n+mid]*CONJ(Hpxr[n+mid])));
		eff->eff_r[n-Nmin_super] = Pz_r/sumPzi;
		eff->N_eff_r[n-Nmin_super] = (double) n;
		eff->theta_r[n-Nmin_super] =  SIGN(sigma[n+mid])*acos(kz_super[n+mid]/k_super) * 180.0/PI;
		eff->phi_r[n-Nmin_super] = asin(ky_0/(k_super*sin(eff->theta_r[n-Nmin_super]))) * 180.0/PI;
	}
/*printf("\nRe eff_R :\n");
SaveDbleTab2file (eff->eff_r, Nmax_super-Nmin_super+1, "stdout", " ");*/
	/* Transmission : efficiencies and directions */
	for (n=Nmin_sub; n<=Nmax_sub; n++) {
/*		Pz_t = cabs(Ext[n+mid]*Hpyt[n+mid] - Eyt[n+mid]*Hpxt[n+mid]);*/
		Pz_t = fabs(creal(Ext[n+mid]*CONJ(Hpyt[n+mid]) - Eyt[n+mid]*CONJ(Hpxt[n+mid])));
		eff->eff_t[n-Nmin_sub] = Pz_t/sumPzi;
		eff->N_eff_t[n-Nmin_sub] = (double) n;
		eff->theta_t[n-Nmin_sub] = SIGN(sigma[n+mid])*acos(kz_sub[n+mid]/k_sub) * 180.0/PI;
		eff->phi_t[n-Nmin_sub] = asin(ky_0/(k_sub*sin(eff->theta_t[n-Nmin_sub]))) * 180.0/PI;
	}
/*printf("\nRe eff_T :\n");
SaveDbleTab2file (eff->eff_t, Nmax_sub-Nmin_sub+1, "stdout", " ");*/
	
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
/*!	\fn		int md2D_amplitudes(complex *Ai, complex *Ar, complex *At, complex **S12, 
 *                               complex **S22, struct Param_struct *par)
 *
 *	\brief	Field amplitudes calculations
 */
/*-------------------------------------------------------------------------------------*/
int md2D_amplitudes(complex *Ai, complex *Ar, complex *At, complex **S12, complex **S22, struct Param_struct *par)
{
	int n;
	int vec_size = par->vec_size;
	complex *kz_super;
	kz_super = par->kz_super;
		
	/*	Vi = Ai*cexp(-I*kz_super*h) */
	for (n=0;n<=vec_size-1;n++){
		par->Vi[n]          = Ai[n]         *cexp(-I*kz_super[n]*par->h);
		par->Vi[n+vec_size] = Ai[n+vec_size]*cexp(-I*kz_super[n]*par->h);
	}
				
	/* Vt = S22*Vi */
	M_x_V (par->Vt, S22, par->Vi, 2*vec_size, 2*vec_size);

	/* Vr = S12*Vi */	
	M_x_V (par->Vr, S12, par->Vi, 2*vec_size, 2*vec_size);

	/* Ar = Vr*cexp(-I*kz_super*h) */
	for (n=0;n<=vec_size-1;n++){
		Ar[n]          = par->Vr[n]         *cexp(-I*kz_super[n]*par->h);
		Ar[n+vec_size] = par->Vr[n+vec_size]*cexp(-I*kz_super[n]*par->h);
	}
	
	/* At = Vt */
	for (n=0;n<=2*vec_size-1;n++){
		At[n] = par->Vt[n];
	}

	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn	int S_matrix(struct Param_struct *par)
 *
 *		\brief	S matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int S_matrix(struct Param_struct *par)
{

	if (par->verbosity) fprintf(stdout,"S-Matrix calculation\n");

	int i, j, nS;
	int vec_size = par->vec_size;
	int NS = par->NS;
	complex **S12, **S22, **T11, **T12, **T21, **T22, **Z, **T_tmp, **T_tmp2;
	S12 = par->S12;
	S22 = par->S22;
	T11 = par->T11;
	T12 = par->T12;
	T21 = par->T21;
	T22 = par->T22;

	Z      = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	T_tmp  = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	T_tmp2 = allocate_CplxMatrix(2*par->vec_size,2*par->vec_size);
	
	/* Initialisations */
	for (i=0; i<=2*vec_size-1; i++) {
		for (j=0; j<=2*vec_size-1; j++) {
			S12[i][j] = 0;
			S22[i][j] = 0;
		}
		S22[i][i] = 1;
	}

	/* Iterations */
	for (nS=0; nS<=NS-1; nS++) {

		/* T-Matrix calculation */
		(*par->T_Matrix)(T11, T12, T21, T22, nS, par);

		/* Z = inv(T11 + T12*S12) */
		invM(Z, add_M(T_tmp2, 
			T11, M_x_M(T_tmp,
				T12,S12,	2*vec_size, 2*vec_size), 2*vec_size, 2*vec_size), 2*vec_size);
		/* S12 = (T21 +T22*S12)*Z */
		M_x_M(S12,
			add_M(T_tmp2, T21, M_x_M(T_tmp,
					T22,S12, 2*vec_size, 2*vec_size), 2*vec_size, 2*vec_size),
			Z, 2*vec_size, 2*vec_size);
		/* S22 = S22*Z */
		M_equals(T_tmp,S22, 2*vec_size, 2*vec_size);
		M_x_M(S22,T_tmp,Z, 2*vec_size, 2*vec_size);

		/* if NEAR_FIELD, field components are saved in the Near_field_matrix */
/*		if (par->SAVE_NEAR_FIELD){
			md2D_save_near_field(nS, par);
		}
*/		
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



/*-------------------------------------------------------------------------------------*/
/*!	\fn	int dm_T_Matrix(complex **T11, complex **T12, complex **T21, complex **T22,
											int nS, struct Param_struct *par)
 *
 *		\brief T_matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int dm_T_Matrix(complex **T11, complex **T12, complex **T21, complex **T22,
				int nS, struct Param_struct *par)
{
	int j, n, np;
	int N = par->N;
	int NS = par->NS;
	int vec_size = par->vec_size;
	double h = par->h;
	complex k2__ky2, kz, k0h, k0h2, ky_sigma, k_super, kzk2, ky_02, ky_0;
		
	double hmin, hmax;
	complex *sigma, *kz0h, *Ex1, *Ey1, *Hpx1, *Hpy1, *Ex2, *Ey2, *Hpx2, *Hpy2, *VE_m2, *VE_p2, *VH_m2, *VH_p2, *F1, *F2, *V2;
	complex VE_m1, VE_p1, VH_m1, VH_p1;
	sigma = par->sigma;
	ky_0 = par->ky_0;
	k_super = par->k_super;

	F1 = (complex *) malloc(sizeof(complex)*4*par->vec_size);
	F2 = (complex *) malloc(sizeof(complex)*4*par->vec_size);
	V2 = (complex *) malloc(sizeof(complex)*4*par->vec_size);
	
	/* [F] vectors are defined by 	 	[V] vectors by
	[F] = |[Ex ]|								[V] = |[VE-]|
			|[Ey ]|								      |[VH-]|
			|[H'x]|								      |[VE+]|
			|[H'y]|  							      |[VH+]|  
	with H'=omega mu H											*/
	
	/* Pointer alignments, (for ease of use and readability) */	
	Ex1  = F1;
	Ey1  = F1 + vec_size;
	Hpx1 = F1 + 2*vec_size;
	Hpy1 = F1 + 3*vec_size;
	Ex2  = F2;
	Ey2  = F2 + vec_size;
	Hpx2 = F2 + 2*vec_size;
	Hpy2 = F2 + 3*vec_size;
	VE_m2 = V2; 
	VH_m2 = V2 + vec_size;
	VE_p2 = V2 + 2*vec_size;
	VH_p2 = V2 + 3*vec_size;
		
	for (n=0; n<=4*vec_size-1; n++){

		/* [F1] vector initialisation */
		for(j=0;j<=4*vec_size-1;j++){
			F1[j]  = 0;
		}

		/* [V1] vector initialisation */
		VE_m1 = 0; /* (since only one component is not nul at a time no need of an array) */
		VE_p1 = 0; 
		VH_m1 = 0; 
		VH_p1 = 0;
		
		/* [V2] vector initialisation */	
		for(j=0;j<=4*vec_size-1;j++){
			V2[j]  = 0;
		}

		/* k0h & kz0h values */
		if (nS==0){ /* 1st S-Matrix iteration : we are in the substrat */
			k0h  = par->k_sub;
			kz0h = par->kz_sub;
		}else{     /* Following iterations : we are in the superstrat */
			k0h  = par->k_super;
			kz0h = par->kz_super;
		}
		k0h2 = k0h*k0h;
		ky_02 = ky_0*ky_0;
		k2__ky2 = k0h2 - ky_02;
				
		/* Shooting method : determination of the initial [V1] vectors*/
		if(n>=0 && n<=vec_size-1){
			/* VE- = 1 */
			VE_m1 = 1;
			np = n;
		}else if(n>=vec_size && n<=2*vec_size-1){
			/* VH- = 1 */
			VH_m1 = 1;
			np = n - vec_size;
		}else if(n>=2*vec_size && n<=3*vec_size-1){
			/* VE+ = 1 */
			VE_p1 = 1;
			np = n - 2*vec_size;
		}else if(n>=3*vec_size && n<=4*vec_size-1){
			/* VH+ = 1 */
			VH_p1 = 1;
			np = n - 3*vec_size;
		}else{
			fprintf(stderr, "%s, line %d : ERROR, n out of normal range\n",__FILE__,__LINE__);
			np = 0;
			exit(EXIT_FAILURE);

		}
		/* Calculation of the corresponding [F1] vector components */
		Ex1 [np] = (1/k2__ky2)*(-ky_0*sigma[np]*VE_m1 - kz0h[np]*VH_m1 - ky_0*sigma[np]*VE_p1 + kz0h[np]*VH_p1);
		Ey1 [np] = VE_m1 + VE_p1;
		Hpx1[np] = (1/k2__ky2)*(kz0h[np]*k0h2*VE_m1 - ky_0*sigma[np]*VH_m1 - kz0h[np]*k0h2*VE_p1 - ky_0*sigma[np]*VH_p1);
		Hpy1[np] = VH_m1 + VH_p1;
		
  		/* Intégration des grandeurs pour une couche */
		hmin = h*nS/NS;
		hmax = h*(nS+1)/NS;
		
		/* (alternative RK4 algorithm */
		/*eq_diff((double *)F1,  (double *)F2, 8*vec_size, hmin, hmax, par->Nstep_S, fun, (void *) par);
*/
		/* Calling the RK4 algorithm */
		ode_solve((double *)F1, (double *)F2, 8*vec_size, hmin, hmax, par->Nstep_S, fun, (void *) par);

		/* Translating resulting [F2] vector in [V2] vector */
		k2__ky2 = k_super*k_super - ky_0*ky_0;
		for(j=0;j<=vec_size-1;j++){
/*			sigma = (j-half_vec_size)*delta_sigma + sigma0;
*/			ky_sigma = ky_0*sigma[j];
/*			kz = csqrt(k0h2- sigma^2 - ky_0^2);
*//*		kzk2 = kz0h[j]*k_super*k_super;
*/			kz = par->kz_super[j];
			kzk2 = kz*k_super*k_super;
			VE_m2[j] = 0.5*(Ey2[j] + (k2__ky2/kzk2)*Hpx2[j] + (ky_sigma/kzk2)*Hpy2[j]);
			VH_m2[j] = 0.5*(-(k2__ky2/kz)*Ex2[j] - (ky_sigma/kz)*Ey2[j] + Hpy2[j]);
			VE_p2[j] = 0.5*(Ey2[j] - (k2__ky2/kzk2)*Hpx2[j] - (ky_sigma/kzk2)*Hpy2[j]);
			VH_p2[j] = 0.5*((k2__ky2/kz)*Ex2[j] + (ky_sigma/kz)*Ey2[j] + Hpy2[j]);
		}
		
		/* Constructing T matrices */
		if(n>=0 && n<=2*vec_size-1){
			np = n;
			for(j=0;j<=2*vec_size-1;j++){
				T11[j][np] = VE_m2[j]; /* notice that VE_m2[vec_size+j] corresponds to VH_m2[j] ...*/
				T21[j][np] = VE_p2[j]; /* VE_p2[vec_size+j] == VH_p2[j] ...*/
			}
		}else{
			np = n - 2*vec_size;
			for(j=0;j<=2*vec_size-1;j++){
				T12[j][np] = VE_m2[j]; /* notice that VE_m2[vec_size+j] corresponds to VH_m2[j] ...*/
				T22[j][np] = VE_p2[j]; /* VE_p2[vec_size+j] == VH_p2[j] ...*/
			}
		}

		/*Affichage du temps restant à l'écran */
		md2D_affichTemps(n,N,nS,NS,par->ni,par->Ni,par);
	}
	
	free(F1);
	free(F2);
	free(V2);
	
	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn	int rcwa_T_Matrix(complex **T11, complex **T12, complex **T21, complex **T22,
					int nS, struct Param_struct *par)
 *
 *		\brief T_matrix calculation for RCWA
 */
/*-------------------------------------------------------------------------------------*/
int rcwa_T_Matrix(complex **T11, complex **T12, complex **T21, complex **T22,
				int nS, struct Param_struct *par)
{
	int i,j;
	int vec_size = par->vec_size;
	
	double hmin, hmax, z, Delta_z;
	complex **EigVectors = par->EigVectors;
	complex **invEigVec = par->invEigVec;
	complex **Psi;

	/* [F] vectors are defined by 	 	[V] vectors by
	[F] = |[Ex ]|								[V] = |[VE-]|
			|[Ey ]|								      |[VH-]|
			|[H'x]|								      |[VE+]|
			|[H'y]|  							      |[VH+]|  
	with H'=omega mu H											*/

	/* Psi matrix */
	if (nS==0){ /* 1st S-Matrix iteration : we are in the substrat */
		Psi = par->Psi_sub;
	}else{     /* Following iterations : we are in the superstrat */
		Psi = par->Psi_super;
	}

	/* z of the considered T matirx slice */
	hmin = par->h*nS/par->NS;
	hmax = par->h*(nS+1)/par->NS;
	if (par->tab_NS_ENABLED){
		hmin = par->tab_NS[nS];
		hmax = par->tab_NS[nS+1];
	}
	Delta_z = hmax - hmin;
	z = (hmax + hmin)/2;

	/* M matrix calculation */
	rcwa_M_Matrix(par->M, z, par);
/*printf("\nIm(par->M) :\n");
SaveMatrix2file (par->M, 4*par->vec_size, 4*par->vec_size, "Im", "stdout");
*/
	/* Diagonalisation of M */
	eigen_values(par->M, par->eig_values, EigVectors, par->eig_buffer, 4*vec_size);

	/* Solution of the diagonalized system */
	for (i=0;i<=4*vec_size-1;i++){
		par->M_sol[i] = cexp(par->eig_values[i]*Delta_z);
	}
	
	/* T matrix : T = inv(Psi_super) * EigVec * M_diag * inv(EigVec) * Psi */
	/* invEigVec = inv(EigVectors) */
	invM(invEigVec, EigVectors, 4*vec_size);

	/* invVec_Psi = invEigVec * Psi*/
	M_x_M(par->invVec_Psi, invEigVec, Psi, 4*vec_size, 4*vec_size);

	/* M_invVec_Psi = M_sol * invVec_Psi */
	for (i=0;i<=4*vec_size-1;i++){
		for (j=0;j<=4*vec_size-1;j++){
			par->M_invVec_Psi[i][j] = par->M_sol[i] * par->invVec_Psi[i][j];
		}
	}
	/* T = invPsi_super * EigVectors * M_invVec_Psi */
	M_x_M(par->T,
			par->invPsi_super, M_x_M(par->M_buffer_4vecsize, 
										EigVectors, par->M_invVec_Psi, 4*vec_size, 4*vec_size), 4*vec_size, 4*vec_size);

	/*Affichage du temps restant à l'écran */
	md2D_affichTemps(par->N,par->N,nS,par->NS,par->ni,par->Ni,par);
	
	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int fun (double z, const double *F_real, double *dF_real, void *param_void)
 *
 *		\brief	Calculation of dF, function called by the ODE solver
 */
/*-------------------------------------------------------------------------------------*/
int fun (double z, const double *F_real, double *dF_real, void *param_void)
{
	int j, vec_size;
	struct Param_struct *par = (struct Param_struct *) param_void;
	complex sigm, *sigma;
	complex ky_0 = par->ky_0;	
	complex ky_02 = ky_0*ky_0;
	complex *F, *dF;
	complex *Ex, *Ey, *Hpx, *Hpy, *dEx, *dEy, *dHpx, *dHpy;
	complex **Qxx, **Qyy, **Qxz, **Qzz, **Qzz_1;
	complex *QxzEx, *Qzz_1QxzEx, *Qzz_1Hpx, *V_tmp1, *sigmaHpy; 
	complex *ky0Qzz_1Hpx, *Qzz_1sigmaHpy, *QxxEx, *QyyEy, *QxzVtmp1;

	vec_size = par->vec_size;
	sigma = par->sigma;
	Qxx = par->Qxx;
	Qyy = par->Qyy;
	Qxz = par->Qxz;
	Qzz = par->Qzz;
	Qzz_1 = par->Qzz_1;

	QxzEx = par->QxzEx;
	Qzz_1QxzEx = par->Qzz_1QxzEx;
	Qzz_1Hpx = par->Qzz_1Hpx;
	V_tmp1 = par->V_tmp1;
	sigmaHpy = par->sigmaHpy;
	ky0Qzz_1Hpx = par->ky0Qzz_1Hpx;
	Qzz_1sigmaHpy = par->Qzz_1sigmaHpy;
	QxxEx = par->QxxEx;
	QyyEy = par->QyyEy;
	QxzVtmp1 = par->QxzVtmp1;
	
	/* Cast F_real & dF_real into complex 
	   - Rem 1 : most ODE solver algorithms need datas in 'double' format
		          and it is more convenient and readable to make calculations in complex
		- Rem 2 : this simple cast does the job since complex are stored as 
		          [complex] = [real,imag] in memory, thus the solver actually see a system 
					 twice bigger composed of Re() & Im() equations of the complex system */
	F  = (complex *) F_real;
	dF = (complex *) dF_real;

	/* Pointer alignments, (for ease of use and readability) */	
	Ex  = F;
	Ey  = F + vec_size;
	Hpx = F + 2*vec_size;
	Hpy = F + 3*vec_size;
	dEx  = dF;
	dEy  = dF + vec_size;
	dHpx = dF + 2*vec_size;
	dHpy = dF + 3*vec_size;

	/* Toeplitz matrices calculations */
	md2D_QMatrix(z, Qxx, Qyy, Qxz, Qzz, Qzz_1, par);

	/* Stored Toeplitz matrices */
	/* ESSAYER AVEC POINTEUR DE FONCTION POUR GAGNER DU TEMPS*/
/*	if (par->STOCKER_TF) {
		FFT_k2_et_invk2_stockee(z, &TF_k2, &TF_invk2, par);
	}else{
		TF_k2 = FFT_k2_directe(z, TF_k2, par);
		TF_invk2 = FFT_invk2_directe(z, TF_invk2, par);
	}*/
			
	/* Precalculating some datas used in the differential system */
	M_x_V(QxzEx, Qxz, Ex, vec_size, vec_size); 				/* QxzEx = Qxz x Ex */
	M_x_V(Qzz_1QxzEx, Qzz_1, QxzEx, vec_size, vec_size);  /* Qzz_1QxzEx = Qzz_1 x QxzEx */
	M_x_V(Qzz_1Hpx, Qzz_1, Hpx, vec_size, vec_size); 		/* Qzz_1Hpx = Qzz_1 x Hpx */
	for(j=0;j<=vec_size-1;j++){
		sigmaHpy[j] = sigma[j]*Hpy[j];
		ky0Qzz_1Hpx[j] = ky_0*Qzz_1Hpx[j];
	}
	M_x_V(Qzz_1sigmaHpy, Qzz_1, sigmaHpy, vec_size, vec_size);  /* Qzz_1sigmaHpy = Qzz_1 x sigmaHpy */
	M_x_V(QxxEx, Qxx, Ex, vec_size, vec_size); 						/* QxxEx = Qxx x Ex */
	M_x_V(QyyEy, Qyy, Ey, vec_size, vec_size); 						/* QyyEy = Qyy x Ey */
	for(j=0;j<=vec_size-1;j++){
		V_tmp1[j] = -Qzz_1QxzEx[j] + ky0Qzz_1Hpx[j] - Qzz_1sigmaHpy[j];
	}
	M_x_V(QxzVtmp1, Qxz, V_tmp1, vec_size, vec_size);  			/* QxzVtmp1 = Qxz x V_tmp1 */
	/* Differential system writting */
	for(j=0;j<=vec_size-1;j++){
		sigm = sigma[j];
		dEx[j]  = I*(sigm*V_tmp1[j] + Hpy[j]);
		dEy[j]  = I*(ky_0*V_tmp1[j] - Hpx[j]);
		dHpx[j] = I*(-ky_0*sigm*Ex[j] + sigm*sigm*Ey[j] - QyyEy[j]);
		dHpy[j] = I*(QxzVtmp1[j] +QxxEx[j] - ky_02*Ex[j] + ky_0*sigm*Ey[j]);
	}

	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_QMatrix(double z, complex **Qxx, complex **Qyy, complex **Qxz, complex **Qzz, complex **Qzz_1, struct Param_struct *par)
 *
 *		\brief	Q matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int md2D_QMatrix(double z, complex **Qxx, complex **Qyy, complex **Qxz, complex **Qzz, complex **Qzz_1, struct Param_struct *par)
{
	int i,j,var_tmp01;
	fftw_plan plan_TFk2, plan_TFinvk2, plan_TFNx2, plan_TFNz2, plan_TFNxNz;
	int N_x = par->N_x;
	int N_tf = 2*par->N;
	int vec_size = par->vec_size;
	double coefnorm;
	complex *tmp_k2, *tmp_invk2, *tmp_Nx2, *tmp_Nz2, *tmp_NxNz;
	tmp_k2    = par->tmp_tf_k2;
	tmp_invk2 = par->tmp_tf_invk2;
	tmp_Nx2   = par->tmp_tf_Nx2;
	tmp_Nz2   = par->tmp_tf_Nz2;
	tmp_NxNz  = par->tmp_tf_NxNz;

	complex *TF_k2, *TF_invk2, *TF_Nx2, *TF_Nz2, *TF_NxNz;
	TF_k2    = par->TF_k2;
	TF_invk2 = par->TF_invk2;
	TF_Nx2   = par->TF_Nx2;
	TF_Nz2   = par->TF_Nz2;
	TF_NxNz  = par->TF_NxNz;
	
	complex **Toep_k2, **Toep_invk2, **Toep_Nx2, **Toep_NxNz, **Toep_Nz2, **invToep_invk2, **M_tmp1, **M_tmp2;
	Toep_k2       = par->Toep_k2;
	Toep_invk2    = par->Toep_invk2;
	Toep_Nx2      = par->Toep_Nx2;
	Toep_NxNz     = par->Toep_NxNz;
	Toep_Nz2      = par->Toep_Nz2;
	invToep_invk2 = par->invToep_invk2;
	M_tmp1 = par->M_tmp1;
	M_tmp2 = par->M_tmp2;
		
	/*--- Calculations of the Fourier transforms of the necessary 'grandeurs' ---*/

	/* Calculating the values in the direct space, using
	function pointors to adapt to surface definition type */
	(*par->k_2)(par, par->k2, z);
	(*par->invk_2)(par, par->invk2, z);
/*printf("\ninvk2_1D : \n");SaveCplxTab2file (par->invk2, par->N_x, "Re", "stdout"," ");
*/	(*par->Normal_function)(par, par->Nx2, par->NxNz, par->Nz2, z);
	
	/* Creating the 'plans' for the FFTW */
	plan_TFk2    = fftw_plan_dft_1d(N_x, (fftw_complex *)par->k2,    (fftw_complex *)tmp_k2,    FFTW_FORWARD, FFTW_ESTIMATE);	
	plan_TFinvk2 = fftw_plan_dft_1d(N_x, (fftw_complex *)par->invk2, (fftw_complex *)tmp_invk2, FFTW_FORWARD, FFTW_ESTIMATE);	
	plan_TFNx2   = fftw_plan_dft_1d(N_x, (fftw_complex *)par->Nx2,   (fftw_complex *)tmp_Nx2,   FFTW_FORWARD, FFTW_ESTIMATE);	
	plan_TFNz2   = fftw_plan_dft_1d(N_x, (fftw_complex *)par->Nz2,   (fftw_complex *)tmp_Nz2,   FFTW_FORWARD, FFTW_ESTIMATE);	
	plan_TFNxNz  = fftw_plan_dft_1d(N_x, (fftw_complex *)par->NxNz,  (fftw_complex *)tmp_NxNz,  FFTW_FORWARD, FFTW_ESTIMATE);	

	/* Calculating the FFT */
	/* (NOTICE : Real DFT could be used for Nx2, Nz2, etc, which would slightly increase speed but also code complexity) */	
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
/*!	\fn		int rcwa_M_Matrix(complex **M, double z, struct Param_struct *par)
 *
 *		\brief	M matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int rcwa_M_Matrix(complex **M, double z, struct Param_struct *par)
{
	int i,j,var_tmp01;
	fftw_plan plan_TFk2, plan_TFinvk2;
	int N_x = par->N_x;
	int N_tf = 2*par->N;
	int vec_size = par->vec_size;
	double coefnorm;
	complex *tmp_k2, *tmp_invk2;
	tmp_k2    = par->tmp_tf_k2;
	tmp_invk2 = par->tmp_tf_invk2;
	
	complex ky_0 = par->ky_0;
	complex ky_02 = ky_0*ky_0;
	complex *sigma = par->sigma;
	complex *TF_k2    = par->TF_k2;
	complex *TF_invk2 = par->TF_invk2;	
	complex **Toep_k2       = par->Toep_k2;
	complex **Toep_invk2    = par->Toep_invk2;
	complex **invToep_invk2 = par->invToep_invk2;
	complex **invToep_k2    = par->invToep_k2;

	/*--- Calculations of the Fourier transforms of the necessary 'grandeurs' ---*/

	/* Calculating the values in the direct space, using
	function pointors to adapt to surface definition type */
	(*par->k_2)(par, par->k2, z);
	(*par->invk_2)(par, par->invk2, z);
	
	/* Creating the 'plans' for the FFTW */
	plan_TFk2    = fftw_plan_dft_1d(N_x, (fftw_complex *)par->k2,    (fftw_complex *)tmp_k2,    FFTW_FORWARD, FFTW_ESTIMATE);	
	plan_TFinvk2 = fftw_plan_dft_1d(N_x, (fftw_complex *)par->invk2, (fftw_complex *)tmp_invk2, FFTW_FORWARD, FFTW_ESTIMATE);	

	/* Calculating the FFT */
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
	/* Inverting Toep_invk2 & Toep_k2 */
	invM(invToep_invk2, Toep_invk2, vec_size);
	invM(invToep_k2, Toep_k2, vec_size);

	for(i=0;i<=4*vec_size-1;i++){
		for(j=0;j<=4*vec_size-1;j++){
			M[i][j] = 0;
		}
	}
	for(i=0;i<=vec_size-1;i++){
		/* M31 = -sigma * ky_0 */
		M[i+2*vec_size][i] = -I*ky_0*sigma[i];
		/* M42 = sigma * ky_0 */
		M[i+3*vec_size][i+  vec_size] = I*ky_0*sigma[i];
		for(j=0;j<=vec_size-1;j++){
			/* M13 = sigma * K2_1 * ky_0 */
			M[i][j+2*vec_size] =  I*sigma[i] * invToep_k2[i][j] * ky_0;
			/* M14 = -sigma * K2_1 * sigma + Id*/
			M[i][j+3*vec_size] = -I*sigma[i] * invToep_k2[i][j] * sigma[j];
			/* M23 = K2_1 * ky_02 - Id */
			M[i+vec_size][j+2*vec_size] = I*invToep_k2[i][j] * ky_02;
			/* M24 = -ky_0 * K2_1 * sigma */
			M[i+vec_size][j+3*vec_size] = -I*ky_0 * invToep_k2[i][j] * sigma[j];
			/* M32 = sigma^2 - K2 */
			M[i+2*vec_size][j+vec_size] = -I*Toep_k2[i][j];
			/* M41 = inv_K2_1 - ky02 * Id */
			M[i+3*vec_size][j] = I*invToep_invk2[i][j];
		}
		/* M14 = -sigma * K2_1 * sigma + Id */
		M[i][i+3*vec_size] +=  I;
		/* M23 = K2_1 * ky_02 - Id */
		M[i+vec_size][i+2*vec_size] -= I;
		/* M32 = sigma^2 - K2 */
		M[i+2*vec_size][i+vec_size] += I*sigma[i]*sigma[i];
		/* M41 = inv_K2_1 - ky02 * Id */
		M[i+3*vec_size][i] -= I*ky_02;
	}


	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md2D_PsiMatrix(complex **Psi, complex k, complex *kz, struct Param_struct *par)
 *
 *		\brief	Psi matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int md2D_PsiMatrix(complex **Psi, complex k, complex *kz, struct Param_struct *par)
{
	int i,j;
	int vec_size = par->vec_size;
	complex p, qe, qh;
	complex ky_0 = par->ky_0;
	complex k2 = k*k;
	complex ky02 = ky_0*ky_0;

	for(i=0;i<=4*vec_size-1;i++){
		for(j=0;j<=4*vec_size-1;j++){
			Psi[i][j] = 0;
		}
	}
	/*	Psi = [ p    qe   p   -qe   ; ...
   	        Id   zero Id   zero ; ...
      	     qh   p   -qh   p    ; ...
         	  zero Id   zero Id  ];
  */
	for(i=0;i<=vec_size-1;i++){
   	p  = -ky_0*par->sigma[i]/(k2-ky02);
   	qe = -kz[i]/(k2-ky02);
   	qh = k2*kz[i]/(k2-ky02);
	
		Psi[i][i]            =  p;
		Psi[i][i+  vec_size] =  qe;
		Psi[i][i+2*vec_size] =  p;
		Psi[i][i+3*vec_size] = -qe;
			
		Psi[i+vec_size][i]            = 1;
		Psi[i+vec_size][i+2*vec_size] = 1;

		Psi[i+2*vec_size][i]            =  qh;
		Psi[i+2*vec_size][i+  vec_size] =  p;
		Psi[i+2*vec_size][i+2*vec_size] = -qh;
		Psi[i+2*vec_size][i+3*vec_size] =  p;

		Psi[i+3*vec_size][i+  vec_size] = 1;
		Psi[i+3*vec_size][i+3*vec_size] = 1;
	}

	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int k2_H_X(struct Param_struct *par, complex *k2_1D, double z)
 *
 *	\brief	Détermine le tableau de complexes k^2(x) pour un z donné
 *
 */
/*-------------------------------------------------------------------------------------*/
int k2_H_X(struct Param_struct *par, complex *k2_1D, double z)
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
/*!	\fn			
 *
 *	\brief	Détermine le tableau de complexes k^2(x) pour un z donné, pour un multicouches
 */
/*-------------------------------------------------------------------------------------*/
int k2_MULTI(struct Param_struct *par, complex *k2_1D, double z)
{
	
	int nx, n_layer=0;
	
	for (nx=0; nx<=par->N_x-1; nx++){
		do{
			if (z <= par->profil[n_layer][nx]){
				if (z >= par->profil[n_layer+1][nx]){
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
/*!	\fn		k2_N_XYZ(struct Param_struct *par, complex *k2_1D, double z)	
 *
 *		\brief	Détermine le tableau de complexes k^2(x) pour un z donné
 */
/*-------------------------------------------------------------------------------------*/
int k2_N_XYZ(struct Param_struct *par, complex *k2_1D, double z)
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
/*!	\fn	invk2_N_XYZ(struct Param_struct *par, complex *invk2_1D, double z)	
 *
 *	\brief	Détermine le tableau de complexes invk^2(x) pour un z donné
 */
/*-------------------------------------------------------------------------------------*/
int invk2_N_XYZ(struct Param_struct *par, complex *invk2_1D, double z)
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
/*!	\fn		int invk_2(complex *invk2_1D, double *profil, struct Param_struct *par, double z)
 *
 *	\brief	Détermine le tableau de complexes 1/k^2(x) pour un z donné
 *
 *	\todo	Prend pour l'instant en compte seulement un profil de type h(x)\n
 *			Doit être plus polyvalent : accepter aussi les profils de type n(x,z)
 */
/*-------------------------------------------------------------------------------------*/
int invk2_H_X(struct Param_struct *par, complex *invk2_1D, double z)
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
/*!	\fn			
 *
 *	\brief	Détermine le tableau de complexes 1/k^2(x) pour un z donné, pour un multicouches
 */
/*-------------------------------------------------------------------------------------*/
int invk2_MULTI(struct Param_struct *par, complex *invk2_1D, double z)
{
	
	int nx, n_layer=0;
	
	for (nx=0; nx<=par->N_x-1; nx++){
		do{
			if (z <= par->profil[n_layer][nx]){
				if (z >= par->profil[n_layer+1][nx]){
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
/*!	\fn		int Normal_H_X(struct Param_struct *par, complex *k2_1D, double z)
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
int Normal_H_X(struct Param_struct *par, complex *Nx2, complex *NxNz, complex *Nz2, double z)
{
	int i;
	int Nx = par->N_x; 	/* Caution : this 'Nx', corresponds to the number of points in x              */
								/* while 'Nx2', corresponds to the x component^2 of the normal to the surface */
	double dhdx, norm_x, norm_z;
	double two_dx = 2*par->L/Nx;

	double *profil = par->profil[0];

	if (par->HX_Normal_CALCULATED == 1){
		Nx2  = par->Nx2;
		NxNz = par->NxNz;
		Nz2  = par->Nz2;
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
/*!	\fn		complex *FFT_k2_directe(double z, complex *TF_k2, struct Param_struct *par)
 *
 *		\brief	
 */
/*-------------------------------------------------------------------------------------*/
complex *FFT_k2_directe(double z, complex *TF_k2, struct Param_struct *par)
{

	int i;
	complex *tmp;
	fftw_plan plan_TFk2;

	int N_x = par->N_x;
	int N_tf = 2*par->N;
	double coefnorm = 1.0/N_x;

	/* Calcul de k^2(x) à z fixé, à partir du profil */
	(*par->k_2)(par, par->k2, z);

	/* Calcul de la TF de k2, avec N_x points */
	
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!   ALLOUER A L'EXTERIEUR   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
	tmp = (complex *) malloc(sizeof(complex) * N_x);
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
	
	return TF_k2;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn		complex *FFT_invk2_directe(double z, complex *TF_invk2, struct Param_struct *par)
 *
 *	\brief	Calcule la TF de 1/k^2(x) pour un z donné et la tronque entre -N et +N
 */
/*-------------------------------------------------------------------------------------*/
complex *FFT_invk2_directe(double z, complex *TF_invk2, struct Param_struct *par)
{

	int i;
	complex *tmp;
	fftw_plan plan_TFk2;

	int N_x = par->N_x;
	int N_tf = 2*par->N;
	double coefnorm = 1.0/N_x;

	/* Calcul de 1 / k^2(x) à z fixé, à partir du profil */
	(*par->invk_2)(par, par->invk2, z);

	/* Calcul de la TF de invk2, avec N_x points */
	tmp = (complex *) malloc(sizeof(complex) * N_x); /*TODO : ALLOUER A L'EXTERIEUR */
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
/*!	\fn		int *FFT_k2_stockee(double z, complex **TF_k2, struct Param_struct *par)
 *
 *	\brief	Calcule la TF de k^2(x) pour un z donné. Si le calcul à déja été fait, FFT_k2_stockee \n
 *          retourne l'adresse de la ligne du tableau tab_TF_k2 où les valeurs ont été stockées \n
 *          sinon le calcul est effectué par FFT_k2_directe et stocké dans tab_TF_k2.
 */
/*-------------------------------------------------------------------------------------*/
int FFT_k2_stockee(double z, complex **ptTF_k2, struct Param_struct *par)
{
	long int **z2n = par->z2n;
	
	/* z est convertit en long int afin de permettre la comparaison == entre deux z */
	long int z_int = (long int)((z/par->h)*DBLE_CMP_EXIGEANCE);

	/* Recherche de z dans z2n pour voir si le calcul à déja été fait */
	long int *adr = cherche(z_int, z2n[0], par->N_z2n);
	if (adr != NULL) {
	/* z trouvé : On renvoie la ligne n indiquée par z2n du tableau tab_FFT_k2 */
		*ptTF_k2 = par->tab_TF_k2[z2n[1][adr-z2n[0]]];
		return 0;
	}

	/* z non touvé : On calcule la TF qu'on ajoute au tableau tab_FFT_k2 et on met à jour z2n */
	/* Réallocation éventuelle de mémoire */
	if (par->N_z2n+1 > par->TAILLE_z2n) {
		/* ROUJOUTER : VERIF que pas trop de memoire allouée */
printf("REALLOCATION DE MEMOIRE\n");	
		par->TAILLE_z2n += par->BLOC_TAILLE_z2n;	
		par->tab_TF_k2 = reallocate_CplxMatrix(par->tab_TF_k2,(par->TAILLE_z2n + par->BLOC_TAILLE_z2n),4*par->N+1);
		z2n[0] = (long int *) realloc(z2n[0], sizeof(long int)*(par->TAILLE_z2n + par->BLOC_TAILLE_z2n));
		z2n[1] = (long int *) realloc(z2n[1], sizeof(long int)*(par->TAILLE_z2n + par->BLOC_TAILLE_z2n));
	}
	/* Mise à jour de z2n : insertion de z et de n dans le tableau décroissant en z */
	int k = par->N_z2n;
	par->N_z2n++; 
	while(z_int < z2n[0][k-1] && k > 0){
		z2n[0][k] = z2n[0][k-1];
		z2n[1][k] = z2n[1][k-1];
		k--;
	}
	z2n[0][k] = z_int;
	z2n[1][k] = par->N_z2n-1;
	/* Calcul de la TF de k2(z) et mise à jour de tab_FFT_k2 */
	FFT_k2_directe(z, par->tab_TF_k2[par->N_z2n-1], par);

	*ptTF_k2 = par->tab_TF_k2[par->N_z2n-1];
	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int FFT_k2_et_invk2_stockee(double z, complex **ptTF_k2, complex **ptTF_invk2, struct Param_struct *par)
 *
 *	\brief	Calcule la TF de k^2(x) et de 1 / k^2(x) pour un z donné.\n
 *          Voir FFT_k2_stockee() pour + de détails
 */
/*-------------------------------------------------------------------------------------*/
int FFT_k2_et_invk2_stockee(double z, complex **ptTF_k2, complex **ptTF_invk2, struct Param_struct *par)
{
	long int **z2n = par->z2n;
	
	/* z est convertit en long int afin de permettre la comparaison == entre deux z */
	long int z_int = (long int)((z/par->h)*DBLE_CMP_EXIGEANCE);

	/* Recherche de z dans z2n pour voir si le calcul à déja été fait */
	long int *adr = cherche(z_int, z2n[0], par->N_z2n);
	if (adr != NULL) {
	/* z trouvé : On renvoie la ligne n indiquée par z2n du tableau tab_FFT_k2 */
		*ptTF_k2    = par->tab_TF_k2[z2n[1][adr-z2n[0]]];
		*ptTF_invk2 = par->tab_TF_invk2[z2n[1][adr-z2n[0]]];
		return 0;		
	}

	/* z non touvé : On calcule la TF qu'on ajoute au tableau tab_FFT_k2 et on met à jour z2n */
	/* Réallocation éventuelle de mémoire */
	if (par->N_z2n+1 > par->TAILLE_z2n) {
		/* ROUJOUTER : VERIF que pas trop de memoire allouée */
		printf("REALLOCATION DE MEMOIRE pour tab_TF_invk2\n");	
		par->TAILLE_z2n += par->BLOC_TAILLE_z2n;	
		par->tab_TF_k2    = reallocate_CplxMatrix(par->tab_TF_k2,(par->TAILLE_z2n + par->BLOC_TAILLE_z2n),4*par->N+1);
		par->tab_TF_invk2 = reallocate_CplxMatrix(par->tab_TF_invk2,(par->TAILLE_z2n + par->BLOC_TAILLE_z2n),4*par->N+1);
		z2n[0] = (long int *) realloc(z2n[0], sizeof(long int)*(par->TAILLE_z2n + par->BLOC_TAILLE_z2n));
		z2n[1] = (long int *) realloc(z2n[1], sizeof(long int)*(par->TAILLE_z2n + par->BLOC_TAILLE_z2n));
	}
	/* Mise à jour de z2n : insertion de z et de n dans le tableau décroissant en z */
	int k = par->N_z2n;
	par->N_z2n++; 
	while(z_int < z2n[0][k-1] && k > 0){
		z2n[0][k] = z2n[0][k-1];
		z2n[1][k] = z2n[1][k-1];
		k--;
	}
	z2n[0][k] = z_int;
	z2n[1][k] = par->N_z2n-1;
	/* Calcul de la TF de k2(z) et mise à jour de tab_FFT_k2 */
	FFT_k2_directe(z, par->tab_TF_k2[par->N_z2n-1], par);
	FFT_invk2_directe(z, par->tab_TF_invk2[par->N_z2n-1], par);

	/* Renvoie des résultats */
	*ptTF_k2    = par->tab_TF_k2[par->N_z2n-1];
	*ptTF_invk2 = par->tab_TF_invk2[par->N_z2n-1];
	return 0;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn		long int *cherche(long int z, long int *tab0, int N)
 *
 *	\brief	Recherche dichotomique de l'élément z dans un tableau tab0 de taille N \n
 *			classé dans l'ORDRE CROISSANT. 
 *
 * 	\return	L'adresse correspondant à l'élément trouvé ou NULL si l'élément n'est pas présent
 */
/*-------------------------------------------------------------------------------------*/
long int *cherche(long int z, long int *tab0, int N)
{

	if (N<=1) {
		if (z==tab0[0]) return tab0;
		else            return NULL;
	}
	
	if (tab0[N>>1] > z) return(cherche(z, tab0, N>>1));
	else                return(cherche(z, tab0+(N>>1), N-(N>>1)));
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
	}else /*if (par->verbosity >= 1)*/{
		int i, n_total;
		time(&par->last_time);
		par->last_clock = clock();
		float t_ecoule = difftime(par->last_time,par->time0);
		float t_total; 

		n_total = 4*(2*N+1);
		if (!strcmp(par->calcul_method,"DM")){
			t_total = t_ecoule*( Ni*NS*n_total)/(ni*NS*n_total+nS*n_total+n);
		}else if (!strcmp(par->calcul_method,"RCWA")){
			t_total = t_ecoule*NS/(nS+1);
		}else{
			fprintf(stderr, "%s, line %d : ERROR, unknown calculation method (\"%s\")\n",__FILE__,__LINE__,par->calcul_method);
			exit(EXIT_FAILURE);
		}
		
		float t_restant = t_total - t_ecoule;
		int pourcent = ROUND(100.0*t_ecoule/t_total);
				
		fprintf(stdout,"\r");
		if (par->Ni > 1){
			 fprintf(stdout,"%3d %%, i = %d° [%ds ", pourcent,ROUND(par->theta_i*180/PI),ROUND(t_ecoule));
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
		sprintf(filename,"S_TM_%s_h%f.txt",par->nom_profil,h_partial);
	}else{
		sprintf(filename,"S_TE_%s_h%f.txt",par->nom_profil,h_partial);
	}	
	fp = fopen(filename, "w");

	/* Ecriture de certains parametres */
	fprintf(fp,"mat_S_name = %s_h%f\n",par->nom_profil,h_partial);
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


/*-------------------------------------------------------------------------------------*/
/*!	\fn	int md2D_make_tab_S_steps(struct Param_struct* par)
 *
 *	\brief	N steps tab making, when some S steps are imposed
 */
/*-------------------------------------------------------------------------------------*/
int md2D_make_tab_S_steps(struct Param_struct* par){

	int NSmin, num, ns, ms;
	double eps = 1e-10;	
	/* read imposed_S_steps */
	
	/**/
	NSmin = CEIL(10*par->h/par->lambda);

	num=0;
	for (ns=0;ns<=NSmin;ns++){
		par->tab_NS[num] = ns*(par->h/NSmin);
		num++;
		for (ms=0;ms<=par->N_imposed_S_steps-1;ms++){
			if((par->tab_imposed_S_steps[ms] > ns*(par->h/NSmin)+eps) && (par->tab_imposed_S_steps[ms] < (ns+1)*(par->h/NSmin)-eps)){
					par->tab_NS[num] = par->tab_imposed_S_steps[ms];
					num++;
			}
		}
	}
	par->NS=num-1;

	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn	int md2D_save_near_field(int nS, struct Param_struct* par)
 *
 *	\brief	field components saving for near field maping
 */
/*-------------------------------------------------------------------------------------*/
int md2D_save_near_field(int nS, struct Param_struct* par)
{

	
		
	return 0;
}



