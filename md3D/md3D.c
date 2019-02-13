/*!	\file		md3D.c
 *
 * 	\brief		Differential method \n
 *				3-Dimensions ! \n
 *				S-Matrices algorithm \n
 * 				FFF algorithm \n
 * 
 *
 *	\date		nov 2007
 *	\version	09.10.15
 *	\authors	Laurent ARNAUD
 */


#include "md3D.h"


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md3D_incident_field(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *		\brief	Incident field determination
 */
/*-------------------------------------------------------------------------------------*/
int md3D_incident_field(struct Param_struct *par, struct Efficacites_struct *eff)
{
	int n;
	COMPLEX Eyi, Hpyi;
	double phi_i = par->phi_i;
	double psi = par->psi;
	double theta_i = par->theta_i;
	COMPLEX k_super = par->k_super;
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
		fprintf(stderr, "%s, line %d : ERROR, \"%s\" : unsupported i_field_mode type\n",__FILE__,__LINE__,par->i_field_mode);
		exit(EXIT_FAILURE);
	}

	
	return 0;	
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md3D_propagativ_limits(struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *	\brief		Calculation of Propagatives modes limits
 */
/*-------------------------------------------------------------------------------------*/
int md3D_propagativ_limits(struct Param_struct *par, struct Efficacites_struct *eff)
{
	int Nx = par->Nx;
	int Ny = par->Ny;
	double Delta_sigma_x = par->Delta_sigma_x;
	double Delta_sigma_y = par->Delta_sigma_y;
	COMPLEX k_super2 = par->k_super*par->k_super;
	COMPLEX k_sub2 = par->k_sub*par->k_sub;
	
	double sigma_x0 = par->sigma_x0;
	double sigma_y0 = par->sigma_y0;
	COMPLEX sigma_y02 = par->sigma_y0*par->sigma_y0;
	COMPLEX sigma_x02 = par->sigma_x0*par->sigma_x0;
	int Nxmin_super, Nxmax_super, Nxmin_sub, Nxmax_sub, Nymin_super, Nymax_super, Nymin_sub, Nymax_sub;

	/* Propagativ modes limits */
	Nxmax_super =  FLOOR( ( sqrt(cabs(k_super2 - sigma_y02)) - sigma_x0 )/Delta_sigma_x );
	Nxmin_super = -FLOOR( ( sqrt(cabs(k_super2 - sigma_y02)) + sigma_x0 )/Delta_sigma_x );
	Nymax_super =  FLOOR( ( sqrt(cabs(k_super2 - sigma_x02)) - sigma_y0 )/Delta_sigma_y );
	Nymin_super = -FLOOR( ( sqrt(cabs(k_super2 - sigma_x02)) + sigma_y0 )/Delta_sigma_y );
	Nxmax_sub =  FLOOR( ( sqrt(cabs(k_sub2 - sigma_y02)) - sigma_x0 )/Delta_sigma_x );
	Nxmin_sub = -FLOOR( ( sqrt(cabs(k_sub2 - sigma_y02)) + sigma_x0 )/Delta_sigma_x );
	Nymax_sub =  FLOOR( ( sqrt(cabs(k_sub2 - sigma_x02)) - sigma_y0 )/Delta_sigma_y );
	Nymin_sub = -FLOOR( ( sqrt(cabs(k_sub2 - sigma_x02)) + sigma_y0 )/Delta_sigma_y );

	int Nxlimit = MAX(MAX(Nxmax_super,Nxmax_sub),MAX(-Nxmin_super,-Nxmin_sub));
	if (Nx < Nxlimit) {
		fprintf(stderr,"CAUTION, Nx not sufficient to cover all propagativ orders\n");
		fprintf(stderr,"Minimum value : Nx = %d \n",Nxlimit);
		if (Nx < Nxmax_super)  Nxmax_super =  Nx;
		if (Nx < Nxmax_sub)    Nxmax_sub   =  Nx;
		if (Nxmin_super < -Nx) Nxmin_super = -Nx;
		if (Nxmin_sub   < -Nx) Nxmin_sub   = -Nx;
	}
	int Nylimit = MAX(MAX(Nymax_super,Nymax_sub),MAX(-Nymin_super,-Nymin_sub));
	if (Ny < Nylimit) {
		fprintf(stderr,"CAUTION, Ny not sufficient to cover all propagativ orders\n");
		fprintf(stderr,"Minimum value : Ny = %d \n",Nylimit);
		if (Ny < Nymax_super)  Nymax_super =  Ny;
		if (Ny < Nymax_sub)    Nymax_sub   =  Ny;
		if (Nymin_super < -Ny) Nymin_super = -Ny;
		if (Nymin_sub   < -Ny) Nymin_sub   = -Ny;
	}
	
	eff->Nxmin_super = Nxmin_super;
	eff->Nxmax_super = Nxmax_super;
	eff->Nxmin_sub = Nxmin_sub;
	eff->Nxmax_sub = Nxmax_sub;
	eff->Nymin_super = Nymin_super;
	eff->Nymax_super = Nymax_super;
	eff->Nymin_sub = Nymin_sub;
	eff->Nymax_sub = Nymax_sub;

	return 0;	
}

	
/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md3D_efficiencies(COMPLEX *Ai, COMPLEX *A0, COMPLEX *Ah, struct Param_struct *par,  struct Efficacites_struct *eff)
 *
 *	\brief		Efficiencies calculation
 */
/*-------------------------------------------------------------------------------------*/
int md3D_efficiencies(COMPLEX *Ai, COMPLEX *Ar, COMPLEX *At, struct Param_struct *par,  struct Efficacites_struct *eff)
{
	int n;
	int vec_size = par->vec_size;
	
	COMPLEX k_super = par->k_super;
	COMPLEX k_sub = par->k_sub;
	COMPLEX k_super2 = k_super*k_super;
	COMPLEX k_sub2 = k_sub*k_sub;
	int nx, ny;
	
	COMPLEX *kz_super, *kz_sub, *sigma_x, *sigma_y, C_super, C_sub, sigx, sigy, sigx2, sigy2;
	kz_super = par->kz_super;
	kz_sub = par->kz_sub;
	sigma_x = par->sigma_x;
	sigma_y = par->sigma_y;

	double Pz_i, Pz_r, Pz_t, sumPzi;
	COMPLEX *Eyi, *Hpyi, *Eyr, *Hpyr, *Eyt, *Hpyt;
	COMPLEX *Exi, *Hpxi, *Exr, *Hpxr, *Ext, *Hpxt;
	Exi = par->Exi;
	Exr = par->Exr;
	Ext = par->Ext;
	Hpxi = par->Hpxi;
	Hpxr = par->Hpxr;
	Hpxt = par->Hpxt;
	
	int Nxmin_super = eff->Nxmin_super;
	int Nxmin_sub   = eff->Nxmin_sub;
	int Nymin_super = eff->Nymin_super;
	int Nymin_sub   = eff->Nymin_sub;
	int Nxmax_super = eff->Nxmax_super;
	int Nxmax_sub   = eff->Nxmax_sub;
	int Nymax_super = eff->Nymax_super;
	int Nymax_sub   = eff->Nymax_sub;

	/* y field components */
	Eyi  = par->Ai;
	Hpyi = par->Ai + vec_size;
	Eyr  = par->Ar;
	Hpyr = par->Ar + vec_size;
	Eyt  = par->At;
	Hpyt = par->At + vec_size;

	/* x field components calculations */
	for(n=0;n<=vec_size-1;n++){
		C_super = (1/(k_super2 - sigma_y[n]*sigma_y[n]));
		C_sub   = (1/(k_sub2   - sigma_y[n]*sigma_y[n]));
		Exi[n]  = C_super * (-kz_super[n]*Hpyi[n] - sigma_y[n]*sigma_x[n]*Eyi[n]);
		Exr[n]  = C_super * ( kz_super[n]*Hpyr[n] - sigma_y[n]*sigma_x[n]*Eyr[n]);
		Ext[n]  = C_sub   * (-kz_sub[n]*Hpyt[n]   - sigma_y[n]*sigma_x[n]*Eyt[n]);
		Hpxi[n] = C_super * ( k_super2*kz_super[n]*Eyi[n] - sigma_y[n]*sigma_x[n]*Hpyi[n]);
		Hpxr[n] = C_super * (-k_super2*kz_super[n]*Eyr[n] - sigma_y[n]*sigma_x[n]*Hpyr[n]);
		Hpxt[n] = C_sub   * ( k_sub2*kz_sub[n]*Eyt[n]     - sigma_y[n]*sigma_x[n]*Hpyt[n]);
	}
	/* Initialising to nan */
	for (ny=Nymin_super;ny<=Nymax_super;ny++){
		for (nx=Nxmin_super;nx<=Nxmax_super;nx++){

			eff->eff_r   [ny-Nymin_super][nx-Nxmin_super] = 1/0.0;
			eff->nx_eff_r[ny-Nymin_super][nx-Nxmin_super] = 1/0.0;
			eff->ny_eff_r[ny-Nymin_super][nx-Nxmin_super] = 1/0.0;
		}
	}
	for (ny=Nymin_sub;ny<=Nymax_sub;ny++){
		for (nx=Nxmin_sub;nx<=Nxmax_sub;nx++){
			eff->eff_t   [ny-Nymin_sub][nx-Nxmin_sub] = 1/0.0;
			eff->nx_eff_t[ny-Nymin_sub][nx-Nxmin_sub] = 1/0.0;
			eff->ny_eff_t[ny-Nymin_sub][nx-Nxmin_sub] = 1/0.0;
		}
	}

	/* Total incident energie */
	sumPzi = 0.0;
	for (ny=Nymin_super;ny<=Nymax_super;ny++){
		for (nx=Nxmin_super;nx<=Nxmax_super;nx++){
			sigx = nx*par->Delta_sigma_x + par->sigma_x0;
			sigy = ny*par->Delta_sigma_y + par->sigma_y0;
			sigx2 = sigx*sigx;
			sigy2 = sigy*sigy;
			if (cabs(sigx2+sigy2) < cabs(k_super*k_super)) {
				n = (ny+par->Ny)*(2*par->Nx+1)+(nx+par->Nx);
				Pz_i = fabs(creal(Exi[n]*CONJ(Hpyi[n]) - Eyi[n]*CONJ(Hpxi[n])));
				sumPzi += Pz_i;
			}
		}
	}
/*	for (n=0; n<=vec_size-1; n++) {
		n_x=nx[n];
		n_y=ny[n];
		if (cabs(sigma_x[n]*sigma_x[n]+sigma_y[n]*sigma_y[n]) < cabs(k_super*k_super)) {
			Pz_i = fabs(creal(Exi[n]*CONJ(Hpyi[n]) - Eyi[n]*CONJ(Hpxi[n])));
			sumPzi += Pz_i;
		}
	}
*/
	/* Reflexion : efficiencies and directions */
/*	eff->sum_eff_r = 0.0;
printf("k2=%f\n",cabs(k_super*k_super));
	for (n=0; n<=vec_size-1; n++) {
		n_x=nx[n];
		n_y=ny[n];
printf("sigx=%f, sigy=%f ",cabs(sigma_x[n]),cabs(sigma_y[n]));
		if (cabs(sigma_x[n]*sigma_x[n]+sigma_y[n]*sigma_y[n]) < cabs(k_super*k_super)) {
printf("nx=%d ny=%d sigx2=%f, sigy2=%f ",n_y,n_x,cabs(sigma_x[n]*sigma_x[n]),cabs(sigma_y[n]*sigma_y[n]));
			Pz_r = fabs(creal(Exr[n]*CONJ(Hpyr[n]) - Eyr[n]*CONJ(Hpxr[n])));
printf("eff_r[%d][%d]\n",n_y-Nymin_super,n_x-Nxmin_super);
			eff->eff_r[n_y-Nymin_super][n_x-Nxmin_super] = Pz_r/sumPzi;
			eff->nx_eff_r[n_y-Nymin_super][n_x-Nxmin_super] = (double)n_x ;
			eff->ny_eff_r[n_y-Nymin_super][n_x-Nxmin_super] = (double)n_y ;
			eff->sum_eff_r += eff->eff_r[n_y-Nymin_super][n_x-Nxmin_super];
		}
	}*/
	/* Reflexion : efficiencies and directions */
	eff->sum_eff_r = 0.0;
	for (ny=Nymin_super;ny<=Nymax_super;ny++){
		for (nx=Nxmin_super;nx<=Nxmax_super;nx++){
			sigx = nx*par->Delta_sigma_x + par->sigma_x0;
			sigy = ny*par->Delta_sigma_y + par->sigma_y0;
			sigx2 = sigx*sigx;
			sigy2 = sigy*sigy;
			if (cabs(sigx2+sigy2) < cabs(k_super*k_super)) {
				n = (ny+par->Ny)*(2*par->Nx+1)+(nx+par->Nx);
				Pz_r = fabs(creal(Exr[n]*CONJ(Hpyr[n]) - Eyr[n]*CONJ(Hpxr[n])));
				eff->eff_r   [ny-Nymin_super][nx-Nxmin_super] = Pz_r/sumPzi;
				eff->nx_eff_r[ny-Nymin_super][nx-Nxmin_super] = (double)nx ;
				eff->ny_eff_r[ny-Nymin_super][nx-Nxmin_super] = (double)ny ;
				eff->sum_eff_r += eff->eff_r[ny-Nymin_super][nx-Nxmin_super];
			}
		}
	}

	/* Transmission : efficiencies and directions */
	eff->sum_eff_t = 0.0;
	for (ny=Nymin_sub;ny<=Nymax_sub;ny++){
		for (nx=Nxmin_sub;nx<=Nxmax_sub;nx++){
			sigx = nx*par->Delta_sigma_x + par->sigma_x0;
			sigy = ny*par->Delta_sigma_y + par->sigma_y0;
			sigx2 = sigx*sigx;
			sigy2 = sigy*sigy;
			if (cabs(sigx2+sigy2) < cabs(k_sub*k_sub)) {
				n = (ny+par->Ny)*(2*par->Nx+1)+(nx+par->Nx);
				Pz_t = fabs(creal(Ext[n]*CONJ(Hpyt[n]) - Eyt[n]*CONJ(Hpxt[n])));
				eff->eff_t[ny-Nymin_sub][nx-Nxmin_sub] = Pz_t/sumPzi;
				eff->nx_eff_t[ny-Nymin_sub][nx-Nxmin_sub] = (double) nx;
				eff->ny_eff_t[ny-Nymin_sub][nx-Nxmin_sub] = (double) ny;
				eff->sum_eff_t += eff->eff_t[ny-Nymin_sub][nx-Nxmin_sub];
			}
		}
	}

	
/*printf("\nRe eff_R :\n");
SaveDbleTab2file (eff->eff_r, Nmax_super-Nmin_super+1, "stdout", " ");*/
	/* Transmission : efficiencies and directions */
/*	eff->sum_eff_t = 0.0;
	for (n=0; n<=vec_size-1; n++) {
		n_x=nx[n];
		n_y=ny[n];
		if (cabs(sigma_x[n]*sigma_x[n]+sigma_y[n]*sigma_y[n]) < cabs(k_sub*k_sub)) {
			Pz_t = fabs(creal(Ext[n]*CONJ(Hpyt[n]) - Eyt[n]*CONJ(Hpxt[n])));
			eff->eff_t[n_y-Nymin_sub][n_x-Nxmin_sub] = Pz_t/sumPzi;
			eff->nx_eff_t[n_y-Nymin_sub][n_x-Nxmin_sub] = (double) n_x;
			eff->ny_eff_t[n_y-Nymin_sub][n_x-Nxmin_sub] = (double) n_y;
			eff->sum_eff_t += eff->eff_t[n_y-Nymin_sub][n_x-Nxmin_sub];
		}
	}*/
/*printf("\nRe eff_T :\n");
SaveDbleTab2file (eff->eff_t, Nmax_sub-Nmin_sub+1, "stdout", " ");*/
	
	/* Sum of efficiencies */
	eff->sum_eff = eff->sum_eff_r + eff->sum_eff_t;

	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md3D_amplitudes(COMPLEX *Ai, COMPLEX *Ar, COMPLEX *At, COMPLEX **S12, 
 *                               COMPLEX **S22, struct Param_struct *par)
 *
 *	\brief	Field amplitudes calculations
 */
/*-------------------------------------------------------------------------------------*/
int md3D_amplitudes(COMPLEX *Ai, COMPLEX *Ar, COMPLEX *At, COMPLEX **S12, COMPLEX **S22, struct Param_struct *par)
{
	int n;
	int vec_size = par->vec_size;
	COMPLEX *kz_super;
	kz_super = par->kz_super;
		
	/* Vi = Ai*cexp(-I*kz_super*h) */
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
/*
printf("\n Ai : \n");SaveCplxTab2file (par->Ai, par->vec_size, "Re", "stdout", " ",999999999,"\n");
printf("\n Ar : \n");SaveCplxTab2file (par->Ar, par->vec_size, "Re", "stdout", " ",999999999,"\n");
printf("\n At : \n");SaveCplxTab2file (par->At, par->vec_size, "Re", "stdout", " ",999999999,"\n");
*/
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
	int i,j, nS;
	int vec_size = par->vec_size;
	int NS = par->NS;
	COMPLEX **S11, **S12, **S21, **S22, **T11, **T12, **T21, **T22, **Z, **T_tmp, **T_tmp2, **T_tmp3;
	S11 = par->S11;
	S12 = par->S12;
	S21 = par->S21;
	S22 = par->S22;
	T11 = par->T11;
	T12 = par->T12;
	T21 = par->T21;
	T22 = par->T22;

	Z      = allocate_CplxMatrix(2*vec_size,2*vec_size);
	T_tmp  = allocate_CplxMatrix(2*vec_size,2*vec_size);
	T_tmp2 = allocate_CplxMatrix(2*vec_size,2*vec_size);
	T_tmp3 = allocate_CplxMatrix(2*vec_size,2*vec_size);
	
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
		T_Matrix(T11, T12, T21, T22, nS, par);

		/* Z = inv(T11 + T12*S12) */
		invM(Z, add_M(T_tmp2, 
			T11, M_x_M(T_tmp,
				T12,S12, 2*vec_size, 2*vec_size), 2*vec_size, 2*vec_size), 2*vec_size);

		/* S12 = (T21 +T22*S12)*Z */
		M_x_M(S12,
			add_M(T_tmp2, T21, M_x_M(T_tmp,
					T22,S12, 2*vec_size, 2*vec_size), 2*vec_size, 2*vec_size),
			Z, 2*vec_size, 2*vec_size);
						
		/* S22 = S22*Z */
		M_equals(T_tmp,S22, 2*vec_size, 2*vec_size);
		M_x_M(S22,T_tmp,Z, 2*vec_size, 2*vec_size);
		
		/* Following blocks only usefull when some light is coming from below */
		/* NOT TESTED YET, but carrefully written... (should work !) */
		/* (Matrices should also be initialized ... !) */
		/* The order (S22, S12),S21,S11 must not be changed */
		/* S21 = S21 - S22*T12*S11 */
/*		sub_M(S21, S21, M_x_M(T_tmp, 
			S22, 
			M_x_M(T_tmp2, T12, S11, 2*vec_size, 2*vec_size), 2*vec_size, 2*vec_size), 2*vec_size, 2*vec_size);
*/
		/* S11 = (T22 - S12*T12)*S11 */
/*		M_equals(T_tmp, S11, 2*vec_size, 2*vec_size);
		M_x_M(S11,
			sub_M(T_tmp2, T22, M_x_M(T_tmp3,
					S12,T12, 2*vec_size, 2*vec_size), 2*vec_size, 2*vec_size),
			T_tmp, 2*vec_size, 2*vec_size);
*/
		/* if NEAR_FIELD, field components are saved in the Near_field_matrix */
/*		if (par->SAVE_NEAR_FIELD){
			md3D_save_near_field(nS, par);
		}
*/		
	}
#if 0
printf("\nT11");SaveMatrix2file (T11, 2*vec_size, 2*vec_size, "Re", "stdout");
printf("\nT12");SaveMatrix2file (T12, 2*vec_size, 2*vec_size, "Re", "stdout");
printf("\nT21");SaveMatrix2file (T21, 2*vec_size, 2*vec_size, "Re", "stdout");
printf("\nT22");SaveMatrix2file (T22, 2*vec_size, 2*vec_size, "Re", "stdout");
printf("\nS12");SaveMatrix2file (S12, 2*vec_size, 2*vec_size, "Re", "stdout");
printf("\nS22");SaveMatrix2file (S22, 2*vec_size, 2*vec_size, "Re", "stdout");
#endif	
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

/*-------------------------------------------------------------------------------------*/
/*!	\fn	int S_matrix_stack(struct Param_struct *par)
 *
 *		\brief	S matrix calculation in case of stack with some homogeneous layers
 */
/*-------------------------------------------------------------------------------------*/
#if 0
int S_matrix_stack(struct Param_struct *par)
{
	if (par->verbosity) fprintf(stdout,"S-Matrix calculation. STACK case\n");
	int i,j, nS, n_stack, N_stack, n_patterned_layer;
	int vec_size = par->vec_size;
	int NS = par->NS;
	COMPLEX **S11, **S12, **S21, **S22, **T11, **T12, **T21, **T22, **Z, **T_tmp, **T_tmp2, **T_tmp3, **Psi_sub_tmp, nu_subsuper;
	S11 = par->S11;
	S12 = par->S12;
	S21 = par->S21;
	S22 = par->S22;
	T11 = par->T11;
	T12 = par->T12;
	T21 = par->T21;
	T22 = par->T22;
	N_stack = par->N_stack;
	n_patterned_layer = par->n_patterned_layer;

	Z      = allocate_CplxMatrix(2*vec_size,2*vec_size);
	T_tmp  = allocate_CplxMatrix(2*vec_size,2*vec_size);
	T_tmp2 = allocate_CplxMatrix(2*vec_size,2*vec_size);
	T_tmp3 = allocate_CplxMatrix(2*vec_size,2*vec_size);
	
	/* Initialisations */
	for (i=0; i<=2*vec_size-1; i++) {
		for (j=0; j<=2*vec_size-1; j++) {
			S12[i][j] = 0;
			S22[i][j] = 0;
		}
		S22[i][i] = 1;
	}

	/* Iterations */
	for (n_stack=0; n_stack<=N_stack-1; n_stack++){

		if (n_stack == n_patterned_layer) { /* Determinate the S-matrix of the structured layer (plus layers situated below), using the differential method and the S-matrix algorithm */
			/* Redefinition of nu_sub to adapt to md3D_T_matrix */
			Psi_sub_tmp = par->Psi_sub;
			if (n_patterned_layer != 0){ 
				par->Psi_sub = par->Psi_super;
			}
			for (nS=0; nS<=NS-1; nS++) {
				/* T-Matrix calculation */
				T_Matrix(T11, T12, T21, T22, nS, par);
				/* Z = inv(T11 + T12*S12) */
				invM(Z, add_M(T_tmp2, 
					T11, M_x_M(T_tmp,
						T12,S12, 2*vec_size, 2*vec_size), 2*vec_size, 2*vec_size), 2*vec_size);
				/* S12 = (T21 +T22*S12)*Z */
				M_x_M(S12,
					add_M(T_tmp2, T21, M_x_M(T_tmp,
							T22,S12, 2*vec_size, 2*vec_size), 2*vec_size, 2*vec_size),
					Z, 2*vec_size, 2*vec_size);
				/* S22 = S22*Z */
				M_equals(T_tmp,S22, 2*vec_size, 2*vec_size);
				M_x_M(S22,T_tmp,Z, 2*vec_size, 2*vec_size);
			}
			/* Reattributing right value to nu_sub */
			par->Psi_sub = Psi_sub_tmp;
		}else{
			/* if first iteration, we consider the infinitly thin layer below the homogeneous layer to be the 
                           same material as the substrate, otherwise it is considered to be the same as the superstrate (same as in T_matrix) */
			if (n_stack == 0){ 
				nu_subsuper = par->nu_sub;
			}else{
				nu_subsuper = par->nu_super;
			}
			/* Calculation of the T-Matrix of the homogeneous layer */
			T_Matrix_homog_layer(T11, T12, T21, T22, par->nu_super, nu_subsuper, par->nu_stack[n_stack], par->h_stack[n_stack], par->lambda, par->sigma_x, par->sigma_y, vec_size, par);
#if 0
printf("\nT11");SaveMatrix2file (T11, 2*vec_size, 2*vec_size, "Im", "stdout");
printf("\nT12");SaveMatrix2file (T12, 2*vec_size, 2*vec_size, "Im", "stdout");
printf("\nT21");SaveMatrix2file (T21, 2*vec_size, 2*vec_size, "Im", "stdout");
printf("\nT22");SaveMatrix2file (T22, 2*vec_size, 2*vec_size, "Im", "stdout");
#endif	
			/* Z = inv(T11 + T12*S12) */
			invM(Z, add_M(T_tmp2, 
				T11, M_x_M(T_tmp,
					T12,S12, 2*vec_size, 2*vec_size), 2*vec_size, 2*vec_size), 2*vec_size);
			/* S12 = (T21 +T22*S12)*Z */
			M_x_M(S12,
				add_M(T_tmp2, T21, M_x_M(T_tmp,
						T22,S12, 2*vec_size, 2*vec_size), 2*vec_size, 2*vec_size),
				Z, 2*vec_size, 2*vec_size);
			/* S22 = S22*Z */
			M_equals(T_tmp,S22, 2*vec_size, 2*vec_size);
			M_x_M(S22,T_tmp,Z, 2*vec_size, 2*vec_size);
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
#endif

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int T_Matrix(COMPLEX **T11, COMPLEX **T12, COMPLEX **T21, COMPLEX **T22,
				int nS, struct Param_struct *par)
 *
 *	\brief Calcul de la matrice T
 */
/*-------------------------------------------------------------------------------------*/
int T_Matrix(COMPLEX **T11, COMPLEX **T12, COMPLEX **T21, COMPLEX **T22, int nS, struct Param_struct *par)
{
	int i,j,vec_size = par->vec_size;
	
	double hmin, hmax, z0, Delta_z;
	COMPLEX **Psi, **invPsi_super, **M_buffer_4vecsize, **T;
	M_buffer_4vecsize = allocate_CplxMatrix(4*vec_size,4*vec_size);
	T   = allocate_CplxMatrix(4*vec_size,4*vec_size);

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

	invPsi_super = par->invPsi_super;

	/* z of the considered T matrix slice */
	hmin = par->h*nS/par->NS;
	hmax = par->h*(nS+1)/par->NS;
	if (par->tab_NS_ENABLED){
		hmin = par->tab_NS[nS];
		hmax = par->tab_NS[nS+1];
	}
	Delta_z = hmax - hmin;
	z0 = hmin;

	/* P_matrix calculation */
	(*par->P_matrix)(par->P, z0, Delta_z, par);

	/* T matrix : T = inv(Psi_super) * P_matrix * Psi */
	/*************************************************************/
	/* SHOULD BE DONE MANUALLY SINCE Psi Matrix are quite simple */
	M_x_M(T,
			invPsi_super, M_x_M(M_buffer_4vecsize, par->P, Psi, 4*vec_size, 4*vec_size), 4*vec_size, 4*vec_size);
	/*************************************************************/

	/* putting it block T matrices */
	for(i=0;i<=2*vec_size-1;i++){
		for(j=0;j<=2*vec_size-1;j++){
			T11[i][j] = T[i][j];
			T12[i][j] = T[i][j+2*vec_size];
			T21[i][j] = T[i+2*vec_size][j];
			T22[i][j] = T[i+2*vec_size][j+2*vec_size];
		}
	}





#if 0
printf("\nRe(par->P) :\n");
SaveMatrix2file (par->P, 4*par->vec_size, 4*par->vec_size, "Re", "stdout");
printf("\nIm(par->P) :\n");
SaveMatrix2file (par->P, 4*par->vec_size, 4*par->vec_size, "Im", "stdout");
printf("\nRe(par->M) :\n");
SaveMatrix2file (par->M, 4*par->vec_size, 4*par->vec_size, "Re", "stdout");
printf("\nIm(par->M) :\n");
SaveMatrix2file (par->M, 4*par->vec_size, 4*par->vec_size, "Im", "stdout");
printf("\nRe(Psi) :\n");
SaveMatrix2file (Psi, 4*par->vec_size, 4*par->vec_size, "Re", "stdout");
printf("\nIm(Psi) :\n");
SaveMatrix2file (Psi, 4*par->vec_size, 4*par->vec_size, "Im", "stdout");
printf("\nRe(invPsi_super) :\n");
SaveMatrix2file (invPsi_super, 4*par->vec_size, 4*par->vec_size, "Re", "stdout");
printf("\nIm(invPsi_super) :\n");
SaveMatrix2file (invPsi_super, 4*par->vec_size, 4*par->vec_size, "Im", "stdout");
printf("\nRe(par->T) :\n");
SaveMatrix2file (par->T, 4*par->vec_size, 4*par->vec_size, "Re", "stdout");
printf("\nIm(par->T) :\n");
SaveMatrix2file (par->T, 4*par->vec_size, 4*par->vec_size, "Im", "stdout");
#endif

	/*Affichage du temps restant à l'écran */
	md3D_affichTemps(par->Nx,par->Ny,nS,par->NS,par);
	
	free(T[0]);
	free(T);
	free(M_buffer_4vecsize[0]);
	free(M_buffer_4vecsize);

	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		T_Matrix_homog_layer(COMPLEX **T11, COMPLEX **T12, COMPLEX **T21, COMPLEX **T22, COMPLEX nu_super, COMPLEX nu_sub, COMPLEX nu_layer, double h_layer, double lambda, COMPLEX *sigma_x, COMPLEX *sigma_y, int vec_size, struct Param_struct *par)
 *
 *	\brief T-matrix of a homogeneous layer
 */
/*-------------------------------------------------------------------------------------*/
#if 0
int T_Matrix_homog_layer(COMPLEX **T11, COMPLEX **T12, COMPLEX **T21, COMPLEX **T22, COMPLEX nu_super, COMPLEX nu_sub, COMPLEX nu_layer, double h_layer, double lambda, COMPLEX *sigma_x, COMPLEX *sigma_y, int vec_size, struct Param_struct *par)
{
	int i,j;
	
	double z0;
	COMPLEX **Psi, **invPsi_super, **M_buffer_4vecsize, **T;
	COMPLEX **k2_ptr_tmp;
	COMPLEX **invk2_ptr_tmp;
	M_buffer_4vecsize = allocate_CplxMatrix(4*vec_size,4*vec_size);
	T   = allocate_CplxMatrix(4*vec_size,4*vec_size);

	/* Psi matrix */
	if (nu_super==par->nu_sub){ /* 1st S-Matrix iteration : we are in the substrat */
		Psi = par->Psi_sub;
	}else if (nu_super==par->nu_super){     /* Following iterations : we are in the superstrat */
		Psi = par->Psi_super;
	}else{
		fprintf(stderr, "%s, line %d : ERROR, wrong nusubsuper value\n",__FILE__,__LINE__);
		exit(EXIT_FAILURE);
	}

	invPsi_super = par->invPsi_super;


	/* P_matrix calculation */
	par->nu_this_layer = nu_layer;
	/* special function pointers for k2 and 1/k2 determination */
	k2_ptr_tmp = par->k_2;
	invk2_ptr_tmp = par->invk_2;
	par->k_2 = &k2_homog;
	par->invk_2 = &invk2_homog;
	z0 = 0; fprintf(stderr,"/////////////////////////////\n// T_Matrix_homog_layer: z0 = 0. Check that other value also works\n/////////////////////////////\n");

	/* Calculate the P-matrix */
	(*zinvar_P_matrix)(par->P, z0, h_layer, par);

	/* Restitute original values to k2 and 1/k2 function pointers */
	par->k_2 = k2_ptr_tmp;
	par->invk_2 = invk2_ptr_tmp;

	/* T matrix : T = inv(Psi_super) * P_matrix * Psi */
	/*************************************************************/
	/* SHOULD BE DONE MANUALLY SINCE Psi Matrix are quite simple */
	M_x_M(T,
			invPsi_super, M_x_M(M_buffer_4vecsize, par->P, Psi, 4*vec_size, 4*vec_size), 4*vec_size, 4*vec_size);
	/*************************************************************/

	/* putting it block T matrices */
	for(i=0;i<=2*vec_size-1;i++){
		for(j=0;j<=2*vec_size-1;j++){
			T11[i][j] = T[i][j];
			T12[i][j] = T[i][j+2*vec_size];
			T21[i][j] = T[i+2*vec_size][j];
			T22[i][j] = T[i+2*vec_size][j+2*vec_size];
		}
	}





#if 0
printf("\nRe(par->P) :\n");
SaveMatrix2file (par->P, 4*par->vec_size, 4*par->vec_size, "Re", "stdout");
printf("\nIm(par->P) :\n");
SaveMatrix2file (par->P, 4*par->vec_size, 4*par->vec_size, "Im", "stdout");
printf("\nRe(par->M) :\n");
SaveMatrix2file (par->M, 4*par->vec_size, 4*par->vec_size, "Re", "stdout");
printf("\nIm(par->M) :\n");
SaveMatrix2file (par->M, 4*par->vec_size, 4*par->vec_size, "Im", "stdout");
printf("\nRe(Psi) :\n");
SaveMatrix2file (Psi, 4*par->vec_size, 4*par->vec_size, "Re", "stdout");
printf("\nIm(Psi) :\n");
SaveMatrix2file (Psi, 4*par->vec_size, 4*par->vec_size, "Im", "stdout");
printf("\nRe(invPsi_super) :\n");
SaveMatrix2file (invPsi_super, 4*par->vec_size, 4*par->vec_size, "Re", "stdout");
printf("\nIm(invPsi_super) :\n");
SaveMatrix2file (invPsi_super, 4*par->vec_size, 4*par->vec_size, "Im", "stdout");
printf("\nRe(par->T) :\n");
SaveMatrix2file (par->T, 4*par->vec_size, 4*par->vec_size, "Re", "stdout");
printf("\nIm(par->T) :\n");
SaveMatrix2file (par->T, 4*par->vec_size, 4*par->vec_size, "Im", "stdout");
#endif

	
	free(T[0]);
	free(T);
	free(M_buffer_4vecsize[0]);
	free(M_buffer_4vecsize);

	return 0;
}
#endif

#if 0 /* Has never worked: Not finished */
/*-------------------------------------------------------------------------------------*/
/*!	\fn		T_Matrix_homog_layer(COMPLEX **T11, COMPLEX **T12, COMPLEX **T21, COMPLEX **T22, COMPLEX nu_super, COMPLEX nu_sub, COMPLEX nu_layer, double h_layer, double lambda, COMPLEX *sigma_x, COMPLEX *sigma_y, int vec_size)
 *
 *	\brief T-matrix of a homogeneous layer
 */
/*-------------------------------------------------------------------------------------*/
int T_Matrix_homog_layer(COMPLEX **T11, COMPLEX **T12, COMPLEX **T21, COMPLEX **T22, COMPLEX nu_super, COMPLEX nu_sub, COMPLEX nu_layer, double h_layer, double lambda, COMPLEX *sigma_x, COMPLEX *sigma_y, int vec_size)
{
	int n,m;
	COMPLEX mu_r0 = 1, mu_r1 = 1, mu_r2 = 1; /* nonmagnetic materials */
	COMPLEX eps_r0, eps_r1, eps_r2, k, *gamma0, *gamma1, *gamma2, C1, C0, C;
	COMPLEX *T11s, *T12s, *T21s, *T22s, *T11p, *T12p, *T21p, *T22p, **T_buffer;
	T_buffer= allocate_CplxMatrix(11,vec_size);
	T11s=T_buffer[1]; T12s=T_buffer[2]; T21s=T_buffer[3]; T22s=T_buffer[4];
	T11p=T_buffer[5]; T12p=T_buffer[6]; T21p=T_buffer[7]; T22p=T_buffer[0];
	gamma0=T_buffer[8]; gamma1=T_buffer[9]; gamma2=T_buffer[10];
	double h = h_layer;
	eps_r0 = nu_sub*nu_sub;
	eps_r1 = nu_layer*nu_layer;
	eps_r2 = nu_super*nu_super;
	k=2*PI/lambda;

	for (n=0;n<=vec_size-1;n++){
		gamma0[n] = csqrt (nu_sub * nu_sub * k*k - sigma_x[n]*sigma_x[n] - sigma_y[n]*sigma_y[n]);
		gamma1[n] = csqrt (nu_layer*nu_layer*k*k - sigma_x[n]*sigma_x[n] - sigma_y[n]*sigma_y[n]);
		gamma2[n] = csqrt (nu_super*nu_super*k*k - sigma_x[n]*sigma_x[n] - sigma_y[n]*sigma_y[n]);
		g0=gamma0[n];
		g1=gamma1[n];
		g2=gamma2[n];
		/* TE */
		C1=(g1/g2)*(mu_r2/mu_r1);
		C0=(g0/g1)*(mu_r1/mu_r0);
		T11s[n]=0.25*((1-C1)*(1-C0)*cexp(I*g1*h) + (1+C1)*(1+C0)*cexp(-I*g1*h));
		T12s[n]=0.25*((1-C1)*(1+C0)*cexp(I*g1*h) + (1+C1)*(1-C0)*cexp(-I*g1*h));
		T21s[n]=0.25*((1+C1)*(1-C0)*cexp(I*g1*h) + (1-C1)*(1+C0)*cexp(-I*g1*h));
		T22s[n]=0.25*((1+C1)*(1+C0)*cexp(I*g1*h) + (1-C1)*(1-C0)*cexp(-I*g1*h));
		/* TM */
		C1=(g1/g2)*(eps_r2/eps_r1);
		C0=(g0/g1)*(eps_r1/eps_r0);
		T11p[n]=0.25*((1-C1)*(1-C0)*cexp(I*g1*h) + (1+C1)*(1+C0)*cexp(-I*g1*h));
		T12p[n]=0.25*((1-C1)*(1+C0)*cexp(I*g1*h) + (1+C1)*(1-C0)*cexp(-I*g1*h));
		T21p[n]=0.25*((1+C1)*(1-C0)*cexp(I*g1*h) + (1-C1)*(1+C0)*cexp(-I*g1*h));
		T22p[n]=0.25*((1+C1)*(1+C0)*cexp(I*g1*h) + (1-C1)*(1-C0)*cexp(-I*g1*h));
	}
	/* T11 = [T11s 0   ],  T12 = [T12s 0   ],  T21 = [T21s 0   ],  T22 = [T22s 0   ] 
	         [0    T11p]         [0    T12p]         [0    T21p]         [0    T22p] */
	for (n=0;n<=2*vec_size-1;n++){
		for (m=0;m<=2*vec_size-1;m++){
			T11[n][m]=0;
		}
	}
	for (n=0;n<=vec_size-1;n++){
		g0=gamma0[n];
		k0=k_sub;
		/* iR0 = inv(R0) */
		C = 1/(1+(g0*g0/(k0*k0))*sin(Phi)*sin(Phi));
		iR011 = C*cos(Phi);
		iR012 = C*(g0/k0*k0)*sin(Phi);
		iR021 = C*g0*sin(Phi);
		iR022 = -C*cos(Phi);		
		/* R2 */

		/* T11 = R2*T11*iR0 */
	C11 = T11s[n]*iR011
	C12 = T11s[n]*iR012
	C21 = T11p[n]*iR021
	C22 = T11p[n]*iR022

		T11ss = R211*C11 + R212*C21
		T


		/* T11 */
		T11[n]         [n]         =T11ss;
		T11[n]         [n+vec_size]=T11sp;
		T11[n+vec_size][n]         =T11ps;
		T11[n+vec_size][n+vec_size]=T11pp;
		/* T12 */
		T12[n]         [n]         =T12s[n];
		T12[n]         [n+vec_size]=0;
		T12[n+vec_size][n]         =0;
		T12[n+vec_size][n+vec_size]=T12p[n];
		/* T21 */
		T21[n]         [n]         =T21s[n];
		T21[n]         [n+vec_size]=0;
		T21[n+vec_size][n]         =0;
		T21[n+vec_size][n+vec_size]=T21p[n];
		/* T22 */
		T22[n]         [n]         =T22s[n];
		T22[n]         [n+vec_size]=0;
		T22[n+vec_size][n]         =0;
		T22[n+vec_size][n+vec_size]=T22p[n];
	}

	free(T_buffer[0]);
	free(T_buffer);

	return 0;
}
#endif

/*-------------------------------------------------------------------------------------*/
/*!	int zinvar_P_matrix(COMPLEX **P, double z, double Delta_z, struct Param_struct *par)
 *
 *	\brief P_matrix calculation in the case of z invariance
 */
/*-------------------------------------------------------------------------------------*/
int zinvar_P_matrix(COMPLEX **P, double z0, double Delta_z, struct Param_struct *par)
{
	int i,j, P_size = 4*par->vec_size;
	double z = z0 + 0.5*Delta_z;
	
	COMPLEX *eig_values, *eig_buffer, *M_sol, **EigVectors, **invEigVec, **invVec_Psi, **M_invVec_Psi, **M_buffer_4vecsize;

	eig_values = (COMPLEX *) malloc(sizeof(COMPLEX)*P_size);
	eig_buffer = (COMPLEX *) malloc(sizeof(COMPLEX)*50*P_size);
	M_sol =      (COMPLEX *) malloc(sizeof(COMPLEX)*P_size);
	EigVectors        = allocate_CplxMatrix(P_size,P_size);
	invEigVec         = allocate_CplxMatrix(P_size,P_size);
	invVec_Psi        = allocate_CplxMatrix(P_size,P_size);
	M_invVec_Psi      = allocate_CplxMatrix(P_size,P_size);
	M_buffer_4vecsize = allocate_CplxMatrix(P_size,P_size);

	/* M matrix calculation */
	(*par->M_matrix)(par->M, z, par);

	/* Diagonalisation of M */
	eigen_values(par->M, eig_values, EigVectors, eig_buffer, P_size);

	/* Solution of the diagonalized system */
	for (i=0;i<=P_size-1;i++){
		M_sol[i] = cexp(eig_values[i]*Delta_z);
	}
		
	/* invEigVec = inv(EigVectors) */
	invM(invEigVec, EigVectors, P_size);

	/* M_invVec_Psi = M_sol * invVec_Psi */
	for (i=0;i<=P_size-1;i++){
		for (j=0;j<=P_size-1;j++){
			M_invVec_Psi[i][j] = M_sol[i] * invEigVec[i][j];
		}
	}

	/* P = EigVec * M_diag * inv(EigVec) */
	M_x_M(P, EigVectors, M_invVec_Psi, P_size, P_size);

	par->N_steps++;

	free(eig_values);
	free(eig_buffer);
	free(M_sol);
	free(EigVectors[0]);
	free(invEigVec[0]);
	free(invVec_Psi[0]);
	free(M_invVec_Psi[0]);
	free(M_buffer_4vecsize[0]);
	free(EigVectors);
	free(invEigVec);
	free(invVec_Psi);
	free(M_invVec_Psi);
	free(M_buffer_4vecsize);
		
	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	int rk4_P_matrix(COMPLEX **P, double z0, double Delta_z, struct Param_struct *par)
 *
 *	\brief P_matrix calculation with 4th order Runge Kutta vectorized method
 */
/*-------------------------------------------------------------------------------------*/
int rk4_P_matrix(COMPLEX **P, double z0, double dz, struct Param_struct *par)
{
	int i,j, P_size = 4*par->vec_size;
	COMPLEX **Mz, **Mzd, **Mzdd, **M1, **M2, **M3, **M4, **M_tmp1, **M_tmp2, **M_tmp3;
	double dz_2 = 0.5*dz;
	double dz_3 = dz/3.0;
	double dz_6 = dz/6.0;

	M_tmp1 = allocate_CplxMatrix(P_size,P_size);
	M_tmp2 = allocate_CplxMatrix(P_size,P_size);
	M_tmp3 = allocate_CplxMatrix(P_size,P_size);


	/* M matrix calculation */


	/* M1 = Mz */
	Mz = M_tmp1;
	(*par->M_matrix)(Mz,   z0,         par);
	M1 = Mz;
	/* P = Id + (h/6)M1 */
	for (i=0;i<=P_size-1;i++){
		for (j=0;j<=P_size-1;j++){
			P[i][j] = dz_6*M1[i][j];
		}
		P[i][i] += 1; 
	}

	/* M2 = Mzd*(Id + M1*dz/2) */
	for (i=0;i<=P_size-1;i++){
		for (j=0;j<=P_size-1;j++){
			M_tmp2[i][j] = dz_2 * M1[i][j]; 
		}
		M_tmp2[i][i] += 1; 
	}
	Mzd = M_tmp1;
	(*par->M_matrix)(Mzd,  z0+dz_2	, par);
	M2  = M_tmp3;
	M_x_M(M2, Mzd, M_tmp2, P_size, P_size);
	/* P = Id + (h/6)M1 + (h/3)M2  */
	for (i=0;i<=P_size-1;i++){
		for (j=0;j<=P_size-1;j++){
			P[i][j] += dz_3*M2[i][j];
		}
	}

	/* M3 = Mzd*(Id + M2*dz/2) */
	for (i=0;i<=P_size-1;i++){
		for (j=0;j<=P_size-1;j++){
			M_tmp2[i][j] = dz_2 * M2[i][j]; 
		}
		M_tmp2[i][i] += 1; 
	}
	M3  = M_tmp3;
	M_x_M(M3, Mzd, M_tmp2, P_size, P_size);
	/* P = Id + (h/6)M1 + (h/3)M2 + (h/3)M3 */
	for (i=0;i<=P_size-1;i++){
		for (j=0;j<=P_size-1;j++){
			P[i][j] += dz_3*M3[i][j];
		}
	}
	
	/* M4 = Mzdd*(Id + M3*dz) */
	for (i=0;i<=P_size-1;i++){
		for (j=0;j<=P_size-1;j++){
			M_tmp2[i][j] = dz * M3[i][j]; 
		}
		M_tmp2[i][i] += 1; 
	}
	M4 = M_tmp3;
	Mzdd = M_tmp1;
	(*par->M_matrix)(Mzdd, z0+dz,      par);
	M_x_M(M4, Mzdd, M_tmp2, P_size, P_size);
	/* P = Id + (h/6)M1 + (h/3)M2 + (h/3)M3 + (h/6)M4 */
	for (i=0;i<=P_size-1;i++){
		for (j=0;j<=P_size-1;j++){
			P[i][j] += dz_6*M4[i][j];
		}
	}
	par->N_steps++;

	free(M_tmp1[0]);
	free(M_tmp2[0]);
	free(M_tmp3[0]);
	free(M_tmp1);
	free(M_tmp2);
	free(M_tmp3);

	return 0;
}

#if 0
/* Version of rk4_P_matrix with greater memory allocations */
int rk4_P_matrix(COMPLEX **P, double z0, double dz, struct Param_struct *par)
{
	int i,j, P_size = 4*par->vec_size;
	COMPLEX **Mz, **Mzd, **Mzdd, **M2, **M3, **M4,**M_tmp;

	Mz    = allocate_CplxMatrix(P_size,P_size);
	Mzd   = allocate_CplxMatrix(P_size,P_size);
	Mzdd  = allocate_CplxMatrix(P_size,P_size);
	M2    = allocate_CplxMatrix(P_size,P_size);
	M3    = allocate_CplxMatrix(P_size,P_size);
	M4    = allocate_CplxMatrix(P_size,P_size);
	M_tmp = allocate_CplxMatrix(P_size,P_size);

	double dz_2 = 0.5*dz;

	/* M matrix calculation */
	(*par->M_matrix)(Mz,   z0,         par);
	(*par->M_matrix)(Mzd,  z0+dz_2	, par);
	(*par->M_matrix)(Mzdd, z0+dz,      par);

	/* M1 = Mz */

	/* M2 = Mzd*(Id + M1*dz/2) */
	for (i=0;i<=P_size-1;i++){
		for (j=0;j<=P_size-1;j++){
			M_tmp[i][j] = dz_2 * Mz[i][j]; 
		}
		M_tmp[i][i] += 1; 
	}
	M_x_M(M2, Mzd, M_tmp, P_size, P_size);

	/* M3 = Mzd*(Id + M2*dz/2) */
	for (i=0;i<=P_size-1;i++){
		for (j=0;j<=P_size-1;j++){
			M_tmp[i][j] = dz_2 * M2[i][j]; 
		}
		M_tmp[i][i] += 1; 
	}
	M_x_M(M3, Mzd, M_tmp, P_size, P_size);

	
	/* M4 = Mzdd*(Id + M3*dz) */
	for (i=0;i<=P_size-1;i++){
		for (j=0;j<=P_size-1;j++){
			M_tmp[i][j] = dz * M3[i][j]; 
		}
		M_tmp[i][i] += 1; 
	}
	M_x_M(M4, Mzdd, M_tmp, P_size, P_size);

	/* P = Id + (h/6)M1 + (h/3)M2 + (h/3)M3 + (h/6)M4 */
	double dz_6 = dz/6.0;
	double dz_3 = dz/3.0;
	for (i=0;i<=P_size-1;i++){
		for (j=0;j<=P_size-1;j++){
			P[i][j] = dz_6*(Mz[i][j] + M4[i][j]) + dz_3*(M2[i][j] + M3[i][j]);
		}
		P[i][i] += 1; 
	}
#if 0
printf("\nIm(M2) :\n");
SaveMatrix2file (M2, P_size, P_size, "Im", "stdout");
printf("\nIm(M3) :\n");
SaveMatrix2file (M3, P_size, P_size, "Im", "stdout");
printf("\nIm(M4) :\n");
SaveMatrix2file (M4, P_size, P_size, "Im", "stdout");
printf("\nIm(P) :\n");
SaveMatrix2file (P, P_size, P_size, "Im", "stdout");
#endif
	par->N_steps++;

	free(Mz[0]);
	free(Mzd[0]);
	free(Mzdd[0]);
	free(M2[0]);
	free(M3[0]);
	free(M4[0]);
	free(M_tmp[0]);
	free(Mz);
	free(Mzd);
	free(Mzdd);
	free(M2);
	free(M3);
	free(M4);
	free(M_tmp);
		
	return 0;
}
#endif

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int M_matrix(COMPLEX **M, double z, struct Param_struct *par)
 *
 *		\brief	M matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int M_matrix(COMPLEX **M, double z, struct Param_struct *par)
{
	int i,j;
	int vec_size = par->vec_size;

	COMPLEX *sigma_x = par->sigma_x;
	COMPLEX *sigma_y = par->sigma_y;
	COMPLEX **M11,**M12,**M13,**M14, **M21,**M22,**M23,**M24, **M31,**M32,**mM33,**M34, **M41,**M42,**M43,**mM44;
	COMPLEX **invQzzQzx, **invQzzQzy, **invQzzSigmax, **invQzzSigmay, **Qtmp;
	COMPLEX **Qxx, **Qxy, **Qxz, **Qyx, **Qyy, **Qyz, **Qzx, **Qzy, **Qzz, **invQzz;
	Qxx = par->Qxx;
	Qxy = par->Qxy;
	Qxz = par->Qxz;
	Qyy = par->Qyy;
	Qyx = Qxy;
	Qyz = par->Qyz;
	Qzx = Qxz;
	Qzy = Qyz;
	Qzz = par->Qzz;
	
	invQzz = par->Qzz_1;
	
	M11 = allocate_CplxMatrix(vec_size,vec_size);
	M12 = allocate_CplxMatrix(vec_size,vec_size);
	M13 = allocate_CplxMatrix(vec_size,vec_size);
	M14 = allocate_CplxMatrix(vec_size,vec_size);
	M21 = allocate_CplxMatrix(vec_size,vec_size);
	M22 = allocate_CplxMatrix(vec_size,vec_size);
	M23 = allocate_CplxMatrix(vec_size,vec_size);
	M24 = allocate_CplxMatrix(vec_size,vec_size);
	M31 = allocate_CplxMatrix(vec_size,vec_size);
	M32 = allocate_CplxMatrix(vec_size,vec_size);
	mM33 = allocate_CplxMatrix(vec_size,vec_size);
	M34 = allocate_CplxMatrix(vec_size,vec_size);
	M41 = allocate_CplxMatrix(vec_size,vec_size);
	M42 = allocate_CplxMatrix(vec_size,vec_size);
	M43 = allocate_CplxMatrix(vec_size,vec_size);
	mM44 = allocate_CplxMatrix(vec_size,vec_size);
	invQzzQzx = allocate_CplxMatrix(vec_size,vec_size);
	invQzzQzy = allocate_CplxMatrix(vec_size,vec_size);
	invQzzSigmax = allocate_CplxMatrix(vec_size,vec_size);
	invQzzSigmay = allocate_CplxMatrix(vec_size,vec_size);
	Qtmp = allocate_CplxMatrix(vec_size,vec_size);

	/* Toeplitz matrices calculations */
	md3D_QMatrix(z, Qxx, Qxy, Qxz, Qyy, Qyz, Qzz, invQzz, par);
	/* invQzzQzx = invQzz*Qzx */
	M_x_M(invQzzQzx, invQzz, Qzx, vec_size, vec_size);
	/* invQzzQzy = invQzz*Qzy */
	M_x_M(invQzzQzy, invQzz, Qzy, vec_size, vec_size);
	/* invQzzSigmax = invQzz*sigma_x */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			invQzzSigmax[i][j] = invQzz[i][j]*sigma_x[j];
		}
	}
	/* invQzzSigmay = invQzz*sigma_y */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			invQzzSigmay[i][j] = invQzz[i][j]*sigma_y[j];
		}
	}

	/* M11 = -sigma_x*invQzzQzx */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M11[i][j] = -sigma_x[i]*invQzzQzx[i][j];
		}
	}
	/* M12 = -sigma_x*invQzzQzy */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M12[i][j] = -sigma_x[i]*invQzzQzy[i][j];
		}
	}
	/* M13 = sigma_x*invQzzSigmay */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M13[i][j] = sigma_x[i]*invQzzSigmay[i][j];
		}
	}
	/* M14 = Id - sigma_x*invQzzSigmax */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M14[i][j] = -sigma_x[i]*invQzzSigmax[i][j];
		}
		M14[i][i] += 1;
	}

	/* M21 = -sigma_y*invQzzQzx */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M21[i][j] = -sigma_y[i]*invQzzQzx[i][j];
		}
	}
	/* M22 = -sigma_y*invQzzQzy */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M22[i][j] = -sigma_y[i]*invQzzQzy[i][j];
		}
	}
	/* M23 = -Id + sigma_y*invQzzSigmay */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M23[i][j] = sigma_y[i]*invQzzSigmay[i][j];
		}
		M23[i][i] += -1;
	}
	/* M24 = -sigma_y*invQzzSigmax */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M24[i][j] = -sigma_y[i]*invQzzSigmax[i][j];
		}
	}

	/* M31 = -sigma_x*sigma_y - Qyx + Qyz*invQzzQzx */
	M_x_M(Qtmp,Qyz,invQzzQzx, vec_size, vec_size);
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M31[i][j] = -Qyx[i][j] + Qtmp[i][j];
		}
		M31[i][i] -= sigma_x[i]*sigma_y[i];
	}
	/* M32 = sigma_x^2 - Qyy + Qyz*invQzzQzy */
	M_x_M(Qtmp,Qyz,invQzzQzy, vec_size, vec_size);
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M32[i][j] = -Qyy[i][j] + Qtmp[i][j];
		}
		M32[i][i] += sigma_x[i]*sigma_x[i];
	}
	/* M33 = -Qyz*invQzzSigmay , mM33 = -M33 */
	M_x_M(mM33,Qyz,invQzzSigmay, vec_size, vec_size);
	/* M34 = Qyz*invQzzSigmax */
	M_x_M(M34,Qyz,invQzzSigmax, vec_size, vec_size);

	/* M41 = -sigma_y^2 + Qxx - Qxz*invQzzQzx */
	M_x_M(Qtmp,Qxz,invQzzQzx, vec_size, vec_size);
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M41[i][j] = Qxx[i][j] - Qtmp[i][j];
		}
		M41[i][i] -= sigma_y[i]*sigma_y[i];
	}
	/* M42 = sigma_y*sigma_x + Qxy - Qxz*invQzzQzy */
	M_x_M(Qtmp,Qxz,invQzzQzy, vec_size, vec_size);
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M42[i][j] = Qxy[i][j] - Qtmp[i][j];
		}
		M42[i][i] += sigma_y[i]*sigma_x[i];
	}
	/* M43 = Qxz*invQzzSigmay */
	M_x_M(M43,Qxz,invQzzSigmay, vec_size, vec_size);
	/* M44 = -Qxz*invQzzSigmax , mM44 = -M44 */
	M_x_M(mM44,Qxz,invQzzSigmax, vec_size, vec_size);

	/* putting it all together */
	for(i=0;i<=vec_size-1;i++){
		for(j=0;j<=vec_size-1;j++){
			M[i][j]            = I * M11[i][j];
			M[i][j+vec_size]   = I * M12[i][j];
			M[i][j+2*vec_size] = I * M13[i][j];
			M[i][j+3*vec_size] = I * M14[i][j];

			M[i+vec_size][j]            = I * M21[i][j];
			M[i+vec_size][j+vec_size]   = I * M22[i][j];
			M[i+vec_size][j+2*vec_size] = I * M23[i][j];
			M[i+vec_size][j+3*vec_size] = I * M24[i][j];

			M[i+2*vec_size][j]            =  I * M31[i][j];
			M[i+2*vec_size][j+vec_size]   =  I * M32[i][j];
			M[i+2*vec_size][j+2*vec_size] = -I * mM33[i][j];
			M[i+2*vec_size][j+3*vec_size] =  I * M34[i][j];

			M[i+3*vec_size][j]            =  I * M41[i][j];
			M[i+3*vec_size][j+vec_size]   =  I * M42[i][j];
			M[i+3*vec_size][j+2*vec_size] =  I * M43[i][j];
			M[i+3*vec_size][j+3*vec_size] = -I * mM44[i][j];
		}
	}

#if 0
printf("\nz=%f\n",z);
printf("\nRe(M) :\n");SaveMatrix2file (M, 4*vec_size, 4*vec_size, "Re", "stdout");
printf("\nIm(M) :\n");SaveMatrix2file (M, 4*vec_size, 4*vec_size, "Im", "stdout");
printf("\nRe(M11) :\n");SaveMatrix2file (M11, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M11) :\n");SaveMatrix2file (M11, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M12) :\n");SaveMatrix2file (M12, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M12) :\n");SaveMatrix2file (M12, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M13) :\n");SaveMatrix2file (M13, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M13) :\n");SaveMatrix2file (M13, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M14) :\n");SaveMatrix2file (M14, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M14) :\n");SaveMatrix2file (M14, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M21) :\n");SaveMatrix2file (M21, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M21) :\n");SaveMatrix2file (M21, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M22) :\n");SaveMatrix2file (M22, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M22) :\n");SaveMatrix2file (M22, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M23) :\n");SaveMatrix2file (M23, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M23) :\n");SaveMatrix2file (M23, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M24) :\n");SaveMatrix2file (M24, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M24) :\n");SaveMatrix2file (M24, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M31) :\n");SaveMatrix2file (M31, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M31) :\n");SaveMatrix2file (M31, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M32) :\n");SaveMatrix2file (M32, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M32) :\n");SaveMatrix2file (M32, vec_size, vec_size, "Im", "stdout");
printf("\nRe(mM33) :\n");SaveMatrix2file (mM33, vec_size, vec_size, "Re", "stdout");
printf("\nIm(mM33) :\n");SaveMatrix2file (mM33, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M34) :\n");SaveMatrix2file (M34, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M34) :\n");SaveMatrix2file (M34, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M41) :\n");SaveMatrix2file (M41, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M41) :\n");SaveMatrix2file (M41, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M42) :\n");SaveMatrix2file (M42, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M42) :\n");SaveMatrix2file (M42, vec_size, vec_size, "Im", "stdout");
printf("\nRe(M43) :\n");SaveMatrix2file (M43, vec_size, vec_size, "Re", "stdout");
printf("\nIm(M43) :\n");SaveMatrix2file (M43, vec_size, vec_size, "Im", "stdout");
printf("\nRe(mM44) :\n");SaveMatrix2file (mM44, vec_size, vec_size, "Re", "stdout");
printf("\nIm(mM44) :\n");SaveMatrix2file (mM44, vec_size, vec_size, "Im", "stdout");
#endif

	free(M11[0]);
	free(M12[0]);
	free(M13[0]);
	free(M14[0]);
	free(M21[0]);
	free(M22[0]);
	free(M23[0]);
	free(M24[0]);
	free(M31[0]);
	free(M32[0]);
	free(mM33[0]);
	free(M34[0]);
	free(M41[0]);
	free(M42[0]);
	free(M43[0]);
	free(mM44[0]);
	free(invQzzQzx[0]);
	free(invQzzQzy[0]);
	free(invQzzSigmax[0]);
	free(invQzzSigmay[0]);
	free(Qtmp[0]);
	free(M11);
	free(M12);
	free(M13);
	free(M14);
	free(M21);
	free(M22);
	free(M23);
	free(M24);
	free(M31);
	free(M32);
	free(mM33);
	free(M34);
	free(M41);
	free(M42);
	free(M43);
	free(mM44);
	free(invQzzQzx);
	free(invQzzQzy);
	free(invQzzSigmax);
	free(invQzzSigmay);
	free(Qtmp);
	
	return 0;
}
/*-------------------------------------------------------------------------------------*/
/*!	\fn		int zinvar_M_Matrix(COMPLEX **M, double z, struct Param_struct *par)
 *
 *		\brief	M matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int zinvar_M_matrix(COMPLEX **M, double z, struct Param_struct *par)
{
	printf("\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#\n#* zinvar_M_matrix : Not implemented...\n#* Calling normal M_matrix instead\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#");
	return M_matrix(M, z, par);
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md3D_QMatrix(double z, COMPLEX **Qxx, COMPLEX **Qxy, COMPLEX **Qxz, COMPLEX **Qyy, COMPLEX **Qyz, COMPLEX **Qzz, COMPLEX **Qzz_1, struct Param_struct *par)
 *
 *		\brief	Q matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
int md3D_QMatrix(double z, COMPLEX **Qxx, COMPLEX **Qxy, COMPLEX **Qxz, COMPLEX **Qyy, COMPLEX **Qyz, COMPLEX **Qzz, COMPLEX **Qzz_1, struct Param_struct *par)
{
	int vec_size=par->vec_size, Nprx=par->Nprx, Npry=par->Npry, Nx=par->Nx, Ny=par->Ny;
	COMPLEX **M_buff1, **toepk2, **invtoepk2, **moinsDelta, *buffer_NprxNpry;

	M_buff1 = allocate_CplxMatrix(vec_size,vec_size);
	toepk2 = allocate_CplxMatrix(vec_size,vec_size);
	invtoepk2 = allocate_CplxMatrix(vec_size,vec_size);
	moinsDelta = allocate_CplxMatrix(vec_size,vec_size);

	md3D_toepNorm(z, par);

	buffer_NprxNpry = (COMPLEX *) malloc(sizeof(COMPLEX)*Nprx*Npry);

	/* toepk2 = [[k2]] */
	toeplitz_2D(toepk2, Nx, Ny, (*par->k_2)(par, buffer_NprxNpry, z), Nprx, Npry);

	/* invtoepk2 = inv([[1/k2]]) */
	toeplitz_2D(M_buff1, Nx, Ny, (*par->invk_2)(par, buffer_NprxNpry, z), Nprx, Npry);
	invM(invtoepk2, M_buff1, vec_size);

	/* moinsDelta = inv([[1/k2]]) - [[k2]] */
	sub_M(moinsDelta, invtoepk2, toepk2, vec_size, vec_size);

	/* Qxx = toepk2 -DeltaNxx */
	add_M(Qxx, toepk2, M_x_M(M_buff1, moinsDelta, par->Nxx, vec_size, vec_size), vec_size, vec_size);

	/* Qxy = -DeltaNxy */
	M_x_M(Qxy, moinsDelta, par->Nxy, vec_size, vec_size);

	/* Qxz = -DeltaNxz */
	M_x_M(Qxz,moinsDelta,par->Nxz,vec_size,vec_size);

	/* Qyy = toepk2 -DeltaNyy */
	add_M(Qyy, toepk2, M_x_M(M_buff1,moinsDelta,par->Nyy,vec_size,vec_size), vec_size, vec_size);

	/* Qyz = -DeltaNyz */
	M_x_M(Qyz,moinsDelta,par->Nyz,vec_size,vec_size);

	/* Qzz = toepk2 -DeltaNzz */
	add_M(Qzz, toepk2, M_x_M(M_buff1,moinsDelta,par->Nzz,vec_size,vec_size), vec_size, vec_size);

	/* Qzz_1 = inv(Qzz) */
	invM(Qzz_1, Qzz, vec_size);
	
#if 0
printf("\n Rek2\n");SaveCplxTab2file ((*par->k_2)(par, buffer_NprxNpry, z), Npry*Nprx, "Re", "stdout", " ", Nprx, "\n");
SaveMatrix2file (toepk2, vec_size, vec_size, "Re", "stdout");
printf("\n toepk2\n");SaveMatrix2file (toepk2, vec_size, vec_size, "Re", "stdout");
printf("\n invtoepk2\n");SaveMatrix2file (invtoepk2, vec_size, vec_size, "Re", "stdout");
printf("\n moinsDelta\n");SaveMatrix2file (moinsDelta, vec_size, vec_size, "Re", "stdout");
printf("\n Qxx\n");SaveMatrix2file (Qxx, vec_size, vec_size, "Re", "stdout");
printf("\n Qxy\n");SaveMatrix2file (Qxy, vec_size, vec_size, "Re", "stdout");
printf("\n Qxz\n");SaveMatrix2file (Qyz, vec_size, vec_size, "Re", "stdout");
printf("\n Qyy\n");SaveMatrix2file (Qyy, vec_size, vec_size, "Re", "stdout");
printf("\n Qyz\n");SaveMatrix2file (Qyz, vec_size, vec_size, "Re", "stdout");
printf("\n Qzz\n");SaveMatrix2file (Qzz, vec_size, vec_size, "Re", "stdout");
#endif
	free(buffer_NprxNpry);
	free(M_buff1[0]);
	free(M_buff1);
	free(toepk2[0]);
	free(toepk2);
	free(invtoepk2[0]);
	free(invtoepk2);
	free(moinsDelta[0]);
	free(moinsDelta);
	
	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md3D_toepNorm(double z, struct Param_struct *par)
 *
 *		\brief	Toeplitz matrix normal vector calculation
 */
/*-------------------------------------------------------------------------------------*/
int md3D_toepNorm(double z, struct Param_struct *par)
{
	int i, j, Nprx=par->Nprx, Npry=par->Npry, Nx=par->Nx, Ny=par->Ny;
	COMPLEX **normx, **normy, **normz, **normxx, **normxy, **normxz, **normyy, **normyz, **normzz;

	if (par->toepNorm_CALCULATED == 1){
		return 0;
	}else{		
		normx = allocate_CplxMatrix(Npry,Nprx);
		normy = allocate_CplxMatrix(Npry,Nprx);
		normz = allocate_CplxMatrix(Npry,Nprx);
		normxx = allocate_CplxMatrix(Npry,Nprx);
		normxy = allocate_CplxMatrix(Npry,Nprx);
		normxz = allocate_CplxMatrix(Npry,Nprx);
		normyy = allocate_CplxMatrix(Npry,Nprx);
		normyz = allocate_CplxMatrix(Npry,Nprx);
		normzz = allocate_CplxMatrix(Npry,Nprx);
	
		(*par->Normal_function)(normx, normy, normz, z, par->profil[0], par);	

		for (i=0;i<=Npry-1;i++){
			for (j=0;j<=Nprx-1;j++){
				normxx[i][j] = normx[i][j]*normx[i][j];
				normxy[i][j] = normx[i][j]*normy[i][j];
				normxz[i][j] = normx[i][j]*normz[i][j];
				normyy[i][j] = normy[i][j]*normy[i][j];
				normyz[i][j] = normy[i][j]*normz[i][j];
				normzz[i][j] = normz[i][j]*normz[i][j];
			}
		}

		/* Nxx, Nxy, Nxz, Nyy, Nyz & Nzz */
		toeplitz_2D(par->Nxx, Nx, Ny, normxx[0], Nprx, Npry);
		toeplitz_2D(par->Nxy, Nx, Ny, normxy[0], Nprx, Npry);
		toeplitz_2D(par->Nxz, Nx, Ny, normxz[0], Nprx, Npry);
		toeplitz_2D(par->Nyy, Nx, Ny, normyy[0], Nprx, Npry);
		toeplitz_2D(par->Nyz, Nx, Ny, normyz[0], Nprx, Npry);
		toeplitz_2D(par->Nzz, Nx, Ny, normzz[0], Nprx, Npry);

		free(normx[0]);
		free(normx);
		free(normy[0]);
		free(normy);
		free(normz[0]);
		free(normz);
		free(normxx[0]);
		free(normxx);
		free(normxy[0]);
		free(normxy);
		free(normxz[0]);
		free(normxz);
		free(normyy[0]);
		free(normyy);
		free(normyz[0]);
		free(normyz);
		free(normzz[0]);
		free(normzz);
		if (par->profile_type == H_XY){ /* N have to be calculated only once for this type of profiles */
			par->toepNorm_CALCULATED = 1;
		}else{
			par->toepNorm_CALCULATED = 0;
		}
		
		return 0;
	}
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		COMPLEX **PsiMatrix(COMPLEX **Psi, COMPLEX k, COMPLEX *kz, COMPLEX *sigma_x, COMPLEX *sigma_y, int vec_size)
 *
 *		\brief	Psi matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX **PsiMatrix(COMPLEX **Psi, COMPLEX k, COMPLEX *kz, COMPLEX *sigma_x, COMPLEX *sigma_y, int vec_size)
{
	int i,j;
	COMPLEX p, qe, qh;
	COMPLEX k2 = k*k;
	COMPLEX sigma_y2;

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
		sigma_y2 = sigma_y[i]*sigma_y[i];
   	p  = -sigma_y[i]*sigma_x[i]/(k2-sigma_y2);
   	qe = -kz[i]/(k2-sigma_y2);
   	qh = k2*kz[i]/(k2-sigma_y2);
	
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
	return Psi;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn	COMPLEX **toeplitz_2D(COMPLEX **toep, int Nx, int Ny, COMPLEX *M_in, int Nxin, int Nyin)
 *
 *	\brief	Toeplitz matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX **toeplitz_2D(COMPLEX **toep, int Nx, int Ny, COMPLEX *M_in, int Nxin, int Nyin)
{
	int i,j;
	COMPLEX **tmp, **TF_2D;
	fftw_plan plan_TF2D;
	
	tmp = allocate_CplxMatrix(Nyin,Nxin);

	plan_TF2D = fftw_plan_dft_2d(Nxin, Nyin, M_in, tmp[0], FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_TF2D); 

	int Ntfx = 2*Nx;
	int Ntfy = 2*Ny;

	TF_2D = allocate_CplxMatrix(2*Ntfy+1,2*Ntfx+1);

	double coefnorm = 1.0/(Nxin*Nyin);
	for (i=0;i<=Ntfy-1;i++){
		for (j=0;j<=Ntfx-1;j++){
			TF_2D[i]     [j]      = tmp[Nyin-Ntfy+i][Nxin-Ntfx+j] * coefnorm;
			TF_2D[i+Ntfy][j]      = tmp[i]          [Nxin-Ntfx+j] * coefnorm;
			TF_2D[i]     [j+Ntfx] = tmp[Nyin-Ntfy+i][j] * coefnorm;
			TF_2D[i+Ntfy][j+Ntfx] = tmp[i]          [j] * coefnorm;
		}
		TF_2D[i]     [2*Ntfx] = tmp[Nyin-Ntfy+i][Ntfx] * coefnorm;
		TF_2D[i+Ntfy][2*Ntfx] = tmp[i]          [Ntfx] * coefnorm;
	}
	for (j=0;j<=Ntfx-1;j++){
		TF_2D[2*Ntfy][j]      = tmp[Ntfy][Nxin-Ntfx+j] * coefnorm;
		TF_2D[2*Ntfy][j+Ntfx] = tmp[Ntfy][j] * coefnorm;
	}
	TF_2D[2*Ntfy][2*Ntfx] = tmp[Ntfy][Ntfx] * coefnorm;

	fftw_destroy_plan(plan_TF2D);
	free(tmp[0]);
	free(tmp);

	int v_size = (2*Nx+1)*(2*Ny+1);
	int *nx, *ny;
	nx = (int *) malloc(sizeof(int)*v_size);
	ny = (int *) malloc(sizeof(int)*v_size);

	for (i=-Ny;i<=Ny;i++){
		for (j=-Nx;j<=Nx;j++){
			nx[(i+Ny)*(2*Nx+1)+(j+Nx)] = j;
			ny[(i+Ny)*(2*Nx+1)+(j+Nx)] = i;
		}
	}

/*	int Ncentre = 4*Nx*Ny+3*Nx+2*Ny+1;*/
	for (i=0;i<=v_size-1;i++){
		for (j=0;j<=v_size-1;j++){
			toep[i][j] = TF_2D[ny[i]-ny[j]+Ntfy][nx[i]-nx[j]+Ntfx];
		}
	}
	free(nx);
	free(ny);
	free(TF_2D[0]);
	free(TF_2D);

	return toep;
}
#if 0 /* Error in this function. Replaced by toeplitz_2D */
/*-------------------------------------------------------------------------------------*/
/*!	\fn	COMPLEX **toeplitz_2D(COMPLEX **toep, int Nx, int Ny, COMPLEX *M_in, int Nxin, int Nyin)
 *
 *	\brief	Toeplitz matrix calculation
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX **old_toeplitz_2D(COMPLEX **toep, int Nx, int Ny, COMPLEX *M_in, int Nxin, int Nyin)
{
	int i,j;
	COMPLEX **tmp, **TF_2D, *coef2;
	fftw_plan plan_TF2D;
	
	tmp = allocate_CplxMatrix(Nyin,Nxin);
	TF_2D = allocate_CplxMatrix(4*Ny+3,2*Nx+1);

	plan_TF2D = fftw_plan_dft_2d(Nxin, Nyin, M_in, tmp[0], FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_TF2D); 

	int Ntfx = Nx;
	int Ntfy = 2*Ny+1;

	double coefnorm = 1.0/(Nxin*Nyin);
	for (i=0;i<=Ntfy-1;i++){
		for (j=0;j<=Ntfx-1;j++){
			TF_2D[i]     [j]      = tmp[Nyin-Ntfy+i][Nxin-Ntfx+j] * coefnorm;
			TF_2D[i+Ntfy][j]      = tmp[i]          [Nxin-Ntfx+j] * coefnorm;
			TF_2D[i]     [j+Ntfx] = tmp[Nyin-Ntfy+i][j] * coefnorm;
			TF_2D[i+Ntfy][j+Ntfx] = tmp[i]          [j] * coefnorm;
		}
		TF_2D[i]     [2*Ntfx] = tmp[Nyin-Ntfy+i][Ntfx] * coefnorm;
		TF_2D[i+Ntfy][2*Ntfx] = tmp[i]          [Ntfx] * coefnorm;
	}
	for (j=0;j<=Ntfx-1;j++){
		TF_2D[2*Ntfy][j]      = tmp[Ntfy][Nxin-Ntfx+j] * coefnorm;
		TF_2D[2*Ntfy][j+Ntfx] = tmp[Ntfy][j] * coefnorm;
	}
	TF_2D[2*Ntfy][2*Ntfx] = tmp[Ntfy][Ntfx] * coefnorm;

	fftw_destroy_plan(plan_TF2D);
	free(tmp[0]);
	free(tmp);

	coef2 = TF_2D[0];

	int Ncentre = 4*Nx*Ny+3*Nx+2*Ny+1;
	int v_size = (2*Nx+1)*(2*Ny+1);
	for (i=0;i<=v_size-1;i++){
		for (j=0;j<=v_size-1;j++){
        toep[i][j] = coef2[Ncentre+i-j];
		}
	}
	free(TF_2D[0]);
	free(TF_2D);

	return toep;
}
#endif

/*-------------------------------------------------------------------------------------*/
/*!	\fn		COMPLEX *k2_H_XY(struct Param_struct *par, COMPLEX *k2_2D, double z)
 *
 *	\brief	Détermine le tableau de COMPLEXes k^2(x) pour un z donné
 *
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX *k2_H_XY(struct Param_struct *par, COMPLEX *k2_2D, double z)
{
	int i;
	
	for (i=0;i<=par->Nprx*par->Npry-1;i++){
		if (z > par->profil[0][i]) 
			k2_2D[i] = par->k2_layer[0]; /* Superstrat */
		else
			k2_2D[i] = par->k2_layer[1]; /* Substrat */
	}

	
	return k2_2D;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		COMPLEX *k2_homog(struct Param_struct *par, COMPLEX *k2_2D, double z)
 *
 *	\brief	Determine the COMPLEX array k^2(x) for a homogeneous layer (for compatibility with non homogeneous profiles)
 *
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX *k2_homog(struct Param_struct *par, COMPLEX *k2_2D, double z)
{
	int i;
	/* the index of the layer is placed  in par->nu_this_layer */
	COMPLEX	k = 2*PI*par->nu_this_layer/par->lambda;
	for (i=0;i<=par->Nprx*par->Npry-1;i++){
		k2_2D[i] = k*k;
	}
printf("\nCOUCOU, I am in k2_homog !\n\n");
	return k2_2D;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		COMPLEX *invk2_homog(struct Param_struct *par, COMPLEX *k2_2D, double z)
 *
 *	\brief	Determine the COMPLEX array 1/k^2(x) for a homogeneous layer (for compatibility with non homogeneous profiles)
 *
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX *invk2_homog(struct Param_struct *par, COMPLEX *invk2_2D, double z)
{
	int i;
	/* the index of the layer is placed  in par->nu_this_layer */
	COMPLEX	k = 2*PI*par->nu_this_layer/par->lambda;
	for (i=0;i<=par->Nprx*par->Npry-1;i++){
		invk2_2D[i] = 1/(k*k);
	}
printf("\nCOUCOU, I am in invk2_homog !\n\n");

	return invk2_2D;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn			COMPLEX *k2_MULTI(struct Param_struct *par, COMPLEX *k2_2D, double z)
 *
 *	\brief	Détermine le tableau de COMPLEXes k^2(xy) pour un z donné, pour un multicouches
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX *k2_MULTI(struct Param_struct *par, COMPLEX *k2_2D, double z)
{
	
	int nxy, n_layer=0;
	
	for (nxy=0; nxy<=par->Nprx*par->Npry-1; nxy++){
		do{
			if (z <= par->profil[n_layer][nxy]){
				if (z >= par->profil[n_layer+1][nxy]){
					k2_2D[nxy] = par->k2_layer[n_layer];
					break;
				}else{
					n_layer++;
				}
			}else{
				n_layer--;
			}	
		}while (1);
	}
	
	return k2_2D;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		COMPLEX *k2_N_XYZ(struct Param_struct *par, COMPLEX *k2_1D, double z)	
 *
 *		\brief	Détermine le tableau de COMPLEXes k^2(x) pour un z donné
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX *k2_N_XYZ(struct Param_struct *par, COMPLEX *k2_1D, double z)
{
	int i, nz;
	double DeuxPisurLambda2 = (2*PI/par->lambda)*(2*PI/par->lambda);
	double z_inv = par->h - z;
	
	nz = (int)((z_inv*par->Nprz)/par->h);
	if (nz > par->Nprz-1) nz = par->Nprz-1;
	if (nz < 0)          nz = 0;
	
	for (i=0;i<=par->Nprx*par->Npry-1;i++){
		k2_1D[i] = par->n_xyz[nz][i]*par->n_xyz[nz][i]*DeuxPisurLambda2; 
	}
/*
printf("\nreal(k2)\n");SaveCplxTab2file (k2_1D, par->Nprx*par->Npry, "Re", "stdout", " ",par->Nprx,"\n");
*/	
	return k2_1D;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn	COMPLEX *invk2_N_XYZ(struct Param_struct *par, COMPLEX *invk2_1D, double z)	
 *
 *	\brief	Détermine le tableau de COMPLEXes invk^2(x) pour un z donné
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX *invk2_N_XYZ(struct Param_struct *par, COMPLEX *invk2_1D, double z)
{
	int i, nz;
	double invDeuxPisurLambda2 = 1/(2*PI/par->lambda*2*PI/par->lambda);
	double z_inv = par->h - z;
	
	nz = (int)((z_inv*par->Nprz)/par->h);
	if (nz > par->Nprz-1) nz = par->Nprz-1;
	if (nz < 0)          nz = 0;
	
	for (i=0;i<=par->Nprx*par->Npry-1;i++){
		invk2_1D[i] = invDeuxPisurLambda2/(par->n_xyz[nz][i]*par->n_xyz[nz][i]); 
	}
/*
printf("\nreal(invk2)\n");SaveCplxTab2file (invk2_1D, par->Nprx*par->Npry, "Re", "stdout", " ",par->Nprx,"\n");
*/	
	return invk2_1D;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		COMPLEX *invk2_H_XY(struct Param_struct *par, COMPLEX *invk2_1D, double z)
 *
 *	\brief	Détermine le tableau de COMPLEXes 1/k^2(x) pour un z donné
 *
 *	\todo	Prend pour l'instant en compte seulement un profil de type h(x)\n
 *			Doit être plus polyvalent : accepter aussi les profils de type n(x,z)
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX *invk2_H_XY(struct Param_struct *par, COMPLEX *invk2_1D, double z)
{
	
	int i;
	
	for (i=0;i<=par->Nprx*par->Npry-1;i++){
		if (z > par->profil[0][i]) 
			invk2_1D[i] = par->invk2_layer[0]; /* Superstrat */
		else
			invk2_1D[i] = par->invk2_layer[1]; /* Substrat */
	}
	
	return invk2_1D;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn	COMPLEX *invk2_MULTI(struct Param_struct *par, COMPLEX *invk2_1D, double z)
 *
 *	\brief	Détermine le tableau de COMPLEXes 1/k^2(x) pour un z donné, pour un multicouches
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX *invk2_MULTI(struct Param_struct *par, COMPLEX *invk2_1D, double z)
{
	
	int nx, n_layer=0;
	
	for (nx=0; nx<=par->Nprx*par->Npry-1; nx++){
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
	
	return invk2_1D;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int Normal_H_XY(COMPLEX **normx, COMPLEX **normy, COMPLEX **normz, double z, struct Param_struct *par)
 *
 *		\brief	Determine normx, normy, normz
 *
 * 	\todo 	The PRECISION can be IMPROVED with a HIGHER ORDER calculation
 */
/*-------------------------------------------------------------------------------------*/
int Normal_H_XY(COMPLEX **norm_x, COMPLEX **norm_y, COMPLEX **norm_z, double z, double *profil, struct Param_struct *par)
{
	int i,j;
	int Nprx = par->Nprx;
	int Npry = par->Npry;
	double dhdx, dhdy;
	double two_dx = 2*par->Lx/Nprx;
	double two_dy = 2*par->Ly/Npry;

/*	double *profil = par->profil[0];*/

	if (Npry >= 2){ /* 2D Profile (general case) */
		for (i=1;i<=Npry-2;i++){
			for (j=1;j<=Nprx-2;j++){
				dhdx = (profil[Nprx*i+j+1]-profil[Nprx*i+j-1])/two_dx;
				dhdy = (profil[Nprx*(i+1)+j]-profil[Nprx*(i-1)+j])/two_dy;
				norm_x[i][j] = -dhdx/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
				norm_y[i][j] = -dhdy/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
				norm_z[i][j] = 1.0/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
			}
		}
		for (i=1;i<=Npry-2;i++){
			dhdx = (profil[Nprx*i+1]-profil[Nprx*i+Nprx-1])/two_dx;
			dhdy = (profil[Nprx*(i+1)]-profil[Nprx*(i-1)])/two_dy;
			norm_x[i][0] = -dhdx/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
			norm_y[i][0] = -dhdy/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
			norm_z[i][0] = 1.0/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
			dhdx = (profil[Nprx*i]-profil[Nprx*i+Nprx-2])/two_dx;
			dhdy = (profil[Nprx*(i+1)+Nprx-1]-profil[Nprx*(i-1)+Nprx-1])/two_dy;
			norm_x[i][Nprx-1] = -dhdx/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
			norm_y[i][Nprx-1] = -dhdy/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
			norm_z[i][Nprx-1] = 1.0/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
		}
		for (j=1;j<=Nprx-2;j++){
			dhdx = (profil[j+1]-profil[j-1])/two_dx;
			dhdy = (profil[Nprx+j]-profil[Nprx*(Npry-1)+j])/two_dy;
			norm_x[0][j] = -dhdx/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
			norm_y[0][j] = -dhdy/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
			norm_z[0][j] = 1.0/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
			dhdx = (profil[Nprx*(Npry-1)+j+1]-profil[Nprx*(Npry-1)+j-1])/two_dx;
			dhdy = (profil[j]-profil[Nprx*(Npry-2)+j])/two_dy;
			norm_x[Npry-1][j] = -dhdx/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
			norm_y[Npry-1][j] = -dhdy/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
			norm_z[Npry-1][j] = 1.0/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
		}
		dhdx = (profil[1]-profil[Nprx-1])/two_dx;
		dhdy = (profil[Nprx]-profil[Nprx*(Npry-1)])/two_dy;
		norm_x[0][0] = -dhdx/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
		norm_y[0][0] = -dhdy/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
		norm_z[0][0] = 1.0/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
		dhdx = (profil[(Npry-1)*Nprx]-profil[(Npry-1)*Nprx+Nprx-2])/two_dx;
		dhdy = (profil[Nprx-1]-profil[(Npry-2)*Nprx+Nprx-1])/two_dy;
		norm_x[Npry-1][Nprx-1] = -dhdx/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
		norm_y[Npry-1][Nprx-1] = -dhdy/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
		norm_z[Npry-1][Nprx-1] = 1.0/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
		dhdx = (profil[0]-profil[Nprx-2])/two_dx;
		dhdy = (profil[2*Nprx-1]-profil[(Npry-1)*Nprx+Nprx-1])/two_dy;
		norm_x[0][Nprx-1] = -dhdx/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
		norm_y[0][Nprx-1] = -dhdy/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
		norm_z[0][Nprx-1] = 1.0/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
		dhdx = (profil[(Npry-1)*Nprx+1]-profil[(Npry-1)*Nprx+Nprx-1])/two_dx;
		dhdy = (profil[0]-profil[(Npry-2)*Nprx])/two_dy;
		norm_x[Npry-1][0] = -dhdx/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
		norm_y[Npry-1][0] = -dhdy/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
		norm_z[Npry-1][0] = 1.0/(csqrt(1+dhdx*dhdx+dhdy*dhdy));
	}else if (Npry == 1){ /* Particular case of treatment of 1D profile (useful if you don't have the specific code but slower) */
		for (j=1;j<=Nprx-2;j++){
			dhdx = (profil[j+1]-profil[j-1])/two_dx;
			norm_x[0][j] = -dhdx/(csqrt(1+dhdx*dhdx));
			norm_y[0][j] = 0.0;
			norm_z[0][j] = 1.0/(csqrt(1+dhdx*dhdx));
		}
		dhdx = (profil[1]-profil[Nprx-1])/two_dx;
		norm_x[0][0] = -dhdx/(csqrt(1+dhdx*dhdx));
		norm_y[0][0] = 0.0;
		norm_z[0][0] = 1.0/(csqrt(1+dhdx*dhdx));
		dhdx = (profil[0]-profil[Nprx-2])/two_dx;
		norm_x[0][Nprx-1] = -dhdx/(csqrt(1+dhdx*dhdx));
		norm_y[0][Nprx-1] = 0.0;
		norm_z[0][Nprx-1] = 1.0/(csqrt(1+dhdx*dhdx));
	}else{ /* We should not be here... */
		fprintf(stderr, "%s, line %d : ERROR, Impossible value for Npry.\n",__FILE__,__LINE__);
		exit(EXIT_FAILURE);
	}

#if 0
/********************************************************/
if(par->verbosity == 11){ /*norm_x = 1 DEBUGGING and validation purpose only */
printf("\nCAUTION, norm_x set to 1, DEBUGGING and VALIDATION purpose only !\n\n");
for (i=0;i<=Npry-1;i++){
for (j=0;j<=Nprx-1;j++){
norm_x[i][j] = 1.0;
norm_y[i][j] = 0.0;
norm_z[i][j] = 0.0;
}}}
if(par->verbosity == 12){ /*norm_y = 1 DEBUGGING and validation purpose only */
printf("\nCAUTION, norm_y set to 1, DEBUGGING and VALIDATION purpose only !\n\n");
for (i=0;i<=Npry-1;i++){
for (j=0;j<=Nprx-1;j++){
norm_x[i][j] = 0.0;
norm_y[i][j] = 1.0;
norm_z[i][j] = 0.0;
}}}
if(par->verbosity == 13){ /*norm_z = 1 DEBUGGING and validation purpose only */
printf("\nCAUTION, norm_z set to 1, DEBUGGING and VALIDATION purpose only !\n\n");
for (i=0;i<=Npry-1;i++){
for (j=0;j<=Nprx-1;j++){
norm_x[i][j] = 0.0;
norm_y[i][j] = 0.0;
norm_z[i][j] = 1.0;
}}}
/********************************************************/
#endif


	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int Normal_H_XY(COMPLEX **normx, COMPLEX **normy, COMPLEX **normz, double z, struct Param_struct *par)
 *
 *	\brief	Determine normx, normy, normz. Normal vector is obtained as a linear interpolation beetwen normal vectors of the 2 adjacent interfaces.
 */
/*-------------------------------------------------------------------------------------*/
int Normal_MULTI(COMPLEX **norm_x, COMPLEX **norm_y, COMPLEX **norm_z, double z, double *profil, struct Param_struct *par)
{
	int nx,ny,nxy, n_layer=0;
	double zb,zt,dz,ct,cb, check_norm;

	COMPLEX ***tab_normx, ***tab_normy, ***tab_normz;
	tab_normx = allocate_CplxMatrix_3(par->Npry,par->Nprx, par->N_layers+3);
	tab_normy = allocate_CplxMatrix_3(par->Npry,par->Nprx, par->N_layers+3);
	tab_normz = allocate_CplxMatrix_3(par->Npry,par->Nprx, par->N_layers+3);
	
	/* Calculation of normal vectors of each interface */
	for (n_layer=0; n_layer<=par->N_layers+2; n_layer++){
		Normal_H_XY(tab_normx[n_layer], tab_normy[n_layer], tab_normz[n_layer], 0, par->profil[n_layer], par);
	}

	/* Linear interpolation (since most points are situated beetwen 2 interfaces and a continuous variation is desired) */
	n_layer=0;
	for (ny=0; ny<=par->Npry-1; ny++){
		for (nx=0; nx<=par->Nprx-1; nx++){
			nxy=par->Nprx*ny+nx;
			do{
				if (z <= par->profil[n_layer][nxy]){
					if (z >= par->profil[n_layer+1][nxy]){
						zt=par->profil[n_layer][nxy];
						zb=par->profil[n_layer+1][nxy];
						dz=zt-zb;
						if (dz <= par->h*1e-10){
							ct=0.5;
							cb=0.5;
						}else{
							ct=(zt-z)/dz;
							cb=(z-zb)/dz;
						}
						norm_x[ny][nx]=ct*tab_normx[n_layer][ny][nx]+cb*tab_normx[n_layer+1][ny][nx];
						norm_y[ny][nx]=ct*tab_normy[n_layer][ny][nx]+cb*tab_normy[n_layer+1][ny][nx];
						norm_z[ny][nx]=ct*tab_normz[n_layer][ny][nx]+cb*tab_normz[n_layer+1][ny][nx];
						check_norm=sqrt(cabs(norm_x[ny][nx]*norm_x[ny][nx]+norm_y[ny][nx]*norm_y[ny][nx]+norm_z[ny][nx]*norm_z[ny][nx]));
						norm_x[ny][nx]=norm_x[ny][nx]/check_norm;
						norm_y[ny][nx]=norm_y[ny][nx]/check_norm;
						norm_z[ny][nx]=norm_z[ny][nx]/check_norm;
						break;
					}else{
						n_layer++;
					}
				}else{
					n_layer--;
				}	
			}while (1);
		}
	}

	free(tab_normx[0][0]);
	free(tab_normx[0]);
	free(tab_normx);
	free(tab_normy[0][0]);
	free(tab_normy[0]);
	free(tab_normy);
	free(tab_normz[0][0]);
	free(tab_normz[0]);
	free(tab_normz);

	return 0;
}


/*!-------------------------------------------------------------------------------------
 * \fn int Normal_N_XY_ZINVAR(COMPLEX **normx, COMPLEX **normy, COMPLEX **normz, double z, struct Param_struct *par)
 *
 * \brief Determine normx, normy, normz for N_XY profile
 *
 *-------------------------------------------------------------------------------------*/
int Normal_N_XY_ZINVAR(COMPLEX **norm_x, COMPLEX **norm_y, COMPLEX **norm_z, double z, double *profil, struct Param_struct *par)
{
	int i,j;
	int Nprx = par->Nprx;
	int Npry = par->Npry;

	double x,y;
printf("\nCAUTION, norm set to RADIAL, DEBUGGING and VALIDATION purpose only !\n\n");
	for (i=0;i<=Npry-1;i++){
		for (j=0;j<=Nprx-1;j++){
			x=((double)j-((double)Nprx-1)/2);
			y=((double)i-((double)Npry-1)/2);
			if (cabs(x*x+y*y) == 0.0) x=0.0000001;
			norm_x[i][j] = x/csqrt(x*x+y*y);
			norm_y[i][j] = y/csqrt(x*x+y*y);
			norm_z[i][j] = 0.0;
		}
	}
return 1;


/*
printf("\nCAUTION, norm set norm_x = 1, DEBUGGING and VALIDATION purpose only !\n\n");
	for (i=0;i<=Npry-1;i++){
		for (j=0;j<=Nprx-1;j++){
			norm_x[i][j] = 1.0;
			norm_y[i][j] = 0.0;
			norm_z[i][j] = 0.0;
		}
	}
*/
	double *rnorm_x, *rnorm_y;
	rnorm_x = (double *) malloc(sizeof(double)*Nprx*Npry); 
	rnorm_y = (double *) malloc(sizeof(double)*Nprx*Npry); 
	if (par->verbosity >= 2) fprintf(stdout,"Reading norm_x and norm_y in %s (norm_z = 0)\n",par->profile_file);
	if (lire_tab(par->profile_file, "norm_x", rnorm_x, Nprx*Npry) == 0 && lire_tab(par->profile_file, "norm_y", rnorm_y, Nprx*Npry) == 0) {
		if (par->verbosity >= 2) fprintf(stdout,"OK\n");
	}else{
		fprintf(stderr,"norm reading ERROR\n");
		exit(EXIT_FAILURE);
	}
	for (i=0;i<=Npry-1;i++){
		for (j=0;j<=Nprx-1;j++){
			norm_x[i][j] = rnorm_x[i*Nprx + j];
			norm_y[i][j] = rnorm_y[i*Nprx + j];
			norm_z[i][j] = 0.0;
		}
	}
	free(rnorm_x);
	free(rnorm_y);


	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md3D_affichTemps(int n, int N, int nS, int NS, struct Param_struct *par)
 *
 *	\brief	Affichage du temps restant estimé en cours de calculs
 */
/*-------------------------------------------------------------------------------------*/
int md3D_affichTemps(int n, int N, int nS, int NS, struct Param_struct *par)
{

	/* Si moins de 5 secondes depuis le dernier affichage, on ne change rien */
	if (CHRONO(clock(), par->last_clock) < 5){
		return 0;
	/* Sinon, estimation et affichage de la durée restante */
	}else if (par->verbosity >= 1){
		int i;
		time(&par->last_time);
		par->last_clock = clock();
		float t_ecoule = difftime(par->last_time,par->time0);
		float t_total; 

		t_total = t_ecoule*NS/(nS+1);
		
		float t_restant = t_total - t_ecoule;
		int pourcent = ROUND(100.0*t_ecoule/t_total);
				
		fprintf(stdout,"\r");
		fprintf(stdout,"%3d %% [%ds ", pourcent,ROUND(t_ecoule));
		for (i=0;i<pourcent/5;i++) {fprintf(stdout,">");}
		for (i=pourcent/5;i<20;i++) {fprintf(stdout," ");}
		fprintf(stdout," %ds] %ds     ",ROUND(t_restant), ROUND(t_total));
		fflush(stdout);
	}
	
	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!   \fn     int md3D_make_tab_S_steps(struct Param_struct* par)
 *
 *    \brief  N steps tab making, when some S steps are imposed
 */
/*-------------------------------------------------------------------------------------*/
int md3D_make_tab_S_steps(struct Param_struct* par){

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

