/*!	\file		md1D.c
 *
 * 	\brief		Méthode différentielle 1D \n
 *				cas TE et TM \n
 *				algorithme Matrice-S
 *	\version	0.1
 *
 *	\date		../../2006
 *	\authors	Laurent ARNAUD
 */


#include "md1D.h"



/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_efficacites(complex *Ai, complex *A0, complex *Ah, 
                                 struct Param_struct *par, struct Efficacites_struct *eff)
 *
 *	\brief	Calcul des efficacités
 */
/*-------------------------------------------------------------------------------------*/
int md1D_efficacites(complex *Ai, complex *A0, complex *Ah, struct Param_struct *par,  struct Efficacites_struct *eff)
{
	int n;
	double coef_eff_T = 0;
/*	double K = 2.0*PI/par->L; */
	int Nsigma = par->N_sigma;
	int Nsigma_min = par->N_sigma_min;

	complex alpha_nh, alpha_n0;
	complex k0 = par->k0;
	complex kh = par->kh;
	complex kh2 = kh*kh;
	complex k02 = k0*k0;
	double *sigma=par->sigma;
	double delta_sigma = par->delta_sigma;
	double sigma0 = par->sigma0;
	int Nmin_super, Nmax_super, Nmin_sub, Nmax_sub,Nlimit_max,Nlimit_min;

	/* Calcul des limites des modes propagatifs */
	Nmax_super =  FLOOR((creal(k0) - sigma0)/delta_sigma) - Nsigma_min; /* Indices des limites des sigmas */
	Nmin_super = -FLOOR((creal(k0) + sigma0)/delta_sigma) - Nsigma_min; /* propagatifs dans le tableau de */
	Nmax_sub =  FLOOR((creal(kh) - sigma0)/delta_sigma) - Nsigma_min;    /* sigmas, pour le superstrat et le substrat */
	Nmin_sub = -FLOOR((creal(kh) + sigma0)/delta_sigma) - Nsigma_min;
	Nlimit_max = MAX(Nmax_super,Nmax_sub);
	Nlimit_min = MIN(Nmin_super,Nmin_sub);
	
	if (Nlimit_min < 0) {
		fprintf(stderr,"ATTENTION, N_sigma_min pas assez petit pour représenter l'ensemble des modes propagatifs\n");
		fprintf(stderr,"valeur minimale : %d \n",Nlimit_min+Nsigma_min);
		if (Nmin_super < 0)  Nmin_super = 0;
		if (Nmin_sub < 0)    Nmin_sub = 0;
	}
	if (Nlimit_max > Nsigma-1) {
		fprintf(stderr,"ATTENTION, N_sigma_max trop petit pour représenter l'ensemble des modes propagatifs\n");
		fprintf(stderr,"valeur minimale : %d \n",Nlimit_max+Nsigma_min);/* +1 ????? */
		if (Nmax_super > Nsigma-1) Nmax_super = Nsigma-1;
		if (Nmax_sub > Nsigma-1) Nmax_sub = Nsigma-1;
	}
	eff->Nmin_super = Nmin_super;
	eff->Nmax_super = Nmax_super;
	eff->Nmin_sub = Nmin_sub;
	eff->Nmax_sub = Nmax_sub;


	/* Calcul de l'énergie incidente */
	double Ei = 0.0;
	for (n=0; n<=Nsigma-1; n++) {
		/*alpha_n0 = csqrt(k02 - n*n*delta_sigma2);*/ /* SANS angle_i */
		alpha_n0 = csqrt(k02 - (sigma[n] + sigma0)*(sigma[n] + sigma0));
		Ei += creal(Ai[n]*conj(Ai[n])*alpha_n0);
	}

	/* Calcul des efficacités, ordres et angles en réflexion */
	for (n=Nmin_super; n<=Nmax_super; n++) {
		alpha_n0 = csqrt(k02 - (sigma[n] + sigma0)*(sigma[n] + sigma0));
		eff->eff_R[n-Nmin_super] = creal( A0[n]*conj(A0[n])*alpha_n0/Ei );
		eff->N_eff_R[n-Nmin_super] = (double) (n+Nsigma_min);
		eff->theta_eff_R[n-Nmin_super] = (180.0/PI)*asin((sigma0+sigma[n])/k0);
	}
	
	/* Calcul des efficacités, ordres et angles en transmission */
	if (par->pola == TM) coef_eff_T = creal((par->n_super/par->n_sub)*(par->n_super/par->n_sub));
	if (par->pola == TE) coef_eff_T = 1;
	for (n=Nmin_sub; n<=Nmax_sub; n++) {
		alpha_nh = csqrt(kh2 - (sigma[n] + sigma0)*(sigma[n] + sigma0));
		eff->eff_T[n-Nmin_sub] = creal( Ah[n]*conj(Ah[n])*alpha_nh/Ei )*coef_eff_T;
		eff->N_eff_T[n-Nmin_sub] = (double) (n+Nsigma_min);
		eff->theta_eff_T[n-Nmin_sub] = (180.0/PI)*asin((sigma0 + sigma[n])/kh);
	}

	/* Calcul des sommes des efficacités */
	eff->somm_eff_R = 0;
	eff->somm_eff_T = 0;
	for (n=eff->Nmin_super; n<=eff->Nmax_super; n++) {
		eff->somm_eff_R += eff->eff_R[n-eff->Nmin_super];
	}
	for (n=eff->Nmin_sub; n<=eff->Nmax_sub; n++) {
		eff->somm_eff_T += eff->eff_T[n-eff->Nmin_sub];
	}
	eff->somm_eff = eff->somm_eff_R + eff->somm_eff_T;

	return 0;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_amplitudes(complex *Ai, complex *A0, complex *Ah, complex **S12, 
                                complex **S22, struct Param_struct *par)
 *
 *	\brief	Calcul des amplitudes
 */
/*-------------------------------------------------------------------------------------*/
int md1D_amplitudes(complex *Ai, complex *A0, complex *Ah, complex **S12, complex **S22, struct Param_struct *par)
{
	int Nsigma = par->N_sigma;
	
	/* Calcul de A0 */
	A0 = M_x_V (A0, S12, Ai, Nsigma, Nsigma);

	/* Calcul de Ah */	
	Ah = M_x_V (Ah, S22, Ai, Nsigma, Nsigma);

	return 0;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn		int matrice_S(struct Param_struct *par)
 *
 *	\brief	Calcul de la matrice S
 *
 *	\todo	Passer les allocations en initialisation (?)
 */
/*-------------------------------------------------------------------------------------*/
int matrice_S(struct Param_struct *par){

	if (par->verbose) fprintf(stdout,"Calcul de la matrice S\n");

	int Nsigma = par->N_sigma;
	int NS = par->NS;
	complex **S12, **S22;
	S12 = par->S12;
	S22 = par->S22;

	int i, j, nS;
	complex **T11, **T12, **T21, **T22, **Z, **tmp, **tmp2;
	complex *F_plus, *F_moins, *F_plus2, *F_moins2;
	
	/* Allocations */
	T11 = allocate_CplxMatrix(Nsigma, Nsigma);
	T12 = allocate_CplxMatrix(Nsigma, Nsigma);
	T21 = allocate_CplxMatrix(Nsigma, Nsigma);
	T22 = allocate_CplxMatrix(Nsigma, Nsigma);
	Z   = allocate_CplxMatrix(Nsigma, Nsigma);
	tmp = allocate_CplxMatrix(Nsigma, Nsigma);
	tmp2= allocate_CplxMatrix(Nsigma, Nsigma);
	F_plus   = (complex *) malloc(sizeof(complex)*2*Nsigma);
	F_moins  = (complex *) malloc(sizeof(complex)*2*Nsigma);
	F_plus2  = (complex *) malloc(sizeof(complex)*2*Nsigma);
	F_moins2 = (complex *) malloc(sizeof(complex)*2*Nsigma);
	
	/* Initialisations */
	for (i=0; i<=Nsigma-1; i++) {
		for (j=0; j<=Nsigma-1; j++) {
			S12[i][j] = 0;
			S22[i][j] = (i==j);
		}
	}
	/* Itérations */
	for (nS=NS; nS>=1; nS--) {

			/* Calcul de la matrice T */
			if (par->pola == TE){
				matrice_T_TE(T11, T12, T21, T22, F_plus, F_moins, F_plus2, F_moins2, nS, par);
			}else{
				matrice_T_TM(T11, T12, T21, T22, F_plus, F_moins, F_plus2, F_moins2, nS, par);
			}

		/* Z = inv(T11 + T12*S12) */
		invM(Z, add_M(tmp2, 
			T11, M_x_M(tmp,
				T12,S12,	Nsigma,Nsigma),Nsigma,Nsigma),Nsigma);
		/* S12 = (T21 +T22*S12)*Z */
		M_x_M(S12,
			add_M(tmp2, T21, M_x_M(tmp,
					T22,S12,Nsigma,Nsigma),Nsigma,Nsigma),
			Z,Nsigma, Nsigma);
		/* S22 = S22*Z */
		M_egal(tmp,S22,Nsigma,Nsigma);
		M_x_M(S22,tmp,Z,Nsigma,Nsigma);

		/* Eventuellement, enregistrement de la matrice S en cours d'itérations */
		if (par->mode_extract_S){
			if (fmod((NS+1.0-nS)*par->h/NS,par->h_extract_S) < par->h/NS){
				md1D_save_S_matrix((NS+1.0-nS)*par->h/NS, par);
printf("extract S, h = %f\n",(NS+1.0-nS)*par->h/NS);
			}
		}
	}
	
	if(par->verbose >0) {fprintf(stdout,"\n");}
	
	/* Libération de la mémoire */
	free(T11[0]); free(T11); free(T12[0]); free(T12);
	free(T21[0]); free(T21); free(T22[0]); free(T22);
	free(Z[0]); free(Z); free(tmp[0]); free(tmp);
	free(tmp2[0]); free(tmp2);
	free(F_plus); free(F_moins); free(F_plus2); free(F_moins2);
	
	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn	int matrice_S_aleat_T(struct Param_struct *par)
 *
 *	\brief	Calcul de la matrice S par combinaison aléatoire de matrices T
 *
 */
/*-------------------------------------------------------------------------------------*/
int matrice_S_aleat_T(struct Param_struct *par){

	if (par->verbose) fprintf(stdout,"Calcul de la *** matrice S aleat ***\n");

	int Nsigma = par->N_sigma;
	int NS = par->NS;
	int NS_total = par->NS_total;
		complex **S12, **S22;
	S12 = par->S12;
	S22 = par->S22;

	int i, j, nS;
	complex ***T11, ***T12, ***T21, ***T22, **Z, **tmp, **tmp2;
	complex *F_plus, *F_moins, *F_plus2, *F_moins2;

	/* Allocations */
	T11 = allocate_CplxMatrix_3(Nsigma, Nsigma, NS);
	T12 = allocate_CplxMatrix_3(Nsigma, Nsigma, NS);
	T21 = allocate_CplxMatrix_3(Nsigma, Nsigma, NS);
	T22 = allocate_CplxMatrix_3(Nsigma, Nsigma, NS);
	Z   = allocate_CplxMatrix(Nsigma, Nsigma);
	tmp = allocate_CplxMatrix(Nsigma, Nsigma);
	tmp2= allocate_CplxMatrix(Nsigma, Nsigma);
	F_plus   = (complex *) malloc(sizeof(complex)*2*Nsigma);
	F_moins  = (complex *) malloc(sizeof(complex)*2*Nsigma);
	F_plus2  = (complex *) malloc(sizeof(complex)*2*Nsigma);
	F_moins2 = (complex *) malloc(sizeof(complex)*2*Nsigma);

	
		/*----- Calcul des matrices T de bases -----*/
	/* Itérations */
	for (nS=NS; nS>=1; nS--) {
			/* Calcul de la matrice T */
			if (par->pola == TE){
				matrice_T_TE(T11[nS-1], T12[nS-1], T21[nS-1], T22[nS-1], F_plus, F_moins, F_plus2, F_moins2, nS, par);
			}else{
				matrice_T_TM(T11[nS-1], T12[nS-1], T21[nS-1], T22[nS-1], F_plus, F_moins, F_plus2, F_moins2, nS, par);
			}
	}
	
		/*----- Composition de la matrice S à l'aide de la séquence de matrices T -----*/
	/* Initialisation de l'algorithme */
	for (i=0; i<=Nsigma-1; i++) {
		for (j=0; j<=Nsigma-1; j++) {
			S12[i][j] = 0;
			S22[i][j] = (i==j);
		}
	}
	fprintf(stdout,"\nIncrémentation de la matrice S ...\n");
	/* Itérations */
	for (i=NS_total-1; i>=0; i--) {

		nS = par->sequence_T[i];
				
		/* Z = inv(T11 + T12*S12) */
		invM(Z, add_M(tmp2, 
			T11[nS], M_x_M(tmp,
				T12[nS],S12,	Nsigma,Nsigma),Nsigma,Nsigma),Nsigma);
		/* S12 = (T21 +T22*S12)*Z */
		M_x_M(S12,
			add_M(tmp2, T21[nS], M_x_M(tmp,
					T22[nS],S12,Nsigma,Nsigma),Nsigma,Nsigma),
			Z,Nsigma, Nsigma);
		/* S22 = S22*Z */
		M_egal(tmp,S22,Nsigma,Nsigma);
		M_x_M(S22,tmp,Z,Nsigma,Nsigma);

		/* Eventuellement, enregistrement de la matrice S en cours d'itérations */
		if (par->mode_extract_S){
			if (fmod((NS_total+1.0-i)*par->h_total_aleat_T/NS_total,par->h_extract_S) < par->h_total_aleat_T/NS_total){
				md1D_save_S_matrix((NS_total+1.0-i)*par->h_total_aleat_T/NS_total, par);
printf("extract S, h = %f\n",(NS_total+1.0-i)*par->h_total_aleat_T/NS_total);
			}
		}

		/*Affichage du temps restant à l'écran */
		md1D_affichTemps(0,0,i,par->NS_total,0,1,par);
	
	}
		
	if(par->verbose >0) {fprintf(stdout,"\n");}
	
	/* Libération de la mémoire */
	free(T11[0][0]); free(T11[0]); free(T11);
	free(T12[0][0]); free(T12[0]); free(T12);
	free(T21[0][0]); free(T21[0]); free(T21);
	free(T22[0][0]); free(T22[0]); free(T22);
	free(Z[0]); free(Z);
	free(tmp[0]); free(tmp);
	free(tmp2[0]); free(tmp2);
	free(F_plus);
	free(F_moins);
	free(F_plus2);
	free(F_moins2);
	
	return 0;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn		int matrice_T_TE(complex **T11, complex **T12, complex **T21, complex **T22,
				complex *F_plus, complex *F_moins, complex *F_plus2, complex *F_moins2,
				int nS, struct Param_struct *par)
 *
 *	\brief Calcul de la matrice T en polarisation TE
 */
/*-------------------------------------------------------------------------------------*/
int matrice_T_TE(complex **T11, complex **T12, complex **T21, complex **T22,
				complex *F_plus, complex *F_moins, complex *F_plus2, complex *F_moins2,
				int nS, struct Param_struct *par)
{

	int Nsigma = par->N_sigma;
	int NS = par->NS;
	double h = par->h;
	complex k02 = par->k0*par->k0;
	complex kh2 = par->kh*par->kh;
	/*double delta_sigma2 = par->delta_sigma*par->delta_sigma;*/
	/*double delta_sigma = par->delta_sigma;*/
	double sigma0 = par->sigma0;
	
	int j, n;
	double hmin, hmax;
	complex L0, alpha_n;
	complex *E_plus, *dE_plus, *E_moins, *dE_moins;
	complex *E_plus2, *dE_plus2, *E_moins2, *dE_moins2;
	
	/* Les vecteurs F sont constitués de la superposition des vecteurs E et dE */
	E_plus   =  F_plus;
	dE_plus  = &F_plus[Nsigma];
	E_moins  =  F_moins;
	dE_moins = &F_moins[Nsigma];
	
	E_plus2   =  F_plus2;
	dE_plus2  = &F_plus2[Nsigma];
	E_moins2  =  F_moins2;
	dE_moins2 = &F_moins2[Nsigma];

		
	for (n=0; n<=Nsigma-1; n++){

		/* Valeur de alpha_n */
		if (nS==NS){ /* 1ere itération Matrice S : On est dans le substrat */
			/*alpha_n = csqrt(kh2- n*n*delta_sigma2);*/ /* SANS angle_i */
			alpha_n = csqrt(kh2- (par->sigma[n] + sigma0)*(par->sigma[n] + sigma0));
		}else{     /* Itérations suivantes : Meme materiau que le superstrat */
			/*alpha_n = csqrt(k02- n*n*delta_sigma2);*/ /* SANS angle_i */
			alpha_n = csqrt(k02- (par->sigma[n] + sigma0)*(par->sigma[n] + sigma0));
		}

		/* Construction des vecteurs E_plus et E_moins */	
		for (j=0;j<=Nsigma-1;j++){
			E_plus[j]   = 0;
			dE_plus[j]  = 0;
			E_moins[j]  = 0;
			dE_moins[j] = 0;
		}
		E_plus[n]   =  1;
		dE_plus[n]  =  I*alpha_n;
		E_moins[n]  =  1;
		dE_moins[n] =  -I*alpha_n;

		
   		/* Intégration des grandeurs pour une couche */
		hmin = h*(nS-1)/NS;
		hmax = h*nS/NS;
/*		
eq_diff((double *)F_plus,  (double *)F_plus2,  8*N+4, hmax, hmin, par->Nstep_S, fun_TE, (void *) par);
eq_diff((double *)F_moins, (double *)F_moins2, 8*N+4, hmax, hmin, par->Nstep_S, fun_TE, (void *) par);
*/

ode_solve((double *)F_plus,  (double *)F_plus2,  4*Nsigma, hmax, hmin, par->Nstep_S, fun_TE, (void *) par);
ode_solve((double *)F_moins, (double *)F_moins2, 4*Nsigma, hmax, hmin, par->Nstep_S, fun_TE, (void *) par);

		/* Construction des matrices T11, T21, T12 et T22*/
		for(j=0;j<=Nsigma-1;j++){
			L0 = 1.0/(I*csqrt(k02 - (par->sigma[j] + sigma0)*(par->sigma[j] + sigma0)));
			T11[j][n] = ( E_plus2[j]  + L0*dE_plus2[j] )*0.5;
			T21[j][n] = ( E_plus2[j]  - L0*dE_plus2[j] )*0.5;
			T12[j][n] = ( E_moins2[j] + L0*dE_moins2[j] )*0.5;
			T22[j][n] = ( E_moins2[j] - L0*dE_moins2[j] )*0.5;
		}
		
		/*Affichage du temps restant à l'écran */
		md1D_affichTemps(n,Nsigma,nS,NS,par->ni,par->Ni,par);
	}
	
	
	return 0;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn		int matrice_T_TM(complex **T11, complex **T12, complex **T21, complex **T22,
				complex *F_plus, complex *F_moins, complex *F_plus2, complex *F_moins2,
				int nS, struct Param_struct *par)
 *
 *	\brief Calcul de la matrice T en polarisation TM
 */
/*-------------------------------------------------------------------------------------*/
int matrice_T_TM(complex **T11, complex **T12, complex **T21, complex **T22,
				complex *F_plus, complex *F_moins, complex *F_plus2, complex *F_moins2,
				int nS, struct Param_struct *par)
{
	

	int Nsigma = par->N_sigma;
	int NS = par->NS;
	double h = par->h;
	complex k02 = par->k0*par->k0;
	complex kh2 = par->kh*par->kh;
	complex k0h2;
	/*double delta_sigma2 = par->delta_sigma*par->delta_sigma;*/
	/*double delta_sigma = par->delta_sigma;*/
	double sigma0 = par->sigma0;
		
	int j, n;
	double hmin, hmax;
	complex L0, alpha_n;
	complex *Eb_plus, *H_plus, *Eb_moins, *H_moins;
	complex *Eb_plus2, *H_plus2, *Eb_moins2, *H_moins2;
	
	/* Les vecteurs [F] sont constitués de la superposition des vecteurs [E'] et [H] : [F] = |[E']| */
	/* avec [E'] = [E]/(i.w.mu)                                                              |[H ]| */
	Eb_plus   =  F_plus;
	H_plus  = &F_plus[Nsigma];
	Eb_moins  =  F_moins;
	H_moins = &F_moins[Nsigma];
	
	Eb_plus2   =  F_plus2;
	H_plus2  = &F_plus2[Nsigma];
	Eb_moins2  =  F_moins2;
	H_moins2 = &F_moins2[Nsigma];

		
	for (n=0; n<=Nsigma-1; n++){
		
		/* Valeur de alpha_n */
		if (nS==NS){ /* 1ere itération Matrice S : On est dans le substrat */
			/*alpha_n = csqrt(kh2- n*n*delta_sigma2);*/ /* SANS angle_i */
			alpha_n = csqrt(kh2- (par->sigma[n] + sigma0)*(par->sigma[n] + sigma0));
			k0h2 = kh2;
		}else{     /* Itérations suivantes : Meme materiau que le superstrat */
			/*alpha_n = csqrt(k02- n*n*delta_sigma2);*/ /* SANS angle_i */
			alpha_n = csqrt(k02- (par->sigma[n] + sigma0)*(par->sigma[n] + sigma0));
			k0h2 = k02;
		}

		/* Construction des vecteurs E_plus et E_moins */	
		for (j=0;j<=Nsigma-1;j++){
			Eb_plus[j]  = 0;
			H_plus[j]   = 0;
			Eb_moins[j] = 0;
			H_moins[j]  = 0;
		}
		Eb_plus[n]  = -I*alpha_n/k0h2;
		H_plus[n]   =  1;
		Eb_moins[n] =  I*alpha_n/k0h2;
		H_moins[n]  =  1;

		
   		/* Intégration des grandeurs pour une couche */
		hmin = h*(nS-1)/NS;
		hmax = h*nS/NS;
/*		
eq_diff((double *)F_plus,  (double *)F_plus2,  8*N+4, hmax, hmin, par->Nstep_S, fun_TM, (void *) par);
eq_diff((double *)F_moins, (double *)F_moins2, 8*N+4, hmax, hmin, par->Nstep_S, fun_TM, (void *) par);
*/

ode_solve((double *)F_plus,  (double *)F_plus2,  4*Nsigma, hmax, hmin, par->Nstep_S, fun_TM, (void *) par);
ode_solve((double *)F_moins, (double *)F_moins2, 4*Nsigma, hmax, hmin, par->Nstep_S, fun_TM, (void *) par);

		/* Construction des matrices T11, T21, T12 et T22*/
		for(j=0;j<=Nsigma-1;j++){
			L0 = I*k02/csqrt(k02 - (par->sigma[j] + sigma0)*(par->sigma[j] + sigma0));
			T11[j][n] = ( H_plus2[j]  + L0*Eb_plus2[j] )*0.5;
			T21[j][n] = ( H_plus2[j]  - L0*Eb_plus2[j] )*0.5;
			T12[j][n] = ( H_moins2[j] + L0*Eb_moins2[j] )*0.5;
			T22[j][n] = ( H_moins2[j] - L0*Eb_moins2[j] )*0.5;
		}
		
		/*Affichage du temps restant à l'écran */
		md1D_affichTemps(n,Nsigma,nS,NS,par->ni,par->Ni,par);
	}
	
	
	return 0;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn		int fun_TE (double z, const double *F_reel, double *dF_reel, void *param_void)
 *
 *	\brief	Calcul la dérivée de F dans le cas TE
 */
/*-------------------------------------------------------------------------------------*/
int fun_TE (double z, const double *F_reel, double *dF_reel, void *param_void)
{
	int n, m, sizeE;
	struct Param_struct *par = (struct Param_struct *) param_void;
	int Nsigma = par->N_sigma;
	/*double delta_sigma2 = par->delta_sigma*par->delta_sigma;*/
	double delta_sigma = par->delta_sigma;
	double sigma0 = par->sigma0;
	
	complex *F, *dF, *null;
	complex *TF_k2 = par->TF_k2;

	CAST_COMPLEX(F, F_reel, Nsigma);
	CAST_COMPLEX(dF, dF_reel, Nsigma);

	if (par->STOCKER_TF) {
		FFT_k2_et_invk2_stockee(z, &TF_k2, &null, par);
	}else{
		TF_k2 = TF_k2_directe(z, TF_k2, par);
	}

	
	/* Calcul de [dF] à partir de [F] */
/*	for (n=-N; n<=N; n++) {
		dF[n+N] = F[n+N+sizeE];
		dF[n+N+sizeE] = (n*delta_sigma + sigma0)*(n*delta_sigma + sigma0)*F[n+N];
		for (m=-N; m<=N; m++) {
			dF[n+N+sizeE] -= delta_sigma*TF_k2[n-m+2*N]*F[m+N];
		}
	}
*/
	/* Delta_sigma quelconque */
	int Nmin2=par->N_sig_TF_min;
	double *sigma=par->sigma;
	sizeE = Nsigma;
	
/*	for (n=Nmin; n<=Nmax; n++) {
		dF[n-Nmin] = F[n-Nmin+sizeE];
		dF[n-Nmin+sizeE] = (n*delta_sigma + sigma0)*(n*delta_sigma + sigma0)*F[n-Nmin];
		for (m=Nmin; m<=Nmax; m++) {
			dF[n-Nmin+sizeE] -= delta_sigma*TF_k2[n-m-Nmin2]*F[m-Nmin];
		}
	}*/
	for (n=0; n<=Nsigma-1; n++) {
		dF[n] = F[n+sizeE];
		dF[n+sizeE] = (sigma[n] + sigma0)*(sigma[n] + sigma0)*F[n];
		for (m=0; m<=Nsigma-1; m++) {
			dF[n+sizeE] -= delta_sigma*TF_k2[n-m-Nmin2]*F[m];
		}
	}
	
/*	int register n_p_sizeE, n_p_2N;
	for (n=0; n<sizeE; n++) {
		n_p_sizeE = n + sizeE;
		n_p_2N = n + 2*N;
		dF[n] = F[n+sizeE];
		dF[n_p_sizeE] = ((n-N)*delta_sigma + sigma0)*((n-N)*delta_sigma + sigma0)*F[n];
		for (m=0; m<sizeE; m++) {
			dF[n_p_sizeE] -= TF_k2[n_p_2N-m]*F[m];
		}
	}
*/
	FREE_IF_CPP(F);
	UNCAST_COMPLEX(dF, dF_reel);
		
	return 0;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn		int fun_TM (double z, const double *F_reel, double *dF_reel, void *param_void)
 *
 *	\brief	Calcul la dérivée de F dans le cas TM
 */
/*-------------------------------------------------------------------------------------*/
int fun_TM (double z, const double *F_reel, double *dF_reel, void *param_void)
{
	int n, m, sizeE;
	struct Param_struct *par = (struct Param_struct *) param_void;
	int Nsigma = par->N_sigma;
	/*double delta_sigma2 = par->delta_sigma*par->delta_sigma;*/
	double delta_sigma = par->delta_sigma;
	double sigma0 = par->sigma0;
	
	complex *F, *dF;
	complex *TF_k2    = par->TF_k2;
	complex *TF_invk2 = par->TF_invk2;
	
	CAST_COMPLEX(F, F_reel, Nsigma);
	CAST_COMPLEX(dF, dF_reel, Nsigma);

	/* ESSAYER AVEC POINTEUR DE FONCTION POUR GAGNER DU TEMPS*/
	if (par->STOCKER_TF) {
		FFT_k2_et_invk2_stockee(z, &TF_k2, &TF_invk2, par);
	}else{
		TF_k2 = TF_k2_directe(z, TF_k2, par);
		TF_invk2 = TF_invk2_directe(z, TF_invk2, par);
	}

	
	/* Calcul de [dF] à partir de [F] */
/*	sizeE = 2*N+1;
	for (n=-N; n<=N; n++) {
		dF[n+N] = F[n+N+sizeE];
		dF[n+N+sizeE] = 0;
		for (m=-N; m<=N; m++) {
			dF[n+N] -= (n*delta_sigma + sigma0)*(m*delta_sigma + sigma0)*TF_invk2[n-m+2*N]*F[m+N+sizeE];
			dF[n+N+sizeE] -= TF_k2[n-m+2*N]*F[m+N];
		}
	}
*/
	/* Calcul de [dF] à partir de [F], Delta_sigma quelconque */
	int Nmin=par->N_sigma_min;
	int Nmax=par->N_sigma_max;
	int Nmin2=par->N_sig_TF_min;
	sizeE = Nsigma;
	for (n=Nmin; n<=Nmax; n++) {
		dF[n-Nmin] = F[n-Nmin+sizeE];
		dF[n-Nmin+sizeE] = 0;
		for (m=Nmin; m<=Nmax; m++) {
			dF[n-Nmin] -= delta_sigma*(n*delta_sigma + sigma0)*(m*delta_sigma + sigma0)*TF_invk2[n-m-Nmin2]*F[m-Nmin+sizeE];
			dF[n-Nmin+sizeE] -= delta_sigma*TF_k2[n-m-Nmin2]*F[m-Nmin];
		}
		/* OPTIMISATION : multiplier ici dF[n+N] par n puis ajouter F[n+N+sizeE] */
	}
	

	FREE_IF_CPP(F);
	UNCAST_COMPLEX(dF, dF_reel);
	
	return 0;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn		
 *
 *	\brief	Détermine le tableau de complexes k^2(x) pour un z donné
 *
 *	\todo	Prend pour l'instant en compte seulement un profil de type h(x)\n
 *			Doit être plus polyvalent : accepter aussi les profils de type n(x,z)
 */
/*-------------------------------------------------------------------------------------*/
int k2_H_X(struct Param_struct *par, complex *k2_1D, double z)
{
	int i;
	
	for (i=0;i<=par->N_x-1;i++){
		if (z < par->profil[0][i]) 
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
			if (z >= par->profil[n_layer][nx]){
				if (z <= par->profil[n_layer+1][nx]){
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
/*!	\fn	k2_N_XYZ(struct Param_struct *par, complex *k2_1D, double z)	
 *
 *	\brief	Détermine le tableau de complexes k^2(x) pour un z donné
 */
/*-------------------------------------------------------------------------------------*/
int k2_N_XYZ(struct Param_struct *par, complex *k2_1D, double z)
{
	int i, nz;
	double DeuxPisurLambda2 = (2*PI/par->lambda)*(2*PI/par->lambda);
	
	nz = (int)((z*par->N_z)/par->h);
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
	
	nz = (int)((z*par->N_z)/par->h);
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
		if (z < par->profil[0][i]) 
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
			if (z >= par->profil[n_layer][nx]){
				if (z <= par->profil[n_layer+1][nx]){
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



/*____________________________________________________________________________________*/
/*!	\fn		complex *TF_k2_directe(double z, complex *TF_k2, struct Param_struct *par)
 *
 *	\brief	Calcule la TF de k^2(x) pour un z donné et la tronque entre -N et +N
 *	ATTENTION NOUVELLE VERSION :  R E V O L U T I O N   ! ! !
 *	... prend maintenant un pas Delta_sigma quelconque !
 *	du coup n'utilise plus de FFT ...
 *	mais de betes TF "a la main" !
 *
 * 	\todo MULTIPLIER par Delta_sigma dans TF_k2 et pas dans fun_TE/M
 */
/*_____________________________________________________________________________________*/
complex *TF_k2_directe(double z, complex *TF_k2, struct Param_struct *par)
{

	/* Calcul de k^2(x) à z fixé, à partir du profil */
	(*par->k_2)(par, par->k2, z);

	/* Calcul de la TF de k2, avec N_x points */
	TF(TF_k2, par->k2, par->N_x, par->sigma_TF, par->N_sig_TF, par->L);
	
	return TF_k2;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		complex *TF_invk2_directe(double z, complex *TF_invk2, struct Param_struct *par)
 *
 *	\brief	Calcule la TF de 1/k^2(x) pour un z donné et la tronque entre -N et +N
 */
/*-------------------------------------------------------------------------------------*/
complex *TF_invk2_directe(double z, complex *TF_invk2, struct Param_struct *par)
{

	/* Calcul de 1 / k^2(x) à z fixé, à partir du profil */
	(*par->invk_2)(par, par->invk2, z);

	/* Calcul de la TF de invk2, avec N_x points */
	TF(TF_invk2, par->invk2, par->N_x, par->sigma_TF, par->N_sig_TF, par->L);
	
	return TF_invk2;
}


#if 0
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
#endif


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
		par->tab_TF_k2 = reallocate_CplxMatrix(par->tab_TF_k2,(par->TAILLE_z2n + par->BLOC_TAILLE_z2n),par->N_sig_TF);
		z2n[0] = (long int *) realloc(z2n[0], sizeof(long int)*(par->TAILLE_z2n + par->BLOC_TAILLE_z2n));
		z2n[1] = (long int *) realloc(z2n[1], sizeof(long int)*(par->TAILLE_z2n + par->BLOC_TAILLE_z2n));
	}
	/* Mise à jour de z2n : insertion de z et de n dans le tableau décroissant en z */
	int k = par->N_z2n;
	par->N_z2n++; 
	while(z_int > z2n[0][k-1] && k>0){
		z2n[0][k] = z2n[0][k-1];
		z2n[1][k] = z2n[1][k-1];
		k--;
	}
	z2n[0][k] = z_int;
	z2n[1][k] = par->N_z2n-1;
	/* Calcul de la TF de k2(z) et mise à jour de tab_FFT_k2 */
	TF_k2_directe(z, par->tab_TF_k2[par->N_z2n-1], par);

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
		par->tab_TF_k2    = reallocate_CplxMatrix(par->tab_TF_k2,(par->TAILLE_z2n + par->BLOC_TAILLE_z2n),par->N_sig_TF);
		par->tab_TF_invk2 = reallocate_CplxMatrix(par->tab_TF_invk2,(par->TAILLE_z2n + par->BLOC_TAILLE_z2n),par->N_sig_TF);
		z2n[0] = (long int *) realloc(z2n[0], sizeof(long int)*(par->TAILLE_z2n + par->BLOC_TAILLE_z2n));
		z2n[1] = (long int *) realloc(z2n[1], sizeof(long int)*(par->TAILLE_z2n + par->BLOC_TAILLE_z2n));
	}
	/* Mise à jour de z2n : insertion de z et de n dans le tableau décroissant en z */
	int k = par->N_z2n;
	par->N_z2n++; 
	while(z_int > z2n[0][k-1] && k>0){
		z2n[0][k] = z2n[0][k-1];
		z2n[1][k] = z2n[1][k-1];
		k--;
	}
	z2n[0][k] = z_int;
	z2n[1][k] = par->N_z2n-1;
	/* Calcul de la TF de k2(z) et mise à jour de tab_FFT_k2 */
	TF_k2_directe(z, par->tab_TF_k2[par->N_z2n-1], par);
	TF_invk2_directe(z, par->tab_TF_invk2[par->N_z2n-1], par);

	/* Renvoie des résultats */
	*ptTF_k2    = par->tab_TF_k2[par->N_z2n-1];
	*ptTF_invk2 = par->tab_TF_invk2[par->N_z2n-1];
	return 0;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn		long int *cherche(long int z, long int *tab0, int N)
 *
 *	\brief	Recherche dichotomique de l'élément z dans un tableau tab0 de taille N \n
 *			classé dans l'ORDRE DECROISSANT. 
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
	
	if (tab0[N>>1] < z) return(cherche(z, tab0, N>>1));
	else                return(cherche(z, tab0+(N>>1), N-(N>>1)));
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn		int md1D_affichTemps(int n, int N, int nS, int NS, struct Param_struct *par)
 *
 *	\brief	Affichage du temps restant estimé en cours de calculs
 */
/*-------------------------------------------------------------------------------------*/
int md1D_affichTemps(int n, int Nsigma, int nS, int NS, int ni, int Ni, struct Param_struct *par)
{

	/* Si moins de 5 secondes depuis le dernier affichage, on ne change rien */
	if (CHRONO(clock(), par->last_clock) < 5){
		return 0;
	/* Sinon, estimation et affichage de la durée restante */
	}else /*if (par->VERBOSE >= 1)*/{
		int i;
		time(&par->last_time);
		par->last_clock = clock();
		float t_ecoule = difftime(par->last_time,par->time0);
		float t_total; 
		t_total = t_ecoule*( Ni*NS*Nsigma)/( ni*NS*Nsigma+(NS-nS)*Nsigma+n );
		float t_restant = t_total - t_ecoule;
		int pourcent = ROUND(100.0*t_ecoule/t_total);
				
		fprintf(stdout,"\r");
		if (par->Ni > 1){
			 fprintf(stdout,"%3d %%, i = %d° [%ds ", pourcent,ROUND(par->angle_i*180/PI),ROUND(t_ecoule));
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
/*!	\fn	int md1D_save_S_matrix(double h_partial, struct Param_struct *par)
 *
 *	\brief	Sauvegarde la matrice S dans un fichier
 */
/*-------------------------------------------------------------------------------------*/
int md1D_save_S_matrix(double h_partial, struct Param_struct *par)
{
	int i,j,Nlign,Ncol;
	FILE *fp;
	char filename[SIZE_STR_BUFFER];
	Nlign = par->N_sigma;
	Ncol  = par->N_sigma;
	
	/* Nom de fichier */
	if (par->pola == TM){
		sprintf(filename,"S_TM_%s_h%f.txt",par->nom_profil,h_partial);
	}else{
		sprintf(filename,"S_TE_%s_h%f.txt",par->nom_profil,h_partial);
	}	
	fp = fopen(filename, "w");

	/* Ecriture de certains parametres */
	fprintf(fp,"mat_S_name = %s_h%f\n",par->nom_profil,h_partial);
	fprintf(fp,"Nsigma = %d\n",par->N_sigma);
	fprintf(fp,"h_partial = %1.6e\n",h_partial);
	fprintf(fp,"Re_k0 = %1.6e\n",creal(par->k0));
	fprintf(fp,"Im_k0 = %1.6e\n",cimag(par->k0));
	fprintf(fp,"Re_kh = %1.6e\n",creal(par->kh));
	fprintf(fp,"Im_kh = %1.6e\n",cimag(par->kh));
	fprintf(fp,"L = %1.6e\n",par->L);
	fprintf(fp,"delta_sigma = %1.6e\n",par->delta_sigma);
	fprintf(fp,"sigma0 = %1.6e\n",par->sigma0);
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

int TF(complex *TF, complex *f, int Nf, double *sigma, int Nsigma, double L)
{
	int i,j;
	double x;
	double delta_x = L/Nf;
	double coef = delta_x/(2*PI);	
	for (i=0;i<=Nsigma-1;i++){
		TF[i] = 0;
		for (j=0;j<=Nf-1;j++){
			x = delta_x * j;
			TF[i] += coef * f[j] * cexp(-I*sigma[i]*x);
		}
	}

	return 0;
}
		
		
		
