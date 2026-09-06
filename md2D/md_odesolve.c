/*!	\file		md_odesolve.c
 *
 * 	\brief		Vectorial Ordinary Differential Equation Solvers
 * 
 *
 *	\date		sept 2010
 *	\authors	Laurent ARNAUD
 */


#include "md_odesolve.h"

/*-------------------------------------------------------------------------------------*/
/*!	int zinvar_P_matrix(COMPLEX **P, double z, double Delta_z, struct Param_struct *par)
 *
 *	\brief P_matrix calculation in the case of z invariance. Uses eigen values resolution of the system (ie RCWA or Fourier Modal Method)
 */
/*-------------------------------------------------------------------------------------*/
int zinvar_P_matrix(COMPLEX **P, double z, double Delta_z, struct Param_struct *par)
{
	int i,j, vec_size = par->vec_size;
	COMPLEX *eig_values, *eig_buffer, *M_sol, **EigVectors, **invEigVec, **M_invVec_Psi;

	eig_values = (COMPLEX *) malloc(sizeof(COMPLEX)*2*vec_size);
	eig_buffer = (COMPLEX *) malloc(sizeof(COMPLEX)*50*2*vec_size);
	M_sol      = (COMPLEX *) malloc(sizeof(COMPLEX)*2*vec_size);
	EigVectors   = allocate_CplxMatrix(2*vec_size,2*vec_size);
	invEigVec    = allocate_CplxMatrix(2*vec_size,2*vec_size);
	M_invVec_Psi = allocate_CplxMatrix(2*vec_size,2*vec_size);

	/* M matrix calculation */
	(*par->M_matrix)(par->M, z, par);

	/* Diagonalisation of M */
	eigen_values(par->M, eig_values, EigVectors, eig_buffer, 2*vec_size);

	/* Solution of the diagonalized system */
	for (i=0;i<=2*vec_size-1;i++){
		M_sol[i] = cexp(eig_values[i]*Delta_z);
	}
		
	/* invEigVec = inv(EigVectors) */
	invM(invEigVec, EigVectors, 2*vec_size);

	/* M_invVec_Psi = M_sol * invVec_Psi */
	for (i=0;i<=2*vec_size-1;i++){
		for (j=0;j<=2*vec_size-1;j++){
			M_invVec_Psi[i][j] = M_sol[i] * invEigVec[i][j];
		}
	}

	/* P = EigVec * M_diag * inv(EigVec) */
	M_x_M(P, EigVectors, M_invVec_Psi, 2*vec_size, 2*vec_size);

	par->N_steps++;
		
	free(eig_values);
	free(eig_buffer);
	free(M_sol);
	free(EigVectors[0]);free(EigVectors);
	free(invEigVec[0]);free(invEigVec);
	free(M_invVec_Psi[0]);free(M_invVec_Psi);

	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	int rk4_P_matrix(COMPLEX **p, double z, double Delta_z, struct Param_struct *par)
 *
 *	\brief P_matrix calculation with 4th order Runge Kutta vectorized method
 *        integration occurs between z-dz/2 and z+dz/2
 */
/*-------------------------------------------------------------------------------------*/
int rk4_P_matrix(COMPLEX **P, double z, double dz, struct Param_struct *par)
{
	int i,j, vec_size = par->vec_size;
	COMPLEX **Mz, **Mzd, **Mzdd, **M2, **M3, **M4, **rkM_tmp;
	Mz    = allocate_CplxMatrix(2*vec_size,2*vec_size);
	Mzd   = allocate_CplxMatrix(2*vec_size,2*vec_size);
	Mzdd  = allocate_CplxMatrix(2*vec_size,2*vec_size);
/*	M1    = allocate_CplxMatrix(2*vec_size,2*vec_size);*/
	M2    = allocate_CplxMatrix(2*vec_size,2*vec_size);
	M3    = allocate_CplxMatrix(2*vec_size,2*vec_size);
	M4    = allocate_CplxMatrix(2*vec_size,2*vec_size);
	rkM_tmp = allocate_CplxMatrix(2*vec_size,2*vec_size);

	double dz_2 = 0.5*dz;

	/* M matrix calculation */
	(*par->M_matrix)(Mz,   z-dz_2, par);
	(*par->M_matrix)(Mzd,  z,      par);
	(*par->M_matrix)(Mzdd, z+dz_2, par);

	/* M1 = Mz */

	/* M2 = Mzd*(Id + M1*dz/2) */
	for (i=0;i<=2*vec_size-1;i++){
		for (j=0;j<=2*vec_size-1;j++){
			rkM_tmp[i][j] = dz_2 * Mz[i][j]; 
		}
		rkM_tmp[i][i] += 1; 
	}
	M_x_M(M2, Mzd, rkM_tmp, 2*vec_size, 2*vec_size);

	/* M3 = Mzd*(Id + M2*dz/2) */
	for (i=0;i<=2*vec_size-1;i++){
		for (j=0;j<=2*vec_size-1;j++){
			rkM_tmp[i][j] = dz_2 * M2[i][j]; 
		}
		rkM_tmp[i][i] += 1; 
	}
	M_x_M(M3, Mzd, rkM_tmp, 2*vec_size, 2*vec_size);

	
	/* M4 = Mzdd*(Id + M3*dz) */
	for (i=0;i<=2*vec_size-1;i++){
		for (j=0;j<=2*vec_size-1;j++){
			rkM_tmp[i][j] = dz * M3[i][j]; 
		}
		rkM_tmp[i][i] += 1; 
	}
	M_x_M(M4, Mzdd, rkM_tmp, 2*vec_size, 2*vec_size);

	/* P = Id + (h/6)M1 + (h/3)M2 + (h/3)M3 + (h/6)M4 */
	double dz_6 = dz/6.0;
	double dz_3 = dz/3.0;
	for (i=0;i<=2*vec_size-1;i++){
		for (j=0;j<=2*vec_size-1;j++){
			P[i][j] = dz_6*(Mz[i][j] + M4[i][j]) + dz_3*(M2[i][j] + M3[i][j]);
		}
		P[i][i] += 1; 
	}

	par->N_steps++;

	free(Mz[0]);
	free(Mzd[0]);
	free(Mzdd[0]);
/*	free(M1[0]); free(M1);*/
	free(M2[0]);
	free(M3[0]);
	free(M4[0]);
	free(rkM_tmp[0]);
	free(Mz);
	free(Mzd);
	free(Mzdd);
	free(M2);
	free(M3);
	free(M4);
	free(rkM_tmp);

	return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	int implicit_rk_P_matrix(COMPLEX **p, double z, double Delta_z, struct Param_struct *par)
 *
 *	\brief P_matrix calculation with implicit Runge Kutta vectorized method
 *  cf. Watanabe "Numerical integration schemes used on the differential theory for anisotropic gratings", JOSAA 2002
 */
/*-------------------------------------------------------------------------------------*/
int implicit_rk_P_matrix(COMPLEX **P, double z, double dz, struct Param_struct *par)
{
	int i,j, vec_size = par->vec_size;
	COMPLEX **Mtmp1,**Mtmp2,**Mtmp3,**Mtmp4,**Mtmp5,**Mz3m,**Mz3p,**invMz3m,**invMz3p,**K2,**invK2,**R1a,**R1b,**R1,**R2b,**R2;
	Mtmp1 = allocate_CplxMatrix(2*vec_size,2*vec_size);
	Mtmp2 = allocate_CplxMatrix(2*vec_size,2*vec_size);
	Mtmp3 = allocate_CplxMatrix(2*vec_size,2*vec_size);
	Mtmp4 = allocate_CplxMatrix(2*vec_size,2*vec_size);
	Mtmp5 = allocate_CplxMatrix(2*vec_size,2*vec_size);

	double dz_3m = dz*(3-sqrt(3))/6;
	double dz_3p = dz*(3+sqrt(3))/6;

	/* M matrix calculation */
	Mz3m=Mtmp3;
	Mz3p=Mtmp4;
	(*par->M_matrix)(Mz3m, z+dz_3m, par);
	(*par->M_matrix)(Mz3p, z+dz_3p, par);

	/* invMz3p=inv(Mz3p) */
	invMz3p = Mtmp1;
	invM(invMz3p,Mz3p,2*vec_size);

	/* invMz3m=inv(Mz3m) */
	invMz3m = Mtmp2;
	invM(invMz3m,Mz3m,2*vec_size);

	/* K2 = inv(Mz3m) - h*Id/4 */
	K2 = Mtmp3;
	for (i=0;i<=2*vec_size-1;i++){
		for (j=0;j<=2*vec_size-1;j++){
			K2[i][j] = invMz3p[i][j]; 
		}
		K2[i][i] -= dz/4; 
	}

	/* invK2 = inv(K2) */
	invK2 = Mtmp2;
	invM(invK2,K2,2*vec_size);

	/* R1a = (inv(Mz3p) -Id*h/4) + h^2*invK2/48 */
	R1a=Mtmp3;
	double dz2_48=dz*dz/48;
	for (i=0;i<=2*vec_size-1;i++){
		for (j=0;j<=2*vec_size-1;j++){
			R1a[i][j] = invMz3p[i][j]+dz2_48*invK2[i][j]; 
		}
		R1a[i][i] -= dz/4; 
	}
	/* R1b = Id + [(3+2*sqrt(3))*h/12]*invK2 */
	R1b=Mtmp1;
	double c1=(3+2*sqrt(3))*dz/12;
	for (i=0;i<=2*vec_size-1;i++){
		for (j=0;j<=2*vec_size-1;j++){
			R1b[i][j] = c1*invK2[i][j]; 
		}
	R1b[i][i] += 1;
	}

	/* R1 = inv(R1a) * R1b */
	R1=Mtmp4;
	M_x_M(R1, invM(Mtmp5,R1a,2*vec_size), R1b, 2*vec_size, 2*vec_size);

	/* R2b = Id + [(3-2*sqrt(3))*h/12]*R1 */
	R2b=Mtmp1;
	double c2 = (3-2*sqrt(3))*dz/12;
	for (i=0;i<=2*vec_size-1;i++){
		for (j=0;j<=2*vec_size-1;j++){
			R2b[i][j] = c2*R1[i][j]; 
		}
	R2b[i][i] += 1;
	}

	/* R2 = invK2 * R2b */
	R2=Mtmp3;
	M_x_M(R2, invK2, R2b, 2*vec_size, 2*vec_size);

	/* P = Id + (h/2)*(R1+R2) */
	for (i=0;i<=2*vec_size-1;i++){
		for (j=0;j<=2*vec_size-1;j++){
			P[i][j] = dz*(R1[i][j] + R2[i][j])/2;
		}
		P[i][i] += 1; 
	}

	par->N_steps++;

	free(Mtmp1[0]);free(Mtmp1);
	free(Mtmp2[0]);free(Mtmp2);
	free(Mtmp3[0]);free(Mtmp3);
	free(Mtmp4[0]);free(Mtmp4);
	free(Mtmp5[0]);free(Mtmp5);

	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	int euler_P_matrix(COMPLEX **P, double z, double Delta_z, struct Param_struct *par)
 *
 *	\brief P_matrix calculation with vectorized Euler method
 */
/*-------------------------------------------------------------------------------------*/
int euler_P_matrix(COMPLEX **P, double z, double dz, struct Param_struct *par)
{
	int i,j, vec_size = par->vec_size;
	COMPLEX **Mz;
	Mz    = allocate_CplxMatrix(2*vec_size,2*vec_size);

	/* M matrix calculation */
	(*par->M_matrix)(Mz,  z	, par);

	/* P = Id + h*M */
	for (i=0;i<=2*vec_size-1;i++){
		for (j=0;j<=2*vec_size-1;j++){
			P[i][j] = dz * Mz[i][j]; 
		}
		P[i][i] += 1; 
	}

	par->N_steps++;

	free(Mz[0]);
	free(Mz);
		
	return 0;
}
/*-------------------------------------------------------------------------------------*/
/*!	int euler2_P_matrix(COMPLEX **P, double z, double Delta_z, struct Param_struct *par)
 *
 *	\brief P_matrix calculation with vectorized Euler method (takes the M matrix at the middle of the interval)
 */
/*-------------------------------------------------------------------------------------*/
int euler2_P_matrix(COMPLEX **P, double z, double dz, struct Param_struct *par)
{
	int i,j, vec_size = par->vec_size;
	COMPLEX **Mzd;
	Mzd   = allocate_CplxMatrix(2*vec_size,2*vec_size);

	/* M matrix calculation */
	(*par->M_matrix)(Mzd,  z+0.5*dz	, par);

	/* P = Id + h*M(h+dz/2) */
	for (i=0;i<=2*vec_size-1;i++){
		for (j=0;j<=2*vec_size-1;j++){
			P[i][j] = dz * Mzd[i][j]; 
		}
		P[i][i] += 1; 
	}

	par->N_steps++;

	free(Mzd[0]);
	free(Mzd);
		
	return 0;
}


