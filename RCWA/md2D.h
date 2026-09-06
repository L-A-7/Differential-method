/* 	\file  md2D.h
 *  	\brief Header file for md2D
 */

#ifndef _md2D_H
#define _md2D_H

#include "std_include.h"

#include <fftw3.h>

#include "md2D_io_utils.h"
#include "md2D_in_out.h"
#include "md2D_utils.h"


/* fonctions */
int md2D_efficacites(REAL complex *Ai, REAL complex *A0, REAL complex *Ah, struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_amplitudes(REAL complex *Ai, REAL complex *A0, REAL complex *Ah, REAL complex **S12, REAL complex **S22, struct Param_struct *par);
int T_matrix(REAL complex **T11, REAL complex **T12, REAL complex **T21, REAL complex **T22,
				REAL complex *F1, REAL complex *F2, int nS, struct Param_struct *par);
long int *cherche(long int z, long int *tab, int N);
int invk_2(struct Param_struct *par, REAL complex *invk2_1D, double z);
int k_2(struct Param_struct *par, REAL complex *k2_1D, double z);
REAL complex *FFT_invk2_directe(double z, REAL complex *TF_invk2, struct Param_struct *par);
REAL complex *FFT_k2_directe(double z, REAL complex *TF_k2, struct Param_struct *par);
int FFT_k2_stockee(double z, REAL complex **ptTF_k2, struct Param_struct *par);
int FFT_k2_et_invk2_stockee(double z, REAL complex **ptTF_k2, REAL complex **ptTF_invk2, struct Param_struct *par);
int fun (double z, const double *F_real, double *dF_real, void *param_void);
int md2D_QMatrix(double z, REAL complex **Qxx, REAL complex **Qyy, REAL complex **Qxz, REAL complex **Qzz, REAL complex **Qzz_1, struct Param_struct *par);
int Normal_H_X(struct Param_struct *par, REAL complex *Nx2, REAL complex *NxNz, REAL complex *Nz2, double z);
int ode_solve(const double *y0, double *y, int N, double t0, double t1, int nstep,
	int (*f)(double, const double *, double *, void *), void *param_void);
int eq_diff(const double *y0, double *y, int N, double t0, double t1, int N_step,
	int (*func)(double, const double *, double *, void *), void *param_void);
int md2D_affichTemps(int n, int N, int nS, int NS, int ni, int Ni, struct Param_struct *par);
int md2D_save_S_matrix(double h_partial, struct Param_struct *par);
int md2D_save_near_field(int nS, struct Param_struct* par);
int rcwa_M_Matrix(REAL complex **M, double z, struct Param_struct *par);
int T_matrix_RCWA(REAL complex **T11, REAL complex **T12, REAL complex **T21, REAL complex **T22,
				REAL complex *F1, REAL complex *F2, int nS, struct Param_struct *par);
#endif /* _md2D_H */



