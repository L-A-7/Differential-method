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
int md2D_efficacites(COMPLEX *Ai, COMPLEX *A0, COMPLEX *Ah, struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_amplitudes(COMPLEX *Ai, COMPLEX *A0, COMPLEX *Ah, COMPLEX **S12, COMPLEX **S22, struct Param_struct *par);
int T_Matrix(COMPLEX **T11, COMPLEX **T12, COMPLEX **T21, COMPLEX **T22, int nS, struct Param_struct *par);
int invk_2(struct Param_struct *par, COMPLEX *invk2_1D, double z);
int k_2(struct Param_struct *par, COMPLEX *k2_1D, double z);
COMPLEX *FFT_invk2_directe(double z, COMPLEX *TF_invk2, struct Param_struct *par);
COMPLEX *FFT_k2_directe(double z, COMPLEX *TF_k2, struct Param_struct *par);
int md2D_QMatrix(double z, COMPLEX **Toep_k2, COMPLEX **invToep_invk2, COMPLEX **Qxx, COMPLEX **Qyy, COMPLEX **Qxz, COMPLEX **Qzz, COMPLEX **Qzz_1, struct Param_struct *par);
int md2D_zinvarQMatrix(double z, COMPLEX **Toep_k2, COMPLEX **invToep_invk2, struct Param_struct *par);
int Normal_H_X(struct Param_struct *par, COMPLEX *Nx2, COMPLEX *NxNz, COMPLEX *Nz2, double z);
int md2D_affichTemps(int n, int N, int nS, int NS, int ni, int Ni, struct Param_struct *par);
int md2D_save_S_matrix(double h_partial, struct Param_struct *par);
int rcwa_P_matrix(COMPLEX **P, double z, double Delta_z, struct Param_struct *par);
int fun (double z, const double *F_reel, double *dF_reel, void *param_void);
int md2D_save_near_field(COMPLEX **S12, COMPLEX **Z, int vec_size, int nS, struct Param_struct *par);

/* obsoletes ? */
int ode_solve(const double *y0, double *y, int N, double t0, double t1, int nstep,
	int (*f)(double, const double *, double *, void *), void *param_void);
int eq_diff(const double *y0, double *y, int N, double t0, double t1, int N_step,
	int (*func)(double, const double *, double *, void *), void *param_void);

#endif /* _md2D_H */



