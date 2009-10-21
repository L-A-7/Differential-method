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
int md2D_efficacites(complex *Ai, complex *A0, complex *Ah, struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_amplitudes(complex *Ai, complex *A0, complex *Ah, complex **S12, complex **S22, struct Param_struct *par);
int T_Matrix(complex **T11, complex **T12, complex **T21, complex **T22, int nS, struct Param_struct *par);
int invk_2(struct Param_struct *par, complex *invk2_1D, double z);
int k_2(struct Param_struct *par, complex *k2_1D, double z);
complex *FFT_invk2_directe(double z, complex *TF_invk2, struct Param_struct *par);
complex *FFT_k2_directe(double z, complex *TF_k2, struct Param_struct *par);
int md2D_QMatrix(double z, complex **Toep_k2, complex **invToep_invk2, complex **Qxx, complex **Qyy, complex **Qxz, complex **Qzz, complex **Qzz_1, struct Param_struct *par);
int md2D_zinvarQMatrix(double z, complex **Toep_k2, complex **invToep_invk2, struct Param_struct *par);
int Normal_H_X(struct Param_struct *par, complex *Nx2, complex *NxNz, complex *Nz2, double z);
int md2D_affichTemps(int n, int N, int nS, int NS, int ni, int Ni, struct Param_struct *par);
int md2D_save_S_matrix(double h_partial, struct Param_struct *par);
int rcwa_P_matrix(complex **P, double z, double Delta_z, struct Param_struct *par);
int fun (double z, const double *F_reel, double *dF_reel, void *param_void);
int md2D_save_near_field(complex **S12, complex **Z, int vec_size, int nS, struct Param_struct *par);

/* obsoletes ? */
int ode_solve(const double *y0, double *y, int N, double t0, double t1, int nstep,
	int (*f)(double, const double *, double *, void *), void *param_void);
int eq_diff(const double *y0, double *y, int N, double t0, double t1, int N_step,
	int (*func)(double, const double *, double *, void *), void *param_void);

#endif /* _md2D_H */



