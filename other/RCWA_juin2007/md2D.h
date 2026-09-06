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
int T_matrix(complex **T11, complex **T12, complex **T21, complex **T22,
				complex *F1, complex *F2, int nS, struct Param_struct *par);
long int *cherche(long int z, long int *tab, int N);
int invk_2(struct Param_struct *par, complex *invk2_1D, double z);
int k_2(struct Param_struct *par, complex *k2_1D, double z);
complex *FFT_invk2_directe(double z, complex *TF_invk2, struct Param_struct *par);
complex *FFT_k2_directe(double z, complex *TF_k2, struct Param_struct *par);
int FFT_k2_stockee(double z, complex **ptTF_k2, struct Param_struct *par);
int FFT_k2_et_invk2_stockee(double z, complex **ptTF_k2, complex **ptTF_invk2, struct Param_struct *par);
int fun (double z, const double *F_real, double *dF_real, void *param_void);
int md2D_QMatrix(double z, complex **Qxx, complex **Qyy, complex **Qxz, complex **Qzz, complex **Qzz_1, struct Param_struct *par);
int Normal_H_X(struct Param_struct *par, complex *Nx2, complex *NxNz, complex *Nz2, double z);
int ode_solve(const double *y0, double *y, int N, double t0, double t1, int nstep,
	int (*f)(double, const double *, double *, void *), void *param_void);
int eq_diff(const double *y0, double *y, int N, double t0, double t1, int N_step,
	int (*func)(double, const double *, double *, void *), void *param_void);
int md2D_affichTemps(int n, int N, int nS, int NS, int ni, int Ni, struct Param_struct *par);
int md2D_save_S_matrix(double h_partial, struct Param_struct *par);
int md2D_save_near_field(int nS, struct Param_struct* par);
int rcwa_M_Matrix(complex **M, double z, struct Param_struct *par);
int T_matrix_RCWA(complex **T11, complex **T12, complex **T21, complex **T22,
				complex *F1, complex *F2, int nS, struct Param_struct *par);
#endif /* _md2D_H */



