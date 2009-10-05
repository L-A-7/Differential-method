/* 	\file  md3D.h
 *  	\brief Header file for md3D
 */

#ifndef _md3D_H
#define _md3D_H

#include "std_include.h"

#include <fftw3.h>

#include "md3D_io_utils.h"
#include "md3D_in_out.h"
#include "md3D_utils.h"


/* fonctions */
int md3D_efficacites(complex *Ai, complex *A0, complex *Ah, struct Param_struct *par, struct Efficacites_struct *eff);
int md3D_amplitudes(complex *Ai, complex *A0, complex *Ah, complex **S12, complex **S22, struct Param_struct *par);
int T_Matrix(complex **T11, complex **T12, complex **T21, complex **T22, int nS, struct Param_struct *par);
complex *invk_2(struct Param_struct *par, complex *invk2_1D, double z);
complex *k_2(struct Param_struct *par, complex *k2_1D, double z);
int Normal_H_XY(complex **norm_x, complex **norm_y, complex **norm_z, double z, struct Param_struct *par);
int md3D_toepNorm(double z, struct Param_struct *par);
int md3D_affichTemps(int n, int N, int nS, int NS, struct Param_struct *par);
int md3D_save_S_matrix(double h_partial, struct Param_struct *par);
int md3D_save_near_field(int nS, struct Param_struct* par);
int rcwa_P_matrix(complex **P, double z, double Delta_z, struct Param_struct *par);
int fun (double z, const double *F_reel, double *dF_reel, void *param_void);
complex **toeplitz_2D(complex **toep, int Nx, int Ny, complex *M_in, int Nxin, int Nyin);
int md3D_QMatrix(double z, complex **Qxx, complex **Qxy, complex **Qxz, complex **Qyy, complex **Qyz, complex **Qzz, complex **Qzz_1, struct Param_struct *par);

/* obsoletes ? */
int ode_solve(const double *y0, double *y, int N, double t0, double t1, int nstep,
	int (*f)(double, const double *, double *, void *), void *param_void);
int eq_diff(const double *y0, double *y, int N, double t0, double t1, int N_step,
	int (*func)(double, const double *, double *, void *), void *param_void);

#endif /* _md3D_H */



