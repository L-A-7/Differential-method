/* 	\file  md2D.h
 *  	\brief Header file for md2D
 */

#ifndef _md2D_H
#define _md2D_H

#include "std_include.h"

#include <fftw3.h>

#include "md_io_utils.h"
#include "md2D_in_out.h"
#include "md_utils.h"
#include "md_maths.h"


/* fonctions */
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
int md2D_save_near_field(COMPLEX **S12, COMPLEX **Z, int vec_size, int nS, struct Param_struct *par);
int PsiMatrix(COMPLEX *Psi11, COMPLEX *Psi12, COMPLEX *Psi21, COMPLEX *Psi22, COMPLEX *kz, COMPLEX k, int vec_size, int pola);
int invPsiMatrix(COMPLEX *Psi11, COMPLEX *Psi12, COMPLEX *Psi21, COMPLEX *Psi22, COMPLEX *kz, COMPLEX k, int vec_size, int pola);
int md2D_make_tab_S_steps(struct Param_struct *par);
int read_S_steps_from_profile(struct Param_struct *par);
int md2D_local_field_map_by_T_products(COMPLEX *V0m, struct Param_struct *par);
int md2D_local_field_map(COMPLEX ***tab_S12, COMPLEX ***tab_Z, struct Param_struct *par);
int md2D_energy_flux(struct Param_struct *par, struct Efficacites_struct *eff);
int monolayer_S_matrix(COMPLEX n_guide, double h_guide, struct Param_struct *par);
int swifts_md2D_energy_flux(double h_guide, struct Param_struct *par, struct Efficacites_struct *eff);
int S_matrix_init(struct Param_struct *par);

/* Redondants */
/*
int md2D_amplitudes(COMPLEX *Ai, COMPLEX *A0, COMPLEX *Ah, COMPLEX **S12, COMPLEX **S22, struct Param_struct *par);
*/

/* obsoletes ? */
/*
int ode_solve(const double *y0, double *y, int N, double t0, double t1, int nstep,
	int (*f)(double, const double *, double *, void *), void *param_void);
int eq_diff(const double *y0, double *y, int N, double t0, double t1, int N_step,
	int (*func)(double, const double *, double *, void *), void *param_void);
int rcwa_P_matrix(COMPLEX **P, double z, double Delta_z, struct Param_struct *par);*/

#endif /* _md2D_H */



