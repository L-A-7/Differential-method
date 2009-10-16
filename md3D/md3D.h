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


/* functions md3D.c */
int md3D_propagativ_limits(struct Param_struct *par, struct Efficacites_struct *eff);
int md3D_incident_field(struct Param_struct *par, struct Efficacites_struct *eff);
int md3D_amplitudes(complex *Ai, complex *A0, complex *Ah, complex **S12, complex **S22, struct Param_struct *par);
int md3D_efficiencies(complex *Ai, complex *Ar, complex *At, struct Param_struct *par,  struct Efficacites_struct *eff);
int S_matrix(struct Param_struct *par);
int S_matrix_stack(struct Param_struct *par);
int T_Matrix(complex **T11, complex **T12, complex **T21, complex **T22, int nS, struct Param_struct *par);
int T_Matrix_homog_layer(complex **T11, complex **T12, complex **T21, complex **T22, complex nu_super, complex nu_sub, complex nu_layer, double h_layer, double lambda, complex *sigma_x, complex *sigma_y, int vec_size, struct Param_struct *par);
int zinvar_P_matrix(complex **P, double z, double Delta_z, struct Param_struct *par);
int rk4_P_matrix(complex **P, double z, double Delta_z, struct Param_struct *par);
int M_matrix(complex **M, double z, struct Param_struct *par);
int zinvar_M_matrix(complex **M, double z, struct Param_struct *par);
complex *invk_2(struct Param_struct *par, complex *invk2_1D, double z);
complex *k_2(struct Param_struct *par, complex *k2_1D, double z);
int Normal_H_XY(complex **norm_x, complex **norm_y, complex **norm_z, double z, struct Param_struct *par);
int Normal_N_XY_ZINVAR(complex **norm_x, complex **norm_y, complex **norm_z, double z, struct Param_struct *par);
int md3D_toepNorm(double z, struct Param_struct *par);
int md3D_affichTemps(int n, int N, int nS, int NS, struct Param_struct *par);
int md3D_save_S_matrix(double h_partial, struct Param_struct *par);
int md3D_save_near_field(int nS, struct Param_struct* par);
complex **toeplitz_2D(complex **toep, int Nx, int Ny, complex *M_in, int Nxin, int Nyin);
int md3D_QMatrix(double z, complex **Qxx, complex **Qxy, complex **Qxz, complex **Qyy, complex **Qyz, complex **Qzz, complex **Qzz_1, struct Param_struct *par);
int md3D_make_tab_S_steps(struct Param_struct* par);


/* Obsolete */ /*
int md3D_efficacites(complex *Ai, complex *A0, complex *Ah, struct Param_struct *par, struct Efficacites_struct *eff);
int rcwa_P_matrix(complex **P, double z, double Delta_z, struct Param_struct *par);
int fun (double z, const double *F_reel, double *dF_reel, void *param_void);
*/

/* functions of md3D_pilot.c */
int md3D_std (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md3D_init (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md3D_variables_init(struct Param_struct *par, struct Efficacites_struct *eff);
int md3D_arrays_init(struct Param_struct *par, struct Efficacites_struct *eff);
int md3D_alloc(struct Param_struct *par, struct Efficacites_struct *eff);
int md3D_alloc_init_profil(struct Param_struct *par);
int md3D_free(struct Param_struct *par, struct Efficacites_struct *eff);


/* fonctions attribuées à des pointeurs de fonctions */
complex *k2_homog(struct Param_struct *par, complex *k2_2D, double z);
complex *invk2_homog(struct Param_struct *par, complex *invk2_2D, double z);
complex *k2_H_XY  (struct Param_struct *par, complex *invk2_1D, double z);
complex *k2_MULTI(struct Param_struct *par, complex *invk2_1D, double z);
complex *invk2_H_XY(struct Param_struct *par, complex *invk2_1D, double z);
complex *invk2_MULTI(struct Param_struct *par, complex *invk2_1D, double z);
complex *k2_N_XYZ(struct Param_struct *par, complex *invk2_1D, double z);
complex *invk2_N_XYZ(struct Param_struct *par, complex *invk2_1D, double z);
complex **PsiMatrix(complex **Psi, complex k, complex *kz, complex *sigma_x, complex *sigma_y, int vec_size);

#endif /* _md3D_H */



