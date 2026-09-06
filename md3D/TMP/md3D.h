/* 	\file  md3D.h
 *  	\brief Header file for md3D
 */

#ifndef _md3D_H
#define _md3D_H

#include "std_include.h"

#include <fftw3.h>

#include "md_io_utils.h"
#include "md3D_in_out.h"
#include "md_utils.h"
#include "md_maths.h"


/* functions md3D.c */
int md3D_propagativ_limits(struct Param_struct *par, struct Efficacites_struct *eff);
int md3D_incident_field(struct Param_struct *par, struct Efficacites_struct *eff);
int md3D_amplitudes(COMPLEX *Ai, COMPLEX *A0, COMPLEX *Ah, COMPLEX **S12, COMPLEX **S22, struct Param_struct *par);
int md3D_efficiencies(COMPLEX *Ai, COMPLEX *Ar, COMPLEX *At, struct Param_struct *par,  struct Efficacites_struct *eff);
int S_matrix(struct Param_struct *par);
int S_matrix_stack(struct Param_struct *par);
int T_Matrix(COMPLEX **T11, COMPLEX **T12, COMPLEX **T21, COMPLEX **T22, int nS, struct Param_struct *par);
int T_Matrix_homog_layer(COMPLEX **T11, COMPLEX **T12, COMPLEX **T21, COMPLEX **T22, COMPLEX nu_super, COMPLEX nu_sub, COMPLEX nu_layer, double h_layer, double lambda, COMPLEX *sigma_x, COMPLEX *sigma_y, int vec_size, struct Param_struct *par);
int zinvar_P_matrix(COMPLEX **P, double z, double Delta_z, struct Param_struct *par);
int rk4_P_matrix(COMPLEX **P, double z, double Delta_z, struct Param_struct *par);
int M_matrix(COMPLEX **M, double z, struct Param_struct *par);
int zinvar_M_matrix(COMPLEX **M, double z, struct Param_struct *par);
COMPLEX *invk_2(struct Param_struct *par, COMPLEX *invk2_1D, double z);
COMPLEX *k_2(struct Param_struct *par, COMPLEX *k2_1D, double z);
int md3D_toepNorm(double z, struct Param_struct *par);
int md3D_affichTemps(int n, int N, int nS, int NS, struct Param_struct *par);
int md3D_save_S_matrix(double h_partial, struct Param_struct *par);
int md3D_save_near_field(int nS, struct Param_struct* par);
COMPLEX **toeplitz_2D(COMPLEX **toep, int Nx, int Ny, COMPLEX *M_in, int Nxin, int Nyin);
int md3D_QMatrix(double z, COMPLEX **Qxx, COMPLEX **Qxy, COMPLEX **Qxz, COMPLEX **Qyy, COMPLEX **Qyz, COMPLEX **Qzz, COMPLEX **Qzz_1, struct Param_struct *par);
int md3D_make_tab_S_steps(struct Param_struct* par);


/* Obsolete */ /*
int md3D_efficacites(COMPLEX *Ai, COMPLEX *A0, COMPLEX *Ah, struct Param_struct *par, struct Efficacites_struct *eff);
int rcwa_P_matrix(COMPLEX **P, double z, double Delta_z, struct Param_struct *par);
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
COMPLEX *k2_homog(struct Param_struct *par, COMPLEX *k2_2D, double z);
COMPLEX *invk2_homog(struct Param_struct *par, COMPLEX *invk2_2D, double z);
COMPLEX *k2_H_XY  (struct Param_struct *par, COMPLEX *invk2_1D, double z);
COMPLEX *k2_MULTI(struct Param_struct *par, COMPLEX *invk2_1D, double z);
COMPLEX *invk2_H_XY(struct Param_struct *par, COMPLEX *invk2_1D, double z);
COMPLEX *invk2_MULTI(struct Param_struct *par, COMPLEX *invk2_1D, double z);
COMPLEX *k2_N_XYZ(struct Param_struct *par, COMPLEX *invk2_1D, double z);
COMPLEX *invk2_N_XYZ(struct Param_struct *par, COMPLEX *invk2_1D, double z);
COMPLEX **PsiMatrix(COMPLEX **Psi, COMPLEX k, COMPLEX *kz, COMPLEX *sigma_x, COMPLEX *sigma_y, int vec_size);
int Normal_H_XY(COMPLEX **norm_x, COMPLEX **norm_y, COMPLEX **norm_z, double z, double *profil, struct Param_struct *par);
int Normal_N_XY_ZINVAR(COMPLEX **norm_x, COMPLEX **norm_y, COMPLEX **norm_z, double z, double *profil, struct Param_struct *par);
int Normal_MULTI(COMPLEX **norm_x, COMPLEX **norm_y, COMPLEX **norm_z, double z, double *profil, struct Param_struct *par);

#endif /* _md3D_H */



