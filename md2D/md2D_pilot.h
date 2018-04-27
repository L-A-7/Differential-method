/* \file  md2D_pilot.h
 *  \brief Header file for md2D_pilot
 */

#ifndef _md2D_PILOT_H
#define _md2D_PILOT_H

/* Bibliotheques standards */
#include "std_include.h"

#include <fftw3.h>

/* Bibliotheques specifiques */
#include "md_io_utils.h"
#include "md2D_in_out.h"
#include "md_utils.h"
#include "md_maths.h"


/* fonctions */
int md2D_classical_FFF (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md2D_standard (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md2D_guided (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md2D_swifts (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md2D_var_i (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md2D_conical_FFF_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md2D_near_field (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md2D_var_lambda_ellipso (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier);
int md2D_init (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md2D_variables_init(struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_arrays_init(struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_alloc(struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_alloc_init_profil(struct Param_struct *par);
int md2D_free(struct Param_struct *par, struct Efficacites_struct *eff);

int md2D_efficiencies(COMPLEX *Ai, COMPLEX *Ar, COMPLEX *At, struct Param_struct *par,  struct Efficacites_struct *eff);
int md2D_propagativ_limits(struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_incident_field(struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_amplitudes(COMPLEX *Ai, COMPLEX *A0, COMPLEX *Ah, COMPLEX **S12, COMPLEX **S22, struct Param_struct *par);
int S_matrix(struct Param_struct *par);
int matrice_S_aleat_T(struct Param_struct *par);
int md2D_comb_mat_S(COMPLEX **S11, COMPLEX **S12, COMPLEX **S21, COMPLEX **S22, 
		COMPLEX **S11_1, COMPLEX **S12_1, COMPLEX **S21_1, COMPLEX **S22_1, 
		COMPLEX **S11_2, COMPLEX **S12_2, COMPLEX **S21_2, COMPLEX **S22_2, int taille_matrice);
int md2D_read_mat_S(struct Param_struct *par);
int PsiMatrixTE(COMPLEX **Psi, COMPLEX k, COMPLEX *kz, struct Param_struct *par);
int PsiMatrixTM(COMPLEX **Psi, COMPLEX k, COMPLEX *kz, struct Param_struct *par);
int md2D_make_tab_S_steps(struct Param_struct* par);

/* fonctions attribuées à des pointeurs de fonctions */
int k2_H_X  (struct Param_struct *par, COMPLEX *invk2_1D, double z);
int k2_MULTI(struct Param_struct *par, COMPLEX *invk2_1D, double z);
int invk2_H_X(struct Param_struct *par, COMPLEX *invk2_1D, double z);
int invk2_MULTI(struct Param_struct *par, COMPLEX *invk2_1D, double z);
int k2_N_XYZ(struct Param_struct *par, COMPLEX *invk2_1D, double z);
int invk2_N_XYZ(struct Param_struct *par, COMPLEX *invk2_1D, double z);
int Normal_H_X(struct Param_struct *par, COMPLEX *Nx2, COMPLEX *NxNz, COMPLEX *Nz2, double z);
int M_matrix_TE(COMPLEX **M, double z, struct Param_struct *par);
int M_matrix_TM(COMPLEX **M, double z, struct Param_struct *par);
int zinvar_M_matrix_TM(COMPLEX **M, double z, struct Param_struct *par);
int zinvar_P_matrix(COMPLEX **P, double z, double Delta_z, struct Param_struct *par);
int rk4_P_matrix(COMPLEX **P, double z, double Delta_z, struct Param_struct *par);
int implicit_rk_P_matrix(COMPLEX **P, double z, double dz, struct Param_struct *par);
int shooting_P_matrix(COMPLEX **P, double z, double Delta_z, struct Param_struct *par);

/* For testing purpose */
int euler_P_matrix(COMPLEX **P, double z, double Delta_z, struct Param_struct *par);
int euler2_P_matrix(COMPLEX **P, double z, double Delta_z, struct Param_struct *par);

#endif /* _md2D_H */

