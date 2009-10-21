/* \file  md2D_pilot.h
 *  \brief Fichier d'en-tête pour le programme md2D_pilot
 */

#ifndef _md2D_PILOT_H
#define _md2D_PILOT_H

/* Bibliotheques standards */
#include "std_include.h"

#include <fftw3.h>

/* Bibliotheques specifiques */
#include "md2D_io_utils.h"
#include "md2D_in_out.h"
#include "md2D_utils.h"


/* fonctions */
int md2D_classical_FFF (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);int md2D_standard (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md2D_conical_FFF_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md2D_near_field (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md2D_var_lambda_ellipso (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier);
int md2D_init (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md2D_variables_init(struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_arrays_init(struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_alloc(struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_alloc_init_profil(struct Param_struct *par);
int md2D_free(struct Param_struct *par, struct Efficacites_struct *eff);

int md2D_efficiencies(complex *Ai, complex *Ar, complex *At, struct Param_struct *par,  struct Efficacites_struct *eff);
int md2D_propagativ_limits(struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_incident_field(struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_amplitudes(complex *Ai, complex *A0, complex *Ah, complex **S12, complex **S22, struct Param_struct *par);
int S_matrix(struct Param_struct *par);
int matrice_S_aleat_T(struct Param_struct *par);
int md2D_comb_mat_S(complex **S11, complex **S12, complex **S21, complex **S22, 
		complex **S11_1, complex **S12_1, complex **S21_1, complex **S22_1, 
		complex **S11_2, complex **S12_2, complex **S21_2, complex **S22_2, int taille_matrice);
int md2D_read_mat_S(struct Param_struct *par);
int PsiMatrixTE(complex **Psi, complex k, complex *kz, struct Param_struct *par);
int PsiMatrixTM(complex **Psi, complex k, complex *kz, struct Param_struct *par);
int md2D_make_tab_S_steps(struct Param_struct* par);
int md2D_near_field_map(complex ***tab_S12, complex ***tab_Z, struct Param_struct *par);

/* fonctions attribuées à des pointeurs de fonctions */
int k2_H_X  (struct Param_struct *par, complex *invk2_1D, double z);
int k2_MULTI(struct Param_struct *par, complex *invk2_1D, double z);
int invk2_H_X(struct Param_struct *par, complex *invk2_1D, double z);
int invk2_MULTI(struct Param_struct *par, complex *invk2_1D, double z);
int k2_N_XYZ(struct Param_struct *par, complex *invk2_1D, double z);
int invk2_N_XYZ(struct Param_struct *par, complex *invk2_1D, double z);
int Normal_H_X(struct Param_struct *par, complex *Nx2, complex *NxNz, complex *Nz2, double z);
int M_matrix_TE(complex **M, double z, struct Param_struct *par);
int M_matrix_TM(complex **M, double z, struct Param_struct *par);
int zinvar_M_matrix_TM(complex **M, double z, struct Param_struct *par);
int zinvar_P_matrix(complex **P, double z, double Delta_z, struct Param_struct *par);
int rk4_P_matrix(complex **P, double z, double Delta_z, struct Param_struct *par);
int shooting_P_matrix(complex **P, double z, double Delta_z, struct Param_struct *par);

/* For testing purpose */
int euler_P_matrix(complex **P, double z, double Delta_z, struct Param_struct *par);
int euler2_P_matrix(complex **P, double z, double Delta_z, struct Param_struct *par);

#endif /* _md2D_H */

