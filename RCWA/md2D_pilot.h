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
int md2D_conical_FFF (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);int md2D_standard (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md2D_conical_FFF_ellipso (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md2D_variation_incidence_A (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier);
int md2D_variation_incidence_B (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier);
int md2D_ellipso (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier);
int md2D_aleat_T_ellipso (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier);
int md2D_var_i_ellipso (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier);
int md2D_var_lambda_ellipso (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier);
int md2D_init (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md2D_variables_init(struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_arrays_init(struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_alloc(struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_alloc_init_profil(struct Param_struct *par);
int md2D_free(struct Param_struct *par, struct Efficacites_struct *eff);

int md2D_efficiencies(REAL complex *Ai, REAL complex *Ar, REAL complex *At, struct Param_struct *par,  struct Efficacites_struct *eff);
int md2D_propagativ_limits(struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_incident_field(struct Param_struct *par, struct Efficacites_struct *eff);
int md2D_amplitudes(REAL complex *Ai, REAL complex *A0, REAL complex *Ah, REAL complex **S12, REAL complex **S22, struct Param_struct *par);
int S_matrix(struct Param_struct *par);
int matrice_S_aleat_T(struct Param_struct *par);
int md2D_comb_mat_S(REAL complex **S11, REAL complex **S12, REAL complex **S21, REAL complex **S22, 
		REAL complex **S11_1, REAL complex **S12_1, REAL complex **S21_1, REAL complex **S22_1, 
		REAL complex **S11_2, REAL complex **S12_2, REAL complex **S21_2, REAL complex **S22_2, int taille_matrice);
int md2D_read_mat_S(struct Param_struct *par);
int md2D_PsiMatrix(REAL complex **Psi, REAL complex k, REAL complex *kz, struct Param_struct *par);
int md2D_make_tab_S_steps(struct Param_struct* par);

/* fonctions attribuées à des pointeurs de fonctions */
int k2_H_X  (struct Param_struct *par, REAL complex *invk2_1D, double z);
int k2_MULTI(struct Param_struct *par, REAL complex *invk2_1D, double z);
int invk2_H_X(struct Param_struct *par, REAL complex *invk2_1D, double z);
int invk2_MULTI(struct Param_struct *par, REAL complex *invk2_1D, double z);
int k2_N_XYZ(struct Param_struct *par, REAL complex *invk2_1D, double z);
int invk2_N_XYZ(struct Param_struct *par, REAL complex *invk2_1D, double z);
int Normal_H_X(struct Param_struct *par, REAL complex *Nx2, REAL complex *NxNz, REAL complex *Nz2, double z);
int dm_T_Matrix(REAL complex **T11, REAL complex **T12, REAL complex **T21, REAL complex **T22, int nS, struct Param_struct *par);
int rcwa_T_Matrix(REAL complex **T11, REAL complex **T12, REAL complex **T21, REAL complex **T22, int nS, struct Param_struct *par);

#endif /* _md2D_H */

