/* \file  md3D_pilot.h
 *  \brief Fichier d'en-tête pour le programme md3D_pilot
 */

#ifndef _md3D_PILOT_H
#define _md3D_PILOT_H

/* Bibliotheques standards */
#include "std_include.h"

#include <fftw3.h>

/* Bibliotheques specifiques */
#include "md3D_io_utils.h"
#include "md3D_in_out.h"
#include "md3D_utils.h"


/* fonctions */
int md3D_std (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md3D_init (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md3D_variables_init(struct Param_struct *par, struct Efficacites_struct *eff);
int md3D_arrays_init(struct Param_struct *par, struct Efficacites_struct *eff);
int md3D_alloc(struct Param_struct *par, struct Efficacites_struct *eff);
int md3D_alloc_init_profil(struct Param_struct *par);
int md3D_free(struct Param_struct *par, struct Efficacites_struct *eff);

int md3D_efficiencies(complex *Ai, complex *Ar, complex *At, struct Param_struct *par,  struct Efficacites_struct *eff);
int md3D_propagativ_limits(struct Param_struct *par, struct Efficacites_struct *eff);
int md3D_incident_field(struct Param_struct *par, struct Efficacites_struct *eff);
int md3D_amplitudes(complex *Ai, complex *A0, complex *Ah, complex **S12, complex **S22, struct Param_struct *par);
int S_matrix(struct Param_struct *par);
int md3D_make_tab_S_steps(struct Param_struct* par);

/* fonctions attribuées à des pointeurs de fonctions */
complex *k2_H_XY  (struct Param_struct *par, complex *invk2_1D, double z);
complex *k2_MULTI(struct Param_struct *par, complex *invk2_1D, double z);
complex *invk2_H_XY(struct Param_struct *par, complex *invk2_1D, double z);
complex *invk2_MULTI(struct Param_struct *par, complex *invk2_1D, double z);
complex *k2_N_XYZ(struct Param_struct *par, complex *invk2_1D, double z);
complex *invk2_N_XYZ(struct Param_struct *par, complex *invk2_1D, double z);
int Normal_H_XY(complex **norm_x, complex **norm_y, complex **norm_z, double z, struct Param_struct *par);
int zinvar_P_matrix(complex **P, double z, double Delta_z, struct Param_struct *par);
int rk4_P_matrix(complex **P, double z, double Delta_z, struct Param_struct *par);
int M_matrix(complex **M, double z, struct Param_struct *par);
int zinvar_M_matrix(complex **M, double z, struct Param_struct *par);
complex **PsiMatrix(complex **Psi, complex k, complex *kz, complex *sigma_x, complex *sigma_y, int vec_size);
#endif /* _md3D_H */

