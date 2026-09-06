/* \file  md1D_pilot.h
 *  \brief Fichier d'en-tête pour le programme md1D_pilot
 */

#ifndef _MD1D_PILOT_H
#define _MD1D_PILOT_H

/* Bibliotheques standards */
#include "std_include.h"

#include <fftw3.h>

/* Bibliotheques specifiques */
#include "md1D_in_out.h"
#include "md1D_utils.h"


/* fonctions */
int md1D_standard (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md1D_variation_incidence_A (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier);
int md1D_variation_incidence_B (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier);
int md1D_ellipso (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier);
int md1D_aleat_T_ellipso (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier);
int md1D_var_i_ellipso (struct Param_struct *par, struct Efficacites_struct *eff, struct Noms_fichiers *nomfichier);
int md1D_init (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier);
int md1D_variables_init(struct Param_struct *par, struct Efficacites_struct *eff);
int md1D_alloc(struct Param_struct *par, struct Efficacites_struct *eff);
int md1D_alloc_init_profil(struct Param_struct *par);
int md1D_free(struct Param_struct *par, struct Efficacites_struct *eff);

int md1D_efficacites(complex *Ai, complex *A0, complex *Ah, struct Param_struct *par, struct Efficacites_struct *eff);
int md1D_amplitudes(complex *Ai, complex *A0, complex *Ah, complex **S12, complex **S22, struct Param_struct *par);
int matrice_S(struct Param_struct *par);
int matrice_S_aleat_T(struct Param_struct *par);
int md1D_comb_mat_S(complex **S11, complex **S12, complex **S21, complex **S22, 
		complex **S11_1, complex **S12_1, complex **S21_1, complex **S22_1, 
		complex **S11_2, complex **S12_2, complex **S21_2, complex **S22_2, int taille_matrice);
int md1D_read_mat_S(struct Param_struct *par);

/* fonctions attribuées à des pointeurs de fonctions */


int k2_H_X  (struct Param_struct *par, complex *invk2_1D, double z);
int k2_MULTI(struct Param_struct *par, complex *invk2_1D, double z);
int invk2_H_X(struct Param_struct *par, complex *invk2_1D, double z);
int invk2_MULTI(struct Param_struct *par, complex *invk2_1D, double z);
int k2_N_XYZ(struct Param_struct *par, complex *invk2_1D, double z);
int invk2_N_XYZ(struct Param_struct *par, complex *invk2_1D, double z);

int matrice_T_TE(complex **T11, complex **T12, complex **T21, complex **T22,
				complex *F_plus, complex *F_moins, complex *F_plus2, complex *F_moins2,
				int nS, struct Param_struct *par);
int matrice_T_TM(complex **T11, complex **T12, complex **T21, complex **T22,
				complex *F_plus, complex *F_moins, complex *F_plus2, complex *F_moins2,
				int nS, struct Param_struct *par);


#endif /* _MD1D_H */

