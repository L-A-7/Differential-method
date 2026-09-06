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
int md1D_init (struct Param_struct *par, struct Efficacites_struct *eff,struct Noms_fichiers *nomfichier, int argc, char **argv);
int md1D_variables_init(struct Param_struct *par, struct Efficacites_struct *eff);
int md1D_alloc(struct Param_struct *par, struct Efficacites_struct *eff);
int md1D_free(struct Param_struct *par, struct Efficacites_struct *eff);

int md1D_efficacites(complex *Ai, complex *A0, complex *Ah, struct Param_struct *par, struct Efficacites_struct *eff);
int md1D_amplitudes(complex *Ai, complex *A0, complex *Ah, complex **S12, complex **S22, struct Param_struct *par);
int matrice_S(struct Param_struct *par);
int md1D_lire_args(char *fichier_profil, struct Param_struct *par, int argc, char **argv);
int md1D_affiche_valeurs_lues(struct Param_struct *par, struct Noms_fichiers *nomfichier);
int lire_str_arg(char *dest, char *label, int argc, char **argvcp);
int lire_dble_arg(double *res, char *label, int argc, char **argvcp);
int lire_int_arg(int *res, char *label, int argc, char **argvcp);

/* fonctions attribuées à des pointeurs de fonctions */

int md1D_lire_profil(const char *nom_fichier, struct Param_struct *par);
int md1D_lire_multi_profil(const char *nom_fichier, struct Param_struct *par);

int k2_H_X  (complex *k2_1D, double **profil, int N_profil, complex *k2_layer, double z);
int k2_MULTI(complex *k2_1D, double **profil, int N_profil, complex *k2_layer, double z);
int invk2_H_X  (complex *invk2_1D, double **profil, int N_profil, complex *invk2_layer, double z);
int invk2_MULTI(complex *invk2_1D, double **profil, int N_profil, complex *invk2_layer, double z);

int matrice_T_TE(complex **T11, complex **T12, complex **T21, complex **T22,
				complex *F_plus, complex *F_moins, complex *F_plus2, complex *F_moins2,
				int nS, struct Param_struct *par);
int matrice_T_TM(complex **T11, complex **T12, complex **T21, complex **T22,
				complex *F_plus, complex *F_moins, complex *F_plus2, complex *F_moins2,
				int nS, struct Param_struct *par);


#endif /* _MD1D_H */


