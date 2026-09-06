/* \file  md1D_in_out.h
 *  \brief Fichier d'en-tête pour le programme md1D_in_out
 */

#ifndef _MD1D_IN_OUT_H
#define _MD1D_IN_OUT_H

#include "std_include.h"

/* Constantes */
#define CHAR_COMMENT '#'

/* entrees_sorties*/
int lire_int(FILE *fp, const char *label, int *value);
int lire_double(FILE *fp, const char *label, double *value);
int lire_string(FILE *fp, const char *label, char *value);
int lire_complex(FILE *fp, const char *label, complex *value);
int md1D_lire_profil_H_X(const char *nom_fichier, struct Param_struct *par);
int md1D_lire_profil_MULTI(const char *nom_fichier, struct Param_struct *par);
int md1D_lire_profil_N_XYZ(const char *nom_fichier, struct Param_struct *par);
int md1D_lire_param(struct Noms_fichiers *nomfichier, struct Param_struct *par);
int md1D_lire_config(const char *nom_fichier, char *fichier_param);
int md1D_affiche_valeurs_param(struct Param_struct *par, struct Noms_fichiers *nomfichier);
int md1D_ecrire_results(char *filename, struct Param_struct *par, struct Efficacites_struct *eff);
int md1D_genere_nom_fichier_results(char *nomfichier_results, struct Param_struct *par);
int	md1D_ecrire_config(struct Param_struct *par, struct Noms_fichiers *nomfichier);
int	ecrire_dble_tab(FILE *fp, double *tab, int N, char *separateur1, int Nmax1, char *separateur2);
char *label_search(char *str,const char *label);
int lire_tab(const char *filename, const char *label, double *tab, int N);
int lire_ligne(FILE *fp, char *line);
void skip_comment(char *string);
int lire_str_arg(char *dest, char *label, int argc, char **argvcp);
int lire_dble_arg(double *res, char *label, int argc, char **argvcp);
int lire_int_arg(int *res, char *label, int argc, char **argvcp);

#endif /* _MD1D_IN_OUT_H */

