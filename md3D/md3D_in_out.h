/* \file  md3D_in_out.h
 *  \brief Fichier d'en-tête pour le programme md3D_in_out
 */

#ifndef _md3D_IN_OUT_H
#define _md3D_IN_OUT_H

#include "std_include.h"
#include "md_io_utils.h"
#include "md_utils.h"

/* Constantes */
#define AUTO -987325984

/* entrees_sorties*/
int md3D_lire_profil_H_XY(const char *nom_fichier, struct Param_struct *par);
int md3D_lire_profil_MULTI(const char *nom_fichier, struct Param_struct *par);
int md3D_lire_profil_N_XYZ(const char *nom_fichier, struct Param_struct *par);
int md3D_lire_param(struct Noms_fichiers *nomfichier, struct Param_struct *par);
int md3D_lire_config(const char *nom_fichier, char *param_file);
int md3D_affiche_valeurs_param(struct Param_struct *par, struct Noms_fichiers *nomfichier);
int md3D_ecrire_results(char *filename, struct Param_struct *par, struct Efficacites_struct *eff);
int md3D_genere_nom_fichier_results(char *nomfichier_results, struct Param_struct *par);
int md3D_ecrire_config(struct Param_struct *par, struct Noms_fichiers *nomfichier);
COMPLEX md3D_indice(char* name,double lambda, char* method);
COMPLEX md3D_index(char* name,double lambda, char* method);

#endif /* _md3D_IN_OUT_H */

