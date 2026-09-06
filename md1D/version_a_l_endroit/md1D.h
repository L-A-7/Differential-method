/* \file  md1D.h
 *  \brief Fichier d'en-tête pour le programme md1D
 */

#ifndef _MD1D_H
#define _MD1D_H

/* Bibliotheques standards */
#include "std_include.h"

#include <fftw3.h>

/* Bibliotheques specifiques */
#include "md1D_io_utils.h"
#include "md1D_in_out.h"
#include "md1D_utils.h"


/* fonctions */
int md1D_efficacites(complex *Ai, complex *A0, complex *Ah, struct Param_struct *par, struct Efficacites_struct *eff);
int md1D_amplitudes(complex *Ai, complex *A0, complex *Ah, complex **S12, complex **S22, struct Param_struct *par);
int matrice_T_TE(complex **T11, complex **T12, complex **T21, complex **T22,
				complex *F_plus, complex *F_moins, complex *F_plus2, complex *F_moins2,
				int nS, struct Param_struct *par);
int matrice_T_TM(complex **T11, complex **T12, complex **T21, complex **T22,
				complex *F_plus, complex *F_moins, complex *F_plus2, complex *F_moins2,
				int nS, struct Param_struct *par);
long int *cherche(long int z, long int *tab, int N);
/*int k2_H_X  (struct Param_struct *par, complex *invk2_1D, double z);
int k2_MULTI(struct Param_struct *par, complex *invk2_1D, double z);
int invk2_H_X(struct Param_struct *par, complex *invk2_1D, double z);
int invk2_MULTI(struct Param_struct *par, complex *invk2_1D, double z);*/
int invk_2(struct Param_struct *par, complex *invk2_1D, double z);
int k_2(struct Param_struct *par, complex *k2_1D, double z);
complex *FFT_invk2_directe(double z, complex *TF_invk2, struct Param_struct *par);
complex *FFT_k2_directe(double z, complex *TF_k2, struct Param_struct *par);
int FFT_k2_stockee(double z, complex **ptTF_k2, struct Param_struct *par);
int FFT_k2_et_invk2_stockee(double z, complex **ptTF_k2, complex **ptTF_invk2, struct Param_struct *par);
int fun_TE (double z, const double *F_reel, double *dF_reel, void *param_void);
int fun_TM (double z, const double *F_reel, double *dF_reel, void *param_void);
int ode_solve(const double *y0, double *y, int N, double t0, double t1, int nstep,
	int (*f)(double, const double *, double *, void *), void *param_void);
int eq_diff(const double *y0, double *y, int N, double t0, double t1, int N_step,
	int (*func)(double, const double *, double *, void *), void *param_void);
int md1D_affichTemps(int n, int N, int nS, int NS, int ni, int Ni, struct Param_struct *par);
int md1D_save_S_matrix(double h_partial, struct Param_struct *par);
int md1D_save_near_field(int nS, struct Param_struct* par);

#endif /* _MD1D_H */



