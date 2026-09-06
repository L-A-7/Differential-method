#ifndef _STD_INCLUDE_H
#define _STD_INCLUDE_H

/* Bibliotheques standards */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
/*#include <sys/time.h>*/

/* Complex.h */
#ifndef __cplusplus
#include <complex.h>
#define c_omplex(a,b) (a+I*b) 
#define CAST_COMPLEX(c,r,N) c=(complex *) r 
#define UNCAST_COMPLEX(c,r)
#define FREE_IF_CPP(c)
#else
#include "complex.h"
#define CAST_COMPLEX(c,r,N) { \
	int i##r; \
	c = (complex *) malloc(sizeof(complex)*N); \
	for (i##r=0; i##r<=N-1; i##r++) { \
		c[i##r].real = r[2*i##r]; \
		c[i##r].imag = r[2*i##r+1]; \
	} \
}while(0)
#define UNCAST_COMPLEX(c,r) { \
	int j##r; \
	for (j##r=0; j##r<=N-1; j##r++) { \
		r[2*j##r]   = c[j##r].real; \
		r[2*j##r+1] = c[j##r].imag; \
	} \
	free(c); \
}while(0)
#define FREE_IF_CPP(c) free(c)
#endif /* __cplusplus */


/* Constantes */
/*#define DBLE_CMP_EXIGEANCE 1000000*/
#define DBLE_CMP_EXIGEANCE 100000000.0
#define PI 3.14159265358979323846
#define TE 1
#define TM 2
#define H_X          1
#define MULTICOUCHES 2
#define N_XYZ        3
#define AUTO -987325984
#define GRATING -84674393

#define SIZE_STR_BUFFER 200
#define SIZE_LINE_BUFFER 50000

#define NON_LU "Et_non_c_pas_lu"
/* Macros */
#define CHRONO(t2,t1) ((double)(t2-t1)/CLOCKS_PER_SEC)
#define ROUND(x) ((int)(x<0 ? x-0.5 : x+0.5))
#define FLOOR(x) ((int)(x))
#define CEIL(x)  (x-(int)(x)>0 ? (int)(x)+1 : (int)(x))
#define MAX(a,b) ((a>b)?a:b)
#define MIN(a,b) ((a<b)?a:b)

/* Structures */
struct Param_struct {
	int argc;
	char **argvcp;
	int pola;
	char type_calcul[SIZE_STR_BUFFER];
	complex n_super;
	complex n_sub;
	double L;
	double h;
	double coef_h;
	double lambda;
	double angle_i;
	double k_sin_i;
	complex k0;
	complex kh;
	int type_profil;
	int N_layers;
/*	complex *indice;*/
	complex *k2_layer;
	complex *invk2_layer;
	int N;
	int NS;
	int Nstep;
	int Nstep_S;
	int Ni;
	int ni;
	double *sigma;
	double *sigma_TF;
	double delta_sigma;
	double sigma_min;
	double sigma_max;
	int N_sigma;
	int N_sigma_min;
	int N_sigma_max;
	int N_sig_TF;
	int N_sig_TF_min;
	int N_sig_TF_max;
	double sigma0;
	double delta_h;
	double **profil;
	complex **n_xyz;
	int N_x;
	int N_z;
	char nom_profil[SIZE_STR_BUFFER];
		
	complex **tab_TF_k2;
	complex **tab_TF_invk2;
	
	complex *k2;
	complex *invk2;
	complex *TF_k2;
	complex *TF_invk2;
	
	complex **S12;
	complex **S22;
	complex **S11;
	complex **S21;

	complex *Ai;
	complex *A0;
	complex *Ah;

	double *var_i;
	double *var_i2;
	
	double delta_s;
	double delta_p;
	double delta;

	/* variable spécifiques pour aleat_T_ellipso */
	double L_segment;
	double ecart_type_segment;
	double h_total_aleat_T;
	int *sequence_T;
	int NS_total;

	/* mode extract S */
	int mode_extract_S;
	double h_extract_S;
	
	/* Pointeurs de fonctions */
	int (*k_2)(struct Param_struct *par, complex *k2_1D, double z);
	int (*invk_2)(struct Param_struct *par, complex *invk2_1D, double z);
	int (*md1D_lire_profil)(const char *, struct Param_struct *);
	int (*matrice_T)(complex **, complex **, complex **, complex **, complex *, complex *, complex *,
				 complex *,	int, struct Param_struct *);
	
	long int *z2n[2];	 /* z2n             : "z to n", tableau de conversion de z vers n, utilisé pour  */
	int N_z2n;           /*                   le stockage des valeurs de TF_k2(z), permet de retrouver   */
	int TAILLE_z2n;      /*                   le numero de la ligne du tableau de stockage correspondant */
	int BLOC_TAILLE_z2n; /*                   a un z donné                                               */
	                     /* N_z2n           : Nombre d'éléments dans z2n                                 */
	                     /* TAILLE_z2n      : Taille alouée à z2n ainsi qu'à tab_FFT_k2                  */
	                     /* BLOC_TAILLE_z2n : valeur initiale de TAILLE_z2n et valeur qui lui est ajouté */
						 /*                   en cas de besoin de réallocation de mémoire                */					
	
	clock_t clock0;         /* Stocke le temps de départ                              */
	clock_t last_clock;     /* Durée écoulée depuis le dernier appel à md1D_temps     */
	time_t  time0;			/* clock(): très précis mais cyclique => tps longs pas OK */
	time_t last_time;		/* time() : précision d'1 s, non cyclique => tps longs OK */

	int READ_MAT_S; /* Si READ_MAT_S = 1, la matrice S est lue dans un fichier au lieu d'etre calculee */		
	char mat_S_file[SIZE_STR_BUFFER];
	char mat_S_name[SIZE_STR_BUFFER];
	int STOCKER_TF; /* Si STOCKER_TF = 0, calculs directs des TFs. si STOCKER_TF = 1 on stocke les  */
					/* valeurs des calculs de TFs : plus rapide, mais nécessite plus de mémoire.    */

	int verbose; /* Si verbose = 1, affiche plus d'infos sur le terminal */

	};

struct Efficacites_struct {
	int Nmin_super;
	int Nmax_super;
	int Nmin_sub;
	int Nmax_sub;

	double *eff_R;
	double *eff_T;
	
	double *N_eff_R;	
	double *N_eff_T;
	
	double *theta_eff_R;
	double *theta_eff_T;
	
	double somm_eff_T;
	double somm_eff_R;
	double somm_eff;
	};

struct Noms_fichiers {
	char fichier_config[SIZE_STR_BUFFER];
	char fichier_param[SIZE_STR_BUFFER];
	char fichier_profil[SIZE_STR_BUFFER];
	char fichier_results[SIZE_STR_BUFFER];
	};

#endif /* _STD_INCLUDE_H */


