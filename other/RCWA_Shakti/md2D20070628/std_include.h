#ifndef _STD_INCLUDE_H
#define _STD_INCLUDE_H

/* Replaced by -D in Makefile */
/*
#define BLAS_OPTIMIZATION   0
#define LAPACK_OPTIMIZATION 0
#define ACML_OPTIMIZATION 1
*/

/* Bibliotheques standards */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#if BLAS_OPTIMIZATION
#include <cblas.h>
#endif
#if LAPACK_OPTIMIZATION
#include <clapack.h>
#endif
#if ACML_OPTIMIZATION
#include <acml.h>
#endif

#include <complex.h>

/* Constantes */
#define DBLE_CMP_EXIGEANCE 100000000.0
#define PI 3.14159265358979323846
#define TE 1
#define TM 2
#define RE 1
#define IM 2

#define H_X          1
#define MULTICOUCHES 2
#define N_XYZ        3
#define SIZE_STR_BUFFER 200
#define SIZE_LINE_BUFFER 50000

#define NON_LU "Et_non_c_pas_lu"
/* Macros */
#define c_omplex(a,b) (a+I*b) 
#define CHRONO(t2,t1) ((double)(t2-t1)/CLOCKS_PER_SEC)
#define ROUND(x) ((int)(x<0 ? x-0.5 : x+0.5))
#define FLOOR(x) ((int)(x))
#define CEIL(x)  (x-(int)(x)>0 ? (int)(x)+1 : (int)(x))
#define MAX(a,b) ((a>b)?a:b)
#define MIN(a,b) ((a<b)?a:b)
#define SIGN(x)  (creal(x)/fabs(creal(x)))
#define CONJ(z)  (creal(z)-I*cimag(z))

/* For complex precision */
/* Using complex is not correct ('double complex' or 'float complex' must be used) */
#define REAL double

/* Structures */
struct Param_struct {
	int argc;
	char **argvcp;
	int pola;
	char calcul_type[SIZE_STR_BUFFER];
	char calcul_method[SIZE_STR_BUFFER];
	REAL complex n_super;
	REAL complex n_sub;
	double L;
	double h;
	double coef_h;
	double lambda;
	double theta_i;
	double phi_i;
	double psi;
/*	double k_sin_i;*/
	int type_profil;
	int N_layers;
	REAL complex *k2_layer;
	REAL complex *invk2_layer;
	int N;
	int NS;
	int Nstep;
	int Nstep_S;
	int Ni;
	int ni;
	double Delta_sigma;
	REAL complex sigma0;
	REAL complex ky_0;
	REAL complex k_super;
	REAL complex k_sub;
	REAL complex *sigma;
	REAL complex *kz_super;
	REAL complex *kz_sub;

	int READ_tab_NS;
        int tab_NS_ENABLED;
        double *tab_NS;
        char tab_NS_filename[SIZE_STR_BUFFER];

	int imposed_S_steps;
        int N_imposed_S_steps;
        double *tab_imposed_S_steps;
        char imposed_S_steps_filename[SIZE_STR_BUFFER];

	
	double delta_h;
	double **profil;
	REAL complex **n_xyz;
	int N_x;
	int N_z;
	int vec_size;
	int vec_middle;
	char nom_profil[SIZE_STR_BUFFER];
		
	REAL complex **tab_TF_k2;
	REAL complex **tab_TF_invk2;
	
	REAL complex *k2;
	REAL complex *invk2;
	REAL complex *Nx2;
	REAL complex *Nz2;
	REAL complex *NxNz;
	
	REAL complex *TF_k2;
	REAL complex *TF_invk2;
	REAL complex *TF_Nx2;
	REAL complex *TF_Nz2;
	REAL complex *TF_NxNz;
	int HX_Normal_CALCULATED;

	REAL complex **T;
	REAL complex **M;
	REAL complex *eig_values;
	REAL complex **EigVectors;
	REAL complex *eig_buffer;
	
	REAL complex **Toep_k2;
	REAL complex **Toep_invk2;
	REAL complex **invToep_k2;
	REAL complex **invToep_invk2;
	REAL complex **Toep_Nx2;
	REAL complex **Toep_NxNz;
	REAL complex **Toep_Nz2;
	REAL complex **M_tmp1;
	REAL complex **M_tmp2;
	
	REAL complex *tmp_tf_k2;   
	REAL complex *tmp_tf_invk2;
	REAL complex *tmp_tf_Nx2;
	REAL complex *tmp_tf_Nz2;
	REAL complex *tmp_tf_NxNz;

	REAL complex **Qxx;
	REAL complex **Qyy;
	REAL complex **Qzz;
	REAL complex **Qxz;
	REAL complex **Qzz_1;

	REAL complex *QxzEx;
	REAL complex *Qzz_1QxzEx;
	REAL complex *Qzz_1Hpx;
	REAL complex *V_tmp1;
	REAL complex *sigmaHpy; 
	REAL complex *ky0Qzz_1Hpx;
	REAL complex *Qzz_1sigmaHpy;
	REAL complex *QxxEx;
	REAL complex *QyyEy;
	REAL complex *QxzVtmp1;

	REAL complex **S12;
	REAL complex **S22;
	REAL complex **S11;
	REAL complex **S21;
	
	REAL complex **T11;
	REAL complex **T12;
	REAL complex **T21;
	REAL complex **T22;

	REAL complex **Psi_sub;
	REAL complex **Psi_super;
	REAL complex **invPsi_super;
	REAL complex **invEigVec;
	REAL complex **invVec_Psi;
	REAL complex **M_invVec_Psi;
	REAL complex **M_buffer_4vecsize;
	REAL complex *M_sol;
	
	REAL complex *Ai;
	REAL complex *Ar;
	REAL complex *At;
	
	REAL complex *Vi;
	REAL complex *Vr;
	REAL complex *Vt;

	REAL complex *Exi;
	REAL complex *Exr;
	REAL complex *Ext;
	REAL complex *Hpxi;
	REAL complex *Hpxr;
	REAL complex *Hpxt;
	
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

	/* NEAR_FIELD mode */
	REAL complex **Near_field_matrix;
	REAL complex *bottom_field;
		
	/* mode extract S */
	int mode_extract_S;
	double h_extract_S;
	
	/* Pointeurs de fonctions */
	int (*k_2)(struct Param_struct *par, REAL complex *k2_1D, double z);
	int (*invk_2)(struct Param_struct *par, REAL complex *invk2_1D, double z);
	int (*Normal_function)(struct Param_struct *par, REAL complex *Nx2, REAL complex *NxNz, REAL complex *Nz2, double z);
	int (*md2D_lire_profil)(const char *, struct Param_struct *);
	int (*T_Matrix)(REAL complex **T11, REAL complex **T12, REAL complex **T21, REAL complex **T22,
				int nS, struct Param_struct *par);
	
	long int *z2n[2];	 /* z2n             : "z to n", tableau de conversion de z vers n, utilisé pour  */
	int N_z2n;           /*                   le stockage des valeurs de TF_k2(z), permet de retrouver   */
	int TAILLE_z2n;      /*                   le numero de la ligne du tableau de stockage correspondant */
	int BLOC_TAILLE_z2n; /*                   a un z donné                                               */
	                     /* N_z2n           : Nombre d'éléments dans z2n                                 */
	                     /* TAILLE_z2n      : Taille alouée à z2n ainsi qu'à tab_FFT_k2                  */
	                     /* BLOC_TAILLE_z2n : valeur initiale de TAILLE_z2n et valeur qui lui est ajouté */
						 /*                   en cas de besoin de réallocation de mémoire                */					
	
	clock_t clock0;         /* Stocke le temps de départ                              */
	clock_t last_clock;     /* Durée écoulée depuis le dernier appel à md2D_temps     */
	time_t  time0;			/* clock(): très précis mais cyclique => tps longs pas OK */
	time_t last_time;		/* time() : précision d'1 s, non cyclique => tps longs OK */

	int READ_MAT_S; /* Si READ_MAT_S = 1, la matrice S est lue dans un fichier au lieu d'etre calculee */		
	char mat_S_file[SIZE_STR_BUFFER];
	char mat_S_name[SIZE_STR_BUFFER];
	int STOCKER_TF; /* Si STOCKER_TF = 0, calculs directs des TFs. si STOCKER_TF = 1 on stocke les  */
					/* valeurs des calculs de TFs : plus rapide, mais nécessite plus de mémoire.    */

	int verbosity; /* Si verbose = 1, affiche plus d'infos sur le terminal */
	char i_field_mode[SIZE_STR_BUFFER]; /* incident field (PLANE_WAVE, GAUSSIAN, FROM_BINARY, FROM_ASCII) */

	};

struct Efficacites_struct {
	int Nmin_super;
	int Nmax_super;
	int Nmin_sub;
	int Nmax_sub;

	double *eff_r;
	double *eff_t;
	
	double *N_eff_r;	
	double *N_eff_t;
	
	double *theta_r;
	double *theta_t;
	double *phi_r;
	double *phi_t;
	
	double sum_eff_r;
	double sum_eff_t;
	double sum_eff;
	};

struct Noms_fichiers {
	char fichier_config[SIZE_STR_BUFFER];
	char fichier_param[SIZE_STR_BUFFER];
	char fichier_profil[SIZE_STR_BUFFER];
	char fichier_results[SIZE_STR_BUFFER];
	};

#endif /* _STD_INCLUDE_H */


