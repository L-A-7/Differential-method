#ifndef _STD_INCLUDE_H
#define _STD_INCLUDE_H

/* Bibliotheques standards */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#ifdef _ACML
#undef _BLAS_OPTIMIZATION
#undef _LAPACK_OPTIMIZATION
#endif

#ifdef _BLAS_OPTIMIZATION
#include <cblas.h>
#endif
#ifdef _LAPACK_OPTIMIZATION
#include <clapack.h>
#endif
#ifdef _ACML
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
/*#define free(x) fprintf(stdout,"%s, line%d, freeing %p\n",__FILE__,__LINE__,x);fflush(stdout);free(x)*/
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
	complex n_super;
	complex n_sub;
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
	complex *k2_layer;
	complex *invk2_layer;
	int N;
	int NS;
	int N_steps;
	int Ni;
	int ni;
	double Delta_sigma;
	complex sigma0;
	complex ky_0;
	complex k_super;
	complex k_sub;
	complex *sigma;
	complex *kz_super;
	complex *kz_sub;

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
	complex **n_xyz;
	int N_x;
	int N_z;
	int vec_size;
	int vec_middle;
	char nom_profil[SIZE_STR_BUFFER];
		
	complex **tab_TF_k2;
	complex **tab_TF_invk2;
	
	complex *k2;
	complex *invk2;
	complex *Nx2;
	complex *Nz2;
	complex *NxNz;
	
	complex *TF_k2;
	complex *TF_invk2;
	complex *TF_Nx2;
	complex *TF_Nz2;
	complex *TF_NxNz;
	int HX_Normal_CALCULATED;

	complex **T;
	complex **P;
	complex **M;
	complex *eig_values;
	complex **EigVectors;
	complex *eig_buffer;
	
	complex **Toep_k2;
	complex **Toep_invk2;
	complex **invToep_k2;
	complex **invToep_invk2;
	complex **Toep_Nx2;
	complex **Toep_NxNz;
	complex **Toep_Nz2;
	complex **M_tmp1;
	complex **M_tmp2;
	complex **M_tmp3;
	
	complex *tmp_tf_k2;   
	complex *tmp_tf_invk2;
	complex *tmp_tf_Nx2;
	complex *tmp_tf_Nz2;
	complex *tmp_tf_NxNz;

	complex **Qxx;
	complex **Qyy;
	complex **Qzz;
	complex **Qxz;
	complex **Qzz_1;

	complex *QxzEx;
	complex *Qzz_1QxzEx;
	complex *Qzz_1Hpx;
	complex *V_tmp1;
	complex *sigmaHpy; 
	complex *ky0Qzz_1Hpx;
	complex *Qzz_1sigmaHpy;
	complex *QxxEx;
	complex *QyyEy;
	complex *QxzVtmp1;

	complex **S12;
	complex **S22;
	complex **S11;
	complex **S21;

	complex **T11;
	complex **T12;
	complex **T21;
	complex **T22;
	
	complex **M_tmp11;
	complex **M_tmp12;
	complex **M_tmp21;
	complex **M_tmp22;

	complex **rkMz;
	complex **rkMzd;
	complex **rkMzdd;
	complex **rkM1;
	complex **rkM2;
	complex **rkM3;
	complex **rkM4;
	complex **rkMtmp1;

	complex **Psi_sub_TE;
	complex **Psi_super_TE;
	complex **invPsi_super_TE;
	complex **Psi_sub_TM;
	complex **Psi_super_TM;
	complex **invPsi_super_TM;
	complex **invEigVec;
	complex **invVec_Psi;
	complex **M_invVec_Psi;
	complex **M_buffer_2vecsize;
	complex *M_sol;
	
	complex *Ai;
	complex *Ar;
	complex *At;
	
	complex *Vi;
	complex *Vr;
	complex *Vt;

	complex *Exi;
	complex *Exr;
	complex *Ext;
	complex *Hpxi;
	complex *Hpxr;
	complex *Hpxt;
	
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
/*	complex **Near_field_matrix;
	complex *bottom_field;*/
	complex ***tab_Z;
	complex ***tab_S12;
		
	/* mode extract S */
	int mode_extract_S;
	double h_extract_S;
	
	/* Pointeurs de fonctions */
	int (*k_2)(struct Param_struct *par, complex *k2_1D, double z);
	int (*invk_2)(struct Param_struct *par, complex *invk2_1D, double z);
	int (*Normal_function)(struct Param_struct *par, complex *Nx2, complex *NxNz, complex *Nz2, double z);
	int (*md2D_lire_profil)(const char *, struct Param_struct *);
	int (*M_matrix)(complex **M, double z, struct Param_struct *par);
	int (*P_matrix)(complex **P, double z, double Delta_z, struct Param_struct *par);
	
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


