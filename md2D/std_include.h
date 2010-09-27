#ifndef _STD_INCLUDE_H
#define _STD_INCLUDE_H

/* Bibliotheques standards */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#if defined (_ACML) || defined (_BLAS) || defined (_CBLAS)
#define _LOWLEVEL_MAT_LIB
#else
#define _NO_LOWLEVEL_MAT_LIB
#endif

#ifdef _ACML
#undef _BLAS
#undef _CBLAS
#undef _LAPACK
#endif

#ifdef _BLAS
/*#include <blas.h>*/
#endif
#ifdef _CBLAS
#include <cblas.h>
#endif
#ifdef _LAPACK
#include <clapack.h>
#endif
#ifdef _ACML
#include <acml.h>
#endif

#define COMPLEX double complex
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
/* For COMPLEX precision */
/* Using COMPLEX is not correct ('double COMPLEX' or 'float COMPLEX' must be used) */
#define REAL double

/* Structures */
struct Param_struct {
	int argc;
	char **argvcp;
	int pola;
	char calcul_type[SIZE_STR_BUFFER];
	char calcul_method[SIZE_STR_BUFFER];
	COMPLEX n_super;
	COMPLEX n_sub;
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
	COMPLEX *k2_layer;
	COMPLEX *invk2_layer;
	int N;
	int NS;
	int N_steps;
	int Ni;
	int ni;
	double Delta_sigma;
	COMPLEX sigma0;
	COMPLEX ky_0;
	COMPLEX k_super;
	COMPLEX k_sub;
	COMPLEX *sigma;
	COMPLEX *kz_super;
	COMPLEX *kz_sub;

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
	COMPLEX **n_xyz;
	int N_x;
	int N_z;
	int vec_size;
	int vec_middle;
	char profile_name[SIZE_STR_BUFFER];
		
	COMPLEX **tab_TF_k2;
	COMPLEX **tab_TF_invk2;
	
	COMPLEX *k2;
	COMPLEX *invk2;
	COMPLEX *Nx2;
	COMPLEX *Nz2;
	COMPLEX *NxNz;
	
	COMPLEX *TF_k2;
	COMPLEX *TF_invk2;
	COMPLEX *TF_Nx2;
	COMPLEX *TF_Nz2;
	COMPLEX *TF_NxNz;
	int HX_Normal_CALCULATED;

	COMPLEX **T;
	COMPLEX **P;
	COMPLEX **M;
	COMPLEX *eig_values;
	COMPLEX **EigVectors;
	COMPLEX *eig_buffer;
	
	COMPLEX **Toep_k2;
	COMPLEX **Toep_invk2;
	COMPLEX **invToep_k2;
	COMPLEX **invToep_invk2;
	COMPLEX **Toep_Nx2;
	COMPLEX **Toep_NxNz;
	COMPLEX **Toep_Nz2;
	COMPLEX **M_tmp1;
	COMPLEX **M_tmp2;
	COMPLEX **M_tmp3;
	
	COMPLEX *tmp_tf_k2;   
	COMPLEX *tmp_tf_invk2;
	COMPLEX *tmp_tf_Nx2;
	COMPLEX *tmp_tf_Nz2;
	COMPLEX *tmp_tf_NxNz;

	COMPLEX **Qxx;
	COMPLEX **Qyy;
	COMPLEX **Qzz;
	COMPLEX **Qxz;
	COMPLEX **Qzz_1;

	COMPLEX *QxzEx;
	COMPLEX *Qzz_1QxzEx;
	COMPLEX *Qzz_1Hpx;
	COMPLEX *V_tmp1;
	COMPLEX *sigmaHpy; 
	COMPLEX *ky0Qzz_1Hpx;
	COMPLEX *Qzz_1sigmaHpy;
	COMPLEX *QxxEx;
	COMPLEX *QyyEy;
	COMPLEX *QxzVtmp1;

	COMPLEX **S12;
	COMPLEX **S22;
	COMPLEX **S11;
	COMPLEX **S21;

	COMPLEX **T11;
	COMPLEX **T12;
	COMPLEX **T21;
	COMPLEX **T22;
	
	COMPLEX **M_tmp11;
	COMPLEX **M_tmp12;
	COMPLEX **M_tmp21;
	COMPLEX **M_tmp22;

/*	COMPLEX **Psi_sub_TE;
	COMPLEX **Psi_super_TE;
	COMPLEX **invPsi_super_TE;
	COMPLEX **Psi_sub_TM;
	COMPLEX **Psi_super_TM;
	COMPLEX **invPsi_super_TM;
	COMPLEX **invEigVec;
	COMPLEX **invVec_Psi;
	COMPLEX **M_invVec_Psi;
	COMPLEX **M_buffer_2vecsize;
	COMPLEX *M_sol;*/
	
	COMPLEX *Ai;
	COMPLEX *Ar;
	COMPLEX *At;
	
	COMPLEX *Vi;
	COMPLEX *Vr;
	COMPLEX *Vt;

	COMPLEX *Exi;
	COMPLEX *Exr;
	COMPLEX *Ext;
	COMPLEX *Hpxi;
	COMPLEX *Hpxr;
	COMPLEX *Hpxt;
	
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
/*	COMPLEX **Near_field_matrix;
	COMPLEX *bottom_field;*/
	COMPLEX ***tab_Z;
	COMPLEX ***tab_S12;
		
	/* mode extract S */
	int mode_extract_S;
	double h_extract_S;
	
	/* Pointeurs de fonctions */
	int (*k_2)(struct Param_struct *par, COMPLEX *k2_1D, double z);
	int (*invk_2)(struct Param_struct *par, COMPLEX *invk2_1D, double z);
	int (*Normal_function)(struct Param_struct *par, COMPLEX *Nx2, COMPLEX *NxNz, COMPLEX *Nz2, double z);
	int (*md2D_lire_profil)(const char *, struct Param_struct *);
	int (*M_matrix)(COMPLEX **M, double z, struct Param_struct *par);
	int (*P_matrix)(COMPLEX **P, double z, double Delta_z, struct Param_struct *par);
	
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
	char profile_file[SIZE_STR_BUFFER];
	char fichier_results[SIZE_STR_BUFFER];
	};

#endif /* _STD_INCLUDE_H */


