#ifndef _STD_INCLUDE_H
#define _STD_INCLUDE_H

/* Bibliotheques standards */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#define COMPLEX double complex

#ifdef _ACML
#undef _BLAS
#undef _LAPACK
#endif

#ifdef _BLAS
#include <cblas.h>
#endif
#ifdef _LAPACK
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

#define H_XY          1
#define MULTICOUCHES 2
#define N_XYZ        3
#define N_XY_ZINVAR  4
#define H_XY_plus_STACK 5
#define SIZE_STR_BUFFER 200
#define SIZE_LINE_BUFFER 50000
#define SIZE_INT_BUFFER 1000

#define NON_LU "Et_non_c_pas_lu"
/* Macros */
/*#define c_omplex(a,b) (a+I*b) */
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
/*#define REAL double*/

/* Structures */
struct Param_struct {
	int argc;
	char **argvcp;
	char calcul_type[SIZE_STR_BUFFER];
	char calcul_method[SIZE_STR_BUFFER];
	COMPLEX nu_super;
	COMPLEX nu_sub;
	double Lx;
	double Ly;
	double h;
	double coef_h;
	double lambda;
	double theta_i;
	double phi_i;
	double psi;
	int profile_type;
	COMPLEX **n_xyz;
	int N_layers;
	COMPLEX *k2_layer;
	COMPLEX *invk2_layer;
	int Nx;
	int Ny;
	
	int NS;
	int N_steps;
	int *nx;
	int *ny;
	double Delta_sigma_x;
	double Delta_sigma_y;
	COMPLEX sigma_x0;
	COMPLEX sigma_y0;
	COMPLEX k_super;
	COMPLEX k_sub;
	COMPLEX *sigma_x;
	COMPLEX *sigma_y;
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
	int Nprx;
	int Npry;
	int Nprz;
	int vec_size;
	int vec_middle;
	char profile_name[SIZE_STR_BUFFER];

	int HXY_Normal_CALCULATED;
	int toepNorm_CALCULATED;
	COMPLEX **P;
	COMPLEX **M;

	COMPLEX **Nxx;
	COMPLEX **Nxy;
	COMPLEX **Nxz;
	COMPLEX **Nyy;
	COMPLEX **Nyz;
	COMPLEX **Nzz;

	COMPLEX **Qxx;
	COMPLEX **Qxy;
	COMPLEX **Qxz;
	COMPLEX **Qyy;
	COMPLEX **Qyz;
	COMPLEX **Qzz;
	COMPLEX **Qzz_1;

	COMPLEX **S12;
	COMPLEX **S22;
	COMPLEX **S11;
	COMPLEX **S21;

	COMPLEX **T11;
	COMPLEX **T12;
	COMPLEX **T21;
	COMPLEX **T22;

	COMPLEX **Psi_sub;
	COMPLEX **Psi_super;
	COMPLEX **invPsi_super;

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
	
	
	/* Pointeurs de fonctions */
	COMPLEX* (*k_2)(struct Param_struct *par, COMPLEX *k2_1D, double z);
	COMPLEX* (*invk_2)(struct Param_struct *par, COMPLEX *invk2_1D, double z);
	int (*Normal_function)(COMPLEX **norm_x, COMPLEX **norm_y, COMPLEX **norm_z, double z, double *profil, struct Param_struct *par);
	int (*md3D_lire_profil)(const char *, struct Param_struct *);
	int (*M_matrix)(COMPLEX **M, double z, struct Param_struct *par);
	int (*P_matrix)(COMPLEX **P, double z, double Delta_z, struct Param_struct *par);
	int (*Normal_to_profile)(COMPLEX **norm_x, COMPLEX **norm_y, COMPLEX **norm_z, double z, struct Param_struct *par);

	clock_t clock0;         /* Stocke le temps de départ                              */
	clock_t last_clock;     /* Durée écoulée depuis le dernier appel à md3D_temps     */
	time_t  time0;			/* clock(): très précis mais cyclique => tps longs pas OK */
	time_t last_time;		/* time() : précision d'1 s, non cyclique => tps longs OK */


	int verbosity; /* Si verbose = 1, affiche plus d'infos sur le terminal */
	char i_field_mode[SIZE_STR_BUFFER]; /* incident field (PLANE_WAVE, GAUSSIAN, FROM_BINARY, FROM_ASCII) */

	/* calcul_type = STACK */
	int N_stack;
	int n_patterned_layer;
	double *h_stack;
	COMPLEX *nu_stack;
	COMPLEX nu_this_layer;

	char profile_file[SIZE_STR_BUFFER]; /* is present both in this structure and in fichiers structures ... */
	};

struct Efficacites_struct {
	int Nxmin_super;
	int Nxmax_super;
	int Nxmin_sub;
	int Nxmax_sub;
	int Nymin_super;
	int Nymax_super;
	int Nymin_sub;
	int Nymax_sub;

	double **eff_r;
	double **eff_t;
	
	double **nx_eff_r;	
	double **ny_eff_r;	
	double **nx_eff_t;
	double **ny_eff_t;
	
	double sum_eff_r;
	double sum_eff_t;
	double sum_eff;
	};

struct Noms_fichiers {
	char fichier_config[SIZE_STR_BUFFER];
	char param_file[SIZE_STR_BUFFER];
	char profile_file[SIZE_STR_BUFFER];
	char fichier_results[SIZE_STR_BUFFER];
	};

#endif /* _STD_INCLUDE_H */


