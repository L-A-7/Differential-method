/*!	\file	eq_diff.c
 * 
 *	\brief	Résolution du système d'équations différentielles (utilise la biliothèque GSL)
 */
 
#include "std_include.h"
#include <gsl/gsl_odeiv.h>


/*! \fn		int eq_diff(const double *y0, double *y, int N, double t0, double t1, int nstep,
	int (*func)(double, const double *, double *, void *), void *param_void)

 */
int eq_diff(const double *y0, double *y, int N, double t0, double t1, int nstep,
	int (*func)(double, const double *, double *, void *), void *param_void)
{

	int i;
	double *y_err;
	
	y_err = (double *) malloc(sizeof(double)*N);

	/* Copie des conditions initiales */
	for(i=0;i<=N-1;i++){
		y[i] = y0[i];
	}

	/* Type de methode d'intégration */
	const gsl_odeiv_step_type * step_type = gsl_odeiv_step_rk4;
/*	const gsl_odeiv_step_type * step_type = gsl_odeiv_step_rk8pd;*/

	/* Pointeur sur la methode d'intégration */
	gsl_odeiv_step * step_method = gsl_odeiv_step_alloc (step_type, N);

	/* Définition du système à intégrer */
	gsl_odeiv_system sys = {func, NULL, N, param_void};

	double step = (t1-t0)/nstep;
	double t = t0;

	for(i=0;i<=nstep-1;i++){
		gsl_odeiv_step_apply (step_method, t, step, y, y_err, NULL, NULL, &sys);
		t += step;
	}


	gsl_odeiv_step_free (step_method);
	free(y_err);

	return 0;
}



