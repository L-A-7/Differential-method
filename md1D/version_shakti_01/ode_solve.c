/*!	\file	ode_solve.c
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void rk4_step(double *y, double *dy, int N, double t, double h, double *yout,
	int (*f)(double, const double *, double *, void *),void *param_void);


/*!	\fn		int ode_solve(const double *y0, double *y, int N, double t0, double t1, int nstep,
	void (*f)(double, const double *, double *, void *), void *param_void)
 *
 *	\brief	Routine de résolution d'un systeme d'équation différentielles réelles du 1er ordre \n
 *			défini par  [dy(t)/dt] = f([y(t)]) où  [dy(t)/dt] et [y(t)] sont des vecteurs de taille \n
 *			N
 *
 *
 *	\param	y0    : y(t0), conditions initiales
 *  \param	y     : y(t1)
 *	\param	N     : taille du systeme
 *	\param	nstep : nombre de pas d'intégration
 *	\param	f     : pointeur vers une fonction calculant [dy(t)/dt] à partir de [y(t)]  \n
 *			        de la forme : void f (double t, const double *y, double *dy, void *param)
 */
int ode_solve(const double *y0, double *y, int N, double t0, double t1, int nstep,
	int (*f)(double, const double *, double *, void *), void *param_void)
{

	int i,n;
	double t,h, *dy, *yout;

	dy = (double *) malloc(sizeof(double)*N);
	yout = (double *) malloc(sizeof(double)*N);

	h=(t1-t0)/nstep;

	for (i=0; i<=N-1; i++) {
		y[i] = y0[i];
	}

	for(n=0; n<=nstep-1; n++) {

		t = t0 + n*h; /* ICI ou à la fin de la boucle ? */
		f(t, y, dy, param_void);
		rk4_step(y, dy, N, t, h, yout, f, param_void);
		for (i=0; i<=N-1; i++) {
			y[i] = yout[i];
		}
	}

	free(dy);
	free(yout);

	return 0;
}


void rk4_step(double *y, double *dy, int N, double t, double h, double *yout,
	int (*f)(double, const double *, double *, void *), void *param_void)
{
	int i;
	double th,hh,h6,*dym,*dyt,*yt;

	dym = (double *) malloc(sizeof(double)*N);
	dyt = (double *) malloc(sizeof(double)*N);
	yt  = (double *) malloc(sizeof(double)*N);

	hh = h*0.5;
	h6 = h/6.0;

	th = t+hh;
	for (i=0;i<=N-1;i++) yt[i]=y[i]+hh*dy[i];
	(*f)(th,yt,dyt,param_void);
	for (i=0;i<=N-1;i++) yt[i]=y[i]+hh*dyt[i];
	(*f)(th,yt,dym,param_void);
	for (i=0;i<=N-1;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*f)(t+h,yt,dyt,param_void);
	for (i=0;i<=N-1;i++)
		yout[i]=y[i]+h6*(dy[i]+dyt[i]+2.0*dym[i]);

	free(yt);
	free(dyt);
	free(dym);

}




