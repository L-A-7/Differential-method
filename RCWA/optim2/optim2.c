/*
 *	optim.c : routine for optimisation
 *
 *	usage : optim -funscript name_of_script.sh -resultsfile results_filename.txt -parametersfile param_filename.txt
 *						-N_exp_points Number_of_data_points
 *	funscript : script calling the function to minimize
 * resultsfile : file where results from funscript are stored
 * parametersfile
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include "../md2D_utils.h"
#include "../md2D_io_utils.h"

#define STR_SIZE 5000

void err_message(){
	fprintf(stderr,	"usage : optim -fs scriptname -resfile resfilename -paramfile paramfilename -N_param N1 -N_datas N2 \n");
}

double func(double *x, void *params);
double my_f(const gsl_vector *x, void *param);
void my_df(const gsl_vector *x, void *param, gsl_vector *g);
void my_fdf(const gsl_vector *x, void *param, double *f, gsl_vector *g);
double fundf(double *x, double *df, void *params);

int main(int argc, char *argv[]){
	
	void *params[10];
	int i, N_datas, N_param, read_verif, status, iter = 0, itmax = 10;
	char scriptname[STR_SIZE], resfilename[STR_SIZE], paramfilename[STR_SIZE];
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
	
	/* Copie des arguments de la ligne de commande */
	char **argvcp; 
	argvcp = (char **) malloc(sizeof(char*)*argc);
	argvcp[0] = (char*) malloc(sizeof(char)*STR_SIZE*argc);
	for(i=1;i<=argc-1;i++){
		argvcp[i] = argvcp[i-1] + STR_SIZE;
		strncpy(argvcp[i], argv[i],STR_SIZE);
	}
	
	/* Vérification de la présence du nombre minimal d'options */
	if (argc < 11){
		err_message();
		return 1;
	}
	
	/* Reading the arguments */
	read_verif = 0;
	read_verif += lire_int_arg(&N_datas, "-N_datas", argc, argvcp);
	read_verif += lire_int_arg(&N_param, "-N_param", argc, argvcp);
	read_verif += lire_str_arg(scriptname, "-fs", argc, argvcp);
	read_verif += lire_str_arg(resfilename, "-resfile", argc, argvcp);
	read_verif += lire_str_arg(paramfilename, "-paramfile", argc, argvcp);
	if (read_verif != 0){
		err_message();
		fprintf(stderr,"ERROR, %s can't read all his arguments, exiting\n",__FILE__);
		exit(EXIT_FAILURE);
	}
	
	
	gsl_multimin_function_fdf my_func;
	my_func.f = &my_f;
	my_func.df = &my_df;
	my_func.fdf = &my_fdf;
	my_func.n = N_param;
	my_func.params = params;
			
	params[0] = (void *) &N_param;
	params[1] = (void *) &N_datas;
	params[2] = (void *) &scriptname[0];
	params[3] = (void *) &resfilename[0];
	
	/* Reading initial x in paramfile */
	gsl_vector *x;
	x = gsl_vector_alloc(N_param);
	lire_tab(paramfilename, "", x->data, N_param);	

/*	T = gsl_multimin_fdfminimizer_conjugate_fr;*/
	T = gsl_multimin_fdfminimizer_vector_bfgs2;
	s = gsl_multimin_fdfminimizer_alloc(T,N_param);
	
	gsl_multimin_fdfminimizer_set(s, &my_func, x, 1, 0.5);

	do{
		iter++;
printf("gsl_iterate...\n");fflush(stdout);
		status = gsl_multimin_fdfminimizer_iterate(s);
printf("done.\n");fflush(stdout);
		
		if (status){
			break;
		}
printf("gsl_check_gradient...\n");fflush(stdout);
		status = gsl_multimin_test_gradient(s->gradient,5e-5);
printf("done.\n");fflush(stdout);
		if (status == GSL_SUCCESS){
			printf("Minimum found at:\n");
		}
		
	}while (status == GSL_CONTINUE && iter < itmax);
	
	gsl_multimin_fdfminimizer_free(s);
	gsl_vector_free(x);
	
	return 0;
}


double fun(double *x, void *params)
{
	int i;
	double Sq, MSqDiff;
	double *simu_datas;
	char command[STR_SIZE];
	
	void **par = (void **) params;
	int N_param = *((int *)par[0]);
	int N_data  = *((int *)par[1]);
	char *scriptname = (char *)par[2];
	char *resfilename = (char *)par[3];
		
	simu_datas = (double *) malloc (sizeof(double)*N_data);

	for(i=0;i<=N_param-1;i++){
		fprintf(stdout,"%f ",x[i]);fflush(stdout);
	}
	/* Calling the script */
	sprintf(command,"%s",scriptname);
	for (i=0;i<=N_param-1;i++){
		sprintf(command,"%s %f",command, x[i]);
	}
	if (system(command) != 0){
		fprintf(stderr,"ERROR, %s, failed to execute script \"%s\". Exiting\n",__FILE__,command);
		exit(EXIT_FAILURE);
	}
	/* Reading the results */
	if (lire_tab(resfilename, "", simu_datas, N_data) !=0 ){
		fprintf(stderr,"ERROR, %s, failed to read th results in %s. Exiting\n",__FILE__,resfilename);
		exit(EXIT_FAILURE);
	}	

	/* Calculating results mean square differences */
	for (i=0,Sq=0.0;i<=N_data-1;i++){
		Sq += simu_datas[i]*simu_datas[i];
	}
	MSqDiff = Sq/N_data;

	fprintf(stdout,"%f\n",MSqDiff);fflush(stdout);

	free(simu_datas);
	
	return MSqDiff;
}

double fundf(double *x, double *df, void *params)
{
	int i;
	double f, dx, tmp;
	
	double period = 310;
	double Nprofil = 2048*4;
	
	void **par = (void **) params;
	int N_param = *((int *)par[0]);

	f = fun(x, params);
	for (i=0;i<=N_param-1;i++){
		if (i==0){
			dx=2.0*period/Nprofil;
		}else{
			dx=0.01;
		}
		tmp = x[i];
		x[i] += dx;
printf("    ");fflush(stdout);
		df[i] = (fun(x, params) - f)/dx;
		x[i] = tmp;
	}
	
	return f;
}



double my_f(const gsl_vector *x, void *param)
{
	double *x_dble;
	
	x_dble = x->data;
printf("f   ");fflush(stdout);

	return fun(x_dble, param);
}

void my_fdf(const gsl_vector *x, void *param, double *f, gsl_vector *g)
{
	double *x_dble, *g_dble;
	
	x_dble = x->data;
	g_dble = g->data;
printf("fdf ");fflush(stdout);

	*f = fundf(x_dble, g_dble, param);
}

void my_df(const gsl_vector *x, void *param, gsl_vector *g)
{
	double *x_dble, *g_dble;
	
	x_dble = x->data;
	g_dble = g->data;
printf("df  ");fflush(stdout);

	fundf(x_dble, g_dble, param);
}

