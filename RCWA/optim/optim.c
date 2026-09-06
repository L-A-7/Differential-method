/*
 *	optim.c : routine for optimisation
 *
 *	usage : optim -funscript name_of_script.sh -resultsfile results_filename.txt -parametersfile param_filename.txt
 *						-N_exp_points Number_of_data_points
 *	funscript : script calling the function to minimize
 * resultsfile : file where results from funscript are stored
 * parametersfile
 *
 *	Uses the Levenberg-Marquardt opimisation algorithm
 *
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "lm.h"
#include "../md2D_utils.h"
#include "../md2D_io_utils.h"

#define STR_SIZE 5000

void err_message(){
	fprintf(stderr,	"usage : optim -fs scriptname -resfile resfilename -paramfile paramfilename -N_param N1 -N_datas N2 \n");
}

void func(double *parameters, double *simu_datas, int N_param, int N_data, void *adatas);

int main(int argc, char *argv[]){
	
	char *adatas[10];
	int i, N_datas, N_param, read_verif, itmax = 1000;
	char scriptname[STR_SIZE], resfilename[STR_SIZE], paramfilename[STR_SIZE];
	double *parameters, *null_datas, info[9], opts[5];
	
	opts[0] = 1e-3; /* mu   */
	opts[1] = 1e-5; /* eps1 */
	opts[2] = 1e-5; /* eps2 */
	opts[3] = 1e-5; /* eps3 */
	opts[4] = 1e-5; /* delta, for derivative evaluation */
	
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
	
	/* Mem alloc */
	parameters = (double *) malloc(sizeof(double)*N_param);
	null_datas = (double *) malloc(sizeof(double)*N_datas);
	
	/* Creating fake experimental data (necessary for the routine) */
	for (i=0;i<=N_datas-1;i++){
		null_datas[i] = 0;
	}
	
	/* Reading initial parameters in paramfile */
	lire_tab(paramfilename, "", parameters, N_param);	
	
	/* Calling the optimization algorithm */	
	adatas[0] = scriptname;
	adatas[1] = resfilename;
	dlevmar_dif(func, parameters, null_datas, N_param, N_datas, itmax, opts, info, NULL, NULL, (void *) adatas);

	/* Returning optimum parameters in paramfile */
	SaveDbleTab2file (parameters, N_param, paramfilename, " ");

	/* freeing memory */
	free(parameters);
	free(null_datas);

	return 0;
}


void func(double *parameters, double *simu_datas, int N_param, int N_data, void *adatas)
{
	int i;
	char** tabdata = (char **) adatas;
	char *scriptname = tabdata[0];
	char *resfilename = tabdata[1];
	char command[STR_SIZE];

	fprintf(stdout,"%f %f ",parameters[0],parameters[1]);fflush(stdout);
	
	/* Calling the script */
	sprintf(command,"%s",scriptname);
	for (i=0;i<=N_param-1;i++){
		sprintf(command,"%s %f",command, parameters[i]);
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
	
	fprintf(stdout,"%f\n",Sq/N_data);fflush(stdout);

	
}
