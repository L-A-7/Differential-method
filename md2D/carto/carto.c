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
#include "../md2D_utils.h"
#include "../md2D_io_utils.h"

#define STR_SIZE 5000

void err_message(){
	fprintf(stderr,	"usage : optim -fs scriptname -resfile resfilename -paramfile paramfilename -N_param N1 -N_datas N2 \n");
}

double fun(double *x, void *params);
int SaveDbleMatrix2file (double **M, int Nlign, int Ncol, char *filename);

int main(int argc, char *argv[]){
	
	void *params[10];
	int i, j, N_CD, N_hPoly, N_datas, N_param, read_verif;
	char scriptname[STR_SIZE], resfilename[STR_SIZE], paramfilename[STR_SIZE];
	double x[10], **carto, CD, hPoly;

	double CDmin=100;
	double CDMax=150;
	double CDstep=5;
	double hPolymin=150;
	double hPolyMax=200;
	double hPolystep=5;

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
	
				
	params[0] = (void *) &N_param;
	params[1] = (void *) &N_datas;
	params[2] = (void *) &scriptname[0];
	params[3] = (void *) &resfilename[0];


	N_CD = ROUND((CDMax-CDmin)/CDstep)+1;
	N_hPoly = ROUND((hPolyMax-hPolymin)/hPolystep)+1;
	carto = allocate_DbleMatrix(N_CD,N_hPoly);

	for(i=0;i<=N_CD-1;i++){
		for(j=0;j<=N_hPoly-1;j++){
			carto[i][j] = 0;
		}
	}
	for(i=0;i<=N_CD-1;i++){
		for(j=0;j<=N_hPoly-1;j++){
			CD = CDmin + i*CDstep;
			hPoly = hPolymin +j*hPolystep;
			x[0] = CD;
			x[1] = hPoly;
			carto[i][j] = fun(x, params);
			SaveDbleMatrix2file (carto, N_CD, N_hPoly, "cartographie.txt");
		}
	}

	free(carto[0]);
	free(carto);
		
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

int SaveDbleMatrix2file (double **M, int Nlign, int Ncol, char *filename)
{
  int i,j;
  FILE *fp;
	
  /* Affichage à l'écran si filename = "stdout" */										
  if(!strcmp(filename,"stdout")) fp = stdout;
  else fp = fopen(filename, "w");
	
  /* Mode = "Re"/"Im" : enregistrement de la partie réelle/Imaginaire */
  for (i=0;i<=Nlign-1;i++){
    for (j=0;j<=Ncol-1;j++){
		fprintf(fp,"% 1.6e  ",M[i][j]);
    }
    fprintf(fp,"\n");
  }

  if(strcmp(filename,"stdout")) fclose(fp);

  return 0;


}
