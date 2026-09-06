/*
 *	sub_tanpsicosdelta.c
 *
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../md2D_utils.h"
#include "../md2D_io_utils.h"

#define STR_SIZE 5000

void err_message(){
	fprintf(stderr,	"usage : sub_tanpsicosdelta Nb_lambda SIMU_FILE_NAME EXP_FILE_NAME\n");
}

int main(int argc, char *argv[]){
	
	int i,nl,Nlambda;
	double eps, *lambda_simu, *lambda_exp, *tanPsi_simu, *tanPsi_exp, *cosDelta_simu, *cosDelta_exp, *diff;
	char simu_file[STR_SIZE], exp_file[STR_SIZE], *endptr;
	
	/* Copie des arguments de la ligne de commande */
	char **argvcp; 
	argvcp = (char **) malloc(sizeof(char*)*argc);
	argvcp[0] = (char*) malloc(sizeof(char)*STR_SIZE*argc);
	for(i=1;i<=argc-1;i++){
		argvcp[i] = argvcp[i-1] + STR_SIZE;
		strncpy(argvcp[i], argv[i],STR_SIZE);
	}
	
	/* Vérification de la présence du nombre minimal d'options */
	if (argc < 2){
		err_message();
		return 1;
	}
	
	/* Reading the arguments */
	Nlambda = (int) strtod(argvcp[1], &endptr);
	if (argvcp[1] == endptr){
		fprintf(stderr,"ERROR, %s can't read argument Nlambda. Exiting.",__FILE__);
		exit(EXIT_FAILURE);
	}
	strncpy(simu_file, argvcp[2],STR_SIZE);
	strncpy(exp_file, argvcp[3],STR_SIZE);
	
	/* Memory allocation */
	lambda_simu = (double *) malloc(sizeof(double)*Nlambda);
	lambda_exp = (double *) malloc(sizeof(double)*Nlambda);
	tanPsi_simu = (double *) malloc(sizeof(double)*Nlambda);
	tanPsi_exp = (double *) malloc(sizeof(double)*Nlambda);
	cosDelta_simu = (double *) malloc(sizeof(double)*Nlambda);
	cosDelta_exp = (double *) malloc(sizeof(double)*Nlambda);
	diff = (double *) malloc(sizeof(double)*2*Nlambda);

	/* Reading the data */
	lire_tab(simu_file, "lambda", lambda_simu, Nlambda);	
	lire_tab( exp_file, "lambda", lambda_exp , Nlambda);	
	lire_tab(simu_file, "tan_Psi", tanPsi_simu, Nlambda);	
	lire_tab( exp_file, "tan_Psi", tanPsi_exp , Nlambda);	
	lire_tab(simu_file, "cos_Delta", cosDelta_simu, Nlambda);	
	lire_tab( exp_file, "cos_Delta", cosDelta_exp , Nlambda);	

	/* Checking that the exp & simu values of lambda are the same */
	eps = 1e-5;
	for (nl=0;nl<=Nlambda-1;nl++){
		if (fabs(lambda_simu[nl]-lambda_exp[nl]) > eps){
			fprintf(stderr,"ERROR, %s, lambda_exp and lambda_simu don't have the same values\n",__FILE__);
			exit(EXIT_FAILURE);
		}
	}
	
	/* Calculating the difference between tanPsi + cosDelta exp and simulated */
	for (nl=0;nl<=Nlambda-1;nl++){
		diff[nl] = tanPsi_exp[nl] - tanPsi_simu[nl];
		diff[nl+Nlambda] = cosDelta_exp[nl] - cosDelta_simu[nl];
	}
	SaveDbleTab2file (diff, 2*Nlambda, "stdout", " ");
	
	free(lambda_simu);
	free(lambda_exp);
	free(tanPsi_simu);
	free(tanPsi_exp);
	free(cosDelta_simu);
	free(cosDelta_exp);
	free(diff);
	
	return 0;
}
