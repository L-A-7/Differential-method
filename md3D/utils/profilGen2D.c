/*!	\file	profilGenMULTI.c
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../md3D_utils.h"
#include "../md3D_io_utils.h"

#define MAX(a,b) ((a>b)?a:b)
#define SQUARE(a) ((a)*(a))
#define SIZE_STR 200
#define NPR_DEFAULT 512
#ifndef N_XY_ZINVAR
#define N_XY_ZINVAR 1
#endif
#ifndef PI
#define PI 3.14159265358979323846
#endif /* PI */


void message_erreur(){
	int Npr_default=NPR_DEFAULT;
	fprintf(stderr,	"  profilGen2D PROFILE [options] > file\n"
					"  PROFILE: profile name\n"
					"  Optionnal arguments :\n"
					"  -Nprx Nprx -Npry Npry ; profile number of points (%d by default)\n\n"
					,Npr_default);
}

int main(int argc, char *argv[]){
	
	int i, j, profile_type, nx, ny;
	char pr_name[SIZE_STR];
	double **pr;

	/* Copie des arguments de la ligne de commande */
	char **argvcp; 
	argvcp = (char **) malloc(sizeof(char*)*argc);
	argvcp[0] = (char*) malloc(sizeof(char)*SIZE_STR*argc);
	for(i=1;i<=argc-1;i++){
		argvcp[i] = argvcp[i-1] + SIZE_STR;
		strncpy(argvcp[i], argv[i],SIZE_STR);
	}
	
	/* Vérification de la présence du nombre minimal d'options */
	if (argc <= 1){
		message_erreur();
		return 1;
	}
	
	/* Valeurs par défaut */
	int Nprx = NPR_DEFAULT;
	int Npry = NPR_DEFAULT;

	lire_int_arg(&Nprx, "-Nprx", argc, argvcp);
	lire_int_arg(&Npry, "-Npry", argc, argvcp);
	strncpy(pr_name,argvcp[1],SIZE_STR);

	

	/*----- Profile generation -----*/
	if (!strcmp(pr_name,"CIRCLE")){
	
		profile_type= N_XY_ZINVAR;

		/* Memory allocation */
		pr = allocate_DbleMatrix(Npry, Nprx);

		double R = 0.25;
		double Lx = 1.0;
		double Ly = 1.0;
		lire_dble_arg(&R, "-R", argc, argvcp);
		lire_dble_arg(&Lx, "-Lx", argc, argvcp);
		lire_dble_arg(&Ly, "-Ly", argc, argvcp);
		if (R>MAX(Lx,Ly) || R<0) {
			fprintf(stderr,"CIRCLE: R must be between 0 and %f\n",MAX(Lx,Ly));
			free(pr);
			return 1;
		}

		for (ny=0;ny<=Npry-1;ny++){
			for (nx=0;nx<=Nprx-1;nx++){
/*				if (SQUARE((nx-(Nprx-1)/2)*Lx/Nprx)+SQUARE((ny-(Npry-1)/2)*Ly/Npry) <= SQUARE(R)){
*/				if (SQUARE(((double)nx-(double)(Nprx-1)/2)*Lx/(double)Nprx)+SQUARE(((double)ny-((double)Npry-1)/2)*Ly/(double)Npry) <= SQUARE(R)){
					pr[ny][nx]= 1;
				}else{
					pr[ny][nx]= 2;
				}
			}
		}
		
	}else{
		fprintf(stderr,"%s : Unkown profile type\n",pr_name);
		free(pr[0]);free(pr);
		return 1;
	}
	
	/* Enregistrement du profil */
	if (profile_type == N_XY_ZINVAR){
		fprintf(stdout,"profile_type = N_XY_ZINVAR\n");
		fprintf(stdout,"Nprx = %d\n",Nprx);
		fprintf(stdout,"Npry = %d\n",Npry);
		fprintf(stdout,"Nprz = %d\n",1);
		fprintf(stdout,"n_xyz = \n");
		for(i=0; i<=Npry-1; i++){
			for(j=0; j<=Nprx-1; j++){
				fprintf(stdout,"%1.0f ",pr[i][j]);
			}
			fprintf(stdout,"\n");
		}
	}else{
		fprintf(stderr,"%s : Unkown profile type\n",pr_name);
		return 1;
	}
	
	/* Libération de la mémoire */
	free(pr[0]); free(pr);
	
	return 0;
	
}
