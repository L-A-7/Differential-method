/* test2.c */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include "../md1D_utils.h"
#include "../std_include.h"


int main()
{

	int i,j,a;
	time_t time0, time1;
	
	time(&time0);
	for (i=0;i<7;i++) {
		for (j=0;j<100000000;j++) {
			a=1+2*i+3*i*i+4*i*i*i;
		}
		time(&time1);
		printf("difftime = %f \n",(float)difftime(time1,time0));
	}
	
	
	return 0;
}




int md1D_temps(int n, int N, int nS, int NS, struct Param_struct *par){

	/* Si moins de 5 secondes depuis le dernier affichage, on ne change rien */
	if (CHRONO(clock(), par->last_clock) < .5){
		return 0;
	/* Sinon, estimation et affichage de la durée restante */
	}else{
	int i;
		par->last_clock = clock();
		float t_ecoule = CHRONO(clock(),par->clock0);
		float t_total = t_ecoule*(NS*(2*N+1))/((NS-nS)*(2*N+1)+n+N+1);
		float t_restant = t_total - t_ecoule;
		int pourcent = ROUND(100.0*t_ecoule/t_total);
				
		fprintf(stdout,"\r");
		fprintf(stdout,"%3d %% |", pourcent);
		for (i=0;i<pourcent/5;i++) {fprintf(stdout,"|");}
	/*	fprintf(stdout,">");
	*/	for (i=pourcent/5;i<20;i++) {fprintf(stdout," ");}
		fprintf(stdout,"| reste %d s sur %d s",ROUND(t_restant), ROUND(t_total));
		fflush(stdout);
	}
	
	return 0;
}
