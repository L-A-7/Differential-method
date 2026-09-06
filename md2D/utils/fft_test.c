/* tests.c */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>

#define CONST_N 30


int main()
{
	int i;
	complex *profil, *TF;
	fftw_plan plan;
	time_t t1;

	srand(time(&t1));

	profil = malloc(sizeof(complex)*CONST_N);
	TF = malloc(sizeof(complex)*CONST_N);

	for (i=0; i<=CONST_N-1; i++) {
		profil[i] = (complex) rand() / RAND_MAX;
		printf("%f\n",creal(profil[i]));
	}	
		

	plan = fftw_plan_dft_1d(CONST_N, profil, TF, FFTW_FORWARD, FFTW_ESTIMATE);

	fftw_execute(plan);
	
	printf("TF : \n");
	for (i=0; i<=CONST_N-1; i++) {
		printf("%f\n",creal(TF[i]));
	}
	
	
	fftw_destroy_plan(plan);	



	return 0;
}
