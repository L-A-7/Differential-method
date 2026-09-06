/* tests.c */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include "md1D_utils.h"

complex **invMM(complex **inv, const complex **A, int );

int main()
{
	int N=100;
	int i,j;
	double sommInv=0,sommA=0,sommB=0;
	complex **A, **inv, **B;
	
	A = allocate_CplxMatrix(N,N);
	inv = allocate_CplxMatrix(N,N);
	B = allocate_CplxMatrix(N,N);

	
	/* Création d'une matrice */
	srand(time(NULL));
	for(i=0;i<=N-1;i++){
		for(j=0;j<=N-1;j++){
			A[i][j] = (complex) rand()/RAND_MAX;
		}
	}

	clock_t clock0=clock();
	for (i=0;i<1;i++) {
		invMM(inv, A, N);
	}
	printf("Temps écoulé : %f s\n",CHRONO(clock(),clock0));
	
	M_x_M(B,inv,A,N,N);

	for(i=0;i<=N-1;i++){
		for(j=0;j<=N-1;j++){
			sommA += cabs(A[i][j]);
			sommInv += cabs(inv[i][j]);
			sommB += cabs(B[i][j]);
		}
	}

	printf("inv : %f\nA : %f\nB : %f\n",sommInv/N,sommA/N,sommB/N);


	printf("Temps écoulé : %f s\n",CHRONO(clock(),clock0));
	return 0;
}




complex **invMM(complex **inv, const complex **A, int N)
{

	void partialPivoting(complex **M, complex **b, int k, int N);

	int i,j,k;
	complex **Id, **M, tmp;
	
	M = allocate_CplxMatrix(N, N);
	Id = allocate_CplxMatrix(N, N);

	/* Copie de la matrice */
	for(i=0;i<=N-1;i++){
		for(j=0;j<=N-1;j++){
			M[i][j] = A[i][j]; 
		}
	}
		
	/* Matrice Identité */
	for(i=0;i<=N-1;i++){
		for(j=0;j<=N-1;j++){
			Id[i][j] = (i==j);
		}
	}

	/* 'Triangularisation' de la matrice */
	for(k=0;k<=N-2;k++){
		/* Permutation des lignes pour que le pivot, M[k][k], soit le  */
		/* plus grand élément de la colonne. Minimise la propagation   */
		/* des erreurs d'arrondi, et évite les divisions par zero      */
		partialPivoting(M,Id,k,N);
		
		/* Eliminations des variables */
		for(i=k+1;i<=N-1;i++){
			tmp = M[i][k];
			for(j=k;j<=N-1;j++){
				M[i][j] -= M[k][j]*tmp/M[k][k];
			}
			for(j=0;j<=N-1;j++) {
				Id[i][j] -= Id[k][j]*tmp/M[k][k]; 
			}
		}
	}

	/* Résolution du système A.[inv] = [Id] après 'triangularisation' */		
	for(k=0;k<=N-1;k++) {
		inv[N-1][k] = Id[N-1][k]/M[N-1][N-1];
		for(i=N-2;i>=0;i--){
			inv[i][k] = Id[i][k]/M[i][i];
			for(j=i+1;j<=N-1;j++){
				inv[i][k] -= inv[j][k]*M[i][j]/M[i][i];
			}
		}
	}

		
	free(Id[0]); free(Id);
	free(M[0]); free(M);
	
	return inv;
}

void partialPivoting(complex **M, complex **b, int k, int N)
{

	int i, max=k;
	complex tmp;
	
	/* Recherche du plus grand élément */
	for(i=k;i<=N-1;i++){
		if (cabs(M[i][k]) > cabs(M[max][k])) max = i;
	}
	
	if(M[max][k]==0){
		fprintf (stderr, "%s : Erreur, impossible d'inverser la matrice", __FILE__); 
		exit (EXIT_FAILURE);
	}
	
	/* Permutation des lignes */
	for(i=k;i<=N-1;i++){
		tmp       = M[max][i];
		M[max][i] = M[k][i];
		M[k][i]   = tmp;
	}
	for(i=0;i<=N-1;i++){
		tmp       = b[max][i];
		b[max][i] = b[k][i];
		b[k][i]   = tmp;
	}

}
