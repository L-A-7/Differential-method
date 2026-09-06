/* tests.c */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>

#define ALEAT() (ceil(20*(double)random()/RAND_MAX)/10)
#define ROUND_LONG(x) ((long int)(x<0 ? x-0.5 : x+0.5))
#define FLOOR(x) ((int)(x<0 ? x : x))	
#define DBLE_CMP_EXIGEANCE 1000000
#define CHRONO(t2,t1) ((double)(t2-t1)/CLOCKS_PER_SEC)
#define CONST_N 10


static int Taille;

long int *cherche(long int z, long int *tab, int N);

int main()
{
	int i, j, k,res;
	double z, h=1 ;
	long int tab[2][CONST_N], z_int;
	clock_t clock0=clock();
	time_t t1;
	long int *adr;	
	
	Taille = 0;

srandom(time(&t1));

	for (i=0; i<=CONST_N-1; i++) {
		
		/* Affichage de la table */
		for(k=0; k<=Taille-1; k++){
			printf("tab[0][%d] = %d, tab[1][%d] = %d\n",k, tab[0][k],k,tab[1][k]);
		}
				
		/* Génération de z */
		z = ALEAT();
		z_int = ROUND_LONG((z/h)*DBLE_CMP_EXIGEANCE);
		printf("z = %d\n", z_int);
			
		/* Recherche de z */
		adr = cherche(z_int, tab[0], Taille);
		if (adr != NULL) {
			/* z trouvé */
			printf("z trouvé en position %d\n",adr-tab[0]);
		}else{
			/* Ajout de z */
			Taille++; /* AJOUTER vérifications et allocation eventuelle */
			k = Taille-2;
			while(z_int > tab[0][k] && k>=0){
				tab[0][k+1] = tab[0][k];
				tab[1][k+1] = tab[1][k];
				k--;
			}
			tab[0][k+1] = z_int;
			tab[1][k+1] = Taille-1;
		}
		
	}
double x = -3.2;
printf("FLOOR(%f) = %d\n",x,FLOOR(x));
/*	z = 1.6;
	z_int = ROUND_LONG((z/h)*DBLE_CMP_EXIGEANCE);
	for (i=0; i<=1000000; i++){
		res = cherche(z_int, tab[0], CONST_N);
	}
	printf("z_int  = %d\n",z_int);
	printf("res = %d\n",res);
	if (res) printf("z trouvé en position %d\n",adr-tab[0]);
*/

	printf("Temps écoulé : %f s\n",CHRONO(clock(),clock0));

	printf("time : %d s\n",time(&t1));


	return 0;
}



/*!	\fn		long int *cherche(long int z, long int *tab0, int N)
 *
 *	\brief	Recherche dichotomique de l'élément z dans un tableau tab0 de taille N \n
 *			classé dans l'ORDRE DECROISSANT. 
 *
 * 	\return	L'adresse correspondant à l'élément trouvé ou NULL si l'élément n'est pas présent
 */
long int *cherche(long int z, long int *tab0, int N)
{

	if (N<=1) {
		if (z==tab0[0]) return tab0;
		else            return NULL;
	}
	
	if (tab0[N>>1] < z) return(cherche(z, tab0, N>>1));
	else                return(cherche(z, tab0+(N>>1), N-(N>>1)));
}



