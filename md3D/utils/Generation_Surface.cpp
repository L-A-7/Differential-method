#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <conio.h>

// Nombre de points :
#define N	1000
// Taille du filtre :
#define TAILLE_FILTRE	300

int main()
{	   
	FILE *Fichier;

	// GENERATION DE N NOMBRES ALEATOIRES COMPRIS ENTRE -0.05 ET 0.05
	double RandomValue[N + 2 * TAILLE_FILTRE + 1];

	Fichier = fopen("FichierNbresAleatoires.xls", "w+");
	//srand(time(0));
	
	// Initialisation de RandomValue à 0
	for (int i = - TAILLE_FILTRE; i < N + TAILLE_FILTRE; i++)
		RandomValue[i] = 0;

	for (i = 0; i < N + TAILLE_FILTRE/2; i++)
	{
		RandomValue[i] = ((rand() % 10000) / 10000.0) - 0.5;
		fprintf(Fichier, "%d\t%lf\n", i, RandomValue[i]);
	}

	fclose(Fichier);

	// Ag = Amplitude de la gaussienne, Lg = Largeur de la gaussienne
	// Ae = Amplitude de l'exponentielle, Le = Largeur de l'exponentielle
	double Ag = 1000.0;
	double Lg = 1;
	double Ae = 0.0;
	double Le = 0;

	// CONVOLUTION DES NOMBRES ALEATOIRES PAR UNE GAUSSIENNE + UNE EXPONENTIELLE
	double Result[N + 2 * TAILLE_FILTRE + 1];
	// Initialisation de Result à 0
	for (i = - TAILLE_FILTRE; i < N + TAILLE_FILTRE; i++)
		Result[i] = 0;
		  
	// filtre = somme de la gausienne et de l'exponentielle
	Fichier = fopen("FichierFiltre.xls", "w+");
	double filtre[TAILLE_FILTRE + 1];
	for (i = 0; i < 1 + TAILLE_FILTRE; i++)
		filtre[i] = 0.0;
	for (i = -TAILLE_FILTRE / 2; i < TAILLE_FILTRE / 2; i++)
	{
		filtre[i + TAILLE_FILTRE/2] = Ag * exp(-((i/Lg)*(i/Lg))) + Ae * exp(-(abs(i)/Le));
	}
	filtre[TAILLE_FILTRE/2] = Ag;

	for (i = 0; i < TAILLE_FILTRE; i++)
	{
		fprintf(Fichier, "%lf\n", filtre[i]);
	}

	fclose(Fichier);

	// Convolution
	Fichier = fopen("FichierSurface.xls", "w+");
	RandomValue[-4] = 0;
	for (int m = 0; m < N + TAILLE_FILTRE/2; m++)
	{
		for (i = 0; i <= TAILLE_FILTRE; i++)
		{
			Result[m] += (RandomValue[m-i] * filtre[i]);
		}
	}
	for (m = 0; m < N; m++)
		fprintf(Fichier, "%lf\n", Result[m + TAILLE_FILTRE/2]);

	fclose(Fichier);

	return 0;
}