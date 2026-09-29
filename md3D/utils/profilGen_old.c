/*!	\file	profilGen.c
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define SIZE_STR_BUFFER 200
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif /* M_PI */


void message_erreur(){
	fprintf(stderr,	"  profilGen -f file -n pr_name\n"
					"  file    : fichier destination\n"
					"  pr_name : nom du profil\n"
					"  Arguments optionnels :\n"
					"  -N N_profil ; nombre de points du profil, par defaut 256\n\n"
					);
}

int main(int argc, char *argv[]){
	
	int i, cpt_options = 0;
	const int MIN_OPTIONS = 2;
	char filename[SIZE_STR_BUFFER], pr_name[SIZE_STR_BUFFER];
	double *pr, min, max;
	FILE *fp;
	
	/* Valeurs par défaut */
	int N_profil = 256;
	
	/* Lecture des options de la ligne de commande */
	for(i=1;i<=argc-1;i++){
		if (!strcmp(argv[i],"help")){
			message_erreur();
			return 0; 
		}
		if (!strcmp(argv[i],"-f") && (i+1<=argc-1)){
			strncpy(filename, argv[i+1],SIZE_STR_BUFFER);
			cpt_options++; 
		}
		if (!strcmp(argv[i],"-n") && (i+1<=argc-1)){
			strncpy(pr_name, argv[i+1],SIZE_STR_BUFFER);
			cpt_options++; 
		}
		if (!strcmp(argv[i],"-N") && (i+1<=argc-1)){
			N_profil = (int) strtol(argv[i+1], (char **)NULL, 10);/* Ajouter gestion d'erreur */
		}
	}
	
	/* Vérification de la présence du nombre minimal d'options */
	if (cpt_options<MIN_OPTIONS){
		message_erreur();
		return 1;
	}
	
	/* Allocation de mémoire pour le profil*/
	pr = (double *) malloc(N_profil*sizeof(double));
	
	/* Génération d'un profil */
	if (!strcmp(pr_name,"COS")){
		for(i=0;i<=N_profil-1;i++){
			pr[i] = cos(2*M_PI*i/N_profil);
		}
	}else if(!strcmp(pr_name,"SIN")){
		for(i=0;i<=N_profil-1;i++){
			pr[i] = sin(2*M_PI*i/N_profil);
		}
	}else if(!strcmp(pr_name,"LINEAIRE")){
		for(i=0;i<=N_profil-1;i++){
			pr[i] = i;
		}
	}else if(!strcmp(pr_name,"ALEAT")){
		srand(time(NULL));
		for(i=0;i<=N_profil-1;i++){
			pr[i] = (double) rand();
		}
	}else if(!strcmp(pr_name,"CARRE")){
		srand(time(NULL));
		double a, n_period = 3.0;
		for(i=0;i<=N_profil-1;i++){
			a = (double)i/N_profil*n_period;
			pr[i] = (double)(int)(a-(int)a+0.5);
		}
	}else if(!strcmp(pr_name,"FLAT")){
		for(i=0;i<=N_profil-1;i++){
			pr[i] = 0.5;
		}
		goto ENREGISTREMENT;
	}else{
		fprintf(stderr,"%s : Nom de profil inconnu\n",pr_name);
		free(pr);
		return 1;
	}
	
	/* Normalisation entre 0 et 1 */
	min = max = pr[0];
	for (i=0;i<=N_profil-1;i++){
		if (pr[i]>max) max = pr[i];
		if (pr[i]<min) min = pr[i];
	}
	
	for (i=0;i<=N_profil-1;i++){
		pr[i] = (pr[i]-min)/(max-min);
	}
	
ENREGISTREMENT:	
	/* Enregistrement du profil */
	if (!(fp = fopen(filename,"w"))){
		fprintf(stderr, "%s ligne %d : Erreur, impossible d'ouvrir %s\n",__FILE__, __LINE__,filename);
		return 1;
	}
	for(i=0;i<=N_profil-1;i++){
			fprintf(fp,"%f\n",pr[i]);
		} 
	fclose(fp);
	
	/* Libération de la mémoire */
	free(pr);
	
	return 0;
	
}


