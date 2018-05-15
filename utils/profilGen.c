/*!	\file	profilGenMULTI.c
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#define MAX(a,b) ((a>b)?a:b)
#define ROUND(x) ((int)(x<0 ? x-0.5 : x+0.5))
#define MAX(a,b) ((a>b)?a:b)
#define MIN(a,b) ((a<b)?a:b)
#define SIZE_STR 200
#define SIZE_STR_BUFFER 200
#define SIZE_LINE_BUFFER 50000
#define CHAR_COMMENT '#'
#ifndef PI
#define PI 3.14159265358979323846
#endif /* PI */

int lire_int(FILE *fp, const char *label, int *value);
int lire_tab(const char *nom_fichier, const char *label, double *tab, int N);
char *label_search(char *str,const char *label);
int lire_str_arg(char *dest, char *label, int argc, char **argvcp);
int lire_dble_arg(double *res, char *label, int argc, char **argvcp);
int lire_int_arg(int *res, char *label, int argc, char **argvcp);
double **allocate_DbleMatrix(int nlign,int ncol);
void skip_comment(char *str_in_out);
int lire_ligne(FILE *fp, char *line);

void message_erreur(){
	fprintf(stderr,	"  profilGen PROFIL [options] > fichier\n"
					"  PROFIL : nom du profil\n"
					"  Arguments optionnels :\n"
					"  -N_profil N; nombre de points du profil, par defaut 256\n\n"
					);
}

int main(int argc, char *argv[]){
	
	int i, j, n, N_layers, N_pr;
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
	int N_profil = 512;
	lire_int_arg(&N_profil, "-N_profil", argc, argvcp);
	strncpy(pr_name,argvcp[1],SIZE_STR);



	/*----- Génération d'un profil -----*/
	if (!strcmp(pr_name,"CMINCES")){
	
	/*     ---------------------------------  -
	                                          |
	                                          | h1
	   pr1 ---------------------------------  -
	                                          | h2 = 1 - h1
	   pr2 ---------------------------------  -	
	*/
fprintf(stderr,"ERREUR : Ancienne convention des axes, reprogrammer la fonction\n");exit(EXIT_FAILURE);
		/* Allocation de mémoire pour le profil*/
		N_layers = 0;
		pr = allocate_DbleMatrix(N_layers+1, N_profil);

		double h1 = 0;
		lire_dble_arg(&h1, "-h1", argc, argvcp);
	
		if (h1>1) {
			fprintf(stderr,"CMINCES : Problème de paramètre(s)");
			free(pr);
			return 1;
		}
	
		for(i=0;i<=N_profil-1;i++){
			pr[0][i] = h1;
			pr[1][i] = 1;
		}		

	}else if (!strcmp(pr_name,"CARRE01")){
	
	/*     <---- 1 - L1 ----><--- L1 ---> 
	                         +----------+     -
	                         |          |     |
	                         |          |     | h
	   pr  ------------------+          +---  -

	*/

		/* Allocation de mémoire pour le profil*/
		N_layers = 0;
		pr = allocate_DbleMatrix(N_layers+1, N_profil);


		double l, L1 = 0.50;
		double h1 = 1;
		int n_period, N_periods = 1;
		

		lire_dble_arg(&L1, "-L1", argc, argvcp);
		lire_dble_arg(&h1, "-h1", argc, argvcp);
		lire_int_arg(&N_periods, "-N_periods", argc, argvcp);
	
	
		if (L1>1 || L1<0) {
			fprintf(stderr,"CARRE01 : Problème de paramètre(s), L1 doit être compris entre 0 et 1");
			free(pr);
			return 1;
		}

		for(i=0;i<=N_profil-1;i++){
			n_period = (i/(N_profil/N_periods));
			l=((double)(i-n_period*(N_profil/N_periods)))/((double)(N_profil/N_periods));
			pr[0][i] = ((l<L1)?h1:0);
		}	
	
	
/*		for(i=0;i<=N_profil*(1-L1)-1;i++){
			pr[0][i] = h1;
		}
		for(i=N_profil*(1-L1);i<=N_profil-1;i++){
			pr[0][i] = 0;
		}
*/
	}else if (!strcmp(pr_name,"SINUS01")){
	
	/*                
	    --.           .---.   
	       `.        '     `.
	        `.      '       `.      '
	         `.___,'         `.___,'  
	*/
fprintf(stderr,"ERREUR : Ancienne convention des axes, reprogrammer la fonction\n");exit(EXIT_FAILURE);
		/* Allocation de mémoire pour le profil*/
		N_layers = 0;
		pr = allocate_DbleMatrix(N_layers+1, N_profil);
		
		double l;
		int N_periods = 16;
		
		lire_int_arg(&N_periods, "-N_periods", argc, argvcp);
		
		for(i=0;i<=N_profil-1;i++){
			l = 2*PI*((double)(N_periods*i)) / ((double) N_profil);
			pr[0][i] = 0.5*sin(l) + 0.5;
		}
	}else if (!strcmp(pr_name,"CARRE02")){
	
	/*     <---- 1 - L1 ----><--- L1 ---> 
	                         +----------+     -
	                         |          |     |
	                         |          |     | h1
	   pr1 ------------------+          +---  -
	                                          | h2 = 1 - h1
	   pr2 ---------------------------------  -	
	*/

		/* Allocation de mémoire pour le profil*/
		N_layers = 1;
		pr = allocate_DbleMatrix(N_layers+1, N_profil);

		double L1 = 0.50;
		double h1 = 1;

		lire_dble_arg(&L1, "-L1", argc, argvcp);
		lire_dble_arg(&h1, "-h1", argc, argvcp);
	
	
		if (h1>1 || L1>1) {
			fprintf(stderr,"CARRE02 : Problème de paramètre(s), L1 et h1 doivent être entre 0 et 1");
			free(pr);
			return 1;
		}
	
		for(i=0;i<=N_profil*(1-L1)-1;i++){
			pr[0][i] = (1-h1);
			pr[1][i] = 0;
		}
		for(i=N_profil*(1-L1);i<=N_profil-1;i++){
			pr[0][i] = 1;
			pr[1][i] = 0;
		}

	}else if (!strcmp(pr_name,"CARRE02_B")){
	
	/* Similar to CARRE02 but the lines are inside the substrate instead of being inside the superstrate
	   pr1 ---------------------------------  -	
	                                          | h2
	   pr2 ------------------+          +---  -
	                         |          |     |
	                         |          |     | h1
	                         +----------+     -
	       <---- 1 - L1 ----><--- L1 ---> 

	*/

		/* Allocation de mémoire pour le profil*/
		N_layers = 1;
		pr = allocate_DbleMatrix(N_layers+1, N_profil);

		double L1 = 0.50;
		double h1 = 1;
		double h2 = 0;

		lire_dble_arg(&L1, "-L1", argc, argvcp);
		lire_dble_arg(&h1, "-h1", argc, argvcp);
		lire_dble_arg(&h2, "-h2", argc, argvcp);
	
	
		if (h1>1 || L1>1) {
			fprintf(stderr,"CARRE02_B : Problème de paramètre(s), L1 et h1 doivent être entre 0 et 1");
			free(pr);
			return 1;
		}
	
		for(i=0;i<=N_profil*(1-L1)-1;i++){
			pr[0][i] = h1+h2;
			pr[1][i] = h1;
		}
		for(i=N_profil*(1-L1);i<=N_profil-1;i++){
			pr[0][i] = h1+h2;
			pr[1][i] = 0;
		}

	}else if (!strcmp(pr_name,"CARRE03")){
	
	/*     <---- 1 - L1 ----><--- L1 ---> 
	                         +----------+     -
	                         |          |     |  h1
	                         +----------+     -     
	                         |          |     |  h2
	   pr1 ------------------+          +---  -
	                                          |  h3 
	   pr2 ---------------------------------  -	
	*/

		/* Allocation de mémoire pour le profil*/
		N_layers = 2;
		pr = allocate_DbleMatrix(N_layers+1, N_profil);

		double L1 = 0.50;
		double h1 = 1;
		double h2 = 0;
		double h3 = 0;

		lire_dble_arg(&L1, "-L1", argc, argvcp);
		lire_dble_arg(&h1, "-h1", argc, argvcp);
		lire_dble_arg(&h2, "-h2", argc, argvcp);
		lire_dble_arg(&h3, "-h3", argc, argvcp);
	
	
		if (h1<0 || h2<0 || h3<0 || L1>1 || L1<0) {
			fprintf(stderr,"MULTI01 : Problème de paramètre(s), L1 doit être entre 0 et 1, h1, h2, h3 doivent etre > 0");
			free(pr);
			return 1;
		}
			
		for(i=0;i<=N_profil*(1-L1)-1;i++){
			pr[0][i] = h3;
			pr[1][i] = h3;
			pr[2][i] = 0;
		}
		for(i=N_profil*(1-L1);i<=N_profil-1;i++){
			pr[0][i] = h1+h2+h3;
			pr[1][i] = h2+h3;
			pr[2][i] = 0;
		}

	}else if (!strcmp(pr_name,"RECTANGLE04")){
	
	/*
	       <---- 1 - L1 ----><--- L1 ---> 
	                   (pr1) +----------+     -
	                         |          |     | h1
                       (pr2) +----------+     -
	                         |          |     | h2
	   pr1/2-----------------+          +---  -
	                                          | h3
       pr3 ---------------------------------  -	
                                              | h4
                                              |
	   pr4 ---------------------------------  -	
	*/

		/* Memory alocation for the profile */
		N_layers = 3;
		pr = allocate_DbleMatrix(N_layers+1, N_profil);

		double L1 = 0.50;
		double h1 = 1;
		double h2 = 0;
		double h3 = 0;
		double h4 = 0;

		lire_dble_arg(&L1, "-L1", argc, argvcp);
		lire_dble_arg(&h1, "-h1", argc, argvcp);
		lire_dble_arg(&h2, "-h2", argc, argvcp);
		lire_dble_arg(&h3, "-h3", argc, argvcp);
		lire_dble_arg(&h4, "-h4", argc, argvcp);
	
		if (h1<0 || h2<0 || h3<0 || h4<0 || L1>1 || L1<0) {
			fprintf(stderr,"%s line %d, RECTANGLE04: parameter(s) problem, L1 must be between 0 and 1 and h1 to h4 must be at least 0.\n",__FILE__, __LINE__);
			free(pr);
			return 1;
		}
	
		for(i=0;i<=N_profil*(1-L1)-1;i++){
			pr[0][i] = h3+h4;
			pr[1][i] = h3+h4;
			pr[2][i] = h4;
			pr[3][i] = 0;
		}
		for(i=N_profil*(1-L1);i<=N_profil-1;i++){
			pr[0][i] = h1+h2+h3+h4;
			pr[1][i] = h2+h3+h4;
			pr[2][i] = h4;
			pr[3][i] = 0;
		}

	}else if (!strcmp(pr_name,"RECTANGLE04_B")){
	
	/* pr1 ---------------------------------  -	
                                              | h4
                                              |
	   pr2 ---------------------------------  -	
	                                          | h3
	   pr3/4-----------------+          +---  -
	                         |          |     | h2
                       (pr3) +----------+     -
	                         |          |     | h1
	                   (pr4) +----------+     -
	       <---- 1 - L1 ----><--- L1 ---> 
	*/

		/* Memory alocation for the profile */
		N_layers = 3;
		pr = allocate_DbleMatrix(N_layers+1, N_profil);

		double L1 = 0.50;
		double h1 = 1;
		double h2 = 0;
		double h3 = 0;
		double h4 = 0;

		lire_dble_arg(&L1, "-L1", argc, argvcp);
		lire_dble_arg(&h1, "-h1", argc, argvcp);
		lire_dble_arg(&h2, "-h2", argc, argvcp);
		lire_dble_arg(&h3, "-h3", argc, argvcp);
		lire_dble_arg(&h4, "-h4", argc, argvcp);
	
		if (h1<0 || h2<0 || h3<0 || h4<0 || L1>1 || L1<0) {
			fprintf(stderr,"%s line %d, RECTANGLE04_B: parameter(s) problem, L1 must be between 0 and 1 and h1 to h4 must be at least 0.\n",__FILE__, __LINE__);
			free(pr);
			return 1;
		}
	
		for(i=0;i<=N_profil*(1-L1)-1;i++){
			pr[0][i] = h1+h2+h3+h4;
			pr[1][i] = h1+h2+h3;
			pr[2][i] = h1+h2;
			pr[3][i] = h1+h2;
		}
		for(i=N_profil*(1-L1);i<=N_profil-1;i++){
			pr[0][i] = h1+h2+h3+h4;
			pr[1][i] = h1+h2+h3;
			pr[2][i] = h1;
			pr[3][i] = 0;
		}

	}else if (!strcmp(pr_name,"FINITE_GRATING_B")){
	
	/* Similar to FINITE_GRATING but the lines are inside the substrate instead of being inside the superstrate
	   pr1 ------------------------------------ -	
	                                            | h2 = 1 - h1
	   pr2 ------------------+  +---+  +---+  + -
	                         |  |   |  |   |  | | h1
	                         +--+   +--+   +--+ -
          <--  "voids"  --><--   patterns   -->

      The patterns are repeated N_patterns times, and voids (with the length of a pattern) N_voids times
	*/

		double L1 = 0.50;
		double h1 = 1;
		int N_patterns = 1;
		int N_voids = 0;

		lire_dble_arg(&L1, "-L1", argc, argvcp);
		lire_dble_arg(&h1, "-h1", argc, argvcp);
		lire_int_arg(&N_patterns, "-N_patterns", argc, argvcp);
		lire_int_arg(&N_voids, "-N_voids", argc, argvcp);
	
		/* Allocation de mémoire pour le profil*/
		N_pr = N_profil;
		N_profil = N_pr*(N_voids+N_patterns);
		N_layers = 1;
		pr = allocate_DbleMatrix(N_layers+1, N_profil);

		if (h1>1 || L1>1) {
			fprintf(stderr,"FINITE_GRATING_B: Wrong parameters, L1 and h1 must be beetwen 0 and 1");
			free(pr);
			return 1;
		}
		int id=0;
		for(n=1;n<=N_voids;n++){
			for(i=0;i<=N_pr-1;i++){
				pr[0][id] = 1;
				pr[1][id] = h1;
				id++;
			}
		}
		for(n=1;n<=N_patterns;n++){
			for(i=0;i<=N_pr*(1-L1)-1;i++){
				pr[0][id] = 1;
				pr[1][id] = h1;
				id++;
			}
			for(i=N_pr*(1-L1);i<=N_pr-1;i++){
				pr[0][id] = 1;
				pr[1][id] = 0;
				id++;
			}
		}

		if (id != N_profil) {
			fprintf(stderr,"FINITE_GRATING_B: Error...\nid = %d and N_profil = %d (Should be equal).\nCheck the source file.",id,N_profil);
			free(pr);
			return 1;
		}

	}else if (!strcmp(pr_name,"TRAPEZ01")){
	
	/*     <--- 1 - L1 ---><--- L1 ----> 
	                      .   +----+   .      -
	                      .  /      \  .      |  
	                      . /        \ .      | h1   
	                      ./          \.      |  
	   pr  ---------------+alpha  alpha+---  -

				
	         /|      |\
           / |      | \
		    /--|      |--\
	------+   |      |   +---

		  alpha      alpha

	*/

		/* Allocation de mémoire pour le profil*/
		N_layers = 0;
		pr = allocate_DbleMatrix(N_layers+1, N_profil);

		int i0,i1,i2;
		double x;
		double L1 = 0.50;
		double Ltotal = 1.0;
		double h1 = 1;
		double alpha = 0;
		

		lire_dble_arg(&L1, "-L1", argc, argvcp);
		lire_dble_arg(&h1, "-h1", argc, argvcp);
		lire_dble_arg(&Ltotal, "-Ltotal", argc, argvcp);
		lire_dble_arg(&alpha, "-alpha", argc, argvcp);
		double tan_alpha = tan(alpha*PI/180.0);
	
		if (h1<0 || L1<0 || L1>=Ltotal) {
			fprintf(stderr,"TRAPEZ01 : Problème de paramètre(s) h1<0 || L1<0 || L1>=Ltotal");
			free(pr);
			return 1;
		}
		if (alpha <=0) {
			fprintf(stderr,"TRAPEZ01 : Problème de paramètre(s) alpha <=0");
			free(pr);
			return 1;
		}
		if (h1*tan_alpha>L1/2) {
			fprintf(stderr,"TRAPEZ01 : Problème de paramètre(s) h1*tan_alpha>L1/2");
			free(pr);
			return 1;
		}
		/* Bottom part ...*/	
		for(i=0;i<=N_profil*(1-L1)-1;i++){
			pr[0][i] = 0;
		}
		/* ascending part  */
		i0 = N_profil*(Ltotal-L1)/Ltotal-1;
		i1 = i0 + ROUND(N_profil*h1*tan_alpha/Ltotal);
printf("i1 = %d\n",i1);
		for(i=i0;i<=i1-1;i++){
			x=(i-i0)*Ltotal/N_profil;
				pr[0][i] = x/tan_alpha;
		}
		/* top part  */
		i2 = ROUND(N_profil*(Ltotal -h1*tan_alpha)/Ltotal);
printf("i2 = %d\n",i2);fflush(stdout);
		for(i=i1;i<=i2-1;i++){
				pr[0][i] = h1;
		}
		/* descending part  */
		for(i=i2;i<=N_profil-1;i++){
			x=(i-i2)*Ltotal/N_profil;
			pr[0][i] = h1-x/tan_alpha;
		}
	}else if (!strcmp(pr_name,"TRAPEZ02")){
	
	/*     <--- 1 - L1 ---><--- L1 ----> 
	                      .   +----+   .      -
	                      .  /      \  .      |  h1
	                      . +--------+ .      -     
	                      ./          \.      |  h2
	   pr1 ---------------+alpha1 alpha2+---  -
	                                          |  h3 
	   pr2 ---------------------------------  -	
				
	         /|      |\
           / |      | \
		    /--|      |--\
	------+   |      |   +---
		  alpha 1    alpha2
	*/
fprintf(stderr,"ERREUR : Ancienne convention des axes, reprogrammer la fonction\n");exit(EXIT_FAILURE);
		/* Allocation de mémoire pour le profil*/
		N_layers = 2;
		pr = allocate_DbleMatrix(N_layers+1, N_profil);
		int i0, imid;
		double x;
		double L1 = 0.50;
		double Ltotal = 1.0;
		double h1 = 1;
		double h2 = 0;
		double h3 = 0;
		double alpha = 0;
		

		lire_dble_arg(&L1, "-L1", argc, argvcp);
		lire_dble_arg(&h1, "-h1", argc, argvcp);
		lire_dble_arg(&h2, "-h2", argc, argvcp);
		lire_dble_arg(&h3, "-h3", argc, argvcp);
		lire_dble_arg(&Ltotal, "-Ltotal", argc, argvcp);
		lire_dble_arg(&alpha, "-alpha", argc, argvcp);
		double tan_alpha = tan(alpha*PI/180.0);
	
		if (h1<0 || h2<0 || h3<0 || L1>1 || L1<0) {
			fprintf(stderr,"TRAPEZ01 : Problème de paramètre(s), L1 doit être entre 0 et 1, h1, h2, h3 doivent etre > 0");
			free(pr);
			return 1;
		}
		/* Bottom part ...*/	
		for(i=0;i<=N_profil*(1-L1)-1;i++){
			pr[0][i] = h1+h2;
			pr[1][i] = h1+h2;
			pr[2][i] = h1+h2+h3;
		}
		/* _/¨ part  */
		i0 = N_profil*(1-L1);
		imid = N_profil*(1-0.5*L1);
		for(i=i0;i<=imid-1;i++){
			pr[2][i] = h1+h2+h3;
			x=(i-i0)*Ltotal/N_profil;
			if (x < h2*tan_alpha){
				pr[1][i] = h1+h2-x/tan_alpha;
			}else{
				pr[1][i] = h1;
			}
			if (x < (h1+h2)*tan_alpha){
				pr[0][i] = h1+h2-x/tan_alpha;
			}else{
				pr[0][i] = 0;
			}
		}
		/* ¨\_ part  */
		for(i=imid;i<=N_profil-1;i++){
			pr[2][i] = h1+h2+h3;
			x=(i-N_profil)*Ltotal/N_profil;
			if (x > -h2*tan_alpha){
				pr[1][i] = h1+h2+x/tan_alpha;
			}else{
				pr[1][i] = h1;
			}
			if (x > -(h1+h2)*tan_alpha){
				pr[0][i] = h1+h2+x/tan_alpha;
			}else{
				pr[0][i] = 0;
			}

		}
	}else if (!strcmp(pr_name,"ORIGAMI")){
	
	/* 

   _.-°¨°-._                   _.-°¨°-._
-°¨         ¨°-._         _.-°¨         ¨°-._         _.-
   _.-°¨°-._     ¨°-._.-°¨     _.-°¨°-._     ¨°-._.-°¨
-°¨         ¨°-._         _.-°¨         ¨°-._         _.-
                 ¨°-._.-°¨                   ¨°-._.-°¨



	                .          +
	               / \         |
	              /   \        |     
	             /     \       | h2 (height of one triangle)
 	       +    /   .   \      |
	       |   /   / \   \     |
	       |  /   /   \   \    +    
	    h2 |     /     \       |
 	       |    /       \      | h1 (height between both shapes)
	       |   /         \     |
          +  /           \    +

	*/
		/* Allocating memory */
		N_layers = 2;
		pr = allocate_DbleMatrix(N_layers+1, N_profil);
		double h1 = 0.5;
		double h2 = 0.5;

		lire_dble_arg(&h1, "-h1", argc, argvcp);
		lire_dble_arg(&h2, "-h2", argc, argvcp);

		double hstep = h2/ROUND(0.5*N_profil);

		int hlf_Nprofil = ROUND(0.5*N_profil);
fprintf(stderr, "hlf_Nprofil=%d \n", hlf_Nprofil);

		if (2*(int)(0.5*N_profil) != N_profil) {
			fprintf(stderr, "%s line %d: ERROR, ORIGAMI: N_profil=%d. Should be an even number\n",__FILE__, __LINE__, N_profil);
			free(pr[0]);free(pr);
			return 1;
		}

		/* Origami profile */	
		for(i=0;i<=hlf_Nprofil;i++){
			pr[0][i] = i*hstep + h1;
			pr[1][i] = i*hstep;
		}
		for(i=hlf_Nprofil+1;i<=N_profil-1;i++){
			pr[0][i] = h2-(i-hlf_Nprofil)*hstep + h1;
			pr[1][i] = h2-(i-hlf_Nprofil)*hstep;
		}
	}else if (!strcmp(pr_name,"ALEAT01")){
	/*    
	       __                           _
	  _/\_/  \        _  /\_            |
	_/        \_/\_  / \/   \_/\   _/   | h
	               \/           \_/     _
	*/

		double Ag = 0.50; /* Amplitude de la Gaussienne */
		double Lg = 3; /* Largeur de la Gaussienne */
		double Ae = 0.50; /* Amplitude de l'exponentielle */
		double Le = 10; /* Largeur de l'exponentielle */
		double h = 1;
		lire_dble_arg(&Ag, "-Ag", argc, argvcp);
		lire_dble_arg(&Lg, "-Lg", argc, argvcp);
		lire_dble_arg(&Ae, "-Ae", argc, argvcp);
		lire_dble_arg(&Le, "-Le", argc, argvcp);
		if(!lire_dble_arg(&h, "-h", argc, argvcp)) fprintf(stdout,"h = %f\n",h);

		/* Allocation de mémoire pour le profil*/
		N_layers = 0;
		int marge = (int)MAX(2*Lg,2*Le);
		int N_profil_plus = N_profil + 2*marge; /* On prend des marges qu'on enlève après convolution */
		pr = allocate_DbleMatrix(2, N_profil_plus); /* Allocation d'1 ligne en +, pr calculs temporaires */
		

		/* Vérification des paramètres */
		
		/* Génération de nombres aléatoires */
		srand(time(NULL));
		for(i=0;i<=N_profil_plus-1;i++){
			pr[1][i] = (double) rand()/RAND_MAX;
		}

		double x,t;
		/* Lissage avec Gaussienne plus exponentielle */		
		for(x=0; x<=N_profil_plus-1; x++){
			pr[0][(int)x] = 0;
			for(t=-N_profil_plus; t<=N_profil_plus-1; t++){
				if ( (x+t) >= 0 && (x+t) <= N_profil_plus-1){
					pr[0][(int)x] += (Ag * exp(-((t/Lg)*(t/Lg))) + Ae * exp(-(abs(t)/Le))) * pr[1][(int)(x+t)];
				}
			}
		}
		
		/* Ebavurage : on enlève les marges */
		for (i=0;i<=N_profil;i++){
			pr[0][i] = pr[0][i+marge];
		}
		
		/* Correction de la pente */
		double correc = (pr[0][N_profil] - pr[0][0])/N_profil;
		for (i=0;i<=N_profil;i++){
			pr[0][i] -= i*correc;
		}
		
		/* Normalisation entre 0 et 1 */
		double min = pr[0][0];
		double max = pr[0][0];
		for (i=0;i<=N_profil;i++){
			min = MIN(min,pr[0][i]);
			max = MAX(max,pr[0][i]);
		}
		double norm = max - min;
		for (i=0;i<=N_profil;i++){
			pr[0][i] -= min;
			pr[0][i] /= norm;
		}

		/* Normalisation entre 0 et h */
		for (i=0;i<=N_profil;i++){
			pr[0][i] *= h;
		}
		
	}else if (!strcmp(pr_name,"ADD")){
	
		/*
		profilGen ADD -f1 pr01.txt -f2 pr02.txt -fact1 0.2 -fact2 0.8 > result.txt
		
		         _________
	    ________|           pr01 x fact1
		         +
		 _   _   _   _   _  
		/ \_/ \_/ \_/ \_/   pr02 x fact2
		         =
		          _   _   _
		 _   _   / \_/ \_/   
		/ \_/ \_/  
		
		*/
fprintf(stderr,"ERREUR : Ancienne convention des axes, reprogrammer la fonction\n");exit(EXIT_FAILURE);
		char file01[SIZE_STR], file02[SIZE_STR];
		int N_layers02;
		FILE *fp01, *fp02;
		double **pr01, **pr02, *pr01_tmp, *pr02_tmp, fact1 = 1.0, fact2 = 1.0;

		lire_str_arg(file01, "-f1", argc, argvcp);
		lire_str_arg(file02, "-f2", argc, argvcp);
		lire_dble_arg(&fact1, "-fact1", argc, argvcp);
		lire_dble_arg(&fact2, "-fact2", argc, argvcp);

		/* Vérification des arguments */

		/* Ouverture des profils d'entrée */
		if (!(fp01 = fopen(file01,"r"))){
			fprintf(stderr, "%s ligne %d : ERREUR, impossible d'ouvrir %s\n",__FILE__, __LINE__,file01);
			return 1;
		}
		if (!(fp02 = fopen(file02,"r"))){
			fprintf(stderr, "%s ligne %d : ERREUR, impossible d'ouvrir %s\n",__FILE__, __LINE__,file02);
			return 1;
		}
		
		if(lire_int(fp01,"N_layers", &N_layers  ) != 0) N_layers  = 0;
		if(lire_int(fp02,"N_layers", &N_layers02) != 0) N_layers02 = 0;
/*	N_layers = 0; N_layers02 = 0;
*/		if (N_layers != N_layers02) {
			fprintf(stderr, "%s ligne %d : ERREUR, les 2 profils n'ont pas le meme nombre de couches\n",__FILE__, __LINE__);
			return 2;;
		}
		fclose(fp01);
		fclose(fp02);
		pr   = allocate_DbleMatrix(N_layers+1, N_profil);
		pr01 = allocate_DbleMatrix(N_layers+1, N_profil);
		pr02 = allocate_DbleMatrix(N_layers+1, N_profil);
		pr01_tmp = *allocate_DbleMatrix(1,(N_layers+1) * N_profil);
		pr02_tmp = *allocate_DbleMatrix(1,(N_layers+1) * N_profil);

		/* Lecture des profils d'entrée */
		lire_tab(file01, "profil", pr01_tmp, N_profil);
		lire_tab(file02, "profil", pr02_tmp, N_profil);

		/* Séparation en plsieurs interfaces */
		for (i=0; i<=N_profil-1;i++){
			for (j=0; j<=N_layers; j++){
				pr01[j][i] = pr01_tmp[i*(N_layers+1)+j];
				pr02[j][i] = pr02_tmp[i*(N_layers+1)+j];
			}
		}

		/* Somme des 2 profils */
		for (i=0;i<=N_profil-1;i++){
			for (j=0;j<=N_layers;j++){
				pr[j][i] = fact1*pr01[j][i] + fact2*pr02[j][i];
			}
		}

		free(pr01); free(pr02); free(pr01_tmp); free(pr02_tmp);

	}else{
		fprintf(stderr,"%s : Nom de profil inconnu\n",pr_name);
		free(pr[0]);free(pr);
		return 1;
	}
	
	/* Enregistrement du profil */
	fprintf(stdout,"N_x = %d\n",N_profil);			
	if (N_layers > 0){
		fprintf(stdout,"type_profil = MULTICOUCHES\n");
		fprintf(stdout,"N_layers = %d\n",N_layers);
	}else{
		fprintf(stdout,"type_profil = H_X\n");
	}
	fprintf(stdout,"profil = \n");
	for(i=0; i<=N_profil-1; i++){
		for(j=0; j<=N_layers; j++){
			fprintf(stdout,"%1.6e\t",pr[j][i]);
		}
		fprintf(stdout,"\n");
	}
	
	/* Libération de la mémoire */
	free(pr[0]); free(pr);
	
	return 0;
	
}

/*---------------------------------------------------------------------------------------------*/
/*! \fn    int lire_int(FILE *fp, char *label, int *value)
 *
 *  \brief	Lit dans le fichier pointé par *fp la valeur entiere 'value' indiquée par 'label' \n
 *			sous la forme label = value (ex.: N2 = 10)
 *  \return	0 si lecture réussie 1 sinon
 */
/*---------------------------------------------------------------------------------------------*/
int lire_int(FILE *fp, const char *label, int *value){

	char stmp[SIZE_LINE_BUFFER];
	char *pos, *pos2;

	rewind(fp);
	while(!feof(fp)){
		if (lire_ligne(fp,stmp) != 0) break;
		skip_comment(stmp); /* Vire les commentaires */
		if ((pos = label_search(stmp,label)) != NULL){ /* Recherche le label */	
			if((pos2 = strchr(pos+strlen(label),'=')) != NULL) *pos2 = ' '; /* remplace '=' par un espace */
			if (sscanf(pos+strlen(label)," %d", value) == 1){ /* Lit la valeur */
				return 0;
			}
		}
	}
	return 1;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int lire_tab(char *nom_fichier, const char *label, double *tab, int N)
 *
 *	\brief	Lit N valeurs de format double dans un fichier et les stocke dans un tableau   \n 
 *			Les valeurs doivent être séparées par un ou plusieurs espaces, tabulations     \n
 *			ou sauts de lignes et précédées d'un label éventuellement suivi d'un signe '='.\n
 *          ex. : (...) tab1 = 3.4  4.5e-3  +46  -7.6e+2 ...                               \n
 *          Remarque : Pour lire un tableau sans label, donner "" comme label.             \n
 *
 *	\return	0 si succès, 1 si le nombre d'éléments lus diffère de N ou si le label n'a pas été trouvé.
 *
 *	\todo	RENDRE PLUS ROBUSTE : PAS DE BUFFER OVERFLOW AU CAS OU IL Y A PLUS DE N LIGNES
 *
 */
/*---------------------------------------------------------------------------------------------*/
int lire_tab(const char *nom_fichier, const char *label, double *tab, int N)
{
	FILE *fp;
	char *pos, *pos2, *endptr, line[SIZE_LINE_BUFFER];
	int cpt=0, line_cpt=0;
	double tmp;
	
	if (!(fp = fopen(nom_fichier,"r"))){
		fprintf(stderr, "%s line %d: ERROR, can't open %s\n",__FILE__, __LINE__,nom_fichier);
		return 1;
	}

	/* Recherche du label */
	while(!feof(fp)){ 
		if (lire_ligne(fp,line) != 0) goto LECTURE_FINIE;
		line_cpt++;
		skip_comment(line);
		if ((pos = label_search(line,label)) != NULL){ /* on cherche le label */	
			pos += strlen(label); 
			if((pos2 = strchr(pos,'=')) != NULL) *pos2 = ' '; /* on remplace '=' par ' ' */
			goto LABEL_TROUVE;
		}
	}

	fprintf(stderr,"%s line %d: ERROR, label '%s' was not found in %s \n", __FILE__, __LINE__, label, nom_fichier);
	fclose(fp);
	return 1;

	/* Lecture des valeurs */
	while(!feof(fp)){
		if (lire_ligne(fp,line) != 0) {goto LECTURE_FINIE;}
		line_cpt++;
		skip_comment(line); 
		pos = line;
		
	LABEL_TROUVE :
		while(isspace(*pos)) pos++; /* on élimine les espaces */
		while(pos < line+strlen(line)) { /* Tant qu'on est pas à la fin de la ligne */
			tmp = strtod(pos, &endptr); /* on lit le 'double' */
			if (pos == endptr) {goto LECTURE_FINIE;} /* conversion ratée */
			if (cpt <= N+1) tab[cpt] = tmp;
			cpt ++;
			pos = endptr;
			while(isspace(*pos)) pos++; /* on élimine les espaces */
		}
	}

	LECTURE_FINIE:
	if (cpt != N) {
		fprintf(stderr, "%s line %d: ERROR, %s contains %d elements instead of %d in %s (line %d)\n", __FILE__, __LINE__, label, cpt, N, nom_fichier, line_cpt);
		fclose(fp);
		return 2;
	}
	fclose(fp);
	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int lire_dble_arg(double *res, char *label, int argc, char **argvcp)
 *
 *	\brief	Lit la valeur de l'argument de la ligne de commande indiqué sous la forme "-label valeur"
 */
/*---------------------------------------------------------------------------------------------*/
int lire_dble_arg(double *res, char *label, int argc, char **argvcp)
{
	int i;
	char strtmp[SIZE_STR_BUFFER], *endptr;
	double tmp;

	for(i=1;i<=argc-2;i++){
		if (!strcmp(argvcp[i],label)){
			strncpy(strtmp, argvcp[i+1], SIZE_STR_BUFFER);
			tmp = strtod(strtmp, &endptr);
			if (strtmp != endptr){
				*res = tmp;
				return 0; 
			}
		}
	}
	return 1;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn	int lire_int_arg(int *res, char *label, int argc, char **argvcp)
 *
 *	\brief	Lit la valeur de l'argument de la ligne de commande indiqué sous la forme "-label valeur"
 */
/*---------------------------------------------------------------------------------------------*/
int lire_int_arg(int *res, char *label, int argc, char **argvcp)
{
	int i, tmp;
	char strtmp[SIZE_STR_BUFFER], *endptr;

	for(i=1;i<=argc-2;i++){
		if (!strcmp(argvcp[i],label)){
			strncpy(strtmp, argvcp[i+1], SIZE_STR_BUFFER);
			tmp = (int) strtod(strtmp, &endptr);
			if (strtmp != endptr){
				*res = tmp;
				return 0;
			}
		}
	}
	return 1;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		int lire_str_arg(char *dest, char *label, int argc, char **argvcp)
 *
 *	\brief	Lit la valeur de l'argument de la ligne de commande indiqué sous la forme "-label valeur"
 */
/*---------------------------------------------------------------------------------------------*/
int lire_str_arg(char *dest, char *label, int argc, char **argvcp)
{
	int i;

	for(i=1;i<=argc-2;i++){
		if (!strcmp(argvcp[i],label)){
			strncpy(dest, argvcp[i+1], SIZE_STR_BUFFER);
			return 0; 
		}
	}
	return 1;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		double **allocate_DbleMatrix(int nlign,int ncol)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
double **allocate_DbleMatrix(int nlign,int ncol)
{
  int i;
  double **tabl;
	
  tabl = (double **) malloc (nlign * sizeof (double *));
  if (tabl == NULL) {
    fprintf (stderr, "%s : Error, allocate_DbleMatrix() can't allocate memory\n", __FILE__); 
    exit (EXIT_FAILURE);
  }
  tabl[0] = (double *) malloc (ncol*nlign * sizeof (double));
  if (tabl[0] == NULL) {
    free (tabl); 
    fprintf (stderr, "%s : Error, allocate_DbleMatrix() can't allocate memory\n", __FILE__); 
    exit (EXIT_FAILURE);
  }
  for(i = 1; i < nlign; i++){
    tabl[i] = tabl[i-1] + ncol;
  }	
  return tabl;
}

/*---------------------------------------------------------------------------------------------*/
/*! \fn		void skip_comment(char *str_in_out)
 *
 *  \brief	Elimine tout ce qui se trouve après un commentaire '#' dans str_in_out
 */
/*---------------------------------------------------------------------------------------------*/
void skip_comment(char *str_in_out){

	char *pos;
	/* Cherche CHAR_COMMENT et le remplace par le charactère nul '\0' */
	if((pos = strchr(str_in_out,CHAR_COMMENT)) != NULL) {
		*pos = '\0';
	}
 }
/*---------------------------------------------------------------------------------------------*/
/*! \fn		void lire_ligne(FILE *fp, char *line)
 *
 *  \brief	Lit une ligne dans un fichier et la stocke dans une chaine de charactères 
 */
/*---------------------------------------------------------------------------------------------*/
int lire_ligne(FILE *fp, char *line)
{
	if (fgets(line, SIZE_LINE_BUFFER, fp) == NULL) return 1;
	if (strlen(line) == SIZE_LINE_BUFFER-1) {
		fprintf(stderr, "%s line %d: ERROR, insufficient buffer size. Can't read more than %d characters per line.\n",__FILE__, __LINE__,SIZE_LINE_BUFFER-1);
		exit(EXIT_FAILURE);
	}
	return 0;
}

/*---------------------------------------------------------------------------------------------*/
/*!	\fn		char *label_search(char *str,const char *label)
 *
 *	\brief	Cherche un label dans une chaine de caractères, le label doit être isolé, c.a.d, \n
 *          en début de ligne ou précédé d'un espace au sens de isspace() et suivi d'un espace \n
 *          ou d'un signe '='
 *
 *	\return	La position de la 1ere occurence du label dans la chaine ou NULL si le label n'a pas été trouvé
 */
/*---------------------------------------------------------------------------------------------*/
char *label_search(char *str,const char *label)
{
	char *pos;

	/* Cas du label vide */
	if (strlen(label) == 0) return str;

	/* Recherche du label */
	while((pos=strstr(str,label)) != NULL) {
		/* Vérification que le label est en début de ligne ou précédé par un espace */
		if (pos != str && !isspace(*(pos-1))) {
			str = pos + strlen(label);
			continue;
		}
		/* Vérification que le label est suivi par un espace, saut de ligne ou signe '=' */
		if (strlen(pos) > strlen(label)) {
			if (!isspace(*(pos+strlen(label))) && *(pos+strlen(label)) != '=') {
				str = pos + strlen(label);
				continue;
			}
		}
		break;
	}
	return pos;
}


