/*!	\file		md1D_utils.c
 *
 * 	\brief		fonctions de manipulation de vecteurs et de matrices
 */



#include "md1D_utils.h"


/*-------------------------------------------------------------------------------------*/
/*!	\fn		double md1D_chrono(struct Param_struct *par)
 *
 *	\brief	mesure du temps écoulé depuis le début de l'execution du programme
 */
/*-------------------------------------------------------------------------------------*/
double md1D_chrono(struct Param_struct *par)
{

	time_t time1;
	double chrono_fin = CHRONO(clock(),par->clock0);
	time(&time1);
	double chrono_long = difftime(time1,par->time0);

	/* time()  permet un accès au temps à la seconde près sur une longue période  */
	/* clock() permet un accès fin au temps (<< 1s), mais est cyclique, ne marche */
	/*         sur une longue période.                                            */
	if (fabs(chrono_fin-chrono_long < 1)) {
		return chrono_fin;
	}
	
	return chrono_long;

}



/********************************/
/*	M_x_V		*/
/********************************/

complex *M_x_V(complex *vector_out, complex **matrix, complex *vector_in, int nlign, int ncol)
{
	int i,j;
	
	for (i=0;i<=nlign-1;i++){
		vector_out[i] = 0;
		for (j=0;j<=ncol-1;j++){
			vector_out[i] += matrix[i][j]*vector_in[j];		
		}
	}
	
	return vector_out;
}

/********************************/
/*	Number_x_Vector		*/
/********************************/

complex *Number_x_Vector(complex *vector_out, complex number, complex *vector_in, int nlign)
{
	int i;
	
	for (i=0;i<=nlign-1;i++){
		vector_out[i] = number * vector_in[i];		
	}
	
	return vector_out;
}

/********************************/
/*	add_Vectors		*/
/********************************/

complex *add_Vectors(complex *vector_out, complex *vector_1, complex *vector_2, int nlign)
{
	int i;
	
	for (i=0;i<=nlign-1;i++){
		vector_out[i] = vector_1[i] + vector_2[i];		
	}
	
	return vector_out;
}

/********************************/
/*	square_Vector		*/
/********************************/

complex *square_Vector(complex *vector_out, complex *vector_in, int nlign)
{
	int i;
	
	for (i=0;i<=nlign-1;i++){
		vector_out[i] = vector_in[i] * vector_in[i];		
	}
	
	return vector_out;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		complex **M_x_M(complex **matrix_out, complex **matrix_1, complex **matrix_2, int nlign, int ncol)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
complex **M_x_M(complex **matrix_out, complex **matrix_1, complex **matrix_2, int nlign, int ncol)
{
	int i,j,k;
	
	for (i=0;i<=nlign-1;i++){
		for (j=0;j<=nlign-1;j++){
			matrix_out[i][j] = 0;
			for (k=0;k<=ncol-1;k++){
				matrix_out[i][j] += matrix_1[i][k] * matrix_2[k][j];		
			}		
		}
	}
	
	return matrix_out;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		complex **Number_x_Matrix(complex **matrix_out, complex number, complex **matrix_in, int nlign, int ncol)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/ 
complex **Number_x_Matrix(complex **matrix_out, complex number, complex **matrix_in, int nlign, int ncol)
{
	int i,j;
	
	for (i=0;i<=nlign-1;i++){
		for (j=0;j<=ncol-1;j++){
			matrix_out[i][j] = number * matrix_in[i][j];		
		}
	}
	
	return matrix_out;
}

/********************************/
/*	copy_M		*/
/********************************/

complex **copy_M(complex **matrix_out, complex **matrix_in, int nlign, int ncol)
{
	int i,j;
	
	for (i=0;i<=nlign-1;i++){
		for (j=0;j<=ncol-1;j++){
			matrix_out[i][j] = matrix_in[i][j];		
		}
	}
	
	return matrix_out;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		complex **add_M(complex **matrix_out, complex **matrix1, complex **matrix2, int nlign, int ncol)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
complex **add_M(complex **matrix_out, complex **matrix1, complex **matrix2, int nlign, int ncol)
{
	int i,j;
	
	for (i=0;i<=nlign-1;i++){
		for (j=0;j<=ncol-1;j++){
			matrix_out[i][j] = matrix1[i][j] + matrix2[i][j];		
		}
	}
	
	return matrix_out;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		complex **sub_Matrices(complex **matrix_out, complex **matrix1, complex **matrix2, int nlign, int ncol)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
complex **sub_Matrices(complex **matrix_out, complex **matrix1, complex **matrix2, int nlign, int ncol)
{
	int i,j;
	
	for (i=0;i<=nlign-1;i++){
		for (j=0;j<=ncol-1;j++){
			matrix_out[i][j] = matrix1[i][j] - matrix2[i][j];		
		}
	}
	
	return matrix_out;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		complex **invM(complex **inv, complex **A, int N)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
complex **invM(complex **inv, complex **A, int N)
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


/*-------------------------------------------------------------------------------------*/
/*!	\fn		void partialPivoting(complex **M, complex **b, int k, int N)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
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





/********************************/
/*    SetMatrixCol_to_Vector    */
/********************************/

complex **SetMatrixCol_to_Vector(complex **matrix, int col_j, complex *vector, int nlign, int ncol)
{
	int i;
	
	for (i=0;i<=nlign-1;i++){
		matrix[i][col_j] = vector[i];		
	}
	
	return matrix;
}

/*************************/
/*       CopyTab         */
/*************************/

complex *CopyTab(complex *tab1, complex *tab2, int N)
{
	int i;
	for (i=0;i<=N-1;i++){
		tab1[i] = tab2[i];
	}
	
	return tab1;
}


/****************************************/
/*	Re_tab1D(cplex_tab,real_tab,N)	*/
/****************************************/

double *Re_tab1D(complex *cplex_tab, double *real_tab,int N)
{
	int i;
	for (i=0;i<=N-1;i++){
		real_tab[i] = creal(cplex_tab[i]);
	}
	
	return real_tab;
} 

/****************************************/
/*	Im_tab1D(cplex_tab,real_tab,N)	*/
/****************************************/

double *Im_tab1D(complex *cplex_tab, double *real_tab,int N)
{
	int i;
	for (i=0;i<=N-1;i++){
		real_tab[i] = cimag(cplex_tab[i]);
	}
	
	return real_tab;
} 

/************************************************/
/*	SavePlot2file(x,y,N,"path/filename")	*/
/************************************************/

int SavePlot2file (double *x, double *y, int N, char *filename)
{
	int i;
	FILE *fp;
	
	fp = fopen(filename, "w");
	
	for (i=0;i<=N-1;i++){
		fprintf(fp,"%1.6e\t%1.6e\n",x[i], y[i]);
	}

	fclose(fp);

	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int SaveDbleTab2file (char *filename, double *tab, int N, char *separateur)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
int SaveDbleTab2file (char *filename, double *tab, int N, char *separateur)
{
	int i;
	FILE *fp;
	
	/* Affichage à l'écran si filename = "stdout" */										
	if(!strcmp(filename,"stdout")) fp = stdout;
	else fp = fopen(filename, "w");
	
	for (i=0;i<=N-1;i++){
		fprintf(fp,"%1.6e%s",tab[i],separateur);
	}

	if(strcmp(filename,"stdout")) fclose(fp);

	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int SaveCplxTab2file (complex *tab, int Nlign, char *mode, char *filename)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
int SaveCplxTab2file (complex *tab, int Nlign, char *mode, char *filename)
{
	int i;
	FILE *fp;
	
	/* Affichage à l'écran si filename = "stdout" */										
	if(!strcmp(filename,"stdout")) fp = stdout;
	else fp = fopen(filename, "w");
	
	/* Mode = "Re"/"Im" : enregistrement de la partie réelle/Imaginaire */
	for (i=0;i<=Nlign-1;i++){
		if(!strcmp(mode,"Re")){
			fprintf(fp,"% 1.6e \n",creal(tab[i]));
		}else if(!strcmp(mode,"Im")){
			fprintf(fp,"% 1.6e \n",cimag(tab[i]));
		}
	}
	

	if(strcmp(filename,"stdout")) fclose(fp);

	return 0;


}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int SaveMatrix2file (complex **M, int Nlign, int Ncol, char *mode, char *filename)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
int SaveMatrix2file (complex **M, int Nlign, int Ncol, char *mode, char *filename)
{
	int i,j;
	FILE *fp;
	
	/* Affichage à l'écran si filename = "stdout" */										
	if(!strcmp(filename,"stdout")) fp = stdout;
	else fp = fopen(filename, "w");
	
	/* Mode = "Re"/"Im" : enregistrement de la partie réelle/Imaginaire */
	for (i=0;i<=Nlign-1;i++){
		for (j=0;j<=Ncol-1;j++){
			if(!strcmp(mode,"Re")){
				fprintf(fp,"% 1.6e  ",creal(M[i][j]));
			}else if(!strcmp(mode,"Im")){
				fprintf(fp,"% 1.6e  ",cimag(M[i][j]));
			}
		}
		fprintf(fp,"\n");
	}

	if(strcmp(filename,"stdout")) fclose(fp);

	return 0;


}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		complex **allocate_CplxMatrix(int ncol,int nlign)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
complex **allocate_CplxMatrix(int ncol,int nlign)
{
	int i;
	complex **tabl;
	
	tabl = (complex **) malloc (ncol * sizeof (complex *));
	if (tabl == NULL) {
		fprintf (stderr, "%s : Error, allocate_CplxMatrix() can't allocate memory", __FILE__); 
		exit (EXIT_FAILURE);
	}
	tabl[0] = (complex *) malloc (ncol*nlign * sizeof (complex));
	if (tabl[0] == NULL) {
		free (tabl); 
		fprintf (stderr, "%s : Error, allocate_CplxMatrix() can't allocate memory", __FILE__); 
		exit (EXIT_FAILURE);
	}
	for(i = 1; i < ncol; i++){
		tabl[i] = tabl[i-1] + nlign; /* nlign*sizeof(complex *) */
	}	
	return tabl;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		complex **reallocate_CplxMatrix (complex **tabl, int ncol, int nlign)
 *
 *
 *	\todo	Ne marche pas, fait planter le prog
 *
 */
/*-------------------------------------------------------------------------------------*/
complex **reallocate_CplxMatrix (complex **tabl, int ncol, int nlign)
{
	int i;
	complex *tmp, **ret;
	
	
	tmp = (complex *) realloc (tabl[0], ncol*nlign*sizeof(complex));
	if (tmp == NULL) { 
		fprintf (stderr, "%s : Error, reallocate_CplxMatrix() can't allocate memory", __FILE__); 
		exit (EXIT_FAILURE);
	}
	
	ret = (complex **) realloc (tabl, ncol*sizeof(complex *));
	if (ret == NULL) {
		fprintf (stderr, "%s : Error, reallocate_CplxMatrix() can't allocate memory", __FILE__); 
		exit (EXIT_FAILURE);
	}
	
	for(i = 0; i < ncol; i++){
		ret[i] = tmp + i*nlign*sizeof(complex *);
	}	
	return ret;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn		double **allocate_DbleMatrix(int ncol,int nlign)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
double **allocate_DbleMatrix(int ncol,int nlign)
{
	int i;
	double **tabl;
	
	tabl = (double **) malloc (ncol * sizeof (double *));
	if (tabl == NULL) {
		fprintf (stderr, "%s : Error, allocate_DbleMatrix() can't allocate memory", __FILE__); 
		exit (EXIT_FAILURE);
	}
	tabl[0] = (double *) malloc (ncol*nlign * sizeof (double));
	if (tabl[0] == NULL) {
		free (tabl); 
		fprintf (stderr, "%s : Error, allocate_DbleMatrix() can't allocate memory", __FILE__); 
		exit (EXIT_FAILURE);
	}
	for(i = 1; i < ncol; i++){
		tabl[i] = tabl[i-1] + nlign;
	}	
	return tabl;
}

