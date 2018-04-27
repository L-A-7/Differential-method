/*!	\file		md_utils.c
 *
 * 	\brief		fonctions de manipulation de vecteurs et de matrices
 */



#include "md_utils.h"



/*!-------------------------------------------------------------------------------------
 *	\fn		double md_chrono(struct Param_struct *par)
 *
 *	\brief	mesure du temps écoulé depuis le début de l'execution du programme
 *
 *--------------------------------------------------------------------------------------*/
double md_chrono(struct Param_struct *par)
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

/************************************************/
/*	SavePlot2file(x,y,N,"path/filename")	*/
/************************************************/

int SavePlot2file (double *x, double *y, int N, char *filename)
{
  int i;
  FILE *fp;
	
  fp = fopen(filename, "w");
	
  for (i=0;i<=N-1;i++){
    fprintf(fp,"%1.8e\t%1.8e\n",x[i], y[i]);
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
int SaveDbleTab2file (double *tab, int N, char *filename, char *separateur, int Nmax1, char *separateur2)
{
  int i,cpt=0;
  FILE *fp;
	
  /* Affichage à l'écran si filename = "stdout" */										
  if(!strcmp(filename,"stdout")) fp = stdout;
  else fp = fopen(filename, "w");
	
  for (i=0;i<=N-1;i++){
    fprintf(fp,"%1.8e%s",tab[i],separateur);
    if (++cpt >= Nmax1){
    	fprintf(fp,"%s",separateur2);
	cpt=0;
    }
  }

  if(strcmp(filename,"stdout")) fclose(fp);

  return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int SaveCplxTab2file (COMPLEX *tab, int Nlign, char *mode, char *filename, char *separateur, int Nmax1, char *separateur2)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
int SaveCplxTab2file (COMPLEX *tab, int Nlign, char *mode, char *filename, char *separateur, int Nmax1, char *separateur2)
{
  int i,cpt=0;
  FILE *fp;
	
  /* Affichage à l'écran si filename = "stdout" */										
  if(!strcmp(filename,"stdout")) fp = stdout;
  else fp = fopen(filename, "w");
	
  /* Mode = "Re"/"Im" : enregistrement de la partie réelle/Imaginaire */
  for (i=0;i<=Nlign-1;i++){
    if(!strcmp(mode,"Re")){
      fprintf(fp,"% 1.8e%s",creal(tab[i]),separateur);
    }else if(!strcmp(mode,"Im")){
      fprintf(fp,"% 1.8e%s",cimag(tab[i]),separateur);
    }
    if (++cpt >= Nmax1){
    	fprintf(fp,"%s",separateur2);
	cpt=0;
    }
  }
	

  if(strcmp(filename,"stdout")) fclose(fp);

  return 0;


}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int SaveMatrix2file (COMPLEX **M, int Nlign, int Ncol, char *mode, char *filename)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
int SaveMatrix2file (COMPLEX **M, int Nlign, int Ncol, char *mode, char *filename)
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
/*	fprintf(fp,"% 1.8e  ",creal(M[i][j]));*/
	fprintf(fp,"%g  ",creal(M[i][j]));
      }else if(!strcmp(mode,"Im")){
/*	fprintf(fp,"% 1.8e  ",cimag(M[i][j]));*/
	fprintf(fp,"%g  ",cimag(M[i][j]));
      }
    }
    fprintf(fp,"\n");
  }

  if(strcmp(filename,"stdout")) fclose(fp);

  return 0;


}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		COMPLEX **allocate_CplxMatrix(int nlin,int ncol)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX **allocate_CplxMatrix(int nlin,int ncol)
{
  int i;
  COMPLEX **tabl;
	
  tabl = (COMPLEX **) malloc (nlin * sizeof (COMPLEX *));
  if (tabl == NULL) {
    fprintf (stderr, "%s : Error, allocate_CplxMatrix() can't allocate memory\n", __FILE__); 
    exit (EXIT_FAILURE);
  }
  tabl[0] = (COMPLEX *) malloc (ncol*nlin * sizeof (COMPLEX));
  if (tabl[0] == NULL) {
    free (tabl); 
    fprintf (stderr, "%s : Error, allocate_CplxMatrix() can't allocate memory\n", __FILE__); 
    exit (EXIT_FAILURE);
  }
  for(i = 1; i < nlin; i++){
    tabl[i] = tabl[i-1] + ncol; /* nlign*sizeof(COMPLEX *) */
  }	
  return tabl;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn	COMPLEX ***allocate_CplxMatrix_3(int nlign, int ncol, int ntab)
 *
 *	\brief	accès avec mat_3[i_tab][i_col][i_ligne]
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX ***allocate_CplxMatrix_3(int nlign, int ncol, int ntab)
{
  int i;
  COMPLEX ***tabl;
	
  tabl = (COMPLEX ***) malloc (ntab * sizeof (COMPLEX **));
  if (tabl == NULL) {
    fprintf (stderr, "%s : Error, allocate_CplxMatrix_3() can't allocate memory\n", __FILE__); 
    exit (EXIT_FAILURE);
  }
  tabl[0] = (COMPLEX **) malloc (ntab*ncol * sizeof (COMPLEX *));
  if (tabl[0] == NULL) {
    free (tabl); 
    fprintf (stderr, "%s : Error, allocate_CplxMatrix_3() can't allocate memory\n", __FILE__); 
    exit (EXIT_FAILURE);
  }
  for(i = 1; i < ntab; i++){
    tabl[i] = tabl[i-1] + ncol; /* ncol*sizeof(COMPLEX **) */
  }	
	
  tabl[0][0] = (COMPLEX *) malloc (ntab*ncol*nlign * sizeof (COMPLEX));
  if (tabl[0][0] == NULL) {
    free(tabl[0]);
    free (tabl); 
    fprintf (stderr, "%s : Error, allocate_CplxMatrix_3() can't allocate memory\n", __FILE__); 
    exit (EXIT_FAILURE);
  }
  for(i = 1; i < ncol*ntab; i++){
    tabl[0][i] = tabl[0][i-1] + nlign; /* nlign*sizeof(COMPLEX *) */
  }	
  return tabl;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		COMPLEX **reallocate_CplxMatrix (COMPLEX **tabl, int ncol, int nlign)
 *
 *
 *	\todo	Ne marche pas, fait planter le prog. Permutter ncol et nlign
 *
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX **reallocate_CplxMatrix (COMPLEX **tabl, int ncol, int nlign)
{
  int i;
  COMPLEX *tmp, **ret;
	
	
  tmp = (COMPLEX *) realloc (tabl[0], ncol*nlign*sizeof(COMPLEX));
  if (tmp == NULL) { 
    fprintf (stderr, "%s : Error, reallocate_CplxMatrix() can't allocate memory\n", __FILE__); 
    exit (EXIT_FAILURE);
  }
	
  ret = (COMPLEX **) realloc (tabl, ncol*sizeof(COMPLEX *));
  if (ret == NULL) {
    fprintf (stderr, "%s : Error, reallocate_CplxMatrix() can't allocate memory\n", __FILE__); 
    exit (EXIT_FAILURE);
  }
	
  for(i = 0; i < ncol; i++){
    ret[i] = tmp + i*nlign*sizeof(COMPLEX *);
  }	
  return ret;
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

/*-------------------------------------------------------------------------------------*/
/*!	\fn	double ***allocate_DbleMatrix_3(int nlign, int ncol, int ntab)
 *
 *	\brief	acces avec mat_3[i_tab][i_col][i_ligne]
 */
/*-------------------------------------------------------------------------------------*/
double ***allocate_DbleMatrix_3(int nlign, int ncol, int ntab)
{
  int i;
  double ***tabl;
	
  tabl = (double ***) malloc (ntab * sizeof (double **));
  if (tabl == NULL) {
    fprintf (stderr, "%s : Error, allocate_DbleMatrix_3() can't allocate memory\n", __FILE__); 
    exit (EXIT_FAILURE);
  }
  tabl[0] = (double **) malloc (ntab*ncol * sizeof (double *));
  if (tabl[0] == NULL) {
    free (tabl); 
    fprintf (stderr, "%s : Error, allocate_DbleMatrix_3() can't allocate memory\n", __FILE__); 
    exit (EXIT_FAILURE);
  }
  for(i = 1; i < ntab; i++){
    tabl[i] = tabl[i-1] + ncol; /* ncol*sizeof(double **) */
  }	
	
  tabl[0][0] = (double *) malloc (ntab*ncol*nlign * sizeof (double));
  if (tabl[0][0] == NULL) {
    free(tabl[0]);
    free (tabl); 
    fprintf (stderr, "%s : Error, allocate_DbleMatrix_3() can't allocate memory\n", __FILE__); 
    exit (EXIT_FAILURE);
  }
  for(i = 1; i < ncol*ntab; i++){
    tabl[0][i] = tabl[0][i-1] + nlign; /* nlign*sizeof(double *) */
  }	
  return tabl;
}


