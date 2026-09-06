/*!	\file		md2D_utils.c
 *
 * 	\brief		fonctions de manipulation de vecteurs et de matrices
 */



#include "md2D_utils.h"

/*-------------------------------------------------------------------------------------*/
/*!	\fn		double md2D_chrono(struct Param_struct *par)
 *
 *	\brief	mesure du temps écoulé depuis le début de l'execution du programme
 */
/*-------------------------------------------------------------------------------------*/
double md2D_chrono(struct Param_struct *par)
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
/*	Number_x_Vector		*/
/********************************/

int Number_x_Vector(REAL complex *vector_out, REAL complex number, REAL complex *vector_in, int nlign)
{
  int i;
	
  for (i=0;i<=nlign-1;i++){
    vector_out[i] = number * vector_in[i];		
  }
	
  return 0;
}

/********************************/
/*	add_Vectors		*/
/********************************/

REAL complex *add_Vectors(REAL complex *vector_out, REAL complex *vector_1, REAL complex *vector_2, int nlign)
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

REAL complex *square_Vector(REAL complex *vector_out, REAL complex *vector_in, int nlign)
{
  int i;
	
  for (i=0;i<=nlign-1;i++){
    vector_out[i] = vector_in[i] * vector_in[i];		
  }
	
  return vector_out;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		complex **Number_x_Matrix(complex **matrix_out, complex number, complex **matrix_in, int nlign, int ncol)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/ 
REAL complex **Number_x_Matrix(REAL complex **matrix_out, REAL complex number, REAL complex **matrix_in, int nlign, int ncol)
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
/*	M_equals		*/
/********************************/

REAL complex **M_equals(REAL complex **matrix_out, REAL complex **matrix_in, int nlign, int ncol)
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
REAL complex **add_M(REAL complex **matrix_out, REAL complex **matrix1, REAL complex **matrix2, int nlign, int ncol)
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
/*!	\fn		complex **M_Id(complex **matrix_Id, int nlign)
 *
 *	\brief	retourne la matrice Identité
 */
/*-------------------------------------------------------------------------------------*/
REAL complex **M_Id(REAL complex **matrix_Id, int nlign){
  int i,j;
	
  for (i=0;i<=nlign-1;i++){
    for (j=0;j<=nlign-1;j++){
      matrix_Id[i][j] = (i==j);		
    }
  }
	
  return matrix_Id;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		complex **M_zero(complex **matrix_zero, int nlign)
 *
 *	\brief	retourne la matrice nulle
 */
/*-------------------------------------------------------------------------------------*/
REAL complex **M_zero(REAL complex **matrix_zero, int nlign){
  int i,j;
	
  for (i=0;i<=nlign-1;i++){
    for (j=0;j<=nlign-1;j++){
      matrix_zero[i][j] = 0;		
    }
  }
	
  return matrix_zero;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		complex **sub_M(complex **matrix_out, complex **matrix1, complex **matrix2, int nlign, int ncol)
 *
 *	\brief 	matrix_out = matrix1 - matrix2
 */
/*-------------------------------------------------------------------------------------*/
REAL complex **sub_M(REAL complex **matrix_out, REAL complex **matrix1, REAL complex **matrix2, int nlign, int ncol)
{
  int i,j;
	
  for (i=0;i<=nlign-1;i++){
    for (j=0;j<=ncol-1;j++){
      matrix_out[j][i] = matrix1[j][i] - matrix2[j][i];		
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
REAL complex **sub_Matrices(REAL complex **matrix_out, REAL complex **matrix1, REAL complex **matrix2, int nlign, int ncol)
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
/*!	\fn		complex *M_x_V(complex *vector_out, complex **matrix, complex *vector_in, int nlign, int ncol)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
REAL complex *M_x_V(REAL complex *vector_out, REAL complex **matrix, REAL complex *vector_in, int nlign, int ncol)
{
#if BLAS_OPTIMIZATION
  if (nlign != ncol){
    fprintf (stderr, "%s : Error, M_x_V nline must equals ncol (edit the source to change this)", __FILE__); 
    exit (EXIT_FAILURE);}
  blas_MxV(vector_out, matrix, vector_in, ncol);
#elif ACML_OPTIMIZATION
  int i, j;
  doublecomplex alpha;
  doublecomplex beta;
  doublecomplex *m1;

  alpha.real = 1.0;
  alpha.imag = 0.0;

  beta.real = 0.0;
  beta.imag = 0.0;

  m1 = malloc(sizeof(doublecomplex) * nlign * ncol);

  /* Copie de la matrice */
  for(i=0;i<nlign;i++){
    for(j=0;j<ncol;j++){
      m1[i + ncol * j].real = creal(matrix[i][j]);
      m1[i + ncol * j].imag = cimag(matrix[i][j]);
    }
  }

  zgemv('N', nlign, ncol, &alpha, m1, nlign, (doublecomplex *)vector_in, 1, &beta, (doublecomplex *)vector_out, 1);

  free(m1);
#else
  int i,j;
  for (i=0;i<=nlign-1;i++){
    vector_out[i] = 0;
    for (j=0;j<=ncol-1;j++){
      vector_out[i] += matrix[i][j]*vector_in[j];		
    }
  }

#endif
  return vector_out;
}

#if BLAS_OPTIMIZATION
/*-------------------------------------------------------------------------------------*/
/*!	\fn		int blas_MxV(complex *v_out, complex **A, complex *v_in, int N)
 *
 *		\brief
 */
/*-------------------------------------------------------------------------------------*/
int blas_MxV(REAL complex *v_out, REAL complex **A, REAL complex *v_in, int N)
{
  /* Produit matrice.vecteur */
  /*
    void cblas_zgemv (
    const enum CBLAS_ORDER order, 
    const enum CBLAS_TRANSPOSE TransA, 
    const int M, const int N, 
    const void * alpha, 
    const void * A, const int lda, 
    const void * x, const int incx, 
    const void * beta, 
    void * y, const int incy) */

  REAL complex alpha = 1.0;
  REAL complex beta = 0.0;

  cblas_zgemv (
	       CblasRowMajor,
	       CblasNoTrans,
	       N, N,
	       &alpha,
	       &A[0][0], N, 
	       &v_in[0], 1,
	       &beta,
	       &v_out[0], 1);

  return 0;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int blas_MxM(complex **M_out, complex **A, complex **B, int N)
 *
 *	\brief	Matrix product using BLAS Library
 */
/*-------------------------------------------------------------------------------------*/
int blas_MxM(REAL complex **M_out, REAL complex **A, REAL complex **B, int N)
{
  REAL complex alpha = 1.0;
  REAL complex beta = 0.0;

  fprintf(stdout, "blas_MxM, there is an unsolved problem... use manual product matrix instead");
  exit(EXIT_FAILURE);

  /*void cblas_zgemm (
    const enum CBLAS_ORDER Order, 
    const enum CBLAS_TRANSPOSE TransA, 
    const enum CBLAS_TRANSPOSE TransB, 
    const int M, const int N, const int K, 
    const void * alpha, 
    const void * A, const int lda, const void * B, const int ldb,
    const void * beta, void * C, const int ldc)*/
  cblas_zgemm (
	       CblasRowMajor,
	       CblasNoTrans,
	       CblasNoTrans,
	       N, N, N, 
	       &alpha, 
	       &A[0][0], N, 
	       &B[0][0], N, 
	       &beta, 
	       &M_out[0][0], N);

  return 0;
}
#endif


/*-------------------------------------------------------------------------------------*/
/*!	\fn	int nolib_MxV(complex *vector_out, complex **matrix, complex *vector_in, int nlign, int ncol)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
int nolib_MxV(REAL complex *vector_out, REAL complex **matrix, REAL complex *vector_in, int nlign, int ncol)
{
  int i,j;
  for (i=0;i<=nlign-1;i++){
    vector_out[i] = 0;
    for (j=0;j<=ncol-1;j++){
      vector_out[i] += matrix[i][j]*vector_in[j];		
    }
  }

  return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		complex **M_x_M(complex **matrix_out, complex **matrix_1, complex **matrix_2, int nlign, int ncol)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
REAL complex **M_x_M(REAL complex **matrix_out, REAL complex **matrix_1, REAL complex **matrix_2, int nlign, int ncol)
{
#if BLAS_OPTIMIZATION
  if (nlign != ncol){
    fprintf (stderr, "%s : Error, M_x_M nline must equals ncol (edit the source to change this)", __FILE__); 
    exit (EXIT_FAILURE);}
  blas_MxM(matrix_out, matrix_1, matrix_2, ncol);
#elif ACML_OPTIMIZATION
  if (nlign != ncol){
    fprintf (stderr, "%s : Error, M_x_M nline must equals ncol (edit the source to change this)", __FILE__); 
    exit (EXIT_FAILURE);}
  acml_MxM(matrix_out, matrix_1, matrix_2, ncol);
#else
  int i,j,k;
	
  for (i=0;i<=nlign-1;i++){
    for (j=0;j<=nlign-1;j++){
      matrix_out[i][j] = 0;
      for (k=0;k<=ncol-1;k++){
	matrix_out[i][j] += matrix_1[i][k] * matrix_2[k][j];		
      }		
    }
  }
#endif	
  return matrix_out;
}



/*-------------------------------------------------------------------------------------*/
/*!	\fn	int acml_MxM(complex **M_out, complex **A, complex **B, int N)
 *
 *	\brief	Matrix product using ACML Library
 */
/*-------------------------------------------------------------------------------------*/
# if ACML_OPTIMIZATION
int acml_MxM(REAL complex **M_out, REAL complex **A, REAL complex **B, int N)
{
  int i,j;
  REAL complex *A_tmp, *B_tmp, *M_out_tmp;
  doublecomplex alpha;
  doublecomplex beta;

  alpha.real = 1.0;
  alpha.imag = 0.0;
  beta.real = 0.0;
  beta.imag = 0.0;

  A_tmp = malloc(sizeof(REAL complex) * N * N);
  B_tmp = malloc(sizeof(REAL complex) * N * N);
  M_out_tmp = malloc(sizeof(doublecomplex) * N * N);

  /* Copying A and B into col major (Fortran style) */
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      A_tmp[i + N * j] = A[i][j];
      B_tmp[i + N * j] = B[i][j];
    }
  }

  /* Computing the matrix product */
  zgemm('N', 'N', N, N, N, &alpha, (doublecomplex *)A_tmp, N, 
	(doublecomplex *)B_tmp, N, &beta, (doublecomplex *)M_out_tmp, N);

  /* Transforming M_out_tmp in row major (C style) */
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      M_out[i][j] = M_out_tmp[j * N + i];
    }
  }
  /*
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
    {
      printf("M[%d][%d]: %f + %fI\n", i, j, creal(M_out[i][j]), cimag(M_out[i][j]));
    }
  */

  free(A_tmp);
  free(B_tmp);
  free(M_out_tmp);

  return 0;
}
#endif

/*-------------------------------------------------------------------------------------*/
/*!	\fn	int eigen_values(complex **A, complex *eig_values, complex **EigVectors, complex *eig_buffer, int N)	
 *
 *	\brief	Eigen values & eigen vector of a complex matrix A
 * 			(calls the appropriate function according to the library at disposal)	
 */
/*-------------------------------------------------------------------------------------*/
int eigen_values(REAL complex **A, REAL complex *eig_values, REAL complex **EigVectors, REAL complex *eig_buffer, int N)
{
#if LAPACK_OPTIMIZATION
  lapack_eigen_values(A, eig_values, EigVectors, eig_buffer, N);
#elif ACML_OPTIMIZATION
  acml_eigen_values(A, eig_values, EigVectors, eig_buffer, N);
#else
  printf("%s, line %d : ERROR, function not available, either lapack or acml must be implemented",__FILE__,__LINE__);
  exit(EXIT_FAILURE);
#endif
  return 0;

}


/*-------------------------------------------------------------------------------------*/
/*!	\fn	lapack_int eigen_values(complex **A, complex *eig_values, complex **EigVectors, complex *eig_buffer, int N)	
 *
 *	\brief	Eigen values & eigen vector of a complex matrix A
 * 			uses LAPACK zgeev function
 * 			eig_values : eigen values tab
 *				EigVectors : matrix containing the eigen vectors
 *				eig_buffer : complex memory buffer, must be of size 50 N
 *				N : size of the matrix
 */
/*-------------------------------------------------------------------------------------*/
int lapack_eigen_values(REAL complex **A, REAL complex *eig_values, REAL complex **EigVectors, REAL complex *eig_buffer, int N)
{
#if LAPACK_OPTIMIZATION
  extern void zgeev_(char *jobvl, char *jobEigVectors, int *n, REAL complex *a, int *lda, REAL complex *w, REAL complex *vl, 
		     int *ldvl, REAL complex *EigVectors, int *ldEigVectors, REAL complex *work, int *lwork, double *rwork, int *info);

  int i,j,info,lwork;
  REAL complex  *work, tmp;
  double *rwork;
  work = eig_buffer;
  lwork = 49*N;	
  rwork = (double *) work + 98*N;

  /* Transforming A in col major (Fortran style) */
  for(i=0;i<=N-1;i++){
    for(j=i;j<=N-1;j++){
      tmp = A[i][j];
      A[i][j] = A[j][i];
      A[j][i] = tmp;
    }
  }
	
  char *jobvl = "N";
  char *jobEigVectors = "V";
	
  /* Computing A * v(j) = lambda(j) * v(j) */	
  zgeev_(
	 jobvl, /* char *jobvl */
	 jobEigVectors, /* char *jobEigVectors */
	 &N, /* int *n */
	 &A[0][0], /* complex *a */ 
	 &N, /* int *lda */
	 &eig_values[0], /* complex *w */ 
	 NULL, /* complex *vl */
	 &N, /* int *ldvl */
	 &EigVectors[0][0], /* complex *EigVectors */
	 &N, /* int *ldEigVectors */
	 &work[0], /* complex *work */ 
	 &lwork, /* int *lwork */
	 &rwork[0], /* double *rwork */
	 &info /* int *info*/
	 );

  if (info != 0){
    printf("ERROR, %s, line %d, can't compute eigen values, error_code info = %d\n",__FILE__,__LINE__,info);
    exit(EXIT_FAILURE);
  }

  /* Transforming EigVectors in row major (C style) */
  for(i=0;i<=N-1;i++){
    for(j=i;j<=N-1;j++){
      tmp = EigVectors[i][j];
      EigVectors[i][j] = EigVectors[j][i];
      EigVectors[j][i] = tmp;
    }
  }

  /* Constructing Lambda Matrix from lambda vector */
  /*for(i=0;i<=N-1;i++){
    for(j=0;j<=N-1;j++){
    Lambda[i][j]= 0;
    }
    Lambda[i][i] = lambda[i];
    }*/
  /* It can be checked that EigVectors*Lambda*inv(EigVectors) = A*/
#else
  printf("%s, line %d : ERROR, function not available, LAPACK not implemented",__FILE__,__LINE__);
  exit(EXIT_FAILURE);
#endif								
  return 0;



  /*  -- LAPACK driver routine (version 3.0) --   
      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
      Courant Institute, Argonne National Lab, and Rice University   
      June 30, 1999   


      Purpose   
      =======   

      ZGEEV computes for an N-by-N complex nonsymmetric matrix A, the   
      eigenvalues and, optionally, the left and/or right eigenvectors.   

      The right eigenvector v(j) of A satisfies   
      A * v(j) = lambda(j) * v(j)   
      where lambda(j) is its eigenvalue.   
      The left eigenvector u(j) of A satisfies   
      u(j)**H * A = lambda(j) * u(j)**H   
      where u(j)**H denotes the conjugate transpose of u(j).   

      The computed eigenvectors are normalized to have Euclidean norm   
      equal to 1 and largest component real.   

      Arguments   
      =========   

      JOBVL   (input) CHARACTER*1   
      = 'N': left eigenvectors of A are not computed;   
      = 'V': left eigenvectors of are computed.   

      JOBEigVectors   (input) CHARACTER*1   
      = 'N': right eigenvectors of A are not computed;   
      = 'V': right eigenvectors of A are computed.   

      N       (input) INTEGER   
      The order of the matrix A. N >= 0.   

      A       (input/output) COMPLEX*16 array, dimension (LDA,N)   
      On entry, the N-by-N matrix A.   
      On exit, A has been overwritten.   

      LDA     (input) INTEGER   
      The leading dimension of the array A.  LDA >= max(1,N).   

      W       (output) COMPLEX*16 array, dimension (N)   
      W contains the computed eigenvalues.   

      VL      (output) COMPLEX*16 array, dimension (LDVL,N)   
      If JOBVL = 'V', the left eigenvectors u(j) are stored one   
      after another in the columns of VL, in the same order   
      as their eigenvalues.   
      If JOBVL = 'N', VL is not referenced.   
      u(j) = VL(:,j), the j-th column of VL.   

      LDVL    (input) INTEGER   
      The leading dimension of the array VL.  LDVL >= 1; if   
      JOBVL = 'V', LDVL >= N.   

      EigVectors      (output) COMPLEX*16 array, dimension (LDEigVectors,N)   
      If JOBEigVectors = 'V', the right eigenvectors v(j) are stored one   
      after another in the columns of EigVectors, in the same order   
      as their eigenvalues.   
      If JOBEigVectors = 'N', EigVectors is not referenced.   
      v(j) = EigVectors(:,j), the j-th column of EigVectors.   

      LDEigVectors    (input) INTEGER   
      The leading dimension of the array EigVectors.  LDEigVectors >= 1; if   
      JOBEigVectors = 'V', LDEigVectors >= N.   

      WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)   
      On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

      LWORK   (input) INTEGER   
      The dimension of the array WORK.  LWORK >= max(1,2*N).   
      For good performance, LWORK must generally be larger.   

      If LWORK = -1, then a workspace query is assumed; the routine   
      only calculates the optimal size of the WORK array, returns   
      this value as the first entry of the WORK array, and no error   
      message related to LWORK is issued by XERBLA.   

      RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N)   

      INFO    (output) INTEGER   
      = 0:  successful exit   
      < 0:  if INFO = -i, the i-th argument had an illegal value.   
      > 0:  if INFO = i, the QR algorithm failed to compute all the   
      eigenvalues, and no eigenvectors have been computed;   
      elements and i+1:N of W contain eigenvalues which have   
      converged.   

  */
}					


/*-------------------------------------------------------------------------------------*/
/*!	\fn	int acml_eigen_values(complex **A, complex *eig_values, complex **EigVectors, complex *eig_buffer, int N)	
 *
 *	\brief	Eigen values & eigen vector of a complex matrix A
 * 			uses ACML zgeev function
 * 			eig_values : eigen values tab
 *				EigVectors : matrix containing the eigen vectors
 *				eig_buffer : (not used in this function)
 *				N : size of the matrix
 */
/*-------------------------------------------------------------------------------------*/
#if ACML_OPTIMIZATION
int acml_eigen_values(REAL complex **A, REAL complex *eig_values, REAL complex **EigVectors, REAL complex *eig_buffer, int N)
{
  REAL complex tmp;
  int i,j,info;

  /* Transforming A in col major (Fortran style) */
  for(i=0;i<N;i++){
    for(j=i;j<N;j++){
      tmp = A[i][j];
      A[i][j] = A[j][i];
      A[j][i] = tmp;
    }
  }

  /* Computing A * v(j) = lambda(j) * v(j) */	
  zgeev(
	'N', 
	'V', 
	N, 
	(doublecomplex *) A[0], 
	N, 
	(doublecomplex *) eig_values, 
	NULL, 
	N,
	(doublecomplex *) EigVectors[0],
	N, 
	&info);

  if (info != 0){
    printf("ERROR, %s, line %d, can't compute eigen values, error_code info = %d\n",__FILE__,__LINE__,info);
    exit(EXIT_FAILURE);
  }

  /* Transforming EigVectors in row major (C style) */
  for(i=0;i<=N-1;i++){
    for(j=i;j<=N-1;j++){
      tmp = EigVectors[i][j];
      EigVectors[i][j] = EigVectors[j][i];
      EigVectors[j][i] = tmp;
    }
  }
			
  return 0;


}					
#endif


/*-------------------------------------------------------------------------------------*/
/*!	\fn		complex **invM(complex **inv, complex **A, int N)
 *
 *		\brief	matrix inversion
 */
/*-------------------------------------------------------------------------------------*/
REAL complex **invM(REAL complex **inv, REAL complex **A, int N)
{
#if LAPACK_OPTIMIZATION
  lapack_invM(inv, A, N);
#elif ACML_OPTIMIZATION
  acml_invM(inv, A, N);
#else
  nolib_invM(inv, A, N);
#endif

  return inv;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn	int lapack_invM(complex **inv, complex **A, int N)
 *
 *		\brief	Inverse Matrix calculation using LU factorisation
 * 			uses LAPACK zgetrf & zgetri
 *				A : complex N x N matrix to invert (input)
 *				inv : complex N x N inverse matrix (output)  
 *				N : size of the matrix
 */
/*-------------------------------------------------------------------------------------*/
#if LAPACK_OPTIMIZATION
int lapack_invM(REAL complex **inv, REAL complex **A, int N)
{
  extern void zgetri_(int *n, REAL complex *a, int *lda, int *ipiv, REAL complex *work, int *lwork, int *info);
  extern void zgetrf_(int *m, int *n, REAL complex *a,	int *lda, int *ipiv, int *info);
	
  int i,j,info,lwork, *ipiv;
  REAL complex  *work, tmp;
  lwork = 100*N;	

  work = (REAL complex *) malloc(sizeof(REAL complex)*lwork);
  ipiv = (int *) malloc(sizeof(int)*N);

  /* Copying A to inv and transforming it in col major (Fortran style) */
  for(i=0;i<=N-1;i++){
    for(j=0;j<=N-1;j++){
      inv[i][j] = A[j][i];
    }
  }
	
  /* LU factorization */
  zgetrf_(
	  &N, /* int *M */
	  &N, /* int *N */
	  &inv[0][0], /* complex *a */
	  &N, /* int *lda */
	  &ipiv[0], /* int *ipiv */
	  &info); /* int *info */

  if (info != 0){
    printf("ERROR, %s, line %d, can't compute LU factorization, error_code info = %d\n",__FILE__,__LINE__,info);
    exit(EXIT_FAILURE);
  }
	
  /* Matrix inversion */
  zgetri_(
	  &N, /* int *N */
	  &inv[0][0], /* complex *A */ 
	  &N, /* int *LDA */ 
	  &ipiv[0], /* int *ipiv */ 
	  &work[0], /* complex *work */ 
	  &lwork, /* int *lwork */ 
	  &info); /* int *info */

  if (info != 0){
    printf("ERROR, %s, line %d, can't invert matrix, error_code info = %d\n",__FILE__,__LINE__,info);
    exit(EXIT_FAILURE);
  }

  /* Transforming inv in row major (C style) */
  for(i=0;i<=N-1;i++){
    for(j=i;j<=N-1;j++){
      tmp = inv[i][j];
      inv[i][j] = inv[j][i];
      inv[j][i] = tmp;
    }
  }

  free(work);
  free(ipiv);
	
  return 0;
}					
#endif

/*-------------------------------------------------------------------------------------*/
/*!	\fn	int acml_invM(complex **inv, complex **A, int N)
 *
 *		\brief	Inverse Matrix calculation using LU factorisation
 * 			uses ACML zgetrf & zgetri
 *				A : complex N x N matrix to invert (input)
 *				inv : complex N x N inverse matrix (output)  
 *				N : size of the matrix
 */
/*-------------------------------------------------------------------------------------*/
#if ACML_OPTIMIZATION
int acml_invM(REAL complex **inv, REAL complex **A, int N)
{
  /*	extern void zgetri_(int *n, complex *a, int *lda, int *ipiv, complex *work, int *lwork, int *info);
	extern void zgetrf_(int *m, int *n, complex *a,	int *lda, int *ipiv, int *info);
  */	
  int i,j,info, *ipiv;
  REAL complex tmp;

  ipiv = (int *) malloc(sizeof(int)*N);

  /* Copying A to inv and transforming it in col major (Fortran style) */
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      inv[i][j] = A[j][i];
    }
  }
	
  /* LU factorization */
  zgetrf(N, N, (doublecomplex *) inv[0], N, ipiv, &info);

  if (info != 0){
    printf("ERROR, %s, line %d, can't compute LU factorization, error_code info = %d\n",__FILE__,__LINE__,info);
    exit(EXIT_FAILURE);
  }
	
  /* Matrix inversion */
  zgetri(N, (doublecomplex *) inv[0], N, ipiv, &info);

  if (info != 0){
    printf("ERROR, %s, line %d, can't invert matrix, error_code info = %d\n",__FILE__,__LINE__,info);
    exit(EXIT_FAILURE);
  }

  /* Transforming inv in row major (C style) */
  for(i=0;i<N;i++){
    for(j=i;j<N;j++){
      tmp = inv[i][j];
      inv[i][j] = inv[j][i];
      inv[j][i] = tmp;
    }
  }

  free(ipiv);
	
  return 0;
}					
#endif

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int nolib_invM(complex **inv, complex **A, int N)
 *
 *		\brief	manual matrix inversion, to use only when there's no optimized lib at disposal
 */
/*-------------------------------------------------------------------------------------*/
int nolib_invM(REAL complex **inv, REAL complex **A, int N)
{

  void partialPivoting(REAL complex **M, REAL complex **b, int k, int N);

  int i,j,k;
  REAL complex **Id, **M, tmp;
	
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
	
  return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		void partialPivoting(complex **M, complex **b, int k, int N)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
void partialPivoting(REAL complex **M, REAL complex **b, int k, int N)
{

  int i, max=k;
  REAL complex tmp;
	
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

REAL complex **SetMatrixCol_to_Vector(REAL complex **matrix, int col_j, REAL complex *vector, int nlign, int ncol)
{
  int i;
	
  for (i=0;i<=nlign-1;i++){
    matrix[i][col_j] = vector[i];		
  }
	
  return matrix;
}

/*************************/
/*       CopyCplxTab         */
/*************************/

REAL complex *CopyCplxTab(REAL complex *tab1, REAL complex *tab2, int N)
{
  int i;
  for (i=0;i<=N-1;i++){
    tab1[i] = tab2[i];
  }
	
  return tab1;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn	int *CopyDbleTab(double *tab1, double *tab2, int N)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
int CopyDbleTab(double *tab1, double *tab2, int N)
{
  int i;
  for (i=0;i<=N-1;i++){
    tab1[i] = tab2[i];
  }
	
  return 0;
}


/****************************************/
/*	Re_tab1D(cplex_tab,real_tab,N)	*/
/****************************************/

double *Re_tab1D(REAL complex *cplex_tab, double *real_tab,int N)
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

double *Im_tab1D(REAL complex *cplex_tab, double *real_tab,int N)
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
int SaveDbleTab2file (double *tab, int N, char *filename, char *separateur)
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
int SaveCplxTab2file (REAL complex *tab, int Nlign, char *mode, char *filename, char *separateur)
{
  int i;
  FILE *fp;
	
  /* Affichage à l'écran si filename = "stdout" */										
  if(!strcmp(filename,"stdout")) fp = stdout;
  else fp = fopen(filename, "w");
	
  /* Mode = "Re"/"Im" : enregistrement de la partie réelle/Imaginaire */
  for (i=0;i<=Nlign-1;i++){
    if(!strcmp(mode,"Re")){
      fprintf(fp,"% 1.6e%s",creal(tab[i]),separateur);
    }else if(!strcmp(mode,"Im")){
      fprintf(fp,"% 1.6e%s",cimag(tab[i]),separateur);
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
int SaveMatrix2file (REAL complex **M, int Nlign, int Ncol, char *mode, char *filename)
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
REAL complex **allocate_CplxMatrix(int ncol,int nlign)
{
  int i;
  REAL complex **tabl;
	
  tabl = (REAL complex **) malloc (ncol * sizeof (REAL complex *));
  if (tabl == NULL) {
    fprintf (stderr, "%s : Error, allocate_CplxMatrix() can't allocate memory", __FILE__); 
    exit (EXIT_FAILURE);
  }
  tabl[0] = (REAL complex *) malloc (ncol*nlign * sizeof (REAL complex));
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
/*!	\fn	complex ***allocate_CplxMatrix_3(int nlign, int ncol, int ntab)
 *
 *	\brief	accès avec mat_3[i_tab][i_col][i_ligne]
 */
/*-------------------------------------------------------------------------------------*/
REAL complex ***allocate_CplxMatrix_3(int nlign, int ncol, int ntab)
{
  int i;
  REAL complex ***tabl;
	
  tabl = (REAL complex ***) malloc (ntab * sizeof (REAL complex **));
  if (tabl == NULL) {
    fprintf (stderr, "%s : Error, allocate_CplxMatrix_3() can't allocate memory", __FILE__); 
    exit (EXIT_FAILURE);
  }
  tabl[0] = (REAL complex **) malloc (ntab*ncol * sizeof (REAL complex *));
  if (tabl[0] == NULL) {
    free (tabl); 
    fprintf (stderr, "%s : Error, allocate_CplxMatrix_3() can't allocate memory", __FILE__); 
    exit (EXIT_FAILURE);
  }
  for(i = 1; i < ntab; i++){
    tabl[i] = tabl[i-1] + ncol; /* ncol*sizeof(complex **) */
  }	
	
  tabl[0][0] = (REAL complex *) malloc (ntab*ncol*nlign * sizeof (REAL complex));
  if (tabl[0][0] == NULL) {
    free(tabl[0]);
    free (tabl); 
    fprintf (stderr, "%s : Error, allocate_CplxMatrix_3() can't allocate memory", __FILE__); 
    exit (EXIT_FAILURE);
  }
  for(i = 1; i < ncol*ntab; i++){
    tabl[0][i] = tabl[0][i-1] + nlign; /* nlign*sizeof(complex *) */
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
REAL complex **reallocate_CplxMatrix (REAL complex **tabl, int ncol, int nlign)
{
  int i;
  REAL complex *tmp, **ret;
	
	
  tmp = (REAL complex *) realloc (tabl[0], ncol*nlign*sizeof(REAL complex));
  if (tmp == NULL) { 
    fprintf (stderr, "%s : Error, reallocate_CplxMatrix() can't allocate memory", __FILE__); 
    exit (EXIT_FAILURE);
  }
	
  ret = (REAL complex **) realloc (tabl, ncol*sizeof(REAL complex *));
  if (ret == NULL) {
    fprintf (stderr, "%s : Error, reallocate_CplxMatrix() can't allocate memory", __FILE__); 
    exit (EXIT_FAILURE);
  }
	
  for(i = 0; i < ncol; i++){
    ret[i] = tmp + i*nlign*sizeof(REAL complex *);
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
