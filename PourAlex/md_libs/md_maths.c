/*!	\file		md_maths.c
 *
 * 	\brief		Mathematics operations
 *				Mainly matrices operations and FFT, using libraries, BLAS, LAPACK, FFTW, ACML
 */

#include "md_maths.h"
#include "md_io_utils.h"
#include "md_utils.h"

/*-------------------------------------------------------------------------------------*/
/*!	\fn		COMPLEX **M_x_M(COMPLEX **matrix_out, COMPLEX **matrix_1, COMPLEX **matrix_2, int nlign, int ncol)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX **M_x_M(COMPLEX **matrix_out, COMPLEX **matrix_1, COMPLEX **matrix_2, int nlign, int ncol)
{
	if (nlign != ncol){
		fprintf (stderr, "%s : Error, M_x_M. nline must equals ncol\n", __FILE__); 
		exit (EXIT_FAILURE);
	}
#ifdef _ACML
	acml_MxM(matrix_out, matrix_1, matrix_2, ncol);
#elif defined _BLAS
	blas_MxM(matrix_out, matrix_1, matrix_2, ncol);
#elif defined _CBLAS
	cblas_MxM(matrix_out, matrix_1, matrix_2, ncol);
#elif defined _NO_LOWLEVEL_MAT_LIB
	nolib_MxM(matrix_out, matrix_1, matrix_2, ncol);
#else
	fprintf (stderr, "%s, line %d: ERROR, we should not be there... Preprocessor library definition error.\n", __FILE__,__LINE__);
	exit(EXIT_FAILURE);
#endif
	return matrix_out;
}


#ifdef _ACML
/* Shakti */
int acml_MxM(COMPLEX **M_out, COMPLEX **A, COMPLEX **B, int N)
{
  int i,j;
  COMPLEX *A_tmp, *B_tmp, tmp1;
  doublecomplex alpha;
  doublecomplex beta;

  alpha.real = 1.0;
  alpha.imag = 0.0;
  beta.real = 0.0;
  beta.imag = 0.0;

  A_tmp = malloc(sizeof(COMPLEX) * N * N);
  B_tmp = malloc(sizeof(COMPLEX) * N * N);

  /* Copying A and B into col major (Fortran style) */
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      A_tmp[i + N * j] = A[i][j];
      B_tmp[i + N * j] = B[i][j];
    }
  }

  /* Computing the matrix product */
  zgemm('N', 'N', N, N, N, &alpha, (doublecomplex *)A_tmp, N, 
	(doublecomplex *)B_tmp, N, &beta, (doublecomplex *)M_out[0], N);

  /* Transforming M_out in row major (C style) */
  for(i=0;i<N-1;i++){
    for(j=i+1;j<N;j++){
      tmp1 = M_out[i][j];
      M_out[i][j] = M_out[j][i];
      M_out[j][i] = tmp1;
    }
  }

  free(A_tmp);
  free(B_tmp);

  return 0;
}
#endif


/*!------------------------------------------------------------------------------------
 *	\fn		int blas_MxM(COMPLEX **M_out, COMPLEX **A, COMPLEX **B, int N)
 *
 *	\brief	Matrix product using BLAS Library
 *
 *-------------------------------------------------------------------------------------*/
#ifdef _BLAS
int blas_MxM(COMPLEX **M_out, COMPLEX **A, COMPLEX **B, int N)
{
  int i,j;
  COMPLEX tmp1, *A_tmp, *B_tmp;
  COMPLEX alpha = 1.0;
  COMPLEX beta = 0.0;
  char *TransA = "N";
  char *TransB = "N";
  extern void zgemm_(char *TRANSA, char *TRANSB, int *, int *, int *, COMPLEX *ALPHA, COMPLEX *a, int *, COMPLEX *, int *, COMPLEX *, COMPLEX *, int *);

  /* Making a copy of A and B (some blas implementations modify their values...) */
  /* The copies are transposed matrices of A and B (because of the different matrices memory storing conventions in Fortran and C) */
  /* (transposition could also be performed by seting TRANSA and TRANSB to "T") */
  A_tmp = malloc(sizeof(COMPLEX) * N * N);
  B_tmp = malloc(sizeof(COMPLEX) * N * N);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      A_tmp[i + N * j] = A[i][j];
      B_tmp[i + N * j] = B[i][j];
    }
  }

  zgemm_(TransA, TransB, &N, &N, &N, &alpha, &A_tmp[0], &N, &B_tmp[0], &N, &beta, &M_out[0][0], &N);

  /* Transposition of resulting matrice */
  for(i=0;i<N-1;i++){
    for(j=i+1;j<N;j++){
      tmp1 = M_out[i][j];
      M_out[i][j] = M_out[j][i];
      M_out[j][i] = tmp1;
    }
  }

free(A_tmp);
free(B_tmp);

  return 0;
}
#endif

/*!------------------------------------------------------------------------------------
 *	\fn		int cblas_MxM(COMPLEX **M_out, COMPLEX **A, COMPLEX **B, int N)
 *
 *	\brief	Matrix product using CBLAS Library
 *
 *-------------------------------------------------------------------------------------*/
#ifdef _CBLAS
int cblas_MxM(COMPLEX **M_out, COMPLEX **A, COMPLEX **B, int N)
{
  int i,j;
  COMPLEX *A_tmp, *B_tmp;
  COMPLEX alpha = 1.0;
  COMPLEX beta = 0.0;

  /* Making a copy of A and B (some cblas implementations modify their values...) */
  A_tmp = malloc(sizeof(COMPLEX) * N * N);
  B_tmp = malloc(sizeof(COMPLEX) * N * N);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      A_tmp[i + N * j] = A[j][i];
      B_tmp[i + N * j] = B[j][i];
    }
  }
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
	       &A_tmp[0], N, 
	       &B_tmp[0], N, 
	       &beta, 
	       &M_out[0][0], N);

free(A_tmp);
free(B_tmp);

  return 0;
}
#endif

/*!------------------------------------------------------------------------------------
 *	\fn		int nolib_MxM(COMPLEX **M_out, COMPLEX **A, COMPLEX **B, int N)
 *
 *	\brief	Matrix product without use of a library
 *
 *-------------------------------------------------------------------------------------*/
int nolib_MxM(COMPLEX **M_out, COMPLEX **A, COMPLEX **B, int N)
{
	int i,j,k;
	int ncol = N;
	int nlign = N;
	printf("nolib_MxM");
	for (i=0;i<=nlign-1;i++){
		for (j=0;j<=nlign-1;j++){
			M_out[i][j] = 0;
			for (k=0;k<=ncol-1;k++){
				M_out[i][j] += A[i][k] * B[k][j];		
			}		
		}
	}
	return 0;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		COMPLEX *M_x_V(COMPLEX *vector_out, COMPLEX **matrix, COMPLEX *vector_in, int nlign, int ncol)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX *M_x_V(COMPLEX *vector_out, COMPLEX **matrix, COMPLEX *vector_in, int nlign, int ncol)
{
#ifdef _CBLAS
  if (nlign != ncol){
    fprintf (stderr, "%s : Error, M_x_V nline must equals ncol (edit the source to change this)", __FILE__); 
    exit (EXIT_FAILURE);
  }
  cblas_MxV(vector_out, matrix, vector_in, ncol);
#else
	#ifdef _ACML
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
printf("MxV");
  int i,j;
  for (i=0;i<=nlign-1;i++){
    vector_out[i] = 0;
    for (j=0;j<=ncol-1;j++){
      vector_out[i] += matrix[i][j]*vector_in[j];		
    }
  }

	#endif
#endif	
  return vector_out;
}

#ifdef _CBLAS
/*-------------------------------------------------------------------------------------*/
/*!	\fn		int cblas_MxV(COMPLEX *v_out, COMPLEX **A, COMPLEX *v_in, int N)
 *
 *		\brief
 */
/*-------------------------------------------------------------------------------------*/
int cblas_MxV(COMPLEX *v_out, COMPLEX **A, COMPLEX *v_in, int N)
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

  COMPLEX alpha = 1.0;
  COMPLEX beta = 0.0;

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
#endif


/*-------------------------------------------------------------------------------------*/
/*!	\fn	int nolib_MxV(COMPLEX *vector_out, COMPLEX **matrix, COMPLEX *vector_in, int nlign, int ncol)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
int nolib_MxV(COMPLEX *vector_out, COMPLEX **matrix, COMPLEX *vector_in, int nlign, int ncol)
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
/*!	\fn	int eigen_values(COMPLEX **A, COMPLEX *eig_values, COMPLEX **EigVectors, COMPLEX *eig_buffer, int N)	
 *
 *	\brief	Eigen values & eigen vector of a COMPLEX matrix A
 * 			(calls the appropriate function according to the library at disposal)	
 */
/*-------------------------------------------------------------------------------------*/
int eigen_values(COMPLEX **A, COMPLEX *eig_values, COMPLEX **EigVectors, COMPLEX *eig_buffer, int N)
{
#ifdef _LAPACK
  lapack_eigen_values(A, eig_values, EigVectors, eig_buffer, N);
#else
	#ifdef _ACML
  acml_eigen_values(A, eig_values, EigVectors, eig_buffer, N);
	#else
  printf("%s, line %d : ERROR, function not available, either lapack or acml must be implemented",__FILE__,__LINE__);
  exit(EXIT_FAILURE);
	#endif
#endif
  return 0;

}


/*-------------------------------------------------------------------------------------*/
/*!	\fn	lapack_eigen_values(COMPLEX **A, COMPLEX *eig_values, COMPLEX **EigVectors, COMPLEX *eig_buffer, int N)	
 *
 *	\brief	Eigen values & eigen vector of a COMPLEX matrix A
 * 			uses LAPACK zgeev function
 * 			eig_values : eigen values tab
 *				EigVectors : matrix containing the eigen vectors
 *				eig_buffer : COMPLEX memory buffer, must be of size 50 N
 *				N : size of the matrix
 */
/*-------------------------------------------------------------------------------------*/
int lapack_eigen_values(COMPLEX **A, COMPLEX *eig_values, COMPLEX **EigVectors, COMPLEX *eig_buffer, int N)
{
#ifdef _LAPACK
  extern void zgeev_(char *jobvl, char *jobEigVectors, int *n, COMPLEX *a, int *lda, COMPLEX *w, COMPLEX *vl, 
		     int *ldvl, COMPLEX *EigVectors, int *ldEigVectors, COMPLEX *work, int *lwork, double *rwork, int *info);

  int i,j,info,lwork;
  COMPLEX  *work, tmp;
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
	 &A[0][0], /* COMPLEX *a */ 
	 &N, /* int *lda */
	 &eig_values[0], /* COMPLEX *w */ 
	 NULL, /* COMPLEX *vl */
	 &N, /* int *ldvl */
	 &EigVectors[0][0], /* COMPLEX *EigVectors */
	 &N, /* int *ldEigVectors */
	 &work[0], /* COMPLEX *work */ 
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

      ZGEEV computes for an N-by-N COMPLEX nonsymmetric matrix A, the   
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
/*!	\fn	int acml_eigen_values(COMPLEX **A, COMPLEX *eig_values, COMPLEX **EigVectors, COMPLEX *eig_buffer, int N)	
 *
 *	\brief	Eigen values & eigen vector of a COMPLEX matrix A
 * 			uses ACML zgeev function
 * 			eig_values : eigen values tab
 *				EigVectors : matrix containing the eigen vectors
 *				eig_buffer : (not used in this function)
 *				N : size of the matrix
 */
/*-------------------------------------------------------------------------------------*/
#ifdef _ACML
int acml_eigen_values(COMPLEX **A, COMPLEX *eig_values, COMPLEX **EigVectors, COMPLEX *eig_buffer, int N)
{
  COMPLEX tmp;
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
/*!	\fn		COMPLEX **invM(COMPLEX **inv, COMPLEX **A, int N)
 *
 *		\brief	matrix inversion
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX **invM(COMPLEX **inv, COMPLEX **A, int N)
{
#ifdef _LAPACK
  lapack_invM(inv, A, N);
#else	
	#ifdef _ACML
  acml_invM(inv, A, N);
	#else
  nolib_invM(inv, A, N);
	#endif
#endif

  return inv;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn	int lapack_invM(COMPLEX **inv, COMPLEX **A, int N)
 *
 *		\brief	Inverse Matrix calculation using LU factorisation
 * 			uses LAPACK zgetrf & zgetri
 *				A : COMPLEX N x N matrix to invert (input)
 *				inv : COMPLEX N x N inverse matrix (output)  
 *				N : size of the matrix
 */
/*-------------------------------------------------------------------------------------*/
#ifdef _LAPACK
int lapack_invM(COMPLEX **inv, COMPLEX **A, int N)
{
  extern void zgetri_(int *n, COMPLEX *a, int *lda, int *ipiv, COMPLEX *work, int *lwork, int *info);
  extern void zgetrf_(int *m, int *n, COMPLEX *a,	int *lda, int *ipiv, int *info);
	
  int i,j,info,lwork, *ipiv;
  COMPLEX  *work, tmp;
  lwork = 100*N;	

  work = (COMPLEX *) malloc(sizeof(COMPLEX)*lwork);
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
	  &inv[0][0], /* COMPLEX *a */
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
	  &inv[0][0], /* COMPLEX *A */ 
	  &N, /* int *LDA */ 
	  &ipiv[0], /* int *ipiv */ 
	  &work[0], /* COMPLEX *work */ 
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
/*!	\fn	int acml_invM(COMPLEX **inv, COMPLEX **A, int N)
 *
 *		\brief	Inverse Matrix calculation using LU factorisation
 * 			uses ACML zgetrf & zgetri
 *				A : COMPLEX N x N matrix to invert (input)
 *				inv : COMPLEX N x N inverse matrix (output)  
 *				N : size of the matrix
 */
/*-------------------------------------------------------------------------------------*/
#ifdef _ACML
int acml_invM(COMPLEX **inv, COMPLEX **A, int N)
{
  /*	extern void zgetri_(int *n, COMPLEX *a, int *lda, int *ipiv, COMPLEX *work, int *lwork, int *info);
	extern void zgetrf_(int *m, int *n, COMPLEX *a,	int *lda, int *ipiv, int *info);
  */	
  int i,j,info, *ipiv;
  COMPLEX tmp;

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
/*!	\fn		int nolib_invM(COMPLEX **inv, COMPLEX **A, int N)
 *
 *		\brief	manual matrix inversion, to use only when there's no optimized lib at disposal
 */
/*-------------------------------------------------------------------------------------*/
int nolib_invM(COMPLEX **inv, COMPLEX **A, int N)
{

  void partialPivoting(COMPLEX **M, COMPLEX **b, int k, int N);

  int i,j,k;
  COMPLEX **Id, **M, tmp;
	
  M = allocate_CplxMatrix(N, N);
  Id = allocate_CplxMatrix(N, N);

  /* Copie de la matrice */
  for(i=0;i<=N-1;i++){
    for(j=0;j<=N-1;j++){
      M[i][j] = A[i][j]; 
    }
  }
		
  /* Matrice Identite */
  for(i=0;i<=N-1;i++){
    for(j=0;j<=N-1;j++){
      Id[i][j] = (i==j);
    }
  }

  /* 'Triangularisation' de la matrice */
  for(k=0;k<=N-2;k++){
    /* Permutation des lignes pour que le pivot, M[k][k], soit le  */
    /* plus grand element de la colonne. Minimise la propagation   */
    /* des erreurs d'arrondi, et evite les divisions par zero      */
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

  /* Resolution du système A.[inv] = [Id] après 'triangularisation' */		
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
/*!	\fn		void partialPivoting(COMPLEX **M, COMPLEX **b, int k, int N)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
void partialPivoting(COMPLEX **M, COMPLEX **b, int k, int N)
{

  int i, max=k;
  COMPLEX tmp;
	
  /* Recherche du plus grand element */
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
/*	Number_x_Vector		*/
/********************************/

int Number_x_Vector(COMPLEX *vector_out, COMPLEX number, COMPLEX *vector_in, int nlign)
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

COMPLEX *add_Vectors(COMPLEX *vector_out, COMPLEX *vector_1, COMPLEX *vector_2, int nlign)
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

COMPLEX *square_Vector(COMPLEX *vector_out, COMPLEX *vector_in, int nlign)
{
  int i;
	
  for (i=0;i<=nlign-1;i++){
    vector_out[i] = vector_in[i] * vector_in[i];		
  }
	
  return vector_out;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		COMPLEX **Number_x_Matrix(COMPLEX **matrix_out, COMPLEX number, COMPLEX **matrix_in, int nlign, int ncol)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/ 
COMPLEX **Number_x_Matrix(COMPLEX **matrix_out, COMPLEX number, COMPLEX **matrix_in, int nlign, int ncol)
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
COMPLEX **M_equals(COMPLEX **matrix_out, COMPLEX **matrix_in, int nlign, int ncol)
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
/*!	\fn		COMPLEX **add_M(COMPLEX **matrix_out, COMPLEX **matrix1, COMPLEX **matrix2, int nlign, int ncol)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX **add_M(COMPLEX **matrix_out, COMPLEX **matrix1, COMPLEX **matrix2, int nlign, int ncol)
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
/*!	\fn		COMPLEX **M_Id(COMPLEX **matrix_Id, int nlign)
 *
 *	\brief	retourne la matrice Identite
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX **M_Id(COMPLEX **matrix_Id, int nlign){
  int i,j;
	
  for (i=0;i<=nlign-1;i++){
    for (j=0;j<=nlign-1;j++){
      matrix_Id[i][j] = (i==j);		
    }
  }
	
  return matrix_Id;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		COMPLEX **minus_M(COMPLEX **matrix_out, COMPLEX **matrix_in, int nlign)
 *
 *	\brief	retourne la matrice Identite
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX **minus_M(COMPLEX **matrix_out, COMPLEX **matrix_in, int nlign, int ncol)
{
  int i,j;
	
  for (i=0;i<=nlign-1;i++){
    for (j=0;j<=ncol-1;j++){
      matrix_out[i][j] = -matrix_in[i][j];		
    }
  }
	
  return matrix_out;
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		COMPLEX **M_zero(COMPLEX **matrix_zero, int nlign)
 *
 *	\brief	retourne la matrice nulle
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX **M_zero(COMPLEX **matrix_zero, int nlign){
  int i,j;
	
  for (i=0;i<=nlign-1;i++){
    for (j=0;j<=nlign-1;j++){
      matrix_zero[i][j] = 0;		
    }
  }
	
  return matrix_zero;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		COMPLEX **sub_M(COMPLEX **matrix_out, COMPLEX **matrix1, COMPLEX **matrix2, int nlign, int ncol)
 *
 *	\brief 	matrix_out = matrix1 - matrix2
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX **sub_M(COMPLEX **matrix_out, COMPLEX **matrix1, COMPLEX **matrix2, int nlign, int ncol)
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
/*!	\fn		COMPLEX **sub_Matrices(COMPLEX **matrix_out, COMPLEX **matrix1, COMPLEX **matrix2, int nlign, int ncol)
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
COMPLEX **sub_Matrices(COMPLEX **matrix_out, COMPLEX **matrix1, COMPLEX **matrix2, int nlign, int ncol)
{
  int i,j;
	
  for (i=0;i<=nlign-1;i++){
    for (j=0;j<=ncol-1;j++){
      matrix_out[i][j] = matrix1[i][j] - matrix2[i][j];		
    }
  }
	
  return matrix_out;
}






/********************************/
/*    SetMatrixCol_to_Vector    */
/********************************/

COMPLEX **SetMatrixCol_to_Vector(COMPLEX **matrix, int col_j, COMPLEX *vector, int nlign, int ncol)
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

COMPLEX *CopyCplxTab(COMPLEX *tab1, COMPLEX *tab2, int N)
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

double *Re_tab1D(COMPLEX *cplex_tab, double *real_tab,int N)
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

double *Im_tab1D(COMPLEX *cplex_tab, double *real_tab,int N)
{
  int i;
  for (i=0;i<=N-1;i++){
    real_tab[i] = cimag(cplex_tab[i]);
  }
	
  return real_tab;
} 


/*!-------------------------------------------------------------------------------------
 *	\fn	COMPLEX **transpose_M(COMPLEX **M, int N)
 *
 *	\brief Calculate the transpose (non conjugated) of a complex matrix
 *
 *--------------------------------------------------------------------------------------*/
COMPLEX **transpose_M(COMPLEX **M, int N)
{
	int i,j;
	COMPLEX tmp;

	for(i=0;i<N-1;i++){
		for(j=i+1;j<N;j++){
			tmp = M[i][j];
			M[i][j] = M[j][i];
			M[j][i] = tmp;
		}
	}
	
	return M;
}


/*-------------------------------------------------------------------------------------*/
/*!	\fn		int check_complex()
 *
 *	\brief	Check that complex real and imag part are stored consecutively in memory for both COMPLEX type and struct{REAL real;REAL imag} typedef
 */
/*-------------------------------------------------------------------------------------*/
int check_complex(int verbosity)
{
	COMPLEX *A;
	typedef struct {double real; double imag;} doublecomplex;
	doublecomplex *B;
	double error,*a,*b;
	int n;
	A=(COMPLEX *) malloc(sizeof(COMPLEX)*3);
	B=(doublecomplex *) malloc(sizeof(doublecomplex)*3);
	a=(double *)A;
	b=(double *)B;
	A[0]=10+I*1;
	A[1]=20+I*2;
	A[2]=30+I*3;
	B[0].real=10;
	B[0].imag=1;
	B[1].real=20;
	B[1].imag=2;
	B[2].real=30;
	B[2].imag=3;

	error=0;
	for (n=0;n<6;n++){
		error+=fabs(a[n]-b[n]);
	}
	free(A);
	free(B);

	fprintf(stdout,"check_complex: ");
	if (error < 1e-18){
		fprintf(stdout,"OK (error = %.1e)\n",error);
		return 0;
	}else{
		fprintf(stdout,"ERROR. Complex number don't seem properly aligned in the memory (error = %.1e)\n",error);
		return 1;
	}
}

/*-------------------------------------------------------------------------------------*/
/*!	\fn		int matrix_operations_check()
 *
 *	\brief
 */
/*-------------------------------------------------------------------------------------*/
int matrix_operations_check()
{
	int i,j,N=30;
	double error;
	COMPLEX **A,**B,**C, **Id, **A1, **A2, **A3, **A4, **A5, **Mtmp1, **Mtmp2, **Mtmp3;
	A = allocate_CplxMatrix(N, N);
	B = allocate_CplxMatrix(N, N);
	C = allocate_CplxMatrix(N, N);
	Id = allocate_CplxMatrix(N, N);
	A1 = allocate_CplxMatrix(3, 3);
	A2 = allocate_CplxMatrix(3, 3);
	A3 = allocate_CplxMatrix(3, 3);
	A4 = allocate_CplxMatrix(3, 3);
	A5 = allocate_CplxMatrix(3, 3);
	Mtmp1 = allocate_CplxMatrix(N, N);
	Mtmp2 = allocate_CplxMatrix(N, N);
	Mtmp3 = allocate_CplxMatrix(N, N);

	srand(time(NULL));
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			A[j][i] = (double) rand()/RAND_MAX + I*(double) rand()/RAND_MAX;
			B[j][i] = (double) rand()/RAND_MAX + I*(double) rand()/RAND_MAX;
			Id[j][i] = 0;
		}
		Id[i][i] = 1;
	}


	fprintf(stdout,"\nMatrix operations checking:\n");
	fprintf(stdout,"(Errors are calculated as the sum of the absolute values of all elements of the resulting matrix, divided by the number of elements)\n\n");

	/* Checking that A - M_equals(A) = 0 */
	fprintf(stdout,"EQUALITY: Checking that A - M_equals(A) = 0... ");
	sub_M(C,A,M_equals(Mtmp1, A, N, N),N,N);
	error = M_norm(C,N)/(N*N);
	fprintf(stdout,"Error = %e\n",error);

	/* Checking that A - ((A)t)t = 0 */
	fprintf(stdout,"TRANSPOSE: Checking that A - ((A)t)t = 0... ");
	M_equals(Mtmp1, A, N, N);
	sub_M(C,A,transpose_M(transpose_M(Mtmp1,N),N),N,N);
	error = M_norm(C,N)/(N*N);
	fprintf(stdout,"Error = %e\n",error);

	/* Checking that A*Id - A = 0 */
	fprintf(stdout,"PRODUCT: Checking that Id*A - A = 0... ");
	sub_M(C, M_x_M(Mtmp1, Id, A, N, N), A, N, N);
	error = M_norm(C,N)/(N*N);
	fprintf(stdout,"Error = %e\n",error);

	/* Checking that A*Id - A = 0 */
	fprintf(stdout,"PRODUCT: Checking that A*Id - A = 0... ");
	sub_M(C, M_x_M(Mtmp1, A, Id, N, N), A, N, N);
/*fprintf(stdout,"\n");
SaveMatrix2file (Mtmp1, N, N, "Re", "stdout");fprintf(stdout,"\n");
SaveMatrix2file (A, N, N, "Re", "stdout");fprintf(stdout,"\n");
SaveMatrix2file (Mtmp1, N, N, "Im", "stdout");fprintf(stdout,"\n");
SaveMatrix2file (A, N, N, "Im", "stdout");fprintf(stdout,"\n");*/
	error = M_norm(C,N)/(N*N);
	fprintf(stdout,"Error = %e\n",error);

	/* Checking that (A)t*(B)t - (B*A)t = 0 */
	fprintf(stdout,"PRODUCT: Checking that (A)t*(B)t - (B*A)t = 0... ");
	M_equals(Mtmp1, A, N, N);
	M_equals(Mtmp2, B, N, N);
	M_x_M(Mtmp3, transpose_M(Mtmp1,N), transpose_M(Mtmp2,N), N, N);
	M_x_M(Mtmp1, B, A, N, N);
	transpose_M(Mtmp1, N);
	sub_M(C, Mtmp3, Mtmp1, N, N);
	error = M_norm(C,N)/(N*N);
	fprintf(stdout,"Error = %e\n",error);

	/* Checking that A1*A2 - A3 = 0 */
	A1[0][0]=1*(2+I);A1[0][1]=2*(2+I);A1[0][2]=3*(2+I);
	A1[1][0]=4*(2+I);A1[1][1]=5*(2+I);A1[1][2]=6*(2+I);
	A1[2][0]=7*(2+I);A1[2][1]=8*(2+I);A1[2][2]=9*(2+I);
	A2[0][0]=3*(2+I);A2[0][1]=7*(2+I);A2[0][2]=8*(2+I);
	A2[1][0]=9*(2+I);A2[1][1]=6*(2+I);A2[1][2]=5*(2+I);
	A2[2][0]=4*(2+I);A2[2][1]=2*(2+I);A2[2][2]=1*(2+I);
    A3[0][0]= 99+I*132;A3[0][1]= 75+I*100;A3[0][2]= 63+I*84;
    A3[1][0]=243+I*324;A3[1][1]=210+I*280;A3[1][2]=189+I*252;
    A3[2][0]=387+I*516;A3[2][1]=345+I*460;A3[2][2]=315+I*420;
	fprintf(stdout,"PRODUCT: Checking that A1*A2 - A3 = 0... ");
	sub_M(A5, M_x_M(A4, A1, A2, 3, 3), A3, 3, 3);
	error = M_norm(A5,3)/(3*3);
	fprintf(stdout,"Error = %e\n",error);
/*	SaveMatrix2file (A1, 3, 3, "Re", "stdout");fprintf(stdout,"\n");
	SaveMatrix2file (A2, 3, 3, "Re", "stdout");fprintf(stdout,"\n");
	SaveMatrix2file (A3, 3, 3, "Re", "stdout");fprintf(stdout,"\n");
	SaveMatrix2file (Mtmp1, 3, 3, "Re", "stdout");fprintf(stdout,"\n");*/
/*	error = M_norm(A,N)/(N*N);
	fprintf(stdout,"Error A= %e\n",error);
	error = M_norm(B,N)/(N*N);
	fprintf(stdout,"Error B= %e\n",error);
	error = M_norm(Mtmp3,N)/(N*N);
	fprintf(stdout,"Error tmp3= %e\n",error);
	error = M_norm(Mtmp1,N)/(N*N);
	fprintf(stdout,"Error tmp1= %e\n",error);*/

	/* Checking that initial matrix values are preserved after a matrix product */
	fprintf(stdout,"PRODUCT: Checking that initial matrix values are preserved after a matrix product... ");
	M_equals(Mtmp1, A, N, N);
	M_equals(Mtmp2, B, N, N);
	M_x_M(Mtmp3, A, B, N, N);
	error = M_norm(sub_M(C,A,Mtmp1,N,N),N)/(N*N);
	fprintf(stdout,"Error A = %e\n",error);
	error = M_norm(sub_M(C,B,Mtmp2,N,N),N)/(N*N);
	fprintf(stdout,"Error B = %e\n",error);


	free(A[0]); free(A);
	free(B[0]); free(B);
	free(C[0]); free(C);
	free(Id[0]); free(Id);
	free(A1[0]); free(A1);
	free(A2[0]); free(A2);
	free(A3[0]); free(A3);
	free(A4[0]); free(A4);
	free(A5[0]); free(A5);
	free(Mtmp1[0]); free(Mtmp1);
	free(Mtmp2[0]); free(Mtmp2);
	free(Mtmp3[0]); free(Mtmp3);
	
	return 0;
}

double M_norm(COMPLEX **M, int N)
{
	double norm = 0.0;
	int i,j;

	for(i=0;i<N;i++){
		for(j=i;j<N;j++){
			norm += cabs(M[j][i]);
		}
	}

	return norm;
}

