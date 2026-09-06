#ifndef _MD2D_MATHS_H
#define _MD2D_MATHS_H


#include "std_include.h"

COMPLEX **M_equals(COMPLEX **matrix_out, COMPLEX **matrix_in, int nlign, int ncol);
COMPLEX *M_x_V(COMPLEX *, COMPLEX **, COMPLEX *, int, int);
int Number_x_Vector(COMPLEX *vector_out, COMPLEX number, COMPLEX *vector_in, int nlign);
COMPLEX *add_Vectors(COMPLEX *, COMPLEX *, COMPLEX *, int);
COMPLEX *square_Vector(COMPLEX *, COMPLEX *, int);
COMPLEX **M_x_M(COMPLEX **, COMPLEX **, COMPLEX **, int, int);
COMPLEX **Number_x_Matrix(COMPLEX **, COMPLEX, COMPLEX **, int, int);
COMPLEX **add_M(COMPLEX **, COMPLEX **, COMPLEX **, int, int);
COMPLEX **sub_M(COMPLEX **, COMPLEX **, COMPLEX **, int, int);
COMPLEX **M_Id(COMPLEX **matrix_Id, int nlign);
COMPLEX **minus_M(COMPLEX **matrix_out, COMPLEX **matrix_in, int, int);
COMPLEX **M_zero(COMPLEX **matrix_zero, int nlign);
COMPLEX **sub_Matrices(COMPLEX **, COMPLEX **, COMPLEX **, int, int);
COMPLEX **SetMatrixCol_to_Vector(COMPLEX **, int, COMPLEX *, int, int);
COMPLEX **transpose_M(COMPLEX **M, int N);
double *Re_tab1D(COMPLEX *, double *, int);
double *Im_tab1D(COMPLEX *, double *, int);
COMPLEX **invM(COMPLEX **inv, COMPLEX **A, int N);
int nolib_MxM(COMPLEX **M_out, COMPLEX **A, COMPLEX **B, int N);
int blas_MxM(COMPLEX **M_out, COMPLEX **A, COMPLEX **B, int N);
int cblas_MxM(COMPLEX **M_out, COMPLEX **A, COMPLEX **B, int N);
int blas_MxV(COMPLEX *v_out, COMPLEX **A, COMPLEX *v_in, int N);
int lapack_invM(COMPLEX **inv, COMPLEX **A, int N);
int nolib_invM(COMPLEX **inv, COMPLEX **A, int N);
int eigen_values(COMPLEX **A, COMPLEX *eig_values, COMPLEX **EigVectors, COMPLEX *eig_buffer, int N);
int lapack_eigen_values(COMPLEX **A, COMPLEX *eig_values, COMPLEX **EigVectors, COMPLEX *eig_buffer, int N);
int check_complex(int verbosity);
int matrix_operations_check();
double M_norm(COMPLEX **M, int N);
#endif /* _MD2D_MATHS_H */


