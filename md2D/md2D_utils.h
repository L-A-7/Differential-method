#ifndef _MD2D_UTILS_H
#define _MD2D_UTILS_H


#include "std_include.h"

double md2D_chrono(struct Param_struct *par);
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
COMPLEX *CopyCplxTab(COMPLEX *, COMPLEX *, int);
COMPLEX **allocate_CplxMatrix(int, int);
COMPLEX ***allocate_CplxMatrix_3(int, int, int);
COMPLEX **reallocate_CplxMatrix(COMPLEX **tabl, int ncol,int nlign);
int CopyDbleTab(double *, double *, int);
double **allocate_DbleMatrix(int ncol,int nlign);
double *Re_tab1D(COMPLEX *, double *, int);
double *Im_tab1D(COMPLEX *, double *, int);
int SavePlot2file (double *, double *, int, char *);
int SaveDbleTab2file (double *tab, int N, char *filename, char *separateur, int Nmax1, char *separateur2);
int SaveCplxTab2file (COMPLEX *tab, int Nlign, char *mode, char *filename, char *separateur);
int SaveMatrix2file (COMPLEX **, int, int, char *, char *);
COMPLEX **invM(COMPLEX **inv, COMPLEX **A, int N);
int blas_MxM(COMPLEX **M_out, COMPLEX **A, COMPLEX **B, int N);
int acml_MxM(COMPLEX **M_out, COMPLEX **A, COMPLEX **B, int N);
int blas_MxV(COMPLEX *v_out, COMPLEX **A, COMPLEX *v_in, int N);
int lapack_invM(COMPLEX **inv, COMPLEX **A, int N);
int acml_invM(COMPLEX **inv, COMPLEX **A, int N);
int nolib_invM(COMPLEX **inv, COMPLEX **A, int N);
int eigen_values(COMPLEX **A, COMPLEX *eig_values, COMPLEX **EigVectors, COMPLEX *eig_buffer, int N);
int lapack_eigen_values(COMPLEX **A, COMPLEX *eig_values, COMPLEX **EigVectors, COMPLEX *eig_buffer, int N);
int acml_eigen_values(COMPLEX **A, COMPLEX *eig_values, COMPLEX **EigVectors, COMPLEX *eig_buffer, int N);
#endif /* _MD2D_UTILS_H */


