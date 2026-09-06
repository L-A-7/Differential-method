#ifndef _MD2D_UTILS_H
#define _MD2D_UTILS_H


#include "std_include.h"

double md2D_chrono(struct Param_struct *par);
REAL complex **M_equals(REAL complex **matrix_out, REAL complex **matrix_in, int nlign, int ncol);
REAL complex *M_x_V(REAL complex *, REAL complex **, REAL complex *, int, int);
int Number_x_Vector(REAL complex *vector_out, REAL complex number, REAL complex *vector_in, int nlign);
REAL complex *add_Vectors(REAL complex *, REAL complex *, REAL complex *, int);
REAL complex *square_Vector(REAL complex *, REAL complex *, int);
REAL complex **M_x_M(REAL complex **, REAL complex **, REAL complex **, int, int);
REAL complex **Number_x_Matrix(REAL complex **, REAL complex, REAL complex **, int, int);
REAL complex **add_M(REAL complex **, REAL complex **, REAL complex **, int, int);
REAL complex **sub_M(REAL complex **, REAL complex **, REAL complex **, int, int);
REAL complex **M_Id(REAL complex **matrix_Id, int nlign);
REAL complex **M_zero(REAL complex **matrix_zero, int nlign);
REAL complex **sub_Matrices(REAL complex **, REAL complex **, REAL complex **, int, int);
REAL complex **SetMatrixCol_to_Vector(REAL complex **, int, REAL complex *, int, int);
REAL complex *CopyCplxTab(REAL complex *, REAL complex *, int);
REAL complex **allocate_CplxMatrix(int, int);
REAL complex ***allocate_CplxMatrix_3(int, int, int);
REAL complex **reallocate_CplxMatrix(REAL complex **tabl, int ncol,int nlign);
int CopyDbleTab(double *, double *, int);
double **allocate_DbleMatrix(int ncol,int nlign);
double *Re_tab1D(REAL complex *, double *, int);
double *Im_tab1D(REAL complex *, double *, int);
int SavePlot2file (double *, double *, int, char *);
int SaveDbleTab2file (double *tab, int N, char *filename, char *separateur);
int SaveCplxTab2file (REAL complex *tab, int Nlign, char *mode, char *filename, char *separateur);
int SaveMatrix2file (REAL complex **, int, int, char *, char *);
REAL complex **invM(REAL complex **inv, REAL complex **A, int N);
int blas_MxM(REAL complex **M_out, REAL complex **A, REAL complex **B, int N);
int acml_MxM(REAL complex **M_out, REAL complex **A, REAL complex **B, int N);
int blas_MxV(REAL complex *v_out, REAL complex **A, REAL complex *v_in, int N);
int lapack_invM(REAL complex **inv, REAL complex **A, int N);
int acml_invM(REAL complex **inv, REAL complex **A, int N);
int nolib_invM(REAL complex **inv, REAL complex **A, int N);
int eigen_values(REAL complex **A, REAL complex *eig_values, REAL complex **EigVectors, REAL complex *eig_buffer, int N);
int lapack_eigen_values(REAL complex **A, REAL complex *eig_values, REAL complex **EigVectors, REAL complex *eig_buffer, int N);
int acml_eigen_values(REAL complex **A, REAL complex *eig_values, REAL complex **EigVectors, REAL complex *eig_buffer, int N);
#endif /* _MD2D_UTILS_H */


