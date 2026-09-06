#ifndef _MD1D_UTILS_H
#define _MD1D_UTILS_H


#include "std_include.h"

double md1D_chrono(struct Param_struct *par);
complex **invM(complex **inv, complex **A, int N);
complex **M_egal(complex **matrix_out, complex **matrix_in, int nlign, int ncol);
complex *M_x_V(complex *, complex **, complex *, int, int);
complex *Number_x_Vector(complex *, complex, complex *, int);
complex *add_Vectors(complex *, complex *, complex *, int);
complex *square_Vector(complex *, complex *, int);
complex **M_x_M(complex **, complex **, complex **, int, int);
complex **Number_x_Matrix(complex **, complex, complex **, int, int);
complex **add_M(complex **, complex **, complex **, int, int);
complex **sub_M(complex **, complex **, complex **, int, int);
complex **M_Id(complex **matrix_Id, int nlign);
complex **M_zero(complex **matrix_zero, int nlign);
complex **sub_Matrices(complex **, complex **, complex **, int, int);
complex **SetMatrixCol_to_Vector(complex **, int, complex *, int, int);
complex *CopyCplxTab(complex *, complex *, int);
complex **allocate_CplxMatrix(int, int);
complex ***allocate_CplxMatrix_3(int, int, int);
complex **reallocate_CplxMatrix(complex **tabl, int ncol,int nlign);
int CopyDbleTab(double *, double *, int);
double **allocate_DbleMatrix(int ncol,int nlign);
double *Re_tab1D(complex *, double *, int);
double *Im_tab1D(complex *, double *, int);
int SavePlot2file (double *, double *, int, char *);
int SaveDbleTab2file (char *filename, double *tab, int N, char *separateur);
int SaveCplxTab2file (complex *tab, int Nlign, char *mode, char *filename, char *separateur);
int SaveMatrix2file (complex **, int, int, char *, char *);
#if 0
int blas_MxM(complex **M_out, complex **A, complex **B, int N);
int blas_MxV(complex *v_out, complex **A, complex *v_in, int N);
#endif

#endif /* _MD1D_UTILS_H */


