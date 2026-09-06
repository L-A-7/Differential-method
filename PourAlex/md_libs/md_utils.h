#ifndef _MD2D_UTILS_H
#define _MD2D_UTILS_H


#include "std_include.h"
double md_chrono(struct Param_struct *par);
int SavePlot2file (double *, double *, int, char *);
int SaveDbleTab2file (double *tab, int N, char *filename, char *separateur, int Nmax1, char *separateur2);
int SaveCplxTab2file (COMPLEX *tab, int Nlign, char *mode, char *filename, char *separateur, int Nmax1, char *separateur2);
int SaveMatrix2file (COMPLEX **, int, int, char *, char *);
COMPLEX *CopyCplxTab(COMPLEX *, COMPLEX *, int);
COMPLEX **allocate_CplxMatrix(int, int);
COMPLEX ***allocate_CplxMatrix_3(int, int, int);
COMPLEX **reallocate_CplxMatrix(COMPLEX **tabl, int ncol,int nlign);
int CopyDbleTab(double *, double *, int);
double **allocate_DbleMatrix(int ncol,int nlign);
double ***allocate_DbleMatrix_3(int nlign, int ncol, int ntab);
#endif /* _MD2D_UTILS_H */


