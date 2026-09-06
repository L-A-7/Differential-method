#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <complex.h>
#include "lm.h"
#define SIZE_LINE_BUFFER 50000
#define SIZE_STR_BUFFER 5000
#define CHAR_COMMENT '#'
#define ROUND(x) ((int)(x<0 ? x-0.5 : x+0.5))
#define CEIL(x)  (x-(int)(x)>0 ? (int)(x)+1 : (int)(x))
#define MIN(a,b) ((a<b)?a:b)

int tanPsi_cosDelta_lambda(double *tanPsicosDelta, double *var_lambda, int Nb_lambda, double *param);
void skip_comment(char *str_in_out);
char *label_search(char *str,const char *label);
int lire_ligne(FILE *fp, char *line);
int lire_tab(const char *nom_fichier, const char *label, double *tab, int N);
void func(double *p, double *hx, int m, int n, void *adata);
int ecrire_dble_tab(FILE *fp, double *tab, int N, char *separateur1, int Nmax1, char *separateur2);
int SaveDbleTab2file (double *tab, int N, char *filename, char *separateur);
complex refractive_index(char* name,double lambda, char* method);

