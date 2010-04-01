/* \file  md2D_io_utils.h
 *  \brief Fichier d'en-tête pour md2D_io_utils.c
 */

#ifndef _md2D_I0_UTILS_H
#define _md2D_I0_UTILS_H

#include "std_include.h"

/* Constantes */
#define CHAR_COMMENT '#'

/* functions*/
int lire_int(FILE *fp, const char *label, int *value);
int lire_double(FILE *fp, const char *label, double *value);
int lire_string(FILE *fp, const char *label, char *value);
int lire_COMPLEX(FILE *fp, const char *label, COMPLEX *value);
int ecrire_dble_tab(FILE *fp, double *tab, int N, char *separateur1, int Nmax1, char *separateur2);
int ecrire_cplx_tab(FILE *fp, COMPLEX *tab, int N, int mode, char *separateur1, int Nmax1, char *separateur2);
char *label_search(char *str,const char *label);
int lire_tab(const char *filename, const char *label, double *tab, int N);
int lire_ligne(FILE *fp, char *line);
void skip_comment(char *string);
int lire_str_arg(char *dest, char *label, int argc, char **argvcp);
int lire_dble_arg(double *res, char *label, int argc, char **argvcp);
int lire_int_arg(int *res, char *label, int argc, char **argvcp);

#endif /* _md2D_IO_UTILS_H */

