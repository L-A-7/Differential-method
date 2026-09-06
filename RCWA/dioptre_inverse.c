/*!	\file	reflexion_invert.c
 * 	Calculate the incident field that produce a given field after
 * 	reflexion on a plane interface
 *
 * 	\usage reflexion_invert args [options]		
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../md1D_utils.h"
#include "../md1D_io_utils.h"


#ifndef PI
#define PI 3.14159265358979323846
#endif /* PI */

int main(int argc, char *argv[])
{
	/* Reading command-line arguments */

	/* Loading the field to invert */
	re_AI = (double *) malloc(sizeof(double)*(2*N+1));
	im_AI = (double *) malloc(sizeof(double)*(2*N+1));
	AI   = (complex *) malloc(sizeof(complex)*(2*N+1));
	if (lire_tab(filename, "Re_Ai", re_Ai, 2*N+1) != 0 ){ /* Real part */
		fprintf(stderr,"%s line %d : ERROR, can't read Re_Ai in %s\n",__FILE__, __LINE__, filename);
		exit(EXIT_FAILURE);
	}
	if (lire_tab(filename, "Im_Ai", im_Ai, 2*N+1) != 0 ){ /* Imaginary part */
		fprintf(stderr,"%s line %d : ERROR, can't read Im_Ai in %s\n",__FILE__, __LINE__, filename);
		exit(EXIT_FAILURE);
	}
	for (i=-N;i<=N+1;i++){ /* Complex reconstitution */
		Ai[i+N] = re_Ai[i+N] + I*im_Ai[i+N];
	}
	free(re_Ai);
	free(im_Ai);

	
	
	
	
	
	
	
	
	
	}



#! /usr/bin/env python

# load system and math module:
import sys, os, math, myUtils       
#from scipy.optimize import *
from math import sqrt
from scipy import *

filename = "result.txt"
wavelength = 0.6328
L = 12.0
sigma_0 = 0.0
delta_sigma = 2.0*pi/L
nsub = 1.5 + 0.0j
nsuper = 1.0 + 0.0j
ksub   = 2*pi*nsub/wavelength
ksuper = 2*pi*nsuper/wavelength

# load real and imaginary components of the field to transform
re_Ai = myUtils.read_array(filename,"Re_Ai")
im_Ai = myUtils.read_array(filename,"Im_Ai")
print re_Ai
print im_Ai
if (len(re_Ai) != len(im_Ai)):
	print "ERROR : Re_Ai and Im_Ai don't have the same size. Exiting"
	sys.exit()
if (len(re_Ai)%2 != 1):
	print "ERROR : Re_Ai and Im_Ai must have 2N+1 elements. Exiting"
	sys.exit()
N=int( (len(re_Ai)-1.0) / 2.0 )	

print N

# Creating the reflexion matrix according to Fresnel laws
r_TE = mat(complex(zeros((2*N+1,2*N+1))))
for i in range(-N,N+1):
	sigma = i*delta_sigma
	print i
	b_sub   = (ksub**2   - sigma**2)**0.5
	print b_sub
	b_super = (ksuper**2 - sigma**2)**0.5
	print b_super
	r_TE[i+N,i+N] = (b_sub - b_super)
	#(b_sub - b_super)/(b_sub + b_super)

print r_TE

# 
# Calculate 



