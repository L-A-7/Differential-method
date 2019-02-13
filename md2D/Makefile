#
#
#

#HOME = /home/arnaud
ATLAS_INCLUDE_DIR = /usr/include/atlas
ATLAS_LIB_DIR = /usr/lib/atlas-base
LAPACK_INCLUDE_DIR = /usr/include/atlas
ICC_DIR = /opt/intel/compilers_and_libraries/linux/bin/intel64/

MKL = /opt/intel/mkl
MKL_LIB = $(MKL)/lib/em64t
MKL_INC = $(MKL)/include
#FFTW_PATH = /opt/scilibs/FFTW/fftw-3.2.1_Bull.9010
#FFTW_LIB = $(FFTW_PATH)/lib
#FFTW_INC = $(FFTW_PATH)/include
#LOCAL_INC = $(HOME)/include
#LOCAL_LIB = $(HOME)/lib

HAVE_MKL=0
HAVE_BLAS=1
HAVE_CBLAS=0
HAVE_LAPACK=1
HAVE_ATLAS=0
#COMPIL=ICC
COMPIL=GCC

NAME = md2D
STD_H = std_include.h
IN_OUT = $(NAME)_in_out
UTILS = md_utils
IO_UTILS = md_io_utils
PILOT = $(NAME)_pilot
MATHS = md_maths
ODE = md_odesolve

######## FFTW ###########
#LIB1 = -static -lfftw3
#LIB= -L$(LOCAL_LIB) $(LOCAL_LIB)/libfftw3.a
#INC= -I$(LOCAL_INC) -I$(FFTW_INC)
LIB= -lfftw3
DEFS=

ifeq ($(HAVE_MKL),1)
INC+= -I$(MKL_INC) -I$(LAPACK_INCLUDE_DIR)
#LIB+= -L$(MKL_LIB) $(MKL_LIB)/libmkl_solver_lp64.a -Wl,--start-group $(MKL_LIB)/libmkl_intel_lp64.a $(MKL_LIB)/libmkl_intel_thread.a $(MKL_LIB)/libmkl_core.a -Wl,--end-group -openmp -lpthread
LIB+=  -Wl,--start-group /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/mkl/lib/intel64/libmkl_gnu_thread.a /opt/intel/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
#LIB+= -mkl
DEFS+= -D_BLAS -D_LAPACK  -m64

endif

ifeq ($(HAVE_BLAS),1)
#INC+= -I$(BLAS_INCLUDE_DIR)
#LIB+= -L$(BLAS_LIB_DIR) -lblas
LIB+= -lblas
DEFS+= -D_BLAS
endif

ifeq ($(HAVE_CBLAS),1)
#INC+= -I$(BLAS_INCLUDE_DIR)
#LIB+= -L$(BLAS_LIB_DIR) -lcblas
LIB+= -lcblas
DEFS+= -D_CBLAS
endif

ifeq ($(HAVE_LAPACK),1)
INC+= -I$(LAPACK_INCLUDE_DIR)
#LIB+= -L$(LAPACK_LIB_DIR) -llapack
LIB+= -llapack
DEFS+= -D_LAPACK
endif

ifeq ($(HAVE_ATLAS),1)
INC+= -I$(ATLAS_INCLUDE_DIR)
LIB+= -L$(ATLAS_LIB_DIR) -latlas -llapack -lblas 
XERBLA=xerbla.o
DEFS+= -D_BLAS -D_LAPACK
endif

ifeq ($(COMPIL),ICC)
#CC=$(ICC_DIR)/icc -wd981
CC = $(ICC_DIR)/icc
#LIB3 = -L/opt/intel/cce/9.1.047/lib
else
CC=gcc
endif

LIB+= -lm

all : $(NAME)

# Options de compilation :
#CFLAGS = -g -Wall -ansi $(DEFS)
CFLAGS = -O3 -ansi $(DEFS)


# Profiling
#PG = -pg
#PG = 

# Compilation
$(NAME).o : $(NAME).c $(NAME).h $(STD_H)
	$(CC) -c $(CFLAGS) $(INC) $(NAME).c $(PG)

#
$(IO_UTILS).o : $(IO_UTILS).c $(IO_UTILS).h $(STD_H)
	$(CC) -c $(CFLAGS) $(INC) $(IO_UTILS).c $(PG)

#
$(IN_OUT).o : $(IN_OUT).c $(IN_OUT).h $(STD_H)
	$(CC) -c $(CFLAGS) $(INC) $(IN_OUT).c $(PG) 

#
$(UTILS).o : $(UTILS).c $(UTILS).h $(STD_H)
	$(CC) -c $(CFLAGS) $(INC) $(UTILS).c $(PG)

#
$(MATHS).o : $(MATHS).c $(MATHS).h $(STD_H)
	$(CC) -c $(CFLAGS) $(INC) $(MATHS).c $(PG)

#
$(ODE).o : $(ODE).c $(ODE).h $(STD_H)
	$(CC) -c $(CFLAGS) $(INC) $(ODE).c $(PG)

#
$(PILOT).o : $(PILOT).c $(NAME).h $(STD_H)
	$(CC) -c $(CFLAGS) $(INC) $(PILOT).c $(PG)

# Link the object files :
$(NAME) : $(NAME).o $(IN_OUT).o $(IO_UTILS).o $(UTILS).o $(PILOT).o $(MATHS).o $(ODE).o
	$(CC) -o $(NAME) $(NAME).o $(IO_UTILS).o $(IN_OUT).o $(UTILS).o $(PILOT).o $(MATHS).o $(ODE).o $(LIB) $(PG)

# Archivage du projet en .tar
tar :  
#	tar -cvf $(NAME).tar ./*.c ./*.h  ./Makefile ./*.txt
	tar -cvf $(NAME).tar ./*.c ./*.h 

# Archivage et Compression du projet en .tar.bz2
bz2 : tar
	rm -f $(NAME).tar.bz2
	bzip2 $(NAME).tar

# Archivage et Compression du projet en .tar.gz
targz : tar
	rm -f $(NAME).tar.gz
	gzip $(NAME).tar

# Archivage et Compression du projet en .zip
zip :
	zip $(NAME) ./*.c ./*.h ./Makefile ./md1D_param.txt ./utils/* ./TESTS/* ./doc/*


# Nettoyage du répertoire de travail	
clean :
	rm -f ./*.o
	rm -f ./*~



