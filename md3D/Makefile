#
#
#
ACML_INCLUDE_DIR = /opt/acml4.3.0/gfortran64/include 
ACML_LIB_DIR = /opt/acml4.3.0/gfortran64/lib
ATLAS_INCLUDE_DIR = /usr/include/atlas
ATLAS_LIB_DIR = /usr/lib64/atlas

#BUILD=ACML
#BUILD=Lapack
BUILD=Lapack_ATLAS

NAME = md3D
STD_H = std_include.h
IN_OUT = $(NAME)_in_out
UTILS = $(NAME)_utils
IO_UTILS = $(NAME)_io_utils
PILOT = $(NAME)_pilot


#INC1 = -I/home/lau/Programmes/md1D/fftw
#LIB1 = -static -lfftw3 -lgsl 
#LIB1 = -static -lfftw3 -lgsl -lgslcblas
LIB1 = -lfftw3 -lgsl -lgslcblas

ifeq ($(BUILD),Lapack)
INC2 = -I$(ATLAS_INCLUDE_DIR)
LIB2 = -llapack -lblas 
XERBLA=xerbla.o
OPTIM_FLAGS=-D_LAPACK
endif

ifeq ($(BUILD),Lapack_ATLAS)
INC2 = -I$(ATLAS_INCLUDE_DIR)
LIB2 = -L$(ATLAS_LIB_DIR) -latlas -llapack -lblas 
XERBLA=xerbla.o
OPTIM_FLAGS=-D_BLAS -D_LAPACK
endif

ifeq ($(BUILD),ACML)
INC2= -I$(ACML_INCLUDE_DIR)
LIB2 = -L$(ACML_LIB_DIR) -lacml
OPTIM_FLAGS= -D_ACML
endif

ifeq ($(COMPIL),ICC)
CC=/opt/intel/cce/9.1.047/bin/icc
LIB3 = -L/opt/intel/cce/9.1.047/lib
else
CC=gcc
endif

ifeq ($(SCHED),FIFO)
OPTIM_FLAGS+=-DSCHED_FIFO=1
endif

INC =$(INC1) $(INC2) $(INC3)
#LIB =$(LIB1) $(LIB2) -lf2c -lm
LIB =$(LIB1) $(LIB2) -lm

all : $(NAME)

# Options de compilation :
#CFLAGS = -g -Wall -ansi $(OPTIM_FLAGS)
CFLAGS = -O3 -Wall -ansi $(OPTIM_FLAGS)

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
$(PILOT).o : $(PILOT).c $(PILOT).h $(STD_H)
	$(CC) -c $(CFLAGS) $(INC) $(PILOT).c $(PG)

# Link the object files :
$(NAME) : $(NAME).o $(IN_OUT).o $(IO_UTILS).o $(UTILS).o $(PILOT).o
	$(CC) -o $(NAME) $(NAME).o $(IO_UTILS).o $(IN_OUT).o $(UTILS).o $(PILOT).o $(LIB) $(PG)

# Archivage du projet en .tar
tar :  
	tar -cvf $(NAME).tar ./*.c ./*.h  ./Makefile ./*.txt

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

