#
#
#
ACML_INCLUDE_DIR = /opt/acml4.4.0/gfortran32/include 
ACML_LIB_DIR = /opt/acml4.4.0/gfortran32/lib
ATLAS_INCLUDE_DIR = /usr/include/atlas
ATLAS_LIB_DIR = /usr/lib64/atlas

#BUILD=ACML
#BUILD=Lapack
#BUILD=Lapack_ATLAS
HAVE_ACML=1
HAVE_BLAS=0
HAVE_CBLAS=0
HAVE_LAPACK=0
HAVE_ATLAS=0


NAME = md3D
STD_H = std_include.h
IN_OUT = $(NAME)_in_out
UTILS = md_utils
IO_UTILS = md_io_utils
PILOT = $(NAME)_pilot
MATHS = md_maths

LIB1 = -lfftw3

ifeq ($(HAVE_ACML),1)
INC2= -I$(ACML_INCLUDE_DIR)
#LIB2 = -L$(ACML_LIB_DIR) -lacml
LIB2 = $(ACML_LIB_DIR)/libacml.a
OPTIM_FLAGS= -D_ACML
endif

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
LIB =$(LIB1) $(LIB2) -lm -lgfortran

all : $(NAME)

# Options de compilation :
CFLAGS = -g -Wall -ansi $(OPTIM_FLAGS)
#CFLAGS = -O3 -Wall -ansi $(OPTIM_FLAGS)

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
$(PILOT).o : $(PILOT).c $(NAME).h $(STD_H)
	$(CC) -c $(CFLAGS) $(INC) $(PILOT).c $(PG)

# Link the object files :
$(NAME) : $(NAME).o $(IN_OUT).o $(IO_UTILS).o $(UTILS).o $(PILOT).o $(MATHS).o
	$(CC) -o $(NAME) $(NAME).o $(IO_UTILS).o $(IN_OUT).o $(UTILS).o $(PILOT).o $(MATHS).o $(LIB) $(PG)

# Archivage du projet en .tar
tar :  
#	tar -cvf $(NAME).tar ./*.c ./*.h  ./Makefile ./*.txt
	tar -cvf $(NAME).tar ./*.c ./*.h ./*.txt

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

