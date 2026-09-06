#
#
#
#ACML_INCLUDE_DIR = /opt/acml3.6.0/gnu64/include 
#ACML_LIB_DIR = /opt/acml3.6.0/gnu64/lib 
ACML_INCLUDE_DIR = /opt/acml3.6.0/gnu32/include 
ACML_LIB_DIR = /opt/acml3.6.0/gnu32/lib 
ATLAS_INCLUDE_DIR = /usr/include/atlas
ATLAS_LIB_DIR = /usr/lib64/atlas

BUILD=Lapack_ATLAS

ifeq ($(BUILD),Original)
CC=gcc
CFLAGS = -O3 -Wall -ansi $(OPTIM_FLAGS)
PATH1 = -I/usr/local/include -I/home/lau/Programmes/md1D/fftw -I$(ACML_INCLUDE_DIR)
PATH2 = -L$(ACML_LIB_DIR)
LIB = -lfftw3 -llapack -lblas -lacml_mv -lacml -lf2c -lm
OPTIM_FLAGS=-DBLAS_OPTIMIZATION=0 -DLAPACK_OPTIMIZATION=0 -DACML_OPTIMIZATION=1
else
ifeq ($(BUILD),Lapack_ATLAS)
CC=gcc
CFLAGS = -O3 -Wall -ansi $(OPTIM_FLAGS)
PATH1 = -I/usr/local/include -I$(ATLAS_INCLUDE_DIR)
# -static: to simplify profiling (oprofile, gprof)
LIB = -lfftw3 -L$(ATLAS_LIB_DIR) -latlas -llapack -lblas -lm 
XERBLA=xerbla.o
OPTIM_FLAGS=-DBLAS_OPTIMIZATION=0 -DLAPACK_OPTIMIZATION=1 -DACML_OPTIMIZATION=0
else
ifeq ($(BUILD),ACML)
CC=gcc
CFLAGS = -O3 -Wall -ansi $(OPTIM_FLAGS)
# -static: to simplify profiling (oprofile, gprof)
PATH1=-I$(ACML_INCLUDE_DIR)
#LIB = -static -lfftw3 -lm -L$(ACML_LIB_DIR) -lacml -lacml_mv -lg2c
LIB = -lfftw3 -L$(ACML_LIB_DIR) -lacml -lm
OPTIM_FLAGS=-DBLAS_OPTIMIZATION=0 -DLAPACK_OPTIMIZATION=0 -DACML_OPTIMIZATION=1
else
ifeq ($(BUILD),ICC)
CC=/opt/intel/cce/9.1.047/bin/icc
CFLAGS = -O3 -Wall -ansi $(OPTIM_FLAGS)
# -static: to simplify profiling (oprofile, gprof)
PATH1=-I$(ACML_INCLUDE_DIR)
LIB = -lfftw3 -lm -L$(ACML_LIB_DIR) -lacml -lacml_mv -L/opt/intel/cce/9.1.047/lib
OPTIM_FLAGS:=-DBLAS_OPTIMIZATION=0 -DLAPACK_OPTIMIZATION=0 -DACML_OPTIMIZATION=1
else
default: help
endif
endif
endif
endif

ifeq ($(SCHED),FIFO)
OPTIM_FLAGS+=-DSCHED_FIFO=1
endif

NAME = md2D
STD_H = std_include.h
IN_OUT = $(NAME)_in_out
UTILS = $(NAME)_utils
IO_UTILS = $(NAME)_io_utils
PILOT = $(NAME)_pilot
ODE_SOLVE = ode_solve
EQ_DIFF = eq_diff

#all : $(NAME)
all : no_eq_diff

help:
	@echo "To compile choose the mathematical library: "
	@echo "  BUILD=Original make clean all"
	@echo "  BUILD=Lapack_ATLAS make clean all"
	@echo "  BUILD=ACML make clean all"
	@echo "  BUILD=ICC make clean all"
	@echo
	@echo "  Add SCHED=FIFO to give max priority to the md2D processus."



# Profilage
#PG = -pg
#PG = 

# Compilation
$(NAME).o : $(NAME).c $(NAME).h $(STD_H)
	$(CC) -c $(CFLAGS) $(PATH1) $(NAME).c $(PG)

#
$(IO_UTILS).o : $(IO_UTILS).c $(IO_UTILS).h $(STD_H)
	$(CC) -c $(CFLAGS) $(PATH1) $(IO_UTILS).c $(PG)

#
$(IN_OUT).o : $(IN_OUT).c $(IN_OUT).h $(STD_H)
	$(CC) -c $(CFLAGS) $(PATH1) $(IN_OUT).c $(PG) 

#
$(UTILS).o : $(UTILS).c $(UTILS).h $(STD_H)
	$(CC) -c $(CFLAGS) $(PATH1) $(UTILS).c $(PG)

#
$(PILOT).o : $(PILOT).c $(PILOT).h $(STD_H)
	$(CC) -c $(CFLAGS) $(PATH1) $(PILOT).c $(PG)

#
$(ODE_SOLVE).o : $(ODE_SOLVE).c $(STD_H)
	$(CC) -c $(CFLAGS) $(PATH1) $(ODE_SOLVE).c $(PG)

#
$(EQ_DIFF).o : $(EQ_DIFF).c $(STD_H)
	$(CC) -c $(CFLAGS) $(PATH1) $(EQ_DIFF).c $(PG)

#xerbla.o: xerbla.f
#	g77 -O3 -o xerbla.o -c xerbla.f

# Link the object files :
$(NAME) : $(NAME).o $(IN_OUT).o $(IO_UTILS).o $(UTILS).o $(PILOT).o $(ODE_SOLVE).o $(EQ_DIFF).o
	$(CC) -o $(NAME) $(NAME).o $(IO_UTILS).o $(IN_OUT).o $(UTILS).o $(PILOT).o $(ODE_SOLVE).o $(EQ_DIFF).o $(LIB) $(PG)

# Cas particulier : évite un conflit entre gsl et "complex.h" quand compilation avec g++ :
no_eq_diff : $(NAME).o $(IN_OUT).o  $(IO_UTILS).o $(UTILS).o $(PILOT).o $(ODE_SOLVE).o
	$(CC) -o $(NAME) $(NAME).o $(IO_UTILS).o $(IN_OUT).o $(UTILS).o $(PILOT).o $(ODE_SOLVE).o $(LIB) $(PG) $(PATH2)


# Archivage du projet en .tar
tar :  
#	tar -cvf $(NAME).tar ./*.c ./*.h ./Makefile ./md1D_param.txt ./utils/* ./TESTS/* ./results/*
	tar -cvf $(NAME).tar ./*.c ./*.h ./*param* ./Makefile
#	tar -cvf $(NAME).tar ./*.c ./*.h ./*.m ./md1D_param.txt ./utils/*  
#	tar -cvf $(NAME).tar ./*.c ./*.h  

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

