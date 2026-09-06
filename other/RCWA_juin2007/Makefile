#
#
#


NAME = md2D
STD_H = std_include.h
IN_OUT = $(NAME)_in_out
UTILS = $(NAME)_utils
IO_UTILS = $(NAME)_io_utils
PILOT = $(NAME)_pilot
ODE_SOLVE = ode_solve
EQ_DIFF = eq_diff

# Uncomment the following depending on the available libraries
#BLAS_FLAG = -D_BLAS_OPTIMIZATION
#LAPACK_FLAG = -D_LAPACK_OPTIMIZATION
ACML_FLAG = -D_ACML_OPTIMIZATION

#all : $(NAME)
all : no_eq_diff

# Compilateur :
CC = gcc

# Options pour le precompilateur
PREPROC = $(BLAS_FLAG) $(LAPACK_FLAG) $(ACML_FLAG)

# Options de compilation :
#CFLAGS = -g -Wall -ansi
CFLAGS = -O3 -Wall -ansi $(PREPROC)

# Chemin des libraries :
PATH1 = -I/usr/local/include -I/home/lau/Programmes/md1D/fftw -I/opt/acml3.6.0/gnu64/include
PATH2 = -L/opt/acml3.6.0/gnu64/lib

# Librairies à inclure :
#LIB = -lfftw3 -lgsl -lgslcblas -llapack -lblas -lg2c -lm
LIB = -lfftw3 -llapack -lblas -lacml -lg2c -lm
#LIB = -lfftw3 -lacml -lg2c -lm
#LIB = -lfftw3 -lm

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

# Link the object files :
$(NAME) : $(NAME).o $(IN_OUT).o $(IO_UTILS).o $(UTILS).o $(PILOT).o $(ODE_SOLVE).o $(EQ_DIFF).o
	$(CC) -o $(NAME) $(NAME).o $(IO_UTILS).o $(IN_OUT).o $(UTILS).o $(PILOT).o $(ODE_SOLVE).o $(EQ_DIFF).o $(LIB) $(PG)

# Cas particulier : évite un conflit entre gsl et "complex.h" quand compilation avec g++ :
no_eq_diff : $(NAME).o $(IN_OUT).o  $(IO_UTILS).o $(UTILS).o $(PILOT).o $(ODE_SOLVE).o
	$(CC) -o $(NAME) $(NAME).o $(IO_UTILS).o $(IN_OUT).o $(UTILS).o $(PILOT).o $(ODE_SOLVE).o $(LIB) $(PG) $(PATH2)


# Archivage du projet en .tar
tar :  
#	tar -cvf $(NAME).tar ./*.c ./*.h ./Makefile ./md1D_param.txt ./utils/* ./TESTS/* ./results/*
	tar -cvf $(NAME).tar ./*.c ./*.h ./*.txt ./utils/* ./Makefile
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

