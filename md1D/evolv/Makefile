#
#
#


NAME = md1D
STD_H = std_include.h
IN_OUT = $(NAME)_in_out
UTILS = $(NAME)_utils
PILOT = $(NAME)_pilot
ODE_SOLVE = ode_solve
EQ_DIFF = eq_diff

all : $(NAME)
#all : no_eq_diff

# Compilateur :
CC = gcc

# Options de compilation :
CFLAGS = -g -Wall -ansi
#CFLAGS = -O3 -Wall -ansi

# Chemin des libraries :
PATH1 = -I/usr/local/include

# Librairies à inclure :
LIB = -lfftw3 -lgsl -lgslcblas -lm
#LIB = -lfftw3 -lm

# Profilage
#PG = -pg
PG = 

# Compilation
$(NAME).o : $(NAME).c $(NAME).h $(STD_H)
	$(CC) -c $(CFLAGS) $(PATH1) $(NAME).c $(PG)

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
$(NAME) : $(NAME).o $(IN_OUT).o $(UTILS).o $(PILOT).o $(ODE_SOLVE).o $(EQ_DIFF).o
	$(CC) -o $(NAME) $(NAME).o $(IN_OUT).o $(UTILS).o $(PILOT).o $(ODE_SOLVE).o $(EQ_DIFF).o $(LIB) $(PG)

# Cas particulier : évite un conflit entre gsl et "complex.h" quand compilation avec g++ :
no_eq_diff : $(NAME).o $(IN_OUT).o $(UTILS).o $(PILOT).o $(ODE_SOLVE).o
	$(CC) -o $(NAME) $(NAME).o $(IN_OUT).o $(UTILS).o $(PILOT).o $(ODE_SOLVE).o $(LIB) $(PG)


# Archivage du projet en .tar
tar :  
#	tar -cvf $(NAME).tar ./*.c ./*.h ./Makefile ./md1D_param.txt ./utils/* ./TESTS/* 
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

