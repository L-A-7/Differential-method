#
#
#


NAME = md1D
STD_H = std_include.h
NAME2 = $(NAME)_in_out
NAME3 = $(NAME)_utils
NAME4 = $(NAME)_pilot
NAME5 = ode_solve
NAME6 = eq_diff

all : $(NAME)

# Compilateur :
CC = gcc

# Options de compilation :
#CFLAGS = -g -Wall -ansi
CFLAGS = -O3 -Wall -ansi

# Chemin des libraries :
#PATH1 = -L/usr/local/lib -I/usr/local/include

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
$(NAME2).o : $(NAME2).c $(NAME2).h $(STD_H)
	$(CC) -c $(CFLAGS) $(PATH1) $(NAME2).c $(PG) 

#
$(NAME3).o : $(NAME3).c $(NAME3).h $(STD_H)
	$(CC) -c $(CFLAGS) $(PATH1) $(NAME3).c $(PG)

#
$(NAME4).o : $(NAME4).c $(NAME4).h $(STD_H)
	$(CC) -c $(CFLAGS) $(PATH1) $(NAME4).c $(PG)

#
$(NAME5).o : $(NAME5).c $(STD_H)
	$(CC) -c $(CFLAGS) $(PATH1) $(NAME5).c $(PG)

#
$(NAME6).o : $(NAME6).c $(STD_H)
	$(CC) -c $(CFLAGS) $(PATH1) $(NAME6).c $(PG)

# Link the object files :
$(NAME) : $(NAME).o $(NAME2).o $(NAME3).o $(NAME4).o $(NAME5).o $(NAME6).o
	$(CC) -o $(NAME) $(NAME).o $(NAME2).o $(NAME3).o $(NAME4).o $(NAME5).o $(NAME6).o $(LIB) $(PG)

# Cas particulier : évite un conflit entre gsl et "complex.h" quand compilation avec g++ :
no_eq_diff : $(NAME).o $(NAME2).o $(NAME3).o $(NAME4).o $(NAME5).o
	$(CC) -o $(NAME) $(NAME).o $(NAME2).o $(NAME3).o $(NAME4).o $(NAME5).o $(LIB) $(PG)


# Compression du projet en tar.bz2
bz2 :
	tar -cvf $(NAME).tar ./*.c ./*.h ./Makefile ./md1D_param.txt ./utils/* ./TESTS/*
	bzip2 $(NAME).tar
	

# Nettoyage du répertoire de travail	
clean :
	rm -f ./*.o
	rm -f ./*~

