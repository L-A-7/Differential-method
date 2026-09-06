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
PROFIL_GEN = profilGen
RK4_TEST = rk4_test
ODESOLVE = ode_solve
TESTS = tests
FFT_TEST = fft_test

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



# ProfilGen (indépendant) :
$(PROFIL_GEN) : $(PROFIL_GEN).c $(STD_H)
	$(CC) -c $(CFLAGS) $(PATH1) $(PROFIL_GEN).c
	$(CC) -o $(PROFIL_GEN) $(PROFIL_GEN).o $(LIB)
	
# rk4_test (indépendant)
$(RK4_TEST) : $(RK4_TEST).c $(ODESOLVE).c
	$(CC) -o $(RK4_TEST) $(CFLAGS) $(RK4_TEST).c $(ODESOLVE).c -lm

# fft_test (indépendant)
$(FFT_TEST) : $(FFT_TEST).c 
	$(CC) -o $(FFT_TEST) $(CFLAGS) $(FFT_TEST).c -lm -lfftw3
	
# lire_donnees (indépendant)
lire_donnees : lire_donnees.c 
	$(CC) -o lire_donnees $(CFLAGS) lire_donnees.c

# lire_plots (indépendant)
lire_plot : lire_plot.c 
	$(CC) -o lire_plot $(CFLAGS) lire_plot.c

# lire_string (indépendant)
lire_string : lire_string.c 
	$(CC) -o lire_string $(CFLAGS) lire_string.c
	
#tests (indépendant)
$(TESTS) : $(TESTS).c 
	$(CC) -o $(TESTS) $(CFLAGS) $(TESTS).c md1D_utils.c -lm -lfftw3
	
clean :
	rm -f ./*.o
	rm -f ./*~

