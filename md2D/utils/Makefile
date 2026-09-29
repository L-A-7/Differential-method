#
#
#


NAME = md2D
STD_H = std_include.h
NAME2 = $(NAME)_io_utils
NAME3 = $(NAME)_utils
NAME4 = $(NAME)_pilot
NAME5 = ode_solve
NAME6 = eq_diff
PROFIL_GEN = profilGen
RK4_TEST = rk4_test
ODESOLVE = ode_solve
TESTS = test2
FFT_TEST = fft_test

# Compilateur :
CC = gcc

# Options de compilation :
CFLAGS = -g -Wall -ansi
#CFLAGS = -O3 -Wall -ansi

# Chemin des libraries :
PATH1 = -I/usr/local/include -I/home/lau/Programmes/md1D/fftw -I/opt/acml3.6.0/gnu64/include
PATH2 = -L/opt/acml3.6.0/gnu64/lib

# Librairies à inclure :
#LIB = -lfftw3 -lgsl -lgslcblas -lm
#LIB = -lm
#LIB = -lfftw3 -llapack -lblas -lacml -lg2c -lm
LIB = -lfftw3 -llapack -lblas -lacml -lm

# Profilage
#PG = -pg
PG = 

utils : $(PROFIL_GEN) $(FFT_TEST) lire_tab lire_ligne lire_N1N2 add_col lire_plot lire_string diffsum $(TESTS)

# ProfilGen (indépendant) :
$(PROFIL_GEN) : $(PROFIL_GEN).c
	$(CC) -c $(CFLAGS) $(PATH1) $(PROFIL_GEN).c ../md2D_utils.c ../md2D_io_utils.c
	$(CC) -o $(PROFIL_GEN) $(PROFIL_GEN).o ../md2D_utils.o ../md2D_io_utils.o $(LIB) $(PATH2)
	
# fft_test (indépendant)
$(FFT_TEST) : $(FFT_TEST).c 
	$(CC) -o $(FFT_TEST) $(CFLAGS) $(FFT_TEST).c -lm -lfftw3
	
# lire_tab (indépendant)
lire_tab : lire_tab.c 
	$(CC) -o lire_tab $(CFLAGS) lire_tab.c

# lire_tab (indépendant)
lire_tab2 : lire_tab2.c 
	$(CC) -o lire_tab2 $(CFLAGS) lire_tab2.c

# lire_ligne (indépendant)
lire_ligne : lire_ligne.c 
	$(CC) -o lire_ligne $(CFLAGS) lire_ligne.c

# lire_N1N2 (indépendant)
lire_N1N2 : lire_N1N2.c 
	$(CC) -o lire_N1N2 $(CFLAGS) lire_N1N2.c

# add_col (indépendant)
add_col : add_col.c 
	$(CC) -o add_col $(CFLAGS) add_col.c

# lire_plots (indépendant)
lire_plot : lire_plot.c 
	$(CC) -o lire_plot $(CFLAGS) lire_plot.c

# lire_string (indépendant)
lire_string : lire_string.c 
	$(CC) -o lire_string $(CFLAGS) lire_string.c

# diffsum (indépendant)
diffsum : diffsum.c 
	$(CC) -o diffsum $(CFLAGS) diffsum.c -lm
	
#tests (indépendant)
$(TESTS) : $(TESTS).c 
	$(CC) -o $(TESTS) $(CFLAGS) $(TESTS).c -lm -lfftw3
	
clean :
	rm -f ./*.o
	rm -f ./*~

