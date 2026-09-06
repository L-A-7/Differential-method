NAME = abcde
CARTO = carto
OPTIM = optim
OPTIM_GSL = optim2
CC = gcc
CFLAGS = -O3 -Wall -ansi
LIB= -lm -L/opt/acml3.6.0/gnu64/lib -L/home/arnaud/lib
INCLUDE= -I/home/arnaud/include -I/opt/acml3.6.0/gnu64/include

#LIB_GSL= -lgsl -lgslcblas -lm 
MD2D = ../md2D
STD_H = std_include.h
IN_OUT = $(NAME)_in_out
UTILS = $(NAME)_utils
IO_UTILS = $(NAME)_io_utils
R_INDEX = refractive_index
SUB = sub_tanpsicosdelta

all : $(NAME)

#
$(NAME).o : $(NAME).c $(NAME).h
	$(CC) -c $(CFLAGS) $(PATH1) $(NAME).c

# Link the object files :
$(NAME) : $(NAME).o
	$(CC) -o $(NAME) $(NAME).o $(LIB)

tar :
	tar -cvf $(NAME).tar *.c *.h *.txt *.sh Makefile

# carto
$(CARTO) : $(CARTO).c $(OPTIM)_utils.o
	$(CC) -c $(CFLAGS) $(INCLUDE) $(CARTO).c
	$(CC) -o $(CARTO) $(CARTO).o $(OPTIM)_utils.o $(OPTIM)_io_utils.o $(LIB)


# optim :
$(OPTIM) : $(OPTIM).c $(OPTIM)_utils.o
	$(CC) -c $(CFLAGS) $(PATH1) $(OPTIM).c
	$(CC) -o $(OPTIM) $(OPTIM).o $(OPTIM)_utils.o $(OPTIM)_io_utils.o $(LIB) $(PATH2)

$(OPTIM_GSL) : $(OPTIM_GSL).c $(OPTIM)_utils.o
	$(CC) -c $(CFLAGS) $(PATH1) $(OPTIM_GSL).c
	$(CC) -o $(OPTIM_GSL) $(OPTIM_GSL).o $(OPTIM)_utils.o $(OPTIM)_io_utils.o $(LIB_GSL) $(PATH2)

$(OPTIM)_utils.o : $(MD2D)_utils.c $(MD2D)_io_utils.c
	$(CC) -c -o $(OPTIM)_utils.o $(CFLAGS)  $(INCLUDE) $(MD2D)_utils.c 
	$(CC) -c -o $(OPTIM)_io_utils.o $(CFLAGS) $(INCLUDE) $(MD2D)_io_utils.c

# refractive_index
$(R_INDEX) : $(R_INDEX).c $(OPTIM)_utils.o
	$(CC) -c $(CFLAGS) $(PATH1) $(R_INDEX).c
	$(CC) -o $(R_INDEX) $(R_INDEX).o $(OPTIM)_utils.o $(OPTIM)_io_utils.o $(LIB) $(PATH2)

# sub_tanpsicosdelta
$(SUB) : $(SUB).c $(OPTIM)_utils.o
	$(CC) -c $(CFLAGS) $(PATH1) $(SUB).c
	$(CC) -o $(SUB) $(SUB).o $(OPTIM)_utils.o $(OPTIM)_io_utils.o $(LIB) $(PATH2)

