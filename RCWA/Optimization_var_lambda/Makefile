NAME = optimiz
OPTIM = optim
CC = gcc
CFLAGS = -O3 -Wall -ansi
LIB= liblevmar.a -lm
MD2D = ../md2D
STD_H = std_include.h
IN_OUT = $(NAME)_in_out
UTILS = $(NAME)_utils
IO_UTILS = $(NAME)_io_utils



all : $(NAME)

#
$(NAME).o : $(NAME).c $(NAME).h
	$(CC) -c $(CFLAGS) $(PATH1) $(NAME).c

# Link the object files :
$(NAME) : $(NAME).o
	$(CC) -o $(NAME) $(NAME).o $(LIB)

tar :
	tar -cvf $(NAME).tar *.c *.h *.txt

# carto
$(CARTO) : carto.c $(NAME).h
	$(CC) -c $(CFLAGS) $(PATH1) carto.c
	$(CC) -o carto carto.o $(LIBS)


# optim :
$(OPTIM) : $(OPTIM).c $(OPTIM)_utils.o
	$(CC) -c $(CFLAGS) $(PATH1) $(OPTIM).c ../md2D_utils.c ../md2D_io_utils.c
	$(CC) -o $(OPTIM) $(OPTIM).o $(OPTIM)_utils.o $(OPTIM)_io_utils.o $(LIB) $(PATH2)

$(OPTIM)_utils.o : $(MD2D)_utils.c $(MD2D)_io_utils.c
	$(CC) -c -o $(OPTIM)_utils.o $(CFLAGS) $(PATH1) $(MD2D)_utils.c 
	$(CC) -c -o $(OPTIM)_io_utils.o $(CFLAGS) $(PATH1) $(MD2D)_io_utils.c


