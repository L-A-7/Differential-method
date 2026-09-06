NAME = test_optim2
CC = gcc
CFLAGS = -O3 -Wall -ansi
LIBS= liblevmar.a -lm

all : $(NAME)

#
$(NAME).o : $(NAME).c
	$(CC) -c $(CFLAGS) $(PATH1) $(NAME).c

# Link the object files :
$(NAME) : $(NAME).o
	$(CC) -o $(NAME) $(NAME).o $(LIBS)

tar :
	tar -cvf $(NAME).tar *.c *.h Makefile *.txt
