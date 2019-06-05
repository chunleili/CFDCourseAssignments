vpath %.H include
vpath %.C src
CC=g++
SRC = ${wildcard src/*.C}
OBJ = ${patsubst %.C, %.o, $(SRC)}
CFLAGS = -g -Wall -Iincludes 

all : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o exe
%.o: %.C
	$(CC) $(CFLAGS) -c  $< -o $@
move :
	mv $(OBJ)  ./objs
clean:
	-rm $(OBJ) objs/*.o data/*.dat exe

.PHONY: clean move
