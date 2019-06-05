vpath %.H include
vpath %.C src
CC=g++
SRC = ${wildcard src/*.C}
OBJ = ${patsubst %.C, %.o, $(SRC)}
CFLAGS = -g -Wall -Iinclude

all : $(OBJ)
	 $(CC) $(CFLAGS) $(OBJ) -o main
	 #mv $(OBJ)  ./objs

$(OBJ):%.o : %.C
	$(CC) $(CFLAGS) -c $^ -o $@ 

move :
	mv $(OBJ)  ./objs

clean:
	rm $(OBJ) objs/*.o data/*.dat main 

.PHONY:all clean move
