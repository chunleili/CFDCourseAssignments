#SRCS=assignment4.C
#OBJS:=$(wildcard *.o)
OBJS=assignment4.o solve.o mesh.o utility.o Roe.o 
CFLAG=-g -O2 -Wall -std=c++11 
all:$(OBJS) 
	g++ $(CFLAG) -lm $(OBJS)

assignment4.o: assignment4.C main.H
	g++ -c $(CFLAG) assignment4.C 

solve.o: solve.C mesh.C main.H
	g++ -c $(CFLAG) solve.C

mesh.o: mesh.C main.H
	g++ -c $(CFLAG) mesh.C

utility.o: utility.C main.H
	g++ -c $(CFLAG) utility.C

Roe.o: Roe.C main.H
	g++ -c $(CFLAG) utility.C

clean:
	rm -rf *.o a.out *.dat
