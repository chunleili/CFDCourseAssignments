#SRCS=assignment4.C
OBJS=assignment4.o solve.o
CFLAG=-g -O2 -Wall -std=c++11 
all:$(OBJS) 
	g++ $(CFLAG) -lm $(OBJS)
assignment4.o: assignment4.C main.H
	g++ -c assignment4.C 
solve.o: solve.C main.H
	g++ -c solve.C