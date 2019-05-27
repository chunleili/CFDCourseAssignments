SOL  = scalarJSTSolver
#SOL  = MacCormackSolver
OBJS = init.o print.o ${SOL}.o main.o UToF.o CFL.o
SRCS = init.C print.C ${SOL}.C main.C UToF.C CFL.C
LIBS = -lm
EXE  = out
CFLG = -g -Wall -std=c++11	
start   : $(OBJS) main.H
	g++ -o ${EXE} ${CFLG} ${OBJS} ${LIBS}

CFL.o :	CFL.C main.H
	g++ ${CFLG} -c CFL.C

print.o : print.C main.H
	g++  ${CFLG} -c print.C

init.o : init.C main.H
	g++  ${CFLG} -c init.C

UToF.o : UToF.C main.H
	g++  ${CFLG} -c UToF.C

${SOL}.o  : ${SOL}.C UToF.C main.H 
	g++  ${CFLG} -c ${SOL}.C

main.o    : $(SRCS) main.H 
	g++  ${CFLG} -c main.C

.PHONY:clean  clean_data clean_log clean_all
clean      :
	rm -rf *.o 
clean_data :
	rm -rf data/*
clean_log  :
	rm -rf log/*
clean_all  :
	rm -rf data/* log/* *.o ${EXE}