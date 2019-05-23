OBJS = init.o print.o MacCormackSolver.o main.o UToF.o
SRCS = init.C print.C MacCormackSolver.C main.C UToF.C
LIBS = -lm
EXE  = out
CFLG = -g -Wall -std=c++11	
main       : $(OBJS) main.H
	g++ -o ${EXE} ${CFLG} ${OBJS} ${LIBS}

print.o   : print.C main.H
	g++ -g ${CFLG} -c print.C

init.o : init.C main.H
	g++ -g ${CFLG} -c init.C

UToF.o : UToF.C main.H
	g++ -g ${CFLG} -c UToF.C
	
MacCormackSolver.o  : MacCormackSolver.C UToF.C main.H
	g++ -g ${CFLG} -c MacCormackSolver.C
	
main.o     : $(SRCS) main.H
	g++ -g ${CFLG} -c main.C

.PHONY:clean  clean_data clean_log clean_all
clean      :
	rm -rf *.o 
clean_data :
	rm -rf data/*
clean_log  :
	rm -rf log/*
clean_all  :
	rm -rf data/* log/* *.o ${EXE}