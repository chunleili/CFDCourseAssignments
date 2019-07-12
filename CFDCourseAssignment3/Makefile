OBJS = $(wildcard *.o)
LIBS = -lm
EXE  = a.out
CFLG = -g -Wall -std=c++11	
start   : $(OBJS) main.H
	g++ -o ${EXE} ${CFLG} ${OBJS} ${LIBS}

main.o    : $(SRCS) main.H 
	g++  ${CFLG} -c main.C

solver1_JST.o: solver1_JST.C main.H
	g++  ${CFLG} -c solver1_JST.C

solver2_MacCormack.o: solver2_MacCormack.C main.H
	g++  ${CFLG} -c solver2_MacCormack.C

solver3_Roe.o: solver3_Roe.C main.H
	g++  ${CFLG} -c solver3_Roe.C

.PHONY:clean  clean_data clean_log clean_all
clean      :
	rm -rf *.o 
clean_data :
	rm -rf data/*
clean_log  :
	rm -rf log/*
clean_all  :
	rm -rf data/* log/* *.o ${EXE}