CONV  = scalarJSTConv
#CONV = MacCormackConv
#CONV = RoeConv
OBJS = init.o print.o ${CONV}.o main.o WToF.o localTime.o BC.o EulerFTime.o
SRCS = init.C print.C ${CONV}.C main.C WToF.C localTime.C BC.C EulerFTime.C
LIBS = -lm
EXE  = out
CFLG = -g -Wall -std=c++11	
start   : $(OBJS) main.H
	g++ -o ${EXE} ${CFLG} ${OBJS} ${LIBS}

localTime.o :	localTime.C main.H
	g++ ${CFLG} -c localTime.C

print.o : print.C main.H
	g++  ${CFLG} -c print.C

init.o : init.C main.H
	g++  ${CFLG} -c init.C

WToF.o : WToF.C main.H
	g++  ${CFLG} -c WToF.C

${CONV}.o  : ${CONV}.C WToF.C main.H 
	g++  ${CFLG} -c ${CONV}.C

BC.o : BC.C main.H
	g++ ${CFLG} -c BC.C
EulerFTime.o : EulerFTime.C main.H
	g++ ${CFLG}	-c EulerFTime.C

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