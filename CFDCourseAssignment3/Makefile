CONV  = scalarJSTConv
#CONV = MacCormackConv
#CONV = RoeConv
TIME_DIS1 = RungeKuttaTime
TIME_DIS2 = EulerFTime
OBJS = init.o print.o ${CONV}.o main.o WToF.o localTime.o BC.o ${TIME_DIS1}.o ${TIME_DIS2}.o calDtGlobal.o
SRCS = init.C print.C ${CONV}.C main.C WToF.C localTime.C BC.C ${TIME_DIS2}.C ${TIME_DIS2}.C calDtGlobal.C
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

calDtGlobal.o : calDtGlobal.C main.H
	g++  ${CFLG} -c calDtGlobal.C

${CONV}.o  : ${CONV}.C WToF.C main.H 
	g++  ${CFLG} -c ${CONV}.C

BC.o : BC.C main.H
	g++ ${CFLG} -c BC.C

minOverScalarField.o : minOverScalarField.C main.H
	g++ ${CFLG} -c minOverScalarField.C

${TIME_DIS1}.o : ${TIME_DIS1}.C main.H
	g++ ${CFLG}	-c ${TIME_DIS1}.C

${TIME_DIS2}.o : ${TIME_DIS2}.C main.H
	g++ ${CFLG}	-c ${TIME_DIS2}.C

main.o    : $(SRCS) main.H 
	g++  ${CFLG} -c main.C

.PHONY:clean_o  clean_data clean_log clean_all
clean_o      :
	rm -rf *.o 
clean_data :
	rm -rf data/*
clean_log  :
	rm -rf log/*
clean_all  :
	rm -rf data/* log/* *.o ${EXE}