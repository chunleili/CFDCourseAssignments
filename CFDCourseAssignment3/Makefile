OBJS = main.o solver3_Roe.o #solver1_JST.o solver2_MacCormack.o 
#$(wildcard *.o)
LIBS = -lm
EXE  = a.out
CFLG = -g -Wall -std=c++11 
main   : $(OBJS) main.H
	g++ -o ${EXE} ${CFLG} ${OBJS} ${LIBS}

main.o    : $(SRCS) main.H 
	g++  ${CFLG} -c main.C

#solver1_JST.o: solver1_JST.C main.H
#	g++  ${CFLG} -c solver1_JST.C
#
#solver2_MacCormack.o: solver2_MacCormack.C main.H
#	g++  ${CFLG} -c solver2_MacCormack.C
#
solver3_Roe.o: solver3_Roe.C main.H
	g++  ${CFLG} -c solver3_Roe.C

.PHONY:clean
clean      :
	rm -rf *.dat log/* *.o ${EXE}