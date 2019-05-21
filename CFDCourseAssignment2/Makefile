OBJS = setField.o print.o iterate.o main.o 
SRCS = setField.C print.C iterate.C main.C 
LIBS = -lm
EXE  = out	
main       : $(OBJS) main.H
	g++ -o ${EXE} $(OBJS) $(LIBS)

print.o   : print.C main.H
	g++ -g -c print.C

setField.o : setField.C main.H
	g++ -g -c setField.C
	
iterate.o  : iterate.C setField.C main.H
	g++ -g -c iterate.C
	
main.o     : $(SRCS) main.H
	g++ -g -c main.C

.PHONY:clean  clean_data clean_log clean_all
clean      :
	rm -rf *.o 
clean_data :
	rm -rf data/*
clean_log  :
	rm -rf log/*
clean_all  :
	rm -rf data/* log/* *.o ${EXE}
