CPLEX_HOME = /home/jarvis/ibm/ILOG/CPLEX_Studio1210/cplex
LIBS = -L${CPLEX_HOME}/lib/x86-64_linux/static_pic -L. -lcplex -lm -lpthread -ldl
INC = -I. -I${CPLEX_HOME}/include/ilcplex
ARGS = -input data/att48.tsp -time_limit 10. -model_type 10 -max_nodes 1000
EXE = mainTSP
OBJ = mainTSP.o


run: debug $(EXE)
	./mainTSP $(ARGS)  >> log/`date +%y%m%d-%H%M%S`.log

debug:
	@echo "compiling with debugging options..."
	gcc $(INC) -g mainTSP.c -o mainTSP.x

$(EXE): $(OBJ)
	gcc -o $(EXE) $(OBJ) $(LIBS)

$(OBJ): mainTSP.c
	@echo "compiling..."
	gcc $(INC) -c mainTSP.c

clean:
	rm *.o *.x mainTSP

cleanlog:
	rm log/*
