CPLEX_HOME = /home/jarvis/ibm/ILOG/CPLEX_Studio1210/cplex
LIBS = -L${CPLEX_HOME}/lib/x86-64_linux/static_pic -L. -lcplex -lm -lpthread -ldl
INC = -I. -I${CPLEX_HOME}/include/ilcplex
ARGS = -input data/att12.tsp -model_type 0 -max_nodes 10000
EXE = mainTSP
OBJ = mainTSP.o
MODELS = 0
TEST_FILES = ls data/*.tsp

test: $(EXE)
	for file in $(TEST_FILES); do \
		for model in $(MODELS); do \
			./$(EXE) -input $$file -model_type $$model -v 100; \
		done \
	done

run: $(EXE)
	./mainTSP $(ARGS)  >> log/`date +%y%m%d-%H%M%S`.log

debug:
	@echo "compiling with debugging options..."
	gcc $(INC) -g mainTSP.c -o mainTSP.x $(LIBS)

$(EXE): $(OBJ)
	gcc -o $(EXE) $(OBJ) $(LIBS)

$(OBJ): mainTSP.c tsp.c tsp.h
	@echo "compiling..."
	gcc $(INC) -c mainTSP.c

clean:
	rm *.o *.x mainTSP

cleanlog:
	rm log/*

cleanmodel:
	rm *.lp

cleanplot:
	rm *.png

cleanall: clean cleanlog cleanmodel cleanplot
