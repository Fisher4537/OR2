CPLEX_HOME = /home/jarvis/ibm/ILOG/CPLEX_Studio1210/cplex
LIBS = -L${CPLEX_HOME}/lib/x86-64_linux/static_pic -L. -lcplex -lm -lpthread -ldl
INC = -I. -I${CPLEX_HOME}/include/ilcplex
ARGS = -input data/att12.tsp -model_type 0 -max_nodes 10000
EXE = mainTSP
OBJ = tsp.o chrono.o mainTSP.o
MODELS = 2
TEST_FILES = `ls data/*.tsp`
TRAINSET_FILES = data/att12.tsp data/att48.tsp data/pr76.tsp data/rat99.tsp data/kroB100.tsp data/a280.tsp
SEED = 0 123456 666 777 1995

test: $(EXE)
	for file in $(TEST_FILES); do \
		for model in $(MODELS); do \
			for seed in $(SEED); do \
				./$(EXE) -input $$file -model_type $$model -randomseed $$seed -v 1; \
			done \
		done \
	done

trainset: $(EXE)
	for file in $(TRAINSET_FILES); do \
		for model in $(MODELS); do \
			for seed in $(SEED); do \
				./$(EXE) -input $$file -model_type $$model -randomseed $$seed -v 10; \
			done \
		done \
	done

run: $(EXE)
	./mainTSP $(ARGS)  >> log/`date +%y%m%d-%H%M%S`.log

debug:
	@echo "compiling with debugging options..."
	gcc $(INC) -g mainTSP.c -o mainTSP.x $(LIBS)

$(EXE): $(OBJ)
	gcc -o $(EXE) $(OBJ) $(LIBS)

$(OBJ): mainTSP.c tsp.c chrono.c
	@echo "compiling..."
	gcc $(INC) -c chrono.c
	gcc $(INC) -c tsp.c
	gcc $(INC) -c mainTSP.c

clean:
	rm *.o *.x mainTSP

cleanlog:
	rm log/*

cleanmodel:
	rm model/*.lp

cleanplot:
	rm plot/*.png

cleanall: clean cleanlog cleanmodel cleanplot
