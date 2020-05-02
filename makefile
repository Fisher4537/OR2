CPLEX_HOME = /home/jarvis/ibm/ILOG/CPLEX_Studio1210/cplex
LIBS = -L${CPLEX_HOME}/lib/x86-64_linux/static_pic -L. -lcplex -lm -lpthread -ldl
INC = -I. -I${CPLEX_HOME}/include/ilcplex
ARGS = -input data_heavy/ali535.tsp -model_type 8 -v 1000 -nthread 4 -time_limit 10.0
EXE = mainTSP
OBJ = tsp.o chrono.o mainTSP.o
MODELS = 0
TEST_LIGHT = `ls data_light/*.tsp`
TEST_AVERAGE = `ls data_average/*.tsp`
TEST_SET = `ls Test_Set/*.tsp`
TRAINSET_FILES = `ls data_trainset/*.tsp`
SEED = 0 123456 666 777 1995

testlight: $(EXE)
	for file in $(TEST_LIGHT); do \
		for model in $(MODELS); do \
			for seed in $(SEED); do \
				./$(EXE) -input $$file -model_type $$model -randomseed $$seed -v 1; \
			done \
		done \
	done

testaverage: $(EXE)
	for file in $(TEST_AVERAGE); do \
		for model in $(MODELS); do \
			for seed in $(SEED); do \
				./$(EXE) -input $$file -model_type $$model -randomseed $$seed -v 1; \
			done \
		done \
	done

testall: $(EXE)
	for file in $(TEST_SET); do \
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
	gcc $(INC) -g mainTSP.c -o mainTSP.x $(LIBS)

$(EXE): $(OBJ)
	gcc -o $(EXE) $(OBJ) $(LIBS)

$(OBJ): mainTSP.c tsp.c chrono.c
	gcc $(INC) -Wall -g -c chrono.c
	gcc $(INC) -Wall -g -c tsp.c
	gcc $(INC) -Wall -g -c mainTSP.c

clean:
	rm *.o *.x mainTSP

cleanlog:
	rm log/*

cleanmodel:
	rm model/*.lp

cleanplot:
	rm plot/*.png

cleanall: clean cleanlog cleanmodel cleanplot
