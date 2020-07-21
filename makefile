CPLEX_HOME = /home/jarvis/ibm/ILOG/CPLEX_Studio1210/cplex
LIBS = -L${CPLEX_HOME}/lib/x86-64_linux/static_pic -L. -lcplex -lm -lpthread -ldl
INC = -I. -I${CPLEX_HOME}/include/ilcplex
ARGS = -input data_heavy/ali535.tsp -setup_model 8 -v 1000 -nthread 4 -time_limit 600
EXE = mainTSP
OBJ = tsp.o chrono.o mainTSP.o
MODELS = 12
ALL_MODELS = 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
TEST_LIGHT = `ls data_light/*.tsp`
TEST_AVERAGE = `ls data_average/*.tsp`
TEST_SET = `ls Test_Set/*.tsp`
FAST_TEST = `ls fast_test/*.tsp`
TRAINSET_FILES = `ls data_trainset/*.tsp`
TRAINSET_HEAVY = data_heavy/att532.tsp data_heavy/fl417.tsp data_heavy/d657.tsp
SEED = 0 123456 666 777 1995

testlight: $(EXE)
	for file in $(TEST_LIGHT); do \
		for model in $(MODELS); do \
			for seed in $(SEED); do \
				./$(EXE) -input $$file -setup_model $$model -randomseed $$seed -v 1 -nthread 4 -time_limit 600; \
			done \
		done \
	done

testaverage: $(EXE)
	for file in $(TEST_AVERAGE); do \
		for model in $(MODELS); do \
			for seed in $(SEED); do \
				./$(EXE) -input $$file -setup_model $$model -randomseed $$seed -v 1 -nthread 4 -time_limit 600; \
			done \
		done \
	done

testall: $(EXE)
	for file in $(TEST_SET); do \
		for model in $(MODELS); do \
			for seed in $(SEED); do \
				./$(EXE) -input $$file -setup_model $$model -randomseed $$seed -v 1 -nthread 4 -time_limit 600; \
			done \
		done \
	done

fasttest: $(EXE)
	for file in $(FAST_TEST); do \
		for model in $(ALL_MODELS); do \
			./$(EXE) -input $$file -setup_model $$model -randomseed 0 -v 1 -nthread 4 -time_limit 4; \
		done \
	done

trainset: $(EXE)
	for file in $(TRAINSET_FILES); do \
		for model in $(MODELS); do \
			for seed in $(SEED); do \
				./$(EXE) -input $$file -setup_model $$model -randomseed $$seed -v 10 -nthread 4 -time_limit 600; \
			done \
		done \
	done

trainsetheavy: $(EXE)
	for file in $(TRAINSET_HEAVY); do \
		for model in $(MODELS); do \
			for seed in $(SEED); do \
				./$(EXE) -input $$file -setup_model $$model -randomseed $$seed -v 10 -nthread 4 -time_limit 600; \
			done \
		done \
	done

test_methods: $(EXE)
	for model in $(ALL_MODELS); do \
			./$(EXE) -input data_light/att48.tsp -setup_model $$model -v 10 -nthread 4 -time_limit 100; \
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
