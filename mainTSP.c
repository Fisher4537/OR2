#include "tsp.h"
#include <stdio.h>
#include <time.h>


#define  VERBOSE 2


int main(int argc, char **argv)
{


	
	if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }
	if ( VERBOSE >= 2 ) {
     for (int a = 0; a < argc; a++) printf("%s ", argv[a]);
		 //printf("\n");
		 fflush(stdout);
	}

	int seed[5] = { 0, 123456, 666, 777, 1995 };
	/*
	char files[45][15] = {  "att48.tsp",
							"berlin52.tsp",
							"bier127.tsp",
							"burma14.tsp",
							"ch130.tsp",
							"ch150.tsp",
							"eil101.tsp",
							"eil51.tsp",
							"eil76.tsp",
							"gr137.tsp",
							"gr96.tsp",
							"kroA100.tsp",
							"kroA150.tsp",
							"kroB100.tsp",
							"kroB150.tsp",
							"kroC100.tsp",
							"kroD100.tsp",
							"kroE100.tsp",
							"lin105.tsp",
							"pr107.tsp",
							"pr124.tsp",
							"pr136.tsp",
							"pr144.tsp",
							"pr76.tsp",
							"rat99.tsp",
							"rd100.tsp",
							"st70.tsp",
							"ulysses16.tsp",
							"ulysses22.tsp"
	};
	*/
	char files[21][15] = {  "att48.tsp",
							"berlin52.tsp",
							"bier127.tsp",
							"burma14.tsp",
							"ch130.tsp",
							"ch150.tsp",
							"eil101.tsp",
							"eil51.tsp",
							"eil76.tsp",
							"gr137.tsp",
							"gr96.tsp",
							"kroA100.tsp",
							"kroB100.tsp",
							"kroC100.tsp",
							"lin105.tsp",
							"pr107.tsp",
							"pr124.tsp",
							"pr76.tsp",
							"st70.tsp",
							"ulysses16.tsp",
							"ulysses22.tsp"
	};
	for (int i = 12; i < 21; i++) {
		for (int j = 0; j < 2; j++) {
			tspinstance inst;

			// parse input args
			parse_command_line(argc, argv, &inst);

			(&inst)->randomseed = seed[j];
			char name[100] = "Test_Set\\";
			concatenate_string(name, files[i]);
			strcpy((&inst)->input_file, name);

			// read input files
			read_input(&inst);

			// TSP optimization
			if (TSPopt(&inst)) print_error(" error within TSPopt()");
			if ((&inst)->verbose >= 1) printf(" ... Executed in %.2lf s, best_lb: %.1lf\n", (&inst)->opt_time, (&inst)->best_lb);

			if ((&inst)->verbose >= 1) plot_instance(&inst);
			if ((&inst)->verbose < 10) save_results(&inst, "res.csv");
			else if ((&inst)->verbose < 100) save_results(&inst, "res_trainset.csv");

			free_instance(&inst);
		}
	}
	// printf("END OF mainTSP, press Return to exit...");
	// getchar(); // pause execution
	return 0;
}

void concatenate_string(char* original, char* add)
{
	while (*original)
		original++;

	while (*add)
	{
		*original = *add;
		add++;
		original++;
	}
	*original = '\0';
}
