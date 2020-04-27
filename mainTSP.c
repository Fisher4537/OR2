#include "pch.h"
#include "tsp.c"
#include <stdio.h>
#include <time.h>

#define  VERBOSE 2

int main(int argc, char **argv)
{
	if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }
	if ( VERBOSE >= 2 ) {
     for (int a = 0; a < argc; a++) printf("%s ", argv[a]);
		 //printf("\n");
  }

	tspinstance inst;

	// parse input args
	parse_command_line(argc, argv, &inst);

	// read input files
	read_input(&inst);

	// TSP optimization
	if ( TSPopt(&inst) ) print_error(" error within VRPopt()");
	if ( (&inst)->verbose >= 1) printf(" ... Executed in %6.3lf s\n", (&inst)->opt_time);

	if ( (&inst)->verbose >= 1 ) plot_instance(&inst);
	if ( (&inst)->verbose <= 1) save_results(&inst, "res.csv");
	else if ( (&inst)->verbose <= 10 ) save_results(&inst, "res_trainset.csv");

	free_instance(&inst);

	// printf("END OF mainTSP, press Return to exit...");
	// getchar(); // pause execution
	return 0;
}
