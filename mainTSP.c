#include "tsp.c"
#include <stdio.h>
#include <time.h>


int main(int argc, char **argv)
{
	if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }
	if ( VERBOSE >= 2 ) {
     for (int a = 0; a < argc; a++) printf("%s ", argv[a]);
		 printf("\n");
  }

	time_t t1 = time(NULL);
	tspinstance inst;

	// parse input args
	parse_command_line(argc, argv, &inst);

	// read input files
	read_input(&inst);

	// TSP optimization
	if ( TSPopt(&inst) ) print_error(" error within VRPopt()");
	if ( VERBOSE >= 10 ) plot_problem_input(&inst);

	time_t t2 = time(NULL);
	if ( VERBOSE >= 1 ) printf("... TSP solved in %lf s.\n", difftime(t2,t1));
	save_results(&inst, "res.csv");
	free_instance(&inst);

	printf("END OF mainTSP, press Return to exit...");
  // getchar(); // pause execution
	return 0;
}
