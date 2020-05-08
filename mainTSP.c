#include "tsp.h"
#include <stdio.h>
#include <time.h>

#define  VERBOSE 2

int load_point(char* pathFileTSP);
int order_by_dis(int firstPoint);

int main(int argc, char **argv)
{
	
	
	if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }
	if ( VERBOSE >= 2 ) {
     for (int a = 0; a < argc; a++) printf("%s ", argv[a]);
		 //printf("\n");
		 fflush(stdout);
	}

	tspinstance inst;

	// parse input args
	parse_command_line(argc, argv, &inst);

	// read input files
	read_input(&inst);

	char* pathFileTSP = (&inst)->input_file;
	load_point(pathFileTSP);
	order_by_dis(0);
	return;

	// TSP optimization
	if ( TSPopt(&inst) ) print_error(" error within TSPopt()");
	if ( (&inst)->verbose >= 1) printf(" ... Executed in %.2lf s, best_lb: %.1lf\n", (&inst)->opt_time, (&inst)->best_lb); 

	if ( (&inst)->verbose >= 1 ) plot_instance(&inst);
	if ( (&inst)->verbose < 10) save_results(&inst, "res.csv");
	else if ( (&inst)->verbose < 100 ) save_results(&inst, "res_trainset.csv");

	free_instance(&inst);

	// printf("END OF mainTSP, press Return to exit...");
	// getchar(); // pause execution
	return 0;
}
