#include "pch.h"         
#include <stdio.h>
#include <time.h>

#include "tsp.h"  

int main(int argc, char **argv){

	printf("Parameter Used : \n");
	for (int i = 0; i < argc; i++)
		printf("%s\n", argv[i]);
	printf("\n\n");

	if (argc < 2) {		
		printf("Usage: -help for help\n"); 
		exit(1);
	}
	
	tspInstance inst;
	result res;

	// parse input args
	parse_command_line(argc, argv, &inst);

	// read input files
	read_input(&inst, &res);

	// TSP optimization
	if (TSPopt(&inst, &res)) 
		print_error(" error within TSPopt()");
	
	get_pipe(&inst, &res);

	printf("TSP solved in %lf s. (CPS=%lf)\n\n", (&res)->opt_time, (double)CLOCKS_PER_SEC);

	save_results(&inst, &res, "Sol/res.csv");
	free_instance(&inst);

	printf("Finish, press Return to exit...");

	// getchar();
	
	return 0;
}

