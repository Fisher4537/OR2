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

#ifdef _WIN32
	// TEST code for Windows:

	int folders = 4;

	/* folder =
		0: data_light
		1: data_average
		2: data_heavy
		3: data_light + data_average
		4: data_light + data_average + data_heavy
	*/

	int seed[5] = { 0, 123456, 666, 777, 1995 };

	folders = folders == 0 ? 29 : folders == 1 ? 16 : folders == 2 ? 43 : folders == 3 ? 45 : 88;

	char** files = (char**)calloc(folders, sizeof(char*));
	for (int i = 0; i < folders; i++)
		files[i] = (char*)calloc(25, sizeof(char));

	FILE* fp;
	fp = fopen("./data/filename.txt", "r");

	int idx = 0;
	while (folders > 0 && fscanf(fp, "%s", files[idx]) == 1) {
		//printf("%s\n", files[idx]);
		idx++;
		folders--;
	}

	// END code for testing

	for (int i = 0; i < idx; i++) {
		for (int j = 0; j < 5; j++) {

			tspinstance inst;

			// parse input args
			parse_command_line(argc, argv, &inst);

			(&inst)->randomseed = seed[j];

			printf("\n____ FILE: %s ____ SEED: %d ____\n", files[i], seed[j]);
			strcpy((&inst)->input_file, files[i]);


			// read input files
			read_input(&inst);

			// TSP optimization
			if (TSPopt(&inst)) print_error(" error within TSPopt()");
			if ((&inst)->verbose >= 1) printf(" ... Executed in %.2lf s, best_lb: %.1lf\n", (&inst)->opt_time, (&inst)->best_lb);

			if ((&inst)->verbose >= 1) plot_instance(&inst);
			if ((&inst)->verbose < 10) save_results(&inst, "res.csv");
			else if ((&inst)->verbose < 100) save_results(&inst, "res_trainset.csv");

			free_instance(&inst);

			// getchar(); // pause execution
		}
	}

	#else
		tspinstance inst;

		// parse input args
		parse_command_line(argc, argv, &inst);

		// read input files
		read_input(&inst);

		// TSP optimization
		if (TSPopt(&inst)) print_error(" error within TSPopt()");
		if ((&inst)->verbose >= 1) printf(" ... Executed in %.2lf s, best_lb: %.1lf\n", (&inst)->opt_time, (&inst)->best_lb);

		if ((&inst)->verbose >= 1) plot_instance(&inst);
		if ((&inst)->verbose < 10) save_results(&inst, "res.csv");
		else if ((&inst)->verbose < 100) save_results(&inst, "res_trainset.csv");

		free_instance(&inst);

		// getchar(); // pause execution
	#endif

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
