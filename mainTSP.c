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

	parse_command_line(argc, argv, &inst);

	//printf(" file %s has %d non-empty lines\n", inst.input_file, number_of_nonempty_lines(inst.input_file)); exit(1);

	read_input(&inst);
	if ( TSPopt(&inst) ) print_error(" error within VRPopt()");
	time_t t2 = time(NULL);

	if ( VERBOSE >= 1 )
	{
		printf("... VRP problem solved in %lf sec.s\n", difftime(t2,t1));
	}

	free_instance(&inst);

	printf("END OF mainTSP, press Return to exit...");
  getchar();
	return 0;
}
