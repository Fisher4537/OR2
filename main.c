#include "pch.h"
#include "tsp.h"           

void debug(const char *err);
void print_error(const char *err);
void read_input(instance *inst);
void parse_command_line(int argc, char** argv, instance *inst);
void free_instance(instance *inst);

void debug(const char *err) { 
	printf("\nDEBUG: %s \n", err); 
	fflush(NULL);
}
void print_error(const char *err) { 
	printf("\n\n ERROR: %s \n\n", err); 
	fflush(NULL); 
	exit(1); 
}


int VRPopt(instance *inst);

int main(int argc, char **argv){

	printf("Parameter Used : \n");
	for (int i = 0; i < argc; i++)
		printf("%s\n", argv[i]);
	printf("\n\n");

	if (argc < 2) {
		
		printf("Usage: -help for help\n"); 
		exit(1);
	}

	if (VERBOSE >= 2) { 
		for (int a = 0; a < argc; a++) 
			printf("%s ", argv[a]); printf("\n"); 
	}

	

	instance inst;

	parse_command_line(argc, argv, &inst);

	read_input(&inst);

	/*
	if (VRPopt(&inst)) 
		print_error(" error within VRPopt()");
	*/

	free_instance(&inst);
	return 0;
}


void parse_command_line(int argc, char** argv, instance *inst){

	if (VERBOSE >= 100) 
		printf(" running %s with %d parameters \n", argv[0], argc - 1);

	// default   
	inst->model_type = 0;
	strcpy(inst->input_file, "NULL");
	strcpy(inst->node_file, "NULL");
	inst->randomseed = 0;
	inst->num_threads = 0;
	inst->timelimit = CPX_INFBOUND;			// CPX_INFBOUND = "infinito" di CPLEX


	inst->available_memory = 12000;   			// available memory, in MB, for Cplex execution (e.g., 12000)
	inst->max_nodes = -1; 						// max n. of branching nodes in the final run (-1 unlimited)        

	int help = 0; 
	if (argc < 1) 
		help = 1;
	
	for (int i = 1; i < argc; i++){

		if (strcmp(argv[i], "-file") == 0) { strcpy(inst->input_file, argv[++i]); continue; } 			// input file
		if (strcmp(argv[i], "-input") == 0) { strcpy(inst->input_file, argv[++i]); continue; } 			// input file
		if (strcmp(argv[i], "-f") == 0) { strcpy(inst->input_file, argv[++i]); continue; } 				// input file
		if (strcmp(argv[i], "-time_limit") == 0) { inst->timelimit = atof(argv[++i]); continue; }		// total time limit
		if (strcmp(argv[i], "-model_type") == 0) { inst->model_type = atoi(argv[++i]); continue; } 		// model type
		if (strcmp(argv[i], "-model") == 0) { inst->model_type = atoi(argv[++i]); continue; } 			// model type
		if (strcmp(argv[i], "-seed") == 0) { inst->randomseed = abs(atoi(argv[++i])); continue; } 		// random seed
		if (strcmp(argv[i], "-threads") == 0) { inst->num_threads = atoi(argv[++i]); continue; } 		// n. threads
		if (strcmp(argv[i], "-memory") == 0) { inst->available_memory = atoi(argv[++i]); continue; }	// available memory (in MB)
		if (strcmp(argv[i], "-node_file") == 0) { strcpy(inst->node_file, argv[++i]); continue; }		// cplex's node file
		if (strcmp(argv[i], "-max_nodes") == 0) { inst->max_nodes = atoi(argv[++i]); continue; } 		// max n. of nodes
		if (strcmp(argv[i], "-help") == 0) { help = 1; continue; } 										// help
		if (strcmp(argv[i], "--help") == 0) { help = 1; continue; } 									// help
		help = 1;
	}

	// print current parameters
	if (help || (VERBOSE >= 10)) {		
		printf("\n\navailable parameters (vers. 15-03-2020) --------------------------------------------------\n");
		printf("-file %s\n", inst->input_file);
		printf("-time_limit %lf\n", inst->timelimit);
		printf("-model_type %d\n", inst->model_type);
		printf("-seed %d\n", inst->randomseed);
		printf("-threads %d\n", inst->num_threads);
		printf("-max_nodes %d\n", inst->max_nodes);
		printf("-memory %d\n", inst->available_memory);
		printf("-node_file %s\n", inst->node_file);
		printf("\nenter -help or --help for help\n");
		printf("----------------------------------------------------------------------------------------------\n\n");
	}

	if (help) 
		exit(1);

}

// simplified CVRP parser, not all SECTIONs detected 
void read_input(instance *inst) {

	FILE *fin = fopen(inst->input_file, "r");
	if (fin == NULL)
		print_error(" input file not found!");

	inst->nnodes = -1;

	char line[180];
	char *par_name;
	char *token1;
	char *token2;

	int active_section = 0;		// = 1 NODE_COORD_SECTION, = 2 DEMAND_SECTION, = 3 DEPOT_SECTION 

	int do_print = (VERBOSE >= 10);

	while (fgets(line, sizeof(line), fin) != NULL) {

		if (VERBOSE >= 2000) {
			printf("%s", line);
			fflush(NULL);
		}

		if (strlen(line) <= 1)
			continue; // skip empty lines

		par_name = strtok(line, " :");

		if (VERBOSE >= 3000) {
			printf("parameter \"%s\" ", par_name);
			fflush(NULL);
		}

		if (strncmp(par_name, "NAME", 4) == 0) {
			active_section = 0;
			continue;
		}

		if (strncmp(par_name, "COMMENT", 7) == 0) {
			active_section = 0;
			token1 = strtok(NULL, "");
			if (VERBOSE >= 10) 
				printf(" ... solving instance %s with model %d\n\n", token1, inst->model_type);
			continue;
		}

		if (strncmp(par_name, "TYPE", 4) == 0) {
			token1 = strtok(NULL, " :");
			if (strncmp(token1, "TSP", 3) != 0)
				print_error(" format error:  only TYPE == TSP implemented so far!");
			active_section = 0;
			continue;
		}


		if (strncmp(par_name, "DIMENSION", 9) == 0) {
			if (inst->nnodes >= 0)
				print_error(" repeated DIMENSION section in input file");

			token1 = strtok(NULL, " :");
			inst->nnodes = atoi(token1);

			if (do_print) 
				printf(" ... nnodes %d\n", inst->nnodes);
				
			inst->xcoord = (double *)calloc(inst->nnodes, sizeof(double));
			inst->ycoord = (double *)calloc(inst->nnodes, sizeof(double));

			active_section = 0;
			continue;
		}

		if (strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0) {
			token1 = strtok(NULL, " :");
			if (strncmp(token1, "ATT", 3) != 0)
				print_error(" format error:  only EDGE_WEIGHT_TYPE == ATT implemented so far!");
			active_section = 0;
			continue;
		}

		if (strncmp(par_name, "NODE_COORD_SECTION", 18) == 0) {
			if (inst->nnodes <= 0)
				print_error(" ... DIMENSION section should appear before NODE_COORD_SECTION section");
			active_section = 1;
			continue;
		}

		if (strncmp(par_name, "EOF", 3) == 0) {
			active_section = 0;
			break;
		}

		// within NODE_COORD_SECTION
		if (active_section == 1) {
			int i = atoi(par_name) - 1;

			if (i < 0 || i >= inst->nnodes)
				print_error(" ... unknown node in NODE_COORD_SECTION section");

			token1 = strtok(NULL, " :,");
			token2 = strtok(NULL, " :,");
			inst->xcoord[i] = atof(token1);
			inst->ycoord[i] = atof(token2);
			
			if (do_print) 
				printf(" ... node %4d at coordinates ( %15.7lf , %15.7lf )\n", i + 1, inst->xcoord[i], inst->ycoord[i]);
			
			continue;
		}

		printf(" final active section %d\n", active_section);
		print_error(" ... wrong format for the current simplified parser!!!!!!!!!");
	}

	fclose(fin);

}

void free_instance(instance *inst) {
	free(inst->xcoord);
	free(inst->ycoord);
}
