#include "tsp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int TSPopt(tspinstance *inst);
int xpos(int i, int j, tspinstance *inst);
void build_model(tspinstance *inst, CPXENVptr env, CPXLPptr lp);
void parse_command_line(int argc, char** argv, tspinstance *inst);
void free_instance(tspinstance *inst);
void read_input(tspinstance *inst);
void print_error(const char *err);
double dist(int i, int j, tspinstance *inst);


int TSPopt(tspinstance *inst)
{

	// open cplex model
	int status;
	CPXENVptr env = CPXopenCPLEX(&status);
	CPXLPptr lp = CPXcreateprob(env, &status, "TSP");

	// set input data in CPX structure
	build_model(inst, env, lp);

	// setup struct to save solution
	inst->nedges = CPXgetnumcols(env, lp);
	inst->best_sol = (double *) calloc(inst->nedges, sizeof(double)); 	// all entries to zero
	inst->zbest = CPX_INFBOUND;

	// Cplex's parameter setting

	// compute cplex
	CPXmipopt(env,lp);

	// get best solution
	// CPXgetbestobjval(env, lp, &inst->best_lb);
	CPXsolution(env, lp, &status, &inst->best_lb, inst->best_sol, NULL, NULL, NULL);

	// use the optimal solution found by CPLEX

	// free and close cplex model
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return !status; // status 0 is ok
}

int xpos(int i, int j, tspinstance *inst)
{
	if ( i == j ) print_error(" i == j in xpos" );
	if ( i > j ) return xpos(j,i,inst);								// simplify returned formula
	return i*inst->nnodes + j - ((i + 1)*(i + 2))/2; 	// default case
}

void build_model(tspinstance *inst, CPXENVptr env, CPXLPptr lp)
{

	// double zero = 0.0;
	char binary = 'B';

	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));

	// add binary var.s x(i,j) for i < j
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);
			double obj = dist(i,j,inst); // cost == distance
			double lb = 0.0;
			double ub = 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) )
				print_error(" wrong CPXnewcols on x var.s");
  		if ( CPXgetnumcols(env,lp)-1 != xpos(i,j, inst) )
				print_error(" wrong position for x var.s");
		}
	}

	// add the degree constraints
	for ( int h = 0; h < inst->nnodes; h++ )  // degree constraints
	{
		int lastrow = CPXgetnumrows(env,lp);
		double rhs = 2.0;
		char sense = 'E';											// 'E' for equality constraint
		sprintf(cname[0], "degree(%d)", h+1);
		if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) )
			print_error(" wrong CPXnewrows [degree]");
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			if ( i == h ) continue;
			if ( CPXchgcoef(env, lp, lastrow, xpos(i,h, inst), 1.0) )
				print_error(" wrong CPXchgcoef [degree]");
		}
	}

	if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model/model.lp", NULL);

	free(cname[0]);
	free(cname);
}

double dist(int i, int j, tspinstance *inst)
{
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j];
	if ( !inst->integer_costs ) return sqrt(dx*dx+dy*dy);
	int dis = sqrt(dx*dx+dy*dy) + 0.499999999; // nearest integer
	return dis+0.0;
}

void print_error(const char *err)
{ printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1); }

void free_instance(tspinstance *inst)
{

	free(inst->xcoord);
	free(inst->ycoord);
	// free(inst->load_min);
	// free(inst->load_max);
}

void read_input(tspinstance *inst) // simplified CVRP parser, not all SECTIONs detected
{

	FILE *fin = fopen(inst->input_file, "r");
	if ( fin == NULL ) print_error(" input file not found!");

	inst->nnodes = -1;
	inst->nedges = -1;

	char line[180];
	char *par_name;
	char *token1;
	char *token2;

	int active_section = 0; // =1 NODE_COORD_SECTION, =2 DEMAND_SECTION, =3 DEPOT_SECTION

	while ( fgets(line, sizeof(line), fin) != NULL )
	{
		if ( VERBOSE >= 2000 ) { printf("[2000] line: %s",line); fflush(NULL); }
		if ( strlen(line) <= 1 ) continue; // skip empty lines

	  par_name = strtok(line, " :");
		if ( VERBOSE >= 3000 ) { printf("[3000] parameter: \"%s\" ",par_name); fflush(NULL); }

		if ( strncmp(par_name, "NAME", 4) == 0 )
		{
			active_section = 0;
			continue;
		}

		if ( strncmp(par_name, "COMMENT", 7) == 0 )
		{
			active_section = 0;
			token1 = strtok(NULL, "");
			if ( VERBOSE >= 10 ) printf(" ... solving instance %s with model %d\n\n", token1, inst->model_type);
			continue;
		}

		if ( strncmp(par_name, "TYPE", 4) == 0 )
		{
			token1 = strtok(NULL, " :");
			if ( strncmp(token1, "TSP",3) != 0 ) print_error(" format error:  only TYPE == CVRP implemented so far!!!!!!");
			active_section = 0;
			continue;
		}

		if ( strncmp(par_name, "DIMENSION", 9) == 0 )
		{
			if ( inst->nnodes >= 0 ) print_error(" repeated DIMENSION section in input file");
			token1 = strtok(NULL, " :");
			inst->nnodes = atoi(token1);
			if ( VERBOSE >= 1000 ) printf("\tnnodes: %d\n", inst->nnodes);
			inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double)*inst->nnodes);
			inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double)*inst->nnodes);
			active_section = 0;
			continue;
		}

		if ( strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0 )
		{
			token1 = strtok(NULL, " :");
			if ( strncmp(token1, "ATT", 3) != 0 ) print_error(" format error:  only EDGE_WEIGHT_TYPE == EUC_2D implemented so far!!!!!!");
			active_section = 0;
			continue;
		}

		if ( strncmp(par_name, "NODE_COORD_SECTION", 18) == 0 )
		{
			if ( inst->nnodes <= 0 ) print_error(" ... DIMENSION section should appear before NODE_COORD_SECTION section");
			active_section = 1;
			continue;
		}

		if ( strncmp(par_name, "EOF", 3) == 0 )
		{
			active_section = 0;
			break;
		}


		if ( active_section == 1 ) // within NODE_COORD_SECTION
		{
			int i = atoi(par_name) - 1;
			if ( i < 0 || i >= inst->nnodes ) print_error(" ... unknown node in NODE_COORD_SECTION section");
			token1 = strtok(NULL, " :,");
			token2 = strtok(NULL, " :,");
			inst->xcoord[i] = atof(token1);
			inst->ycoord[i] = atof(token2);
			if ( VERBOSE >= 1000 )
				printf("\tnode %4d = (%.2f, %.2f)\n", i+1, inst->xcoord[i], inst->ycoord[i]);
			continue;
		}

		printf(" final active section %d\n", active_section);
		print_error(" ... wrong format for the current simplified parser!!!!");
	}

	fclose(fin);
}

void parse_command_line(int argc, char** argv, tspinstance *inst)
{

	if ( VERBOSE >= 100 )
    printf(" running %s with %d parameters \n", argv[0], argc-1);

	// default
	inst->model_type = 0;
	inst->randomseed = 0;
	inst->timelimit = 20.; // CPX_INFBOUND;
	strcpy(inst->input_file, "NULL");

	inst->available_memory = 12000;   // available memory, in MB, for Cplex execution (e.g., 12000)
	inst->max_nodes = -1; 						// max n. of branching nodes in the final run (-1 unlimited)
  inst->integer_costs = 0;

  int help = 0; if ( argc < 1 ) help = 1;
	for ( int i = 1; i < argc; i++ )
  {
		if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-input") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-f") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 				// input file
		if ( strcmp(argv[i],"-time_limit") == 0 ) { inst->timelimit = atof(argv[++i]); continue; }		// total time limit
		if ( strcmp(argv[i],"-model_type") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 	// model type
		if ( strcmp(argv[i],"-model") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 			// model type
		if ( strcmp(argv[i],"-seed") == 0 ) { inst->randomseed = abs(atoi(argv[++i])); continue; } 		// random seed
		if ( strcmp(argv[i],"-memory") == 0 ) { inst->available_memory = atoi(argv[++i]); continue; }	// available memory (in MB)
		if ( strcmp(argv[i],"-node_file") == 0 ) { strcpy(inst->node_file,argv[++i]); continue; }		// cplex's node file
		if ( strcmp(argv[i],"-max_nodes") == 0 ) { inst->max_nodes = atoi(argv[++i]); continue; } 		// max n. of nodes
		if ( strcmp(argv[i],"-int") == 0 ) { inst->integer_costs = 1; continue; } 						// inteher costs
		if ( strcmp(argv[i],"-help") == 0 ) { help = 1; continue; } 									// help
		if ( strcmp(argv[i],"--help") == 0 ) { help = 1; continue; } 									// help
		help = 1;
  }

	if ( help || (VERBOSE >= 10) )		// print current parameters
	{
		printf("[VERBOSE 10] available parameters:\n");
		printf("-file %s\n", inst->input_file);
		printf("-time_limit %f\n", inst->timelimit);
		printf("-model_type %d\n", inst->model_type);
		printf("-seed %d\n", inst->randomseed);
		printf("-max_nodes %d\n", inst->max_nodes);
		printf("-memory %d\n", inst->available_memory);
		printf("-int %d\n", inst->integer_costs);
		printf("-node_file %s\n", inst->node_file);
		printf("---------------------------------\n\n");
	}

	if ( help ) exit(1);

}

void plot_problem_input(tspinstance *inst)
{
	// TODO: consider to use a gnuplot interface if things become complicated
	// open gnuplot process
  FILE *gnuplot = popen("gnuplot", "w");
	char *pngname = "plot/fig.png";

	// set up line and point style
	fprintf(gnuplot,"set style line 1 \
							    lc rgb '#0060ad' \
									pointtype 7 \
									pointsize 1.0\n");

	fprintf(gnuplot,"set style line 2 \
							    lc rgb '#0060ad' \
									linetype 1 linewidth 2 \
									pointtype 7 pointsize 1.0\n");

	// set title
	fprintf(gnuplot, "set title \"Input data\"\n");

	// start plotting points
  fprintf(gnuplot, "plot '-' w p ls 1\n");
  for (size_t i = 0; i < inst->nnodes; i++)
    fprintf(gnuplot, "%f %f\n", inst->xcoord[i], inst->ycoord[i]);
	fprintf(gnuplot, "e\n");

	// start plotting edges
	if (inst->nedges > 0) // check if solution is available
	{

		// save png into FILE
		if (VERBOSE < 1000) // save plot in file
			fprintf(gnuplot, "set terminal png\n\
												set output '%s'\n", pngname);

		fprintf(gnuplot, "plot '-' w linespoints linestyle 2\n");
		for (int i = 0; i < inst->nnodes; i++)
			for (int j = i+1; j < inst->nnodes; j++)
				if (0.5 < inst->best_sol[xpos(i,j,inst)] ) // && inst->best_sol[i] < 1.00001)
					fprintf(gnuplot, "%f %f\n%f %f\n\n",
								 		inst->xcoord[i], inst->ycoord[i],
										inst->xcoord[j], inst->ycoord[j]);
		fprintf(gnuplot, "e\n");
	}

	// show plot or save in file and close
  fflush(gnuplot);
	if (VERBOSE >= 1000) getchar(); // pause execution to see the plot
	fclose(gnuplot);
}
