#include "tsp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int TSPopt(tspinstance *inst);
int xpos(int i, int j, tspinstance *inst);
void build_model(tspinstance *inst, CPXENVptr env, CPXLPptr lp); // interface
void build_sym_std(tspinstance *inst, CPXENVptr env, CPXLPptr lp); // sym, std
void build_mtz(tspinstance *inst, CPXENVptr env, CPXLPptr lp); // asym, MTZ
void build_sol(tspinstance *inst, int *succ, int *comp, int *ncomp);
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
	CPXsetintparam(env,CPXPARAM_Threads, 3);		// allow executing N parallel threads

	// set input data in CPX structure
	build_model(inst, env, lp);

	// setup struct to save solution
	inst->nedges = CPXgetnumcols(env, lp);
	inst->best_sol = (double *) calloc(inst->nedges, sizeof(double)); 	// all entries to zero
	inst->zbest = CPX_INFBOUND;

	// Cplex's parameter setting
	printf("build model succesfully.\n");
	getchar();
	printf("optimizing model...\n");

	// compute cplex
	time_t init_time = time(NULL);
	CPXmipopt(env,lp);
	inst->opt_time = init_time - time(NULL);

	// get best solution
	// CPXgetbestobjval(env, lp, &inst->best_lb);
	CPXsolution(env, lp, &status, &inst->best_lb, inst->best_sol, NULL, NULL, NULL);

	// use the optimal solution found by CPLEX
	int *succ = (int*) calloc(inst->nnodes, sizeof(int));
	int *comp = (int*) calloc(inst->nnodes, sizeof(int));
	int *ncomp = (int*) calloc(1, sizeof(int));
	build_sol(inst, succ, comp, ncomp);

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

int asym_xpos(int i, int j, tspinstance *inst)
{
	if ( i == j ) print_error(" i == j in asym_upos" );
	return i*(inst->nnodes - 1) + ( i < j ? j-1 : j );
}

int asym_upos(int i, tspinstance *inst)
{
	if ( i < 1 ) print_error(" i < 1 in asym_upos" );
	return inst->nnodes*(inst->nnodes - 1) + i - 1;
}

void build_model(tspinstance *inst, CPXENVptr env, CPXLPptr lp)
{

	switch (inst->model_type)
	{
		case 0 :		// basic model with asymmetric x and q
			build_sym_std(inst, env,lp);
			break;

		case 1 :		// MTZ contraints
		 	build_mtz(inst, env,lp);
			break;

		default:
			print_error(" model type unknown!!");
			break;
	}
}

void build_sym_std(tspinstance *inst, CPXENVptr env, CPXLPptr lp)
{

	// double zero = 0.0;
	char xctype = CPX_BINARY;

	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));

	// add xctype var.s x(i,j) for i < j
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);
			double obj = dist(i,j,inst); // cost == distance
			double lb = 0.0;
			double ub = 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &xctype, cname) )
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

void build_mtz(tspinstance *inst, CPXENVptr env, CPXLPptr lp)
{

	char xctype = CPX_BINARY;	// type of variable
	double obj; // objective function constant
	double lb;	// lower bound
	double ub;	// upper bound

	char **cname = (char **) calloc(1, sizeof(char *));	// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));			// name of the variable

	// add binary constraints and objective const x(i,j) for i < j
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = 0; j < inst->nnodes; j++ )
		{
			if (i != j)
			{
				sprintf(cname[0], "x(%d,%d)", i+1,j+1);
				obj = dist(i,j,inst); // cost == distance
				lb = 0.0;
				ub = i == j ? 0.0 : 1.0;
				if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &xctype, cname) )
					print_error(" wrong CPXnewcols on x var.s");
				if ( CPXgetnumcols(env,lp)-1 != asym_xpos(i,j, inst) )
					print_error(" wrong position for x var.s");
			}
		}
	}

	// add nodes index in the circuits
	xctype = CPX_INTEGER;		// maybe not necessary
	obj = 0.0;
	lb = 0.0;
	ub = inst->nnodes-2;

	for ( int i = 1; i < inst->nnodes; i++)
	{
		sprintf(cname[0], "u(%d)", i+1);
		if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &xctype, cname) )
			print_error(" wrong CPXnewcols on u var.s");
		if ( CPXgetnumcols(env,lp)-1 != asym_upos(i, inst) )
			print_error(" wrong position for u var.s");
	}


	// add the degree constraints
	int lastrow;	// the number of rows
	double rhs;		// right head size
	char sense;		// 'L', 'E' or 'G'
	int big_M = inst->nnodes - 1;
	for ( int h = 0; h < inst->nnodes; h++ )  // degree constraints
	{
		lastrow = CPXgetnumrows(env,lp);
		rhs = 1.0;
		sense = 'E';	// 'E' for equality constraint
		sprintf(cname[0], "income_d(%d)", h+1);

		// create a new row
		if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) )
			print_error(" wrong CPXnewrows [degree]");

		// income vertex
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			if ( i == h ) continue;
			if ( CPXchgcoef(env, lp, lastrow, asym_xpos(i,h, inst), 1.0) ) // income vertex
				print_error(" wrong CPXchgcoef [degree]");
		}


		// outcome vertex
		lastrow = CPXgetnumrows(env,lp);
		rhs = 1.0;
		sense = 'E';	// 'E' for equality constraint
		sprintf(cname[0], "outcome_d(%d)", h+1);

		// create a new row
		if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) )
			print_error(" wrong CPXnewrows [degree]");
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			if ( i == h ) continue;
			if ( CPXchgcoef(env, lp, lastrow, asym_xpos(h,i, inst), 1.0) ) // income vertex
				print_error(" wrong CPXchgcoef [degree]");
		}


		// u constraints: nodes index
		if (h == 0) continue;	// skip when i == 0
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			if ( i != h && i != 0 )
			{
				lastrow = CPXgetnumrows(env,lp);
				rhs = big_M -2;
				sense = 'L';	// 'L' for lower of equal
				sprintf(cname[0], "mtz_i(%d, %d)", h+1, i+1);

				// create a new row
				if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) )	// new row
					print_error(" wrong CPXnewrows [degree]");

				if ( CPXchgcoef(env, lp, lastrow, asym_upos(h, inst), 1.0) ) // u constraints
					print_error(" wrong CPXchgcoef [degree]");
				if ( CPXchgcoef(env, lp, lastrow, asym_upos(i, inst), -1.0) ) // u constraints
					print_error(" wrong CPXchgcoef [degree]");
				if ( CPXchgcoef(env, lp, lastrow, asym_xpos(h,i, inst), big_M) ) // u constraints
					print_error(" wrong CPXchgcoef [degree]");
			}
		}
	}

	if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model/asym_model.lp", NULL);

	free(cname[0]);
	free(cname);
}


// #define DEBUG    // da commentare se non ci vuole il debugging
#define EPS 1e-5

void build_sol(tspinstance *inst, int *succ, int *comp, int *ncomp) // build succ() and comp() wrt xstar()...
{

	// check if nodes degree is 2 for each node
#ifdef DEBUG

	int *degree = (int *) calloc(inst->nnodes, sizeof(int));
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			int k = xpos(i,j,inst);
			if ( fabs(inst->best_sol[k]) > EPS && fabs(inst->best_sol[k]-1.0) > EPS ) print_error(" wrong inst->best_sol in build_sol()");
			if ( inst->best_sol[k] > 0.5 )
			{
				++degree[i];
				++degree[j];
			}
		}
	}
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		if ( degree[i] != 2 ) print_error("wrong degree in build_sol()");
	}
	free(degree);
#endif

	// initialization of succ, comp and ncomp
	*ncomp = 0;
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		succ[i] = -1;
		comp[i] = -1;
	}

	// tour
	for ( int start = 0; start < inst->nnodes; start++ )
	{
		if ( comp[start] >= 0 ) continue;  // node "start" has already been setted

		// a new component is found
		(*ncomp)++;
		comp[start] = *ncomp;

		int i = start;
		//int done = 0;
		while ( succ[i] == -1  )//!done )  // go and visit the current component
		{
			comp[i] = *ncomp;
			// done = 1;
			for ( int j = 0; j < inst->nnodes; j++ )
			{
				if (j == i) continue;

				if ( inst->best_sol[xpos(i,j,inst)] > 0.5) // the edge [i,j] is selected in inst->best_sol and j was not visited before
				{
					// intern edge of the cycle
					if (comp[j] == -1)
					{
						succ[i] = j;
						i = j;
						break;
					}
					// last edge of the cycle
					if (start == j)
					{
						succ[i] = j;
					}
				}
			}
		}	// while
	// go to the next component...
	}

	// print succ, comp and ncomp
#ifdef DEBUG
	printf("\ni:      ");
	for (int i = 0; i < inst->nnodes; i++) printf("%6d", i);
	printf("\nsucc:   ");
	for (int i = 0; i < inst->nnodes; i++) printf("%6d", succ[i]);
	printf("\ncomp:   ");
	for (int i = 0; i < inst->nnodes; i++) printf("%6d", comp[i]);
	printf("\n");
#endif
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

int save_results(tspinstance *inst, char *f_name)
{
	FILE *outfile;
	outfile = fopen(f_name, "a");
	char dataToAppend[200];
	snprintf(dataToAppend, sizeof dataToAppend, "%s; %d; %d; %d; %lf ;\n",
	 			inst->input_file, inst->nnodes, inst->model_type, 3, inst->opt_time ); // 3 is number of thread
	/* fopen() return NULL if unable to open file in given mode. */
	if (outfile == NULL)
	{
		/* Unable to open file hence exit */
		printf("\nUnable to open '%s' file.\n", f_name);
		printf("Please check whether file exists and you have write privilege.\n");
		exit(EXIT_FAILURE);
	}

	/* Append data to file */
	fputs(dataToAppend, outfile);

	fclose(outfile);
}
