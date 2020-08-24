#include "pch.h"
#include "tsp.h"
#include "chrono.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>

#include <omp.h>

#ifdef _WIN32
	#define DIR_DELIM '\\'
#else
	#define DIR_DELIM '/'
	#include <unistd.h>
#endif

/**
	TODO:	- a way to get CPLEX error code: status?
			- PLZ change all code with index-value, change_coef is really bad (said by Fischetti)
			- check all array and if free is used everytime is possible	(check all cname for example)

			- Latex: how to use generic callback (if it's write good he can take the part for next year courses)
	PACCO	- User_cut di Concorde (difficile implementarle e non è possibile utilizzarle ogni volta bensì solo ogni tot..)

			13/07/2020
	RISOLTO	- Aggiungere limite di tempo in tutti i metodi che lo richiedono (io lo metterei in tutti se il timelimit viene settato)
			  e dobbiamo stare attenti, dovremmo tenere un remaining_time globale in modo che se utilizziamo modelli con più metodi scalino tutti dallo stesso remaining_time

			- Tabu search: we can add variable length of list and add the check also for tabu solution

	RISOLTO - Errore nei modelli 19-20-21-22 da sistemare! RIGA : 883 => CPXaddmipstarts genera l'errore... l'euristica ritorna un tour unico, dopo aver fatto CPXaddmipstarts il tour viene diviso
			  in più tour.. già provato a cambiare i parametri del metodo

			17/08/2020
			- modelli 19-20-21-22 aggiungiamo anche le callback o lasciamo senza?

*/

char * model_name(int i) {
	switch (i) {
		case 0: return "subtour";							// basic model with asymmetric x and q
		case 1: return "mtz";								// MTZ contraints
		case 2: return "flow1_n-2";							// FLOW 1 with y_0j <= x_0j*(n-2) if i != 0
		case 3: return "flow1_n-1";							// FLOW 1 with y_0j <= x_0j*(n-1)
		case 4: return "mtz_lazy";							// MTZ with LAZY
		case 5: return "subtour_heur";						// Subtour with HEUR
		case 6: return "subtour_callback_lazy";				// Subtour_callback_lazy
		case 7: return "subtour_callback_general";			// Subtour_callback_general
		case 8: return "hard_fixing";						// Hard-Fixing
		case 9: return "local_branching";					// Local-Branching
		case 10: return "heuristic_greedy";					// Greedy (no CPLEX)
		case 11: return "heuristic_greedy_cgal";			// Greedy with CGAL (no CPLEX)
		case 12: return "heuristic_grasp";					// GRASP (no CPLEX)
		case 13: return "heuristic_insertion";				// Heuristic Insertion (no CPLEX)
		case 14: return "grasp";							// GRASP + best_two_opt
		case 15: return "patching";							// Patching
		case 16: return "vns";
		case 17: return "n_greedy_tabu_search";				// Greedy + TABU' SEARCH (linked list version)
		case 18: return "n_greedy_tabu_search_array";		// Greedy + TABU' SEARCH (array version)
		case 19: return "heuristic_greedy_cplex";			// Greedy (Warm Start for CPLEX)
		case 20: return "heuristic_greedy_cgal_cplex";		// Greedy (Warm Start for CPLEX)
		case 21: return "heuristic_grasp_cplex";			// GRASP (Warm Start for CPLEX)
		case 22: return "heuristic_insertion_cplex";		// Heuristic Insertion (Warm Start for CPLEX)
		case 23: return "simulating_annealing";				// GRASP + Simulating Annealing
		case 24: return "genetic_algorithm";				// Genetic Algorithm
		case 25: return "greedy_best_two_opt";
		case 26: return "insertion_best_two_opt";
		case 27: return "n_greedy";							// Greedy of a single tour
		case 28: return "n_grasp";							// GRASP N_TIMES
		case 29: return "n_grasp_best_two_opt";
		case 30: return "n_greedy_best_two_opt";
		case 31: return "vns_n_greedy";
		case 32: return "vns_n_grasp";
		default: return "not_supported";
	}
}

char * setup_model(tspinstance* inst) {
/**
NUM			model_type				warm_start					heuristic						mip_opt							callback

 0		build_sym_std												mip_optimization		subtour_iter_opt
 1		  build_mtz					heur_greedy				hard_fixing			subtour_heur_iter_opt				lazy
 2		 build_flow1			heur_greedy_cgal	local_branching				 CPXmipopt					 generic CAND
 3		build_mtz_lazy			heur_grasp				best_two_opt													generic CAND, GLOBAL
 4											heur_insertion				patching
 5																							vns
 6																				tabu_search (LL)
 7																				tabu_search (Ar)
 8																			simulating_annealing
 9																				genetic_algorithm
*/

	inst->callback = -1;		// callback if needed (generic/lazy etc.)
	inst->heuristic = 0;  	// type of optimization
	inst->warm_start = -1;	// constructive heuristic
	inst->mip_opt = -1;			// MIP optimization
	inst->useCplex = 0;			// 1 for method that use CPLEX, otherwise 0

	switch (inst->setup_model) {
		case 0:
			inst->model_type = 0;
			inst->mip_opt = 0;
			inst->useCplex = 1;
			return "subtour";							// loop model for sTSP
		case 1:
			inst->model_type = 1;
			inst->mip_opt = 2;
			inst->useCplex = 1;
			return "mtz";								// MTZ contraints
		case 2:
			inst->model_type = 2;
			inst->mip_opt = 2;
			inst->useCplex = 1;
			return "flow1_n-2";							// FLOW 1 with y_0j <= x_0j*(n-2) if i != 0
		case 3:
			inst->model_type = 2;
			inst->mip_opt = 2;
			inst->useCplex = 1;
			return "flow1_n-1";							// FLOW 1 with y_0j <= x_0j*(n-1)
		case 4:
			inst->model_type = 3;
			inst->mip_opt = 2;
			inst->useCplex = 1;
			return "mtz_lazy";							// MTZ with LAZY constraints
		case 5:
			inst->model_type = 0;
			inst->mip_opt = 1;
			inst->useCplex = 1;
			return "subtour_ffi";						// Subtour with fast first incumb
		case 6:
			inst->model_type = 0;
			inst->callback = 1;
			inst->mip_opt = 2;
			inst->useCplex = 1;
			return "subtour_callback_lazy";				// Subtour_callback_lazy
		case 7:
			inst->model_type = 0;
			inst->callback = 2;
			inst->mip_opt = 2;
			inst->useCplex = 1;
			return "subtour_callback_general";			// Subtour_callback_general
		case 8:
			inst->model_type = 0;
			inst->heuristic = 1;
			inst->callback = 2;
			inst->mip_opt = 2;
			inst->useCplex = 1;
			return "hard_fixing";						// Hard-Fixing
		case 9:
			inst->model_type = 0;
			inst->heuristic = 2;
			inst->callback = 2;
			inst->mip_opt = 2;
			inst->useCplex = 1;
			return "local_branching";					// Soft-Fixing => Local Branching
		case 10:
			inst->model_type = 0;
			inst->warm_start = 1;
			return "heuristic_greedy";					// Heuristic Greedy (no CPLEX)
		case 11:
			inst->model_type = 0;
			inst->warm_start = 2;
			return "heuristic_greedy_cgal";				// Heuristic Greedy CGAL (no CPLEX)
		case 12:
			inst->model_type = 0;
			inst->warm_start = 3;
			return "heuristic_grasp";					// Heuristic GRASP (no CPLEX)
		case 13:
			inst->model_type = 0;
			inst->warm_start = 4;
			return "heuristic_insertion";				// Heuristic Insertion (no CPLEX)
		case 14:
			inst->model_type = 0;
			inst->warm_start = 3;
			inst->heuristic = 3;
			return "grasp_best_two_opt";								// GRASP + best_two_opt
		case 15:
			inst->model_type = 0;
			inst->warm_start = 0;
			inst->heuristic = 4;
			return "patching";							// patching, no warm_start, no CPLEX.
		case 16:
			inst->model_type = 0;
			inst->warm_start = 3;
			inst->heuristic = 5;
			return "vns";	 							// VSN: GRASP + 2opt and random5opt
		case 17:
			inst->model_type = 0;
			inst->warm_start = 5;
			inst->heuristic = 6;
			return "n_greedy_tabu_search";				// N_Greedy + TABU' SEARCH (Linked List version)
		case 18:
			inst->model_type = 0;
			inst->warm_start = 5;
			inst->heuristic = 7;
			return "n_greedy_tabu_search_array";		// N_Greedy + TABU' SEARCH (Array version)
		case 19:
			inst->model_type = 0;
			inst->warm_start = 1;
			inst->mip_opt = 0;
			inst->useCplex = 1;
			return "heuristic_greedy_cplex";					// Heuristic Greedy (Warm Start for CPLEX)
		case 20:
			inst->model_type = 0;
			inst->warm_start = 2;
			inst->mip_opt = 0;
			inst->useCplex = 1;
			return "heuristic_greedy_cgal_cplex";				// Heuristic Greedy CGAL (Warm Start for CPLEX)
		case 21:
			inst->model_type = 0;
			inst->warm_start = 3;
			inst->mip_opt = 0;
			inst->useCplex = 1;
			return "heuristic_grasp_cplex";					// Heuristic GRASP (Warm Start for CPLEX)
		case 22:
			inst->model_type = 0;
			inst->warm_start = 4;
			inst->mip_opt = 0;
			inst->useCplex = 1;
			return "heuristic_insertion_cplex";				// Heuristic Insertion (Warm Start for CPLEX)
		case 23:
			inst->model_type = 0;
			inst->warm_start = 6;
			inst->heuristic = 8;
			return "simulating_annealing";				// n_greedy + Simulating Annealing
		case 24:
			inst->model_type = 0;
			inst->heuristic = 9;
			return "genetic_algorithm";					// Genetic Algorithm
		case 25:
			inst->model_type = 0;
			inst->warm_start = 1;
			inst->heuristic = 3;
			return "greedy_best_two_opt";				// greedy + best_two_opt
		case 26:
			inst->model_type = 0;
			inst->warm_start = 4;
			inst->heuristic = 3;
			return "insertion_best_two_opt";			// Insertion + best_two_opt
		case 27:
			inst->model_type = 0;
			inst->warm_start = 5;
			return "n_greedy";						// Heuristic Greedy N_TIMES
		case 28:
			inst->model_type = 0;
			inst->warm_start = 6;
			return "n_grasp";						// Heuristic GRASP N_TIMES
		case 29:
			inst->model_type = 0;
			inst->warm_start = 6;
			inst->heuristic = 3;
			return "n_grasp_best_two_opt";  // n_grasp_best_two_opt
		case 30:
			inst->model_type = 0;
			inst->warm_start = 5;
			inst->heuristic = 3;
			return "n_greedy_best_two_opt";  // n_greedy_best_two_opt
		case 31:
			inst->model_type = 0;
			inst->warm_start = 5;
			inst->heuristic = 5;
			return "vns_n_greedy";	 							// VSN: n_greedy + 2opt and random_n_opt
		case 32:
			inst->model_type = 0;
			inst->warm_start = 6;
			inst->heuristic = 5;
			return "vns_n_grasp";	 							// VSN: n_grasp + 2opt and random_n_opt
		default: return "not_supported";
	}
}

int TSPopt(tspinstance *inst) {

	inst->init_time = second();  // initial time, used for time measurement

	// set all the parameters of model chosen
	setup_model(inst);
	CPXENVptr env = NULL;
	CPXLPptr lp = NULL;
	int status;
	double init_opt_time;
	srand(inst->randomseed);

	// open cplex model
	if (inst->useCplex) {
		env = CPXopenCPLEX(&status);
		lp = CPXcreateprob(env, &status, "TSP");

		// Cplex's parameter setting
		CPXsetintparam(env,CPX_PARAM_THREADS, inst->nthread);		// allow executing N parallel threads
		CPXsetintparam(env,CPX_PARAM_RANDOMSEED, inst->randomseed);		// avoid performace variability
		CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit);		// set time limit
		if (inst->verbose >= 100) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);	// show CPLEX log

		// CPX_PARAM_MIPEMPHASIS: it balances optimality and integer feasibility.
		//	(CPX_MIPEMPHASIS_BALANCED, CPX_MIPEMPHASIS_FEASIBILITY, CPX_MIPEMPHASIS_OPTIMALITY,
		//	CPX_MIPEMPHASIS_BESTBOUND, CPX_MIPEMPHASIS_HIDDENFEAS)
		// CPX_PARAM_MIPSEARCH: Dynamic search or B&C ?
		// CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 4);			// Display new incumbents, and display a log line every n nodes

		// set input data in CPX structure
		build_model(inst, env, lp);
		if (inst->verbose >= 100) printf("\nbuild model succesfully.\n");
	}

	init_opt_time = second();  // used to calculate optimization time

	// setup struct to save solution
	inst->nedges = inst->model_type == 0 ? inst->nnodes*(inst->nnodes -1)/2 : inst->nnodes*(inst->nnodes -1);
	int best_sol_size = 0;
	if (inst->model_type == 0) {
		best_sol_size = inst->nnodes*(inst->nnodes-1)/2;
	} else if (inst->model_type == 1 || inst->model_type == 3) {
		best_sol_size = inst->nnodes*inst->nnodes;
	} else if (inst->model_type == 2) {
		best_sol_size = inst->nnodes*(inst->nnodes-1)*2;
	} else {
		print_error("invalid model_type");
	}
	inst->best_sol = (double *) calloc(best_sol_size, sizeof(double));
	inst->zbest = CPX_INFBOUND;

	// set callback if selected
	switch_callback(inst, env, lp);

	// set warm start if used
	switch_warm_start(inst, env, lp, &status);


	// compute cplex and calculate opt_time w.r.t. OS used
	if (inst->verbose >= 100) printf("optimizing model...\n");
	optimization(env, lp, inst, &status);
	inst->opt_time = (double)(second() - init_opt_time);
	if (inst->verbose >= 100) printf("optimization complete!\n");

	if(inst->useCplex) {
		CPXsolution(env, lp, &status, &inst->best_lb, inst->best_sol, NULL, NULL, NULL);
		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(&env);

		// free and close cplex model
		if (inst->verbose >= 100) printf("free instance object...\n");
	}

	return 0; // status 0 is ok
}

int xpos(int i, int j, tspinstance *inst) {
	if ( i == j ) print_error(" i == j in xpos" );
	if ( i > j ) return xpos(j,i,inst);								// simplify returned formula
	return i*inst->nnodes + j - ((i + 1)*(i + 2))/2; 				// default case
}

int* invers_xpos(int pos, tspinstance* inst) {
	if (pos < 0) print_error(" position cannot be negative");
	if (pos > (inst->nnodes * inst->nnodes - inst->nnodes - 2) / 2) print_error(" position exceeds max value");

	int* ij = (int*)calloc(2, sizeof(int));
	int temp_pos = 0;
	for (int i = 0; i < inst->nnodes - 1; i++){
		for (int j = i + 1; j < inst->nnodes; j++) {
			if (temp_pos == pos) {
				ij[0] = i;
				ij[1] = j;
				return ij;
			}
			temp_pos++;
		}
	}
	free(ij);
}

int asym_xpos(int i, int j, tspinstance *inst) {
	if ( i == j ) print_error(" i == j in asym_upos" );
	return i*(inst->nnodes - 1) + ( i < j ? j-1 : j );
}

int asym_upos(int i, tspinstance *inst) {
	if ( i < 1 ) print_error(" i < 1 in asym_upos" );
	return inst->nnodes*(inst->nnodes - 1) + i - 1;
}

int asym_ypos(int i, int j, tspinstance* inst) {
	if (i == j) print_error(" i == j in asym_ypos");
	return (inst->nnodes * (inst->nnodes - 1)) + i * (inst->nnodes - 1) + (i < j ? j - 1 : j);
}


// CPLEX: build_model methods add the constraints to OPTIMIZER structures
void build_model(tspinstance *inst, CPXENVptr env, CPXLPptr lp) {

	switch (inst->model_type)
	{
		case 0 :		// basic model with asymmetric x and q => Subtour
			build_sym_std(inst, env,lp);
			break;

		case 1 :		// MTZ contraints
		 	build_mtz(inst, env,lp);
			break;

		case 2 : 		// FLOW 1
			build_flow1(inst, env, lp);
			break;

		case 3: 		// MTZ with LAZY
			build_mtz_lazy(inst, env, lp);
			break;

		default:
			print_error(" model type unknown!!");
			break;
	}

	// save model in .lp file
	if ( inst->verbose >= 1 ){
		char lpname[sizeof(inst->input_file)+20+sizeof(inst->setup_model)];
		char name[sizeof(inst->input_file)];
		snprintf(lpname, sizeof(lpname),
								"%s%c%s_%d.lp",
								"model", DIR_DELIM,
								get_file_name(inst->input_file, name),
								inst->setup_model);  // TODO: input_file check
		// printf("saving %s\n", lpname);
		CPXwriteprob(env, lp, lpname, NULL);
	}
}

void build_sym_std(tspinstance *inst, CPXENVptr env, CPXLPptr lp) {

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

	free(cname[0]);
	free(cname);
}

void build_mtz(tspinstance *inst, CPXENVptr env, CPXLPptr lp) {

	char xctype = 'B';	// type of variable
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
	xctype = 'I';		// maybe not necessary
	obj = 0.0;
	lb = 0.0;
	ub = inst->nnodes - 2.0;

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
	int big_M = inst->nnodes;
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
				rhs = big_M - 1.0;
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

	free(cname[0]);
	free(cname);
}

void build_flow1(tspinstance* inst, CPXENVptr env, CPXLPptr lp) {  // auto adapt for model 2 and 3
	char xctype = CPX_BINARY;		// Binary 0,1
	double obj;						// objective function constant
	double lb;						// lower bound
	double ub,ub1;					// upper bound

	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));			// name of the variable

	// add binary constraints and objective const x(i,j)
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = 0; j < inst->nnodes; j++) {
			if (i != j) {
				sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
				obj = dist(i, j, inst);							// cost == distance
				lb = 0.0;
				ub = i == j ? 0.0 : 1.0;
				if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &xctype, cname))
					print_error(" wrong CPXnewcols on x var.s");
				if (CPXgetnumcols(env, lp) - 1 != asym_xpos(i, j, inst))
					print_error(" wrong position for x var.s");
			}
		}
	}

	// add nodes index 'yij' in the circuits
	xctype = CPX_INTEGER;						// Integer values
	obj = 0.0;
	lb = 0.0;
	ub = (double)inst->nnodes - 1.0;
	ub1 = (double)inst->nnodes - 2.0;

	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = 0; j < inst->nnodes; j++) {
			if (i != j) {
				sprintf(cname[0], "y(%d,%d)", i + 1, j + 1);
				if (j == 0) {														// yi0 = 0
					if (CPXnewcols(env, lp, 1, &obj, &lb, &lb, &xctype, cname))
						print_error(" wrong CPXnewcols on y var.s");
				} else if (i == 0) {													// y0j <= n-1
					if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &xctype, cname))
						print_error(" wrong CPXnewcols on y var.s");
				} else {																// yij <= n-2
					if (CPXnewcols(env, lp, 1, &obj, &lb, &ub1, &xctype, cname))
						print_error(" wrong CPXnewcols on y var.s");
				}

				if (CPXgetnumcols(env, lp) - 1 != asym_ypos(i, j, inst))
					print_error(" wrong position for y var.s");

			}
		}
	}

	// add the degree constraints
	int lastrow;	// the number of rows
	double rhs;		// right head size
	char sense;		// 'L', 'E' or 'G'
	for (int h = 0; h < inst->nnodes; h++) {
		lastrow = CPXgetnumrows(env, lp);
		rhs = 1.0;
		sense = 'E';	// 'E' for equality constraint
		sprintf(cname[0], "income_d(%d)", h + 1);

		// create a new row
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
			print_error(" wrong CPXnewrows [degree]");

		// income vertex
		for (int i = 0; i < inst->nnodes; i++) {
			if (i == h)
				continue;
			if (CPXchgcoef(env, lp, lastrow, asym_xpos(i, h, inst), 1.0)) // income vertex
				print_error(" wrong CPXchgcoef [degree]");
		}

		lastrow = CPXgetnumrows(env, lp);
		sprintf(cname[0], "outcome_d(%d)", h + 1);

		// create a new row
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
			print_error(" wrong CPXnewrows [degree]");

		// outcome vertex
		for (int i = 0; i < inst->nnodes; i++) {
			if (i == h)
				continue;
			if (CPXchgcoef(env, lp, lastrow, asym_xpos(h, i, inst), 1.0)) // outcome vertex
				print_error(" wrong CPXchgcoef [degree]");
		}
	}

	// sum_{h in V} y_1h  = n - 1
	// sum_{j in V} y_hj = sum_{i in V} y_ih -1
	for (int h = 0; h < inst->nnodes; h++) {
		// y constraints: nodes index
		if (h == 0 ) {
			rhs = (double)inst->nnodes - 1.0;
			sense = 'E';
			lastrow = CPXgetnumrows(env, lp);

			sprintf(cname[0], "flow(%d)", h+1);

			// create a new row
			if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))	// new row
				print_error(" wrong CPXnewrows [flow(1)]");

			for (int i = 1; i < inst->nnodes; i++) { 	// sum_{h in V} y_1h  = n - 1

				if (CPXchgcoef(env, lp, lastrow, asym_ypos(h, i, inst), 1.0))		// outcome vertex from 0
					print_error(" wrong CPXchgcoef [flow(1)]");
			}
		} else {
			rhs = 1.0;
			sense = 'E';
			lastrow = CPXgetnumrows(env, lp);

			sprintf(cname[0], "flow(%d)", h+1);

			// create a new row
			if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))	// new row
				print_error(" wrong CPXnewrows [flow(1)]");

			for (int i = 0; i < inst->nnodes; i++) {		// sum_{j in V} y_hj
				if (i != h) {
					if (CPXchgcoef(env, lp, lastrow, asym_ypos(h, i, inst), -1.0))	// outcome
						print_error(" wrong CPXchgcoef [flow(i)]");
				}
			}

			for (int i = 0; i < inst->nnodes; i++) {		// sum_{i in V} y_ih
				if (i != h) {
					if (CPXchgcoef(env, lp, lastrow, asym_ypos(i, h, inst), 1.0))		// income
						print_error(" wrong CPXchgcoef [flow(i)]");
				}
			}
		}
	}

	// y_ij <= x_ij*(n-1)
	rhs = 0.0;
	sense = 'L';
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = 1; j < inst->nnodes; j++) {
			if (i == j) continue;
			lastrow = CPXgetnumrows(env, lp);

			sprintf(cname[0], "y_cut(%d,%d)", i + 1, j + 1);

			// create a new row
			if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))	// new row
				print_error(" wrong CPXnewrows [y_cut()]");

			if (CPXchgcoef(env, lp, lastrow, asym_ypos(i, j, inst), 1.0)) // y_ij
				print_error(" wrong CPXchgcoef [y_cut()]");

			if (i == 0 || inst->setup_model == 3) {		// y_0j <= x_0j*(n-1)
				if (CPXchgcoef(env, lp, lastrow, asym_xpos(i, j, inst), -(double)inst->nnodes + 1.0))
					print_error(" wrong CPXchgcoef [y_cut()]");
			}	else {  																// y_ij <= x_ij*(n-2)
				if (CPXchgcoef(env, lp, lastrow, asym_xpos(i, j, inst), -(double)inst->nnodes + 2.0))
					print_error(" wrong CPXchgcoef [y_cut()]");
			}
		}
	}
}

void build_mtz_lazy(tspinstance* inst, CPXENVptr env, CPXLPptr lp) {

	char xctype = 'B';	// type of variable
	double obj; // objective function constant
	double lb;	// lower bound
	double ub;	// upper bound

	char** cname = (char**)calloc(1, sizeof(char*));	// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));			// name of the variable

	// add binary constraints and objective const x(i,j) for i < j
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = 0; j < inst->nnodes; j++)
		{
			if (i != j)
			{
				sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
				obj = dist(i, j, inst); // cost == distance
				lb = 0.0;
				ub = i == j ? 0.0 : 1.0;
				if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &xctype, cname))
					print_error(" wrong CPXnewcols on x var.s");
				if (CPXgetnumcols(env, lp) - 1 != asym_xpos(i, j, inst))
					print_error(" wrong position for x var.s");
			}
		}
	}

	// add nodes index in the circuits
	xctype = 'I';		// maybe not necessary
	obj = 0.0;
	lb = 0.0;
	ub = inst->nnodes - 2.0;

	for (int i = 1; i < inst->nnodes; i++)
	{
		sprintf(cname[0], "u(%d)", i + 1);
		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &xctype, cname))
			print_error(" wrong CPXnewcols on u var.s");
		if (CPXgetnumcols(env, lp) - 1 != asym_upos(i, inst))
			print_error(" wrong position for u var.s");
	}


	// add the degree constraints
	int lastrow;	// the number of rows
	double rhs;		// right head size
	char sense;		// 'L', 'E' or 'G'
	for (int h = 0; h < inst->nnodes; h++)  // degree constraints
	{
		lastrow = CPXgetnumrows(env, lp);
		rhs = 1.0;
		sense = 'E';	// 'E' for equality constraint
		sprintf(cname[0], "income_d(%d)", h + 1);

		// create a new row
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
			print_error(" wrong CPXnewrows [degree]");

		// income vertex
		for (int i = 0; i < inst->nnodes; i++)
		{
			if (i == h) continue;
			if (CPXchgcoef(env, lp, lastrow, asym_xpos(i, h, inst), 1.0)) // income vertex
				print_error(" wrong CPXchgcoef [degree]");
		}


		// outcome vertex
		lastrow = CPXgetnumrows(env, lp);
		rhs = 1.0;
		sense = 'E';	// 'E' for equality constraint
		sprintf(cname[0], "outcome_d(%d)", h + 1);

		// create a new row
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
			print_error(" wrong CPXnewrows [degree]");
		for (int i = 0; i < inst->nnodes; i++)
		{
			if (i == h) continue;
			if (CPXchgcoef(env, lp, lastrow, asym_xpos(h, i, inst), 1.0)) // income vertex
				print_error(" wrong CPXchgcoef [degree]");
		}

	}

	add_lazy_mtz(inst, env, lp);

	if (inst->verbose >= -100) CPXwriteprob(env, lp, "model/asym_lazy_model.lp", NULL);

	free(cname[0]);
	free(cname);
}

void add_lazy_mtz(tspinstance* inst, CPXENVptr env, CPXLPptr lp) {

	int izero = 0;
	int index[3];
	double value[3];

	char** cname = (char**)calloc(1, sizeof(char*));
	cname[0] = (char*)calloc(100, sizeof(char));

	// add lazy constraints  u_i - u_j + M * x_ij <= M - 1, for each arc (i,j) not touching node 0
	double big_M = inst->nnodes - 1.0;
	double rhs = big_M - 1.0;
	char sense = 'L';
	int nnz = 3;
	for (int i = 1; i < inst->nnodes; i++){
		for (int j = 1; j < inst->nnodes; j++){
			if (i != j) {
				sprintf(cname[0], "u-consistency for arc (%d,%d)", i + 1, j + 1);
				index[0] = asym_upos(i, inst);
				value[0] = 1.0;
				index[1] = asym_upos(j, inst);
				value[1] = -1.0;
				index[2] = asym_xpos(i, j, inst);
				value[2] = big_M;
				if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname))
					print_error("wrong CPXlazyconstraints() for u-consistency");
			}
		}
	}

	// add lazy constraints x_ij + x_ji <= 1, for each arc (i,j) with i < j
	rhs = 1.0;
	nnz = 2;
	for (int i = 0; i < inst->nnodes; i++){
		for (int j = i + 1; j < inst->nnodes; j++){
			sprintf(cname[0], "SEC on node pair (%d,%d)", i + 1, j + 1);
			index[0] = asym_xpos(i, j, inst);
			value[0] = 1.0;
			index[1] = asym_xpos(j, i, inst);
			value[1] = 1.0;
			if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname))
				print_error("wrong CPXlazyconstraints on 2-node SECs");
		}
	}
}

// heuristic warm start
int switch_warm_start(tspinstance* inst, CPXENVptr env, CPXLPptr lp, int* status) {

	switch (inst->warm_start) {

		case -1:
			break;
		case 0:
		break;

		case 1:													// Heuristic greedy
			heur_greedy(env, lp, inst, status);
		break;

		case 2:													// Heuristic greedy CGAL
			heur_greedy_cgal(env, lp, inst, status);
		break;

		case 3:													// Heuristic GRASP
			heur_grasp(inst, status);
		break;

		case 4:													// Heuristic Insertion
			heur_insertion(env, lp, inst, status);
		break;

		case 5:													// Heuristic Greedy Single + best_two_opt
			n_greedy(env, lp, inst, status, (inst->nnodes < 1000) ? 10 : (inst->nnodes > 10000) ? 2 : 5);
		break;

		case 6:													// Heuristic GRASP N_TIMES
			n_grasp(inst, status, (inst->nnodes < 1000) ? 10 : (inst->nnodes > 10000) ? 2 : 5, .95, .03);			// nnodes: < 1000 = 10 times - 1000:5000 = 5 times - >5000 = 2 times
		break;

		default:
			print_error(" model type unknown!!");
		break;
	}
	if (inst->useCplex && inst->warm_start > 0) {

		CPXsetintparam(env, CPX_PARAM_ADVIND, 1);

		/*
			CPX_MIPSTART_REPAIR;			CPLEX attempts to repair the MIP start if it is infeasible, according to
											the parameter that sets the frequency to try to repair an infeasible MIP start
											number of attempts to repair infeasible MIP start = CPXPARAM_MIP_Limits_RepairTries. Default => CPLEX choose
			CPX_MIPSTART_SOLVEMIP;			CPLEX solves a subMIP.
			CPX_MIPSTART_SOLVEFIXED;		CPLEX solves the fixed problem specified by the MIP start (requires to provide values for all discrete variables)

		******** Add a mip start

			if (CPXaddmipstarts(env, lp, 1, inst->nnodes, &izero, best_sol, &val, &nocheck_warmstart, NULL)) {
				print_error("Error during warm start: adding new start, check CPXaddmipstarts\n");
				return status;
			}

			mip_optimization(env, lp, inst, &status);
			CPXsolution(env, lp, &status, &inst->best_lb, inst->best_sol, NULL, NULL, NULL);

		******** Delete a mip start :

			if (CPXdelmipstarts(env, lp, 0, 0)) {
				print_error("Error post warm start: deleting previous start, check CPXdelmipstarts\n");
				return status;
			}

		******** Add multi-mip start:

			int* all_sol = (int*)calloc((int)(inst->nnodes * inst->nnodes), sizeof(int));
			int* all_izero = (int*)calloc(inst->nnodes, sizeof(int));
			for (int j = 0; j < inst->nnodes - 1; j++) {
				all_sol[i * inst->nnodes + j] = xpos(sol[j], sol[j + 1], inst);

			}
			all_sol[(i+1)*inst->nnodes - 1] = xpos(sol[inst->nnodes - 1], i, inst);
			all_izero[i] = i * inst->nnodes;
			}

							   (env, lp, #mip_start, #elem foreach mip_start, start_point, xstar, value, effort, name)
			if (CPXaddmipstarts(env, lp, inst->nnodes, &inst->nnodes, &all_izero, all_sol, &val, &nocheck_warmstart, NULL)) {
				print_error("Error during warm start: adding new start, check CPXaddmipstarts\n");
				return status;
			}

		******** Try if Out of Memory :

			CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);
			CPXsetintparam(env, CPX_PARAM_NODEFILEIND, 3);
			CPXsetintparam(env, CPX_PARAM_WORKMEM, 6144);

		******** Stop ad first incumbent found :

			CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);

		*/
		int nocheck_warmstart = CPX_MIPSTART_NOCHECK;		// CPLEX accepts the MIP start without any checks. The MIP start needs to be complete.

		int* best_sol_complete = (int*)calloc(inst->nedges, sizeof(int));
		for (int i = 0; i < inst->nedges; i++) {
			best_sol_complete[i] = i;
		}

		int izero = 0;

		if (*status = CPXaddmipstarts(env, lp, 1, inst->nedges, &izero, best_sol_complete, inst->best_sol, &nocheck_warmstart, NULL)) {
			print_error("Error during warm start: adding new start, check CPXaddmipstarts\n");
		}
		free(best_sol_complete);
	}

	return *status;

}

int heur_greedy_cgal(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status) {

	#ifdef _WIN32

		double best_lb = CPX_INFBOUND;
		inst->best_lb = CPX_INFBOUND;
		double val = 1.0;
		int* best_sol = (int*)calloc(inst->nnodes, sizeof(int));
		int izero = 0;

		set_verbose(inst->verbose);

		load_point(inst->input_file);

		greedy_alg();

		for (int i = 0; i < inst->nnodes; i++) {
			int* sol = (int*)calloc(inst->nnodes, sizeof(int));
			for (int k = 0; k < inst->nnodes; k++) {
				sol[k] = -1;
			}

			sol = get_greedy_sol(i);

			best_lb = 0.0;
			for (int j = 0; j < inst->nnodes - 1; j++) {
				if(inst->verbose > 100)
					printf("%d,%d\n", sol[j], sol[j + 1]);
				best_lb += dist(sol[j], sol[j + 1], inst);
				sol[j] = xpos(sol[j], sol[j + 1], inst);
			}
			if (inst->verbose > 100)
				printf("%d,%d\n\n", sol[inst->nnodes - 1], i);
			best_lb += dist(sol[inst->nnodes - 1], i, inst);
			sol[inst->nnodes - 1] = xpos(sol[inst->nnodes - 1], i, inst);

			if (inst->verbose >= 100)
				printf("Solution: %d\tBEST_LB found: [%f]\n", i, best_lb);
			if (best_lb < inst->best_lb) {
				if (inst->verbose >= 100)
					printf("BEST_LB update from -> to : [%f] -> [%f]\n", inst->best_lb, best_lb);
				inst->best_lb = best_lb;
				for (int k = 0; k < inst->nnodes; k++)
					best_sol[k] = sol[k];
			}
			//free(sol);
		}
		free_cgal();
		if (inst->verbose >= 10)
			printf("BEST_LB Greedy Heuristic CGAL found: [%f]\n", inst->best_lb);

		// Copy and convert to double sol in best_sol
		for (int i = 0; i < inst->nnodes; i++) {
			inst->best_sol[best_sol[i]] = 1.0;
		}

		return *status;

	#endif

	return *status;
}

int heur_greedy(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status) {

	double best_lb = CPX_INFBOUND;
	double val = 1.0;
	inst->best_lb = CPX_INFBOUND;
	int* best_sol = (int*)calloc(inst->nnodes, sizeof(int));
	int izero = 0;

	for (int i = 0; i < inst->nnodes; i++) {

		int* sol = (int*)calloc(inst->nnodes, sizeof(int));
		for (int k = 0; k < inst->nnodes; k++) {
			sol[k] = -1;
		}

		int succ = succ_not_contained(i, sol, inst);
		sol[0] = i;
		sol[1] = succ;
		if (inst->verbose > 1000)
			printf("%d\n%d\n", sol[0], sol[1]);
		int idx_pred = succ;

		for (int j = 2; j < inst->nnodes; j++) {
			succ = succ_not_contained(idx_pred, sol, inst);
			sol[j] = succ;
			idx_pred = succ;
			if (inst->verbose > 1000)
				printf("%d\n", sol[j]);

		}

		best_lb = 0.0;
		for (int j = 0; j < inst->nnodes - 1; j++) {
			best_lb += dist(sol[j], sol[j + 1], inst);
			sol[j] = xpos(sol[j], sol[j + 1], inst);
		}

		best_lb += dist(sol[inst->nnodes - 1], i, inst);
		sol[inst->nnodes - 1] = xpos(sol[inst->nnodes - 1], i, inst);

		if (best_lb < inst->best_lb) {

			inst->best_lb = best_lb;
			for (int k = 0; k < inst->nnodes; k++)
				best_sol[k] = sol[k];

		}
		free(sol);
	}

	if (inst->verbose >= 100)
		printf("BEST_LB Greedy Heuristic found: [%f]\n", inst->best_lb);

	// Copy and convert to double sol in best_sol
	for (int i = 0; i < inst->nnodes; i++) {
		inst->best_sol[best_sol[i]] = 1.0;
	}
	free(best_sol);

	return *status;
}
int n_greedy(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status, int times) {
	if (inst->verbose >= 100) printf("Heuristic GREEDY %d_TIMES\n", times);

	double best_lb = 0.0;  		// the cost of the solution
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	int cur_node;
	int first_node;
	int cur_nearest;  		// index of the nearest from current node
	double cur_dist;  		// distance from current node
	double d_ij;

	for (int i = 0; i < times; i++) {

		int* best_sol = (int*)calloc(inst->nnodes, sizeof(int));  // list index of selected edges (x_ij = 1)
		// get first node, randomly selected
		for (int i = 0; i < inst->nnodes; i++) succ[i] = -1;

		//int cur_node = round(((double)rand() / RAND_MAX) * (inst->nnodes - 1.0));
		cur_node = rand() % inst->nnodes;
		if (inst->verbose >= 100)
			printf("\n%d - Starting Nodes: %d", i + 1, cur_node);

		first_node = cur_node;

		// find the closer node
		best_lb = 0.0;
		for (int tour_length = 0; tour_length < inst->nnodes; tour_length++) {

			if (tour_length != inst->nnodes-1) {

				// initialization
				cur_dist = INT_MAX; cur_nearest = -1;

				// get the closer node
				for (int j = 0; j < inst->nnodes; j++) {	// for each other node
					if (cur_node == j) continue;  					// except cur_node == j
					if (succ[j] != -1) continue;						// the node is in the tour

					// get distances from last added node and j
					d_ij = dist(cur_node, j, inst);

					// check insertion condition
					if (d_ij < cur_dist) {
						cur_dist = d_ij;
						cur_nearest = j;
					}
				}

			} else {   // add the last node
				cur_dist = dist(cur_node, first_node, inst);
				cur_nearest = first_node;
			}

			// set succ
			succ[cur_node] = cur_nearest;
			best_sol[tour_length] = xpos(cur_node, cur_nearest, inst);
			best_lb += cur_dist;
			cur_node = cur_nearest;
		}

		// save the tour and the cost
		if (inst->verbose >= 101) print_succ(succ, inst);

		if (best_lb < inst->best_lb) {
			clear_sol(inst);
			for (int i = 0; i < inst->nnodes; i++) {
				inst->best_sol[best_sol[i]] = 1.0;
			}
			inst->best_lb = best_lb;
		}
		free(best_sol);

		if (inst->verbose >= 100) printf("GREEDY %d_TIMES BEST_LB: %lf\n", times, inst->best_lb);
		fflush(stdout);

	}
	free(succ);

	return *status;
}
int succ_not_contained(int node, int* sol, tspinstance* inst) {
	double d = INT_MAX;
	int succ = -1;
	int contained;

	for (int i = 0; i < inst->nnodes; i++) {
		if (node != i) {
			contained = 0;
			for (int j = 0; j < inst->nnodes; j++) {
				if (i == sol[j]) {
					contained = 1;
					break;
				}
				else if (sol[j] == -1) {
					break;
				}
			}
			if (!contained) {
				if (d > dist(node, i, inst)) {
					d = dist(node, i, inst);
					succ = i;
				}
			}
		}
	}
	return succ;
}

int heur_grasp(tspinstance* inst, int* status){

	if (inst->verbose >= 100) printf("Heuristic GRASP\n");

	int* best_sol = (int*)calloc(inst->nnodes, sizeof(int));  // list index of selected edges (x_ij = 1)
	// get first node, randomly selected
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	for (int i = 0; i < inst->nnodes; i++) succ[i] = -1;


	int cur_node = round(((double)rand() / RAND_MAX) * (inst->nnodes - 1.0));
	int first_node = cur_node;

	// find the three nearest node of each node
	int cur_nearest[3];  		// index of the nearest from current node
	double cur_dist[3];  		// distance from current node
	double best_lb = 0.0;  	// the cost of the solution


	for (int tour_length = 0; tour_length < inst->nnodes; tour_length++) {

		// initialization
		for (int k = 0; k < 3; k++) {
			cur_dist[k] = INT_MAX; cur_nearest[k] = -1;
		}

		// get the three nearest node
		for (int j = 0; j < inst->nnodes; j++) {	// for each other node
			if (cur_node == j) continue;  					// except cur_node == j
			if (succ[j] != -1) continue;						// the node is in the tour

			// get distances of last added node and j
			double d_ij = dist(cur_node, j, inst);

			// check insertion condition
			for (int k = 0; k < 3; k++) {
				if (d_ij < cur_dist[k]) {
					// insert the new candidate in the vector
					for(int l = k; l < 3-1; l++) {
						cur_dist[l+1] = cur_dist[l];
 						cur_nearest[l+1] = cur_nearest[l];
					}
					cur_dist[k] = d_ij;
					cur_nearest[k] = j;
					break;
				}
			}

		}

		// select succ randomly

		int rand_node = rand() % 3;
		int next_node;
		double next_dist;

		if (cur_nearest[0] != -1) {
			next_node = cur_nearest[rand_node] != -1 ? cur_nearest[rand_node] : cur_nearest[0];
			next_dist = cur_nearest[rand_node] != -1 ? cur_dist[rand_node] : cur_dist[0];
		} else {  // last node, close the tour
			next_node = first_node;
			next_dist = dist(cur_node, first_node, inst);
		}
		succ[cur_node] = next_node;
		best_sol[tour_length] = xpos(cur_node, next_node, inst);
		best_lb += next_dist;
		cur_node = next_node;
	}

	// save the tour and the cost
	if (inst->verbose >= 101) print_succ(succ, inst);
	free(succ);
	for (int i = 0; i < inst->nnodes; i++) {
		inst->best_sol[best_sol[i]] = 1.0;
	}
	inst->best_lb = best_lb;
	free(best_sol);
	if (inst->verbose >= 100) printf("GRASP BEST_LB: %lf\n", inst->best_lb);
	fflush(stdout);

	return 0;
}
int n_grasp(tspinstance* inst, int* status, int times, double x1, double x2) {
	/**

		x1: percentage of the times that select the closer node (from 0.0 to 1.0)
		x2: percentage of times that select the second closer node (from 0.0 to 1.0)
		1-x1-x2 = percentage of time that select the third closer node.
	*/

	if (x1 > 1. || x1 < 0. || x2 > 1. || x2 < 0.) print_error("bad input: n_grasp");

	if (inst->verbose >= 100) printf("Heuristic GRASP %d_TIMES\n", times);

	int  temp_node = -1;
	double best_lb = INT_MAX;  		// the cost of the solution
	inst->best_lb = INT_MAX;
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	int cur_node;
	int first_node;
	int cur_nearest[3];  		// index of the nearest from current node
	double cur_dist[3];  		// distance from current node
	double d_ij;
	int rand_node;
	double rand_double;
	int next_node;
	double next_dist;

	for (int i = 0; i < times; i++) {

		int* best_sol = (int*)calloc(inst->nnodes, sizeof(int));  // list index of selected edges (x_ij = 1)
		// get first node, randomly selected
		for (int i = 0; i < inst->nnodes; i++) succ[i] = -1;

		//int cur_node = round(((double)rand() / RAND_MAX) * (inst->nnodes - 1.0));
		cur_node = rand() % inst->nnodes;
		if (inst->verbose >= 100)
			printf("\n%d - Starting Nodes: %d", i + 1, cur_node);

		first_node = cur_node;

		// find the three nearest node of each node
		best_lb = 0.0;
		for (int tour_length = 0; tour_length < inst->nnodes; tour_length++) {

			// initialization
			for (int k = 0; k < 3; k++) { cur_dist[k] = INT_MAX; cur_nearest[k] = -1; }

			// get the three nearest node
			for (int j = 0; j < inst->nnodes; j++) {	// for each other node
				if (cur_node == j) continue;  					// except cur_node == j
				if (succ[j] != -1) continue;						// the node is in the tour

				// get distances of last added node and j
				d_ij = dist(cur_node, j, inst);

				// check insertion condition
				for (int k = 0; k < 3; k++) {
					if (d_ij < cur_dist[k]) {
						// insert the new candidate in the vector
						for (int l = k; l < 3 - 1; l++) {
							cur_dist[l + 1] = cur_dist[l];
							cur_nearest[l + 1] = cur_nearest[l];
						}
						cur_dist[k] = d_ij;
						cur_nearest[k] = j;
						break;
					}
				}

			}

			// select succ randomly
			rand_node = -1; // initialization
			rand_double = (double)rand()/RAND_MAX;
			if (rand_double <= x1) rand_node = 0;
			else if (rand_double <= x2) rand_node = 1;
			else rand_node = 2;

			if (cur_nearest[0] != -1) {
				next_node = cur_nearest[rand_node] != -1 ? cur_nearest[rand_node] : cur_nearest[0];
				next_dist = cur_nearest[rand_node] != -1 ? cur_dist[rand_node] : cur_dist[0];
			}
			else {  // last node, close the tour
				next_node = first_node;
				next_dist = dist(cur_node, first_node, inst);
			}
			succ[cur_node] = next_node;
			best_sol[tour_length] = xpos(cur_node, next_node, inst);
			best_lb += next_dist;
			cur_node = next_node;
		}

		// save the tour and the cost
		if (inst->verbose >= 101) print_succ(succ, inst);

		if (best_lb < inst->best_lb) {
			clear_sol(inst);
			for (int i = 0; i < inst->nnodes; i++) {
				inst->best_sol[best_sol[i]] = 1.0;
			}
			inst->best_lb = best_lb;
		}
		free(best_sol);

		if (inst->verbose >= 100) printf("GRASP %d_TIMES BEST_LB: %lf\n", times, inst->best_lb);
		fflush(stdout);

	}
	free(succ);

	return *status;
}

int heur_insertion(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status) {
	int* best_sol = (int*)calloc(inst->nnodes, sizeof(int));
	for (int k = 0; k < inst->nnodes; k++)
		best_sol[k] = -1;
	inst->best_lb = 0.0;

	// Creare a random 3-vertex circuit
	int rand1 = -1, rand2 = -1, rand3 = -1;

	while (rand1 == -1 || rand2 == -1 || rand3 == -1) {

		if (rand1 == -1) {
			rand1 = rand() % (inst->nnodes - 1);
			if (rand1 == rand2 || rand1 == rand3)
				rand1 = -1;
		}
		if (rand2 == -1) {
			rand2 = rand() % (inst->nnodes - 1);
			if (rand2 == rand1 || rand2 == rand3)
				rand2 = -1;
		}
		if (rand3 == -1) {
			rand3 = rand() % (inst->nnodes - 1);
			if (rand3 == rand1 || rand3 == rand2)
				rand3 = -1;
		}
	}
	best_sol[0] = xpos(rand1, rand2, inst);
	best_sol[1] = xpos(rand2, rand3, inst);
	best_sol[2] = xpos(rand3, rand1, inst);
	if (inst->verbose > 10)
		printf("\nFirst 3 vertices: %d, %d, %d\t- Their side: %d, %d, %d\n", rand1, rand2, rand3, best_sol[0], best_sol[1], best_sol[2]);
	int count_sol = 3;

	for (int i = 0; i < inst->nnodes; i++) {
		if (i != rand1 && i != rand2 && i != rand3) {
			insertion_move(inst, best_sol, count_sol, i);
			count_sol++;
		}
	}

	// Copy and convert to double sol in best_sol
	for (int i = 0; i < inst->nnodes; i++) {
		inst->best_sol[best_sol[i]] = 1.0;
		inst->best_lb += dist(invers_xpos(best_sol[i], inst)[0], invers_xpos(best_sol[i], inst)[1], inst);
	}

	if(inst->useCplex)
		CPXsetintparam(env, CPX_PARAM_ADVIND, 1);

	return *status;
}
int insertion_move(tspinstance* inst, int* best_sol, int count_sol, int vertex) {
	double extra_mileage = INT_MAX;
	double temp_min = INT_MAX;
	int best_i = -1, best_j = -1, pos, replace_pos = -1;


	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = i; j < inst->nnodes; j++) {
			if (i == j)
				continue;

			pos = contained_in_index(best_sol, count_sol, xpos(i, j, inst));
			if (pos != -1) {
				if (inst->verbose > 100)
					printf("Side %d in pos %d is contained in best_sol => calculate extra_mileage...", best_sol[pos], pos);
				temp_min = dist(i, vertex, inst) + dist(vertex, j, inst) - dist(i, j, inst);
				if (temp_min < extra_mileage) {
					extra_mileage = temp_min;
					replace_pos = pos;
					best_i = i;
					best_j = j;
					if (inst->verbose > 100)
						printf("\t\t=>Extra_mileage upload: %f\n", extra_mileage);
				}
				else if (inst->verbose > 100)
					printf("\n");
			}
		}
	}
	if (inst->verbose > 10)
		printf("***** Replace %d ", best_sol[replace_pos]);
	best_sol[replace_pos] = xpos(vertex, best_i, inst);
	best_sol[count_sol] = xpos(best_j, vertex, inst);
	if (inst->verbose > 10)
		printf(" with %d=[%d,%d] - %d=[%d,%d] added in best_sol (actual length %d) *****\n\n", best_sol[replace_pos], vertex, best_i, best_sol[count_sol], vertex, best_j, count_sol + 1);
	return 0;
}
int contained_in_index(int* vector, int count_sol, int elem) {
	for (int i = 0; i < count_sol; i++) {
		if (vector[i] == elem)
			return i;
	}
	return -1;
}


// optimization
void optimization(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status) {
	switch (inst->heuristic){

		case 0:													// No Heuristic used
			*status = mip_optimization(env, lp, inst, status);
		break;

		case 1:													// Hard-Fixing
			*status = hard_fixing(env, lp, inst, status);
		break;

		case 2:													// Local-Branching
			*status = local_branching(env, lp, inst, status);
		break;

		case 3:
			*status = best_two_opt(inst);
		break;

		case 4:
			patching(inst);
		break;

		case 5:
			*status = vns(inst);
		break;

		case 6:
			*status = tabu_search(env, inst, status);
		break;

		case 7:
			*status = tabu_search_array(env, inst, status);
		break;

		case 8:
			*status = simulating_annealing(env, inst, status);
			break;

		case 9:
			*status = genetic_algorithm(env, inst, status);
			break;

		default:
			print_error("model ì_type not implemented in optimization method");
		break;
	}
}

int vns(tspinstance* inst) {

	if (inst->verbose >= 90) printf("VNS\n");
	fflush(stdout);

	// struct to keep solution
	double* best_sol = (double *) calloc(inst->nedges, sizeof(double));

	// copy current solution
	for (int i = 0; i < inst->nedges; i++) best_sol[i] = inst->best_sol[i];
	double best_lb = inst->best_lb;

	// stop condition: time_limit, number of iteration (i)
	int i = 0;  // number of iteration
	int max_iteration = 100;	// max number of iteration before stop
	int n;
	if (inst->nnodes < 100) n = 5;
	else if (inst->nnodes < 1000) n = 7;
	else if (inst->nnodes < 10000) n = 10;
	else if (inst->nnodes < 100000) n = 15;
	else n = 20;
	// int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	// int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	// int* ncomp = (int*)calloc(1, sizeof(int));
	// build_sol(inst, succ, comp, ncomp);
	// if (*ncomp != 1) {print_error("VNS need to start with a single tour solution.");}

	while (inst->timelimit > second() - inst->init_time && i < max_iteration) {  // stop condition

		best_two_opt(inst);

		// save the solution
		if (best_lb > inst->best_lb) {
			for (int i = 0; i < inst->nedges; i++) best_sol[i] = inst->best_sol[i];
			best_lb = inst->best_lb;
			if (inst->verbose >= 90)
				printf("VNS: iteration = %7d, cur best_lb = %10.2lf, time = %10.2lf UPDATE!\n", i, inst->best_lb, second() - inst->init_time);
		} else if (inst->verbose >= 90)
			printf("VNS: iteration = %7d, cur best_lb = %10.2lf, time = %10.2lf\n", i, inst->best_lb, second() - inst->init_time);
		// move to a random 5 opt solution
		random_n_opt(inst, n);
		i++;  // one iteration compleate
	}


	for (int i = 0; i < inst->nedges; i++) inst->best_sol[i] = best_sol[i];
	inst->best_lb = best_lb;
	return 0;
}

void patching(tspinstance* inst) {

	if (inst->verbose >= 100) printf("PATCHING\n"); fflush(stdout);

	// check if current solution has only one tour
	int *succ = (int*) calloc(inst->nnodes, sizeof(int));
	int *comp = (int*) calloc(inst->nnodes, sizeof(int));
	int *ncomp = (int*) calloc(1, sizeof(int));
	build_sol(inst, succ, comp, ncomp);
	if (inst->verbose >= 100) print_succ(succ, inst);
	if (*ncomp == 1) {
		if(inst->verbose > 10) printf("WARNING: solution already has 1 tour, patching has no effect.\n");
		free(succ);
		free(comp);
		free(ncomp);
		return;
	}

	while (*ncomp > 1) {
		single_patch(inst, succ, comp, ncomp);
		if (inst->verbose >= 100) {
			print_succ(succ, inst);
			printf("comp:   "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", comp[i]);
			printf("\n");
		}
		plot_instance(inst);
	}
}

int hard_fixing(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status) {

	double internal_time_limit = 50.0;
	double next_time_limit = internal_time_limit;
	double objval_p = 0.0;
	double gap = 1.0;						// best_lb - objval_p / best_lb
	double fr = 0.9;				// fixing_ratio

	// set next internal time limit
	double remaining_time = inst->timelimit - (second()-inst->init_time);
	if (remaining_time > internal_time_limit*2)
		next_time_limit = internal_time_limit;
	else
		next_time_limit = remaining_time;

	CPXsetdblparam(env, CPX_PARAM_TILIM, next_time_limit);
	CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 10);

	// first optimization step
	mip_optimization(env, lp, inst, status);

	// get first step solution gap
	CPXsolution(env, lp, status, &inst->best_lb, inst->best_sol, NULL, NULL, NULL);
	CPXgetbestobjval(env, lp, &objval_p);
	gap = (inst->best_lb - objval_p) / inst->best_lb;

	// build_sol(inst, succ, comp, ncomp);
	// if (inst->verbose >= 100) printf("Partial solution, ncomp = %d\n",*ncomp );
	if (inst->verbose >= 1000) plot_instance(inst);
	CPXsetintparam(env, CPX_PARAM_INTSOLLIM, INT_MAX);


	while ( (second()-inst->init_time) < inst->timelimit &&
	 (fr > 0.0 || gap > (inst->optimal_gap)))
	{
		// printf("BEST SOLUTION: %lf\n", inst->best_lb);

		// fix a % of bounds
		fix_bound(env, lp, inst, status, fr);

		// set next internal time limit
		remaining_time = inst->timelimit - (second()-inst->init_time);
		if (remaining_time > internal_time_limit*2)
			next_time_limit = internal_time_limit;
		else
			next_time_limit = remaining_time;
		CPXsetdblparam(env, CPX_PARAM_TILIM, next_time_limit);

		// run again with fixed bound
		mip_optimization(env, lp, inst, status); // this might not be the best solution

		// update gap
		CPXsolution(env, lp, status, &(inst->best_lb), inst->best_sol, NULL, NULL, NULL);
		CPXgetbestobjval(env, lp, &objval_p);
		gap = (inst->best_lb - objval_p) / (inst->best_lb);

		// update fixing_ratio
		if (gap < inst->good_gap) // the solution is good, relax fixing_ration
			fr = fr - inst->decr_fr >= 0.0 ? fr - inst->decr_fr : 0.0;
		else	// the solution is not good, increase frn to get
			fr = fr + inst->incr_fr <= inst->max_fr ? fr + inst->incr_fr : inst->max_fr;

		// build_sol(inst, succ, comp, ncomp);
		// if (inst->verbose >= 100) printf("Partial solution, ncomp = %d\n", *ncomp);
		if (inst->verbose >= 1000) plot_instance(inst);
		if (inst->verbose >= 80) 
			printf("%10.1lf,%.3lf\n", inst->best_lb, second() - inst->init_time);  // used to plot time vs cost
	}
	// if (inst->verbose >= 100) printf("best solution found. ncomp = %d\n", *ncomp);
	// free(succ);
	// free(comp);
	// free(ncomp);
	return 0;
}
void fix_bound(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status, double fixing_ratio) {

	if (inst->verbose >= 100) printf("FIXING %.5lf %%\n", fixing_ratio);

	int k = 0; // position of the ij-th arch in best_sol: xpos(i, j, inst);
	int* indices = (int*) calloc(inst->nnodes,sizeof(int));
	char* lu = (char*)calloc(inst->nnodes, sizeof(char));
	double* bd = (double*)calloc(inst->nnodes, sizeof(double));

	double random;
	int cnt = 0;  // from 0 to inst->nnodes = nedges

	switch (inst->model_type) {
		case -1:
		  break;
		case 0 :
			for (int i = 0; i < inst->nnodes; i++) {
				for (int j = i+1; j < inst->nnodes; j++) {
					k = xpos(i, j, inst);
					if (inst->best_sol[k] > 0.5) {

						random = (double) rand()/RAND_MAX;
						if (random <= fixing_ratio) { 	// fix bound to 1.0
							indices[cnt] = k;
							lu[cnt] = 'B';
							bd[cnt] = 1.0;
						} else {												// relax the bound to be in (0, 1)
							indices[cnt] = k;
							lu[cnt] = 'L';
							bd[cnt] = 0.0;
						}
						cnt++;
					}
				}
			}
			if (cnt > inst->nnodes) print_error("unaspected value of cnt in fix_bound!");
			*status = CPXchgbds(env, lp, cnt, indices, lu, bd);

			break;
		default:
			print_error(" model type unknown!!");
			break;
	}
	free(indices);
	free(lu);
	free(bd);
}

int local_branching(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status) {

	double ini = second();
	int k_index = 0;
	double k[5] = { 3.0, 5.0, 10.0, 15.0, 20.0};
	double timelimit = 300;									// internal timelimit
	double temp_timelimit = timelimit;
	double remaining_time = inst->timelimit;

	CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);			// abort Cplex after the first incument update
	CPXsetdblparam(env, CPX_PARAM_TILIM, remaining_time);

	double* best_sol = (double*)calloc(inst->nedges, sizeof(double));
	double best_lb = CPX_INFBOUND;

	mip_optimization(env, lp, inst, status);
	CPXsolution(env, lp, status, &inst->best_lb, inst->best_sol, NULL, NULL, NULL);

	CPXsetintparam(env, CPX_PARAM_INTSOLLIM, INT_MAX);

	for (int h = 0; remaining_time > 0.0; h++) {

		if( h != 0)
			ini = second();

		int nnz = 0;
		for (int i = 0; i < inst->nnodes * (inst->nnodes - 1) / 2; i++) {
			if (inst->best_sol[i] > 0.5) {
				nnz++;
			}
		}

		int izero = 0;
		int* index = (int*)calloc(nnz, sizeof(int));
		double* value = (double*)calloc(nnz, sizeof(double));

		nnz = 0;

		char** cname = (char**)calloc(1, sizeof(char*));
		cname[0] = (char*)calloc(100, sizeof(char));

		double rhs;
		double gap = 1.0;
		gap = ((inst->best_lb - inst->best_int) / inst->best_lb) * 100;
		if (inst->verbose >= 1) printf("%.1lf,%lf\n", inst->best_lb, second() - inst->init_time);
		if (best_lb > inst->best_lb &&  gap > 0.09) {
			if (inst->verbose >= 100) {
				printf("BEST_LB update from -> to : [%f] -> [%f]\tBEST_INT : %f\tGAP : %f\n",
									best_lb, inst->best_lb, inst->best_int, gap);

			}

			best_lb = inst->best_lb;
			if(k_index < 5)
				rhs = (double)inst->nnodes - k[k_index];
			else
				rhs = (double)inst->nnodes - k[4] > (double)inst->nnodes ? (double)inst->nnodes : k[4];

		} else {
			if (k_index < 4) {
				k_index++;
				rhs = (double)inst->nnodes - k[k_index];
				temp_timelimit = timelimit * (k_index + 1.0);
			} else {
				k[4] = (k[4] * 2.0 > (double)inst->nnodes) ? (double)inst->nnodes : k[4] * 2.0;
				rhs = (double)inst->nnodes - k[4];
			}
		}

		if (temp_timelimit > remaining_time || k[4] == inst->nnodes)
			CPXsetdblparam(env, CPX_PARAM_TILIM, remaining_time);
		else
			CPXsetdblparam(env, CPX_PARAM_TILIM, temp_timelimit);

		sprintf(cname[0], "local-branching constraint, k_index = %f", k_index < 5 ? k[k_index] : k[4]);

		if (inst->verbose >= 100) {
			printf("\n********** Round = %d  -  K = %f  -  Remaining_time = %6.3lf **********\n\n", h, k[k_index], remaining_time);
		}

		for (int i = 0; i < inst->nnodes; i++){
			for (int j = i + 1; j < inst->nnodes; j++) {
				if (inst->best_sol[xpos(i, j, inst)] == 1.0) {
					index[nnz] = xpos(i, j, inst);
					value[nnz] = 1.0;
					nnz++;
				}
			}
		}

		char sense = 'G';
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname))
			print_error("wrong CPXpreaddrows() for adding local-branching constraint\n");

		free(index);
		free(value);

		if (mip_optimization(env, lp, inst, status)) {
			printf("Error in CPXmipopt\n");
		}

		remaining_time -= (second() - ini);

		if (remaining_time <= 0.01) {
			if (CPXgetstat(env, lp) == 101 || CPXgetstat(env, lp) == 102) {
				if (best_lb > inst->best_lb)
					for (int i = 0; i < inst->nedges; i++)
						best_sol[i] = inst->best_sol[i];
				if (inst->verbose >= 1) printf("%.1lf,%lf\n", inst->best_lb, second() - inst->init_time);
				if(inst->verbose >= 100) printf("*** Better soluzion found! ***\n");
			}else {

				if (inst->verbose >= 100)printf("\nUpdate inst->best_sol...: \n");
				for (int i = 0; i < inst->nedges; i++)
					inst->best_sol[i] = best_sol[i];
			}
			if (CPXdelrows(env, lp, CPXgetnumrows(env, lp) - 1, CPXgetnumrows(env, lp) - 1))
				print_error("wrong CPXdelrows() for deleting local-branching constraint\n");
			return 0;
		}

		if (k[4] == inst->nnodes) {
			int status = CPXgetstat(env, lp);
			if (CPXgetstat(env, lp) == 101 || CPXgetstat(env, lp) == 102) {
				if (inst->verbose >= 100)
					printf("CPXgetstat: %s", status == 101 ? "Optimal integer solution found\n" :
														"Optimal sol. within epgap or epagap tolerance found\n");
				return 0;
			}
			if (inst->verbose >= 100)
				printf("CPXgetstat: %s", (status == 107) ? "Time limit exceeded, integer solution exists" : "debug this to see");
		}

		CPXsolution(env, lp, status, &inst->best_lb, inst->best_sol, NULL, NULL, NULL);

		int status = CPXgetstat(env, lp);
		if (inst->verbose >= 100)
			printf("CPXgetstat: %s",
							(status==101) ? "Optimal integer solution found\n" :
							(status==102) ? "Optimal sol. within epgap or epagap tolerance found\n " :
							(status==107) ? "Time limit exceeded, integer solution exists" : "debug this to see");

		if (CPXgetstat(env, lp) == 101 || CPXgetstat(env, lp) == 102) {
			if(inst->verbose >= 100)
				printf("\nNew solution found... update solution: \n");
			for (int i = 0; i < inst->nedges; i++)
				best_sol[i] = inst->best_sol[i];
		}

		if(CPXdelrows(env, lp, CPXgetnumrows(env, lp) - 1, CPXgetnumrows(env, lp) - 1))
			print_error("wrong CPXdelrows() for deleting local-branching constraint\n");

	}

	return 0;
}

int tabu_search(CPXENVptr env, tspinstance* inst, int* status){
	if (inst->verbose >= 100) printf("Tabu Search\n");

	// check if current solution has only one tour
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int* ncomp = (int*)calloc(1, sizeof(int));
	build_sol(inst, succ, comp, ncomp);
	if (*ncomp != 1) {
		printf("Error: tabu_search is called with best_sol with multiple tour!");
		return 1;
	}

	//CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit);
	int* best_sol = (int*)calloc(inst->nnodes, sizeof(int));
	double best_lb = inst->best_lb;
	double best_temp_lb = best_lb;
	double remaining_time = inst->timelimit;
	int isImprovement = 1;
	int best_pre_i = -1;
	int best_pre_j = -1;

	// tabù list (linked list)
	tabu_list* head = NULL;
	tabu_list* head_save = NULL;
	int countListSize = 0;

	for (int times = 0; remaining_time > 0.0 /*&& times < 1000000*/; times++) {
		double ini = second();

		if (inst->verbose > 100) printf("*** Calc new 2_opt sol ***\n");

		// search in 2opt if a better solution is found
		int i = 0;		// start from node 0
		int j = succ[succ[0]];
		int best_i = i;
		int best_j = i;

		double best_improve = 0.0;
		if (!isImprovement)
			best_improve = -CPX_INFBOUND;
		double d_i1_i2;
		double d_j1_j2;
		double d_i1_j1;
		double d_i2_j2;
		double cur_improve;


		// Best 2-opt
		for (int ti = 0; ti < inst->nnodes - 3; ti++) {
			if (!isImprovement && head != NULL && contained_in_posix(&head, xpos(i, succ[i], inst)) != -1 ) {
				i = succ[i];
				j = succ[succ[i]];
			}
			else {
				d_i1_i2 = dist(i, succ[i], inst);
				for (int tj = ti + 2; tj < inst->nnodes; tj++) {			// no 2opt with consequence arches
					if (i == succ[j] || contained_in_posix(&head, xpos(j, succ[j], inst)) != -1) break;								// should happen only when i = 0,
					d_j1_j2 = dist(j, succ[j], inst);
					d_i1_j1 = dist(i, j, inst);
					d_i2_j2 = dist(succ[i], succ[j], inst);

					if (inst->verbose > 100) printf("*** i - j: %d - %d\n", i, j);

					cur_improve = (d_i1_i2 + d_j1_j2) - (d_i1_j1 + d_i2_j2);
					if (cur_improve > best_improve) {						// cross is better
						best_i = i;
						best_j = j;
						best_improve = cur_improve;

						if (inst->verbose > 99) printf("Best_improve: %f\n", best_improve);
					}
					j = succ[j];
				}
				i = succ[i];
				j = succ[succ[i]];
			}
		}
		// if find new best solution change the arches
		if (best_i != best_j) {
			int i2 = succ[best_i];
			int j2 = succ[best_j];
			int pre_node = succ[i2];
			succ[best_i] = best_j;
			succ[i2] = j2;
			int cur_node = i2;
			int suc_node = j2;
			while (cur_node != best_j) {  // reverse the succ
				suc_node = cur_node;
				cur_node = pre_node;
				pre_node = succ[pre_node];
				succ[cur_node] = suc_node;
				if (inst->verbose >= 101) print_succ(succ, inst);
			}
			best_pre_i = best_i;
			best_pre_j = best_j;
		}else {
			if (contained_in_posix(&head, xpos(best_pre_i, best_pre_j, inst)) == -1) {
				push(&head, xpos(best_pre_i, best_pre_j, inst), 1);
				countListSize++;
			}
			if (contained_in_posix(&head, xpos(succ[best_pre_i], succ[best_pre_j], inst)) == -1) {
				push(&head, xpos(succ[best_pre_i], succ[best_pre_j], inst), 1);
				countListSize++;
			}
		}

		if (!isImprovement) {
			if (contained_in_posix(&head, xpos(best_pre_i, best_pre_j, inst)) == -1) {
				push(&head, xpos(best_pre_i, best_pre_j, inst), 1);
				countListSize++;
			}
			if (contained_in_posix(&head, xpos(succ[best_pre_i], succ[best_pre_j], inst)) == -1) {
				push(&head, xpos(succ[best_pre_i], succ[best_pre_j], inst), 1);
				countListSize++;
			}

			while (countListSize - inst->nnodes / 2 > 0) {
				pop_last(head);
				countListSize --;
			}

			if(inst->verbose > 100) print_list(head);
		}

		best_lb -= best_improve;
		if (inst->verbose >= 1) printf("%.1lf,%lf\n", best_lb, second() - inst->init_time);
		if (best_temp_lb == best_lb){
			// Local/Global minimum found
			isImprovement = 0;
			if(inst->verbose > 100) printf("Local (possible global) Minimum found! Start to rise again\n");
			if (best_lb < inst->best_lb) {
				inst->best_lb = best_lb;
				best_temp_lb = best_lb;
				times = 0;
				clear_sol(inst);
				int first_node = 0;
				int second_node = succ[first_node];
				for (int i = 0; i < inst->nnodes; i++) {
					inst->best_sol[xpos(first_node, second_node, inst)] = 1.0;
					first_node = second_node;
					second_node = succ[second_node];
				}
				if (inst->verbose >= 100) printf("BEST_LB GLOBAL update to : [%f]\n", inst->best_lb);
				//plot_instance(inst);

			}

		} else {
			if (best_improve > 0) {
				isImprovement = 1;

				if (inst->verbose >= 100) printf("BEST_LB update from -> to : [%f] -> [%f]\n", best_temp_lb, best_lb);
				best_temp_lb = best_lb;

			}
		}

		//push(&head_save, best_lb, 0);

		remaining_time -= second() - ini;
	}

	if (inst->verbose > 100) print_succ(succ, inst);
	free(succ);
	free(comp);
	free(ncomp);
	free(best_sol);
	if (inst->verbose >= 100) printf("BEST FINAL GLOBAL LB found: [%f]\n", inst->best_lb);

	if (inst->verbose >= 100) write_list_lb(head_save);

}
void push(tabu_list** head_ref, int arc, int isArc) {
	tabu_list* new_node = (tabu_list*)malloc(sizeof(tabu_list));

	new_node->arc = arc;

	new_node->next = *head_ref;
	*head_ref = new_node;

	// used to distinguish tabu_list for arcs and list for LBs found
	// if(isArc)	printf("++++ ADD element %d to Tabu List\n", arc);
}
int pop_first(tabu_list** head) {
	int retval = -1;
	tabu_list* next_node = NULL;

	if (*head == NULL) {
		return retval;
	}

	next_node = (*head)->next;
	retval = (*head)->arc;
	free(*head);
	*head = next_node;

	//printf("---- REMOVE element %d to Tabu List\n", retval);
	return retval;
}
int pop_last(tabu_list* head) {
	int retval = -1;
	/* if there is only one item in the list, remove it */
	if (head->next == NULL) {
		retval = head->arc;
		free(head);
		return retval;
	}

	/* get to the second to last node in the list */
	tabu_list* current = head;
	while (current->next->next != NULL) {
		current = current->next;
	}

	/* now current points to the second to last item of the list, so let's remove current->next */
	retval = current->next->arc;
	free(current->next);
	current->next = NULL;

	//printf("---- REMOVE element %d to Tabu List\n", retval);
	return retval;
}
int remove_by_index(tabu_list** head, int n) {
	int i = 0;
	int retval = -1;
	tabu_list* current = *head;
	tabu_list* temp_node = NULL;

	if (n == 0) {
		return pop_first(head);
	}

	for (i = 0; i < n; i++) {
		if (current->next == NULL) {
			return retval;
		}
		current = current->next;
	}

	temp_node = current->next;
	retval = temp_node->arc;

	current->next = temp_node->next;
	free(temp_node);

	return retval;

}
int contained_in_posix(tabu_list** head, int arc) {
	int retval = -1, i = 0;
	tabu_list* current = *head;
	int found = -1;

	while (current != NULL) {

		found = current->arc;

		if (found == arc) {
			retval = i;
			return retval;
		}

		i++;
		current = current->next;

	}

	return retval;
}
void print_list(tabu_list* head) {
	tabu_list* current = head;
	int printval;
	//printf("--- Tabu list elements: ---\n");
	while (current != NULL) {
		printval = current->arc;
		if(verbose_without_inst) printf("Arc:\t%d\n", printval);
		current = current->next;
	}
}
void delete_list(tabu_list* head, tabu_list** head_ref) {
	tabu_list* current = head, *next = head;

	while (current) {
		next = current->next;
		free(current);
		current = next;
	}

	*head_ref = NULL;
}
void write_list_lb(tabu_list* head) {
	FILE* fptr;

	#ifdef _WIN32
		fptr = fopen("plot\\lb_result.txt", "w");
	#elif __linux__
		fptr = fopen("plot/lb_result.txt", "w");
	#endif

	tabu_list* current = head;
	int printval;
	printf("--- Writing LB list elements ---\n");
	while (current != NULL) {
		printval = current->arc;
		fprintf(fptr, "%d,\n", printval);
		current = current->next;
	}

	fclose(fptr);
}

int tabu_search_array(CPXENVptr env, tspinstance* inst, int* status) {
	if (inst->verbose >= 100) printf("Tabu Search\n");

	// check if current solution has only one tour
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int* ncomp = (int*)calloc(1, sizeof(int));
	build_sol(inst, succ, comp, ncomp);
	if (*ncomp != 1) {
		printf("Error: tabu_search is called with best_sol with multiple tour!");
		return 1;
	}

	//CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit);
	int* best_sol = (int*)calloc(inst->nnodes, sizeof(int));
	double best_lb = inst->best_lb;
	double best_temp_lb = best_lb;
	double remaining_time = inst->timelimit;
	int isImprovement = 1;
	int best_pre_i = -1;
	int best_pre_j = -1;

	// tabù list (array)
	tabu_list* head_save = NULL;

	// Arrotonda per eccesso nnodes / 2. +1 viene usato per tenere conto di dove inserire l'arco proibito
	int tabu_array_size = (inst->nnodes / 2) % 2 == 0 ? (inst->nnodes / 2) + 1 : ((inst->nnodes + 1) / 2) + 1;
	int* tabu_array = (int*)calloc(tabu_array_size, sizeof(int));
	for (int i = 0; i < tabu_array_size; i++)
		tabu_array[i] = -1;



	for (int times = 0; remaining_time > 0.0 && times < 1000000; times++) {
		double ini = second();

		if (inst->verbose > 100) printf("*** Calc new 2_opt sol ***\n");

		// search in 2opt if a better solution is found
		int i = 0;		// start from node 0
		int j = succ[succ[0]];
		int best_i = i;
		int best_j = i;

		double best_improve = 0.0;
		if (!isImprovement)
			best_improve = -CPX_INFBOUND;
		double d_i1_i2;
		double d_j1_j2;
		double d_i1_j1;
		double d_i2_j2;
		double cur_improve;


		// Best 2-opt
		for (int ti = 0; ti < inst->nnodes - 3; ti++) {
			if (!isImprovement && contained_in_posix_array(tabu_array_size, tabu_array, xpos(i, succ[i], inst)) != -1) {
				i = succ[i];
				j = succ[succ[i]];
			}
			else {
				d_i1_i2 = dist(i, succ[i], inst);
				for (int tj = ti + 2; tj < inst->nnodes; tj++) {			// no 2opt with consequence arches
					if (i == succ[j] || contained_in_posix_array(tabu_array_size, tabu_array, xpos(j, succ[j], inst)) != -1) break;								// should happen only when i = 0,
					d_j1_j2 = dist(j, succ[j], inst);
					d_i1_j1 = dist(i, j, inst);
					d_i2_j2 = dist(succ[i], succ[j], inst);

					if (inst->verbose > 100) printf("*** i - j: %d - %d\n", i, j);

					cur_improve = (d_i1_i2 + d_j1_j2) - (d_i1_j1 + d_i2_j2);
					if (cur_improve > best_improve) {						// cross is better
						best_i = i;
						best_j = j;
						best_improve = cur_improve;

						if (inst->verbose > 99) printf("Best_improve: %f\n", best_improve);
					}
					j = succ[j];
				}
				i = succ[i];
				j = succ[succ[i]];
			}
		}
		// if find new best solution change the arches
		if (best_i != best_j) {
			int i2 = succ[best_i];
			int j2 = succ[best_j];
			int pre_node = succ[i2];
			succ[best_i] = best_j;
			succ[i2] = j2;
			int cur_node = i2;
			int suc_node = j2;
			while (cur_node != best_j) {  // reverse the succ
				suc_node = cur_node;
				cur_node = pre_node;
				pre_node = succ[pre_node];
				succ[cur_node] = suc_node;
				if (inst->verbose >= 101) print_succ(succ, inst);
			}
			best_pre_i = best_i;
			best_pre_j = best_j;
		}
		else {
			if (contained_in_posix_array(tabu_array_size, tabu_array, xpos(best_pre_i, best_pre_j, inst)) == -1) {
				for(int i = 0; i < tabu_array_size; i++)
					if (tabu_array[i] == -1) {
						tabu_array[i] = xpos(best_pre_i, best_pre_j, inst);
						if (inst->verbose >= 100) printf("++++ ADD element %d to Tabu List\n", tabu_array[i]);
						if (inst->verbose >= 100) printf("---- REPLACE element %d with -1 to Tabu List\n", tabu_array[(i + 1) == tabu_array_size ? 0 : (i + 1)]);
						tabu_array[(i + 1) == tabu_array_size ? 0 : (i + 1)] = -1;
						break;
					}
			}
			if (contained_in_posix_array(tabu_array_size, tabu_array, xpos(succ[best_pre_i], succ[best_pre_j], inst)) == -1) {
				for (int i = 0; i < tabu_array_size; i++)
					if (tabu_array[i] == -1) {
						tabu_array[i] = xpos(succ[best_pre_i], succ[best_pre_j], inst);
						if (inst->verbose >= 100) printf("++++ ADD element %d to Tabu List\n", tabu_array[i]);
						if (inst->verbose >= 100) printf("---- REPLACE element %d with -1 to Tabu List\n", tabu_array[(i + 1) == tabu_array_size ? 0 : (i + 1)]);
						tabu_array[(i + 1) == tabu_array_size ? 0 : (i + 1)] = -1;
						break;
					}
			}
		}

		if (!isImprovement) {
			if (contained_in_posix_array(tabu_array_size, tabu_array, xpos(best_pre_i, best_pre_j, inst)) == -1) {
				for (int i = 0; i < tabu_array_size; i++)
					if (tabu_array[i] == -1) {
						tabu_array[i] = xpos(best_pre_i, best_pre_j, inst);
						if (inst->verbose >= 100) printf("++++ ADD element %d to Tabu List\n", tabu_array[i]);
						if (inst->verbose >= 100) printf("---- REPLACE element %d with -1 to Tabu List\n", tabu_array[(i + 1) == tabu_array_size ? 0 : (i + 1)]);
						tabu_array[(i + 1) == tabu_array_size ? 0 : (i + 1)] = -1;
						break;
					}
			}
			if (contained_in_posix_array(tabu_array_size, tabu_array, xpos(succ[best_pre_i], succ[best_pre_j], inst)) == -1) {
				for (int i = 0; i < tabu_array_size; i++)
					if (tabu_array[i] == -1) {
						tabu_array[i] = xpos(succ[best_pre_i], succ[best_pre_j], inst);
						if (inst->verbose >= 100) printf("++++ ADD element %d to Tabu List\n", tabu_array[i]);
						if (inst->verbose >= 100) printf("---- REPLACE element %d with -1 to Tabu List\n", tabu_array[(i + 1) == tabu_array_size ? 0 : (i + 1)]);
						tabu_array[(i + 1) == tabu_array_size ? 0 : (i + 1)] = -1;
						break;
					}
			}

			if (inst->verbose > 99) {
				printf("Tabu List: \n");
				for (int k = 0; k < tabu_array_size; k++)
					printf("%d ", tabu_array[k]);
				printf("\n");
			}
		}

		best_lb -= best_improve;
		if (inst->verbose >= 1) printf("%.1lf,%lf\n", best_lb, second() - inst->init_time);
		fflush(stdout);
		if (best_temp_lb == best_lb) {
			// Local/Global minimum found
			isImprovement = 0;
			if (inst->verbose >= 100) printf("Local (possible global) Minimum found! Start to rise again\n");
			if (best_lb < inst->best_lb) {
				inst->best_lb = best_lb;
				best_temp_lb = best_lb;
				times = 0;
				clear_sol(inst);
				int first_node = 0;
				int second_node = succ[first_node];
				for (int i = 0; i < inst->nnodes; i++) {
					inst->best_sol[xpos(first_node, second_node, inst)] = 1.0;
					first_node = second_node;
					second_node = succ[second_node];
				}
				if (inst->verbose >= 100) printf("BEST_LB GLOBAL update to : [%f]\n", inst->best_lb);

			}

		}
		else {
			if (best_improve > 0) {
				isImprovement = 1;

				if (inst->verbose >= 100) printf("BEST_LB update from -> to : [%f] -> [%f]\n", best_temp_lb, best_lb);
				best_temp_lb = best_lb;

			}
		}

		//push(&head_save, best_lb, 0);

		remaining_time -= second() - ini;
	}

	if (inst->verbose >= 100) print_succ(succ, inst);
	free(succ);
	free(comp);
	free(ncomp);
	free(best_sol);
	if (inst->verbose >= 100) printf("BEST FINAL GLOBAL LB found: [%f]\n", inst->best_lb);

	if (inst->verbose >= 10000) write_list_lb(head_save);

}
int contained_in_posix_array(int tabu_array_size, int* tabu_array, int arc) {
	for (int i = 0; i < tabu_array_size; i++)
		if (tabu_array[i] == arc)
			return i;
	return -1;
}

int simulating_annealing(CPXENVptr env, tspinstance* inst, int* status) {

	if (inst->verbose >= 100) printf("Simulating Annealing\n");

	// check if current solution has only one tour
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int* ncomp = (int*)calloc(1, sizeof(int));
	build_sol(inst, succ, comp, ncomp);
	if (*ncomp != 1) {
		printf("Error: simulating_annealing is called with best_sol with multiple tour!");
		return 1;
	}

	//CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit);

	inst->timelimit -= second() - inst->init_time;

	double* best_sol = (double*)calloc(inst->nedges, sizeof(double));
	double* best_global_sol = (double*)calloc(inst->nedges, sizeof(double));
	double best_lb = inst->best_lb;
	double best_global_lb = inst->best_lb;
	double extra_time = inst->timelimit * 5 / 100;
	double remaining_time = inst->timelimit;
	double temp_time = remaining_time;

	double temperature = INT_MAX;
	double temperature_perc = 1.0;
	double decrease = temperature * 99 / 100;	// -1 = 100%
	double maxTemp = decrease == -1 ? (double)INT_MAX + 1e8 : temperature - decrease;

	double delta = -1;
	double delta_perc = 1.0;
	double K = 1.0;					// Boltzmann constant
	double single_step = 0.0;
	double prob = 0.0;

	int dont_print_same_number = 100;		// used to have a slim verbose


	tabu_list* head_save = NULL;

	if (inst->verbose >= 100)
		printf("\nSimulating Annealing uses a maxTemperature with an extra value of\
		 				1e8. This allow to do more iterations when temperature goes to 0.0!\n");

	double max_dist = max_dist_couple_nodes(inst) * 2.0;
	if (max_dist == -1)
		print_error("Probably error during max_dist_couple_nodes excecution\n");


	while (remaining_time > 0.0) {
		double ini = second();

		for (int i = 0; i < inst->nedges; i++)
			if (inst->best_sol[i] == 1.0)
				best_sol[i] = 1.0;
			else
				best_sol[i] = 0.0;

		if (decrease != -1) {
			if (decrease >= temperature)
				temperature = 0;
			else
				temperature -= decrease;
		}
		temperature_perc = temperature / INT_MAX;
		random_two_opt(inst);


		delta = inst->best_lb - best_lb;

		if (delta > 0.0)  {
			delta_perc = delta / max_dist;																	// Can accept bad move

			prob = (double)rand() / (double)RAND_MAX;

			// Used to print some results during run
			if (((int)(temperature_perc * 100) % 10 == 0 && dont_print_same_number != (int)(temperature_perc * 100)) || inst->verbose > 100) {
				dont_print_same_number = (int)(temperature_perc * 100);
				if (inst->verbose >= 100)
					printf("**** Temp_perc, Delta_perc = [%0.2f,%0.4f] --- Probability [%0.3f] >= [%0.3f] ?\n",
									temperature_perc, delta_perc, prob, exp(-delta_perc / (K * temperature_perc)));
			}

			if (prob >= exp(- delta_perc / (K * temperature_perc))) {			// Rejected bad move
				for (int i = 0; i < inst->nedges; i++)
					if (best_sol[i] == 1.0)
						inst->best_sol[i] = 1.0;
					else
						inst->best_sol[i] = 0.0;

				if (inst->verbose > 100) {
					printf("BAD Moves found but refused!!\n");
					printf("LB reset from -> to : [%f] -> [%f]\n", inst->best_lb, best_lb);
				}
				inst->best_lb = best_lb;
				//push(&head_save, inst->best_lb, 0);
				if (inst->verbose >= 1) printf("%.1lf,%lf\n", best_lb, second() - inst->init_time);
			}else {
				if (inst->verbose > 100) {
					printf("BAD Moves found and Accepted\n");													// Accepted bad move
					printf("LB update from -> to : [%f] -> [%f]\n", best_lb, inst->best_lb);
				}
				best_lb = inst->best_lb;
				//push(&head_save, best_lb, 0);
				if (inst->verbose >= 1) printf("%.1lf,%lf\n", best_lb, second() - inst->init_time);
			}
		}else {	 // good move always accepted
			if (inst->verbose > 100) {
				printf("GOOD Moves found and Accepted\n");
				printf("LB update from -> to : [%f] -> [%f]\n", best_lb, inst->best_lb);
			}
			best_lb = inst->best_lb;
			//push(&head_save, best_lb, 0);
			if (inst->verbose >= 1) printf("%.1lf,%lf\n", best_lb, second() - inst->init_time);

			if (inst->best_lb < best_global_lb) {
				if (inst->verbose >= 100) printf("BEST_LB update from -> to : [%f] -> [%f]\n", best_global_lb, inst->best_lb);
				best_global_lb = inst->best_lb;
				for (int i = 0; i < inst->nedges; i++)
					if (inst->best_sol[i] == 1.0)
						best_global_sol[i] = 1.0;
					else
						best_global_sol[i] = 0.0;
			}
		}
		fflush(stdout);

		single_step = (temp_time - remaining_time) / (inst->timelimit);
		decrease = maxTemp * single_step;

		if (inst->verbose > 100)
			printf("Single Step [%f] -- Decrease [%.0f] => Step_time [%f] -- Temperature [%f]\n", single_step, decrease, temp_time - remaining_time, temperature);

		temp_time = remaining_time;
		remaining_time -= second() - ini;
	}

	inst->best_lb = best_global_lb;
	for (int i = 0; i < inst->nedges; i++)
		if (best_global_sol[i] == 1.0)
			inst->best_sol[i] = 1.0;
		else
			inst->best_sol[i] = 0.0;

	//if (inst->verbose > 100) write_list_lb(head_save);

	best_two_opt(inst);

	if (inst->verbose > 100)
		print_succ(succ, inst);

	free(succ);
	free(comp);
	free(ncomp);

	if (inst->verbose >= 100) printf("BEST FINAL GLOBAL LB found: [%f]\n", inst->best_lb);
}
int max_dist_couple_nodes(tspinstance* inst) {
	int max = -1;

	for (int i = 0; i < inst->nnodes; i++)
		for (int j = i; j < inst->nnodes; j++)
			dist(i, j, inst) > max ? max = dist(i, j, inst) : 0;
	return max;
}

int genetic_algorithm(CPXENVptr env, tspinstance* inst, int* status) {

	/*	popolazione di taglia fissa di individui inizializzati random
		repair = salta nodi doppi e fai una costruzione per espansione per i nodi non visitati

	#pragma omp parallel
		{
			printf("\nHello from process : % d\n\n", omp_get_thread_num());
		}
		printf("\nHello from process : % d\n\n", omp_get_thread_num());

	*/

	double global_best_lb = CPX_INFBOUND;
	int updateLB = 0;
	int maxIter = 500;		// 10 * number generation without update LB

	int nPop = 10, nKids = 10;			// GA-EAX/Stage1 & Parallel GA-EAX/Stage1
	int sChunk = 10, nChunk = 30;		// Only for Parallel GA-EAX/Stage1
	int nStag;							// Termination Criterion of GA-EAX/Stage1

	double** population = (double**)calloc(nPop, sizeof(double*));
	double* best_global_sol = (double*)calloc(inst->nedges, sizeof(double));

	int* frequencyTable = (int*)calloc(inst->nedges, sizeof(int));

	int pA, pB;	// Indexes of parents: pA as acceptor and pB as donor

	init_population(inst, population, nPop);
	// save best population LB found
	for (int i = 0; i < nPop; i++) {
		double cost = 0.0;
		for (int j = 0; j < inst->nedges; j++) {
			if (population[i][j] != 0.0) {
				int* ij = invers_xpos(j, inst);
				cost += dist(ij[0], ij[1], inst);
				free(ij);
			}
		}
		if (cost < global_best_lb) {
			global_best_lb = cost;
			for (int k = 0; k < inst->nedges; k++)
				if (population[i][k] == 1.0)
					inst->best_sol[k] = 1.0;
				else
					inst->best_sol[k] = 0.0;
		}
	}
	if (inst->verbose > 1) {
		printf("\n**** INDIVIDUALS POPULATION ****");
		print_population(inst, population, nPop);
	}
	init_frequency_edges(inst, population, frequencyTable, nPop);
	if (inst->verbose > 101) print_frequency_table(inst, frequencyTable);

	double remaining_time = inst->timelimit - (second() - inst->init_time);
	for (int g = 0; updateLB < maxIter && remaining_time > 0.0; g++) {

		double ini = second();

		shuffle_individuals(inst, population, nPop);
		if (inst->verbose > 100) {
			printf("\n**** INDIVIDUALS POPULATION AFTER SHUFFLING ****");
			print_population(inst, population, nPop);
		}

		for (int i = 0; i < nPop && remaining_time > 0.0; i++) {
			pA = i;
			pB = (i == nPop - 1) ? 0 : i + 1;


			double** kids = (double**)calloc(nKids, sizeof(double*));

			int real_nKids = EAX_Single(inst, population, kids, pA, pB, nKids);						// Generate nKids offspring solutions from pA, pB using differents AB-Cycles
			if(real_nKids != 0)
				survival_selection(inst, population, nPop, frequencyTable, real_nKids, pA, kids);

			// free kids matrix
			for (int j = 0; j < real_nKids; j++)
				free(kids[j]);
			free(kids);

			if (inst->verbose > 100) {
				printf("\n**** INDIVIDUALS POPULATION ****");
				print_population(inst, population, nPop);
				if(inst->verbose > 1000) print_frequency_table(inst, frequencyTable);
			}

			// Check if this kid update the global_LB 
			updateLB++;
			double cost = 0.0;
			for (int j = 0; j < inst->nedges; j++) {
				if (population[i][j] != 0.0) {
					int* ij = invers_xpos(j, inst);
					cost += dist(ij[0], ij[1], inst);
					free(ij);
				}
			}
			if (inst->verbose >= 10) printf("Individual %d: %0.f\n", i, cost);
			if (inst->verbose >= 1 && inst->verbose < 10) printf("%.1lf,%lf\n", cost, second() - inst->init_time);
			if (cost < global_best_lb) {
				updateLB = 0;
				global_best_lb = cost;
				for (int k = 0; k < inst->nedges; k++)
					if (population[i][k] == 1.0)
						inst->best_sol[k] = 1.0;
					else
						inst->best_sol[k] = 0.0;
			}

			
		}
		if (inst->verbose > 10) {
			printf("\n**** INDIVIDUALS POPULATION ****");
			print_population(inst, population, nPop);
		}
		//plot_population(inst, population, nPop);
		
		if(inst->verbose >= 10) printf("\n*********** FINISH GENERATION %d ***********\n\n", g);
		
		// print population LBs
		/*
		if (inst->verbose >= 10) printf("\nEvaluate LBs:\n");
		for (int i = 0; i < nPop; i++) {
			double cost = 0.0;
			for (int j = 0; j < inst->nedges; j++) {
				if (population[i][j] != 0.0) {
					int* ij = invers_xpos(j, inst);
					cost += dist(ij[0], ij[1], inst);
					free(ij);
				}
			}
			if (inst->verbose >= 10) printf("Individual %d: %0.f\n", i, cost);
			if (inst->verbose >= 1 && inst->verbose < 10) printf("%.1lf,%lf\n", cost, second() - inst->init_time);
			if (cost < global_best_lb) {
				global_best_lb = cost;
				for (int k = 0; k < inst->nedges; k++)
					if (population[i][k] == 1.0)
						inst->best_sol[k] = 1.0;
					else
						inst->best_sol[k] = 0.0;
			}
		}
		*/

		if (inst->verbose >= 10) {
			printf("\nBEST LB of Population: %.1lf \n", global_best_lb);
			printf("\n");
		}
		fflush(stdout);
		remaining_time -= second() - ini;
	}
	inst->best_lb = global_best_lb;
	free_ga(population, frequencyTable, nPop);


}
void init_population(tspinstance* inst, double** population, int nPop) {
	int status = 1;
	for (int i = 0; i < nPop; i++) {
		n_grasp(inst, &status, (inst->nnodes < 1000) ? 10 : (inst->nnodes > 10000) ? 2 : 5, .33, .33);
		best_two_opt(inst);
		population[i] = (double*)calloc(inst->nedges, sizeof(double));
		for (int j = 0; j < inst->nedges; j++) {
			population[i][j] = inst->best_sol[j];
		}
		clear_sol(inst);
	}
}
void init_frequency_edges(tspinstance* inst, double** population, int* frequencyTable, int nPop) {
	for (int i = 0; i < nPop; i++) {
		for (int j = 0; j < inst->nedges; j++) {
			if (population[i][j] == 1.0)
				frequencyTable[j]++;
		}
	}
}
void shuffle_individuals(tspinstance* inst, double** population, int nPop) {
	for (int i = nPop - 1; i > 0; i--)
	{
		// Pick a random index from 0 to i

		int j = rand() % (i + 1);

		// Swap arr[i] with the element at random index
		swap(&population[i], &population[j]);
	}
}
void swap(double* a, double* b)
{
	double temp = *a;
	*a = *b;
	*b = temp;
}
int EAX_Single(tspinstance* inst, double** population, double** kids, int pA, int pB, int nKids) {
	/*
		- Generate an undirected multigraph defined as GAB = (V; EA U EB).
		- Extract AB-cycles from GAB by repeating a procedure of traversing edges of EA and ones of EB alternatively in
		GAB until an AB-cycle is obtained. Here, an AB-cycle is defined as a cycle in GAB, such that edges of EA and ones of EB are alternately linked.
		- Copy parent pA to intermediate solution y.
		- Choose an AB-cycle from the set of AB-cycles randomly.
		- Generate an intermediate solution by removing the edges of EA in the chosen AB-cycle from the intermediate
		solution y and adding those of EB in the AB-cycle to y.
		- Connect all the subtours in the intermediate solution y to generate an offspring solution y(a single tour) by using a local search heuristic.
	*/
	double** ABcycles = (double**)calloc(nKids, sizeof(double*));			// nKids * 3.0 =  3 * nKids set as maximum number of AB-cycles
	for (int i = 0; i < nKids; i++) {
		ABcycles[i] = (double*)calloc(inst->nedges, sizeof(double));
	}

	double* graph_AB = (double*)calloc(inst->nedges, sizeof(double));
	for (int i = 0; i < inst->nedges; i++) {
		if (population[pA][i] == 1.0)
			graph_AB[i]++;
		if (population[pB][i] == 1.0)
			graph_AB[i]++;
		if (inst->verbose > 1000)
			printf("%.0f %.0f %.0f - ", graph_AB[i], population[pA][i], population[pB][i]);
	}
	if (inst->verbose > 1000) printf("\n");

	// NOOO!! <= nnodes * 2 : perchè potrei fare tutto il giro dei nodi tranne l'ultimo edge che chiude il ciclo e invece tornare indietro per gli stessi e chiuderlo a ritroso
	int** edges_cycles_EA = (int**)calloc(nKids, sizeof(int*));		// Edges di EA che appartengono ai ABcycles
	for (int i = 0; i < nKids; i++) {
		edges_cycles_EA[i] = (int*)calloc(inst->nnodes, sizeof(int));
		for (int j = 0; j < inst->nnodes; j++)
			edges_cycles_EA[i][j] = -1;
	}

	int countCycle = 0;
	extract_ABcycles(inst, population, pA, pB, ABcycles, graph_AB, &countCycle, nKids, edges_cycles_EA);
	if (inst->verbose >= 100) {
		printf("\n**** INDIVIDUALS POPULATION ****");
		print_population(inst, ABcycles, countCycle);
		printf("\n**** EDGES CYCLES EA ****");
		for (int i = 0; i < countCycle; i++) {
			printf("\nCycle EA %d: ", i);
			for (int j = 0; j < inst->nnodes; j++) {
				if (edges_cycles_EA[i][j] != -1)
					printf("%d ", edges_cycles_EA[i][j]);
			}
		}
		printf("\n");
	}
	//plot_population(inst, ABcycles, countCycle);
	//plot_population(inst, population, pB+1);

	int real_nKids = nKids > countCycle ? countCycle : nKids;

	/*
		- We define the size of an AB-cycle as the number of edges of EA (or EB) included in it.
		Note that some of the AB-cycles might consist of two overlapping edges, one from EA and one from EB.
		We call such an ABcycle “ineffective” because the inclusion of ineffective AB-cycles in an E-set does not affect the resulting intermediate solution.
		We call an AB-cycle consisting of more than four edges “effective”
		In Step 3, we select only effective AB-cycles for constructing E-sets unless stated otherwise. We define the size of an E-set as the number of edges of EA (or EB) included in it.
		- According to the definition of an E-set and the procedure in Step 4, EAX generates an intermediate solution
		from EA by replacing edges with the same number of edges selected from EB, under the condition that every vertex is linked by just two edges.
		- An intermediate solution therefore consists of one or more subtours.

		Ec = (Ea \ (E-set intersected EA)) united (E-set intersected Eb)
	*/

	int* idx_effective = (int*)calloc(real_nKids, sizeof(int));
	int i_eff = 0;
	for (int i = 0; i < real_nKids; i++) {

		int* countN = (int*)calloc(inst->nnodes, sizeof(int));
		for (int k = 0; k < inst->nedges; k++) {
			if (ABcycles[i][k] >= 1.0) {
				countN[invers_xpos(k, inst)[0]]++;
				countN[invers_xpos(k, inst)[1]]++;
			}
		}
		int effective = 0;
		for (int k = 0; k < inst->nnodes; k++) {
			if (countN[k] >= 1) {
				effective++;
			}
		}

		free(countN);

		if (effective > 3) {
			idx_effective[i_eff] = i;
			i_eff++;
		}

	}

	int avoid = 0;
	for (int i = 0, temp_i = 0; i < i_eff; i++, temp_i++) {

		int countN = 0;
		if (inst->verbose > 10) {
			printf("\n\n*** Current Effective AB_CYCLE %d: ", i);
			for (int j = 0; j < inst->nedges; j++) {
				if (ABcycles[idx_effective[i]][j] >= 1.0)
					printf("%d ", j);
			}
		}


		double* y = (double*)calloc(inst->nedges, sizeof(double));
		if (inst->verbose > 10) printf("\nEdges of y (EA): \n");
		for (int j = 0; j < inst->nedges; j++) {
			y[j] = population[pA][j];
			if (y[j] == 1.0 && inst->verbose > 10)
				printf("%d ", j);
		}
		if (inst->verbose > 101)
			plot_single(inst, y);

		if (inst->verbose > 10) {
			printf("\nEdges of EB: \n");
			for (int j = 0; j < inst->nedges; j++) {
				if (population[pB][j] == 1.0)
					printf("%d ", j);
			}
		}
		if (inst->verbose > 101)
			plot_single(inst, population[pB]);

		if (inst->verbose > 10) printf("\nEdges of EA removed from y: \n");
		for (int j = 0; j < inst->nnodes; j++)
			if (edges_cycles_EA[idx_effective[i]][j] != -1)
				if (population[pA][edges_cycles_EA[idx_effective[i]][j]] == 1.0)
					if(ABcycles[idx_effective[i]][edges_cycles_EA[idx_effective[i]][j]] >= 1.0) {
						if (inst->verbose > 10) printf("%d ", edges_cycles_EA[idx_effective[i]][j]);
						y[edges_cycles_EA[idx_effective[i]][j]]--;
						countN++;
					}
		if (inst->verbose > 10){
			printf("\nEdges of y After Removing: \n");
			for (int j = 0; j < inst->nedges; j++) {
				if (y[j] == 1.0)
					printf("%d ", j);
			}
		}
		if (inst->verbose > 101)
			plot_single(inst, y);

		if (inst->verbose > 10) printf("\nEdges of EB added to y: \n");
		for (int j = 0; j < inst->nedges; j++) {
			if (population[pB][j] == 1.0 && ABcycles[idx_effective[i]][j] >= 1.0) {
				int found = 0;															// trovato tra gli edges di EA e l'arco non è doppio
				for (int w = 0; w < inst->nnodes; w++) {
					if (edges_cycles_EA[idx_effective[i]][w] == j && ABcycles[idx_effective[i]][j] <= 1.0) {
						found = 1;
						if (inst->verbose > 1000) printf("%d(NO) ", j);
					}
				}
				if (!found) {
					y[j]++;
					countN--;
					if (inst->verbose > 10) printf("%d ", j);
				}
			}
		}
		if (inst->verbose > 10) {
			printf("\nEdges of y After Adding: \n");
			for (int j = 0; j < inst->nedges; j++) {
				if (y[j] == 1.0)
					printf("%d ", j);
				else if (y[j] == 2.0)
					printf("%d(2) ", j);
			}
		}
		if (inst->verbose > 101)
			plot_single(inst, y);



		// Patching
		for (int j = 0; j < inst->nedges; j++) {
			if (y[j] != 0.0 && inst->verbose >= 1000)
				printf("%d <- %d, %d\n", j, invers_xpos(j, inst)[0], invers_xpos(j, inst)[1]);
			inst->best_sol[j] = y[j];
		}

		if (countN != 0) {
			if (inst->verbose >= 1) printf("*********************** Cycle y, before patching, created uncorrecty!! (row 3003) ***********************");
			temp_i--;
			avoid++;
			free(y);
			continue;
		}

		int res = patching_two_edges(inst, y);
		if (res) {
			temp_i--;
			avoid++;
			if (inst->verbose > 10) printf("\nAvoid wrong tours patching!\n");
		}
		else {

			kids[temp_i] = (double*)calloc(inst->nedges, sizeof(double));

			for (int j = 0; j < inst->nedges; j++) {
				if (inst->best_sol[j] != 0.0 && inst->verbose >= 1000)
					printf("%d <- %d, %d\n", j, invers_xpos(j, inst)[0], invers_xpos(j, inst)[1]);
				if (inst->best_sol[j] == 2.0)
					kids[temp_i][j] = 1.0;
				else
					kids[temp_i][j] = inst->best_sol[j];
			}

			//plot_single(inst, kids[i]);
		}


		free(y);
	}

	if (inst->verbose > 100) {
		printf("\nKIDS OBTAINED:");
		print_population(inst, kids, i_eff - avoid);
	}

	for (int i = 0; i < nKids; i++)
		free(ABcycles[i]);
	free(ABcycles);
	free(graph_AB);

	return i_eff - avoid;
}

void extract_ABcycles(tspinstance* inst, double** population, int pA, int pB, double** ABcycles, double* graph_AB, int* idxCycle, int maxNcycles, int** edges_cycles_EA) {
	/*
		- The procedure is started by randomly selecting a vertex.
		- Starting from the selected vertex, trace the edges of EA and EB in GAB in turn until an AB-cycle is found in the traced path,
		where the edge to be traced next is randomly selected (if two candidates exist) and the traced edges are immediately removed from GAB.
		- If an AB-cycle is found in the traced path (a portion of the traced path including the end may form an ABcycle),
		store it and remove the edges constituting it from the traced path.
		- If the current traced path is not empty, start the tracing process again from the end of the current traced path.
		Otherwise, start the tracing process by randomly selecting a vertex from among those linked by at least one edge in GAB.
		- If there is no edge in GAB, iterations of the tracing process are terminated.


	int* succA = (int*)calloc(inst->nnodes, sizeof(int));
	int* prevA = (int*)calloc(inst->nnodes, sizeof(int));
	build_sol_ga(inst, population[pA], succA, prevA);

	int* succB = (int*)calloc(inst->nnodes, sizeof(int));
	int* prevB = (int*)calloc(inst->nnodes, sizeof(int));
	build_sol_ga(inst, population[pB], succB, prevB);

	int init_rand_vertex = -1;
	double* traced_AB = (double*)calloc(inst->nedges, sizeof(double));

	int A_or_B = 1;
	int nextRight = -1, nextLeft = -1, nextEdgeR = -1, nextEdgeL = -1;
	int tourFound = 0;
	int EdgeInGAB = 1;
	while (EdgeInGAB) {

		if (init_rand_vertex == -1) {
			init_rand_vertex = rand() % inst->nnodes;
		}

		// Find next vertex parentA or parentB and so find EA or EB. Remove it from GAB, add it to traced_GAB
		if (A_or_B) {

			nextLeft = prevA[init_rand_vertex];
			nextEdgeL = xpos(nextLeft, init_rand_vertex, inst);

			nextRight = succA[init_rand_vertex];
			nextEdgeR = xpos(init_rand_vertex, nextRight, inst);

			if (graph_AB[nextEdgeL] == 0.0 && graph_AB[nextEdgeR] == 0.0) {
				init_rand_vertex = -1;
			}
			else if (graph_AB[nextEdgeL] == 1.0 && graph_AB[nextEdgeR] == 1.0) {
				if (rand() % 2 == 0) {
					graph_AB[nextEdgeL] = 0.0;
					traced_AB[nextEdgeL] = 1.0;
					init_rand_vertex = nextLeft;
				}
				else {
					graph_AB[nextEdgeR] = 0.0;
					traced_AB[nextEdgeR] = 1.0;
					init_rand_vertex = nextRight;
				}
			}
			else if (graph_AB[nextEdgeL] == 1.0) {
				graph_AB[nextEdgeL] = 0.0;
				traced_AB[nextEdgeL] = 1.0;
				init_rand_vertex = nextLeft;
			}
			else {
				graph_AB[nextEdgeR] = 0.0;
				traced_AB[nextEdgeR] = 1.0;
				init_rand_vertex = nextRight;
			}

			A_or_B = 0;
		}
		else {

			nextLeft = prevB[init_rand_vertex];
			nextEdgeL = xpos(nextLeft, init_rand_vertex, inst);

			nextRight = succB[init_rand_vertex];
			nextEdgeR = xpos(init_rand_vertex, nextRight, inst);

			if (graph_AB[nextEdgeL] == 0.0 && graph_AB[nextEdgeR] == 0.0) {
				init_rand_vertex = -1;
			}
			else if (graph_AB[nextEdgeL] == 1.0 && graph_AB[nextEdgeR] == 1.0) {
				if (rand() % 2 == 0) {
					graph_AB[nextEdgeL] = 0.0;
					traced_AB[nextEdgeL] = 1.0;
					init_rand_vertex = nextLeft;
				}
				else {
					graph_AB[nextEdgeR] = 0.0;
					traced_AB[nextEdgeR] = 1.0;
					init_rand_vertex = nextRight;
				}
			}
			else if (graph_AB[nextEdgeL] == 1.0) {
				graph_AB[nextEdgeL] = 0.0;
				traced_AB[nextEdgeL] = 1.0;
				init_rand_vertex = nextLeft;
			}
			else {
				graph_AB[nextEdgeR] = 0.0;
				traced_AB[nextEdgeR] = 1.0;
				init_rand_vertex = nextRight;
			}

			A_or_B = 1;
		}

		if (init_rand_vertex != -1) {
			evaluate_traced_ABcycle(inst, traced_AB, ABcycles, idxCycle, &tourFound);
			if (tourFound) {
				int* countN = (int*)calloc(inst->nnodes, sizeof(int));
				for (int k = 0; k < inst->nedges; k++) {
					if (traced_AB[k] == 1.0) {
						countN[invers_xpos(k, inst)[0]]++;
						countN[invers_xpos(k, inst)[1]]++;
					}
				}

				if (inst->verbose >= 1000) {
					printf("\ni:   "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", w);
					printf("\nc:   "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", countN[w]);
					printf("\n");
				}

				printf("Node with 1 edge remained in traced path\n");
				int* nodes = (int*)calloc(inst->nnodes, sizeof(int));
				int idx_nodes = 0;
				for (int k = 0; k < inst->nnodes; k++) {
					if (countN[k] == 1) {
						nodes[idx_nodes] = k;
						printf("%d ", nodes[idx_nodes]);
						idx_nodes++;
					}
				}
				printf("\n");

				if (idx_nodes > 0) {
					//init_rand_vertex = nodes[idx_nodes-1];
					init_rand_vertex = rand() % (idx_nodes-1);
					printf("Node choosen: %d\n", init_rand_vertex);
				}else {
					init_rand_vertex = -1;
				}

				free(countN);
				free(nodes);
			}
		}

		printf("\nEdges graph_AB:   ");
		for (int w = 0; w < inst->nedges; w++)
			if (graph_AB[w] == 1.0)
				printf("%2d ", w);
		printf("\n");

		// Termination criterion
		EdgeInGAB = 0;
		for (int w = 0; w < inst->nedges; w++) {
			if (graph_AB[w] == 1.0) {
				EdgeInGAB = 1;
				break;
			}
		}
		if(idxCycle == maxNcycles - 1)
			EdgeInGAB = 1;
	}
	printf("Graph_AB doesn't have any edge or max number [%d] of cycles reatched! ***\n", maxNcycles);
	free(succA);
	free(prevA);
	free(succB);
	free(prevB);
	free(traced_AB);

	*/

	int* succA = (int*)calloc(inst->nnodes, sizeof(int));
	int* prevA = (int*)calloc(inst->nnodes, sizeof(int));
	build_sol_ga(inst, population[pA], succA, prevA, NULL, NULL);
	double* copy_EA = (double*)calloc(inst->nedges, sizeof(double));

	int* succB = (int*)calloc(inst->nnodes, sizeof(int));
	int* prevB = (int*)calloc(inst->nnodes, sizeof(int));
	build_sol_ga(inst, population[pB], succB, prevB, NULL, NULL);
	double* copy_EB = (double*)calloc(inst->nedges, sizeof(double));

	for (int i = 0; i < inst->nedges; i++) {
		copy_EA[i] = population[pA][i];
		copy_EB[i] = population[pB][i];
	}

	int init_rand_vertex = -1;
	double* traced_AB = (double*)calloc(inst->nedges, sizeof(double));

	int A_or_B = 1;
	int nextRight = -1, nextLeft = -1, nextEdgeR = -1, nextEdgeL = -1;
	int tourFound = 0;
	int EdgeInGAB = 1;
	int countEdges = 0;
	int countLoop = 0;

	while (EdgeInGAB && second() - inst->init_time < inst->timelimit) {

		if (init_rand_vertex == -1) {

			init_rand_vertex = rand() % inst->nnodes;
		}

		// Find next vertex parentA or parentB and so find EA or EB. Remove it from GAB, add it to traced_GAB
		if (A_or_B) {

			nextLeft = prevA[init_rand_vertex];
			nextEdgeL = xpos(nextLeft, init_rand_vertex, inst);

			nextRight = succA[init_rand_vertex];
			nextEdgeR = xpos(init_rand_vertex, nextRight, inst);

			if (graph_AB[nextEdgeL] == 0.0 && graph_AB[nextEdgeR] == 0.0) {
				init_rand_vertex = -1;
			}
			else if (graph_AB[nextEdgeL] == 2.0 && graph_AB[nextEdgeR] == 2.0) {

				if (rand() % 2 == 0) {
					graph_AB[nextEdgeL]--;
					traced_AB[nextEdgeL]++;
					edges_cycles_EA[*idxCycle][countEdges] = nextEdgeL;
					countEdges++;
					init_rand_vertex = nextLeft;
					A_or_B = 0;
					copy_EA[nextEdgeL]--;
				}
				else {
					graph_AB[nextEdgeR]--;
					traced_AB[nextEdgeR]++;
					edges_cycles_EA[*idxCycle][countEdges] = nextEdgeR;
					countEdges++;
					init_rand_vertex = nextRight;
					A_or_B = 0;
					copy_EA[nextEdgeR]--;
				}
			}
			else if (graph_AB[nextEdgeL] == 2.0) {
				graph_AB[nextEdgeL]--;
				traced_AB[nextEdgeL]++;
				edges_cycles_EA[*idxCycle][countEdges] = nextEdgeL;
				countEdges++;
				init_rand_vertex = nextLeft;
				A_or_B = 0;
				copy_EA[nextEdgeL]--;
			}
			else if (graph_AB[nextEdgeR] == 2.0) {
				graph_AB[nextEdgeR]--;
				traced_AB[nextEdgeR]++;
				edges_cycles_EA[*idxCycle][countEdges] = nextEdgeR;
				countEdges++;
				init_rand_vertex = nextRight;
				A_or_B = 0;
				copy_EA[nextEdgeR]--;
			}
			else if (graph_AB[nextEdgeL] == 1.0 && copy_EA[nextEdgeL] != 0.0){
				graph_AB[nextEdgeL]--;
				traced_AB[nextEdgeL]++;
				edges_cycles_EA[*idxCycle][countEdges] = nextEdgeL;
				countEdges++;
				init_rand_vertex = nextLeft;
				A_or_B = 0;
				copy_EA[nextEdgeL]--;
			}
			else if (graph_AB[nextEdgeR] == 1.0 && copy_EA[nextEdgeR] != 0.0){
				graph_AB[nextEdgeR]--;
				traced_AB[nextEdgeR]++;
				edges_cycles_EA[*idxCycle][countEdges] = nextEdgeR;
				countEdges++;
				init_rand_vertex = nextRight;
				A_or_B = 0;
				copy_EA[nextEdgeR]--;
			}
			else {
				init_rand_vertex = -1;
			}
		}
		else {

			nextLeft = prevB[init_rand_vertex];
			nextEdgeL = xpos(nextLeft, init_rand_vertex, inst);

			nextRight = succB[init_rand_vertex];
			nextEdgeR = xpos(init_rand_vertex, nextRight, inst);

			if (graph_AB[nextEdgeL] == 0.0 && graph_AB[nextEdgeR] == 0.0) {
				init_rand_vertex = -1;
			}
			else if (graph_AB[nextEdgeL] == 2.0 && graph_AB[nextEdgeR] == 2.0) {

				if (rand() % 2 == 0) {
					graph_AB[nextEdgeL]--;
					traced_AB[nextEdgeL]++;
					init_rand_vertex = nextLeft;
					A_or_B = 1;
					copy_EB[nextEdgeL]--;
				}
				else {
					graph_AB[nextEdgeR]--;
					traced_AB[nextEdgeR]++;
					init_rand_vertex = nextRight;
					A_or_B = 1;
					copy_EB[nextEdgeR]--;
				}
			}
			else if (graph_AB[nextEdgeL] == 2.0) {
				graph_AB[nextEdgeL]--;
				traced_AB[nextEdgeL]++;
				init_rand_vertex = nextLeft;
				A_or_B = 1;
				copy_EB[nextEdgeL]--;
			}
			else if (graph_AB[nextEdgeR] == 2.0) {
				graph_AB[nextEdgeR]--;
				traced_AB[nextEdgeR]++;
				init_rand_vertex = nextRight;
				A_or_B = 1;
				copy_EB[nextEdgeR]--;
			}
			else if (graph_AB[nextEdgeL] == 1.0 && copy_EB[nextEdgeL] != 0.0) {
				graph_AB[nextEdgeL]--;
				traced_AB[nextEdgeL]++;
				init_rand_vertex = nextLeft;
				A_or_B = 1;
				copy_EB[nextEdgeL]--;
			}
			else if (graph_AB[nextEdgeR] == 1.0 && copy_EB[nextEdgeR] != 0.0) {
				graph_AB[nextEdgeR]--;
				traced_AB[nextEdgeR]++;
				init_rand_vertex = nextRight;
				A_or_B = 1;
				copy_EB[nextEdgeR]--;
			}
			else {
				init_rand_vertex = -1;
			}
		}

		if (init_rand_vertex != -1) {
			evaluate_traced_ABcycle(inst, traced_AB, ABcycles, idxCycle, &tourFound, edges_cycles_EA[*idxCycle]);
			if (tourFound) {
				countLoop = 0;
				tourFound = 0;

				int* countN = (int*)calloc(inst->nnodes, sizeof(int));
				for (int k = 0; k < inst->nedges; k++) {
					if (traced_AB[k] == 1.0) {
						countN[invers_xpos(k, inst)[0]]++;
						countN[invers_xpos(k, inst)[1]]++;
					}else if (traced_AB[k] == 2.0) {
						countN[invers_xpos(k, inst)[0]] += 2;
						countN[invers_xpos(k, inst)[1]] += 2;
					}
				}

				if (inst->verbose >= 10000) {
					printf("\ni:   "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", w);
					printf("\nc:   "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", countN[w]);
					printf("\n");
				}

				if (inst->verbose > 10) printf("Node with at least 1 edge remained in traced path\n");
				int* nodes = (int*)calloc(inst->nnodes, sizeof(int));
				int idx_nodes = 0;
				for (int k = 0; k < inst->nnodes; k++) {
					if (countN[k] >= 1) {
						nodes[idx_nodes] = k;
						if(inst->verbose > 10) printf("%d ", nodes[idx_nodes]);
						idx_nodes++;
					}
				}
				if (inst->verbose > 10) printf("\n");

				if (idx_nodes > 0) {
					//init_rand_vertex = nodes[idx_nodes-1];

					init_rand_vertex = idx_nodes - 1 == 0? 0 : nodes[rand() % (idx_nodes - 1)];
					if (inst->verbose > 10) printf("Node choosen: %d\n", init_rand_vertex);
				}
				else {
					init_rand_vertex = -1;
				}

				free(countN);
				free(nodes);

				// remove edges of EA from edges_cycles_EA that do not belong to the corresponding tour of ABcycles
				int temp = -1;
				countEdges = 0;
				for (int w = 0; w < inst->nnodes; w++) {
					if (ABcycles[*idxCycle - 1][edges_cycles_EA[*idxCycle - 1][w]] == 0.0) {
						if(*idxCycle != maxNcycles) 
							edges_cycles_EA[*idxCycle][countEdges] = edges_cycles_EA[*idxCycle - 1][w];
						edges_cycles_EA[*idxCycle - 1][w] = -1;
						countEdges++;
					}
				}
				if (inst->verbose > 10) {
					printf("\nAfter tour removing:\n - ABcycles: %d\n", (*idxCycle) - 1);
					for (int w = 0; w < inst->nedges; w++) {
						if (ABcycles[*idxCycle - 1][w] == 1.0)
							printf("%d ", w);
						else if (ABcycles[*idxCycle - 1][w] == 2.0)
							printf("%d(2) ", w);
					}
					printf("\n - edges_cycles_EA: %d\n", (*idxCycle) - 1);
					for (int w = 0; w < inst->nnodes; w++) {
						if (edges_cycles_EA[*idxCycle - 1][w] != -1 && ABcycles[*idxCycle - 1][edges_cycles_EA[*idxCycle - 1][w]] != 0.0)
							printf("%d ", edges_cycles_EA[*idxCycle - 1][w]);
					}
					if (*idxCycle != maxNcycles) {
						printf("\n => Remaining values on next edges_cycles_EA: %d\n", *idxCycle);
						for (int w = 0; w < inst->nnodes; w++) {
							if (edges_cycles_EA[*idxCycle][w] != -1)
								printf("%d ", edges_cycles_EA[*idxCycle][w]);
						}
					}
					printf("\n");
					/*
					printf("\nedges_cycles_EA: %d\n", (*idxCycle));
					for (int w = 0; w < inst->nnodes; w++) {
						if (edges_cycles_EA[*idxCycle][w] != -1)
							printf("%d ", edges_cycles_EA[*idxCycle][w]);
					}
					*/
				}
			}
		}

		if (inst->verbose > 100) {
			printf("\nEdges graph_AB:   ");
			for (int w = 0; w < inst->nedges; w++) {
				if (graph_AB[w] == 1.0)
					printf("%2d ", w);
				else if (graph_AB[w] == 2.0)
					printf("%2d(%.0f) ", w, graph_AB[w]);
			}
			printf("\nEdges_EA:   ");
			for (int w = 0; w < inst->nedges; w++) {
				if (copy_EA[w] == 1.0)
					printf("%2d ", w);
			}
			printf("\nEdges_EB:   ");
			for (int w = 0; w < inst->nedges; w++) {
				if (copy_EB[w] == 1.0)
					printf("%2d ", w);
			}
		}
		if (inst->verbose > 100) printf("\n********************************* END TOUR SEARCH *********************************\n");

		// Termination criterion
		EdgeInGAB = 0;
		for (int w = 0; w < inst->nedges; w++) {
			if (graph_AB[w] >= 1.0) {
				EdgeInGAB = 1;
				if (countLoop > 2e5) {
					if (inst->verbose > 10) printf("\nSomething goes wrong! Break this search!\n");
					free(succA);
					free(prevA);
					free(succB);
					free(prevB);
					free(traced_AB);
					free(copy_EA);
					free(copy_EB);
					return;
				}
				break;
			}
		}
		if (*idxCycle == maxNcycles)
			EdgeInGAB = 0;

		countLoop++;
	}
	if (inst->verbose > 10) printf("Graph_AB doesn't have any edge or max number [%d] of cycles reatched! Actual: %d ***\n", maxNcycles, *idxCycle);
	free(succA);
	free(prevA);
	free(succB);
	free(prevB);
	free(traced_AB);
	free(copy_EA);
	free(copy_EB);
}
int build_sol_ga(tspinstance* inst, const double* sol, int* succ, int* prev, int* comp, int* ncomp) {

	if (comp != NULL && ncomp != NULL) {

		// initialization of succ, comp and ncomp
		*ncomp = 0;
		for (int i = 0; i < inst->nnodes; i++) {
			succ[i] = -1;
			comp[i] = -1;
		}

		// tour
		for (int start = 0; start < inst->nnodes; start++) {
			if (comp[start] >= 0)
				continue;

			(*ncomp)++;						// the tour id number
			int prv = -1;					// previous node of i in the tour, used to keep j != i
			int i = start;				// iterate over nodes to complete the tour
			int found_succ = 0;		// 1 when found a succ of i
			while (comp[start] == -1) {
				found_succ = 0;
				for (int j = 0; j < inst->nnodes; j++) {	// j iterate to be the subsequent node of i in the tour
					//if (i != j && sol[xpos(i, j, inst)] > 0.5 && j != prv) {  // found subsequent of i (j)
					if (i == j)
						continue;
					int pox = sol[xpos(i, j, inst)];
					if (pox == 1.0 && j != prv) {
						succ[i] = j;
						comp[j] = *ncomp;
						prv = i;
						i = j;
						found_succ = 1;
						break;
					}else if (pox == 2.0 && j != prv) {
						succ[i] = j;
						comp[j] = *ncomp;
						prv = i;
						succ[j] = i;
						comp[i] = *ncomp;
						i = j;
						found_succ = 1;
						break;
					}
				}
				if (!found_succ) {  // no succ found, i is isolated
					if (prv == -1) {	// if no prv found either
						comp[start] = *ncomp;
						break;
					}
					else {  // it's not a tour, i has no subsequent, but is connected to prv
						return 1;
					}
				}
			}
		}

		// print succ and comp
		if (inst->verbose >= 1000) {
			printf("\ni:      "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", i);
			printf("\nsucc:   "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", succ[i]);
			printf("\ncomp:   "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", comp[i]);
			printf("\n");
		}
	}
	else if (prev != NULL) {

		// initialization of succ and prev
		for (int k = 0; k < inst->nnodes; k++) {
			succ[k] = -1;
			prev[k] = -1;
		}

		int* i = (int*)calloc(inst->nnodes, sizeof(int));
		int* j = (int*)calloc(inst->nnodes, sizeof(int));
		int t = 0;

		for (int e = 0; e < inst->nedges; e++) {

			if (sol[e] == 1.0) {
				i[t] = invers_xpos(e, inst)[0];
				j[t] = invers_xpos(e, inst)[1];
				if (inst->verbose > 1000)
					printf("%d <- [%d, %d]\n", e, i[t], j[t]);
				t++;

			}
		}

		if (inst->verbose >= 10000) {
			printf("\i:   "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", i[w]);
			printf("\nj:   "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", j[w]);
			printf("\n");
		}

		int current = i[0];
		succ[current] = j[0];
		prev[succ[current]] = current;
		i[0] = -1;
		j[0] = -1;
		if (inst->verbose >= 10000)
			printf(" prev -> current -> succ: [ %d -> %d -> %d ]\n", prev[current], current, succ[current]);

		for (int count = 0; count < inst->nnodes; count++) {

			current = succ[current];
			for (int n = 0; n < inst->nnodes; n++) {
				if (i[n] == current) {
					succ[current] = j[n];
					prev[succ[current]] = current;
					i[n] = -1;
					j[n] = -1;
					break;
				}
				if (j[n] == current) {
					succ[current] = i[n];
					prev[succ[current]] = current;
					i[n] = -1;
					j[n] = -1;
					break;
				}
			}
			if (inst->verbose >= 10000) {
				printf("\i:   "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", i[w]);
				printf("\nj:   "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", j[w]);
				printf("\n");
				printf(" prev -> current -> succ: [ %d -> %d -> %d ]\n", prev[current], current, succ[current]);
			}
		}

		// print succ and prev
		if (inst->verbose >= 10000) {
			printf("\ni:      "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", i);
			printf("\nsucc:   "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", succ[i]);
			printf("\nprev:   "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", prev[i]);
			printf("\n");
		}

		free(i);
		free(j);
	}
}
void evaluate_traced_ABcycle(tspinstance* inst, double* traced_AB, double** ABcycles, int* idxCycle, int* tourFound, int* edges_cycles_EA_current) {

	/*	Thiella's idea:

		- For each edge in traced path, find relative index i-j
		- Count number of edges for each node
		- Delete all edges that have at least one node with one single edge attached
		- Recursively eliminate all the edges that are conditioned by this
		- If there are remaining arches they will form a cycle. Save and then delete this cycle.

		- By using this sequence of operations, possible errors are avoided (e.g. double cycles)

		=> not works if the condition is to have alternate edges EA and EB :(


	int* i = (int*)calloc(inst->nnodes, sizeof(int));
	int* j = (int*)calloc(inst->nnodes, sizeof(int));
	int* tour = (int*)calloc(inst->nnodes, sizeof(int));
	for (int k = 0; k < inst->nnodes; k++) {
		i[k] = -1;
		j[k] = -1;
		tour[k] = -1;
	}


	for (int e = 0, t = 0; e < inst->nedges; e++) {

		if (traced_AB[e] == 1.0) {
			i[t] = invers_xpos(e, inst)[0];
			j[t] = invers_xpos(e, inst)[1];
			if (inst->verbose >= 100)
				printf("%d <- [%d, %d]\n", e, i[t], j[t]);
			t++;

		}else if (traced_AB[e] == 2.0) {
			i[t] = invers_xpos(e, inst)[0];
			j[t] = invers_xpos(e, inst)[1];
			if (inst->verbose >= 100)
				printf("%d <- [%d, %d]\n", e, i[t], j[t]);
			t++;
			i[t] = invers_xpos(e, inst)[0];
			j[t] = invers_xpos(e, inst)[1];
			if (inst->verbose >= 100)
				printf("%d <- [%d, %d]\n", e, i[t], j[t]);
			t++;

		}
	}

	if (inst->verbose >= 100) {
		printf("i:   "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", i[w]);
		printf("\nj:   "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", j[w]);
		printf("\n");
	}


	int repeat = 0;
	do{
		repeat = 0;
		int* countN = (int*)calloc(inst->nnodes, sizeof(int));
		for (int k = 0; k < inst->nnodes; k++) {
			if (i[k] != -1)
				countN[i[k]]++;
			if (j[k] != -1)
				countN[j[k]]++;
		}

		for (int k = 0; k < inst->nnodes; k++) {
			if (countN[k] == 1) {
				for (int h = 0; h < inst->nnodes; h++) {
					if (i[h] == k) {
						i[h] = -1;
						j[h] = -1;
						repeat = 1;
						break;
					}
					if (j[h] == k) {
						i[h] = -1;
						j[h] = -1;
						repeat = 1;
						break;
					}
				}
			}
		}
		if (inst->verbose >= 1000) {
			printf("\ni:   "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", w);
			printf("\nc:   "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", countN[w]);
			printf("\n");
		}
	} while (repeat);

	// check in i-j if a tour is remained!
	int start = -1;
	int current = -1;
	int next = -1;
	int idx = 0;
	for (int k = 0; k < inst->nnodes; k++) {
		if (i[k] != -1) {
			current = i[k];
			start = current;
			next = j[k];
			i[k] = -1;
			j[k] = -1;
			tour[idx] = current;
			idx++;
			printf("current, next -> [%d, %d]\n", current, next);
			break;
		}
		else if (j[k] != -1) {
			current = j[k];
			start = current;
			next = i[k];
			i[k] = -1;
			j[k] = -1;
			tour[idx] = current;
			idx++;
			printf("current, next -> [%d, %d]\n", current, next);
			break;
		}
	}

	if (inst->verbose >= 1000) {
		printf("\ni:    "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", i[w]);
		printf("\nj:    "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", j[w]);
		printf("\n");
	}

	*tourFound = 0;
	if (current != -1) {
		do {
			current = next;
			tour[idx] = current;
			idx++;
			for (int k = 0; k < inst->nnodes; k++) {
				if (i[k] == current) {
					next = j[k];
					i[k] = -1;
					j[k] = -1;
					printf("current, next -> [%d, %d]\n", current, next);
					break;
				}
				else if (j[k] == current) {
					next = i[k];
					i[k] = -1;
					j[k] = -1;
					printf("current, next -> [%d, %d]\n", current, next);
					break;
				}
			}
			if (inst->verbose >= 1000) {
				printf("\ni:    "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", i[w]);
				printf("\nj:    "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", j[w]);
				printf("\ntour: "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", tour[w]);
				printf("\n");
			}
			if (next == start) {
				*tourFound = 1;
				break;
			}
		} while (current != -1);
	}

	if (inst->verbose >= 100) {
		printf("\ni:    "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", i[w]);
		printf("\nj:    "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", j[w]);
		printf("\ntour: "); for (int w = 0; w < inst->nnodes; w++) printf("%6d", tour[w]);
		printf("\n");
	}

	if (*tourFound) {
		for (int k = 0; k < inst->nnodes; k++) {
			if (tour[k] == -1) {
				break;
			}else{
				int edge = xpos(tour[k], tour[k + 1] == -1 ? tour[0] : tour[k + 1], inst);
				ABcycles[*idxCycle][edge]++;
				traced_AB[edge]--;
			}
		}
	}

	if (inst->verbose >= 99) {
		printf("\nABcycles:   ");
		for (int w = 0; w < inst->nedges; w++) {
			if(ABcycles[*idxCycle][w] == 1.0)
				printf("%2d ", w);
			else if (ABcycles[*idxCycle][w] == 2.0)
				printf("%2d(%.0f) ", w, ABcycles[*idxCycle][w]);
		}
		printf("\ntraced_AB:  ");
		for (int w = 0; w < inst->nedges; w++) {
			if(traced_AB[w] == 1.0)
				printf("%2d ", w);
			else if (traced_AB[w] == 2.0)
				printf("%2d(%.0f) ", w, traced_AB[w]);
		}
		printf("\n");
	}
	if (tour[0] != -1)
		(*idxCycle)++;

	free(i);
	free(j);
	free(tour);
	*/

	int* i = (int*)calloc(inst->nnodes * 2.0, sizeof(int));
	int* j = (int*)calloc(inst->nnodes * 2.0, sizeof(int));
	int* tour = (int*)calloc(inst->nnodes * 2.0, sizeof(int));
	for (int k = 0; k < inst->nnodes * 2; k++) {
		i[k] = -1;
		j[k] = -1;
		tour[k] = -1;
	}

	for (int e = 0, t = 0; e < inst->nedges; e++) {

		if (traced_AB[e] == 1.0) {
			i[t] = invers_xpos(e, inst)[0];
			j[t] = invers_xpos(e, inst)[1];
			if (inst->verbose > 100)
				printf("%d <- [%d, %d]\n", e, i[t], j[t]);
			t++;

		}
		else if (traced_AB[e] == 2.0) {
			i[t] = invers_xpos(e, inst)[0];
			j[t] = invers_xpos(e, inst)[1];
			if (inst->verbose > 100)
				printf("%d <- [%d, %d]\n", e, i[t], j[t]);
			t++;
			i[t] = invers_xpos(e, inst)[0];
			j[t] = invers_xpos(e, inst)[1];
			if (inst->verbose > 100)
				printf("%d <- [%d, %d]\n", e, i[t], j[t]);
			t++;

		}
	}

	if (inst->verbose > 100) {
		printf("\ni:   "); for (int w = 0; w < inst->nnodes * 2; w++) printf("%2d ", i[w]);
		printf("\nj:   "); for (int w = 0; w < inst->nnodes * 2; w++) printf("%2d ", j[w]);
		printf("\n");
	}

	tour_list* tours = NULL;

	tours = grapth_to_tree(inst, i, j, tours, edges_cycles_EA_current);

	if (tours != NULL) {
		*tourFound = 1;

		int size_tours = print_list_of_list(tours);

		int rand_tour = rand() % size_tours;

		int idx_k = 0;
		while (tours != NULL) {

			if (idx_k == rand_tour) {
				tabu_list* current = tours->entry, * del = tours->entry;
				int idx_t = 0;
				while (current != NULL) {
					tour[idx_t] = current->arc;
					current = current->next;
					idx_t++;
				}
				delete_list(del, &del);
				break;
			}

			tours = tours->next;
			idx_k++;
		}

		delete_list(tours, &tours);

		if (inst->verbose > 10) {
			printf("\nVector Tour:\n");
			for (int k = 0; k < inst->nnodes * 2; k++) {
				if(tour[k] != -1)
					printf("%2d ", tour[k]);
			}
			printf("\n");
		}

		for (int k = 0; k < inst->nnodes; k++) {
			if (tour[k] == -1) {
				break;
			}
			else {
				ABcycles[*idxCycle][tour[k]]++;
				traced_AB[tour[k]]--;
			}
		}

		if (inst->verbose >= 100) {
			printf("\nABcycles:   ");
			for (int w = 0; w < inst->nedges; w++) {
				if (ABcycles[*idxCycle][w] == 1.0)
					printf("%2d ", w);
				else if (ABcycles[*idxCycle][w] == 2.0)
					printf("%2d(%.0f) ", w, ABcycles[*idxCycle][w]);
			}
			printf("\ntraced_AB:  ");
			for (int w = 0; w < inst->nedges; w++) {
				if (traced_AB[w] == 1.0)
					printf("%2d ", w);
				else if (traced_AB[w] == 2.0)
					printf("%2d(%.0f) ", w, traced_AB[w]);
			}
			printf("\n");
		}
		if (tour[0] != -1)
			(*idxCycle)++;
	}

	free(i);
	free(j);
	free(tour);


}

tour_list* grapth_to_tree(tspinstance* inst, int* nodes_one, int* nodes_two, tour_list* tours, int* edges_cycles_EA_current) {

	/*
	int root = -1;
	for (int k = 0; k < inst->nnodes * 2; k++) {
		if (nodes_one[k] == -1)
			break;
		else {
			root = rand() % 2 == 0 ? nodes_one[k] : nodes_two[k];
			break;
		}
	}
	if (root != -1) {													// ERROR!! MUST CHECK ALL OF NODES IN CASE WE HAVE MULTI-SUBTOURS!!!!
	*/

	int found = 0;
	for(int start_node = 0; start_node < inst->nnodes; start_node++){

		if (!found) {
			int pred = -1;
			for (int k = 0; k < inst->nnodes * 2; k++) {
				if (!found) {
					if (nodes_one[k] == start_node && nodes_two[k] != pred) {
						int* i = (int*)calloc(inst->nnodes * 2.0, sizeof(int));
						int* j = (int*)calloc(inst->nnodes * 2.0, sizeof(int));
						copy_in_i_j(inst, nodes_one, nodes_two, i, j);

						int current = j[k];
						pred = current;									// if there is, let the other edge from the same nodes free to be able to use it later
						i[k] = -1;
						j[k] = -1;

						if (inst->verbose >= 10000) {
							printf("\ni:   "); for (int w = 0; w < inst->nnodes * 2; w++) printf("%2d ", i[w]);
							printf("\nj:   "); for (int w = 0; w < inst->nnodes * 2; w++) printf("%2d ", j[w]);
							printf("\n");
						}

						tabu_list* pathlist = NULL;
						tabu_list* visited_nodes = NULL;

						push(&pathlist, xpos(start_node, current, inst), 1);
						push(&visited_nodes, start_node, 1);

						if (inst->verbose >= 10000) {
							printf("\n**** Path List *****\n");
							print_list(pathlist);
							printf("**** Visited Node List *****\n");
							print_list(visited_nodes);
							printf("\n");
						}

						tours = Tree_recursive(inst, current, i, j, &found, pathlist, visited_nodes, tours, edges_cycles_EA_current);

						free(i);
						free(j);
						delete_list(pathlist, &pathlist);
						delete_list(visited_nodes, &visited_nodes);

					}
					else if (nodes_two[k] == start_node && nodes_one[k] != pred) {
						int* i = (int*)calloc(inst->nnodes * 2.0, sizeof(int));
						int* j = (int*)calloc(inst->nnodes * 2.0, sizeof(int));
						copy_in_i_j(inst, nodes_one, nodes_two, i, j);

						int current = i[k];
						pred = current;									// if there is, let the other edge from the same nodes free to be able to use it later
						i[k] = -1;
						j[k] = -1;

						if (inst->verbose >= 10000) {
							printf("\ni:   "); for (int w = 0; w < inst->nnodes * 2; w++) printf("%2d ", i[w]);
							printf("\nj:   "); for (int w = 0; w < inst->nnodes * 2; w++) printf("%2d ", j[w]);
							printf("\n");
						}

						tabu_list* pathlist = NULL;
						tabu_list* visited_nodes = NULL;

						push(&pathlist, xpos(start_node, current, inst), 1);
						push(&visited_nodes, start_node, 1);

						if (inst->verbose >= 10000) {
							printf("\n**** Path List *****\n");
							print_list(pathlist);
							printf("**** Visited Node List *****\n");
							print_list(visited_nodes);
							printf("\n");
						}

						tours = Tree_recursive(inst, current, i, j, &found, pathlist, visited_nodes, tours, edges_cycles_EA_current);

						free(i);
						free(j);
						delete_list(pathlist, &pathlist);
						delete_list(visited_nodes, &visited_nodes);

					}
					else if (nodes_one[k] == -1) {
						break;
					}
				}
			}
		}
	}
	return tours;
}
tour_list* Tree_recursive(tspinstance* inst, int current, int* nodes_one, int* nodes_two, int* found, tabu_list* pathlist, tabu_list* visited_nodes, tour_list* tours, int* edges_cycles_EA_current) {


	int pos = contained_in_posix(&visited_nodes, current);
	if (pos == -1) {
		if (second() - inst->init_time > inst->timelimit) {
			if (tours == NULL) {
				return NULL;
			}
			else {
				return tours;
			}
		}
		int pred = -1;
		for (int k = 0; k < inst->nnodes * 2; k++) {
			if (!(*found)) {
				if (nodes_one[k] == current && nodes_two[k] != pred) {
					int* i = (int*)calloc(inst->nnodes * 2.0, sizeof(int));
					int* j = (int*)calloc(inst->nnodes * 2.0, sizeof(int));
					copy_in_i_j(inst, nodes_one, nodes_two, i, j);

					int next = j[k];
					pred = current;
					i[k] = -1;
					j[k] = -1;

					if (inst->verbose >= 10000) {
						printf("\ni:   "); for (int w = 0; w < inst->nnodes * 2; w++) printf("%2d ", i[w]);
						printf("\nj:   "); for (int w = 0; w < inst->nnodes * 2; w++) printf("%2d ", j[w]);
						printf("\n");
					}

					tabu_list* copy_pathlist = NULL;
					tabu_list* copy_visited_nodes = NULL;

					copy_pathlist = copy(pathlist);
					copy_visited_nodes = copy(visited_nodes);

					push(&copy_pathlist, xpos(current, next, inst), 1);
					push(&copy_visited_nodes, current, 1);

					if (inst->verbose >= 10000) {
						printf("\n**** Path List *****\n");
						print_list(copy_pathlist);
						printf("**** Visited Node List *****\n");
						print_list(copy_visited_nodes);
						printf("\n");
					}

					tours = Tree_recursive(inst, next, i, j, found, copy_pathlist, copy_visited_nodes, tours, edges_cycles_EA_current);

					free(i);
					free(j);
					delete_list(copy_pathlist, &copy_pathlist);
					delete_list(copy_visited_nodes, &copy_visited_nodes);

				}
				else if (nodes_two[k] == current && nodes_one[k] != pred) {
					int* i = (int*)calloc(inst->nnodes * 2.0, sizeof(int));
					int* j = (int*)calloc(inst->nnodes * 2.0, sizeof(int));
					copy_in_i_j(inst, nodes_one, nodes_two, i, j);

					int next = i[k];
					i[k] = -1;
					j[k] = -1;

					if (inst->verbose >= 10000) {
						printf("\ni:   "); for (int w = 0; w < inst->nnodes * 2; w++) printf("%2d ", i[w]);
						printf("\nj:   "); for (int w = 0; w < inst->nnodes * 2; w++) printf("%2d ", j[w]);
						printf("\n");
					}

					tabu_list* copy_pathlist = NULL;
					tabu_list* copy_visited_nodes = NULL;

					copy_pathlist = copy(pathlist);
					copy_visited_nodes = copy(visited_nodes);

					push(&copy_pathlist, xpos(current, next, inst), 1);
					push(&copy_visited_nodes, current, 1);

					if (inst->verbose >= 10000) {
						printf("\n**** Path List *****\n");
						print_list(copy_pathlist);
						printf("**** Visited Node List *****\n");
						print_list(copy_visited_nodes);
						printf("\n");
					}

					tours = Tree_recursive(inst, next, i, j, found, copy_pathlist, copy_visited_nodes, tours, edges_cycles_EA_current);

					free(i);
					free(j);
					delete_list(copy_pathlist, &copy_pathlist);
					delete_list(copy_visited_nodes, &copy_visited_nodes);
				}
			}
		}
		if (tours == NULL) {
			return NULL;
		}else {
			return tours;
		}
	}
	else {

		push_list_on_list(inst, &tours, &pathlist, pos, edges_cycles_EA_current);

		if (tours != NULL) {
			(*found) = 1;
		}else
			(*found) = 0;

		return tours;
	}
}
void copy_in_i_j(tspinstance* inst, int* nodes_one, int* nodes_two, int* i, int* j) {
	for (int k = 0; k < inst->nnodes * 2; k++) {
		i[k] = nodes_one[k];
		j[k] = nodes_two[k];
	}
}
void push_list_on_list(tspinstance* inst, tour_list** head_ref, tour_list** pathlist, int pos, int* edges_cycles_EA_current) {

	tabu_list* current_path = *pathlist;
	tabu_list* test_path = *pathlist;

	int* copy_edges_cycles_EA_current = (int*)calloc(inst->nnodes, sizeof(int));

	if(inst->verbose >= 1000)
		printf("\nedges_cycles_EA:\n");
	for (int w = 0; w < inst->nnodes; w++) {
		copy_edges_cycles_EA_current[w] = edges_cycles_EA_current[w];
		if (inst->verbose >= 1000 && copy_edges_cycles_EA_current[w] != -1)
			printf("%d ", copy_edges_cycles_EA_current[w]);
	}
	if (inst->verbose >= 1000)
		printf("\n");

	if (inst->verbose > 1000) printf("\nTOUR to be considered:\n");
	int i = 0, count = 0, pos_EA = 0, alternate = -1, notAlternate = 0;
	while (test_path != NULL) {

		if (i <= pos) {
			count++;

			if (inst->verbose > 1000) printf("%d ", test_path->arc);

			if (alternate == -1) {
				pos_EA = contained_in_posix_array(inst->nnodes, copy_edges_cycles_EA_current, test_path->arc);
				if (pos_EA != -1) {
					alternate = 1;		// contained
					copy_edges_cycles_EA_current[pos_EA] = -1;
				}
				else
					alternate = 0;		// not contained
			}
			else {
				if (alternate) {
					pos_EA = contained_in_posix_array(inst->nnodes, copy_edges_cycles_EA_current, test_path->arc);
					if (pos_EA != -1) {
						notAlternate = 1;
						break;
					}
					else {
						alternate = 0;
					}
				}
				else {
					pos_EA = contained_in_posix_array(inst->nnodes, copy_edges_cycles_EA_current, test_path->arc);
					if (pos_EA != -1) {
						alternate = 1;
						copy_edges_cycles_EA_current[pos_EA] = -1;
					}
					else {
						notAlternate = 1;
						break;
					}
				}
			}
		}
		else
			break;

		test_path = test_path->next;
		i++;
	}
	if (inst->verbose > 1000) printf("\n");

	if (count % 2 != 0 || notAlternate) {							// Odd number of edges certainly not respected the alternating conditions between EA and EB
		if (inst->verbose > 1000) printf("Tour doesn't respect the alternation. Skipped!\n");
		free(copy_edges_cycles_EA_current);
		return;
	}
	if (inst->verbose > 10) printf("\nCYCLE FOUND in current list from pos %d to the last!\n", pos);

	tabu_list* new_tour = NULL;
	tour_list* tours = (tour_list*)malloc(sizeof(tour_list));

	i = 0;
	while (current_path != NULL) {
		if (i <= pos) {
			push(&new_tour, current_path->arc, 1);
			if (inst->verbose > 10) printf("%d ", current_path->arc);
		}
		current_path = current_path->next;
		i++;
	}

	tours->entry = new_tour;

	tours->next = *head_ref;
	(*head_ref) = tours;

	if (verbose_without_inst) {
		printf("\nTour Added:\n");
		print_list(new_tour);
	}
	free(copy_edges_cycles_EA_current);
}
int print_list_of_list(tour_list* tours) {
	tour_list* print_tours = tours;
	int idx_k = 0;
	while (print_tours != NULL) {

		tabu_list* current = print_tours->entry;
		int printval;
		if (verbose_without_inst) printf("\n--- TOUR: %d ---\n", idx_k);
		while (current != NULL) {
			printval = current->arc;
			if (verbose_without_inst) printf("%d ", printval);
			current = current->next;
		}

		print_tours = print_tours->next;
		idx_k++;
	}
	return idx_k;
}
tabu_list* copy(tabu_list* org) {

	tabu_list* new = NULL, ** tail = &new;

	for (; org; org = org->next) {
		*tail = malloc(sizeof * *tail);
		(*tail)->arc = org->arc;
		(*tail)->next = NULL;
		tail = &(*tail)->next;
	}
	return new;
}
int patching_two_edges(tspinstance* inst, double* tour)	 {

	// check if current solution has only one tour
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int* ncomp = (int*)calloc(1, sizeof(int));
	int res = build_sol_ga(inst, tour, succ, NULL, comp, ncomp);
	if (res > 1)
		return 1;
	if (*ncomp == 1) {
		if (inst->verbose > 10) printf("WARNING: solution already has 1 tour, patching has no effect.\n");
		free(succ);
		free(comp);
		free(ncomp);
		return 0;
	}

	while (*ncomp > 1) {
		single_patch(inst, succ, comp, ncomp);
		if (inst->verbose > 101) {
			print_succ(succ, inst);
			if (inst->verbose >= 10000) {
				printf("comp:   "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", comp[i]);
				printf("\n");
			}
			plot_instance(inst);
		}

	}
	free(succ);
	free(comp);
	free(ncomp);
	return 0;
}

void survival_selection(tspinstance* inst, double** population, int nPop, int* frequencyTable, int nKids, int pA, double** kids) {

	double L_parents = calc_L(inst, population, nPop);					// Average Tour Length of Population
	double H_parents = calc_H(inst, frequencyTable, nPop);				// Edge Entropy of Population H = - Σ(e in E) [ F(e) / Npop (log(F(e) / Npop)) ]
	int y = -1;															// index of best kid
	double best_eval = INT_MIN;
	double epsilon = 0.1;												// TODO: best value?

	for (int i = 0; i < nKids; i++) {
		double** new_population = (double**)calloc(nPop, sizeof(double*));
		int* new_frequencyTable = (int*)calloc(inst->nedges, sizeof(int));

		for (int j = 0; j < nPop; j++) {
			new_population[j] = (double*)calloc(inst->nedges, sizeof(double));
			for (int k = 0; k < inst->nedges; k++) {
				new_population[j][k] = population[j][k];
				if (j == 0)
					new_frequencyTable[k] = frequencyTable[k];
			}
		}

		update_frequency_table(inst, new_frequencyTable, new_population[pA], kids[i]);
		for (int j = 0; j < inst->nedges; j++) {
			new_population[pA][j] = kids[i][j];
		}

		double delta_L = L_parents - calc_L(inst, new_population, nPop);
		double delta_H = H_parents - calc_H(inst, new_frequencyTable, nPop);

		double eval = INT_MIN;
		if (delta_L < 0 && delta_H < 0) {
			eval = delta_L / delta_H;
		}
		else if (delta_L < 0 && delta_H >= 0) {
			eval = -delta_L / epsilon;
		}
		else if (delta_L >= 0) {
			eval = -delta_L;
		}
		if (eval > 0 && eval > best_eval) {
			best_eval = eval;
			y = i;
		}

		for (int j = 0; j < nPop; j++)
			free(new_population[j]);
		free(new_population);
		free(new_frequencyTable);
	}

	if (y != -1) {
		if (inst->verbose > 10) printf("\n------ BEST EVAL FOUND: %.2f ------\n", best_eval);
		update_frequency_table(inst, frequencyTable, population[y], kids[y]);
		for (int k = 0; k < inst->nedges; k++) {
			population[y][k] = kids[y][k];
		}
	}
	else {
		if (inst->verbose > 10) printf("\n------ NO BETTER KID FOUND ------\n");
	}
}
void update_frequency_table(tspinstance* inst, int* frequencyTable, double* pA, double* kid) {
	for (int i = 0; i < inst->nedges; i++) {
		if (pA[i] == 1.0)
			frequencyTable[i]--;
		if (kid[i] == 1.0)
			frequencyTable[i]++;
	}
}
double calc_L(tspinstance* inst, double** population, int nPop) {
	double L = 0.0;
	int* ij = (int*)calloc(2, sizeof(int));
	for (int k = 0; k < nPop; k++) {

		double indiv_cost = 0.0;
		int i = 0, j = 0;
		for (int w = 0; w < inst->nedges; w++) {
			if (population[k][w] == 1.0) {
				ij = invers_xpos(w, inst);
				indiv_cost += dist(ij[0], ij[1], inst);
			}
		}
		L += indiv_cost;

	}
	L = L / nPop;
	return L;
}
double calc_H(tspinstance* inst, int* frequencyTable, int nPop){
	double H = 0.0;
	for (int i = 0; i < inst->nedges; i++) {
		if (frequencyTable[i] != 0) {
			H += (double)frequencyTable[i] / nPop * (log((double)frequencyTable[i]) / nPop);
		}
	}
	H *= -1;
	return H;
}

void print_population(tspinstance* inst, double** population, int nPop) {
	if (inst->verbose > 1) {
		for (int i = 0; i < nPop; i++) {
			printf("\nIndividual %d: ", i);
			for (int j = 0; j < inst->nedges; j++) {
				if (population[i][j] == 1.0)
					printf("%d ", j);
				else if (population[i][j] == 2.0)
					printf("%d(%.0f) ", j, population[i][j]);
			}
		}
		printf("\n");
	}
}
void plot_population(tspinstance* inst, double** population, int nPop) {
	printf("\n**** PLOT POPULATION ****\n");
	for (int i = 0; i < nPop; i++) {
		for (int j = 0; j < inst->nedges; j++) {
			inst->best_sol[j] = population[i][j];
		}
		plot_instance(inst);
	}
}
void plot_single(tspinstance* inst, double* individual) {
	printf("\n**** PLOT SINGLE ****\n");
	for (int j = 0; j < inst->nedges; j++) {
		inst->best_sol[j] = individual[j];
	}
	plot_instance(inst);
}
void print_frequency_table(tspinstance* inst, int* frequency_table) {
	if (inst->verbose > 10) {
		printf("\n**** FREQUENCY TABLE ****\n");
		for (int i = 0; i < inst->nedges; i++) {
			printf("Edge - Value: %d - %d\n", i, frequency_table[i]);
		}
	}
}
void free_ga(double** population, int* frequency_table, int nPop) {
	for (int i = 0; i < nPop; i++)
		free(population[i]);
	free(population);
	free(frequency_table);
}

int best_two_opt(tspinstance* inst) {

	double last_best_lb;

	// execute two_opt while local minimum is found
	last_best_lb = inst->best_lb;

	while (inst->timelimit > second() - inst->init_time) {  // iterate two_opt
		two_opt(inst);
		if (inst->verbose >= 80)
			printf("%10.1lf,%.3lf\n", inst->best_lb, second() - inst->init_time);  // used to plot time vs cost
		if (last_best_lb > inst->best_lb) {  // found a new best, continue iteration
			last_best_lb = inst->best_lb;
		} else {		// no new best_lb found, stop iteration
			return 0;
		}
	}
	return 0;
}

// optimization methods run the problem optimization
int mip_optimization(CPXENVptr env, CPXLPptr lp, tspinstance *inst, int *status) {

	switch (inst->mip_opt)
	{
		case -1: // Empty. Used to evaluate heuristic alone
			break;

		case 0:			// subtour_iter_opt, symmetric, without callback
			*status = subtour_iter_opt(env, lp, inst, status);
			break;

		case 1:			// Subtour with HEUR
			*status = subtour_heur_iter_opt(env, lp, inst, status, 0);
			break;

		case 2:			// Simply CPLEX optimization
			*status = CPXmipopt(env,lp);
			break;

		default:
			print_error("model ì_type not implemented in optimization method");
			break;
	}
	return 0;
}

int subtour_iter_opt(CPXENVptr env, CPXLPptr lp, tspinstance *inst, int *status) {

	// structure init
	int *succ = (int*) calloc(inst->nnodes, sizeof(int));
	int *comp = (int*) calloc(inst->nnodes, sizeof(int));
	int *ncomp = (int*) calloc(1, sizeof(int));
	char sense = 'L';
	double rhs;
	int lastrow;

	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));			// name of the variable
	*ncomp = inst->nnodes;
	int subtour_counter = 0;

	CPXmipopt(env,lp);
	CPXsolution(env, lp, status, &inst->best_lb, inst->best_sol, NULL, NULL, NULL);
	build_sol(inst, succ, comp, ncomp);
	if (inst->verbose >= 100) printf("Iter %3d partial solution, ncomp = %d\n", subtour_counter,*ncomp );
	if (inst->verbose >= 1000) plot_instance(inst);
	while (*ncomp >= 2)
	{

		// subtour elimination constraints, one for each component
		for ( int comp_i = 1; comp_i <= *ncomp; comp_i++)
		{
			// create a new row
			sprintf(cname[0], "subtour_comp(%d)", comp_i); // constraint name

			// calculate rhs in O(nedges)
			rhs = 0.0;
			for (int i = 0; i < inst->nnodes; i++)
				if (comp_i == comp[i])
					rhs++;

			rhs--; // rhs = |S| -1

			lastrow = CPXgetnumrows(env,lp);
			if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) )
				print_error(" wrong CPXnewrows [degree]");
			// printf("added row %d with constraint %s %f\n", CPXgetnumrows(env,lp)-1, "<=",rhs);

			for (int i = 0; i < inst->nnodes; i++)
			{
				if (comp[i] == comp_i)
				{
					for (int j = i+1; j < inst->nnodes; j++)
					{
						if (comp[j] == comp_i)
						{
							if ( (*status = CPXchgcoef(env, lp, lastrow, xpos(i, j, inst), 1.0)) ) // income vertex
								print_error(" wrong CPXchgcoef [degree]");
						}
					}
				}
			}
		}
		subtour_counter++; //update subtour counter
		CPXmipopt(env,lp);
		CPXsolution(env, lp, status, &inst->best_lb, inst->best_sol, NULL, NULL, NULL);
		build_sol(inst, succ, comp, ncomp);
		if (inst->verbose >= 100) printf("Iter %3d partial solution, ncomp = %d\n", subtour_counter, *ncomp);
		if (inst->verbose >= 1000) plot_instance(inst);
	}
	if (inst->verbose >= 100) printf("best solution found. ncomp = %d\n", *ncomp);
	free(cname);
	free(succ);
	free(comp);
	free(ncomp);
	return 0;
}

int subtour_heur_iter_opt(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status, int heuristic) {

	// structure init
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int* ncomp = (int*)calloc(1, sizeof(int));
	char sense = 'L';
	double rhs;
	int lastrow;

	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));			// name of the variable
	*ncomp = inst->nnodes;
	int subtour_counter = 0;

	if (heuristic == 0) {
		CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);	// abort Cplex after the first incument update
		//CPXsetintparam(env, CPX_PARAM_NODELIM, 0);	// Solve only root node
		//CPXsetdblparam(env, CPX_PARAM_EPGAP, 0.01);  	// abort Cplex when gap below 10%

		CPXmipopt(env, lp);

		CPXsolution(env, lp, status, &inst->best_lb, inst->best_sol, NULL, NULL, NULL);
		build_sol(inst, succ, comp, ncomp);
		if (inst->verbose >= 100) printf("Iter %3d partial solution on Root node, ncomp = %d\n", subtour_counter, *ncomp);
		if (inst->verbose >= 1000) plot_instance(inst);
		while (*ncomp >= 2)
		{

			// subtour elimination constraints, one for each component
			for (int comp_i = 1; comp_i <= *ncomp; comp_i++)
			{
				// create a new row
				sprintf(cname[0], "subtour_comp(%d)", comp_i); // constraint name

				// calculate rhs in O(nedges)
				rhs = 0.0;
				for (int i = 0; i < inst->nnodes; i++)
					if (comp_i == comp[i])
						rhs++;

				rhs--; // rhs = |S| -1

				lastrow = CPXgetnumrows(env, lp);
				if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
					print_error(" wrong CPXnewrows [degree]");
				// printf("added row %d with constraint %s %f\n", CPXgetnumrows(env,lp)-1, "<=",rhs);

				for (int i = 0; i < inst->nnodes; i++)
				{
					if (comp[i] == comp_i)
					{
						for (int j = i + 1; j < inst->nnodes; j++)
						{
							if (comp[j] == comp_i)
							{
								if ( (*status = CPXchgcoef(env, lp, lastrow, xpos(i, j, inst), 1.0))) // income vertex
									print_error(" wrong CPXchgcoef [degree]");
							}
						}
					}
				}
			}
			subtour_counter++; //update subtour counter
			CPXmipopt(env, lp);
			CPXsolution(env, lp, status, &inst->best_lb, inst->best_sol, NULL, NULL, NULL);
			build_sol(inst, succ, comp, ncomp);
			if (inst->verbose >= 100) printf("Iter %3d partial solution, ncomp = %d\n", subtour_counter, *ncomp);
			if (inst->verbose >= 1000) plot_instance(inst);
		}
		if (inst->verbose >= 100) printf("best solution found in Root node. ncomp = %d\n", *ncomp);
		subtour_heur_iter_opt(env, lp, inst, status, 1);
	}
	else {
		CPXsetintparam(env, CPX_PARAM_INTSOLLIM, INT_MAX);

		CPXmipopt(env, lp);
		CPXsolution(env, lp, status, &inst->best_lb, inst->best_sol, NULL, NULL, NULL);
		build_sol(inst, succ, comp, ncomp);
		if (inst->verbose >= 100) printf("Iter %3d partial solution, ncomp = %d\n", subtour_counter, *ncomp);
		if (inst->verbose >= 1000) plot_instance(inst);
		while (*ncomp >= 2)
		{

			// subtour elimination constraints, one for each component
			for (int comp_i = 1; comp_i <= *ncomp; comp_i++)
			{
				// create a new row
				sprintf(cname[0], "subtour_comp(%d)", comp_i); // constraint name

				// calculate rhs in O(nedges)
				rhs = 0.0;
				for (int i = 0; i < inst->nnodes; i++)
					if (comp_i == comp[i])
						rhs++;

				rhs--; // rhs = |S| -1

				lastrow = CPXgetnumrows(env, lp);
				if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
					print_error(" wrong CPXnewrows [degree]");
				// printf("added row %d with constraint %s %f\n", CPXgetnumrows(env,lp)-1, "<=",rhs);

				for (int i = 0; i < inst->nnodes; i++)
				{
					if (comp[i] == comp_i)
					{
						for (int j = i + 1; j < inst->nnodes; j++)
						{
							if (comp[j] == comp_i)
							{
								if ( (*status = CPXchgcoef(env, lp, lastrow, xpos(i, j, inst), 1.0))) // income vertex
									print_error(" wrong CPXchgcoef [degree]");
							}
						}
					}
				}
			}
			subtour_counter++; //update subtour counter
			CPXmipopt(env, lp);
			CPXsolution(env, lp, status, &inst->best_lb, inst->best_sol, NULL, NULL, NULL);
			build_sol(inst, succ, comp, ncomp);
			if (inst->verbose >= 100) printf("Iter %3d partial solution, ncomp = %d\n", subtour_counter, *ncomp);
			if (inst->verbose >= 1000) plot_instance(inst);
		}
		if (inst->verbose >= 100) printf("best solution found. ncomp = %d\n", *ncomp);
	}

	free(cname);
	free(succ);
	free(comp);
	free(ncomp);
	return 0;
}


// callback
void switch_callback(tspinstance* inst, CPXENVptr env, CPXLPptr lp) {

	if (inst->callback == 1) {										// Lazy Constraint Callback
		CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);				// let MIP callbacks work on the original model
		if (CPXsetlazyconstraintcallbackfunc(env, lazycallback, inst)) {
			print_error(" Error in setLazyCallback()\n");
			CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_ON);			// reset to allow CPX working on its (probably reduce) model
			return;
		}
	} else if (inst->callback == 2) {									// Generic Callback
		CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
		if (CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE, genericcallback, inst)) {
			print_error(" Error in setGenericCallback2()\n");
			CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_ON);
			return;
		}
	} else if (inst->callback == 3) {
		CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
		if (CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_GLOBAL_PROGRESS, genericcallback, inst)) {
			print_error(" Error in setGenericCallback3()\n");
			CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_ON);
			return;
		}
	} else {																// Callback not used
		return;
	}

	CPXsetintparam(env, CPX_PARAM_THREADS, inst->nthread);			// it was reset after callback
	inst->ncols = CPXgetnumcols(env, lp);							// è necessario ricontrollare il numero di colonne?
}

int CPXPUBLIC lazycallback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p){			// Dynamic Search is disabled
	*useraction_p = CPX_CALLBACK_DEFAULT;					// solution is ok, don't add cuts
	tspinstance* inst = (tspinstance*)cbhandle; 			// casting of cbhandle (which is pointing to the above parameter inst)

	// get solution xstar
	double* xstar = (double*)malloc(inst->ncols * sizeof(double));
	if (CPXgetcallbacknodex(env, cbdata, wherefrom, xstar, 0, inst->ncols - 1))			// xstar = current x from CPLEX-- xstar starts from position 0 (getx is not defined)
		return 1;

	// get some random information at the node (as an example)
	double objval = CPX_INFBOUND; CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);					// Valore rilassamento continuo (lower bound) al nodo corrente
	int mythread = -1; CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &mythread);	// Thread chiamante
	double zbest; CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &zbest);			// Valore incumbent al nodo corrente

	//apply cut separator and possibly add violated cuts
	int ncuts = mylazy_separation(inst, xstar, env, cbdata, wherefrom);
	free(xstar);

	if (ncuts >= 1)
		*useraction_p = CPX_CALLBACK_SET; 					// tell CPLEX that cuts have been created
	return 0; 												// return 1 would mean error --> abort Cplex's execution
}

int CPXPUBLIC genericcallback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* cbhandle) {			// Dynamic Search is used anyway

	if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE) {
		tspinstance* inst = (tspinstance*)cbhandle; 			// casting of cbhandle (which is pointing to the above parameter inst)
				// get solution xstar
		double* xstar = (double*)malloc(inst->ncols * sizeof(double));
		double objval = CPX_INFBOUND;
		if (CPXcallbackgetcandidatepoint(context, xstar, 0, inst->ncols - 1, &objval))			// xstar = current x from CPLEX-- xstar starts from position 0 (getx is not defined)
			return 1;

		int mythread = -1; CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADS, &mythread);
		double zbest; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &zbest);				//valore incumbent al nodo corrente
		double best_int = -CPX_INFBOUND; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_BND, &best_int);
		inst->best_int = best_int;

		//apply cut separator and possibly add violated cuts
		mygeneric_separation(inst, xstar, context);
		free(xstar);											//avoid memory leak

		return 0; 												// return 1 would mean error --> abort Cplex's execution
	}
	if (contextid == CPX_CALLBACKCONTEXT_GLOBAL_PROGRESS) {
		tspinstance* inst = (tspinstance*)cbhandle;
		double best_int = CPX_INFBOUND; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &best_int);			// WARNING : it's not globally if callback isn't CPX_CALLBACKCONTEXT_GLOBAL_PROGRESS
		double best_lb = -CPX_INFBOUND; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_BND, &best_lb);

		//printf("\t\t____BEST SOL___ : %f, ____BEST BOUND___ : %f, ____Actual GAP___ : %f\n", best_sol, best_bound, best_sol - best_bound);
		//printf("Actual Gap : [%f]\n", best_sol - best_bound);

		if (inst->best_int > best_int) {
			if (inst->verbose >= 100)
				printf("BEST_INT update from -> to : [%f] -> [%f]\n", inst->best_int, best_int);
			inst->best_int = best_int;
		}

		if (inst->best_lb < best_lb) {
			if (inst->verbose >= 100)
				printf("BEST_LB update from -> to : [%f] -> [%f]\n", inst->best_lb, best_lb);
			inst->best_lb = best_lb;
		}

		return 0;
	}
	return 0;
}

int mylazy_separation(tspinstance* inst, const double* xstar, CPXCENVptr env, void* cbdata, int wherefrom) {
	// structure init
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int ncomp = INT_MAX;
	char sense = 'L';
	double rhs;
	int nnz;
	int val_length = 0;

	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));			// name of the variable

	build_sol_lazy_std(inst, xstar, succ, comp, &ncomp);
	if (ncomp >= 2){

		// subtour elimination constraints, one for each component
		for (int comp_i = 1; comp_i <= ncomp; comp_i++){

			// create a new row
			sprintf(cname[0], "lazy_subtour_comp(%d)", comp_i); // constraint name

			// calculate rhs in O(nedges)
			rhs = 0.0;
			for (int i = 0; i < inst->nnodes; i++)
				if (comp_i == comp[i])
					rhs++;

			rhs--; // rhs = |S| -1

			for (int i = 0; i < inst->nnodes; i++){
				if (comp[i] == comp_i){
					for (int j = i + 1; j < inst->nnodes; j++){
						if (comp[j] == comp_i){
							val_length++;
						}
					}
				}
			}
			nnz = 0;
			int* index = (int*)calloc(val_length, sizeof(int));
			double* value = (double*)calloc(val_length, sizeof(double));
			for (int i = 0; i < inst->nnodes; i++) {
				if (comp[i] == comp_i) {
					for (int j = i + 1; j < inst->nnodes; j++) {
						if (comp[j] == comp_i) {
							index[nnz] = xpos(i, j, inst);
							value[nnz] = 1.0;
							nnz++;
						}
					}
				}
			}
			if (CPXcutcallbackadd(env, cbdata, wherefrom, nnz, rhs, sense, index, value, 0))					// 0 means that the cut is managed by CPX
				print_error("User_Separation: CPXcutcallbackadd error");
				//if ( CPXcutcallbackaddlocal(env, cbdata, wherefrom, nnz, rhs, sense, index, value) )			// Usato per tagli validi localmente (ex. Gomory)
				//print_error("USER_separation: CPXcutcallbackaddlocal error");
			free(index);
			free(value);
		}
		if (inst->verbose >= 100) printf("LAZY partial solution, ncomp = %d\n", ncomp);
		free(cname);
		return ncomp == 1 ? 0 : ncomp;
	}
	if (inst->verbose >= 100) printf("best solution found in LAZY. ncomp = %d\n", ncomp);
	free(cname);
	return 0;
}

int mygeneric_separation(tspinstance* inst, const double* xstar, CPXCALLBACKCONTEXTptr context) {
	// structure init
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int ncomp = INT_MAX;
	char sense = 'L';
	double rhs;
	int nnz;
	int val_length = 0;
	int izero = 0;

	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));			// name of the variable

	build_sol_lazy_std(inst, xstar, succ, comp, &ncomp);
	if (ncomp >= 2) {

		// subtour elimination constraints, one for each component
		for (int comp_i = 1; comp_i <= ncomp; comp_i++) {

			// create a new row
			sprintf(cname[0], "GenericLazy_subtour_comp(%d)", comp_i); // constraint name

			// calculate rhs in O(nedges)
			rhs = 0.0;
			for (int i = 0; i < inst->nnodes; i++)
				if (comp_i == comp[i])
					rhs++;

			rhs--; // rhs = |S| -1

			for (int i = 0; i < inst->nnodes; i++) {
				if (comp[i] == comp_i) {
					for (int j = i + 1; j < inst->nnodes; j++) {
						if (comp[j] == comp_i) {
							val_length++;
						}
					}
				}
			}
			nnz = 0;
			int* index = (int*)calloc(val_length, sizeof(int));
			double* value = (double*)calloc(val_length, sizeof(double));
			for (int i = 0; i < inst->nnodes; i++) {
				if (comp[i] == comp_i) {
					for (int j = i + 1; j < inst->nnodes; j++) {
						if (comp[j] == comp_i) {
							index[nnz] = xpos(i, j, inst);
							value[nnz] = 1.0;
							nnz++;
						}
					}
				}
			}
			if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, value))					// 0 means that the cut is managed by CPX
				print_error("User_Separation: CPXcallbackrejectcandidate error");
			free(index);
			free(value);
		}
		if (inst->verbose >= 100) printf("GenericLAZY partial solution, ncomp = %d\n", ncomp);
		free(cname);
		return ncomp == 1 ? 0 : ncomp;
	}
	if (inst->verbose >= 100) printf("best solution found in GenericLAZY. ncomp = %d\n", ncomp);
	free(cname);
	return 0;
}


// Refining
void two_opt(tspinstance *inst) {
	// improving the best_sol if possibile in 2opt set

	if (inst->verbose >= 100) printf("two_opt\n");
	fflush(stdout);

	// check if current solution has only one tour
	int *succ = (int*) calloc(inst->nnodes, sizeof(int));
	int *comp = (int*) calloc(inst->nnodes, sizeof(int));
	int *ncomp = (int*) calloc(1, sizeof(int));
	build_sol(inst, succ, comp, ncomp);
	if (inst->verbose >= 100) print_succ(succ, inst);
	if (*ncomp != 1) print_error("call two_opt with best_sol with multiple tour");

	// search in 2opt if a better solution is found
	int i = 0;					// start from node 0
	int j = succ[succ[0]];
	double d_i1_i2;			// distance from i and succ[i], -
	double d_j1_j2;			// distance from j and succ[j], -
	double d_i1_j1;			// distance from i and j, +
	double d_i2_j2;			// dist from succ[i] and succ[j], +
	double cur_improve;	// + d_i1_i2 + d_j1_j2 - d_i1_j1 - d_i2_j2
	int best_i = i;			// i that perform the best improve
	int best_j = i;			// j that perform the best improve
	double best_improve = 0.0;

	for (int ti = 0; ti < inst->nnodes - 3; ti++) {
		d_i1_i2 = dist(i, succ[i], inst);
		for (int tj = ti+2; tj < inst->nnodes; tj++) {  // no 2opt with consequence arches
			if (i == succ[j]) break;  // should happen only when i = 0,
			d_j1_j2 = dist(j, succ[j], inst);
			d_i1_j1 = dist(i, j, inst);
			d_i2_j2 = dist(succ[i], succ[j], inst);
			if (inst->verbose > 100) printf("*** i - j: %d - %d\n", i,j);
			cur_improve = (d_i1_i2 + d_j1_j2) - (d_i1_j1 + d_i2_j2);
			if ( cur_improve > best_improve ) {  // cross is better
				best_i = i;
				best_j = j;
				best_improve = cur_improve;
				if (inst->verbose >= 100) printf("BTO best_improve: %7.1f\n", best_improve);
			}
			j = succ[j];
		}
		i = succ[i];
		j = succ[succ[i]];
	}

	// if find new best solution change the arches
	if (best_i != best_j) {

		int i2 = succ[best_i];
		int j2 = succ[best_j];
		int pre_node = succ[i2];
		inst->best_sol[xpos(best_i, i2, inst)] = 0.0;
		inst->best_sol[xpos(best_j, j2, inst)] = 0.0;

		succ[best_i] = best_j;
		succ[i2] = j2;
		inst->best_sol[xpos(best_i, best_j, inst)] = 1.0;
		inst->best_sol[xpos(i2, j2, inst)] = 1.0;

		int cur_node = i2;
		int suc_node = j2;
		while (cur_node != best_j) {  // reverse the succ
			suc_node = cur_node;
			cur_node = pre_node;
			pre_node = succ[pre_node];
			succ[cur_node] = suc_node;
		}
		if (inst->verbose >= 100) print_succ(succ, inst);
	}
	if (inst->verbose >= 100) printf("BTO i j: %3d %3d\n", i,j);
	if (inst->verbose >= 100) printf("BTO Best_improve: %f\n", best_improve);
	inst->best_lb -= best_improve;  // update best_lb
	if (inst->verbose >= 100) printf("BTO new best_lb: %7.1f\n", inst->best_lb);
	if (inst->verbose >= 100) print_succ(succ, inst);
	free(succ);
	free(comp);
	free(ncomp);
}

void random_two_opt(tspinstance* inst) {

	// return random solution in inst->best_sol
	// the improvement of best_lb is already applied to inst->best_lb

	if (inst->verbose > 110) printf("Random two_opt\n");

	// check if current solution has only one tour
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int* ncomp = (int*)calloc(1, sizeof(int));
	build_sol(inst, succ, comp, ncomp);
	if (inst->verbose > 110) print_succ(succ, inst);
	if (*ncomp != 1) print_error("call random_two_opt with best_sol with multiple tour");

	// random solution in 2opt

	int i = rand() % inst->nnodes;
	int j = rand() % inst->nnodes;
	while (j == i || j == succ[i]) {
		j = rand() % inst->nnodes;
	}

	double d_i1_i2 = dist(i, succ[i], inst);
	double d_j1_j2 = dist(j, succ[j], inst);
	double d_i1_j1 = dist(i, j, inst);
	double d_i2_j2 = dist(succ[i], succ[j], inst);
	double cur_improve = (d_i1_i2 + d_j1_j2) - (d_i1_j1 + d_i2_j2);

	// Change the arches
	if (i != j) {

		int i2 = succ[i];
		int j2 = succ[j];
		int pre_node = succ[i2];
		succ[i] = j;
		succ[i2] = j2;
		int cur_node = i2;
		int suc_node = j2;
		while (cur_node != j) {  // reverse the succ
			suc_node = cur_node;
			cur_node = pre_node;
			pre_node = succ[pre_node];
			succ[cur_node] = suc_node;
			if (inst->verbose > 110) print_succ(succ, inst);
		}
	}
	inst->best_lb -= cur_improve;  // update best_lb

	// Assign random solution found in inst->best_sol
	clear_sol(inst);
	int first_node = 0;
	int second_node = succ[first_node];
	for (int i = 0; i < inst->nnodes; i++) {
		inst->best_sol[xpos(first_node, second_node, inst)] = 1.0;
		first_node = second_node;
		second_node = succ[second_node];
	}

	if (inst->verbose > 110) print_succ(succ, inst);
	free(succ);
	free(comp);
	free(ncomp);
}

void random_n_opt(tspinstance* inst, int n) {

	if (inst->verbose >= 110) printf("RANDOM_N_OPT\n");
	if (inst->nnodes < n) print_error("n should be lower than the number of nodes of the problem.");
	fflush(stdout);
	// check if current solution has only one tour
	int *succ = (int*) calloc(inst->nnodes, sizeof(int));
	int *comp = (int*) calloc(inst->nnodes, sizeof(int));
	int *ncomp = (int*) calloc(1, sizeof(int));
	build_sol(inst, succ, comp, ncomp);
	if (inst->verbose >= 110) print_succ(succ, inst);
	if (*ncomp != 1) print_error("call random_n_opt with best_sol with multiple tour");

	// ramdoly select n nodes and save n successor
	int *c_1 = (int*) calloc(n, sizeof(int));		// this will be used as array
	int *c_2 = (int*) calloc(n, sizeof(int));		// of pairs

	// 	initialization to -1
	for (int i = 0; i < n; i++) {c_1[i] = -1; c_2[i] = -1;}

		// set random seed
	int i;  // tour iterator, assume index nodes value
	int c_iter = 0;		// iterate over n pairs, each time a random index is selected
	int already_selected;  // tell if i has already been added in the c_2 struct
	while (c_iter < n) {
		i = rand() % inst->nnodes; 	// random node
		already_selected = 0;

		// check if i have been already selected
		for (int j = 0; j < c_iter; j++) if (c_2[j] == i) already_selected = 1;

		if (!already_selected) {  // add to selected
			c_2[c_iter] = i;

			if (succ[i] == -1)  // isolated node, should never happen
				print_error("Selected isolated node. Random_n_opt error");
			else
				c_1[c_iter] = succ[i];

			c_iter++;		// update c_iter
			inst->best_lb -= dist(i, succ[i], inst);  // update distances
			inst->best_sol[xpos(i, succ[i], inst)] = 0.0;
			succ[i] = -1;  // break the edges
		} else {
			// continue and choose another index
		}
	}

	// for n times, merge a random node with random successor. (no need to reverse the orientation)

	int i_1 = rand() % n;
	int first_node = c_1[i_1];
	int i_2 = -1;
	c_iter = 0;
	while (c_iter < n) {

		// select one of the c_1 at random and set as first nodes
		int i = c_1[i_1];
		while (succ[i] != -1) i = succ[i];

		// select a c_1[i_1] not already selected and != from first_node
		int already_selected = 0;
		int is_first = 0;
		for (i_1 = 0; i_1 < n; i_1++){

			is_first = c_1[i_1] == first_node;  // the node is the first_node, skip
			if (is_first) continue;

			// check if c_1[i_1] is already selected
			for (int j = 0; j < n; j++) {
				already_selected = succ[c_2[j]] == c_1[i_1];
				if (already_selected) {  // already selected c_1[i]
					break;		// no matter what, it can not be selected, change c_1[i]
				}
			}
			// if not already seelected here, it is selected
			if (!already_selected) break;
		}

		// check if it is the last edge or not
		if (i_1 == n) { // close the tour and exit
			succ[i] = first_node;
			inst->best_lb += dist(i, first_node, inst);
			inst->best_sol[xpos(i, first_node, inst)] = 1.0;
			c_iter++;
			break;
		} else {  	// add new arch i -> c_1[i_1]
			succ[i] = c_1[i_1];
			inst->best_lb += dist(i, c_1[i_1], inst);
			inst->best_sol[xpos(i, c_1[i_1], inst)] = 1.0;
			c_iter++;
		}
	}

	if (inst->verbose >= 110) print_succ(succ, inst);
	fflush(stdout);
	free(succ);
	free(comp);
	free(ncomp);
	free(c_1);
	free(c_2);
}


// Repair
void single_patch(tspinstance* inst, int* succ, int* comp, int* ncomp) {

	if (inst->verbose >= 100) printf("\nSINGLE_PATCHING"); fflush(stdout);
	// single patch merge the isolated nodes or tours,
	// if already merged, exit the method
	if (*ncomp == 1) return;

	int initial_i = rand() % inst->nnodes;  // randomize the first node

	// get the closer node which is in another tour
	int closer_j = initial_i;
	double cn_dist = INT_MAX;		// dist(initial_i, closer_j, inst)
	double d_ij = 0.;
	for (int j = 0; j < inst->nnodes; j++) {
		if (comp[initial_i] == comp[j]) continue; // skip if it's in the same components

		d_ij = dist(initial_i, j, inst);
		if (d_ij < cn_dist) {
			closer_j = j;
			cn_dist = d_ij;
		}
	}

	int i = initial_i;
	int i_tour_size = 1;
	while (succ[i] != -1 && succ[i] != initial_i) { i = succ[i]; i_tour_size++; } // components of i_tour

	int j = closer_j;
	int j_tour_size = 1;
	while (succ[j] != -1 && succ[j] != closer_j) { j = succ[j]; j_tour_size++; } // components of j_tour

	// find the node of the i-tour which is closer to closer_j: closer_i
	i = initial_i;
	int closer_i = i;
	int counter = 1;
	if (succ[i] > 0) { // only if i is not isolated

		double best_improve;
		if (succ[closer_j] > 0) { // j is not isolated
			best_improve = dist(i, succ[closer_j], inst) + dist(succ[i], closer_j, inst)
													- dist(i, succ[i], inst) - dist(closer_j, succ[closer_j], inst);
		} else {									// j is isolated
			best_improve = dist(i, closer_j, inst) + dist(closer_j, succ[i], inst)
													- dist(i, succ[i], inst);
		}
		double cur_improve = best_improve;
		i = succ[i];

		while (i != initial_i) {  // for each node of the i-tour
			if (succ[closer_j] > 0) {  	// j is not isolated
				cur_improve = dist(i, succ[closer_j], inst) + dist(succ[i], closer_j, inst)
													- dist(i, succ[i], inst) - dist(closer_j, succ[closer_j], inst);
			} else  {  									// j is isolated
				cur_improve = dist(i, closer_j, inst) + dist(closer_j, succ[i], inst)
														- dist(i, succ[i], inst);
			}

			if (cur_improve < best_improve) {
				best_improve = cur_improve;
				closer_i = i;
			}
			counter++;
			i = succ[i];
		}
	}
	cn_dist = dist(closer_i, closer_j, inst);

	// merge the tours: update best_sol, succ, comp, ncomp, best_lb
	if (succ[closer_j] < 0 && succ[closer_i] < 0) { 	// closer_j and closer_i are isolated
		// update inst->best_sol
		inst->best_sol[xpos(closer_i, closer_j, inst)] = 1.;

		// update succ
		succ[closer_j] = closer_i;
		succ[closer_i] = closer_j;

		// update comp
		comp[closer_j] = comp[closer_i];
	} else if (succ[closer_j] < 0) { 									// closer_j isolated node
		// update inst->best_sol
		if (closer_i != succ[succ[closer_i]])  // i_tour_size == 2
			inst->best_sol[xpos(closer_i, succ[closer_i], inst)] = 0.;
		inst->best_sol[xpos(closer_i, closer_j, inst)] = 1.;
		inst->best_sol[xpos(closer_j, succ[closer_i], inst)] = 1.;

		// update succ
		if (closer_i == succ[succ[closer_i]]) {
			if ( !is_clockwise(inst, closer_i, closer_j, succ[closer_i]) ) {
				closer_i = succ[closer_i];	// reverse the orientation
			}
		}
		succ[closer_j] = succ[closer_i];
		succ[closer_i] = closer_j;

		// update comp
		comp[closer_j] = comp[closer_i];

	} else if (succ[closer_i] < 0) {  								// i is isolated
		// update inst->best_sol
		if (closer_j != succ[succ[closer_j]])
			inst->best_sol[xpos(closer_j, succ[closer_j], inst)] = 0.;
		inst->best_sol[xpos(closer_i, closer_j, inst)] = 1.;
		inst->best_sol[xpos(closer_i, succ[closer_j], inst)] = 1.;

		// update succ
		if (closer_j == succ[succ[closer_j]]) {
			if ( !is_clockwise(inst, closer_j, closer_i, succ[closer_j]) ) {
				closer_j = succ[closer_j];	// reverse the orientation
			}
		}
		succ[closer_i] = succ[closer_j];
		succ[closer_j] = closer_i;

		// update comp
		comp[closer_i] = comp[closer_j];

	} else {  																				// i and closer_j are not isolated
		// update inst->best_sol
		if (closer_j != succ[succ[closer_j]])
			inst->best_sol[xpos(closer_j, succ[closer_j], inst)] = 0.;
		if (closer_i != succ[succ[closer_i]])
			inst->best_sol[xpos(closer_i, succ[closer_i], inst)] = 0.;
		inst->best_sol[xpos(closer_i, succ[closer_j], inst)] = 1.;
		inst->best_sol[xpos(closer_j, succ[closer_i], inst)] = 1.;

		// update succ
		int tmp = succ[closer_i];
		succ[closer_i] = succ[closer_j];
		succ[closer_j] = tmp;

		// update comp
		comp[closer_i] = comp[closer_j];
		while (tmp != closer_i) {
			comp[tmp] = comp[closer_j];
			tmp = succ[tmp];
		}
	}

	// update inst->best_lb
	inst->best_lb = 0.0;
	for (int i = 0; i < inst->nnodes; i++)
	 	for (int j = i+1; j < inst->nnodes; j++)
		 	if (inst->best_sol[xpos(i, j, inst)] > 0.5)
				inst->best_lb += dist(i, j, inst);
	(*ncomp)--;
}

int is_clockwise(tspinstance *inst, int x1, int x2, int x3) {
	return (inst->xcoord[x3] - inst->xcoord[x1])*(inst->ycoord[x2] - inst->ycoord[x1]) <
				 (inst->ycoord[x3] - inst->ycoord[x1])*(inst->xcoord[x2] - inst->xcoord[x1]);
}

void clear_sol(tspinstance* inst) {
	for (int i = 0; i < inst->nedges; i++)
		inst->best_sol[i] = 0.0;
}

// build_sol methods use the optimized solution to plot
void build_sol(tspinstance *inst, int *succ, int *comp, int *ncomp) {

	switch (inst->model_type)
	{
		case -1:
		  break;

		case 0 :		// basic model with asymmetric x and q
			build_sol_sym(inst, succ, comp, ncomp);
			break;

		case 1 :		// MTZ contraints
		 	build_sol_mtz(inst, succ, comp, ncomp);
			break;

		case 2:			// FLOW 1
			build_sol_flow1(inst, succ, comp, ncomp);
			break;

		case 3:			// MTZ contraints
			build_sol_mtz(inst, succ, comp, ncomp);
			break;

		default:
			print_error(" model type unknown!!");
			break;
	}
}

void build_sol_sym(tspinstance *inst, int *succ, int *comp, int *ncomp) {	// build succ() and comp() wrt xstar()...

	// check if nodes degree is 2 or 0 (isolated node) for each node
	if (inst->verbose >= 2000) {
		int *degree = (int *) calloc(inst->nnodes, sizeof(int));
		printf("nnodes=%d\n", inst->nnodes);
		for (int i = 0; i < inst->nnodes; i++) {
			for (int j = i+1; j < inst->nnodes; j++) {
				int k = xpos(i,j,inst);
				if (fabs(inst->best_sol[k]) > EPS && fabs(inst->best_sol[k]-1.0) > EPS ) print_error(" wrong inst->best_sol in build_sol()");
				if (inst->best_sol[k] > 0.5) {
					printf("x[%d,%d] = 1\n", i, j);
					++degree[i];
					++degree[j];
				} else {
					printf("x[%d,%d] = 0\n", i, j);
				}

			}
		}
		for (int i = 0; i < inst->nnodes; i++) {
			if (degree[i] != 2 || degree[i] != 0) {
				char msg[40];
				snprintf(msg, sizeof(msg), "wrong degree[%d] = %d in build_sol_sym", i, degree[i]);
				print_error(msg);
			}
		}
		free(degree);
	}

	// initialization of succ, comp and ncomp
	*ncomp = 0;
	for (int i = 0; i < inst->nnodes; i++) {
		succ[i] = -1;
		comp[i] = -1;
	}

	// tour
	for (int start = 0; start < inst->nnodes; start++) {
		if (comp[start] >= 0)
			continue;

		(*ncomp)++;						// the tour id number
		int prv = -1;					// previous node of i in the tour, used to keep j != i
		int i = start;				// iterate over nodes to complete the tour
		int found_succ = 0;		// 1 when found a succ of i
		while (comp[start] == -1) {
			found_succ = 0;
			for (int j = 0; j < inst->nnodes; j++) {	// j iterate to be the subsequent node of i in the tour
				if (i != j && inst->best_sol[xpos(i, j, inst)] > 0.5  && j != prv) {  // found subsequent of i (j)
					succ[i] = j;
					comp[j] = *ncomp;
					prv = i;
					i = j;
					found_succ = 1;
					break;
				}
			}
			if (!found_succ) {  // no succ found, i is isolated
				if (prv == -1) {	// if no prv found either
					comp[start] = *ncomp;
					break;
				} else {  // it's not a tour, i has no subsequent, but is connected to prv
					print_error("inst->best_sol is not made of tours and isolated nodes");
				}
			}
		}
	}

	// print succ, comp and ncomp
	if (inst->verbose >= 2000) {
		printf("\ni:    "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", i);
		printf("\nsucc: "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", succ[i]);
		printf("\ncomp: "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", comp[i]);
		printf("\n");
		fflush(stdout);
	}
}

void build_sol_lazy_std(tspinstance* inst, const double* xstar, int* succ, int* comp, int* ncomp) {	// build succ() and comp() wrt xstar()...

	// check if nodes degree is 2 for each node
	if (inst->verbose >= 2000)
	{
		int* degree = (int*)calloc(inst->nnodes, sizeof(int));
		printf("nnodes=%d\n", inst->nnodes);
		for (int i = 0; i < inst->nnodes; i++)
		{
			for (int j = i + 1; j < inst->nnodes; j++)
			{
				int k = xpos(i, j, inst);
				if (fabs(xstar[k]) > EPS && fabs(xstar[k] - 1.0) > EPS) print_error(" wrong xstar in build_sol()");
				if (xstar[k] > 0.5)
				{
					printf("x[%d,%d] = 1\n", i, j);
					++degree[i];
					++degree[j];
				}
				else {
					printf("x[%d,%d] = 0\n", i, j);
				}

			}
		}
		for (int i = 0; i < inst->nnodes; i++)
		{
			if (degree[i] != 2)
			{
				char msg[40];
				snprintf(msg, sizeof msg, "wrong degree[%d] = %d in build_sol_sym", i, degree[i]);
				print_error(msg);
			}
		}
		free(degree);
	}

	// initialization of succ, comp and ncomp
	*ncomp = 0;
	for (int i = 0; i < inst->nnodes; i++)
	{
		succ[i] = -1;
		comp[i] = -1;
	}

	// tour
	for (int start = 0; start < inst->nnodes; start++) {
		if (comp[start] >= 0)
			continue;

		(*ncomp)++;
		int prv = -1;
		int i = start;
		int found_succ = 0;
		while (comp[start] == -1) {
			for (int j = 0; j < inst->nnodes; j++) {
				if (i != j && xstar[xpos(i, j, inst)] > 0.5 && j != prv) {

					succ[i] = j;
					comp[j] = *ncomp;
					prv = i;
					i = j;
					found_succ = 1;
					break;
				}
			}
			if (!found_succ) {  // no succ found, i is isolated
				comp[start] = *ncomp;
				break;
			}
		}
	}
	/* SBAGLIATO?
	for (int start = 0; start < inst->nnodes; start++)
	{
		if (comp[start] >= 0) continue;  // node "start" has already been setted

		// a new component is found
		(*ncomp)++;
		comp[start] = *ncomp;

		int i = start;
		//int done = 0;
		while (succ[i] == -1)  // go and visit the current component
		{
			comp[i] = *ncomp;
			// done = 1;
			for (int j = 0; j < inst->nnodes; j++)
			{
				if (j == i) continue;

				if (xstar[xpos(i, j, inst)] > 0.5) // the edge [i,j] is selected in xstar and j was not visited before
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
	*/

	// print succ, comp and ncomp
	if (inst->verbose >= 2000)
	{
		printf("\ni:      "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", i);
		printf("\nsucc:   "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", succ[i]);
		printf("\ncomp:   "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", comp[i]);
		printf("\n");
	}
}

void build_sol_mtz(tspinstance *inst, int *succ, int *comp, int *ncomp) {	// build succ() and comp() wrt xstar()...

	// check if nodes degree is 2 for each node
	if (inst->verbose >= 2000) {
		printf("debugging sym mtz solution...\n");
		if (inst->verbose >= 2000) {
			printf("Solution:\n      ");
			for (int i = 0; i < inst->nnodes; i++) printf("%5d|", i);
			printf("\n");
			for (int i = 0; i < inst->nnodes; i++) {
				printf("%5d)", i);
				for (int j = 0; j < inst->nnodes; j++)
					if (i == j) printf("      ");
					else printf("%6.1f", round(inst->best_sol[asym_xpos(i,j,inst)]) );
				printf("\n");
			}
			printf("\n    u)      ");
			for (int u = 1; u < inst->nnodes; u++)
				printf("%6.1f", round(inst->best_sol[asym_upos(u,inst)]) );
		}
		int *degree = (int *) calloc(inst->nnodes, sizeof(int));
		for ( int i = 0; i < inst->nnodes; i++ ) {
			for ( int j = 0; j < inst->nnodes; j++ ) {
				if (i == j) continue;
				int k = asym_xpos(i,j,inst);
				if ( fabs(inst->best_sol[k]) > EPS && fabs(inst->best_sol[k]-1.0) > EPS ) print_error(" wrong inst->best_sol in build_sol()");
				if ( inst->best_sol[k] > 0.5 ) {
					++degree[i];
					++degree[j];
				}
			}
		}
		for ( int i = 0; i < inst->nnodes; i++ ) {

			if ( degree[i] != 2 ) {
				char msg[100];
				snprintf(msg, sizeof(msg), "wrong degree in build_sol_mtz: degree(%d) = %d", i, degree[i]);
				print_error(msg);
			}
		}
		free(degree);
		printf("\ndebug completed succesfully.\n");
	}

	// initialization of succ, comp and ncomp
	*ncomp = 0;
	for ( int i = 0; i < inst->nnodes; i++ ) { succ[i] = -1; comp[i] = -1; }

	// tour
	for ( int start = 0; start < inst->nnodes; start++ ) {
		if ( comp[start] >= 0 ) continue;  // node "start" has already been setted

		// a new component is found
		(*ncomp)++;
		comp[start] = *ncomp;

		int i = start;
		while ( succ[i] == -1 ) { // go and visit the current component
			comp[i] = *ncomp;
			for ( int j = 0; j < inst->nnodes; j++ ) {
				if (j == i) continue;

				if ( inst->best_sol[asym_xpos(i,j,inst)] > 0.5) { // the edge [i,j] is selected in inst->best_sol and j was not visited before
					// intern edge of the cycle
					if (comp[j] == -1) {
						succ[i] = j;
						i = j;
						break;
					}
					// last edge of the cycle
					if (start == j) succ[i] = j;
				}
			}
		}	// while
	// go to the next component...
	}

	// print succ, comp and ncomp
	if (inst->verbose >= 2000) {
		printf("\ni:   "); for (int i = 0; i < inst->nnodes; i++) printf("%5d", i);
		printf("\nsucc:"); for (int i = 0; i < inst->nnodes; i++) printf("%5d", succ[i]);
		printf("\ncomp:"); for (int i = 0; i < inst->nnodes; i++) printf("%5d", comp[i]);
		printf("\n");
	}
}

void build_sol_flow1(tspinstance *inst, int *succ, int *comp, int *ncomp) {	// build succ() and comp() wrt xstar()...

	// check if nodes degree is 2 for each node
	if (inst->verbose >= 2000) {
		printf("debugging sym flow1 solution...\n");
		if (inst->verbose >= 2000) {
			printf("Solution:\n      ");
			for (int i = 0; i < inst->nnodes; i++) printf("%5d|", i);
			printf("\n");
			for (int i = 0; i < inst->nnodes; i++) {
				printf("%5d)", i);
				for (int j = 0; j < inst->nnodes; j++)
					if (i == j) printf("      ");
					else printf("%6.1f", round(inst->best_sol[asym_xpos(i,j,inst)]) );
				printf("\n");
			}
			// printf("\n    u)      ");
			// for (int u = 1; u < inst->nnodes; u++)
			// 	printf("%6.1f", round(inst->best_sol[asym_upos(u,inst)]) );
		}
		int *degree = (int *) calloc(inst->nnodes, sizeof(int));
		for ( int i = 0; i < inst->nnodes; i++ ) {
			for ( int j = 0; j < inst->nnodes; j++ ) {
				if (i == j) continue;
				int k = asym_xpos(i,j,inst);
				if ( fabs(inst->best_sol[k]) > EPS && fabs(inst->best_sol[k]-1.0) > EPS ) print_error(" wrong inst->best_sol in build_sol()");
				if ( inst->best_sol[k] > 0.5 ) {
					++degree[i];
					++degree[j];
				}
			}
		}
		for ( int i = 0; i < inst->nnodes; i++ ) {

			if ( degree[i] != 2 ) {
				char msg[100];
				snprintf(msg, sizeof(msg), "wrong degree in build_sol_mtz: degree(%d) = %d", i, degree[i]);
				print_error(msg);
			}
		}
		free(degree);
		printf("\ndebug completed succesfully.\n");
	}

	// initialization of succ, comp and ncomp
	*ncomp = 0;
	for ( int i = 0; i < inst->nnodes; i++ ) { succ[i] = -1; comp[i] = -1; }

	// tour
	for ( int start = 0; start < inst->nnodes; start++ ) {
		if ( comp[start] >= 0 ) continue;  // node has already been setted

		// a new component is found
		(*ncomp)++;
		comp[start] = *ncomp;

		int i = start;
		while ( succ[i] == -1 ) { // go and visit the current component
			comp[i] = *ncomp;
			for ( int j = 0; j < inst->nnodes; j++ ) {
				if (j == i) continue;

				if ( inst->best_sol[asym_xpos(i,j,inst)] > 0.5) { // the edge [i,j] is selected in inst->best_sol and j was not visited before
					// intern edge of the cycle
					if (comp[j] == -1) {
						succ[i] = j;
						i = j;
						break;
					}
					// last edge of the cycle
					if (start == j) succ[i] = j;
				}
			}
		}	// while
	// go to the next component...
	}

	// print succ, comp and ncomp
	if (inst->verbose >= 2000) {
		printf("\ni:   "); for (int i = 0; i < inst->nnodes; i++) printf("%5d", i);
		printf("\nsucc:"); for (int i = 0; i < inst->nnodes; i++) printf("%5d", succ[i]);
		printf("\ncomp:"); for (int i = 0; i < inst->nnodes; i++) printf("%5d", comp[i]);
		printf("\n");
	}
}


// distance functions
double dist(int i, int j, tspinstance *inst) {
	switch (inst->edge_weight_type) {
		case 0: {
			double dx = inst->xcoord[i] - inst->xcoord[j];
			double dy = inst->ycoord[i] - inst->ycoord[j];
			double rij = sqrt((dx * dx + dy * dy) / 10.0);
			double tij = round(rij);
			if (tij < rij)
				return tij + 1.0;
			else
				return tij;
			break;
		}
		case 1: {
			double dx = inst->xcoord[i] - inst->xcoord[j];
			double dy = inst->ycoord[i] - inst->ycoord[j];
			if (!inst->integer_costs) return sqrt(dx * dx + dy * dy);
			return round(sqrt(dx * dx + dy * dy)); 			// nearest integer
			break;
		}
		case 2: {
			double PI = 3.141592;
			double RRR = 6378.388;

			double deg = (int)inst->xcoord[i] + 0.0;
			double min = inst->xcoord[i] - deg;								// min è diviso per 100
			double latitude_i = PI * (deg + 5.0 * min / 3.0) / 180.0;		// PI * ( deg + min/60)/180
			deg = (int)inst->ycoord[i] + 0.0;
			min = inst->ycoord[i] - deg;
			double longitude_i = PI * (deg + 5.0 * min / 3.0) / 180.0;

			deg = (int)inst->xcoord[j] + 0.0;
			min = inst->xcoord[j] - deg;
			double latitude_j = PI * (deg + 5.0 * min / 3.0) / 180.0;
			deg = (int)inst->ycoord[j] + 0.0;
			min = inst->ycoord[j] - deg;
			double longitude_j = PI * (deg + 5.0 * min / 3.0) / 180.0;

			double q1 = cos(longitude_i - longitude_j);
			double q2 = cos(latitude_i - latitude_j);
			double q3 = cos(latitude_i + latitude_j);
			return (int)(RRR * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0) + 0.0;
			break;
		}
		default:
			print_error(" model type unknown!!");
			return -1;
		break;
	}

}


// User input
void read_input(tspinstance *inst) { // simplified CVRP parser, not all SECTIONs detected

	FILE *fin = fopen(inst->input_file, "r");
	if ( fin == NULL ) print_error(" input file not found!");

	inst->nnodes = -1;
	inst->nedges = -1;
	inst->best_lb = INT_MAX;

	char line[180];
	char *par_name;
	char *token1;
	char *token2;

	int active_section = 0; // =1 NODE_COORD_SECTION, =2 DEMAND_SECTION, =3 DEPOT_SECTION

	while ( fgets(line, sizeof(line), fin) != NULL )
	{
		if ( inst->verbose >= 2000 ) { printf("read_input: line: %s",line); fflush(NULL); }
		if ( strlen(line) <= 1 ) continue; // skip empty lines

	  par_name = strtok(line, " :");
		if ( inst->verbose >= 2000 ) { printf("read_input: parameter: \"%s\" ",par_name); fflush(NULL); }

		if ( strncmp(par_name, "NAME", 4) == 0 )
		{
			active_section = 0;
			continue;
		}

		if ( strncmp(par_name, "COMMENT", 7) == 0 )
		{
			active_section = 0;
			token1 = strtok(NULL, "");
			if ( inst->verbose >= 2000 ) printf("read_input: solving instance: %s\b with model %d...\n", token1, inst->setup_model);
			continue;
		}

		if ( strncmp(par_name, "TYPE", 4) == 0 )
		{
			token1 = strtok(NULL, " :");
			if ( strncmp(token1, "TSP",3) != 0 ) print_error(" format error:  only TYPE == TSP implemented so far!!!!!!");
			active_section = 0;
			continue;
		}

		if ( strncmp(par_name, "DIMENSION", 9) == 0 )
		{
			if ( inst->nnodes >= 0 ) print_error(" repeated DIMENSION section in input file");
			token1 = strtok(NULL, " :");
			inst->nnodes = atoi(token1);
			if ( inst->verbose >= 2000 ) printf("\tnnodes: %d\n", inst->nnodes);
			inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double)*inst->nnodes);
			inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double)*inst->nnodes);
			active_section = 0;
			continue;
		}

		if ( strncmp(par_name, "DISPLAY_DATA_TYPE", 17) == 0 )
		{
			token1 = strtok(NULL, " :");
			if ( strncmp(token1, "COORD_DISPLAY", 13) != 0 ) print_error(" format error:  only DISPLAY_DATA_TYPE == COORD_DISPLAY implemented so far!!!!!!");
			active_section = 0;
			continue;
		}

		if ( strncmp(par_name, "EDGE_WEIGHT_FORMAT", 18) == 0 )
		{
			token1 = strtok(NULL, " :");
			if ( strncmp(token1, "FUNCTION", 8) != 0 ) print_error(" format error:  only EDGE_WEIGHT_FORMAT == FUNCTION implemented so far!!!!!!");
			active_section = 0;
			continue;
		}



		if ( strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0 )
		{
			token1 = strtok(NULL, " :");
			if (	strncmp(token1, "ATT", 3) != 0 &&
						strncmp(token1, "EUC_2D", 6) != 0 &&
						strncmp(token1, "GEO", 3) != 0 )
				print_error(" format error:  only EDGE_WEIGHT_TYPE == EUC_2D, ATT, GEO implemented so far!!!!!!");
			active_section = 0;
			inst->edge_weight_type = strncmp(token1, "ATT", 3) == 0 ? 0 : strncmp(token1, "EUC_2D", 6) == 0 ? 1 : 2;
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
			if ( inst->verbose >= 2000 )
				printf("\tnode %4d = (%.2f, %.2f)\n", i+1, inst->xcoord[i], inst->ycoord[i]);
			continue;
		}

		printf("read_input: final active section %d\n", active_section);
		print_error("... wrong format for the current simplified parser!!!!");
	}

	fclose(fin);
}

void parse_command_line(int argc, char** argv, tspinstance *inst) {

	// default
	inst->setup_model = 0;
	inst->randomseed = 0;
	inst->nthread = 1;
	inst->timelimit = 300.; 	// CPX_INFBOUND;
	strcpy(inst->input_file, "NULL");

	inst->available_memory = 12000;   // available memory, in MB, for Cplex execution (e.g., 12000)
	inst->max_nodes = -1; 						// max n. of branching nodes in the final run (-1 unlimited)
	inst->integer_costs = 1;
	inst->verbose = 1000;							// VERBOSE

	// hard fixing
	double max_fr = 0.9;		// maximum fixing_ratio
	double incr_fr = 0.1;		// increase of fixing_ratio when good gap
	double decr_fr = 0.1;		// decreasing of fixing_ratio when bad gap
	double good_gap = 0.1;			// a good gap allow to decrease fixing_ratio
	double optimal_gap = 0.05; 	// under this value, solution is optimal

  int help = 0; if ( argc < 1 ) help = 1;
	for ( int i = 1; i < argc; i++ )
  {
		if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-input") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-f") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 				// input file
		if ( strcmp(argv[i],"-time_limit") == 0 ) { inst->timelimit = atof(argv[++i]); continue; }		// total time limit

		// hard fixing
		if ( strcmp(argv[i],"-max_fr") == 0 ) { max_fr = atof(argv[++i]); continue; }
		if ( strcmp(argv[i],"-incr_fr") == 0 ) { incr_fr = atof(argv[++i]); continue; }
		if ( strcmp(argv[i],"-decr_fr") == 0 ) { decr_fr = atof(argv[++i]); continue; }
		if ( strcmp(argv[i],"-good_gap") == 0 ) { good_gap = atof(argv[++i]); continue; }
		if ( strcmp(argv[i],"-optimal_gap") == 0 ) { optimal_gap = atof(argv[++i]); continue; }

		if ( strcmp(argv[i],"-setup_model") == 0) { inst->setup_model = atoi(argv[++i]); continue; } 	// model type
		if ( strcmp(argv[i],"-model") == 0 ) { inst->setup_model = atoi(argv[++i]); continue; } 			// model type
		if ( strcmp(argv[i],"-randomseed") == 0 ) { inst->randomseed = abs(atoi(argv[++i])); continue; } 		// random seed
		if ( strcmp(argv[i],"-nthread") == 0 ) { inst->nthread = abs(atoi(argv[++i])); continue; } 		// random seed
		if ( strcmp(argv[i],"-memory") == 0 ) { inst->available_memory = atoi(argv[++i]); continue; }	// available memory (in MB)
		if ( strcmp(argv[i],"-max_nodes") == 0 ) { inst->max_nodes = atoi(argv[++i]); continue; } 		// max n. of nodesfile
		if ( strcmp(argv[i],"-callback") == 0) { inst->callback = atoi(argv[++i]); continue; }			// 1 = lazy_callback, 2 = generic_callback
		if ( strcmp(argv[i],"-v") == 0 ) { inst->verbose = atoi(argv[++i]); continue; } 		// max n. of nodes
		if ( strcmp(argv[i],"-float") == 0 ) { inst->integer_costs = 0; continue; } 						// inteher costs
		if ( strcmp(argv[i],"-help") == 0 ) { help = 1; continue; } 									// help
		if ( strcmp(argv[i],"--help") == 0 ) { help = 1; continue; } 									// help
		help = 1;
  }

	if (inst->setup_model == 8) {
		inst->max_fr = max_fr;		// maximum fixing_ratio
		inst->incr_fr = incr_fr;		// increase of fixing_ratio when good gap
		inst->decr_fr = decr_fr;		// decreasing of fixing_ratio when bad gap
		inst->good_gap = good_gap;			// a good gap allow to decrease fixing_ratio
		inst->optimal_gap = optimal_gap;
	}


	if ( help || (inst->verbose >= 2000) )		// print current parameters
	{
		printf("\nAvailable parameters:\n");
		printf("-file %s\n", inst->input_file);
		printf("-time_limit %f\n", inst->timelimit);
		printf("-setup_model %d\n", inst->setup_model);
		printf("-randomseed %d\n", inst->randomseed);
		printf("-max_nodes %d\n", inst->max_nodes);
		printf("-memory %d\n", inst->available_memory);
		printf("-float %d\n", !inst->integer_costs);
		printf("-verbose %d\n", inst->verbose);
		printf("---------------------------------\n\n");
	}

	if ( help ) exit(1);

}

void free_instance(tspinstance *inst) {

	free(inst->xcoord);
	free(inst->ycoord);
	free(inst->best_sol);
	// free(inst->load_min);
	// free(inst->load_max);
}


// Plot functions
void plot_instance(tspinstance *inst) {

	// open gnuplot process
	FILE *gnuplot;
	#ifdef _WIN32
			// printf("Windows\n");
		if(inst->verbose < 100)
			gnuplot = _popen("D:\\Programmi\\gnuplot\\bin\\gnuplot.exe", "w");
		else
			gnuplot = _popen("D:\\Programmi\\gnuplot\\bin\\gnuplot.exe -persist", "w");

			if (gnuplot != NULL){
				fprintf(gnuplot, "set term wx\n");    // set the terminal
				fprintf(gnuplot, "plot '-' with lines\n");  // plot type
				for (int i = 0; i < 10; i++)    // loop over the data [0,...,9]
					fprintf(gnuplot, "%d\n", i);    // data terminated with \n
				fprintf(gnuplot, "%s\n", "e");     // termination character
				fflush(gnuplot);        // flush the gnuplot
			}
	#elif __linux__
			// printf("Linux\n");
			gnuplot = popen("gnuplot", "w");
	#else
			gnuplot = popen("gnuplot", "w");
	#endif

	if (gnuplot == NULL) {
		printf("Could not open gnuplot");
		return;
	}

	char pngname[sizeof(inst->input_file)+20+sizeof(inst->setup_model)];
	char name[sizeof(inst->input_file)];
	snprintf(pngname, sizeof(pngname),
		"plot%c%s_%d.png",
		DIR_DELIM,
		get_file_name(inst->input_file, name),
		inst->setup_model);  // TODO: input_file check

	// set up line and point style depending on the model
	setup_style(gnuplot, inst);

	// set title
	fprintf(gnuplot, "set title \"%s\"\n", pngname);


	// start plotting points
	//plot_points(gnuplot, pngname, inst);

	// save png into FILE
	if (inst->verbose < 100) fprintf(gnuplot, "set terminal png\nset output '%s'\n", pngname);

	// start plotting edges
	plot_edges(gnuplot, pngname, inst);

	// show plot or save in file and close
	fflush(gnuplot);

	#ifdef _WIN32
		//Sleep(2000);		//Don't need in windows if gnuplot doesn't have -persist param
		fclose(gnuplot);
	#else
		if (inst->verbose >= 1000)
			sleep(2); // pause execution to see the plot
		fclose(gnuplot);
	#endif
}

void setup_style(FILE *gnuplot, tspinstance *inst) {

	switch (inst->model_type) {
		case -1:
		  break;
		case 17:
		case 0:			// Line
			setup_linestyle2(gnuplot);
			break;

		case 1:			// Arrow
			setup_arrowstyle2(gnuplot);
			break;

		case 2:			// Arrow
			setup_arrowstyle2(gnuplot);
			break;

		case 3:			// Arrow
			setup_arrowstyle2(gnuplot);
			break;

		default:
			fclose(gnuplot);
			print_error(" Model type unknown!!\n");
			break;
	}

}

void setup_linestyle1(FILE *gnuplot) {
	fprintf(gnuplot,"set style line 1 \
									lc rgb '#0060ad' \
									pointtype 7 \
									pointsize 1.0\n");
}

void setup_linestyle2(FILE *gnuplot) {
	fprintf(gnuplot,"set style line 2 \
									lc rgb '#0060ad' \
									linetype 1 linewidth 2 \
									pointtype 7 pointsize 1.0\n");
}

void setup_arrowstyle2(FILE *gnuplot) {
	fprintf(gnuplot,"set style line 2 \
									lc rgb '#0060ad' \
									linetype 1 linewidth 2 \
									pointtype 7 pointsize 1.0\n\
									set style arrow 2 \
									head filled size screen 0.025,30,45 \
									ls 2\n");
}

void plot_points(FILE *gnuplot, char *pngname, tspinstance *inst) {
	fprintf(gnuplot, "plot '-' w p ls 1\n");
	for (size_t i = 0; i < inst->nnodes; i++)
		fprintf(gnuplot, "%f %f\n", inst->xcoord[i], inst->ycoord[i]);
	fprintf(gnuplot, "e\n");
}

void plot_edges(FILE *gnuplot, char *pngname, tspinstance *inst) {
	if (inst->nedges > 0) // check if solution is available
	{

		switch (inst->model_type) {
		case -1:
		  break;
		case 0:			// Line
			plot_lines_sym(gnuplot, pngname, inst);
			break;

		case 1:			// Arrow
			plot_arrow_asym(gnuplot, pngname, inst);
			break;

		case 2:			// Arrow
			plot_arrow_asym(gnuplot, pngname, inst);
			break;

		case 3:			// Arrow
			plot_arrow_asym(gnuplot, pngname, inst);
			break;

		default:
			fclose(gnuplot);
			print_error(" model type unknown!!");
			break;
		}

	}
}

void plot_lines_sym(FILE *gnuplot, char *pngname, tspinstance *inst) {
	fprintf(gnuplot, "plot '-' w linespoints linestyle 2\n");
	for (int i = 0; i < inst->nnodes; i++)
		for (int j = i+1; j < inst->nnodes; j++)
			if (0.5 < inst->best_sol[xpos(i,j,inst)] ) // && inst->best_sol[i] < 1.00001)
				fprintf(gnuplot, "%f %f\n%f %f\n\n",
									inst->xcoord[i], inst->ycoord[i],
									inst->xcoord[j], inst->ycoord[j]);
	fprintf(gnuplot, "e\n");
}

void plot_arrow_asym(FILE *gnuplot, char *pngname, tspinstance *inst) {
	fprintf(gnuplot, "plot '-' using 1:2:3:4 with vectors arrowstyle 2\n");
	for (int i = 0; i < inst->nnodes; i++)
		for (int j = 0; j < inst->nnodes; j++)
			if (i != j)
				if (0.5 < inst->best_sol[asym_xpos(i,j,inst)] ) // && inst->best_sol[i] < 1.00001)
					fprintf(gnuplot, "%f %f %f %f\n",
										inst->xcoord[i], inst->ycoord[i],
										inst->xcoord[j]-inst->xcoord[i], inst->ycoord[j]-inst->ycoord[i]);
	fprintf(gnuplot, "e\n");
}


char * get_file_name(char *path, char *name) {	// data/att48.tsp -> att48
	strcpy(name, path);
	int start_name = 0;
	for (int i = 0; name[i] != '\0'; i++) {
		if (name[i] == DIR_DELIM) start_name = i + 1;
		if (name[i] == '.') {
			name[i] = '\0';
			break;
		}
	}
	name += start_name;
	return name;
}

// Saving data
int save_results(tspinstance *inst, char *f_name) {
	FILE *outfile;
	outfile = fopen(f_name, "a");
	char dataToAppend[sizeof(inst->input_file)+sizeof(inst->nnodes)*4+ sizeof(inst->opt_time) + 20];
	char name[sizeof(inst->input_file)];
	switch (inst->setup_model) {
		case 8:
			snprintf(dataToAppend, sizeof(dataToAppend),
				"%s; %d; %s_%.3lf_%.3lf_%.3lf_%.3lf_%.3lf; %d; %d; %lf; %lf\n",
				get_file_name(inst->input_file, name),
				inst->nnodes,
				model_name(inst->setup_model),
				inst->max_fr, inst->incr_fr, inst->decr_fr,
				inst->good_gap, inst->optimal_gap,
				inst->randomseed, inst->nthread, inst->opt_time, inst->best_lb);
			break;

		default:
			snprintf(dataToAppend, sizeof(dataToAppend),
								"%s; %d; %s; %d; %d; %lf; %lf\n",
								get_file_name(inst->input_file, name), inst->nnodes,
								model_name(inst->setup_model),
								inst->randomseed, inst->nthread, inst->opt_time, inst->best_lb);
			break;
	}

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
	return 0;
}

int save_res(char * input_file, char * user_name, char * method, double opt_time, double best_lb) {
	FILE *outfile;
	char * f_name = "res.csv";
	outfile = fopen(f_name, "a");
	char dataToAppend[sizeof(input_file)+ sizeof(user_name)+ sizeof(method)+ sizeof(opt_time)+sizeof(best_lb) + 20];
	snprintf(dataToAppend, sizeof(dataToAppend),
						"%s; %s; %s; %lf; %lf\n",
						input_file, user_name, method, opt_time, best_lb);
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
	return 0;
}


// DEBUG only
void print_succ(int* succ, tspinstance* inst ) {
	printf("\ni:      "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", i);
	printf("\nsucc:   "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", succ[i]);
	printf("\n");
	fflush(NULL);
}

void pause_execution() {
	printf("Paused execution. Return to restart: ");
	getchar();
	printf("\n");
}

void print_error(const char *err) {
	printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1);
}
