#include "pch.h"
#include "tsp.h"
#include "chrono.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>


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
			- User_cut di Concorde (difficile implementarle e non è possibile utilizzarle ogni volta bensì solo ogni tot..)		-> per la lode

			- funzione della distanza corrette a seconda di EDGE_WEIGHT_TYPE. vedi tsplib/doc.ps

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
		case 8: return "hard_fixing";
		case 9: return "local_branching";
		case 10: return "heuristic_greedy";
		case 11: return "heuristic_greedy_cgal";
		case 12: return "heuristic_grasp";
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
 4																						patching
 5
 6
*/

	inst->callback = 0;
	inst->heuristic = 0;
	inst->warm_start = 0;
	inst->mip_opt = -1;

	switch (inst->setup_model) {
		case 0:
			inst->model_type = 0;
			inst->mip_opt = 0;
			return "subtour";							// basic model with asymmetric x and q
		case 1:
			inst->model_type = 1;
			inst->mip_opt = 2;
			return "mtz";								// MTZ contraints
		case 2:
			inst->model_type = 2;
			inst->mip_opt = 2;
			return "flow1_n-2";							// FLOW 1 with y_0j <= x_0j*(n-2) if i != 0
		case 3:
			inst->model_type = 2;
			inst->mip_opt = 2;
			return "flow1_n-1";							// FLOW 1 with y_0j <= x_0j*(n-1)
		case 4:
			inst->model_type = 3;
			inst->mip_opt = 2;
			return "mtz_lazy";							// MTZ with LAZY constraints
		case 5:
			inst->model_type = 0;
			inst->mip_opt = 1;
			return "subtour_ffi";						// Subtour with fast first incumb
		case 6:
			inst->model_type = 0;
			inst->callback = 1;
			inst->mip_opt = 2;
			return "subtour_callback_lazy";				// Subtour_callback_lazy
		case 7:
			inst->model_type = 0;
			inst->callback = 2;
			inst->mip_opt = 2;
			return "subtour_callback_general";			// Subtour_callback_general
		case 8:
			inst->model_type = 0;
			inst->heuristic = 1;
			inst->callback = 2;
			inst->mip_opt = 2;
			return "hard_fixing";						// Hard-Fixing
		case 9:
			inst->model_type = 0;
			inst->heuristic = 2;
			inst->callback = 2;
			inst->mip_opt = 2;
			return "local_branching";					// Soft-Fixing => Local Branching
		case 10:
			inst->model_type = 0;
			inst->warm_start = 1;
			inst->mip_opt = 100;
			return "heuristic_greedy";					// Heuristic Greedy (no CPLEX)
		case 11:
			inst->model_type = 0;
			inst->warm_start = 2;
			inst->mip_opt = 100;
			return "heuristic_greedy_cgal";				// Heuristic Greedy CGAL (no CPLEX)
		case 12:
			inst->model_type = 0;
			inst->warm_start = 3;
			inst->mip_opt = 100;
			return "heuristic_grasp";					// Heuristic GRASP (no CPLEX)
		case 13:
			inst->model_type = 0;
			inst->warm_start = 4;
			inst->mip_opt = 100;
			return "heuristic_insertion";				// Heuristic Insertion (no CPLEX)
		case 14:
			inst->warm_start = 3;
			inst->model_type = 0;
			inst->heuristic = 3;
			return "grasp_best_two_opt";				// GRASP + best_two_opt
		case 15:
			inst->warm_start = 0;
			inst->model_type = 0;
			inst->heuristic = 4;
			return "patching";

		default: return "not_supported";
	}
}

int TSPopt(tspinstance *inst) {

	// open cplex model
	int status;
	CPXENVptr env = CPXopenCPLEX(&status);
	CPXLPptr lp = CPXcreateprob(env, &status, "TSP");

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

	// set all the parameters of model chosen
	setup_model(inst);

	// set input data in CPX structure
	build_model(inst, env, lp);

	// setup struct to save solution
	inst->nedges = CPXgetnumcols(env, lp);
	inst->best_sol = (double *) calloc(inst->nedges, sizeof(double)); 	// all entries to zero
	inst->zbest = CPX_INFBOUND;

	// set callback if selected
	switch_callback(inst, env, lp);

	// set warm start if used
	switch_warm_start(inst, env, lp, &status);



	if (inst->verbose >= 100) printf("\nbuild model succesfully.\n");
	if (inst->verbose >= 100) printf("optimizing model...\n");

	// compute cplex and calculate opt_time w.r.t. OS used
	double ini = second();
	//tmp_heur_grasp(inst, &status); // TODO save best_val inside
	//best_two_opt(inst);
	optimization(env, lp, inst, &status);
	double fin = second();
	inst->opt_time = (double)(fin - ini);


	if (inst->verbose >= 100) printf("optimization complete!\n");

	if(inst->setup_model != 9 && inst->setup_model != 14)
		CPXsolution(env, lp, &status, &inst->best_lb, inst->best_sol, NULL, NULL, NULL);

	if (inst->verbose >= 100) printf("free instance object...\n");

	// free and close cplex model
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return 0; // status 0 is ok
}

int xpos(int i, int j, tspinstance *inst) {
	if ( i == j ) print_error(" i == j in xpos" );
	if ( i > j ) return xpos(j,i,inst);								// simplify returned formula
	return i*inst->nnodes + j - ((i + 1)*(i + 2))/2; 				// default case
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


// build_model methods add the constraints to OPTIMIZER structures
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
				rhs = big_M -1;
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


void switch_warm_start(tspinstance* inst, CPXENVptr env, CPXLPptr lp, int* status) {

	switch (inst->warm_start) {					// No CPLEX used

		case 0:
		break;

		case 1:													// Heuristic greedy
			heur_greedy(env, lp, inst, status);
		break;

		case 2:													// Heuristic greedy CGAL
			heur_greedy_cgal(env, lp, inst, status);
		break;

		case 3:													// Heuristic GRASP (no CPLEX)
			heur_grasp(inst, status); // do not add the solution to CPLEX
			// int izero = 0;
			// double val = 1.0;
			// int nocheck_warmstart = CPX_MIPSTART_CHECKFEAS;  // CPX_MIPSTART_SOLVEFIXED;
			// if ( (*status = CPXaddmipstarts(env, lp, 1, inst->nnodes, &izero, inst->best_sol, &val, &nocheck_warmstart, NULL)) ) {
			// 	print_error("Error during warm start: adding new start, check CPXaddmipstarts\n");
			// }
		break;

		case 4:
			heur_insertion(env, lp, inst, status);				// Heuristic Insertion
		break;

		default:
			print_error(" model type unknown!!");
		break;
	}
}
void test_warm_start(CPXENVptr env, CPXLPptr lp) {
	CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);

	CPXmipopt(env, lp);
}


// constructive heuristic
int heur_greedy_cgal(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status) {

	#ifdef _WIN32
		CPXsetintparam(env, CPX_PARAM_ADVIND, 1);

		double best_lb = CPX_INFBOUND;
		inst->best_lb = CPX_INFBOUND;
		double val = 1.0;
		int* best_sol = (int*)calloc(inst->nnodes, sizeof(int));
		int izero = 0;

		int nocheck_warmstart = CPX_MIPSTART_SOLVEFIXED;	// CPLEX solves the fixed problem specified by the MIP start (requires to provide values for all discrete variables)
								//CPX_MIPSTART_NOCHECK;		// CPLEX accepts the MIP start without any checks. The MIP start needs to be complete.

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
			free(sol);
		}
		free_cgal();
		if (inst->verbose >= 100)
			printf("BEST_LB Greedy Heuristic CGAL found: [%f]\n", inst->best_lb);

		if (CPXaddmipstarts(env, lp, 1, inst->nnodes, &izero, best_sol, &val, &nocheck_warmstart, NULL)) {
			print_error("Error during warm start: adding new start, check CPXaddmipstarts\n");
			return *status;
		}

		/********* Add a mip start

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
	#endif

	return *status;
}

int heur_greedy(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status) {

	CPXsetintparam(env, CPX_PARAM_ADVIND, 1);

	double best_lb = CPX_INFBOUND;
	double val = 1.0;
	inst->best_lb = CPX_INFBOUND;
	int* best_sol = (int*)calloc(inst->nnodes, sizeof(int));
	int izero = 0;

	int nocheck_warmstart = CPX_MIPSTART_SOLVEFIXED;	// CPLEX solves the fixed problem specified by the MIP start (requires to provide values for all discrete variables)
							//CPX_MIPSTART_NOCHECK;		// CPLEX accepts the MIP start without any checks. The MIP start needs to be complete.


	for (int i = 0; i < inst->nnodes; i++) {

		int* sol = (int*)calloc(inst->nnodes, sizeof(int));
		for (int k = 0; k < inst->nnodes; k++) {
			sol[k] = -1;
		}

		int succ = succ_not_contained(i, sol, inst);
		sol[0] = i;
		sol[1] = succ;
		if(inst->verbose > 1000)
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
			if (inst->verbose > 100)
				printf("%d,%d\n", sol[j], sol[j + 1]);
			best_lb += dist(sol[j], sol[j + 1], inst);
			sol[j] = xpos(sol[j], sol[j + 1], inst);

		}
		if (inst->verbose > 100)
			printf("%d,%d\n\n", sol[inst->nnodes - 1], i);
		best_lb += dist(sol[inst->nnodes - 1], i, inst);
		sol[inst->nnodes - 1] = xpos(sol[inst->nnodes - 1], i, inst);

		if (inst->verbose >= 100)
			printf("Solution: %d\tBEST_LB found: [%f]\n",i, best_lb);
		if (best_lb < inst->best_lb) {
			if (inst->verbose >= 100)
				printf("BEST_LB update from -> to : [%f] -> [%f]\n", inst->best_lb, best_lb);
			inst->best_lb = best_lb;
			for (int k = 0; k < inst->nnodes; k++)
				best_sol[k] = sol[k];
		}
		free(sol);
	}

	if (inst->verbose >= 100)
		printf("BEST_LB Greedy Heuristic found: [%f]\n", inst->best_lb);

	if (CPXaddmipstarts(env, lp, 1, inst->nnodes, &izero, best_sol, &val, &nocheck_warmstart, NULL)) {
		print_error("Error during warm start: adding new start, check CPXaddmipstarts\n");
		return *status;
	}
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

void heur_grasp(tspinstance* inst, int* status) {
	if (inst->verbose >= 100) printf("heur_grasp helloworld!\n");

	double* best_sol = (double*)calloc(inst->nedges, sizeof(double));

	// get first node, randomly selected
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	for (int i = 0; i < inst->nnodes; i++) succ[i] = -1;

	int cur_node = round(((double)rand()/RAND_MAX)*(inst->nnodes-1));
	int first_node = cur_node;

	// find the three nearest node of each node
	int cur_nearest[3];  		// index of the nearest from current node
	double cur_dist[3];  		// distance from current node
	double best_lb = 0.0;  	// the cost of the solution

	int tour_length = 0;
	while (tour_length < inst->nnodes) {

		// initialization
		for (int k = 0; k < 3; k++) { cur_dist[k] = INT_MAX; cur_nearest[k] = -1; }

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
		best_sol[xpos(cur_node, next_node, inst)] = 1.0;
		best_lb += next_dist;
		cur_node = next_node;
		tour_length++;
	}

	// save the tour and the cost // TODO: does it really save the solution?!
	print_succ(succ, inst);
	free(succ);
	for (int i = 0; i < inst->nnodes; i++)
	 	for (int j = i+1; j < inst->nnodes; j++)
			inst->best_sol[xpos(i,j,inst)] = best_sol[xpos(i,j,inst)];
	inst->best_lb = best_lb;
	if (inst->verbose >= 100) printf("GRASP BEST_LB: %lf\n", inst->best_lb);

}

int heur_insertion(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status) {
	int* best_sol = (int*)calloc(inst->nnodes, sizeof(int));
	for (int k = 0; k < inst->nnodes; k++)
		best_sol[k] = -1;

	srand(inst->randomseed);

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

	int izero = 0;
	double val = 1.0;
	int nocheck_warmstart = CPX_MIPSTART_CHECKFEAS;

	if ((*status = CPXaddmipstarts(env, lp, 1, inst->nnodes, &izero, best_sol, &val, &nocheck_warmstart, NULL))) {
		print_error("Error during warm start: adding new start, check CPXaddmipstarts\n");
		return *status;
	}

	CPXsetintparam(env, CPX_PARAM_ADVIND, 1);

	return 0;
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
	best_sol[replace_pos] = xpos(best_i, vertex, inst);
	best_sol[count_sol] = xpos(vertex, best_j, inst);
	if (inst->verbose > 10)
		printf("***** Sides %d, %d added in best_sol (actual length %d) *****\n\n", best_sol[replace_pos], best_sol[count_sol], count_sol + 1);
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
			best_two_opt(inst);
		break;

		case 4:
			patching(inst);
		break;

		default:
			print_error("model ì_type not implemented in optimization method");
		break;
	}
}

int hard_fixing(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status) {


	double external_time_limit = inst->timelimit;
	double internal_time_limit = 50.0;
	double init_time = second();
	double next_time_limit = internal_time_limit;
	double remaining_time = external_time_limit - (second()-init_time);
	double objval_p = 0.0;
	double gap = 1.0;						// best_lb - objval_p / best_lb
	double fr = 0.9;				// fixing_ratio


	// structure init
	// int *succ = (int*) calloc(inst->nnodes, sizeof(int));
	// int *comp = (int*) calloc(inst->nnodes, sizeof(int));
	// int *ncomp = (int*) calloc(1, sizeof(int));

	// set next internal time limit
	remaining_time = external_time_limit - (second()-init_time);
	if (remaining_time > internal_time_limit*2)
		next_time_limit = internal_time_limit;
	else
		next_time_limit = remaining_time;
	CPXsetdblparam(env, CPX_PARAM_TILIM, next_time_limit);

	// first optimization step
	mip_optimization(env, lp, inst, status);

	// get first step solution gap
	CPXsolution(env, lp, status, &inst->best_lb, inst->best_sol, NULL, NULL, NULL);
	CPXgetbestobjval(env, lp, &objval_p);
	gap = (inst->best_lb - objval_p) / inst->best_lb;

	// build_sol(inst, succ, comp, ncomp);
	// if (inst->verbose >= 100) printf("Partial solution, ncomp = %d\n",*ncomp );
	if (inst->verbose >= 1000) plot_instance(inst);


	while ( (second()-init_time) < external_time_limit &&
	 (fr > 0.0 || gap > (inst->optimal_gap)))
	{
		// printf("BEST SOLUTION: %lf\n", inst->best_lb);

		// fix a % of bounds
		fix_bound(env, lp, inst, status, fr);

		// set next internal time limit
		remaining_time = external_time_limit - (second()-init_time);
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

		case 1:
			print_error(" model type unknown!!");
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
	// TODO resolve returning value
	int k_index = 0;
	double k[5] = { 3.0, 5.0, 10.0, 15.0, 20.0};
	double timelimit = 300;									// internal timelimit
	double temp_timelimit = timelimit;
	double remaining_time = inst->timelimit;

	CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);			// abort Cplex after the first incument update
	CPXsetdblparam(env, CPX_PARAM_TILIM, timelimit);

	double* best_sol = (double*)calloc(inst->nedges, sizeof(double));
	double best_lb = CPX_INFBOUND;

	mip_optimization(env, lp, inst, status);
	CPXsolution(env, lp, status, &inst->best_lb, inst->best_sol, NULL, NULL, NULL);

	CPXsetintparam(env, CPX_PARAM_INTSOLLIM, INT_MAX);

	for (int h = 0; remaining_time > 0.0; h++) {
		double ini = second();

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
			if (temp_timelimit > remaining_time)
				CPXsetdblparam(env, CPX_PARAM_TILIM, remaining_time);
			else
				CPXsetdblparam(env, CPX_PARAM_TILIM, temp_timelimit);
		}else {
			if (k_index < 4) {
				k_index++;
				rhs = (double)inst->nnodes - k[k_index];
				temp_timelimit = timelimit * (k_index + 1.0);
			}else {
				k[4] = (k[4] * 2.0 > (double)inst->nnodes) ? (double)inst->nnodes : k[4] * 2.0;
				rhs = (double)inst->nnodes - k[4];
			}
			if(temp_timelimit > remaining_time)
				CPXsetdblparam(env, CPX_PARAM_TILIM, remaining_time);
			else
				CPXsetdblparam(env, CPX_PARAM_TILIM, temp_timelimit);

		}

		sprintf(cname[0], "local-branching constraint, k_index = %f", k_index < 5 ? k[k_index] : k[4]);

		char sense = 'G';

		if (inst->verbose >= 100) {
			printf("\n********** Round = %d  -  K = %f  -  Remaining_time = %6.3lf **********\n\n", h, k[k_index], remaining_time);
		}


		for (int i = 0; i < inst->nnodes; i++){
			for (int j = i + 1; j < inst->nnodes; j++) {
				if (inst->best_sol[xpos(i, j, inst)] == 1) {
					index[nnz] = xpos(i, j, inst);
					value[nnz] = 1.0;
					nnz++;
				}
			}
		}
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname))
			print_error("wrong CPXpreaddrows() for adding local-branching constraint\n");

		free(index);
		free(value);

		if (mip_optimization(env, lp, inst, status)) {
			printf("Error in CPXmipopt\n");
		}

		double fin = second();
		remaining_time -= (fin - ini);

		if (remaining_time <= 0.01) {
			if (CPXgetstat(env, lp) == 101 || CPXgetstat(env, lp) == 102) {
				if(best_lb > inst->best_lb)
					inst->best_sol = best_sol;
			}
			if (CPXdelrows(env, lp, CPXgetnumrows(env, lp) - 1, CPXgetnumrows(env, lp) - 1))
				print_error("wrong CPXdelrows() for deleting local-branching constraint\n");

			printf("*** Better soluzion found! ***\n");
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
			best_sol = inst->best_sol;
		}

		if(CPXdelrows(env, lp, CPXgetnumrows(env, lp) - 1, CPXgetnumrows(env, lp) - 1))
			print_error("wrong CPXdelrows() for deleting local-branching constraint\n");

	}
	return 0;
}


// optimization methods run the problem optimization
int mip_optimization(CPXENVptr env, CPXLPptr lp, tspinstance *inst, int *status) {

	switch (inst->mip_opt)
	{

		case 0:			// subtour_iter_opt, symmetric, without callback
			*status = subtour_iter_opt(env, lp, inst, status);
			break;

		case 1:			// Subtour with HEUR
			*status = subtour_heur_iter_opt(env, lp, inst, status, 0);
			break;

		case 2:
			*status = CPXmipopt(env,lp);
			break;

		case 100:
			test_warm_start(env, lp);
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
							if ( (*status = CPXchgcoef(env, lp, lastrow, xpos(i,j, inst), 1.0)) ) // income vertex
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
	}
	else if (inst->callback == 2) {									// Generic Callback
		CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
		if (CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE, genericcallback, inst)) {
			print_error(" Error in setGenericCallback2()\n");
			CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_ON);
			return;
		}
	}
	else if (inst->callback == 3) {
		CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
		if (CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_GLOBAL_PROGRESS, genericcallback, inst)) {
			print_error(" Error in setGenericCallback3()\n");
			CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_ON);
			return;
		}
	}
	else {																// Callback not used
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
void best_two_opt(tspinstance *inst) {
	// improving the best_sol if possibile in 2opt set

	if (inst->verbose >= 100) printf("BEST_TWO_OPT\n");

	// check if current solution has only one tour
	int *succ = (int*) calloc(inst->nnodes, sizeof(int));
	int *comp = (int*) calloc(inst->nnodes, sizeof(int));
	int *ncomp = (int*) calloc(1, sizeof(int));
	build_sol(inst, succ, comp, ncomp);
	print_succ(succ, inst);
	if (*ncomp != 1) print_error("call best_two_opt with best_sol with multiple tour");

	// search in 2opt if a better solution is found
	int i = 0;		// start from node 0
	int j = succ[succ[0]];
	int best_i = i;
	int best_j = i;
	double best_improve = 0.0;
	double d_i1_i2;
	double d_j1_j2;
	double d_i1_j1;
	double d_i2_j2;
	double cur_improve;

	for (int ti = 0; ti < inst->nnodes; ti++) {
		d_i1_i2 = dist(i, succ[i], inst);
		for (int tj = ti+2; tj < inst->nnodes; tj++) {  // no 2opt with consequence arches
			if (i == succ[j]) break;  // should happen only when i = 0,
			d_j1_j2 = dist(j, succ[j], inst);
			d_i1_j1 = dist(i, j, inst);
			d_i2_j2 = dist(succ[i], succ[j], inst);
			cur_improve = (d_i1_i2 + d_j1_j2) - (d_i1_j1 + d_i2_j2);
			if ( cur_improve > best_improve ) { // cross is better
				best_i = i;
				best_j = j;
				best_improve = cur_improve;
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
		succ[best_i] = best_j;
		succ[i2] = j2;
		int cur_node = i2;
		int suc_node = j2;
		while (cur_node != best_j) {  // reverse the succ
			suc_node = cur_node;
			cur_node = pre_node;
			pre_node = succ[pre_node];
			succ[cur_node] = suc_node;
			if (inst->verbose >= 100) print_succ(succ, inst);
		}
	}
	inst->best_lb -= best_improve;  // update best_lb

	if (inst->verbose >= 100) print_succ(succ, inst);
	free(succ);
	free(comp);
	free(ncomp);
}


// Repair
void patching(tspinstance* inst) {

	if (inst->verbose >= 100) printf("PATCHING\n");

	// check if current solution has only one tour
	int *succ = (int*) calloc(inst->nnodes, sizeof(int));
	int *comp = (int*) calloc(inst->nnodes, sizeof(int));
	int *ncomp = (int*) calloc(1, sizeof(int));
	build_sol(inst, succ, comp, ncomp);
	print_succ(succ, inst);
	if (*ncomp == 1) {
		printf("WARNING: solution already has 1 tour, patching has no effect.\n");
		free(succ);
		free(comp);
		free(ncomp);
		return;
	}

	while (*ncomp > 1) {
		single_patch(inst, succ, comp, ncomp);
		print_succ(succ, inst);
		printf("\ncomp:   "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", comp[i]);
		plot_instance(inst);
	}
}
void single_patch(tspinstance* inst, int* succ, int* comp, int* ncomp) {

	if (*ncomp == 1) {
		return;
	} else {
		int initial_i = rand() % inst->nnodes;  // randomize the first node

		// get the closer node which is in another tour
		int closer_j = initial_i;
		double cn_dist = INT_MAX;
		double d_ij = 0.;
		for (int j = 0; j < inst->nnodes; j++) {
			if (comp[initial_i] == comp[j]) continue; // skip if it's in the same components

			d_ij = dist(initial_i, j, inst);
			if (d_ij < cn_dist) {
				closer_j = j;
				cn_dist = d_ij;
			}
		}
		// find the node of the i-tour which is closer to closer_j: closer_i
		int i = initial_i;
		int closer_i = i;
		int counter = 0;
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


			while (i != initial_i) {  // for each node of the i-tour (the tour which contain i)
				if (succ[closer_j] > 0) {  // j is not isolated
					cur_improve = dist(i, succ[closer_j], inst) + dist(succ[i], closer_j, inst)
											- dist(i, succ[i], inst) - dist(closer_j, succ[closer_j], inst);
				} else  {  // j is isolated
					cur_improve = dist(i, closer_j, inst) + dist(closer_j, succ[i], inst)
															- dist(i, succ[i], inst);
				}
				if (cur_improve < best_improve) {
					best_improve = cur_improve;
					closer_i = i;
					cn_dist = dist(closer_i, closer_j, inst);
				}
				counter++;
				i = succ[i];
			}
		}

		// merge the tours: update best_sol, succ, comp, ncomp, best_lb
 		if (succ[closer_j] < 0 && succ[closer_i] < 0) { // closer_j and closer_i are isolated
			// update inst->best_sol
			(inst->best_sol)[xpos(closer_i, closer_j, inst)] = 1.;

			// update succ
			succ[closer_j] = closer_i;
			succ[closer_i] = closer_j;

			// update comp
			comp[closer_j] = comp[closer_i];

			// update inst->best_lb
			inst->best_lb += dist(closer_i, closer_j, inst);
		} else if (succ[closer_j] < 0) { 		// closer_j isolated node
			// update inst->best_sol
			if (closer_i != succ[succ[closer_i]])
				(inst->best_sol)[xpos(closer_i, succ[closer_i], inst)] = 0.;
			(inst->best_sol)[xpos(closer_i, closer_j, inst)] = 1.;
			(inst->best_sol)[xpos(closer_j, succ[closer_i], inst)] = 1.;

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

			// update inst->best_lb
			inst->best_lb += cn_dist
										+ dist(closer_j, succ[closer_j], inst)
										- dist(closer_i, succ[closer_j], inst);
		} else if (succ[closer_i] < 0) {  		// i is isolated
			// update inst->best_sol
			if (closer_j != succ[succ[closer_j]])
				(inst->best_sol)[xpos(closer_j, succ[closer_j], inst)] = 0.;
			(inst->best_sol)[xpos(closer_i, closer_j, inst)] = 1.;
			(inst->best_sol)[xpos(closer_i, succ[closer_j], inst)] = 1.;

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

			// update inst->best_lb
			inst->best_lb += cn_dist
										+ dist(closer_i, succ[closer_i], inst)
										- dist(closer_j, succ[closer_i], inst);
		} else {  										// i and closer_j are not isolated
			// update inst->best_sol
			if (closer_j != succ[succ[closer_j]])
				(inst->best_sol)[xpos(closer_j, succ[closer_j], inst)] = 0.;
			if (closer_i != succ[succ[closer_i]])
				(inst->best_sol)[xpos(closer_i, succ[closer_i], inst)] = 0.;
			(inst->best_sol)[xpos(closer_i, succ[closer_j], inst)] = 1.;
			(inst->best_sol)[xpos(closer_j, succ[closer_i], inst)] = 1.;

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

			// update inst->best_lb
			inst->best_lb += dist(closer_j, succ[closer_j], inst)
										+ dist(closer_i, succ[closer_i], inst)
										- dist(closer_j, succ[closer_i], inst)
										- dist(closer_i, succ[closer_j], inst);
		}
		(*ncomp)--;
	}
}

int is_clockwise(tspinstance *inst, int x1, int x2, int x3) {
	return (inst->xcoord[x3] - inst->xcoord[x1])*(inst->ycoord[x2] - inst->ycoord[x1]) <
				 (inst->ycoord[x3] - inst->ycoord[x1])*(inst->xcoord[x2] - inst->xcoord[x1]);
}

// build_sol methods use the optimized solution to plot
void build_sol(tspinstance *inst, int *succ, int *comp, int *ncomp) {

	switch (inst->model_type)
	{
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

	// check if nodes degree is 2 for each node
	if (inst->verbose >= 2000)
	{
		int *degree = (int *) calloc(inst->nnodes, sizeof(int));
		printf("nnodes=%d\n", inst->nnodes );
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			for ( int j = i+1; j < inst->nnodes; j++ )
			{
				int k = xpos(i,j,inst);
				if ( fabs(inst->best_sol[k]) > EPS && fabs(inst->best_sol[k]-1.0) > EPS ) print_error(" wrong inst->best_sol in build_sol()");
				if ( inst->best_sol[k] > 0.5 )
				{
					printf("x[%d,%d] = 1\n", i, j );
					++degree[i];
					++degree[j];
				} else {
					printf("x[%d,%d] = 0\n", i, j );
				}

			}
		}
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			if ( degree[i] != 2 )
			{
				char msg[40];
				snprintf(msg, sizeof(msg), "wrong degree[%d] = %d in build_sol_sym", i, degree[i]);
				print_error(msg);
			}
		}
		free(degree);
	}

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
		while ( succ[i] == -1  )  // go and visit the current component
		{
			comp[i] = *ncomp;
			// done = 1;
			for ( int j = 0; j < inst->nnodes; j++ ) {
				if (j == i) continue;

				if ( inst->best_sol[xpos(i,j,inst)] > 0.5) {  // the edge [i,j] is selected in inst->best_sol and j was not visited before
					// intern edge of the cycle
					if (comp[j] == -1) {
						succ[i] = j;
						i = j;
						break;
					}
					// last edge of the cycle
					if (start == j) {
						succ[i] = j;
					}
				}
			}
			if (succ[i] == -1) {
				break;
			}
		}	// while
	// go to the next component...
	}

	// print succ, comp and ncomp
	if (inst->verbose >= 100)
	{
		printf("\ni:      "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", i);
		printf("\nsucc:   "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", succ[i]);
		printf("\ncomp:   "); for (int i = 0; i < inst->nnodes; i++) printf("%6d", comp[i]);
		printf("\n");
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
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j];
	if ( !inst->integer_costs ) return sqrt(dx*dx+dy*dy);
	int dis = sqrt(dx*dx+dy*dy) + 0.499999999; // nearest integer
	return dis+0.0;
}


// User input
void read_input(tspinstance *inst) { // simplified CVRP parser, not all SECTIONs detected

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
	inst->integer_costs = 0;
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
		if ( strcmp(argv[i],"-node_file") == 0 ) { strcpy(inst->node_file,argv[++i]); continue; }		// cplex's node file
		if ( strcmp(argv[i],"-max_nodes") == 0 ) { inst->max_nodes = atoi(argv[++i]); continue; } 		// max n. of nodesfile
    if ( strcmp(argv[i],"-callback") == 0) { inst->callback = atoi(argv[++i]); continue; }			// 1 = lazy_callback, 2 = generic_callback
		if ( strcmp(argv[i],"-v") == 0 ) { inst->verbose = atoi(argv[++i]); continue; } 		// max n. of nodes
		if ( strcmp(argv[i],"-int") == 0 ) { inst->integer_costs = 1; continue; } 						// inteher costs
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
		printf("-int %d\n", inst->integer_costs);
		printf("-verbose %d\n", inst->verbose);
		printf("-node_file %s\n", inst->node_file);
		printf("---------------------------------\n\n");
	}

	if ( help ) exit(1);

}

void free_instance(tspinstance *inst) {

	free(inst->xcoord);
	free(inst->ycoord);
	// free(inst->load_min);
	// free(inst->load_max);
}


// Plot functions
void plot_instance(tspinstance *inst) {

	// open gnuplot process
	FILE *gnuplot;
	#ifdef _WIN32
			// printf("Windows\n");
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
	if (inst->verbose < 1000) fprintf(gnuplot, "set terminal png\nset output '%s'\n", pngname);

	// start plotting edges
	plot_edges(gnuplot, pngname, inst);

	// show plot or save in file and close
	fflush(gnuplot);
	if (inst->verbose >= 1000)
		#ifdef _WIN32
			Sleep(2000); // pause execution to see the plot
		#else
			sleep(2); // pause execution to see the plot
		#endif
	fclose(gnuplot);
}

void setup_style(FILE *gnuplot, tspinstance *inst) {

	switch (inst->model_type) {

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
