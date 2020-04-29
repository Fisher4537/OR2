#include "pch.h"
#include "tsp.h"
#include "chrono.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

/**
	TODO:	- a way to get CPLEX error code: status?
			- subtour_iter_opt better performance with index-value?
*/

int TSPopt(tspinstance *inst) {

	// open cplex model
	int status;
	CPXENVptr env = CPXopenCPLEX(&status);
	CPXLPptr lp = CPXcreateprob(env, &status, "TSP");

	// Cplex's parameter setting
	CPXsetintparam(env,CPXPARAM_Threads, inst->nthread);		// allow executing N parallel threads
	CPXsetintparam(env,CPX_PARAM_RANDOMSEED, inst->randomseed);		// avoid performace variability
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit);		// set time limit
	if (inst->verbose >= 100) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);	// show CPLEX log
    // CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 4);

	// CPX_PARAM_MIPEMPHASIS: it balances optimality and integer feasibility.
	//	(CPX_MIPEMPHASIS_BALANCED, CPX_MIPEMPHASIS_FEASIBILITY, CPX_MIPEMPHASIS_OPTIMALITY,
	//	CPX_MIPEMPHASIS_BESTBOUND, CPX_MIPEMPHASIS_HIDDENFEAS)
	// CPX_PARAM_MIPSEARCH: Dynamic search or B&C ?
	// CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 4);

	// set input data in CPX structure
	build_model(inst, env, lp);

	// setup struct to save solution
	inst->nedges = CPXgetnumcols(env, lp);
	inst->best_sol = (double *) calloc(inst->nedges, sizeof(double)); 	// all entries to zero
	inst->zbest = CPX_INFBOUND;

	if (inst->verbose >= 100) printf("build model succesfully.\n");
	if (inst->verbose >= 1000) pause_execution();
	if (inst->verbose >= 100) printf("optimizing model...\n");

	// compute cplex and calculate opt_time w.r.t. OS used
	#ifdef __linux__
		double ini = second();
		mip_optimization(env, lp, inst, &status);
		double fin = second();
		inst->opt_time = (double)(fin - ini);
	#elif _WIN32

		/* Metodo 1 -> QueryPerformanceCounter*/
		double ini = second();
		mip_optimization(env, lp, inst, &status);
		double fin = second();
		inst->opt_time = (double)(fin - ini);
		
		/* Metodo 2	-> CPXgettime
		struct timespec ts, ts2;
		CPXgettime(env, &ts);
		mip_optimization(env, lp, inst, &status);
		CPXgettime(env, &ts2);
		inst->opt_time = (double)ts2.tv_nsec + 1.0e-9 * ((double)ts2.tv_nsec) - (double)ts.tv_nsec + 1.0e-9 * ((double)ts.tv_nsec);
		*/
	#else
		clock_t init_time = clock();
		mip_optimization(env, lp, inst, &status);
		inst->opt_time = (double)(clock() - init_time) / CLOCKS_PER_SEC;
	#endif

	if (inst->verbose >= 100) printf("optimization complete!\n");

	// get best solution
	if (inst->verbose >= 100) printf("getting succ and comp...\n");
	// CPXgetbestobjval(env, lp, &inst->best_lb);
	CPXsolution(env, lp, &status, &inst->best_lb, inst->best_sol, NULL, NULL, NULL);

	// use the optimal solution found by CPLEX
	int *succ = (int*) calloc(inst->nnodes, sizeof(int));
	int *comp = (int*) calloc(inst->nnodes, sizeof(int));
	int *ncomp = (int*) calloc(1, sizeof(int));
	build_sol(inst, succ, comp, ncomp);
	if (inst->verbose >= 100) printf("free instance object...\n");

	// free and close cplex model
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return !status; // status 0 is ok
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
		case 0 :		// basic model with asymmetric x and q
			build_sym_std(inst, env,lp);
			break;

		case 1 :		// MTZ contraints
		 	build_mtz(inst, env,lp);
			break;

		case 2 : 		// FLOW 1
			build_model_flow1(inst, env, lp);
			break;

		case 3: 		// FLOW 1
			build_model_flow1(inst, env, lp);
			break;

		case 4: 		// MTZ with LAZY
			build_model_mtz_lazy(inst, env, lp);
			break;

		case 5:			// STD with HEUR subtour
			build_sym_std(inst, env, lp);
			break;

		default:
			print_error(" model type unknown!!");
			break;
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

	if ( inst->verbose >= -100 ) CPXwriteprob(env, lp, "model/model.lp", NULL);

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

	if ( inst->verbose >= -100 ) CPXwriteprob(env, lp, "model/asym_model.lp", NULL);

	free(cname[0]);
	free(cname);
}

void build_flow1(tspinstance *inst, CPXENVptr env, CPXLPptr lp) {
	char xctype = 'B';	// type of variable
	double obj; // objective function constant
	double lb;	// lower bound
	double ub;	// upper bound

	char **cname = (char **) calloc(1, sizeof(char *));	// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));			// name of the variable

	// add binary constraints and objective
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
	ub = inst->nnodes-1;

	for ( int i = 0; i < inst->nnodes; i++)
	{
		for (int j = 0; j < inst->nnodes; j++)
		{
			if ( i != j )
			{
				sprintf(cname[0], "y(%d,%d)", i+1, j+1);
				if (j == 0)	// y_i1
				{
					if ( CPXnewcols(env, lp, 1, &obj, &lb, &lb, &xctype, cname) )
						print_error(" wrong CPXnewcols on u var.s");
				}
				else				// y_ij
				{
					if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &xctype, cname) )
						print_error(" wrong CPXnewcols on u var.s");
				}

				if ( CPXgetnumcols(env,lp)-1 != asym_ypos(i, j, inst) )
						print_error(" wrong position for u var.s");
			}
		}
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

		// income vertex
		if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) // create a new row
			print_error(" wrong CPXnewrows [degree]");
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

		if ( CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname) ) // create a new row
			print_error(" wrong CPXnewrows [degree]");
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			if ( i == h ) continue;
			if ( CPXchgcoef(env, lp, lastrow, asym_xpos(h,i, inst), 1.0) ) // income vertex
				print_error(" wrong CPXchgcoef [degree]");
		}
	}
}

void build_model_flow1(tspinstance* inst, CPXENVptr env, CPXLPptr lp) {
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
	double coef;
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

			for (int i = 1; i < inst->nnodes; i++) {

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

			for (int i = 0; i < inst->nnodes; i++) {
				if (i != h) {
					if (CPXchgcoef(env, lp, lastrow, asym_ypos(h, i, inst), -1.0))	// outcome
						print_error(" wrong CPXchgcoef [flow(i)]");
				}
			}

			for (int i = 0; i < inst->nnodes; i++) {
				if (i != h) {
					if (CPXchgcoef(env, lp, lastrow, asym_ypos(i, h, inst), 1.0))		// income
						print_error(" wrong CPXchgcoef [flow(i)]");
				}
			}
		}
	}

	rhs = 0.0;
	sense = 'L';
	for (int i = 0; i < inst->nnodes; i++) {
		if (i == 0) {
			for (int j = 1; j < inst->nnodes; j++) {
				lastrow = CPXgetnumrows(env, lp);

				sprintf(cname[0], "y_cut(%d,%d)", i + 1, j + 1);

				// create a new row
				if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))	// new row
					print_error(" wrong CPXnewrows [y_cut()]");

				if (CPXchgcoef(env, lp, lastrow, asym_ypos(i, j, inst), 1.0))
					print_error(" wrong CPXchgcoef [y_cut()]");

				if (CPXchgcoef(env, lp, lastrow, asym_xpos(i, j, inst), - (double) inst->nnodes + 1.0))
					print_error(" wrong CPXchgcoef [y_cut()]");
			}
		}
		else {
			for (int j = 1; j < inst->nnodes; j++) {
				if (i != j) {
					lastrow = CPXgetnumrows(env, lp);

					sprintf(cname[0], "y_cut(%d,%d)", i + 1, j + 1);

					// create a new row
					if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))	// new row
						print_error(" wrong CPXnewrows [y_cut()]");

					if (CPXchgcoef(env, lp, lastrow, asym_ypos(i, j, inst), 1.0))
						print_error(" wrong CPXchgcoef [y_cut()]");

					if (CPXchgcoef(env, lp, lastrow, asym_xpos(i, j, inst), -(double)inst->nnodes + 2.0))
						print_error(" wrong CPXchgcoef [y_cut()]");
				}
			}
		}
	}
}

void build_model_mtz_lazy(tspinstance* inst, CPXENVptr env, CPXLPptr lp) {

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
	int big_M = inst->nnodes;
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


// optimization methods run the problem optimization
void mip_optimization(CPXENVptr env, CPXLPptr lp, tspinstance *inst, int *status) {

	switch (inst->model_type)
	{

		case 0:
			*status = subtour_iter_opt(env, lp, inst, status);
			break;

		case 1:
			*status = CPXmipopt(env,lp);
			break;

		case 2:
			*status = CPXmipopt(env,lp);
			break;

		case 3:
			*status = CPXmipopt(env, lp);
			break;
		
		case 4:
			*status = CPXmipopt(env, lp);
			break;

		case 5:
			*status = subtour_heur_iter_opt(env, lp, inst, status, 0);
			break;

		default:
			print_error("model ì_type not implemented in optimization method");
			break;
	} // switch
}

int subtour_iter_opt(CPXENVptr env, CPXLPptr lp, tspinstance *inst, int *status) {

	// structure init
	int *succ = (int*) calloc(inst->nnodes, sizeof(int));
	int *comp = (int*) calloc(inst->nnodes, sizeof(int));
	int *ncomp = (int*) calloc(1, sizeof(int));
	char sense = 'L';
	double rhs;
	int first;
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
							if ( *status = CPXchgcoef(env, lp, lastrow, xpos(i,j, inst), 1.0) ) // income vertex
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
	return 0;
}

int subtour_heur_iter_opt(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status, int heuristic) {

	

	// structure init
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int* ncomp = (int*)calloc(1, sizeof(int));
	char sense = 'L';
	double rhs;
	int first;
	int lastrow;

	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));			// name of the variable
	*ncomp = inst->nnodes;
	int subtour_counter = 0;


	if (heuristic == 0) {
		CPXsetintparam(env, CPX_PARAM_NODELIM, 0);
		
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
								if (*status = CPXchgcoef(env, lp, lastrow, xpos(i, j, inst), 1.0)) // income vertex
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
		CPXsetintparam(env, CPX_PARAM_NODELIM, 9999999);	//CPX_INFBOUND

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
								if (*status = CPXchgcoef(env, lp, lastrow, xpos(i, j, inst), 1.0)) // income vertex
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
	return 0;
}


void switch_callback(tspinstance* inst, CPXENVptr env, CPXLPptr lp) {
	CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);				// let MIP callbacks work on the original model

	if (inst->callback == 1) {										// Lazy Constraint Callback
		CPXsetlazyconstraintcallbackfunc(env, lazycallback, inst);

	}else if (inst->callback == 2) {								// Generic Callback
		CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE, genericcallback, inst);		// CPX_CALLBACKCONTEXT_CANDIDATE can be used with other params with |
	}
	int ncores = 1;
	CPXgetnumcores(env, &ncores);
	CPXsetintparam(env, CPX_PARAM_THREADS, ncores); 				// it was reset after callback
	inst->ncols = CPXgetnumcols(env, lp);
}

static int CPXPUBLIC lazycallback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p){
	*useraction_p = CPX_CALLBACK_DEFAULT;
	tspinstance* inst = (tspinstance*)cbhandle; 			// casting of cbhandle (which is pointing to the above parameter inst)

	// get solution xstar
	double* xstar = (double*)malloc(inst->ncols * sizeof(double));
	if (CPXgetcallbacknodex(env, cbdata, wherefrom, xstar, 0, inst->ncols - 1))			// xstar = current x from CPLEX-- xstar starts from position 0 (getx is not defined)
		return 1;

	// get some random information at the node (as an example)
	double objval = CPX_INFBOUND; CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);					//valore rilassamento continuo (LB) al nodo corrente
	int mythread = -1; CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &mythread);
	double zbest; CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &zbest);			//valore incumbent al nodo corrente

	//apply cut separator and possibly add violated cuts
	//int ncuts = mylazy_separation(inst, xstar, env, cbdata, wherefrom);
	free(xstar);													//avoid memory leak

	//if (ncuts >= 1)
		*useraction_p = CPX_CALLBACK_SET; 		// tell CPLEX that cuts have been created
	return 0; 												// return 1 would mean error --> abort Cplex's execution
}

static int CPXPUBLIC genericcallback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* cbhandle) {

	if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE) {
		tspinstance* inst = (tspinstance*)cbhandle; 			// casting of cbhandle (which is pointing to the above parameter inst)
				// get solution xstar
		double* xstar = (double*)malloc(inst->ncols * sizeof(double));
		double objval = CPX_INFBOUND;
		if (CPXcallbackgetcandidatepoint(context, xstar, 0, inst->ncols - 1, &objval))			// xstar = current x from CPLEX-- xstar starts from position 0 (getx is not defined)
			return 1;

		int mythread = -1; CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADS, &mythread);
		double zbest; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &zbest);				//valore incumbent al nodo corrente

		//apply cut separator and possibly add violated cuts
		int ncuts = mygeneric_separation(inst, xstar, context);
		free(xstar);											//avoid memory leak

		return 0; 												// return 1 would mean error --> abort Cplex's execution
	}
	else
		return 0;
}

int mygeneric_separation(tspinstance* inst, const double* xstar, CPXCALLBACKCONTEXTptr context) {
	return 0;
}

int mylazy_separation(tspinstance* inst, const double* xstar, CPXCALLBACKCONTEXTptr context) {
	return 0;
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

		case 2:
			build_sol_flow1(inst, succ, comp, ncomp);
			break;
		
		case 3:
			build_sol_flow1(inst, succ, comp, ncomp);
			break;

		case 4:
			build_sol_mtz(inst, succ, comp, ncomp);
			break;

		case 5:		
			build_sol_sym(inst, succ, comp, ncomp);
			break;

		default:
			print_error(" model type unknown!!");
			break;
	}
}

void build_sol_sym(tspinstance *inst, int *succ, int *comp, int *ncomp) {	// build succ() and comp() wrt xstar()...

	// check if nodes degree is 2 for each node
	if (inst->verbose >= 1000)
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
				snprintf(msg, sizeof msg, "wrong degree[%d] = %d in build_sol_sym", i, degree[i]);
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
	if (inst->verbose >= 1000)
	{
		printf("\ni:      ");
		for (int i = 0; i < inst->nnodes; i++) printf("%6d", i);
		printf("\nsucc:   ");
		for (int i = 0; i < inst->nnodes; i++) printf("%6d", succ[i]);
		printf("\ncomp:   ");
		for (int i = 0; i < inst->nnodes; i++) printf("%6d", comp[i]);
		printf("\n");
	}
}

void build_sol_mtz(tspinstance *inst, int *succ, int *comp, int *ncomp) {	// build succ() and comp() wrt xstar()...

	// check if nodes degree is 2 for each node
	if (inst->verbose >= 1000) {
		printf("debugging sym mtz solution...\n");
		if (inst->verbose >= 1000) {
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
				snprintf(msg, sizeof msg, "wrong degree in build_sol_mtz: degree(%d) = %d", i, degree[i]);
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
	if (inst->verbose >= 1000) {
		printf("\ni:   ");
		for (int i = 0; i < inst->nnodes; i++) printf("%5d", i);
		printf("\nsucc:");
		for (int i = 0; i < inst->nnodes; i++) printf("%5d", succ[i]);
		printf("\ncomp:");
		for (int i = 0; i < inst->nnodes; i++) printf("%5d", comp[i]);
		printf("\n");
	}
}

void build_sol_flow1(tspinstance *inst, int *succ, int *comp, int *ncomp) {	// build succ() and comp() wrt xstar()...

	// check if nodes degree is 2 for each node
	if (inst->verbose >= 1000) {
		printf("debugging sym flow1 solution...\n");
		if (inst->verbose >= 1000) {
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
				snprintf(msg, sizeof msg, "wrong degree in build_sol_mtz: degree(%d) = %d", i, degree[i]);
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
	if (inst->verbose >= 1000) {
		printf("\ni:   ");
		for (int i = 0; i < inst->nnodes; i++) printf("%5d", i);
		printf("\nsucc:");
		for (int i = 0; i < inst->nnodes; i++) printf("%5d", succ[i]);
		printf("\ncomp:");
		for (int i = 0; i < inst->nnodes; i++) printf("%5d", comp[i]);
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
			if ( inst->verbose >= 2000 ) printf("read_input: solving instance: %s\b with model %d...\n", token1, inst->model_type);
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
			if ( inst->verbose >= 2000 ) printf("\tnnodes: %d\n", inst->nnodes);
			inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double)*inst->nnodes);
			inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double)*inst->nnodes);
			active_section = 0;
			continue;
		}

		if ( strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0 )
		{
			token1 = strtok(NULL, " :");
			if ( strncmp(token1, "ATT", 3) != 0 && strncmp(token1, "EUC_2D", 6) != 0 )
				print_error(" format error:  only EDGE_WEIGHT_TYPE == EUC_2D implemented so far!!!!!!");
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
			if ( inst->verbose >= 1000 )
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
	inst->model_type = 0;
	inst->randomseed = 0;
	inst->nthread = 1;
	inst->timelimit = 300.; 	// CPX_INFBOUND;
	strcpy(inst->input_file, "NULL");

	inst->available_memory = 12000;   // available memory, in MB, for Cplex execution (e.g., 12000)
	inst->max_nodes = -1; 						// max n. of branching nodes in the final run (-1 unlimited)
  inst->integer_costs = 0;
	inst->verbose = 1000;							// VERBOSE
    inst->callback = 0;

  int help = 0; if ( argc < 1 ) help = 1;
	for ( int i = 1; i < argc; i++ )
  {
		if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-input") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-f") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 				// input file
		if ( strcmp(argv[i],"-time_limit") == 0 ) { inst->timelimit = atof(argv[++i]); continue; }		// total time limit
		if ( strcmp(argv[i],"-model_type") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 	// model type
		if ( strcmp(argv[i],"-model") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 			// model type
		if ( strcmp(argv[i],"-randomseed") == 0 ) { inst->randomseed = abs(atoi(argv[++i])); continue; } 		// random seed
		if ( strcmp(argv[i],"-nthread") == 0 ) { inst->nthread = abs(atoi(argv[++i])); continue; } 		// random seed
		if ( strcmp(argv[i],"-memory") == 0 ) { inst->available_memory = atoi(argv[++i]); continue; }	// available memory (in MB)
		if ( strcmp(argv[i],"-node_file") == 0 ) { strcpy(inst->node_file,argv[++i]); continue; }		// cplex's node file
		if ( strcmp(argv[i],"-max_nodes") == 0 ) { inst->max_nodes = atoi(argv[++i]); continue; } 		// max n. of nodesfile
        if (strcmp(argv[i], "-callback") == 0) { inst->callback = atoi(argv[++i]); continue; }			// 1 = lazy_callback, 2 = generic_callback
		if ( strcmp(argv[i],"-v") == 0 ) { inst->verbose = atoi(argv[++i]); continue; } 		// max n. of nodes
		if ( strcmp(argv[i],"-int") == 0 ) { inst->integer_costs = 1; continue; } 						// inteher costs
		if ( strcmp(argv[i],"-help") == 0 ) { help = 1; continue; } 									// help
		if ( strcmp(argv[i],"--help") == 0 ) { help = 1; continue; } 									// help
		help = 1;
  }

	if ( help || (inst->verbose >= 1000) )		// print current parameters
	{
		printf("\nAvailable parameters:\n");
		printf("-file %s\n", inst->input_file);
		printf("-time_limit %f\n", inst->timelimit);
		printf("-model_type %d\n", inst->model_type);
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

	char pngname[sizeof(inst->input_file)+20+sizeof(inst->model_type)];
	snprintf(pngname, sizeof pngname, "plot\\%s_%d.png", get_file_name(inst->input_file), inst->model_type);  // TODO: input_file check

	// set up line and point style
	setup_style(gnuplot, inst);

	// set title
	char settitle[sizeof(pngname)+20];
	fprintf(gnuplot, "set title \"%s\"\n", pngname);


	// start plotting points
	//plot_points(gnuplot, pngname, inst);

	// save png into FILE
	if (inst->verbose < 1000) fprintf(gnuplot, "set terminal png\nset output '%s'\n", pngname);

	// start plotting edges
	plot_edges(gnuplot, pngname, inst);

	// show plot or save in file and close
  fflush(gnuplot);
	if (inst->verbose >= 1000) getchar(); // pause execution to see the plot
	fclose(gnuplot);
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

			case 0:
				fprintf(gnuplot, "plot '-' w linespoints linestyle 2\n");
				for (int i = 0; i < inst->nnodes; i++)
					for (int j = i+1; j < inst->nnodes; j++)
						if (0.5 < inst->best_sol[xpos(i,j,inst)] ) // && inst->best_sol[i] < 1.00001)
							fprintf(gnuplot, "%f %f\n%f %f\n\n",
												inst->xcoord[i], inst->ycoord[i],
												inst->xcoord[j], inst->ycoord[j]);
				fprintf(gnuplot, "e\n");
				break;

			case 1:
				fprintf(gnuplot, "plot '-' using 1:2:3:4 with vectors arrowstyle 2\n");
				for (int i = 0; i < inst->nnodes; i++)
					for (int j = 0; j < inst->nnodes; j++)
						if (i != j)
							if (0.5 < inst->best_sol[asym_xpos(i,j,inst)] ) // && inst->best_sol[i] < 1.00001)
								fprintf(gnuplot, "%f %f %f %f\n",
													inst->xcoord[i], inst->ycoord[i],
													inst->xcoord[j]-inst->xcoord[i], inst->ycoord[j]-inst->ycoord[i]);
				fprintf(gnuplot, "e\n");
				break;

			case 2:
				fprintf(gnuplot, "plot '-' using 1:2:3:4 with vectors arrowstyle 2\n");
				for (int i = 0; i < inst->nnodes; i++)
					for (int j = 0; j < inst->nnodes; j++)
						if (i != j)
							if (0.5 < inst->best_sol[asym_xpos(i,j,inst)] ) // && inst->best_sol[i] < 1.00001)
								fprintf(gnuplot, "%f %f %f %f\n",
													inst->xcoord[i], inst->ycoord[i],
													inst->xcoord[j]-inst->xcoord[i], inst->ycoord[j]-inst->ycoord[i]);
				fprintf(gnuplot, "e\n");
				break;

			default:
				fclose(gnuplot);
				print_error(" model type unknown!!");
				break;
		}

	}
}

void setup_style(FILE *gnuplot, tspinstance *inst) {

	fprintf(gnuplot,"set style line 1 \
									lc rgb '#0060ad' \
									pointtype 7 \
									pointsize 1.0\n");

	switch (inst->model_type) {

		case 0:
			fprintf(gnuplot,"set style line 2 \
									    lc rgb '#0060ad' \
											linetype 1 linewidth 2 \
											pointtype 7 pointsize 1.0\n");
			break;

		case 1:
			fprintf(gnuplot,"set style line 2 \
									    lc rgb '#0060ad' \
											linetype 1 linewidth 2 \
											pointtype 7 pointsize 1.0\n\
											set style arrow 2 \
											head filled size screen 0.025,30,45 \
											ls 2\n");
			break;

		case 2:
			fprintf(gnuplot, "set style line 2 \
									    lc rgb '#0060ad' \
											linetype 1 linewidth 2 \
											pointtype 7 pointsize 1.0\n\
											set style arrow 2 \
											head filled size screen 0.025,30,45 \
											ls 2\n");
			break;
		
		default:
			fclose(gnuplot);
			print_error(" Model type unknown!!\n");
			break;
	}

}

char * get_file_name(char *path) {
    int start_name = 0;
	for (int i = 0; path[i] != '\0'; i++) {
		#ifdef _WIN32
			if (path[i] == '\\') start_name = i + 1;
		#else
				if (path[i] == '/') start_name = i + 1;
		#endif
	}
	return path+start_name;
}


// Saving data
int save_results(tspinstance *inst, char *f_name) {
	FILE *outfile;
	outfile = fopen(f_name, "a");
	char dataToAppend[sizeof(inst->input_file)+sizeof(inst->nnodes)*4+ sizeof(inst->opt_time) + 20];
	snprintf(dataToAppend, sizeof dataToAppend, "%s; %d; %d; %d; %d; %lf;\n",
	 												inst->input_file, inst->nnodes,
													inst->model_type, inst->randomseed,
													inst->nthread, inst->opt_time );
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


// DEBUG only
void pause_execution() {
	printf("Paused execution. Return to restart: ");
	getchar();
	printf("\n");
}

void print_error(const char *err) {
	printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1);
}


