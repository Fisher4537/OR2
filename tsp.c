#include "pch.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "tsp.h"

/*** TODO:
		- migliora commenti
		- controlla mtz																=> a posto?
		- sistema plot																=> DONE (grazie fish)

		- non generare vincolo vuoto e cambiare poi i coef 
		  => genera tutto index value e poi passali così
		- look LOOP and performing profile throught LOOP and lazycallback
		- old-lazycallback - new-lazycallback (per non perdere il dynamic search) 						
		  generic-callback															=> da sistemare	

***/


/******************************************************************************************************/
int TSPopt(tspInstance* inst, result* res) {
	// open cplex model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

	if(inst->verbose >= 80)
		CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);

	CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 4);
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->randomseed);

	// build the chosen model
	switch_model(inst, env, lp);

	// setup struct to save solution
	res->nedges = CPXgetnumcols(env, lp);
	res->best_sol = (double*)calloc(res->nedges, sizeof(double)); 	// all entries to zero
	res->zbest = CPX_INFBOUND;

	// Cplex's parameter setting
	printf("build model succesfully.\n");
	printf("optimizing model...\n");

	// compute cplex
	clock_t init_time = clock();
	mip_optimization(env, lp, inst, res, &error);
	res->opt_time = ((double)clock() - (double)init_time) / CLOCKS_PER_SEC;
	printf("optimization complete!\n");

	// get best solution
	CPXsolution(env, lp, &error, &res->best_lb, res->best_sol, NULL, NULL, NULL);

	// use the optimal solution found by CPLEX
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int* ncomp = (int*)calloc(1, sizeof(int));
	build_sol(inst, res, succ, comp, ncomp);


	// free and close cplex model
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return !error;
}

/******************************************************************************************************/
int xpos(int i, int j, tspInstance* inst){
	if (i == j) 
		print_error(" i == j in xpos");
	if (i > j) 
		return xpos(j, i, inst);									// recall xpos with correct indexes
	return i * inst->nnodes + j - ((i + 1) * (i + 2)) / 2; 			// default case
}

/**********************************	********************************************************************/
int asymmetric_xpos(int i, int j, tspInstance* inst){
	if (i == j) 
		print_error(" i == j in asymmetric_xpos");
	return i * (inst->nnodes - 1) + (i < j ? j - 1 : j);
}

/******************************************************************************************************/
int asymmetric_upos(int i, tspInstance* inst){
	if (i < 1) 
		print_error(" i < 1 in asymmetric_upos");
	return inst->nnodes * (inst->nnodes - 1) + i - 1;
}

/******************************************************************************************************/
int asymmetric_ypos(int i, int j, tspInstance* inst) {
	if (i == j)
		print_error(" i == j in asymmetric_ypos");
	return i * (inst->nnodes - 1) + (i < j ? j - 1 : j);
}

/******************************************************************************************************/
void switch_model(tspInstance* inst, CPXENVptr env, CPXLPptr lp) {

	switch (inst->model_type){
		case 0:									// basic model with asymmetric x and q
			build_model_std(inst, env, lp);
			if (inst->callback >= 1) {			// implementation of callbacks
				switch_callback(inst, env, lp);
			}
		break;
		case 1:									// MTZ contraints
			build_model_mtz(inst, env, lp);
		break;
		case 2: 								// FLOW1
			build_model_flow1(inst, env, lp);
		break;
		default:
			print_error(" model type unknown!!");
		break;
	}
}

/******************************************************************************************************/
void build_model_std(tspInstance* inst, CPXENVptr env, CPXLPptr lp){

	char binary = CPX_BINARY;								// Binary 0,1

	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));

	// add binary var.s x(i,j) for i < j
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = i + 1; j < inst->nnodes; j++) {
			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);	// save names of variables
			double obj = dist(i, j, inst);					// cost object function == distance
			double lb = 0.0;
			double ub = 1.0;
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname))
				print_error(" wrong CPXnewcols on x var.s");
			if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, inst))
				print_error(" wrong position for x var.s");
		}
	}

	// add the degree constraints
	for (int h = 0; h < inst->nnodes; h++) {							// degree constraints
		int lastrow = CPXgetnumrows(env, lp);	
		double rhs = 2.0;												// right hand side = termine noto
		char sense = 'E';												// 'E' for equality constraint
		sprintf(cname[0], "degree(%d)", h + 1);
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))			// add empty row
			print_error(" wrong CPXnewrows [degree]");
		for (int i = 0; i < inst->nnodes; i++) {
			if (i == h) 
				continue;
			if (CPXchgcoef(env, lp, lastrow, xpos(i, h, inst), 1.0))	// change coef to 1.0 of the row
				print_error(" wrong CPXchgcoef [degree]");
		}
	}

	if (inst->verbose >= 100) 
		CPXwriteprob(env, lp, "model/std_model.lp", NULL);

	free(cname[0]);
	free(cname);
}

/******************************************************************************************************/
void build_model_mtz(tspInstance* inst, CPXENVptr env, CPXLPptr lp) {

	char xctype = CPX_BINARY;		// Binary 0,1
	double obj;						// objective function constant
	double lb;						// lower bound
	double ub;						// upper bound

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
				if (CPXgetnumcols(env, lp) - 1 != asymmetric_xpos(i, j, inst))
					print_error(" wrong position for x var.s");
			}
		}
	}

	// add nodes index 'ui' in the circuits
	xctype = CPX_INTEGER;						// Integer values
	obj = 0.0;
	lb = 0.0;
	ub = (double) inst->nnodes - 2;

	for (int i = 1; i < inst->nnodes; i++) {			// i=0 -> ui=0	useless 
		sprintf(cname[0], "u(%d)", i + 1);
		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &xctype, cname))
			print_error(" wrong CPXnewcols on u var.s");
		if (CPXgetnumcols(env, lp) - 1 != asymmetric_upos(i, inst))
			print_error(" wrong position for u var.s");
	}

	// add the degree constraints
	int lastrow;	// the number of rows
	double rhs;		// right head size
	char sense;		// 'L', 'E' or 'G'
	double big_M = inst->nnodes - 1.0;				// we can do better? -2.0?	 look pdf compact mtz
 
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
			if (CPXchgcoef(env, lp, lastrow, asymmetric_xpos(i, h, inst), 1.0)) // income vertex
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
			if (CPXchgcoef(env, lp, lastrow, asymmetric_xpos(h, i, inst), 1.0)) // outcome vertex
				print_error(" wrong CPXchgcoef [degree]");
		}


		// u constraints: nodes index
		if (h == 0) 	// skip when i == 0
			continue;
		for (int i = 1; i < inst->nnodes; i++) {
			if (i != h) {
				lastrow = CPXgetnumrows(env, lp);
				rhs = big_M - 1.0;					
				sense = 'L';						// 'L' for lower-equal
				sprintf(cname[0], "mtz_i(%d, %d)", h + 1, i + 1);

				// create a new row
				if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))	// new row
					print_error(" wrong CPXnewrows [degree]");

				if (CPXchgcoef(env, lp, lastrow, asymmetric_upos(h, inst), 1.0)) // u constraints
					print_error(" wrong CPXchgcoef [degree]");
				if (CPXchgcoef(env, lp, lastrow, asymmetric_upos(i, inst), -1.0)) // u constraints
					print_error(" wrong CPXchgcoef [degree]");
				if (CPXchgcoef(env, lp, lastrow, asymmetric_xpos(h, i, inst), big_M)) // u constraints
					print_error(" wrong CPXchgcoef [degree]");
			}
		}
	}

	if (inst->verbose >= 10) 
		CPXwriteprob(env, lp, "model/mtz_model.lp", NULL);

	free(cname[0]);
	free(cname);
}

/******************************************************************************************************/
void build_model_flow1(tspInstance* inst, CPXENVptr env, CPXLPptr lp) {
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
				if (CPXgetnumcols(env, lp) - 1 != asymmetric_xpos(i, j, inst))
					print_error(" wrong position for x var.s");
			}
		}
	}

	// add nodes index 'yij' in the circuits
	xctype = CPX_INTEGER;						// Integer values
	obj = 0.0;
	lb = 0.0;
	ub = (double)inst->nnodes - 1;
	ub1 = (double)inst->nnodes - 2;

	for (int i = 0; i < inst->nnodes; i++) {			
		for (int j = 0; j < inst->nnodes; j++) {
			if (i != j) {
				sprintf(cname[0], "y(%d,%d)", i + 1, j + 1);
				if (j == 0) {														// yi0 = 0
					if (CPXnewcols(env, lp, 1, &obj, &lb, &lb, &xctype, cname))
						print_error(" wrong CPXnewcols on y var.s");
				}else if(i == 0){													// y0j <= n-1			missing * xij
					if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &xctype, cname))
						print_error(" wrong CPXnewcols on y var.s");
				}else {																// yij <= n-2			missing * xij
					if (CPXnewcols(env, lp, 1, &obj, &lb, &ub1, &xctype, cname))
						print_error(" wrong CPXnewcols on y var.s");
				}
				if (CPXgetnumcols(env, lp) - 1 != asymmetric_ypos(i, j, inst))
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
			if (CPXchgcoef(env, lp, lastrow, asymmetric_xpos(i, h, inst), 1.0)) // income vertex
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
			if (CPXchgcoef(env, lp, lastrow, asymmetric_xpos(h, i, inst), 1.0)) // outcome vertex
				print_error(" wrong CPXchgcoef [degree]");
		}

		// y constraints: nodes index
		if (h == 0 ) {
			rhs = (double)inst->nnodes - 1.0;		
			sense = 'E';
			lastrow = CPXgetnumrows(env, lp);

			sprintf(cname[0], "flow(%d)", 1);

			// create a new row
			if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))	// new row
				print_error(" wrong CPXnewrows [flow(0)]");

			for (int i = 1; i < inst->nnodes; i++) {
				
				if (CPXchgcoef(env, lp, lastrow, asymmetric_xpos(h, i, inst), 1.0))		// outcome vertex from 0
					print_error(" wrong CPXchgcoef [degree]");
			}
			// check coloumn
		}else {
			rhs = 1.0;
			sense = 'E';
			lastrow = CPXgetnumrows(env, lp);

			sprintf(cname[0], "flow(%d)", h+1);

			// create a new row
			if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))	// new row
				print_error(" wrong CPXnewrows [flow(1)]");

			for (int i = 0; i < inst->nnodes; i++) {
				if (i != h) {
					if (CPXchgcoef(env, lp, lastrow, asymmetric_xpos(h, i, inst), -1.0))
						print_error(" wrong CPXchgcoef [degree]");
				}
			}

			for (int i = 0; i < inst->nnodes; i++) {
				if (i != h) {
					if (CPXchgcoef(env, lp, lastrow, asymmetric_xpos(i, h, inst), 1.0))
						print_error(" wrong CPXchgcoef [degree]");
				}
			}
			// check coloumn
		}
	}
}

/******************************************************************************************************/
double dist(int i, int j, tspInstance* inst) {
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j];
	if (!inst->integer_costs)
		return sqrt(dx * dx + dy * dy);
	int dis = sqrt(dx * dx + dy * dy) + 0.499999999;			 // nearest integer
	return dis + 0.0;
}

/******************************************************************************************************/
void mip_optimization(CPXENVptr env, CPXLPptr lp, tspInstance* inst, result* res, int* error){
	// optimation methods run the problem optimizationiz

	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));
	int* ncomp = (int*)calloc(1, sizeof(int));
	char sense = 'L';
	double rhs;
	int lastrow;

	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));			// name of the variable

	switch (inst->model_type) {			
		case 0:												// subtour_iter_opt
			*ncomp = inst->nnodes;
			int subtour_counter = 0;

			CPXmipopt(env, lp);
			CPXsolution(env, lp, error, &res->best_lb, res->best_sol, NULL, NULL, NULL);
			build_sol(inst, res, succ, comp, ncomp);
			while (*ncomp >= 2) {

				// subtour elimination constraints
				for (int comp_i = 1; comp_i <= *ncomp; comp_i++) {
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

					for (int i = 0; i < inst->nnodes; i++) {
						if (comp[i] == comp_i) {
							for (int j = i + 1; j < inst->nnodes; j++) {
								if (comp[j] == comp_i) {
									if (*error = CPXchgcoef(env, lp, lastrow, xpos(i, j, inst), 1.0)) // income vertex
										print_error(" wrong CPXchgcoef [degree]");
								}
							}
						}
					}
				}
				subtour_counter++;
				CPXmipopt(env, lp);
				CPXsolution(env, lp, error, &res->best_lb, res->best_sol, NULL, NULL, NULL);
				build_sol(inst, res, succ, comp, ncomp);
				if (inst->verbose >= 50) printf("Iter %3d partial solution, ncomp = %d\n", subtour_counter, *ncomp);
				if (inst->verbose >= 100) get_pipe(&inst, &res);
			}
			if (inst->verbose >= 100)
				printf("best solution found. ncomp = %d\n", *ncomp);
			break;

		case 1:
			*error = CPXmipopt(env, lp);
		break;
		default:
			print_error("model_type not implemented in optimization method");
		break;
	}
}

/******************************************************************************************************/
void build_sol(tspInstance* inst, result* res, int* succ, int* comp, int* ncomp) {
	switch (inst->model_type) {
		case 0 :
			build_sol_std(inst, res, succ, comp, ncomp);
		break;
		case 1 :
			build_sol_mtz(inst, res, succ, comp, ncomp);
		break;
		case 2:
			build_sol_flow1(inst, res, succ, comp, ncomp);
			break;
		default :
			print_error("model_type UNKNOWN!");
		break;
	}
}

/******************************************************************************************************/	
void build_sol_std(tspInstance* inst, const result* res, int* succ, int* comp, int* ncomp){		// build succ() and comp() wrt res->best_sol()

	if (inst->verbose >= 100) {
		int* degree = (int*)calloc(inst->nnodes, sizeof(int));
		printf("nnodes=%d\n", inst->nnodes);

		for (int i = 0; i < inst->nnodes; i++) {
			for (int j = i + 1; j < inst->nnodes; j++) {
				
				int k = xpos(i, j, inst);
				if (fabs(res->best_sol[k]) > EPS && fabs(res->best_sol[k] - 1.0) > EPS)
					print_error(" wrong inst -> best_sol in build_sol()");

				if (res->best_sol[k] > 0.5) {
					printf("x[%d,%d] = 1\n", i, j);
					++degree[i];
					++degree[j];
				}else {
					printf("x[%d,%d] = 0\n", i, j);
				}
			}
		}
		// check if nodes degree is 2 for each node
		for (int i = 0; i < inst->nnodes; i++) {
			if (degree[i] != 2) {
				char msg[40];
				snprintf(msg, sizeof msg, "wrong degree[%d] = %d in build_sol_sym", i, degree[i]);
				print_error(msg);
			}
		}
		free(degree);
	}

	// initialization of succ, comp and ncomp
	*ncomp = 0;
	for (int i = 0; i < inst->nnodes; i++){
		succ[i] = -1;
		comp[i] = -1;
	}

	if (TRUE) {				// FALSE to use alternative method
		// tour
		for (int start = 0; start < inst->nnodes; start++) {
			if (comp[start] >= 0) 								// node "start" has already been setted
				continue;

			(*ncomp)++;											// a new component is found
			comp[start] = *ncomp;

			int i = start;
			while (succ[i] == -1){  							// go and visit the current component
				comp[i] = *ncomp;
				for (int j = 0; j < inst->nnodes; j++) {
					if (j == i) 
						continue;
					if (res->best_sol[xpos(i, j, inst)] > 0.5) { // the edge [i,j] is selected in best_sol and j was not visited before

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
			}
		}
	}else{
		// alternative tour
		for (int start = 0; start < inst->nnodes; start++) {
			if (comp[start] >= 0) 							
				continue;

			(*ncomp)++;										
			int prv = -1;
			int i = start;
			while (comp[start] == -1){													
				for (int j = 0; j < inst->nnodes; j++) {
					if (res->best_sol[xpos(i, j, inst)] > 0.5 && i != j && j != prv) {
					
						succ[i] = j;
						comp[j] = *ncomp;
						prv = i;
						i = j;
						break;
					}
				}
			}
		}
	}

	// print succ, comp and ncomp
	if (inst->verbose >= 100) {
		printf("\ni:      ");
		for (int i = 0; i < inst->nnodes; i++)
			printf("%6d", i);
		printf("\nsucc:   ");
		for (int i = 0; i < inst->nnodes; i++)
			printf("%6d", succ[i]);
		printf("\ncomp:   ");
		for (int i = 0; i < inst->nnodes; i++)
			printf("%6d", comp[i]);
		printf("\n");
	}
}

/******************************************************************************************************/
void build_sol_mtz(tspInstance* inst, result* res, int* succ, int* comp, int* ncomp) {		// build succ() and comp() wrt res->best_sol()

	// check if nodes degree is 2 for each node
	if (inst->verbose >= 100){
		printf("debugging sym mtz solution...\n");

		printf("Solution:\n      ");
		for (int i = 0; i < inst->nnodes; i++) 
			printf("%5d|", i);
		printf("\n");

		for (int i = 0; i < inst->nnodes; i++){
			printf("%5d)", i);
			for (int j = 0; j < inst->nnodes; j++)
				if (i == j) 
					printf("      ");
				else 
					printf("%6.1f", round(res->best_sol[asymmetric_xpos(i, j, inst)]));
			printf("\n");
		}
		printf("\n    u)      ");
		for (int u = 1; u < inst->nnodes; u++)
			printf("%6.1f", round(res->best_sol[asymmetric_upos(u, inst)]));

		int* degree = (int*)calloc(inst->nnodes, sizeof(int));
		for (int i = 0; i < inst->nnodes; i++){
			for (int j = 0; j < inst->nnodes; j++){
				if (i == j) 
					continue;
				int k = asymmetric_xpos(i, j, inst);
				if (fabs(res->best_sol[k]) > EPS && fabs(res->best_sol[k] - 1.0) > EPS) 
					print_error(" wrong inst->best_sol in build_sol()");
				if (res->best_sol[k] > 0.5){
					++degree[i];
					++degree[j];
				}
			}
		}
		for (int i = 0; i < inst->nnodes; i++){

			if (degree[i] != 2){
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
	for (int i = 0; i < inst->nnodes; i++){
		succ[i] = -1;
		comp[i] = -1;
	}

	if (TRUE) {				// FALSE to use alternative method
		// tour
		for (int start = 0; start < inst->nnodes; start++) {
			if (comp[start] >= 0) 								// node "start" has already been setted
				continue;

			(*ncomp)++;											// a new component is found
			comp[start] = *ncomp;

			int i = start;
			while (succ[i] == -1) {  							// go and visit the current component
				comp[i] = *ncomp;
				for (int j = 0; j < inst->nnodes; j++) {
					if (j == i)
						continue;
					if (res->best_sol[asymmetric_xpos(i, j, inst)] > 0.5) { // the edge [i,j] is selected in best_sol and j was not visited before

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
			}
		}
	}
	else {
		// alternative tour
		for (int start = 0; start < inst->nnodes; start++) {
			if (comp[start] >= 0)
				continue;

			(*ncomp)++;
			int prv = -1;
			int i = start;
			while (comp[start] == -1) {
				for (int j = 0; j < inst->nnodes; j++) {
					if (res->best_sol[asymmetric_xpos(i, j, inst)] > 0.5 && i != j && j != prv) {

						succ[i] = j;
						comp[j] = *ncomp;
						prv = i;
						i = j;
						break;
					}
				}
			}
		}
	}

	// print succ, comp and ncomp
	if (inst->verbose >= 100) {
		printf("\ni:      ");
		for (int i = 0; i < inst->nnodes; i++)
			printf("%5d", i);
		printf("\nsucc:   ");
		for (int i = 0; i < inst->nnodes; i++)
			printf("%5d", succ[i]);
		printf("\ncomp:   ");
		for (int i = 0; i < inst->nnodes; i++)
			printf("%5d", comp[i]);
		printf("\n");
	}
}

/******************************************************************************************************/
void build_sol_flow1(tspInstance* inst, result* res, int* succ, int* comp, int* ncomp) {
	
	// check if nodes degree is 2 for each node
	if (inst->verbose >= 100) {
		printf("debugging sym flow1 solution...\n");

		printf("Solution:\n      ");
		for (int i = 0; i < inst->nnodes; i++)
			printf("%5d|", i);
		printf("\n");

		for (int i = 0; i < inst->nnodes; i++) {
			printf("%5d)", i);
			for (int j = 0; j < inst->nnodes; j++)
				if (i == j)
					printf("      ");
				else
					printf("%6.1f", round(res->best_sol[asymmetric_xpos(i, j, inst)]));
			printf("\n");
		}
		printf("\n    y)      ");
		for (int u = 1; u < inst->nnodes; u++)
			printf("%6.1f", round(res->best_sol[asymmetric_upos(u, inst)]));

		int* degree = (int*)calloc(inst->nnodes, sizeof(int));
		for (int i = 0; i < inst->nnodes; i++) {
			for (int j = 0; j < inst->nnodes; j++) {
				if (i == j)
					continue;
				int k = asymmetric_xpos(i, j, inst);
				if (fabs(res->best_sol[k]) > EPS && fabs(res->best_sol[k] - 1.0) > EPS)
					print_error(" wrong inst->best_sol in build_sol()");
				if (res->best_sol[k] > 0.5) {
					++degree[i];
					++degree[j];
				}
			}
		}
		for (int i = 0; i < inst->nnodes; i++) {

			if (degree[i] != 2) {
				char msg[100];
				snprintf(msg, sizeof msg, "wrong degree in build_sol_mtz: degree(%d) = %d", i, degree[i]);
				print_error(msg);
			}
		}
		free(degree);
		printf("\ndebug completed succesfully.\n");
	}
	
	
	if (TRUE) {				// FALSE to use alternative method
		// tour
		for (int start = 0; start < inst->nnodes; start++) {
			if (comp[start] >= 0) 								// node "start" has already been setted
				continue;

			(*ncomp)++;											// a new component is found
			comp[start] = *ncomp;

			int i = start;
			while (succ[i] == -1) {  							// go and visit the current component
				comp[i] = *ncomp;
				for (int j = 0; j < inst->nnodes; j++) {
					if (j == i)
						continue;
					if (res->best_sol[asymmetric_xpos(i, j, inst)] > 0.5) { // the edge [i,j] is selected in best_sol and j was not visited before

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
			}
		}
	}
	else {
		// alternative tour
		for (int start = 0; start < inst->nnodes; start++) {
			if (comp[start] >= 0)
				continue;

			(*ncomp)++;
			int prv = -1;
			int i = start;
			while (comp[start] == -1) {
				for (int j = 0; j < inst->nnodes; j++) {
					if (res->best_sol[asymmetric_xpos(i, j, inst)] > 0.5 && i != j && j != prv) {

						succ[i] = j;
						comp[j] = *ncomp;
						prv = i;
						i = j;
						break;
					}
				}
			}
		}
	}

	// print succ, comp and ncomp
	if (inst->verbose >= 100) {
		printf("\ni:      ");
		for (int i = 0; i < inst->nnodes; i++)
			printf("%5d", i);
		printf("\nsucc:   ");
		for (int i = 0; i < inst->nnodes; i++)
			printf("%5d", succ[i]);
		printf("\ncomp:   ");
		for (int i = 0; i < inst->nnodes; i++)
			printf("%5d", comp[i]);
		printf("\n");
	}
}

/******************************************************************************************************/
void switch_callback(tspInstance* inst, CPXENVptr env, CPXLPptr lp) {
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

/******************************************************************************************************/
static int CPXPUBLIC lazycallback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p){
	*useraction_p = CPX_CALLBACK_DEFAULT;
	tspInstance* inst = (tspInstance*)cbhandle; 			// casting of cbhandle (which is pointing to the above parameter inst)

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

/******************************************************************************************************/
static int CPXPUBLIC genericcallback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* cbhandle) {

	if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE) {
		tspInstance* inst = (tspInstance*)cbhandle; 			// casting of cbhandle (which is pointing to the above parameter inst)
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

/******************************************************************************************************/
int mygeneric_separation(tspInstance* inst, const double* xstar, CPXCALLBACKCONTEXTptr context) {
	int ncomp = 9999;
	int* succ = (int*)calloc(inst->nnodes, sizeof(int));
	int* comp = (int*)calloc(inst->nnodes, sizeof(int));

	build_sol_std(inst, xstar, succ, comp, &ncomp);
	if (ncomp > 1) {
		for (int i = 1; i <= ncomp; i++) {
			int count_nodes = 0;
			char** cname = (char**)calloc(1, sizeof(char*));
			cname[0] = (char**)calloc(1000, sizeof(char));
			strcpy(cname[0], "SEC callback(");
			char str[3];
			sprintf(str, "%d", i);
			strcat(cname[0], str);			//OSCENO cambialo!!!
			strcat(cname[0], ")");
			int* comp_nodes = (int*)calloc(inst->nnodes, sizeof(int));
			if (inst->verbose > 100)
				printf("node in comp %d:\n", i);
			for (int j = 0; j < inst->nnodes; j++) {
				if (comp[j] == i) {
					comp_nodes[count_nodes++] = j;
					if (inst->verbose > 100) {
						printf("- node %d\n", j);
					}
				}
			}
			if (inst->verbose > 100)
				printf("\n");
			int izero = 0;
			int* index = (int*)calloc(count_nodes * count_nodes, sizeof(int));
			double* value = (double*)calloc(count_nodes * count_nodes, sizeof(double));
			int* start_indexes = 0;

			double rhs = (double)(count_nodes - 1);
			char sense = 'L';
			int	nnz = 0;
			for (int j = 0; j < count_nodes; j++) {
				for (int k = 0; k < count_nodes; k++) {
					if (inst->verbose > 100)
						printf("x(%d,%d) +", comp_nodes[j], comp_nodes[k]);
					index[nnz] = xpos(comp_nodes[j], comp_nodes[k], inst);
					value[nnz] = 1.0;
					nnz++;
				}
			}
			if (inst->verbose > 100)
				printf("<=%f\n\n", rhs);
			if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &start_indexes, index, value))
				print_error("errore callback reject");
			free(index);
			free(value);
		}
		if (ncomp > 2 && inst->verbose > 80)
			printf(" %d SEC added using generic callback\n", ncomp);
		return (ncomp == 1 ? 0 : ncomp);
	}
	else
		return 0;
}

/******************************************************************************************************/
int mylazy_separation(tspInstance* inst, const double* xstar, CPXCALLBACKCONTEXTptr context) {
	return 0;
}

/******************************************************************************************************/
void read_input(tspInstance* inst, result* res) {	// simplified CVRP parser, not all SECTIONs detected 

	FILE* fin = fopen(inst->input_file, "r");
	if (fin == NULL)
		print_error(" input file not found!");

	inst->nnodes = -1;
	res->nedges = -1;

	char line[180];
	char* par_name;
	char* token1;
	char* token2;

	int active_section = 0;		// = 1 NODE_COORD_SECTION, = 2 DEMAND_SECTION, = 3 DEPOT_SECTION 

	int do_print = (inst->verbose >= 1000);

	while (fgets(line, sizeof(line), fin) != NULL) {

		if (inst->verbose >= 2000) {
			printf("%s", line);
			fflush(NULL);
		}

		if (strlen(line) <= 1)
			continue; // skip empty lines

		par_name = strtok(line, " :");

		if (inst->verbose >= 3000) {
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
			if (inst->verbose >= 10)
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

			inst->xcoord = (double*)calloc(inst->nnodes, sizeof(double));
			inst->ycoord = (double*)calloc(inst->nnodes, sizeof(double));

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

/******************************************************************************************************/
void parse_command_line(int argc, char** argv, tspInstance* inst) {

	if (inst->verbose >= 100)
		printf(" running %s with %d parameters \n", argv[0], argc - 1);

	// default   
	inst->model_type = 0;
	strcpy(inst->input_file, "NULL");
	strcpy(inst->node_file, "NULL");
	inst->randomseed = 0;
	inst->num_threads = 0;
	inst->timelimit = CPX_INFBOUND;				// CPX_INFBOUND = "infinito" di CPLEX


	inst->available_memory = 12000;   			// available memory, in MB, for Cplex execution (e.g., 12000)
	inst->max_nodes = -1; 						// max n. of branching nodes in the final run (-1 unlimited)        
	inst->integer_costs = 0;
	inst->verbose = 1000;
	inst->callback = 0;

	int help = 0;
	if (argc < 1)
		help = 1;

	for (int i = 1; i < argc; i++) {

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
		if (strcmp(argv[i], "-v") == 0) { inst->verbose = atoi(argv[++i]); continue; }
		if (strcmp(argv[i], "-callback") == 0) { inst->callback = atoi(argv[++i]); continue; }			// 1 = lazy_callback, 2 = generic_callback
		if (strcmp(argv[i], "-int") == 0) { inst->integer_costs = 1; continue; }
		if (strcmp(argv[i], "-help") == 0) { help = 1; continue; } 										// help
		if (strcmp(argv[i], "--help") == 0) { help = 1; continue; } 									// help
		help = 1;
	}

	// print current parameters
	if (help || (inst->verbose >= 10)) {
		printf("\n\navailable parameters (vers. 15-04-2020) --------------------------------------------------\n");			//UPDATE VERSION
		printf("-file %s\n", inst->input_file);
		printf("-time_limit %lf\n", inst->timelimit);
		printf("-model_type %d\n", inst->model_type);
		printf("-seed %d\n", inst->randomseed);
		printf("-threads %d\n", inst->num_threads);
		printf("-max_nodes %d\n", inst->max_nodes);
		printf("-memory %d\n", inst->available_memory);
		printf("-node_file %s\n", inst->node_file);
		printf("-int %d\n", inst->integer_costs);
		printf("-v %d\n", inst->verbose);
		printf("-callback %d\n", inst->callback);
		printf("\nenter -help or --help for help\n");
		printf("----------------------------------------------------------------------------------------------\n\n");
	}

	if (help)
		exit(1);

}

/******************************************************************************************************/
void print_error(const char* err) {
	printf("\n\n ERROR: %s \n\n", err);
	fflush(NULL);
	exit(1);
}

/******************************************************************************************************/
void free_instance(tspInstance* inst) {

	free(inst->xcoord);
	free(inst->ycoord);
	// free(inst->load_min);
	// free(inst->load_max);
}

/******************************************************************************************************/
void get_pipe(tspInstance* inst, result* res) {
	
	FILE* pipe = _popen("D:\\Programmi\\gnuplot\\bin\\gnuplot.exe -persist", "w");

	if (pipe != NULL){
		fprintf(pipe, "set term wx\n");				// set the terminal
		fprintf(pipe, "plot '-' with lines\n");		// plot type
		for (int i = 0; i < 10; i++)				// loop over the data [0,...,9]
			fprintf(pipe, "%d\n", i);				// data terminated with \n
		fprintf(pipe, "%s\n", "e");					// termination character
		fflush(pipe);								// flush the pipe

		getchar();

		plot_problem_input(inst, res, pipe);
	}
	else
		printf("Could not open pipe");
}

void plot_problem_input(tspInstance* inst, result* res, FILE* gnuplot){
	// TODO: consider to use a gnuplot interface if things become complicated
	// open gnuplot process
	char pngname[sizeof(inst->input_file) + 20 + sizeof(inst->model_type)];
	snprintf(pngname, sizeof pngname, "graph/%s_%d.png", inst->input_file, inst->model_type);

	// set up line and point style
	setup_style(gnuplot, inst);

	// set title
	//char settitle[sizeof(pngname) + 20];
	fprintf(gnuplot, "set title \"%s\"\n", pngname);


	// start plotting points
	//plot_points(gnuplot, pngname, inst);

	// save png into FILE
	if (inst->verbose < 1000) // save plot in file
		fprintf(gnuplot, "set terminal png\nset output '%s'\n", pngname);

	// start plotting edges
	plot_edges(gnuplot, pngname, inst, res);

	// show plot or save in file and close
	fflush(gnuplot);
	if (inst->verbose >= 1000)
		getchar(); // pause execution to see the plot
	_pclose(gnuplot);
}

void plot_points(FILE* gnuplot, char* pngname, tspInstance* inst){
	fprintf(gnuplot, "plot '-' w p ls 1\n");
	for (size_t i = 0; i < inst->nnodes; i++)
		fprintf(gnuplot, "%f %f\n", inst->xcoord[i], inst->ycoord[i]);
	fprintf(gnuplot, "e\n");
}

void plot_edges(FILE* gnuplot, char* pngname, tspInstance* inst, result* res){
	if (res->nedges > 0){				// check if solution is available
		switch (inst->model_type) {
			case 0:
				fprintf(gnuplot, "plot '-' w linespoints linestyle 2\n");
				for (int i = 0; i < inst->nnodes; i++)
					for (int j = i + 1; j < inst->nnodes; j++)
						if (0.5 < res->best_sol[xpos(i, j, inst)]) // && inst->best_sol[i] < 1.00001)
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
							if (0.5 < res->best_sol[asymmetric_xpos(i, j, inst)]) // && inst->best_sol[i] < 1.00001)
								fprintf(gnuplot, "%f %f %f %f\n",
									inst->xcoord[i], inst->ycoord[i],
									inst->xcoord[j] - inst->xcoord[i], inst->ycoord[j] - inst->ycoord[i]);
				fprintf(gnuplot, "e\n");
			break;
			default:
				fclose(gnuplot);
				print_error(" model type unknown!!");
			break;
		}

	}
}

void setup_style(FILE* gnuplot, tspInstance* inst){

	fprintf(gnuplot, "set style line 1 \
									lc rgb '#0060ad' \
									pointtype 7 \
									pointsize 1.0\n");

	switch (inst->model_type) {
		case 0:
			fprintf(gnuplot, "set style line 2 \
											lc rgb '#0060ad' \
												linetype 1 linewidth 2 \
												pointtype 7 pointsize 1.0\n");
		break;
		case 1:
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
			print_error(" model type unknown!!");
		break;
	}
}


/******************************************************************************************************/
void save_results(tspInstance* inst, result* res, char* f_name){
	FILE* outfile;
	outfile = fopen(f_name, "a");
	char dataToAppend[200];
	snprintf(dataToAppend, sizeof dataToAppend, "%s; %d; %d; %d; %lf ;\n",
		inst->input_file, inst->nnodes, inst->model_type, 8, res->opt_time); // 8 is number of thread
	/* fopen() return NULL if unable to open file in given mode. */
	if (outfile == NULL){
		/* Unable to open file hence exit */
		printf("\nUnable to open '%s' file.\n", f_name);
		printf("Please check whether file exists and you have write privilege.\n");
		exit(EXIT_FAILURE);
	}

	/* Append data to file */
	fputs(dataToAppend, outfile);

	fclose(outfile);

}
