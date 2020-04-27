#ifndef TSP_H_

#define TSP_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <cplex.h>
#include <pthread.h>
#ifdef _WIN32
    #include <windows.h>
#endif

//hard-wired parameters
#define XSMALL		  		  1e-5 		    // 1e-4*	// tolerance used to decide
                                      //    ingerality of 0-1 var.s
#define EPS 1e-5
#define TICKS_PER_SECOND 	  1000.0  	// cplex's ticks on Intel Core i7 quadcore
                                      //@2.3GHZ

//data structures

typedef struct {

	//input data
	int nnodes;
	double *xcoord;
	double *ycoord;
	int verbose;

	// parameters
	int model_type;         // TSP
	int randomseed;
	int nthread;						// number of threads
	double timelimit;				// overall time limit, in sec.s
	char input_file[100];	// input file
	char node_file[1000];		// cplex node file
	int available_memory;
	int max_nodes; 					// max n. of branching nodes in the final run
                          // (-1 unlimited)
	// double cutoff; 					// cutoff (upper bound) for master
	int integer_costs;

	//global data
	double opt_time;					// optimization time
	double	tstart;
	double zbest;							// best sol. available
	double tbest;							// time for the best sol. available
	double *best_sol;						// best sol. available
	double	best_lb;						// best lower bound available
	int nedges;									// number of edges in best_sol
	// double *load_min;						// minimum load when leaving a node
	// double *load_max;						// maximum load when leaving a node
    
	double ncols;
	int callback;

	// model;
	// int xstart;
	// int qstart;
	// int bigqstart;
	// int sstart;
	// int bigsstart;
	// int ystart;
	// int fstart;
	// int zstart;
} tspinstance;

//inline
inline int imax(int i1, int i2) { return ( i1 > i2 ) ? i1 : i2; }
inline double dmin(double d1, double d2) { return ( d1 < d2 ) ? d1 : d2; }
inline double dmax(double d1, double d2) { return ( d1 > d2 ) ? d1 : d2; }

void switch_callback(tspInstance* inst, CPXENVptr env, CPXLPptr lp);
static int CPXPUBLIC lazycallback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p);
static int CPXPUBLIC genericcallback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* cbhandle);
int mygeneric_separation(tspInstance* inst, const double* xstar, CPXCALLBACKCONTEXTptr context);
int mylazy_separation(tspInstance* inst, const double* xstar, CPXCALLBACKCONTEXTptr context);
void build_sol_flow1(tspInstance* inst, result* res, int* succ, int* comp, int* ncomp);

#endif   /* TSP_H_ */
