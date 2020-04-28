#ifndef TSP_H_

#define TSP_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <cplex.h>

#ifdef _WIN32
#else
	#include <pthread.h>
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


int TSPopt(tspinstance *inst);
int xpos(int i, int j, tspinstance *inst);
int asym_xpos(int i, int j, tspinstance *inst);
int asym_upos(int i, tspinstance *inst);

void build_model(tspinstance *inst, CPXENVptr env, CPXLPptr lp); 		// interface
void build_sym_std(tspinstance *inst, CPXENVptr env, CPXLPptr lp); 	// sym, std
void build_mtz(tspinstance *inst, CPXENVptr env, CPXLPptr lp); 			// asym, MTZ
void build_flow1(tspinstance *inst, CPXENVptr env, CPXLPptr lp);
void build_model_flow1(tspinstance* inst, CPXENVptr env, CPXLPptr lp);

void mip_optimization(CPXENVptr env, CPXLPptr lp, tspinstance *inst, int *status);
int subtour_iter_opt(CPXENVptr env, CPXLPptr lp, tspinstance *inst, int *status);

void switch_callback(tspinstance* inst, CPXENVptr env, CPXLPptr lp);
static int CPXPUBLIC lazycallback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p);
static int CPXPUBLIC genericcallback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* cbhandle);
int mygeneric_separation(tspinstance* inst, const double* xstar, CPXCALLBACKCONTEXTptr context);
int mylazy_separation(tspinstance* inst, const double* xstar, CPXCALLBACKCONTEXTptr context);

void build_sol(tspinstance *inst, int *succ, int *comp, int *ncomp);
void build_sol_sym(tspinstance *inst, int *succ, int *comp, int *ncomp);
void build_sol_mtz(tspinstance *inst, int *succ, int *comp, int *ncomp);

void parse_command_line(int argc, char** argv, tspinstance *inst);
void read_input(tspinstance *inst);
void free_instance(tspinstance *inst);

double dist(int i, int j, tspinstance *inst);		// get distance between two nodes

int save_results(tspinstance *inst, char *f_name);	// Save model performance

void plot_instance(tspinstance *inst);
char * get_file_name(char *path);
void setup_style(FILE *gnuplot, tspinstance *inst);
void plot_points(FILE *gnuplot, char *pngname, tspinstance *inst);
void plot_edges(FILE *gnuplot, char *pngname, tspinstance *inst);

// Debug functions
void pause_execution();
void print_error(const char *err);

#endif   /* TSP_H_ */
