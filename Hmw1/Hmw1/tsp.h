#ifndef TSP_H_  

#define TSP_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>  

#include <cplex.h>  
#include <windows.h>

											// printing level  (=10 only incumbent,
											//    =20 little output, =50-60 good,
											//    =70 verbose, >=100 cplex log)

//hard-wired parameters
#define XSMALL		  		  1e-5 		    // 1e-4*	tolerance used to decide
											//			ingerality of 0-1 var.s
#define EPS					  1e-5
#define EPSILON		  		  1e-9		    // 1e-9		very small numerical tolerance
#define TICKS_PER_SECOND 	  1000.0  		// cplex's ticks on Intel Core i5-9th quadcore
											// @2.4GHZ up to 4.10GHZ

// data structures 
typedef struct {

	//input data
	int nnodes;
	double *xcoord;
	double *ycoord;

	// parameters 
	int model_type;
	int randomseed;
	int num_threads;
	double timelimit;						// overall time limit, in sec.s
	char input_file[1000];		  			// input file
	char node_file[1000];		  			// cplex node file
	int available_memory;
	int max_nodes; 							// max n. of branching nodes in the final run (-1 unlimited)
	int integer_costs;
	int verbose;
	double ncols;
	int callback;

} tspInstance;

// result structures
typedef struct {

	//global data
	double opt_time;						// optimization time
	double	tstart;
	double zbest;							// best sol. available  
	double tbest;							// time for the best sol. available  
	double *best_sol;						// best sol. available    
	double	best_lb;						// best lower bound available  
	//double *load_min;						// minimum load when leaving a node
	//double *load_max;						// maximum load when leaving a node
	int nedges;

} result;

//inline
inline int imax(int i1, int i2) { return (i1 > i2) ? i1 : i2; }
inline double dmin(double d1, double d2) { return (d1 < d2) ? d1 : d2; }
inline double dmax(double d1, double d2) { return (d1 > d2) ? d1 : d2; }


int TSPopt(tspInstance* inst, result* res);
int xpos(int i, int j, tspInstance* inst);
void switch_model(tspInstance* inst, CPXENVptr env, CPXLPptr lp);								// interface
void build_model_std(tspInstance* inst, CPXENVptr env, CPXLPptr lp);							// sym, std
void build_model_mtz(tspInstance* inst, CPXENVptr env, CPXLPptr lp);
void mip_optimization(CPXENVptr env, CPXLPptr lp, tspInstance* inst, result* res, int* error);
void build_sol(tspInstance* inst, result* res, int* succ, int* comp, int* ncomp);
void build_sol_std(tspInstance* inst, result* res, int* succ, int* comp, int* ncomp);
void build_sol_mtz(tspInstance* inst, result* res, int* succ, int* comp, int* ncomp);
void build_model_flow1(tspInstance* inst, CPXENVptr env, CPXLPptr lp);

void switch_callback(tspInstance* inst, CPXENVptr env, CPXLPptr lp);
static int CPXPUBLIC lazycallback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p);
static int CPXPUBLIC genericcallback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* cbhandle);
int mygeneric_separation(tspInstance* inst, const double* xstar, CPXCALLBACKCONTEXTptr context);
int mylazy_separation(tspInstance* inst, const double* xstar, CPXCALLBACKCONTEXTptr context);

void parse_command_line(int argc, char** argv, tspInstance* inst);
void free_instance(tspInstance* inst);
void read_input(tspInstance* inst, result* res);
void print_error(const char* err);
double dist(int i, int j, tspInstance* inst);
void get_pipe(tspInstance* inst, result* res);
void plot_problem_input(tspInstance* inst, result* res, FILE* gnuplot);
void plot_points(FILE* gnuplot, char* pngname, tspInstance* inst);
void plot_edges(FILE* gnuplot, char* pngname, tspInstance* inst, result* res);
void setup_style(FILE* gnuplot, tspInstance* inst);
void save_results(tspInstance* inst, result* res, char* f_name);
#endif   /* TSP_H_ */ 
