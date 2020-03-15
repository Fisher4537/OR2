#ifndef VRP_H_  

#define VRP_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>  

#include <cplex.h>  
#include <windows.h>  

#define VERBOSE				    50		// printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)

//data structures  
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
	
} instance;

//NUOVA STRUTTURA DATI PER I RISULTATI
typedef struct {
	//global data
	double	tstart;
	double zbest;							// best sol. available  
	double tbest;							// time for the best sol. available  
	double *best_sol;						// best sol. available    
	double	best_lb;						// best lower bound available  
	double *load_min;						// minimum load when leaving a node
	double *load_max;						// maximum load when leaving a node
} result;

//inline
inline int imax(int i1, int i2) { return (i1 > i2) ? i1 : i2; }
inline double dmin(double d1, double d2) { return (d1 < d2) ? d1 : d2; }
inline double dmax(double d1, double d2) { return (d1 > d2) ? d1 : d2; }

#endif   /* VRP_H_ */ 
