#ifndef TSP_H_

#define TSP_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <cplex.h>

#ifdef _WIN32
	//TODO : use pthread on Win?
#else
	#include <pthread.h>
#endif

//hard-wired parameters
#define XSMALL		  		  1e-5 		// 1e-4*	// tolerance used to decide // ingerality of 0-1 var.s
#define EPS 1e-5
#define TICKS_PER_SECOND 	  1000.0  	// cplex's ticks on Intel Core i7 quadcore // @2.3GHZ

//data structures
typedef struct {
	double max_fr;		// maximum fixing_ratio
	double incr_fr;		// increase of fixing_ratio when good gap
	double decr_fr;		// decreasing of fixing_ratio when bad gap
	double good_gap;			// a good gap allow to decrease fixing_ratio
	double optimal_gap; 	// under this value, solution is optimal
} hardfix;

typedef struct arches_tuple {
	int arc;
	struct arches_tuple* next;
} tabu_list;

typedef struct tour {
	tabu_list* entry;
	struct tour* next;		//struct to store an array of Linked lists
} tour_list;


typedef struct {

	//input data
	int nnodes;
	double *xcoord;
	double *ycoord;
	int verbose;

	// parameters
	int setup_model;
	int model_type;					// TSP
	int randomseed;
	int nthread;					// number of threads
	double init_time;			// initial execution time
	double timelimit;				// overall time limit, in sec.s
	char input_file[100];			// input file
	int edge_weight_type;
	int available_memory;
	int max_nodes; 					// max n. of branching nodes in the final run
									// (-1 unlimited)
	// double cutoff; 				// cutoff (upper bound) for master
	int integer_costs;

	//global data
	double opt_time;				// optimization time
	double zbest;					// best sol. available
	double tbest;					// time for the best sol. available
	double *best_sol;				// best sol. available
	double	best_lb;				// best lower bound available
	int nedges;						// number of edges in best_sol
	// double *load_min;			// minimum load when leaving a node
	// double *load_max;			// maximum load when leaving a node
	double best_int;

	double ncols;
	int callback;
	int mip_opt;
	int build_sol;
	int warm_start;
	int useCplex;
	int heuristic;
	int plot_style;
	int plot_edge;

	// hard fixing
	double max_fr;		// maximum fixing_ratio
	double incr_fr;		// increase of fixing_ratio when good gap
	double decr_fr;		// decreasing of fixing_ratio when bad gap
	double good_gap;			// a good gap allow to decrease fixing_ratio
	double optimal_gap; 	// under this value, solution is optimal

} tspinstance;


//inline
inline int imax(int i1, int i2) { return ( i1 > i2 ) ? i1 : i2; }
inline double dmin(double d1, double d2) { return ( d1 < d2 ) ? d1 : d2; }
inline double dmax(double d1, double d2) { return ( d1 > d2 ) ? d1 : d2; }


char* model_name(int i);
char * setup_model(tspinstance* inst);
int TSPopt(tspinstance *inst);
int xpos(int i, int j, tspinstance *inst);
int* invers_xpos(int pos, tspinstance* inst);
int asym_xpos(int i, int j, tspinstance *inst);
int asym_upos(int i, tspinstance *inst);
int asym_ypos(int i, int j, tspinstance* inst);

void build_model(tspinstance *inst, CPXENVptr env, CPXLPptr lp); 		// interface
void build_sym_std(tspinstance *inst, CPXENVptr env, CPXLPptr lp); 		// sym, std
void build_mtz(tspinstance *inst, CPXENVptr env, CPXLPptr lp); 			// asym, MTZ
void build_flow1(tspinstance *inst, CPXENVptr env, CPXLPptr lp);
void build_mtz_lazy(tspinstance* inst, CPXENVptr env, CPXLPptr lp);
void add_lazy_mtz(tspinstance* inst, CPXENVptr env, CPXLPptr lp);

void switch_warm_start(tspinstance* inst, CPXENVptr env, CPXLPptr lp, int* status);

void optimization(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status);
int local_branching(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status);
int hard_fixing(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status);
void fix_bound(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status, double fixing_ratio);

int tabu_search(CPXENVptr env, tspinstance* inst, int* status);
void push(tabu_list** head, int arc, int isArc);
int pop_first(tabu_list** head);
int pop_last(tabu_list* head);
int remove_by_index(tabu_list** head, int n);
int contained_in_posix(tabu_list** head, int arc);
void print_list(tabu_list* head);
void write_list_lb(tabu_list* head);
void delete_list(tabu_list* head, tabu_list** head_ref);

int tabu_search_array(CPXENVptr env, tspinstance* inst, int* status);
int contained_in_posix_array(int tabu_array_size, int* tabu_array, int arc);

int simulating_annealing(CPXENVptr env, tspinstance* inst, int* status);
int max_dist_couple_nodes(tspinstance* inst);

int genetic_algorithm(CPXENVptr env, tspinstance* inst, int* status);
void init_population(tspinstance* inst, double** population, int nPop);
void init_frequency_edges(tspinstance* inst, double** population, int* frequencyTable, int nPop);
void shuffle_individuals(tspinstance* inst, double** population, int nPop);
void swap(double* a, double* b);
int EAX_Single(tspinstance* inst, double** population, double** kids, int pA, int pB, int nKids);
void extract_ABcycles(tspinstance* inst, double** population, int pA, int pB, double** ABcycles, double* graph_AB, int* idxCycle, int maxNcycles, int** edges_cycles_EA);
int build_sol_ga(tspinstance* inst, const double* xstar, int* succ, int* prev, int* comp, int* ncomp);
void evaluate_traced_ABcycle(tspinstance* inst, double* traced_AB, double** ABcycles, int* idxCycle, int* tourFound, int* edges_cycles_EA_current);

tour_list* grapth_to_tree(tspinstance* inst, int* nodes_one, int* nodes_two, tour_list* tours, int* edges_cycles_EA_current);
tour_list* Tree_recursive(tspinstance* inst, int current, int* nodes_one, int* nodes_two, int* found, tabu_list* pathlist, tabu_list* visited_nodes, tour_list* tours, int* edges_cycles_EA_current);
void copy_in_i_j(tspinstance* inst, int* nodes_one, int* nodes_two, int* i, int* j);
void push_list_on_list(tspinstance* inst, tour_list** head_ref, tour_list** pathlist, int pos, int* edges_cycles_EA_current);
int print_list_of_list(tour_list* tours);
tabu_list* copy(tabu_list* org);
int patching_two_edges(tspinstance* inst, double* tour);

void survival_selection(tspinstance* inst, double** population, int nPop, int* frequencyTable, int nKids, int pA, double** kids);
void update_frequency_table(tspinstance* inst, int* frequencyTable, double* pA, double* kid);
double calc_L(tspinstance* inst, double** population, int nPop);
double calc_H(tspinstance* inst, int* frequencyTable, int nPop);

void print_population(tspinstance* inst, double** population, int nPop);
void plot_population(tspinstance* inst, double** population, int nPop);
void plot_single(tspinstance* inst, double* individual);
void print_frequency_table(tspinstance* inst, int* frequency_table);
void free_ga(double** population, int* frequency_table, int nPop);

int* heur_greedy_cgal(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status);
int* heur_greedy(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status);
int succ_not_contained(int node, int* sol, tspinstance* inst);
int* heur_grasp(tspinstance* inst, int* status);
int* heur_insertion(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status);
int insertion_move(tspinstance* inst, int* best_sol, int count_sol, int vertex);
int contained_in_index(int* vector, int count_sol, int elem);

void two_opt(tspinstance *inst);
int best_two_opt(tspinstance* inst);
void random_two_opt(tspinstance* inst);
void random_n_opt(tspinstance* inst, int n);

void patching(tspinstance* inst);
void single_patch(tspinstance* inst, int* succ, int* comp, int* ncomp);

int is_clockwise(tspinstance* inst, int x1, int x2, int x3);

void clear_sol(tspinstance* inst);

int mip_optimization(CPXENVptr env, CPXLPptr lp, tspinstance *inst, int *status);
int subtour_iter_opt(CPXENVptr env, CPXLPptr lp, tspinstance *inst, int *status);
int subtour_heur_iter_opt(CPXENVptr env, CPXLPptr lp, tspinstance* inst, int* status, int heuristic);

void switch_callback(tspinstance* inst, CPXENVptr env, CPXLPptr lp);
int CPXPUBLIC lazycallback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p);
int CPXPUBLIC genericcallback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* cbhandle);
int mylazy_separation(tspinstance* inst, const double* xstar, CPXCENVptr env, void* cbdata, int wherefrom);
int mygeneric_separation(tspinstance* inst, const double* xstar, CPXCALLBACKCONTEXTptr context);

void build_sol(tspinstance *inst, int *succ, int *comp, int *ncomp);
void build_sol_sym(tspinstance *inst, int *succ, int *comp, int *ncomp);
void build_sol_lazy_std(tspinstance* inst, const double* xstar, int* succ, int* comp, int* ncomp);
void build_sol_mtz(tspinstance *inst, int *succ, int *comp, int *ncomp);
void build_sol_flow1(tspinstance *inst, int *succ, int *comp, int *ncomp);

void parse_command_line(int argc, char** argv, tspinstance *inst);
void read_input(tspinstance *inst);
void free_instance(tspinstance *inst);

double dist(int i, int j, tspinstance *inst);		// get distance between two nodes

int save_results(tspinstance *inst, char *f_name);	// Save model performance

void plot_instance(tspinstance *inst);
char * get_file_name(char *path, char *name);
void setup_style(FILE *gnuplot, tspinstance *inst);
void setup_linestyle1(FILE *gnuplot);
void setup_linestyle2(FILE *gnuplot);
void setup_arrowstyle2(FILE *gnuplot);
void plot_points(FILE *gnuplot, char *pngname, tspinstance *inst);
void plot_edges(FILE *gnuplot, char *pngname, tspinstance *inst);
void plot_lines_sym(FILE *gnuplot, char *pngname, tspinstance *inst);
void plot_arrow_asym(FILE *gnuplot, char *pngname, tspinstance *inst);

// Debug functions
void print_succ(int* succ, tspinstance* inst);
void pause_execution();
void print_error(const char *err);

//*********************************** CGAL Methods ***********************************
void set_verbose(int v);
int load_point(char* pathFileTSP);
int order_by_dis(int firstPoint, int with_sqrt_distance);
int greedy_alg();
int* get_greedy_sol(int i);
void free_cgal();

#endif   /* TSP_H_ */
