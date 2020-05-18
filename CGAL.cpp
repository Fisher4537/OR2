#include <sstream>
#include <string>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_2.h>
#include <cmath>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Search_traits_2<K> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;

using namespace std;


int nnodes = 0;
int verbose = 0;                    // 100 some details
                                    // 1000 all the accurate details

class varCGAL {

    public :
        vector<K::Point_2> points;
        vector<K::Point_2> points_ordered;
        vector<int> indexes;
        vector<double> sqrt_distance;

        vector<vector<int>> sol;

};

varCGAL* var = new varCGAL();


int status;                         // 0 : all OK
                                    // 1 : error on loading points
                                    // 2 : error, points are not loaded into vector

extern "C" {
    
    void set_verbose(int v) {
        verbose = v;
    }

    int load_point(char* pathFileTSP) {

        var->points.clear();
        ifstream tspFile;
        tspFile.open(pathFileTSP);

        string line;
        bool startPoint = false;
        while (getline(tspFile, line)) {
            if (line == "EOF")
                break;
            if (startPoint) {
                istringstream iss(line);
                double id, x, y;

                if (!(iss >> id >> x >> y)) {
                    printf("Error pushing Point number: %d\n", (int)id);
                    return status = 1;
                }

                if (verbose > 1000)
                    printf("Id:\t%.0f,\tCoordinate:\t[%.0f, %.0f]\n", id, x, y);
                nnodes++;

                var->points.push_back(K::Point_2(x, y));


            }
            if (line == "NODE_COORD_SECTION")
                startPoint = true;
        }
        return status = 0;
    }

    int order_by_dis(int firstPoint, int with_sqrt_distance) {

        if(verbose > 1000)
            printf("\n____ Euclidean Distances Solutions: %d ____\n", firstPoint);

        var->points_ordered.clear();
        var->indexes.clear();

        K::Point_2 first_point(var->points[firstPoint].x(), var->points[firstPoint].y());
        //points.erase(points.begin() + firstPoint);
            
        // Search N points every iter
        const unsigned int N = 1;

        Tree tree(var->points.begin(), var->points.end());

        // Initialize the search structure, and search all points
        Neighbor_search search(tree, first_point, var->points.size());
        
        

        // report the N nearest neighbors and their distance
        // This should sort all N points by increasing distance from origin
        Neighbor_search::iterator it = search.begin();
        var->points_ordered.push_back(it->first);
        var->indexes.push_back(firstPoint);
        if (verbose > 1000)
            printf("%d\t", 0);
            //printf("Point %d [%.0f,%.0f]\tat distance: %.4f\n", 0, it->first.x(), it->first.y(), 0.0);
        ++it;
        
        for (; it != search.end(); ++it) {
            
            var->points_ordered.push_back(it->first);
            
            int id = distance(var->points.begin(), find(var->points.begin(), var->points.end(), it->first));

            bool duplicate = true;
            while (duplicate) {
                if (find(var->indexes.begin(), var->indexes.end(), id) != var->indexes.end()) {
                    id = distance(var->points.begin(), find(var->points.begin() + id + 1, var->points.end(), it->first));
                }
                else
                    duplicate = false;
            }
                
            
            var->indexes.push_back(id);
            if (with_sqrt_distance)
                var->sqrt_distance.push_back(sqrt(it->second));

            if(verbose > 1000)
                //printf("Point %d [%.0f,%.0f]\tat distance: %.4f\n", id, it->first.x(), it->first.y(), sqrt(it->second));
                printf("%d\t", id);
        }
        return 0;
    }

    int greedy_alg() {

        vector<vector<int>> neigh_sol_idx(nnodes);
        
        var->sol = vector<vector<int>>(nnodes);
        
        for (int i = 0; i < nnodes; i++) {

            order_by_dis(i, 1);
            
            if (verbose > 1000)
                printf("\n");
            
            neigh_sol_idx[i] = var->indexes;
            var->sol[i] = vector<int>(nnodes);
            fill(var->sol[i].begin(), var->sol[i].end(), -1);
        }
        if (verbose > 1000)
            printf("\n____ Track Solutions: ____\n");

        for (int i = 0; i < nnodes; i++) {
            if (verbose > 1000)
                printf("\n____ Solutions: %d ____\n", i);
            
            var->sol[i][0] = neigh_sol_idx[i][0];                // first_point ( euclidean distance = 0 from itself)
            var->sol[i][1] = neigh_sol_idx[i][1];
            int succ = neigh_sol_idx[i][1];
            
            if (verbose > 1000)
                printf("%d\t%d\t", var->sol[i][0], var->sol[i][1]);

            for (int j = 2; j < nnodes; j++) {
                
                vector<int>::iterator it = find_if(neigh_sol_idx[succ].begin(), neigh_sol_idx[succ].end(), [&](int neig) {
                    return find(var->sol[i].begin(), var->sol[i].end(), neig) == var->sol[i].end();
                });

                succ = *it;
                var->sol[i][j] = succ;

                if (verbose > 1000)
                    printf("%d\t", succ);
            }
        }
        return 0;
    }

    int* get_greedy_sol(int i) {
        return var->sol[i].data();
    }

    void free_cgal() {
        
        size_t mem = 0;
        var->points.clear(); mem += var->points.capacity();
        var->points_ordered.clear(); mem += var->points_ordered.capacity();
        var->indexes.clear(); mem += var->indexes.capacity();
        var->sqrt_distance.clear(); mem += var->sqrt_distance.capacity();
        
        var->points.shrink_to_fit(); mem -= var->points.capacity();
        var->points_ordered.shrink_to_fit(); mem -= var->points_ordered.capacity();
        var->indexes.shrink_to_fit(); mem -= var->indexes.capacity();
        var->sqrt_distance.shrink_to_fit(); mem -= var->sqrt_distance.capacity();
        
        if (verbose > 100)
            printf("CGAL Memory freed: %zd", mem);
    }
}