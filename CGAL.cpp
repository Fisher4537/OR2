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

int verbose = 0;
int nnodes = 0;

vector<K::Point_2> points;
vector<K::Point_2> points_ordered;
vector<int> indexes;
vector<double> sqrt_distance;

vector<vector<int>> sol;

int status;                         // 0 : all OK
                                    // 1 : error on loading points
                                    // 2 : error, points are not loaded into vector

extern "C" {
    
    void set_verbose(int v) {
        verbose = v;
    }

    int load_point(char* pathFileTSP) {

        points.clear();
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

                if (verbose > 100)
                    printf("Id:\t%.0f,\tCoordinate:\t[%.0f, %.0f]\n", id, x, y);
                nnodes++;

                points.push_back(K::Point_2(x, y));


            }
            if (line == "NODE_COORD_SECTION")
                startPoint = true;
        }
        return status = 0;
    }

    int order_by_dis(int firstPoint, int with_sqrt_distance) {

        if(verbose > 100)
            printf("\n____ Euclidean Distances Solutions: %d ____\n", firstPoint);

        points_ordered.clear();
        indexes.clear();

        K::Point_2 first_point(points[firstPoint].x(), points[firstPoint].y());
        //points.erase(points.begin() + firstPoint);
            
        // Search N points every iter
        const unsigned int N = 1;

        Tree tree(points.begin(), points.end());

        // Initialize the search structure, and search all points
        Neighbor_search search(tree, first_point, points.size());
        
        

        // report the N nearest neighbors and their distance
        // This should sort all N points by increasing distance from origin
        Neighbor_search::iterator it = search.begin();
        points_ordered.push_back(it->first);
        indexes.push_back(firstPoint);
        if (verbose > 100)
            printf("Point %d [%.0f,%.0f]\tat distance: %.4f\n", 0, it->first.x(), it->first.y(), 0);
        ++it;
        
        for (; it != search.end(); ++it) {
            
            points_ordered.push_back(it->first);
            
            int id = distance(points.begin(), find(points.begin(), points.end(), it->first));

            bool duplicate = true;
            while (duplicate) {
                if (find(indexes.begin(), indexes.end(), id) != indexes.end()) {
                    id = distance(points.begin(), find(points.begin() + id + 1, points.end(), it->first));
                }
                else
                    duplicate = false;
            }
                
            
            indexes.push_back(id);
            if (with_sqrt_distance)
                sqrt_distance.push_back(sqrt(it->second));

            if(verbose > 100)
                //printf("Point %d [%.0f,%.0f]\tat distance: %.4f\n", id, it->first.x(), it->first.y(), sqrt(it->second));
                printf("%d\t", id);
        }
        return 0;
    }

    int greedy_alg() {

        vector<vector<int>> neigh_sol_idx(nnodes);
        
        sol = vector<vector<int>>(nnodes);

        if (verbose > 100)
            printf("\n____ Track Solutions: ____\n");
        
        for (int i = 0; i < nnodes; i++) {

            order_by_dis(i, 1);
            
            if (verbose > 100)
                printf("\n");
            
            neigh_sol_idx[i] = indexes;
            sol[i] = vector<int>(nnodes);
            fill(sol[i].begin(), sol[i].end(), -1);
        }
        if (verbose > 100)
            printf("\n");
        
        for (int i = 0; i < nnodes; i++) {
            if (verbose > 100)
                printf("\n");
            
            sol[i][0] = neigh_sol_idx[i][0];                // first_point ( euclidean distance = 0 from itself)
            sol[i][1] = neigh_sol_idx[i][1];
            int succ = neigh_sol_idx[i][1];
            
            if (verbose > 100)
                printf("%d\t%d\t", sol[i][0], sol[i][1]);

            for (int j = 2; j < nnodes; j++) {
                
                vector<int>::iterator it = find_if(neigh_sol_idx[succ].begin(), neigh_sol_idx[succ].end(), [&](int neig) {
                    return find(sol[i].begin(), sol[i].end(), neig) == sol[i].end();
                });

                succ = *it;
                sol[i][j] = succ;

                if (verbose > 100)
                    printf("%d\t", succ);
            }
        }
        return 0;
    }

    int* get_greedy_sol(int i) {
        return sol[i].data();
    }
}