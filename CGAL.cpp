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

vector<K::Point_2> points;
vector<K::Point_2> points_ordered;
vector<int> indexes;
vector<double> sqrt_distance;

int status;                         // 0 : all OK
                                    // 1 : error on loading points
                                    // 2 : error, points are not loaded into vector

extern "C" {
    
    void set_verbose(int v) {
        verbose = v;
    }

    int load_point(char* pathFileTSP)   {

        points.clear();
        ifstream tspFile;
        tspFile.open(pathFileTSP);

        string line;
        bool startPoint = false;
        while (getline(tspFile, line)){
            if (line == "EOF")
                break;
            if (startPoint) {
                istringstream iss(line);
                double id, x, y;
                if (!(iss >> id >> x >> y)) { 
                    printf("Error pushing Point number: %d\n", (int)id);
                    return status = 1;
                }
                points.push_back(K::Point_2(x, y));
                if(verbose > 100)
                    printf("Coordinate: [%f, %f]\n", x, y);
            }
            if (line == "NODE_COORD_SECTION")
                startPoint = true;
        }
        return status = 0;
    }

    int order_by_dis(int firstPoint, int with_sqrt_distance) {

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
        ++it;                                                       //don't consider first_point
        for (int i = 0; it != search.end() && i < 2; ++it, i++) {
            
            //points_ordered.push_back(it->first);
            int id = distance(points.begin(), find(points.begin(), points.end(), it->first));
            indexes.push_back(id);
            if (with_sqrt_distance)
                sqrt_distance.push_back(sqrt(it->second));

            if(verbose > 100)
                printf("Point %d [%.0f,%.0f]\tat distance: %.4f\n", id, it->first.x(), it->first.y(), sqrt(it->second));

        }
        return 0;
    }

    int get_first() {
        return indexes.at(0);
    }

    int get_second() {
        return indexes.at(1);
    }
}