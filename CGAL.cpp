#include <sstream>
#include <string>
#include <fstream>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/Cartesian.h>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Cartesian<double> K;

using namespace std;

vector<K::Point_2> points;
double f_x, f_y;

int status;                         // 0 : all OK
                                    // 1 : error on loading points
                                    // 2 : error, points are not loaded into vector

extern "C" {
    
    int load_point(char* pathFileTSP){

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
                printf("Coordinate: [%f, %f]\n", x, y);
                f_x = x;
                f_y = y;
            }
            if (line == "NODE_COORD_SECTION")
                startPoint = true;
        }
        return status = 0;
    }

    int order_by_dis(int firstPoint) {
        if (status) {

           K::Point_2 firstPoint(f_x, f_y); //points.at(firstPoint).x, points.at(firstPoint).y);

            sort(points.begin(), points.end(), [&firstPoint](const K::Point_2& point1, const K::Point_2& point2 ) {
                return CGAL::squared_distance(firstPoint, point1) < CGAL::squared_distance(firstPoint, point2) ?
                    point1 : point2;
                });

            return 0;
        }
        else return status = 2;
    }
}