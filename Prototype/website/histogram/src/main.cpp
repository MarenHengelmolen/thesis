#include <iostream>
#include "definitions.h"
#include "input.h"
#include "buildings.h"
#include "output.h"

#include <string>

using namespace std::chrono;

int main(int argc, const char * argv[]) {

    const string file = argv[1];
    const string filename(file);

    const string d_file = argv[2];
    const string file_with_distances(d_file);

    const string b = argv[3];
    const string buildings(b);

//    const string filename = (argc > 1) ? argv[1] : "../data/Maastoren.obj";
//    string buildings = "../output/buildings.obj";
//    string file_with_distances = "../output/distances.txt";

    //Input
    string input_file;
    input_file = filename;

    double wd = 90;
    double unit = 1;
    double h_gs = 3;
    double theta_gs = 45;

    int idx = input_file.find_last_of(".");

    if (idx != string::npos){
        string format = input_file.substr(idx+1);

        if (format == "obj") {
            OBJ_data(input_file, wd, unit);
        }

        else if (format == "stl") {
            STL_data(input_file, wd, unit);
        }
    }

    define_buildings(buildings);
    define_ground_surfaces(h_gs, theta_gs);
    distances_for_histogram(file_with_distances);

    return 0;
}
