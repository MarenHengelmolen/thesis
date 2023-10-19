#include <iostream>
#include "definitions.h"
#include "input.h"
#include "output.h"
#include "buildings.h"
//#include "cfd.h"
#include "geometry.h"
#include "meshing.h"

using namespace std::chrono;

int main(int argc, const char * argv[]) {

    const string filename = (argc > 1) ? argv[1] : "../data/Mesh.obj";
    //const string filename = (argc > 1) ? argv[1] : "../data/Mesh_without_terrain.obj";
    //const string filename = (argc > 1) ? argv[1] : "../data/TUDelft_EWISurroundingsStandard.stl";
    //const string filename = (argc > 1) ? argv[1] : "../data/cube.obj";
    //const string filename = (argc > 1) ? argv[1] : "../data/rectangle.obj";
    //const string filename = (argc > 1) ? argv[1] : "../data/TUDelft.stl";
    //const string filename = (argc > 1) ? argv[1] : "../data/Binnenstad.obj";
    //const string filename = (argc > 1) ? argv[1] : "../data/Buitenwatersloot.obj";
    //const string filename = (argc > 1) ? argv[1] : "../data/Begijnhofbuurt.obj";
    //const string filename = (argc > 1) ? argv[1] : "../data/TUDelft_EWISurroundingsStandard_Rotated.stl";
    //const string filename = (argc > 1) ? argv[1] : "../data/trapeze.obj";

    //BASICS
    string input_file;
    input_file = filename;

    clock_t start, end;
    double time_taken;
    start = clock();

    double wd = 0;

    int idx = input_file.find_last_of(".");
    if (idx != string::npos){
        string format = input_file.substr(idx+1);

        if (format == "obj") {
            cout << "This is an OBJ file!" << endl;
            OBJ_data(input_file, wd);
            //OBJ_data(input_file);

            save2obj_terrain("../output/terrain.obj");
            cout << "Terrain successfully saved" << endl;

//            save2obj_iBuildings("../output/buildings.obj");
//            cout << "Buildings successfully saved" << endl;

        }

        else if (format == "stl") {
            cout << "This is an STL file!" << endl;
            STL_data(input_file, wd);
        }
    }

    //BUILDING DEFINITION
    start = clock();
    define_buildings();
    end = clock();

    time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by COMMON VERTICES is : " << fixed << time_taken << setprecision(5);
    cout << " sec " << endl;

//    //COMPUTATIONAL DOMAIN
//    start = clock();
//    computational_domain("../output/computational_domain.obj", "../output/computational_domain.txt", 20);
//    end = clock();
//    time_taken = double(end - start) / double(CLOCKS_PER_SEC);
//    cout << "Time taken by COMPUTATIONAL DOMAIN is : " << fixed << time_taken << setprecision(5);
//    cout << " sec " << endl;

//    //REFINEMENT BOXES
//    start = clock();
//    RefinementBoxes(5);
//    end = clock();
//    time_taken = double(end - start) / double(CLOCKS_PER_SEC);
//    cout << "Time taken by REFINEMENT BOXES is : " << fixed << time_taken << setprecision(5);
//    cout << " sec " << endl;
//
    //SHORT EDGES
    start = clock();
    double th_edges = 3;
    shortedges("../output/short_edges.obj", "../output/short_edges.txt", th_edges);
    end = clock();
    time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by SHORT EDGE IDENTIFICATION is : " << fixed << time_taken << setprecision(5);
    cout << " sec " << endl;

    //SHARP ANGLES
    start = clock();
    double th_angles = 40;
    sharp_angles("../output/sharp_angles.obj", "../output/sharp_angles.txt", th_angles);
    end = clock();
    time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by SHARP ANGLE IDENTIFICATION is : " << fixed << time_taken << setprecision(5);
    cout << " sec " << endl;

    //SLIVER TRIANGLES
    start = clock();
    double th_slivers = 0.1;
    sliver_triangles("../output/sliver_triangles.obj", "../output/sliver_triangles.txt", th_slivers);
    end = clock();
    time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by SLIVER TRIANGLES IDENTIFICATION is : " << fixed << time_taken << setprecision(5);
    cout << " sec " << endl;

//    //SHORT DISTANCE BETWEEN BUILDINGS IDENTIFICATION
//    start = clock();
//    double th_distance = 3;
//    distances_between_buildings("../output/short_distances.obj", "../output/short_distances.txt", th_distance);
//    end = clock();
//    time_taken = double(end - start) / double(CLOCKS_PER_SEC);
//    cout << "Time taken by DISTANCE BETWEEN BUILDINGS is : " << fixed << time_taken << setprecision(5);
//    cout << " sec " << endl;

    define_ground_surfaces();
    save2obj_gs("../output/gs.obj");

    start = clock();
    double threshold = 0.5;
    check_topological_relationship(0.5);
    topological_relationships("../output/topo.obj", "../output/topo.txt",  threshold);

    end = clock();
    time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by TOPOLOGICAL RELATIONSHIP CHECK is : " << fixed << time_taken << setprecision(5);
    cout << " sec " << endl;



//    start = clock();
//    //Liu_region_of_interest(1, 0, 0);
//    region_of_interest("../output/RoI.obj", "../output/RoI.txt", 1, 5, 5);
//    end = clock();
//    time_taken = double(end - start) / double(CLOCKS_PER_SEC);
//    cout << "Time taken by ROI is : " << fixed << time_taken << setprecision(5);
//    cout << " sec " << endl;


//    start = clock();
//    region_of_interest("../output/detailed_buildings.obj", "../output/detailed_buildings.txt", 1, 0, 0);
//    BuildingVolume(1, 0, 0);
//    end = clock();
//    time_taken = double(end - start) / double(CLOCKS_PER_SEC);
//    cout << "Time taken by CYLINDER is : " << fixed << time_taken << setprecision(5);
//    cout << " sec " << endl;
//
//    start = clock();
//    double thresholdBS = 2;
//    BuildingSeparation(1, 0, 0, thresholdBS);
//    end = clock();
//    time_taken = double(end - start) / double(CLOCKS_PER_SEC);
//    cout << "Time taken by BUILDING SEPARATION is : " << fixed << time_taken << setprecision(5);
//    cout << " sec " << endl;
//
//    start = clock();
//    Meshing3(0, 0, 0, 2);
//    end = clock();
//    time_taken = double(end - start) / double(CLOCKS_PER_SEC);
//    cout << "Time taken by MESHING is : " << fixed << time_taken << setprecision(5);
//    cout << " sec " << endl;

//    save2obj_building("../output/zerovolume.obj", 146);
    start = clock();
    meshing();
    end = clock();
    time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by MESHING is : " << fixed << time_taken << setprecision(5);
    cout << " sec " << endl;

    write_blockMeshDict();
    double eh_target = 1.5;
    int idx2 = filename.find_last_of("/");
    string name = input_file.substr(idx2+1);
    write_snappyHexMeshDict(name, eh_target);

    cout << "\nfilename " << filename << endl;
    cout << "wd " << wd << endl;

    cout << "n buildings " << buildings.size() << endl;

//    start = clock();
//    computational_domain("../output/computational_domain.obj", "../output/computational_domain.txt", 20);
//    end = clock();
//    time_taken = double(end - start) / double(CLOCKS_PER_SEC);
//    cout << "Time taken by COMPUTATIONAL DOMAIN is : " << fixed << time_taken << setprecision(5);
//    cout << " sec " << endl;

    //cout << "building size " << buildings.size() << endl;
//    for (unordered_map<int, Building>::iterator bb=buildings.begin(); bb!=buildings.end(); bb++) {
//        string name = str(boost::format{"../output/building_%s.obj"} % bb->first);
//        save2obj_building(name, bb->first);
//    }

//    for (unordered_map<int, Building>::iterator bb=buildings.begin(); bb!=buildings.end(); bb++) {
//
//        double vol = CGAL::Polygon_mesh_processing::volume(bb->second.mesh);
//        cout << "id " << bb->first << " vol " << vol << endl;
//    }

//    buildingseparations("../output/buildingseparations1.obj", "../output/buildingseparations2.obj");
//    building_volume("../output/buildingvolume.obj");

    return 0;
}
