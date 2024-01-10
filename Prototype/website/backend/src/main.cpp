#include <iostream>
#include "definitions.h"
#include "input.h"
#include "output.h"
#include "buildings.h"
#include "geometry.h"
#include "mesh3.h"
#include "mesh2.h"
#include "parameters.h"

#include <string>

using namespace std::chrono;

int main(int argc, const char * argv[]) {

    //I. INPUT FILE

    const string file = argv[1];
    const string filename(file);

    const string up = argv[2];
    const string user_parameters(up);

//    const string filename = (argc > 1) ? argv[1] : "../data/Binnenstad.obj"; // needed when the code is directly executed from backend
    string input_file;
    input_file = filename;

    //II. PARAMETERS

    //II.a. Read input parameters from website

    parameters(user_parameters);
//    parameters("../data/user_parameters.txt"); // needed when the code is directly executed from backend

    //CFD parameters
    string save_buildings = user.buildings; //path to OBJ file including buildings of the input model
    double fd = user.fd; //flow direction
    double unit = user.unit; //unit parameter
    int n_boxes = user.n_boxes; //number of refinement boxes
    double h_user = user.h_user; //target evaluation height
    double z0 = user.z0; //maximum roughness height
    double N_max = user.N_max; //limitation of the number of cells
    double r_min = user.r_min; //required cell ratio
    double d_separation = user.d_separation; //building separation threshold
    int BR = user.BR; //indicates which blockage ratios are used
    int dim = user.dim; //indicates whether hmax (0) or dmax (1) must be used for domain
    int target = user.target; //indicates whether there is a target building
    double x = user.x; //x-coordinate target building
    double y = user.y; //y-coordinate target building
    vector<double> f1 = user.f1; //dimensions refinement box 1
    vector<double> f2 = user.f2; //dimensions refinement box 2
    vector<double> f3 = user.f3; //dimensions refinement box 3

    //Geometric parameters
    double snap_tolerance = user.snap_tolerance; //snap tolerance
    double planarity_tolerance = user.planarity_tolerance; //planarity tolerance
    double overlap_tolerance = user.overlap_tolerance; //overlap tolerance
    double th_topo = user.th_topo; //threshold topological relationships validations
    int ground_level = user.ground_level; //indicates whether a ground level is inserted
    double z_ground = user.z_ground; //z-value of this ground level
    double th_slivers = user.th_slivers; //threshold sliver triangles
    double th_angles = user.th_angles; //threshold sharp angles
    double th_edges = user.th_edges; //threshold short edges
    double h_gs = user.h_gs; //threshold height ground surfaces
    double theta_gs = user.theta_gs; //threshold angle ground surfaces

    //Path output files
    string name = user.name;
    string save_gs = user.save_gs;
    string topocheckOBJ = user.topocheckOBJ;
    string topocheckTXT = user.topocheckTXT;
    string sliversOBJ = user.sliversOBJ;
    string sliversTXT = user.sliversTXT;
    string sharpanglesOBJ = user.sharpanglesOBJ;
    string sharpanglesTXT = user.sharpanglesTXT;
    string shortedgesOBJ = user.shortedgesOBJ;
    string shortedgesTXT = user.shortedgesTXT;
    string overlappingbuildingsOBJ = user.overlappingbuildingsOBJ;
    string overlappingbuildingsTXT = user.overlappingbuildingsTXT;
    string val3dity_report = user.val3dity_report;
    string rotated_modelOBJ = user.rotated_modelOBJ;
    string iblockMeshDict = user.iblockMeshDict;
    string blockMeshDict = user.blockMeshDict;
    string blockMeshDict_noNmax = user.blockMeshDict_noNmax;
    string blockMeshDict_Nmax = user.blockMeshDict_Nmax;
    string isnappyHexMeshDict = user.isnappyHexMeshDict;
    string snappyHexMeshDict = user.snappyHexMeshDict;
    string snappyHexMeshDict_noNmax = user.snappyHexMeshDict_noNmax;
    string snappyHexMeshDict_Nmax = user.snappyHexMeshDict_Nmax;
    string domainOBJ = user.domainOBJ;
    string domainOBJ_noNmax = user.domainOBJ_noNmax;
    string domainOBJ_Nmax = user.domainOBJ_Nmax;
    string b_volumeOBJ = user.b_volumeOBJ;
    string b_volumeTXT = user.b_volumeTXT;
    string b_separations_buildings = user.b_separations_buildings;
    string b_separations_edges = user.b_separations_edges;
    string b_separationsTXT = user.b_separationsTXT;

    string b_volumeOBJ_noNmax = user.b_volumeOBJ_noNmax;
    string b_volumeTXT_noNmax = user.b_volumeTXT_noNmax;
    string b_separations_buildings_noNmax = user.b_separations_buildings_noNmax;
    string b_separations_edges_noNmax = user.b_separations_edges_noNmax;
    string b_separationsTXT_noNmax = user.b_separationsTXT_noNmax;

    string b_volumeOBJ_Nmax = user.b_volumeOBJ_Nmax;
    string b_volumeTXT_Nmax = user.b_volumeTXT_Nmax;
    string b_separations_buildings_Nmax = user.b_separations_buildings_Nmax;
    string b_separations_edges_Nmax = user.b_separations_edges_Nmax;
    string b_separationsTXT_Nmax = user.b_separationsTXT_Nmax;

    string RoIOBJ = user.RoIOBJ;
    string RoITXT = user.RoITXT;
    string RoICylinder = user.RoICylinder;
    string output_for_UI_file = user.output_for_UI;

    string file_with_distances = user.distances;

    //III. Read input file
    int idx = input_file.find_last_of(".");
    if (idx != string::npos){
        string format = input_file.substr(idx+1);

        if (format == "obj") {
            OBJ_data(input_file, fd, unit); //store vertices and faces
            run_val3dity(input_file, val3dity_report, snap_tolerance, planarity_tolerance, overlap_tolerance);
            store_rotated_model(rotated_modelOBJ);
        }

        else if (format == "stl") {
            STL_data(input_file, fd, unit); //store vertices and faces
            run_val3dity(input_file, val3dity_report, snap_tolerance, planarity_tolerance, overlap_tolerance);
            store_rotated_model(rotated_modelOBJ);
        }
    }

    define_buildings(save_buildings); //define and store buildings in "save_buildings"
    cout << "Number of buildings " << buildings.size() << endl;

    define_ground_surfaces(h_gs, theta_gs);
    save2obj_gs(save_gs); //save ground surfaces in save_gs


    //IV. Geometric validations
    //IV.a. Topological relationships
    topological_relationships(topocheckOBJ, topocheckTXT,  th_topo, ground_level, z_ground);

    //IV.b. Sliver triangles
    sliver_triangles(sliversOBJ, sliversTXT, th_slivers);

    //IV.c. Sharp angles
    sharp_angles(sharpanglesOBJ, sharpanglesTXT, th_angles, save_buildings);

    //IV.d. Short edges
    shortedges(shortedgesOBJ, shortedgesTXT, th_edges);

    //IV.e. Overlapping buildings
    overlapping_buildings(overlappingbuildingsOBJ, overlappingbuildingsTXT);

//    int idx2 = filename.find_last_of("/"); // needed when the code is directly executed from backend
//    string name = input_file.substr(idx2+1);

    //V. Mesh parameters definition
    if (n_boxes == 2){ //V.a. with two refinement boxes

        B2_mesh(h_user, z0, N_max, r_min, d_separation, dim, BR, f1, f2, f3, target, x, y); //computes mesh parameters
        B2_write_blockMeshDict(iblockMeshDict, blockMeshDict);
        B2_write_snappyHexMeshDict(name, isnappyHexMeshDict, snappyHexMeshDict);

        computational_domain(domainOBJ); //store computational domain

        building_volume(b_volumeOBJ, b_volumeTXT); //store invalid building volumes

        buildingseparations(b_separations_buildings, b_separations_edges, b_separationsTXT); //store invalid building separations

        if (fTerrain.size() == 0){ //if there is no terrain, ground refinement can be applied

            B2_ground_refinement_noNmax(h_user, z0, N_max, r_min, d_separation, BR, f1, f2, f3, target, x, y); //computes mesh parameters with ground refinement and without number of cells limitation
            B2_write_blockMeshDict(iblockMeshDict, blockMeshDict_noNmax);
            B2_write_snappyHexMeshDict(name, isnappyHexMeshDict, snappyHexMeshDict_noNmax);

            computational_domain(domainOBJ_noNmax); //store computational domain

            building_volume(b_volumeOBJ_noNmax, b_volumeTXT_noNmax); //store invalid building volumes

            buildingseparations(b_separations_buildings_noNmax, b_separations_edges_noNmax, b_separationsTXT_noNmax); //store invalid building separations

            B2_ground_refinement_Nmax(h_user, z0, N_max, r_min, d_separation, BR, f1, f2, f3, target, x, y); //computes mesh parameters with ground refinement and number of cells limitation
            B2_write_blockMeshDict(iblockMeshDict, blockMeshDict_Nmax);
            B2_write_snappyHexMeshDict(name, isnappyHexMeshDict, snappyHexMeshDict_Nmax);

            computational_domain(domainOBJ_Nmax); //store computational domain

            building_volume(b_volumeOBJ_Nmax, b_volumeTXT_Nmax); //store invalid building volumes

            buildingseparations(b_separations_buildings_Nmax, b_separations_edges_Nmax, b_separationsTXT_Nmax); //store invalid building separations
        }

    } else if (n_boxes == 3){ //V.b. with three refinement boxes
        B3_mesh(h_user, z0, N_max, r_min, d_separation, BR, dim, f1, f2, f3, target, x, y); //computes mesh parameters
        B3_write_blockMeshDict(iblockMeshDict, blockMeshDict);
        B3_write_snappyHexMeshDict(name, isnappyHexMeshDict, snappyHexMeshDict);

        computational_domain(domainOBJ);

        //Invalid building volumes
        building_volume(b_volumeOBJ, b_volumeTXT);

        //Invalid building separations
        buildingseparations(b_separations_buildings, b_separations_edges, b_separationsTXT);

        if (fTerrain.size() == 0){
            B3_ground_refinement_noNmax(h_user, z0, N_max, r_min, d_separation, BR, f1, f2, f3, target, x, y);
            B3_write_blockMeshDict(iblockMeshDict, blockMeshDict_noNmax);
            B3_write_snappyHexMeshDict(name, isnappyHexMeshDict, snappyHexMeshDict_noNmax);
            computational_domain(domainOBJ_noNmax);

            //Invalid building volumes
            building_volume(b_volumeOBJ_noNmax, b_volumeTXT_noNmax);

            //Invalid building separations
            buildingseparations(b_separations_buildings_noNmax, b_separations_edges_noNmax,  b_separationsTXT_noNmax);

            B3_ground_refinement_Nmax(h_user, z0, N_max, r_min, d_separation, BR, f1, f2, f3, target, x, y);
            B3_write_blockMeshDict(iblockMeshDict, blockMeshDict_Nmax);
            B3_write_snappyHexMeshDict(name, isnappyHexMeshDict, snappyHexMeshDict_Nmax);

            computational_domain(domainOBJ_Nmax);

            //Invalid building volumes
            building_volume(b_volumeOBJ_Nmax, b_volumeTXT_noNmax);

            //Invalid building separations
            buildingseparations(b_separations_buildings_Nmax, b_separations_edges_Nmax,  b_separationsTXT_Nmax);
        }

    }

    //V.c. Region of Interest (RoI)
    region_of_interest(RoIOBJ, RoITXT, target, x, y); //Save buildings in RoI
    region_of_interest_cylinder(RoICylinder); //Save RoI representing by a cylinder

    //VI. Store results
    output_for_UI(output_for_UI_file);

    return 0;
}
