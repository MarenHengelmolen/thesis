#ifndef BACKEND_PARAMETERS_H
#define BACKEND_PARAMETERS_H

void parameters(string input_file);

struct User { //stores parameters inserted by users in UI

    string buildings;

    double fd; //flow direction
    double unit; //unit parameter
    int n_boxes; //number of refinement boxes
    double h_user; //target evaluation height
    double z0; //roughness height
    double N_max; //number of cells limit
    double r_min; //cell ratio
    double d_separation; //street width threshold
    int BR; //blockage ratio option
    int dim; //dimension to define domain (hmax or dmax)
    int target; //target building
    double x; //coordinates target building
    double y;
    vector<double> f1; //refinement box dimensions
    vector<double> f2;
    vector<double> f3;

    double snap_tolerance; //val3dity parameters
    double planarity_tolerance;
    double overlap_tolerance;
    double th_topo; //threshold topological relationships
    int ground_level; //indicates whether there is a ground level selected
    double z_ground; //ground level
    double th_slivers; //sliver parameter threshold
    double th_angles; //sharp angles threshold
    double th_edges; //short edges threshold
    double h_gs; //height threshold ground surfaces
    double theta_gs; //angle threshold ground surfaces

    string name; //path to output files
    string save_gs;
    string topocheckOBJ, topocheckTXT;
    string sliversOBJ, sliversTXT;
    string sharpanglesOBJ, sharpanglesTXT;
    string shortedgesOBJ, shortedgesTXT;
    string overlappingbuildingsOBJ, overlappingbuildingsTXT;
    string val3dity_report;
    string rotated_modelOBJ;
    string iblockMeshDict, blockMeshDict, blockMeshDict_noNmax, blockMeshDict_Nmax;
    string isnappyHexMeshDict, snappyHexMeshDict, snappyHexMeshDict_noNmax, snappyHexMeshDict_Nmax;
    string domainOBJ, domainOBJ_noNmax, domainOBJ_Nmax;
    string b_volumeOBJ, b_volumeOBJ_noNmax, b_volumeOBJ_Nmax;
    string b_volumeTXT, b_volumeTXT_noNmax, b_volumeTXT_Nmax;
    string b_separations_buildings, b_separations_edges, b_separations_buildings_noNmax, b_separations_edges_noNmax, b_separations_buildings_Nmax, b_separations_edges_Nmax, b_separationsTXT, b_separationsTXT_noNmax, b_separationsTXT_Nmax;
    string RoIOBJ, RoITXT, RoICylinder;
    string output_for_UI;
    string distances;
};

extern User user;

#endif
