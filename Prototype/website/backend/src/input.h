#include "definitions.h"

#ifndef BACKEND_INPUT_H
#define BACKEND_INPUT_H

struct Vertex {
    int id;
    Point3 p;
};

struct Face {
    int id;
    int id_building;
    vector<int> v;
    int gs; //indicate whether it is a ground surface
    //used for approximation number of cells with terrain surfaces:
    int box1; //indicates whether the face is located in refinement box 1, 2, or 3 or the remaining area of the domain
    int box2;
    int box3;
    int domain;

};

struct topo { //store topological errors
    int id_building;
    Point3 pBuilding; //vertex of a building
    Point3 pz; //same vertex at terrain height
    Triangle3 triTerrain; //terrain face, when terrain surfaces are included
    double d; //distance between pBuilding and pz or triTerrain, depending on the presence of terrain features
    int error;
    int warning;
};

struct Building {
    int id;
    SMesh mesh;
    vector<Face> faces;
    vector<Point3> pts;
    vector<Face> gs; //ground surfaces
    vector<Point2> gs2D;
    vector<pair<pair<Point2, Point2>, double>> separations;
    vector<pair<pair<Point2, Point2>, double>> invalid_separations;
    int invalid_volume; //indicates whether the building satisfy at least 10 cells per cube root of the building volume (0: yes, 1: no)
    int invalid_separation; //indicates whether the building has separations that are less than 10 cells (0: no, 1: yes)
};

struct Citymodel { //stores some data of the urban model
    vector<Point3> vb; //building vertices
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double dmax; //dimension defining the building domain dimensions
    double zmin_terrain; //lowest z-value of the terrain, needed to determine the lowest z coordinate of the domain for models with terrain
};

struct MeshCFD {
    CGAL::Bbox_3 domain;
    CGAL::Bbox_3 min_domain; //minimum dimensions of domain
    double nx, ny, nz; //number of cells in x, y, and z directions
    double z0; //roughness value
    double N; //number of cells approximation

    Circle2 RoICircle; //region of interest in 2D
    SMesh RoI; //region of interest in 3D

    double hcell; //cell height
    double wcell; //cell width
    double hmin; //minimum cell height

    CGAL::Bbox_3 box1; //refinement boxes 1, 2, and 3
    CGAL::Bbox_3 box2;
    CGAL::Bbox_3 box3;

    int ground_refinement; //indicates whether ground refinement is required
    pair<double, double> gr_pair; //indicates the height and cell up to which the ground refinement must be applied

    double v_box1; //volume of refinement boxes 1, 2, 3, and remaining area of domain
    double v_box2;
    double v_box3;
    double v_domain;
    double volume;

    int BV; //indicates whether there are enough cells to respect at least 10 cells per building cube root volume (0: false, 1: true)
    double N_refined; //N cells with the highest refinement level
};

struct Output { //stores output data and paths to return output files
    int valid; //validity of the input file based on val3dity

    int n_topo; //number of topological errors
    int n_topo_errors; //number of topological errors
    int n_topo_warnings; //number of topological warnings
    int n_sharp_angles; //number of sharp angles
    int n_short_edges; //number of short edges
    int n_slivers; //number of sliver triangles

    CGAL::Bbox_3 domain; //domain dimensions without ground refinement
    CGAL::Bbox_3 domain_GR_noNmax; //domain dimensions with ground refinement and without number of cell limitation
    CGAL::Bbox_3 domain_GR_Nmax; //domain dimensions with ground refinement and number of cell limitation

    pair<double, double> cell; //cell dimensions without ground refinement
    pair<double, double> cell_GR_noNmax; //cell dimensions with ground refinement and without number of cell limitation
    pair<double, double> cell_GR_Nmax; //cell dimensions with ground refinement and number of cell limitation

    double N; //number of cell approximation without ground refinement
    double N_GR_noNmax; //number of cell approximation with ground refinement and without number of cell limitation
    double N_GR_Nmax; //number of cell approximation with ground refinement and number of cell limitation

    pair<double, double> eh; //min and max evaluation height (without ground refinement)
    double eh_GR_noNmax; //evaluation height (with ground refinement and without number of cell limitation)
    double eh_GR_Nmax; //evaluation height (with ground refinement and number of cell limitation)

    int n_RoI; //number of buildings in RoI
    double per_RoI; //percentage of buildings in RoI

    int BV; //indicates whether there are enough cells to respect at least 10 cells per building cube root volume (0: false, 1: true)
    double BV_valid; //percentage of buildings with at least 10 cells per building cube root volume (without ground refinement)
    double BV_valid_GR_noNmax; // " (with ground refinement and without number of cell limitation)
    double BV_valid_GR_Nmax; // " (with ground refinement and number of cell limitation)

    double S_valid; //percentage of buildings with at least 10 cells per building separation (without ground refinement)
    double S_valid_GR_noNmax; // " (with ground refinement and without number of cell limitation)
    double S_valid_GR_Nmax; // " (with ground refinement and number of cell limitation)

    double BR; //blockage ratio defined by Franke and Baklanov (2007)
    double BRL; //blockage ratio in the lateral horizontal direction defined by Blocken (2015)
    double BRH; //blockage ratio in the lateral vertical direction defined by Blocken (2015)
    double hx; //number of largest building or height dimension in domain along the x-direction
    double hy; // " in the y-direction
    double hz; // " in the z-direction

    double BR_noNmax; //blockage ratios data mesh with ground refinement and no number of cell limitation
    double BRL_noNmax;
    double BRH_noNmax;
    double hx_noNmax;
    double hy_noNmax;
    double hz_noNmax;

    double BR_Nmax; //blockage ratios data mesh with ground refinement and number of cell limitation
    double BRL_Nmax;
    double BRH_Nmax;
    double hx_Nmax;
    double hy_Nmax;
    double hz_Nmax;

    double dmax; //largest building dimension or height, dependent on user parameters

    int overlap; //number of overlapping buildings

};


extern Citymodel city; //stores data of the urban area

extern map<int, Face> fiBuildings; //stores building faces (before building definition)
extern map<int, Face> fTerrain; //stores terrain faces
extern map<int, Vertex> vi;  //stores vertices (before building definition)

extern map<int, Building> buildings; //stores buildings
extern map<int, Face> fBuildings; //stores building faces (after building definition)
extern map<int, Vertex> vertices; //stores vertices (after building definition)

extern MeshCFD meshCFD; //stores mesh parameters

extern Output output; //stores output data and paths to return output files

void OBJ_data(string input_file, double wd, double unit); //read and stores data from OBJ input file
void STL_data(string input_file, double wd, double unit); //read and stores data from STL input file

#endif
