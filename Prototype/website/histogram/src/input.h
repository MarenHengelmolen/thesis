#include "definitions.h"

#ifndef THESIS_V3_INPUT_H
#define THESIS_V3_INPUT_H

struct Vertex {
    int id;
    Point3 p;
};

struct Face {
    int id;
    int type;
    int id_building;
    vector<int> v;
    vector<int> neigh;
    int gs;
    map<int, vector<int>> v_common;
    map<int, vector<int>> v_uncommon;
};

struct Building {
    int id;
    SMesh mesh;
    vector<Face> faces;
    vector<Point3> pts;
    vector<Face> gs;
    vector<Point2> gs2D;
    vector<pair<pair<Point2, Point2>, double>> separations;
    vector<pair<pair<Point2, Point2>, double>> invalid_separations;
    int invalid_volume;
    int invalid_separation;
    vector<int> neighbours;
    vector<Point2> convex_hull;
};

struct Citymodel {
    vector<Point3> vb;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double hmax;

};

struct MeshCFD {
    CGAL::Bbox_3 domain;
    CGAL::Bbox_3 min_domain;
    double nx, ny, nz;
    double z0;
    double N;

    Circle2 RoICircle;
    SMesh RoI;

    double hcell;
    double wcell;
    double hmin;

    CGAL::Bbox_3 box1;
    CGAL::Bbox_3 box2;
    CGAL::Bbox_3 box3;

    int ground_refinement;
    pair<double, double> gr_pair;
};

struct Output {
    int valid;

    int n_topo;
    int n_sharp_angles;
    int n_short_edges;
    int n_slivers;

    CGAL::Bbox_3 domain;
    CGAL::Bbox_3 domain_GR_noNmax;
    CGAL::Bbox_3 domain_GR_Nmax;

    pair<double, double> cell;
    pair<double, double> cell_GR_noNmax;
    pair<double, double> cell_GR_Nmax;

    double N;
    double N_GR_noNmax;
    double N_GR_Nmax;

    pair<double, double> eh;
    double eh_GR_noNmax;
    double eh_GR_Nmax;

    int n_RoI;
    double per_RoI;

    double BV_valid;
    double BV_valid_GR_noNmax;
    double BV_valid_GR_Nmax;

    double S_valid;
    double S_valid_GR_noNmax;
    double S_valid_GR_Nmax;

    double BR;
    double BRL;
    double BRH;
    double hx;
    double hy;
    double hz;

    double BR_noNmax;
    double BRL_noNmax;
    double BRH_noNmax;
    double hx_noNmax;
    double hy_noNmax;
    double hz_noNmax;

    double BR_Nmax;
    double BRL_Nmax;
    double BRH_Nmax;
    double hx_Nmax;
    double hy_Nmax;
    double hz_Nmax;

    double hmax;

};


extern Citymodel city;

extern map<int, Face> fiBuildings;
extern map<int, Face> fTerrain;
extern map<int, Vertex> vi;

extern map<int, Building> buildings;
extern map<int, Face> fBuildings;
extern map<int, Vertex> vertices;

void OBJ_data(string input_file, double wd, double unit);
void STL_data(string input_file, double wd, double unit);

#endif
