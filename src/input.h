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
//    map<int, std::unordered_set<int>> v_common; //map<int, vector<int>> v_common;
//    map<int, std::unordered_set<int>> v_uncommon; //map<int, vector<int>> v_uncommon;
    int gs;
    map<int, vector<int>> v_common;
    map<int, vector<int>> v_uncommon;
//    Triangle3 tri_gs;
};

struct topo {
    int id_building;
    Point3 pBuilding;
    Triangle3 triTerrain;
    double d;
};

struct Building {
    int id;
    SMesh mesh;
    vector<Face> faces;
    vector<Point3> pts;
    vector<Face> gs;
    vector<Point2> gs2D;
    pair<pair<Point2, Point2>, double> d_shortest;
    vector<pair<pair<Point2, Point2>, double>> invalid_separations;
    double d_min;
    int invalid_volume;
    int invalid_separation;

    pair<pair<Point2, Point2>, double> smin;
};

struct Citymodel {
    vector<Point3> vb;

    double xmin, xmax, ymin, ymax, zmin, zmax;
    double hmax;
};

struct MeshCFD {
    CGAL::Bbox_3 domain;
    CGAL::Bbox_3 refinementBox1;
    CGAL::Bbox_3 refinementBox2;
    CGAL::Bbox_3 refinementBox3;
    double nx, ny, nz;
    int nCellsBetweenLevels;
    //double wCell;
    double z0;
    double wmin;
    double min;
//    double wCell1 = wCell/(double)4;
//    double wCell2 = wCell/8;
//    double wCell3 = wCell/16;
//    double wCell4 = wCell/32;

    Circle2 RoICircle;
    SMesh RoI;

    Circle2 RoICircle_intersection;
    SMesh RoI_intersection;

    double hcell;
    double wcell;
    double hmin;

    CGAL::Bbox_3 box1;
    CGAL::Bbox_3 box2;
    CGAL::Bbox_3 box3;



};



extern Citymodel city;

extern map<int, Face> fiBuildings;
extern map<int, Face> fTerrain;
extern map<int, Vertex> vi;

extern map<int, Building> buildings;
extern map<int, Face> fBuildings;
extern map<int, Vertex> vertices;

extern MeshCFD meshCFD;


void OBJ_data(string input_file, double wd);
void STL_data(string input_file, double wd);


#endif
