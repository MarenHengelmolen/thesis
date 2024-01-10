#ifndef BACKEND_GEOMETRY_H
#define BACKEND_GEOMETRY_H

struct FaceSA {
    int id;
    vector<Point3> v;
};

int run_val3dity(string input_file, string val3dity_report, double snap_tolerance, double planarity_tolerance, double overlap_tolerance);
vector<pair<pair<Point3, Point3>, double>> identify_short_edges(double threshold);
vector<pair<pair<int, int>, double>> identify_sharp_angles(double threshold, string input, string OBJfilename, string TXTfilename);
vector<pair<int, double>> identify_sliver_triangles(double threshold);
map<int, topo> check_topological_relationship(double threshold, int ground_level, double z_ground);
vector<pair<int, int>> identify_overlapping_buildings();

#endif
