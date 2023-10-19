#ifndef THESIS_V3_GEOMETRY_H
#define THESIS_V3_GEOMETRY_H

vector<pair<pair<Point3, Point3>, double>> identify_short_edges(double threshold);
vector<pair<pair<int, int>, double>> identify_sharp_angles(double threshold);
vector<pair<int, double>> identify_sliver_triangles(double threshold);
vector<pair<pair<Point2, Point2>, double>> identify_small_distances(double threshold);

map<int, topo> check_topological_relationship(double threshold);

#endif
