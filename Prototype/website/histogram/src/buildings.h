#ifndef THESIS_V3_BUILDINGS_H
#define THESIS_V3_BUILDINGS_H

void define_buildings(string input);
Vector3 normal_vector(Point3 a, Point3 b, Point3 c);
void define_ground_surfaces(double z, double theta);
vector<size_t> convex_hull(vector<Point2> vertices);
void distances_for_histogram(string file_with_distances);

#endif
