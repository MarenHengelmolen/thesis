#ifndef BACKEND_BUILDINGS_H
#define BACKEND_BUILDINGS_H

void define_buildings(string input);
Vector3 normal_vector(Point3 a, Point3 b, Point3 c);
void define_ground_surfaces(double z, double theta);
double building_height(Building b);
vector<size_t> convex_hull(vector<Point2> vertices);
void distances_between_buildings(double threshold);
double perimeter_building(vector<Point2> gsA);
double maximum_dimension(Building b);

#endif
