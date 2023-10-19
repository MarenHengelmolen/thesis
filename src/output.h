#ifndef THESIS_V3_OUTPUT_H
#define THESIS_V3_OUTPUT_H

#include "definitions.h"

void save2obj_terrain(string output_file);
void save2obj_iBuildings(string output_file);
void save2obj_building(string output_file, int b);

void save2obj_gs(string filename);

void computational_domain(string OBJfilename, string TXTfilename, double wCell);

void shortedges(string OBJfilename, string TXTfilename, double threshold);
void sharp_angles(string OBJfilename, string TXTfilename,  double threshold);
void sliver_triangles(string OBJfilename, string TXTfilename,  double threshold);
void distances_between_buildings(string OBJfilename, string TXTfilename, double threshold);

void topological_relationships(string OBJfilename, string TXTfilename,  double threshold);

void region_of_interest(string OBJfilename, string TXTfilename, bool target, double x, double y);

void buildingseparations(string OBJfilename1, string OBJfilename2);
void building_volume(string output_file);
#endif
