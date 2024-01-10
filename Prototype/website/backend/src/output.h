#ifndef BACKEND_OUTPUT_H
#define BACKEND_OUTPUT_H

void store_rotated_model(string output_file);
void save2obj_iBuildings(string output_file);
void save2obj_gs(string filename);
void computational_domain(string OBJfilename);
void shortedges(string OBJfilename, string TXTfilename, double threshold);
void sharp_angles(string OBJfilename, string TXTfilename,  double threshold, string input);
void sliver_triangles(string OBJfilename, string TXTfilename,  double threshold);
void topological_relationships(string OBJfilename, string TXTfilename,  double threshold, int ground_level, double z_ground);
void region_of_interest(string OBJfilename, string TXTfilename, bool target, double x, double y);
void region_of_interest_cylinder(string OBJfilename);
void buildingseparations(string OBJfilename1, string OBJfilename2, string TXTfile);
void building_volume(string output_file, string TXTfile);
void overlapping_buildings(string OBJfilename, string TXTfilename);
void output_for_UI(string TXTfilename);

#endif
