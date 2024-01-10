#ifndef BACKEND_MESH3_H
#define BACKEND_MESH3_H

void B3_mesh(double h_user, double z0, double N_max, double r_min, double d_separation, int dim, int BR, vector<double> f1, vector<double> f2, vector<double> f3, int target, double x, double y);
void B3_write_blockMeshDict(string initial_file, string output_file);
void B3_write_snappyHexMeshDict(string filename, string initial_file, string output_file);
void B3_ground_refinement_noNmax(double h_user, double z0, double N_max, double r_min, double d_separation, int BR, vector<double> f1, vector<double> f2, vector<double> f3, int target, double x, double y);
void B3_ground_refinement_Nmax(double h_user, double z0, double N_max, double r_min, double d_separation, int BR, vector<double> f1, vector<double> f2, vector<double> f3, int target, double x, double y);

#endif

