#ifndef BACKEND_MESH2_H
#define BACKEND_MESH2_H

void B2_mesh(double h_user, double z0, double N_max, double r_min, double d_separation, int dim, int BR, vector<double> f1, vector<double> f2, vector<double> f3, int target, double x, double y);
void B2_write_blockMeshDict(string initial_file, string output_file);
void B2_write_snappyHexMeshDict(string filename, string initial_file, string output_file);
void B2_ground_refinement_noNmax(double h_user, double z0, double N_max, double r_min, double d_separation, int BR, vector<double> f1, vector<double> f2, vector<double> f3, int target, double x, double y);
void B2_ground_refinement_Nmax(double h_user, double z0, double N_max, double r_min, double d_separation, int BR, vector<double> f1, vector<double> f2, vector<double> f3, int target, double x, double y);

#endif
