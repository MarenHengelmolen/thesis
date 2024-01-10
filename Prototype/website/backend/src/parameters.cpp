#include "definitions.h"
#include "parameters.h"

User user;

void parameters(string input_file) {

    //Reads parameters inserted by users in UI

    ifstream input_stream;
    input_stream.open(input_file);

    if (input_stream.is_open()) {
        string line, data;
        while (getline(input_stream, line)) {
            istringstream tmp_stream(line);
            tmp_stream >> data;

            if (data == "name") {
                string name;
                tmp_stream >> name;
                user.name = name;
            }

            if (data == "buildings") {
                string buildings;
                tmp_stream >> buildings;
                user.buildings = buildings;
            }

            if (data == "wd") {
                double fd;
                tmp_stream >> fd;
                user.fd = fd;
            }

            if (data == "unit") {
                double unit;
                tmp_stream >> unit;
                user.unit = unit;
            }

            if (data == "n_boxes") {
                double n_boxes;
                tmp_stream >> n_boxes;
                user.n_boxes = n_boxes;
            }

            if (data == "h_user") {
                double h_user;
                tmp_stream >> h_user;
                user.h_user = h_user;
            }

            if (data == "z0") {
                double z0;
                tmp_stream >> z0;
                user.z0 = z0;
            }

            if (data == "N_max") {
                double N_max;
                tmp_stream >> N_max;
                user.N_max = N_max;
            }

            if (data == "r_min") {
                double r_min;
                tmp_stream >> r_min;
                user.r_min = r_min;
            }

            if (data == "d_separation") {
                double d_separation;
                tmp_stream >> d_separation;
                user.d_separation = d_separation;
            }

            if (data == "BR") {
                int BR;
                tmp_stream >> BR;
                user.BR = BR;
            }

            if (data == "dim") {
                int dim;
                tmp_stream >> dim;
                user.dim = dim;
            }
            if (data == "target") {
                int target;
                tmp_stream >> target;
                user.target = target;
            }

            if (data == "x") {
                double x;
                tmp_stream >> x;
                user.x = x;
            }

            if (data == "y") {
                double y;
                tmp_stream >> y;
                user.y = y;
            }

            if (data == "f1") {
                double fin, fout, flat, ftop;
                tmp_stream >> fin >> fout >> flat >> ftop;
                user.f1 = {fin, fout, flat, ftop};
            }

            if (data == "f2") {
                double fin, fout, flat, ftop;
                tmp_stream >> fin >> fout >> flat >> ftop;
                user.f2 = {fin, fout, flat, ftop};
            }

            if (data == "f3") {
                double fin, fout, flat, ftop;
                tmp_stream >> fin >> fout >> flat >> ftop;
                user.f3 = {fin, fout, flat, ftop};
            }

            if (data == "snap_tolerance") {
                double snap_tolerance;
                tmp_stream >> snap_tolerance;
                user.snap_tolerance = snap_tolerance;
            }

            if (data == "planarity_tolerance") {
                double planarity_tolerance;
                tmp_stream >> planarity_tolerance;
                user.planarity_tolerance = planarity_tolerance;
            }

            if (data == "overlap_tolerance") {
                double overlap_tolerance;
                tmp_stream >> overlap_tolerance;
                user.overlap_tolerance = overlap_tolerance;
            }

            if (data == "th_topo") {
                double th_topo;
                tmp_stream >> th_topo;
                user.th_topo = th_topo;
            }

            if (data == "ground_level") {
                int ground_level;
                tmp_stream >> ground_level;
                user.ground_level = ground_level;
            }

            if (data == "z_ground") {
                double z_ground;
                tmp_stream >> z_ground;
                user.z_ground = z_ground;
            }

            if (data == "th_slivers") {
                double th_slivers;
                tmp_stream >> th_slivers;
                user.th_slivers = th_slivers;
            }

            if (data == "th_angles") {
                double th_angles;
                tmp_stream >> th_angles;
                user.th_angles = th_angles;
            }

            if (data == "th_edges") {
                double th_edges;
                tmp_stream >> th_edges;
                user.th_edges = th_edges;
            }

            if (data == "h_gs") {
                double h_gs;
                tmp_stream >> h_gs;
                user.h_gs = h_gs;
            }

            if (data == "theta_gs") {
                double theta_gs;
                tmp_stream >> theta_gs;
                user.theta_gs = theta_gs;
            }

            if (data == "save_gs") {
                string save_gs;
                tmp_stream >> save_gs;
                user.save_gs = save_gs;
            }

            if (data == "topocheck") {
                string topocheckOBJ, topocheckTXT;
                tmp_stream >> topocheckOBJ >> topocheckTXT;
                user.topocheckOBJ = topocheckOBJ;
                user.topocheckTXT = topocheckTXT;
            }

            if (data == "slivers") {
                string sliversOBJ, sliversTXT;
                tmp_stream >> sliversOBJ >> sliversTXT;
                user.sliversOBJ = sliversOBJ;
                user.sliversTXT = sliversTXT;
            }

            if (data == "sharpangles") {
                string sharpanglesOBJ, sharpanglesTXT;
                tmp_stream >> sharpanglesOBJ >> sharpanglesTXT;
                user.sharpanglesOBJ = sharpanglesOBJ;
                user.sharpanglesTXT = sharpanglesTXT;
            }

            if (data == "shortedges") {
                string shortedgesOBJ, shortedgesTXT;
                tmp_stream >> shortedgesOBJ >> shortedgesTXT;
                user.shortedgesOBJ = shortedgesOBJ;
                user.shortedgesTXT = shortedgesTXT;
            }

            if (data == "overlappingbuildings") {
                string overlappingbuildingsOBJ, overlappingbuildingsTXT;
                tmp_stream >> overlappingbuildingsOBJ >> overlappingbuildingsTXT;
                user.overlappingbuildingsOBJ = overlappingbuildingsOBJ;
                user.overlappingbuildingsTXT = overlappingbuildingsTXT;
            }


            if (data == "val3dity") {
                string val3dity_report;
                tmp_stream >> val3dity_report;
                user.val3dity_report = val3dity_report;
            }

            if (data == "rotated_model") {
                string rotated_modelOBJ;
                tmp_stream >> rotated_modelOBJ;
                user.rotated_modelOBJ = rotated_modelOBJ;
            }

            if (data == "blockMeshDict") {
                string iblockMeshDict, blockMeshDict, blockMeshDict_noNmax, blockMeshDict_Nmax;
                tmp_stream >> iblockMeshDict >> blockMeshDict >> blockMeshDict_noNmax >> blockMeshDict_Nmax;
                user.iblockMeshDict = iblockMeshDict;
                user.blockMeshDict = blockMeshDict;
                user.blockMeshDict_noNmax = blockMeshDict_noNmax;
                user.blockMeshDict_Nmax = blockMeshDict_Nmax;
            }

            if (data == "snappyHexMeshDict") {
                string isnappyHexMeshDict, snappyHexMeshDict, snappyHexMeshDict_noNmax, snappyHexMeshDict_Nmax;
                tmp_stream >> isnappyHexMeshDict >> snappyHexMeshDict >> snappyHexMeshDict_noNmax >> snappyHexMeshDict_Nmax;
                user.isnappyHexMeshDict = isnappyHexMeshDict;
                user.snappyHexMeshDict = snappyHexMeshDict;
                user.snappyHexMeshDict_noNmax = snappyHexMeshDict_noNmax;
                user.snappyHexMeshDict_Nmax = snappyHexMeshDict_Nmax;
            }

            if (data == "domain") {
                string domainOBJ, domainOBJ_noNmax, domainOBJ_Nmax;
                tmp_stream >> domainOBJ >> domainOBJ_noNmax >> domainOBJ_Nmax;
                user.domainOBJ = domainOBJ;
                user.domainOBJ_noNmax = domainOBJ_noNmax;
                user.domainOBJ_Nmax = domainOBJ_Nmax;
            }

            if (data == "b_volume") {
                string b_volumeOBJ, b_volumeOBJ_noNmax, b_volumeOBJ_Nmax, b_volumeTXT, b_volumeTXT_noNmax, b_volumeTXT_Nmax;
                tmp_stream >> b_volumeOBJ >> b_volumeOBJ_noNmax >> b_volumeOBJ_Nmax >>  b_volumeTXT, b_volumeTXT_noNmax, b_volumeTXT_Nmax;
                user.b_volumeOBJ = b_volumeOBJ;
                user.b_volumeOBJ_noNmax = b_volumeOBJ_noNmax;
                user.b_volumeOBJ_Nmax = b_volumeOBJ_Nmax;
                user.b_volumeTXT = b_volumeTXT;
                user.b_volumeTXT_noNmax = b_volumeTXT_noNmax;
                user.b_volumeTXT_Nmax = b_volumeTXT_Nmax;
            }

            if (data == "b_separations") {
                string b_separations_buildings, b_separations_edges, b_separations_buildings_noNmax, b_separations_edges_noNmax, b_separations_buildings_Nmax, b_separations_edges_Nmax, b_separationsTXT, b_separationsTXT_noNmax, b_separationsTXT_Nmax;;
                tmp_stream >> b_separations_buildings >> b_separations_edges >> b_separations_buildings_noNmax >> b_separations_edges_noNmax >> b_separations_buildings_Nmax >> b_separations_edges_Nmax >> b_separationsTXT >> b_separationsTXT_noNmax >> b_separationsTXT_Nmax;;
                user.b_separations_buildings = b_separations_buildings;
                user.b_separations_edges = b_separations_edges;
                user.b_separations_buildings_noNmax = b_separations_buildings_noNmax;
                user.b_separations_edges_noNmax = b_separations_edges_noNmax;
                user.b_separations_buildings_Nmax = b_separations_buildings_Nmax;
                user.b_separations_edges_Nmax= b_separations_edges_Nmax;
                user.b_separationsTXT = b_separationsTXT;
                user.b_separationsTXT_noNmax = b_separationsTXT_noNmax;
                user.b_separationsTXT_Nmax = b_separationsTXT_Nmax;
            }

            if (data == "RoI") {
                string RoIOBJ, RoITXT, RoICylinder;
                tmp_stream >> RoIOBJ >> RoITXT >> RoICylinder;
                user.RoIOBJ = RoIOBJ;
                user.RoITXT = RoITXT;
                user.RoICylinder = RoICylinder;
            }

            if (data == "output_for_UI"){
                string output_for_UI;
                tmp_stream >> output_for_UI;
                user.output_for_UI = output_for_UI;
            }

            if (data == "distances"){
                string distances;
                tmp_stream >> distances;
                user.distances = distances;
            }

        }
    }
};

