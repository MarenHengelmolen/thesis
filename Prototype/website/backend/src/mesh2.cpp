#include "definitions.h"
#include "input.h"
#include "domain.h"
#include "buildings.h"
#include "refinementboxes.h"
#include "meshtools.h"
#include "RoI.h"

void B2_evaluation_height(int BR, vector<double> f1, vector<double> f2, double h_user, int dim){

    //defines initial cell dimensions based on guidelines for evaluation heights

    double eh_min = (double)h_user/(double)2.5; //minimum cell height (level 3) based on user-defined target height
    double eh_max ;
    if (h_user <= 5){ //maximum cell height (level 3) when user-defined target height lower than 5m
        eh_max = (double)5/(double)2.5;
    } else {
        eh_max = eh_min; //maximum cell height (level 3) when user-defined target height higher than 5m
    }

    //store cell dimensions
    meshCFD.hmin = eh_min*8;
    meshCFD.hcell = eh_max*8;
    meshCFD.wcell = meshCFD.hcell;

    //adjustements of domain and refinement boxes dimensions
    computational_domain_I(BR, dim);
    B2_refinement_boxes(f1, f2);

    cout << "Evaluation height " << endl;
    cout << "hmin: " << meshCFD.hmin << endl;
    cout << "hmax: " << meshCFD.hcell << endl;
    cout << "hcellxwcell: " << meshCFD.hcell << "x" << meshCFD.wcell << endl;
    cout << "\n";
}

void B2_roughness(int BR, vector<double> f1, vector<double> f2){

    //defines cell dimensions based on guidelines for the roughness height

    double rh_min = meshCFD.z0*2*8; //minimum cell height (level 3) based on roughness height

    if (rh_min > meshCFD.hmin){
        //when rh_min higher than current minimum cell height (hmin)
        //store new hmin
        meshCFD.hmin = rh_min;
        if (rh_min > meshCFD.hcell){
            //when rh_min than current minimum cell height (hcell)
            //store new hcell and wcell
            meshCFD.hcell = rh_min;
            meshCFD.wcell = meshCFD.hcell ;

            //adjustements of domain and refinement boxes dimensions
            computational_domain_II(BR);
            B2_refinement_boxes(f1, f2);
        }
    }

    cout << "Roughness " << endl;
    cout << "hmin based on z0: " << rh_min << endl;
    cout << "hmin final: " << meshCFD.hmin << endl;
    cout << "hcellxwcell: " << meshCFD.hcell << "x" << meshCFD.wcell << endl;
    cout << "\n";
}

void B2_max_number_of_cells(int BR, vector<double> f1, vector<double> f2, double Nmax, int gr, double h){

    //adjusts cell dimensions based on the limit for the number of cells

    double N = B2_n_cells_snappy(gr, h); //approximate number of cells
    double reduction = 0; //indicates whether a reduction has already been performed
    double increase = 0; //indicates whether an increase has already been performed

    cout << "Maximum number of cells" << endl;
    cout << "N initial: " << meshCFD.N  << endl;
    CGAL::Bbox_3 domain = meshCFD.domain;
    cout << "V domain: " << abs(domain.xmax()-domain.xmin())*abs(domain.ymax()-domain.ymin())*abs(domain.zmax()-domain.zmin()) << endl;

    while (N > Nmax && increase != 1){
        //when the number of cells approximation is greater than cell limitation
        //increase cell dimensions to not exceed the number of cells limitation
        reduction = 1;
        meshCFD.hcell = meshCFD.hcell + 0.1;
        meshCFD.wcell = meshCFD.wcell + 0.1;
        computational_domain_II(BR);
        B2_refinement_boxes(f1, f2);
        N=B2_n_cells_snappy(gr, h);
    }

    while (N<Nmax && meshCFD.hcell >= meshCFD.hmin && reduction != 1){\
        //when the number of cells approximation is lower than cell limitation
        //reduce cell dimensions to maximize the resolution
        increase = 1;
        meshCFD.hcell = meshCFD.hcell - 0.1;
        meshCFD.wcell = meshCFD.wcell - 0.1;
        computational_domain_II(BR);
        B2_refinement_boxes(f1, f2);
        N  = B2_n_cells_snappy(gr, h);
    }

    if (meshCFD.hcell < meshCFD.hmin || N > Nmax){
        meshCFD.hcell = meshCFD.hcell + 0.1;
        meshCFD.wcell = meshCFD.wcell + 0.1;
        computational_domain_II(BR);
        B2_refinement_boxes(f1, f2);
        N  = B2_n_cells_snappy(gr, h);
    }

    meshCFD.N = N;

    cout << "N: " << meshCFD.N << endl;
    cout << "hcellxwcell: " << meshCFD.hcell << "x" << meshCFD.wcell << endl;
    cout << "\n";
}

double B2_percentage_building_volume(){

    //computes the number of buildings with at least 10 cells per building cube root

    double n_score; //percentage of buildings satisfying this condition

    double n_buildings = buildings.size(); //total number of buildings
    double n_valid = 0; //number of buildings satisfying this condition

    double N = meshCFD.N_refined; //number of cells with the highest refinement level

    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {
        SMesh building = b->second.mesh;
        double v = CGAL::Polygon_mesh_processing::volume(building); //building volume

        if (v >= 0){
            double v_cube_root = cbrt(v); //cube root building volume
            N = N-30; //subtract 30 cells from N, per cube root building volume
            if (v_cube_root/(meshCFD.wcell/8)>=10 && v_cube_root/(meshCFD.hcell/8)>=10){
                //if 10 cells fit in v_cube_root, label building as satisfying this condition
                b->second.invalid_volume = 0;
                n_valid += 1;
            }

            else {
                //otherwise, label building as not satisfying this condition
                b->second.invalid_volume = 1;
            }
        } else {
            //indicates that building do not have a volume
            b->second.invalid_volume = 2;
        }
    }

    if (N>0){
        //when N is higher than 0, there are enough cells for 10 cells per building volume
        meshCFD.BV = 1;
    } else {
        //otherwise, there is not
        meshCFD.BV = 0;
    }

    n_score=n_valid/(double)n_buildings;

    return n_score;
}

void B2_min_cells_building_volume(int BR, vector<double> f1, vector<double> f2, double r_min, double N_max, bool target, double x, double y, int gr, double h){

    //reduces cell width when not all buildings satisfy 10 cells per building cube root
    //while satisfying cell ratio and limitation number of cells

    double score = B2_percentage_building_volume(); //percentage of buildings satisfying this condition
    double ratio = meshCFD.hcell/meshCFD.wcell; //cell ratio
    double N = B2_n_cells_snappy(gr, h); //number of cells approximation

    while (N < N_max && score < 1 && ratio <= r_min){
        meshCFD.wcell = meshCFD.wcell - 0.01;
        ratio = meshCFD.hcell/(double)meshCFD.wcell;
        computational_domain_II(BR);
        B2_refinement_boxes(f1, f2);
        N = B2_n_cells_snappy(gr, h);
        score = B2_percentage_building_volume();
    }

    if (N > N_max || ratio > r_min){
        meshCFD.wcell = meshCFD.wcell + 0.01;
        ratio = meshCFD.hcell/(double)meshCFD.wcell;
        computational_domain_II(BR);
        B2_refinement_boxes(f1, f2);
        N = B2_n_cells_snappy(gr, h);
        score = B2_percentage_building_volume();
    }

    cout << "Cells per cube root building volume" << endl;
    cout << "Number of \"valid\" buildings: " << score << "%" << endl;
    cout << "Cell ratio: " << ratio << endl;
    cout << "hcellxwcell: " << meshCFD.hcell << "x" << meshCFD.wcell << endl;
    cout << "N " << B2_n_cells_snappy(gr, h) << endl;
    cout << "\n";
}

double B2_percentage_building_separation(double d_separation, double wcell){

    //computes number of buildings satisfying 10 cells per building separation

    double score; //percentage satisfying this condition

    int nCellsBetweenLevels = 4;
    double wcell3 = meshCFD.wcell/(double)8; //cell width at refinement level 3
    double min = 2*nCellsBetweenLevels*wcell3+2*wcell; //minimum distance between two buildings

    int n_buildings = buildings.size();
    int n_valid = 0;
    distances_between_buildings(d_separation); //compute building separations

    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {

        SMesh building = b->second.mesh;
        double v = CGAL::Polygon_mesh_processing::volume(building); //building volume

        vector<pair<pair<Point2, Point2>, double>> invalid_separations;

        if (v >= 0) {

            for (auto &sep: b->second.separations){
                if (sep.second != -1){ //building separation is higher than a specific distance defined by users
                    if (sep.second < min){
                        //if building separation is lower than minimum distance computed previously,
                        //this separation is considered invalid
                        invalid_separations.push_back(sep);
                    }
                }
            }

            b->second.invalid_separations = invalid_separations; //store invalid separations per building

            if (invalid_separations.empty()) { //if there are no invalid separations, building is labelled as satisfying this condition
                n_valid += 1;
                b->second.invalid_separation = 0;
            } else {
                b->second.invalid_separation = 1;
            }
        }
    }

    score = n_valid/(double)n_buildings;

    return score;
}

void B2_min_cells_building_separation(int BR, vector<double> f1, vector<double> f2, double r_min, double N_max, bool target, double x, double y, int gr, double h, double wcell){

    //reduces cell width when not all buildings satisfy 10 cells per building separation
    //while satisfying cell ratio and limitation number of cells

    double score = B2_percentage_building_separation(N_max, wcell); //number of buildings satisfying 10 cells per building separation
    double ratio = meshCFD.hcell/meshCFD.wcell; //cell ratio
    double N = B2_n_cells_snappy(gr, h); //number of cells approximation

    while (N < N_max && score < 1 && ratio <= r_min){
        meshCFD.wcell = meshCFD.wcell - 0.01;
        ratio = meshCFD.hcell/(double)meshCFD.wcell;
        computational_domain_II(BR);
        B2_refinement_boxes(f1, f2);
        N = B2_n_cells_snappy(gr, h);
        score = B2_percentage_building_separation(N_max, wcell);
    }

    if (N >  N_max || ratio > r_min){
        meshCFD.wcell = meshCFD.wcell + 0.01;
        ratio = meshCFD.hcell/(double)meshCFD.wcell;
        computational_domain_II(BR);
        B2_refinement_boxes(f1, f2);
        N = B2_n_cells_snappy(gr, h);
        score = B2_percentage_building_separation(N_max, wcell);
    }

    cout << "Cells per building separation" << endl;
    cout << "Number of \"valid\" buildings: " << score << "%" << endl;
    cout << "hcellxwcell: " << meshCFD.hcell << "x" << meshCFD.wcell << endl;
    cout << "N " << B2_n_cells_snappy(gr, h) << endl;
    cout << "\n";

}

pair<double, double> B2_ground_refinement(double h_user){

    //defines height and cell at which ground refinement must be applied

    pair<double, double> gr;

    double hcell3 = meshCFD.hcell/8;
    double h;
    double cell;

    if (2*hcell3 <= h_user && 3*hcell3 > h_user){
        h = 3*hcell3;
        cell = 3;
    } else if (3*hcell3 <= h_user && 4*hcell3 > h_user){
        h = 4*hcell3;
        cell = 4;
    } else if (h_user < 2*hcell3){
        h = 3*hcell3;
        cell = trunc(h_user/hcell3);
    }

    gr = make_pair(h, cell);

    return gr;
}

void B2_zmin_terrain(){

    //defines minimum z-value of the terrain
    vector<double> z_terrain;

    for (auto &f: fTerrain){
        for (auto &v: f.second.v){
            z_terrain.push_back(vi[v].p.z());
        }
    }

    auto min = min_element(z_terrain.begin(), z_terrain.end());
    double zmin = *min;

    city.zmin_terrain = zmin;
}

void B2_mesh(double h_user, double z0, double N_max, double r_min, double d_separation, int dim, int BR, vector<double> f1, vector<double> f2, vector<double> f3, int target, double x, double y){

    //defines mesh without ground refinement

    if (fTerrain.size() != 0){ //when input model has terrain features
        B2_zmin_terrain(); //define minimum z-value of the terrain
    }

    meshCFD.z0 = z0;
    B2_evaluation_height(BR, f1, f2, h_user, dim); //define initial cell dimensions based on guidelines for evaluation heights
    B2_roughness(BR, f1, f2); //define cell dimensions based on guidelines for roughness height
    B2_max_number_of_cells(BR, f1, f2, N_max, 0, 0); //adjust cell dimensions based on the limitation on the number of cells

    B2_min_cells_building_volume(BR, f1, f2, r_min, N_max, target, x, y, 0, 0); //adjust cell dimensions when not all buildings satisfy 10 cells per building cube root

    if(buildings.size()>1){ //when there is more than 1 building
        B2_min_cells_building_separation(BR, f1, f2, r_min, d_separation, target, x, y, 0, 0,  meshCFD.wcell/4); //adjust cell dimensions when not all buildings satisfy 10 cells per building separations
    }

    cout << "Final cell dimensions (without ground refinement):" << endl;
    cout << "hcellxwcell: " << meshCFD.hcell << "x" << meshCFD.wcell << endl;
    cout << "Evaluation height: " << meshCFD.hcell*2.5 << "-" << (meshCFD.hcell/8)*2.5 << endl;
    cout << "Roughness height check: z0 is " << meshCFD.z0 << ", minimum distance between ground and first evaluation node is " << (meshCFD.hcell/(double)8)/(double)2 << endl;
    cout << "Number of cells in blockMesh  " << fixed << n_cells_block() << endl;
    cout << "Volume blockMesh " << fixed << (meshCFD.domain.xmax()-meshCFD.domain.xmin())*(meshCFD.domain.ymax()-meshCFD.domain.ymin())*(meshCFD.domain.zmax()-meshCFD.domain.zmin()) << endl;
    double N_final = B2_n_cells_snappy(0, 0);
    cout << "Number of cells in snappy (Nmax: " << N_max << "): " << fixed << N_final << endl;

    double BV_valid = B2_percentage_building_volume()*100;
    cout << "Number of buildings having 10 cells per cube root volume: " << BV_valid << "%" << endl;

    double S_valid = B2_percentage_building_separation(d_separation, meshCFD.wcell/4)*100;
    cout << "Number of buildings having 10 cells per building separation: " <<  S_valid << "%" << endl;

    //Store results
    output.cell = make_pair(meshCFD.hcell, meshCFD.wcell);
    output.N = N_final;
    output.eh = make_pair(meshCFD.hcell*2.5, (meshCFD.hcell/8)*2.5);
    output.BV = meshCFD.BV;
    output.BV_valid = BV_valid;
    output.S_valid = S_valid;
    output.domain = meshCFD.domain;
    output.BR = ((city.ymax - city.ymin)*(city.zmax - city.zmin))/((meshCFD.domain.ymax()-meshCFD.domain.ymin())*(meshCFD.domain.zmax()-meshCFD.domain.zmin()));
    output.BRL = (city.ymax - city.ymin)/(meshCFD.domain.ymax()-meshCFD.domain.ymin());
    output.BRH = (city.zmax-city.zmin)/(meshCFD.domain.zmax()-meshCFD.domain.zmin());
    output.hx = ((meshCFD.domain.xmax()-meshCFD.domain.xmin())-(city.xmax-city.xmin))/city.dmax;
    output.hy = ((meshCFD.domain.ymax()-meshCFD.domain.ymin())-(city.ymax-city.ymin))/city.dmax;
    output.hz = ((meshCFD.domain.zmax()-meshCFD.domain.zmin())-(city.zmax-city.zmin))/city.dmax;

}

void B2_ground_refinement_noNmax(double h_user, double z0, double N_max, double r_min, double d_separation, int BR, vector<double> f1, vector<double> f2, vector<double> f3, int target, double x, double y){

    //defines mesh with ground refinement and no limitation of number of cells

    pair<double, double> gr_pair = B2_ground_refinement(h_user); //defines height and cell at which ground refinement must be applied
    meshCFD.ground_refinement = 1;
    meshCFD.gr_pair = gr_pair;

    cout << "\nFinal cell dimensions (with ground refinement NO N):" << endl;
    cout << "hcellxwcell: " << meshCFD.hcell << "x" << meshCFD.wcell << endl;
    cout << "Evaluation height: " << meshCFD.hcell*2.5 << "-" << (meshCFD.hcell/8)*2.5 << endl;
    cout << "Roughness height check: z0 is " << meshCFD.z0 << ", minimum distance between ground and first evaluation node is " << (meshCFD.hcell/(double)8)/(double)2 << endl;
    cout << "Number of cells in blockMesh  " <<  fixed << n_cells_block() << endl;
    cout << "Volume blockMesh " << fixed << (meshCFD.domain.xmax()-meshCFD.domain.xmin())*(meshCFD.domain.ymax()-meshCFD.domain.ymin())*(meshCFD.domain.zmax()-meshCFD.domain.zmin()) << endl;

    double N_final = B2_n_cells_snappy(1, meshCFD.gr_pair.second);
    cout << "Number of cells in snappy (Nmax: " << N_max << "): " << N_final << endl;

    double invalid_building_volumes =  B2_percentage_building_volume()*100;
    cout << "Number of buildings having 10 cells per cube root volume: " <<  fixed << B2_percentage_building_volume()*100 << "%" << endl;

    double invalid_building_separations = B2_percentage_building_separation(d_separation, meshCFD.wcell/8)*100;
    cout << "Number of buildings having 10 cells per building separation: " <<  fixed <<invalid_building_separations << "%\n" << endl;

    //Store results
    output.cell_GR_noNmax = make_pair(meshCFD.hcell, meshCFD.wcell);
    output.eh_GR_noNmax = (meshCFD.hcell/8)*2.5;
    output.N_GR_noNmax = N_final;
    output.BV_valid_GR_noNmax = invalid_building_volumes;
    output.S_valid_GR_noNmax = invalid_building_separations;
    output.domain_GR_noNmax = meshCFD.domain;
    output.BR_noNmax = ((city.ymax - city.ymin)*(city.zmax - city.zmin))/((meshCFD.domain.ymax()-meshCFD.domain.ymin())*(meshCFD.domain.zmax()-meshCFD.domain.zmin()));
    output.BRL_noNmax = (city.ymax - city.ymin)/(meshCFD.domain.ymax()-meshCFD.domain.ymin());
    output.BRH_noNmax = (city.zmax-city.zmin)/(meshCFD.domain.zmax()-meshCFD.domain.zmin());
    output.hx_noNmax = ((meshCFD.domain.xmax()-meshCFD.domain.xmin())-(city.xmax-city.xmin))/city.dmax;
    output.hy_noNmax = ((meshCFD.domain.ymax()-meshCFD.domain.ymin())-(city.ymax-city.ymin))/city.dmax;
    output.hz_noNmax = ((meshCFD.domain.zmax()-meshCFD.domain.zmin())-(city.zmax-city.zmin))/city.dmax;
    output.dmax = city.dmax;

}

void B2_ground_refinement_Nmax(double h_user, double z0, double N_max, double r_min, double d_separation, int BR, vector<double> f1, vector<double> f2, vector<double> f3, int target, double x, double y){

    //defines mesh with ground refinement and limitation of number of cells

    B2_max_number_of_cells(BR, f1, f2, N_max, 1, meshCFD.gr_pair.first); //adjust cell dimensions based on the limitation on the number of cells
    B2_min_cells_building_volume(BR, f1, f2, r_min, N_max, target, x, y, 1, meshCFD.gr_pair.first); //adjust cell dimensions when not all buildings satisfy 10 cells per building cube root

    if(buildings.size()>1){//when there is more than 1 building
        B2_min_cells_building_separation(BR, f1, f2, r_min, d_separation, target, x, y, 1, meshCFD.gr_pair.first, meshCFD.wcell/8); //adjust cell dimensions when not all buildings satisfy 10 cells per building separations
    }

    cout << "\nFinal cell dimensions (with ground refinement and N):" << endl;
    cout << "hcellxwcell: " << meshCFD.hcell << "x" << meshCFD.wcell << endl;
    cout << "Evaluation height: " << meshCFD.hcell*2.5 << "-" << (meshCFD.hcell/8)*2.5 << endl;
    cout << "Roughness height check: z0 is " << meshCFD.z0 << ",  minimum distance between ground and first evaluation node is " << (meshCFD.hcell/(double)8)/(double)2 << endl;
    cout << "Number of cells in blockMesh  " <<  fixed <<n_cells_block() << endl;
    cout << "Volume blockMesh " <<  fixed <<(meshCFD.domain.xmax()-meshCFD.domain.xmin())*(meshCFD.domain.ymax()-meshCFD.domain.ymin())*(meshCFD.domain.zmax()-meshCFD.domain.zmin()) << endl;

    double N_final =  B2_n_cells_snappy(1, meshCFD.gr_pair.first);
    cout << "Number of cells snappy (Nmax: " << N_max << "): " << N_final << endl;

    double invalid_building_volumes =  B2_percentage_building_volume()*100;
    cout << "Number of buildings having 10 cells per cube root volume: " <<  fixed <<invalid_building_volumes << "%" << endl;

    double invalid_building_separations = B2_percentage_building_separation(d_separation, meshCFD.wcell/8)*100;
    cout << "Number of buildings having 10 cells per building separation: " <<  fixed <<invalid_building_separations << "%\n" << endl;

    //store results
    output.cell_GR_Nmax = make_pair(meshCFD.hcell, meshCFD.wcell);
    output.eh_GR_Nmax = (meshCFD.hcell/8)*2.5;
    output.N_GR_Nmax = N_final;
    output.BV_valid_GR_Nmax = invalid_building_volumes;
    output.S_valid_GR_Nmax = invalid_building_separations;
    output.domain_GR_Nmax = meshCFD.domain;
    output.BR_Nmax = ((city.ymax - city.ymin)*(city.zmax - city.zmin))/((meshCFD.domain.ymax()-meshCFD.domain.ymin())*(meshCFD.domain.zmax()-meshCFD.domain.zmin()));
    output.BRL_Nmax = (city.ymax - city.ymin)/(meshCFD.domain.ymax()-meshCFD.domain.ymin());
    output.BRH_Nmax = (city.zmax-city.zmin)/(meshCFD.domain.zmax()-meshCFD.domain.zmin());
    output.hx_Nmax = ((meshCFD.domain.xmax()-meshCFD.domain.xmin())-(city.xmax-city.xmin))/city.dmax;
    output.hy_Nmax = ((meshCFD.domain.ymax()-meshCFD.domain.ymin())-(city.ymax-city.ymin))/city.dmax;
    output.hz_Nmax = ((meshCFD.domain.zmax()-meshCFD.domain.zmin())-(city.zmax-city.zmin))/city.dmax;

}


void B2_write_blockMeshDict(string initial_file, string output_file) {

    meshCFD.nx = (meshCFD.domain.xmax()-meshCFD.domain.xmin())/(double)meshCFD.wcell;
    meshCFD.ny = (meshCFD.domain.ymax()-meshCFD.domain.ymin())/(double)meshCFD.wcell;
    meshCFD.nz = (meshCFD.domain.zmax()-meshCFD.domain.zmin())/(double)meshCFD.hcell;

    ifstream input_stream;
    input_stream.open(initial_file);
    ofstream ofile(output_file);

    int i = 0;
    if (input_stream.is_open()) {
        string line, data;

        while (getline(input_stream, line)) {
            istringstream tmp_stream(line);
            tmp_stream >> data;

            if (data == "xMin") {
                ofile << "    xMin  " << setprecision(5) << fixed << meshCFD.domain.xmin() << ";" << endl;
            }
            else if (data == "xMax") {
                ofile << "    xMax  " << setprecision(5) << fixed << meshCFD.domain.xmax() << ";" << endl;
            }

            else if(data == "yMin") {
                ofile << "    yMin  " << setprecision(5) << fixed << meshCFD.domain.ymin() << ";" << endl;
            }
            else if(data == "yMax") {
                ofile << "    yMax  " << setprecision(5) << fixed << meshCFD.domain.ymax() << ";" << endl;
            }

            else if(data == "zMin") {
                ofile << "    zMin  " << setprecision(5) << fixed << meshCFD.domain.zmin() << ";" << endl;
            }
            else if(data == "zMax") {
                ofile << "    zMax  " << setprecision(5) << fixed << meshCFD.domain.zmax() << ";" << endl;
            }

            else if(data == "xCells") {
                ofile << "    xCells  " << setprecision(0) << fixed << (meshCFD.nx) << ";" << endl;
            }
            else if(data == "yCells") {
                ofile << "    yCells  " << setprecision(0) << fixed << (meshCFD.ny) << ";" << endl;
            }
            else if(data == "zCells") {
                ofile << "    zCells  " << setprecision(0) << fixed << (meshCFD.nz) << ";" << endl;
            }

            else {
                ofile << line << endl;
            }
        }
    }

    ofile.close();
}

Point3 inside_point(){
    Point3 p (meshCFD.box1.xmin(), meshCFD.box1.ymin(), meshCFD.box1.zmax());
    return p;
}

void B2_write_snappyHexMeshDict(string filename, string initial_file, string output_file){

    ifstream input_stream;
    input_stream.open(initial_file);
    ofstream ofile(output_file);

    int idx = filename.find_last_of(".");
    string name = filename.substr(0, idx);

    int i = 0;
    if (input_stream.is_open()) {
        string line, data;

        while (getline(input_stream, line)) {
            istringstream tmp_stream(line);
            tmp_stream >> data;

            if (data == "input1") {
                ofile << "        file \"" << filename << "\";" << endl;
            }
            else if (data == "input2") {
                ofile << "        min  (" << setprecision(5) << fixed <<meshCFD.box1.xmin() << " " << setprecision(5) << fixed <<meshCFD.box1.ymin() << " " << setprecision(5) << fixed << meshCFD.box1.zmin() << ");" << endl;
            }

            else if(data == "input3") {
                ofile << "        max  (" << setprecision(5) << fixed << meshCFD.box1.xmax() << " " << setprecision(5) << fixed << meshCFD.box1.ymax() << " " << setprecision(5) << fixed << meshCFD.box1.zmax() << ");" << endl;
            }

            else if (data == "input4") {
                ofile << "        min  (" << meshCFD.box2.xmin() << " " << meshCFD.box2.ymin() << " " << meshCFD.box2.zmin() << ");" << endl;
            }

            else if(data == "input5") {
                ofile << "        max  (" << meshCFD.box2.xmax() << " " << meshCFD.box2.ymax() << " " << meshCFD.box2.zmax() << ");" << endl;
            }

            else if(data == "input8") {
                ofile << "        { file " << "\"" << name << ".eMesh" << "\"; level 1; }" << endl;
            }

            else if(data == "input9") {
                Point3 p = inside_point();
                ofile << "    locationInMesh (" << p.x() << " " << p.y() << " " << p.z()<<");"<< endl;
            }

            else if(data == "input10") {
                if (meshCFD.ground_refinement == 1){

                    ofile << " " << endl;
                    ofile << "    ground" << endl;
                    ofile << "    {" << endl;
                    ofile << "        type            searchablePlane;" << endl;
                    ofile << "        planeType       pointAndNormal;" << endl;
                    ofile << " " << endl;
                    ofile << "        pointAndNormalDict" << endl;
                    ofile << "        {" << endl;
                    ofile << "            basePoint       (" << meshCFD.domain.xmin() << " " << meshCFD.domain.ymin() << " " <<  meshCFD.domain.zmin() << ");" << endl;
                    ofile << "            normal          (0 0 1);" << endl;
                    ofile << "        }" << endl;
                    ofile << "    }" << endl;
                    ofile << " " << endl;

                } else {
                    ofile << " " << endl;
                }

            }

            else if(data == "input11") {
                if (meshCFD.ground_refinement == 1){
                    double h = meshCFD.gr_pair.first;

                    ofile << " " << endl;
                    ofile << "        ground" << endl;
                    ofile << "        {" << endl;
                    ofile << "            mode distance;" << endl;
                    ofile << "            levels ((" << h << " 3));" << endl;
                    ofile << "        }" << endl;
                    ofile << " " << endl;

                } else {
                    ofile << " " << endl;
                }
            }

            else {
                ofile << line << endl;
            }
        }
    }

    ofile.close();
}