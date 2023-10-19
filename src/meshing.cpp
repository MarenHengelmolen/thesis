#include "definitions.h"
#include "input.h"

vector<int> Buildings_in_RoI(bool target, double x, double y);

double round5(double d){
    d = round(d*100000)/100000;
    return d;
}

pair<double, double> fit_to_cell_dimensions(double min, double max, double dcell, bool z) {

    double check = abs(max-min)/(double)dcell;

    if (check - trunc(check) > 0){

        min = round5(min);
        max = round5(max);

        double l = abs(max-min);
        int n = trunc(l/(double)dcell);
        double k = round5(dcell*n);

        double c = (min+max)*0.5;

        if (z == 1){
            min = min;
            max = min+k;
        } else {
            min = c-k*0.5;
            max = c+k*0.5;
        }
    }

    return pair(min, max);
}

pair<double, double> fit_to_cell_dimensions_RB(double min, double max, double dcell, bool z) {

    double check = abs(max-min)/(double)dcell;

    if (check - trunc(check) > 0){

        min = round5(min);
        max = round5(max);

        double l = abs(max-min);
        int n = ceil(l/(double)dcell);
        double k = round5(dcell*n);

        double c = (min+max)*0.5;

        if (z == 1){
            min = min;
            max = min+k;
        } else {
            min = c-k*0.5;
            max = c+k*0.5;
        }
    }

    return pair(min, max);
}

void get_city_boundaries() {
    auto x_minmax = minmax_element(city.vb.begin(), city.vb.end(),
                                   [](const Point3 &p1, const Point3 &p2) { return p1.x() < p2.x(); });
    auto y_minmax = minmax_element(city.vb.begin(), city.vb.end(),
                                   [](const Point3 &p1, const Point3 &p2) { return p1.y() < p2.y(); });
    auto z_minmax = minmax_element(city.vb.begin(), city.vb.end(),
                                   [](const Point3 &p1, const Point3 &p2) { return p1.z() < p2.z(); });

    city.xmin = x_minmax.first->x();
    city.xmax = x_minmax.second->x();
    city.ymin = y_minmax.first->y();
    city.ymax = y_minmax.second->y();
    city.zmin = z_minmax.first->z();
    city.zmax = z_minmax.second->z();

//    cout << "xmin " << city.xmin << endl;
//    cout << "xmax " << city.xmax << endl;
//    cout << "ymin " << city.ymin << endl;
//    cout << "ymax " << city.ymax << endl;
//    cout << "zmin " << city.zmin << endl;
//    cout << "zmax " << city.zmax << endl;
}

double building_height(Building b){
    vector<double> z_pts;
    for (auto &v: b.pts){
        double z = v.z();
        z_pts.push_back(z);
    }

    auto z_minmax = minmax_element(z_pts.begin(), z_pts.end(),[](const double& p1, const double& p2) { return p1 < p2; });
    double height = *z_minmax.second-*z_minmax.first;

    return height;
}

void highest_building() {

    vector<double> lengths;
    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++){
        double length = building_height(b->second);
        lengths.push_back(length);
    }

    auto max = max_element(lengths.begin(), lengths.end());
    double hmax = *max;

    city.hmax = hmax;
}

void BR_option() {
//    cout << "BR OPTION" << endl;
    double L_urban = abs(city.ymax - city.ymin);
    double L_domain = abs(meshCFD.domain.ymax()-meshCFD.domain.ymin());
//    cout << "ymin " << meshCFD.domain.ymin() << endl;
//    cout << "ymax " << meshCFD.domain.ymax() << endl;
    double BR_horizontal = L_urban/L_domain;
//    cout << "L_urban " << L_urban << endl;
//    cout << "L_domain " << L_domain << endl;
//    cout << "BR horizontal " << BR_horizontal << endl;

    double H_urban = abs(city.zmax-city.zmin);
    double H_domain = abs(meshCFD.domain.zmax()-meshCFD.domain.zmin());
    double BR_vertical = H_urban/H_domain;
//    cout << "H_urban " << H_urban << endl;
//    cout << "H_domain " << H_domain << endl;
//    cout << "BR vertical " << BR_vertical << endl;

    double xmin_i = meshCFD.domain.xmin();
    double xmax_i = meshCFD.domain.xmax();
    double ymin_i = meshCFD.domain.ymin();
    double ymax_i = meshCFD.domain.ymax();
    double zmin_i = meshCFD.domain.zmin();
    double zmax_i = meshCFD.domain.zmax();

    double xmin_final, xmax_final;
    double ymin_final, ymax_final;
    double zmin_final, zmax_final;

    double ymin_BRh = ymin_i;
    double ymax_BRh = ymax_i;

    while (BR_horizontal>0.17) {
        ymin_BRh = ymin_BRh-0.01;
        ymax_BRh = ymax_BRh+0.01;
        L_domain = abs(ymax_BRh-ymin_BRh);
        BR_horizontal = (double)L_urban/(double)L_domain;
    }
//    cout << "L domain " << L_domain << endl;
//    cout << "BR horizontal " << BR_horizontal << endl;

    if (L_domain <= abs(xmax_i-xmin_i)) {

        double zmax_BRv = zmax_i;

        while (BR_vertical > 0.17) {
            zmax_BRv = zmax_BRv + 0.01;
//            cout << "zmax_BRv " << zmax_BRv << endl;
            H_domain = abs(zmax_BRv - zmin_i);
            BR_vertical = H_urban / H_domain;
        }

//        cout << "BR vertical " << BR_vertical << endl;

        xmin_final = xmin_i;
        xmax_final = xmax_i;
        ymin_final = ymin_BRh;
        ymax_final = ymax_BRh;
        zmin_final = zmin_i;
        zmax_final = zmax_BRv;

    } else {
        double ymin_BR = ymin_i;
        double ymax_BR = ymax_i;
        double zmax_BR = zmax_i;

        double A_urban = abs(city.ymax - city.ymin) * abs(city.zmax - city.zmin);
        double A_domain = abs(ymax_BR - ymin_BR) * abs(zmax_BR - zmin_i);
        double BR = A_urban / A_domain;

//        cout << "BR " << BR << endl;
//        cout << "A_urban " << A_urban << endl;
        while (BR > 0.03) {
            ymin_BR = ymin_BR - 0.01;
            ymax_BR = ymax_BR + 0.01;
            zmax_BR = zmax_BR + 0.01;
            A_domain = abs(ymax_BR- ymin_BR) * abs(zmax_BR - zmin_i);
            BR = A_urban / A_domain;
        }

//        cout << "BR " << BR << endl;

        xmin_final = xmin_i;
        xmax_final = xmax_i;
        ymin_final = ymin_BR;
        ymax_final = ymax_BR;
        zmin_final = zmin_i;
        zmax_final = zmax_BR;
    }

    pair<double, double> xComp = fit_to_cell_dimensions(xmin_final, xmax_final, meshCFD.wcell, 0);
    pair<double, double> yComp = fit_to_cell_dimensions(ymin_final, ymax_final, meshCFD.wcell, 0);
    pair<double, double> zComp = fit_to_cell_dimensions(zmin_final, zmax_final, meshCFD.hcell, 1);

    xmin_final = xComp.first;
    xmax_final = xComp.second;
    ymin_final = yComp.first;
    ymax_final = yComp.second;
    zmin_final = zComp.first;
    zmax_final = zComp.second;

    CGAL::Bbox_3 domain = CGAL::Bbox_3(xmin_final, ymin_final, zmin_final, xmax_final, ymax_final, zmax_final);

//    cout << "xmin_final " << xmin_final << endl;
//    cout << "xmax_final " << xmax_final << endl;
//    cout << "ymin_final " << ymin_final << endl;
//    cout << "ymax_final " << ymax_final << endl;
//    cout << "zmin_final " << zmin_final << endl;
//    cout << "zmax_final " << zmax_final << endl;

    meshCFD.domain = domain;
}


void computational_domain_I(int BR){
    get_city_boundaries();
    highest_building();

    double hcell = meshCFD.hcell;
    double wcell = meshCFD.wcell;

    double hmax = city.hmax;
    double xmin = city.xmin - 5*hmax;
    double xmax = city.xmax + 15*hmax;
    double ymin = city.ymin - 5*hmax;
    double ymax = city.ymax + 5*hmax;
    double zmin = city.zmin;
    double zmax = city.zmax + 5*hmax;

    pair<double, double> xComp = fit_to_cell_dimensions(xmin, xmax, wcell, 0);
    pair<double, double> yComp = fit_to_cell_dimensions(ymin, ymax, wcell, 0);
    pair<double, double> zComp = fit_to_cell_dimensions(zmin, zmax, hcell, 1);

    xmin = xComp.first;
    xmax = xComp.second;
    ymin = yComp.first;
    ymax = yComp.second;
    zmin = zComp.first;
    zmax = zComp.second;

    CGAL::Bbox_3 domain(xmin, ymin, zmin, xmax, ymax, zmax);
//    cout << "Without BR I " << endl;
//    cout << "xmin_final " << xmin << endl;
//    cout << "xmax_final " << xmax << endl;
//    cout << "ymin_final " << ymin << endl;
//    cout << "ymax_final " << ymax << endl;
//    cout << "zmin_final " << zmin << endl;
//    cout << "zmax_final " << zmax << endl;
    meshCFD.domain = domain;

    if (BR == 1){
        BR_option();
    }
}

void computational_domain_II(int BR){

    double hcell = meshCFD.hcell;
    double wcell = meshCFD.wcell;

    double xmin = meshCFD.domain.xmin();
    double xmax = meshCFD.domain.xmax();
    double ymin = meshCFD.domain.ymin();
    double ymax = meshCFD.domain.ymax();
    double zmin = meshCFD.domain.zmin();
    double zmax = meshCFD.domain.zmax();

    pair<double, double> xComp = fit_to_cell_dimensions(xmin, xmax, wcell, 0);
    pair<double, double> yComp = fit_to_cell_dimensions(ymin, ymax, wcell, 0);
    pair<double, double> zComp = fit_to_cell_dimensions(zmin, zmax, hcell, 1);

    xmin = xComp.first;
    xmax = xComp.second;
    ymin = yComp.first;
    ymax = yComp.second;
    zmin = zComp.first;
    zmax = zComp.second;

    CGAL::Bbox_3 domain(xmin, ymin, zmin, xmax, ymax, zmax);

    meshCFD.domain = domain;

//    cout << "Without BR II " << endl;
//    cout << "xmin_final " << xmin << endl;
//    cout << "xmax_final " << xmax << endl;
//    cout << "ymin_final " << ymin << endl;
//    cout << "ymax_final " << ymax << endl;
//    cout << "zmin_final " << zmin << endl;
//    cout << "zmax_final " << zmax << endl;

    if (BR == 1){
        BR_option();
    }
}


CGAL::Bbox_3 refinement_box(double finlet, double foutlet, double flateral, double ftop) {

    double hcell = meshCFD.hcell;
    double wcell = meshCFD.wcell;

    double ymin = city.ymin-flateral*city.hmax;
    double ymax = city.ymax+flateral*city.hmax;
    double xmin = city.xmin-finlet*city.hmax;
    double xmax = city.xmax+foutlet*city.hmax;
    double zmin = city.zmin;
    double zmax = city.zmax+city.hmax*ftop;

    pair<double, double> xComp = fit_to_cell_dimensions_RB(xmin, xmax, wcell, 0);
    pair<double, double> yComp = fit_to_cell_dimensions_RB(ymin, ymax, wcell, 0);
    pair<double, double> zComp = fit_to_cell_dimensions_RB(zmin, zmax, hcell, 1);

    xmin = xComp.first;
    xmax = xComp.second;
    ymin = yComp.first;
    ymax = yComp.second;
    zmin = zComp.first;
    zmax = zComp.second;

    CGAL::Bbox_3 bbox(xmin, ymin, zmin, xmax, ymax, zmax);

    return bbox;
}

void refinement_boxes(vector<double> f1, vector<double> f2, vector<double> f3){

    CGAL::Bbox_3 box1 = refinement_box(f1[0], f1[1], f1[2], f1[3]);
    CGAL::Bbox_3 box2 = refinement_box(f2[0], f2[1], f2[2], f2[3]);
    CGAL::Bbox_3 box3 = refinement_box(f3[0], f3[1], f3[2], f3[3]);

    meshCFD.box1 = box1;
    meshCFD.box2 = box2;
    meshCFD.box3 = box3;

//    cout << "REFINEMENT BOX 1" << endl;
//    cout << box1.xmin() << endl;
//    cout << box1.xmax() << endl;
//    cout << box1.ymin() << endl;
//    cout << box1.ymax() << endl;
//    cout << box1.zmin() << endl;
//    cout << box1.zmax() << endl;
//
//    cout << "REFINEMENT BOX 2" << endl;
//    cout << box2.xmin() << endl;
//    cout << box2.xmax() << endl;
//    cout << box2.ymin() << endl;
//    cout << box2.ymax() << endl;
//    cout << box2.zmin() << endl;
//    cout << box2.zmax() << endl;
//
//    cout << "REFINEMENT BOX 3" << endl;
//    cout << box3.xmin() << endl;
//    cout << box3.xmax() << endl;
//    cout << box3.ymin() << endl;
//    cout << box3.ymax() << endl;
//    cout << box3.zmin() << endl;
//    cout << box3.zmax() << endl;
}

void evaluation_height(int BR, vector<double> f1, vector<double> f2, vector<double> f3){

    double eh_min = (double)1.5/(double)2.5;
    double eh_max = (double)5/(double)2.5;

    meshCFD.hmin = eh_min;
    meshCFD.hcell = eh_max;
    meshCFD.wcell = meshCFD.hcell;

    cout << "Evaluation height " << endl;
    cout << "min: " << meshCFD.hmin << endl;
    cout << "max: " << meshCFD.hcell << endl;
    cout << "hcellxwcell: " << meshCFD.hcell << "x" << meshCFD.wcell << endl;
    cout << "\n";

    computational_domain_I(BR);
    refinement_boxes(f1, f2, f3);
}

void roughness(int BR, vector<double> f1, vector<double> f2, vector<double> f3){
    double rh_min = meshCFD.z0*2*16;

    if (rh_min > meshCFD.hmin){
        meshCFD.hmin = rh_min;
        if (rh_min > meshCFD.hcell){
            meshCFD.hcell = rh_min;
            meshCFD.wcell = meshCFD.hcell ;

            computational_domain_II(BR);
            refinement_boxes(f1, f2, f3);
        }
    }

    cout << "Roughness " << endl;
    cout << "min roughness: " << rh_min << endl;
    cout << "hmin: " << meshCFD.hmin << endl;
    cout << "hcellxwcell: " << meshCFD.hcell << "x" << meshCFD.wcell << endl;
    cout << "\n";
}

double n_cells(){

    CGAL::Bbox_3 domain = meshCFD.domain;
    double V_domain = abs(domain.xmax()-domain.xmin())*abs(domain.ymax()-domain.ymin())*abs(domain.zmax()-domain.zmin());
    double V_cell = meshCFD.wcell*meshCFD.wcell*meshCFD.hcell;
    double N = (double)V_domain/(double)V_cell;

    return N;
}

void max_number_of_cells(int BR, vector<double> f1, vector<double> f2, vector<double> f3, double threshold){

    double N = n_cells();
    double reduction = 0;
    double increase = 0;

    cout << "Maximum number of cells" << endl;
    cout << "N initial: " << N << endl;
    CGAL::Bbox_3 domain = meshCFD.domain;
    cout << "V domain: " << abs(domain.xmax()-domain.xmin())*abs(domain.ymax()-domain.ymin())*abs(domain.zmax()-domain.zmin()) << endl;

    while (N>threshold && increase != 1){
        reduction = 1;
        meshCFD.hcell = meshCFD.hcell + 0.1;
        meshCFD.wcell = meshCFD.wcell + 0.1;
        computational_domain_II(BR);
        refinement_boxes(f1, f2, f3);
        N = n_cells();
    }

    while (N<threshold && meshCFD.hcell >= meshCFD.hmin && reduction != 1){
        increase = 1;
        meshCFD.hcell = meshCFD.hcell - 0.1;
        meshCFD.wcell = meshCFD.wcell - 0.1;
        computational_domain_II(BR);
        refinement_boxes(f1, f2, f3);
        N = n_cells();
    }

    cout << "N: " << N << endl;
    cout << "hcellxwcell: " << meshCFD.hcell << "x" << meshCFD.wcell << endl;
    cout << "\n";
}

pair<double, double> score_cells_per_cube_root_volume(){

    double n_score;
    double v_score;

    double n_buildings = buildings.size();
    double n_valid = 0;

    double v_buildings = 0;
    double v_valid = 0;

    for (map<int, Building>::iterator db=buildings.begin(); db!=buildings.end(); db++) {
        SMesh building = db->second.mesh;
        double v = CGAL::Polygon_mesh_processing::volume(building);
        if (v >= 0){
            v_buildings += v;
            double v_cube_root = cbrt(v);

            if (v_cube_root/(meshCFD.wcell/16)>=10){
                db->second.invalid_volume = 0;
                n_valid += 1;
                v_valid += v;
            }

            if (v_cube_root/(meshCFD.wcell/16)<10){
                db->second.invalid_volume = 1;
            }
        } else {
            db->second.invalid_volume = 2;
        }
    }

    n_score=n_valid/(double)n_buildings;
    v_score=v_valid/(double)v_buildings;


    pair<double, double> score = make_pair(n_score, v_score);

    return score;
}

pair<double, double> score_cells_per_cube_root_volume_RoI(bool target, double x, double y){

    vector<int> RoIBuildings = Buildings_in_RoI(target, x, y);

    double n_score;
    double v_score;

    double n_buildings = RoIBuildings.size();
    double n_valid = 0;

    double v_buildings = 0;
    double v_valid = 0;

    for (vector<int>::iterator b=RoIBuildings.begin(); b!=RoIBuildings.end(); b++) {
        SMesh building = buildings[*b].mesh;
        double v = CGAL::Polygon_mesh_processing::volume(building);
        v_buildings += v;
        double v_cube_root = cbrt(v);

        if (buildings[*b].invalid_volume == 0){
            n_valid += 1;
            v_valid += v;
        }
    }

    n_score=n_valid/n_buildings;
    v_score=v_valid/v_buildings;

    pair<double, double> score = make_pair(n_score, v_score);

    return score;
}

void cells_per_cube_root_building(int BR, vector<double> f1, vector<double> f2, vector<double> f3, double th_ratio, double th_cells, bool target, double x, double y){
    pair<double, double> score = score_cells_per_cube_root_volume();
    double ratio = meshCFD.hcell/meshCFD.wcell;

    while (score.first < 1 && ratio <= th_ratio){
        meshCFD.wcell = meshCFD.wcell - 0.1;
        ratio = meshCFD.hcell/(double)meshCFD.wcell;
        computational_domain_II(BR);
        refinement_boxes(f1, f2, f3);
        max_number_of_cells(BR, f1, f2, f3, th_cells);
        score = score_cells_per_cube_root_volume();
    }

    pair<double, double> RoIScore = score_cells_per_cube_root_volume_RoI(target, x, y);
    cout << "Cells per cube root building" << endl;
    cout << "Score n entire model: " << score.first << endl;
    cout << "Score v entire model: " << score.second << endl;
    cout << "Score n RoI: " << RoIScore.first << endl;
    cout << "Score v RoI: " << RoIScore.second << endl;
    cout << "Ratio: " << ratio << endl;
    cout << "hcellxwcell: " << meshCFD.hcell << "x" << meshCFD.wcell << endl;
    cout << "\n";
}

vector<size_t> convex_hull(vector<Point2> vertices) {
    vector<size_t> indices(vertices.size()), out;
    iota(indices.begin(), indices.end(),0);
    CGAL::convex_hull_2(indices.begin(), indices.end(), back_inserter(out),Convex_hull_traits_2(CGAL::make_property_map(vertices)));
    return out;
}

pair<pair<Point2, Point2>, double> shortest_separation(vector<Point2> gsA, vector<Point2> gsB, double threshold){

    pair<pair<Point2, Point2>, double> sp;

    vector<size_t> a_out = convex_hull(gsA);
    vector<size_t> b_out = convex_hull(gsB);

    vector<double> distances;
    vector<pair<int, int>> pts;
    for (size_t ia: a_out){
        for (size_t ib: b_out){
            double d = sqrt(CGAL::squared_distance(gsA[ia], gsB[ib]));
            pair<int, int> ab = {ia, ib};
            if (d >= threshold){
                distances.push_back(d);
                pts.push_back(ab);
            }
        }
    }

    auto min = min_element(distances.begin(), distances.end());
    int minIndex = distance(distances.begin(), min);


    if (min != distances.end()) {
        double d = *min;
        pair <Point2, Point2> ab = {gsA[pts[minIndex].first], gsB[pts[minIndex].second]};
        sp = make_pair(ab, d);
    } else {
        sp.second = -1;
    }

    return sp;
}

double score_cells_per_building_separation(double threshold){

    double score;

    double wcell4 = meshCFD.wcell/(double)16;
    double wcell3 = meshCFD.wcell/(double)8;
    double min = 2*meshCFD.nCellsBetweenLevels*wcell4+2*wcell3;

    double n_buildings = buildings.size();
    double n_valid = 0;

    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {

        SMesh building = b->second.mesh;
        double v = CGAL::Polygon_mesh_processing::volume(building);

        if (v >= 0) {
            vector < pair < pair < Point2, Point2 >, double >> separation_pairs;
            vector < pair < pair < Point2, Point2 >, double >> invalid_separations;

            for (map<int, Building>::iterator bb = buildings.begin(); bb != buildings.end(); bb++) {
                if (b != bb) {
                    pair<pair<Point2, Point2>, double> smin = shortest_separation(b->second.gs2D, bb->second.gs2D,
                                                                                  threshold);
                    separation_pairs.push_back(smin);

                    if (smin.second < min && smin.second != -1) {
                        invalid_separations.push_back(smin);
                    }

                }
            }

            if (invalid_separations.empty()) {
                n_valid += 1;
                b->second.invalid_separation = 0;
            } else {
                b->second.invalid_separations = invalid_separations;
                b->second.invalid_separation = 1;
            }
        } else {
            b->second.invalid_separation = 2;
        }
    }

    score = n_valid/n_buildings;

    return score;

}

double score_cells_per_building_separation_RoI(double threshold, bool target, double x, double y){

    double score;

    vector<int> RoIBuildings = Buildings_in_RoI(target, x, y);

    double n_buildings = RoIBuildings.size();
    double n_valid = 0;

    for (vector<int>::iterator b=RoIBuildings.begin(); b!=RoIBuildings.end(); b++) {
        if (buildings[*b].invalid_separation == 0){
            n_valid += 1;
        }
    }

    score = n_valid/n_buildings;

    return score;

}


void cells_per_building_separation(int BR, vector<double> f1, vector<double> f2, vector<double> f3, double th_ratio, double th_cells, bool target, double x, double y){

    double score = score_cells_per_building_separation(th_cells);

    double ratio = meshCFD.hcell/meshCFD.wcell;

    while (score < 1 && ratio <= th_ratio){
        meshCFD.wcell = meshCFD.wcell - 0.1;
        ratio = meshCFD.hcell/(double)meshCFD.wcell;
        computational_domain_II(BR);
        refinement_boxes(f1, f2, f3);
        max_number_of_cells(BR, f1, f2, f3, th_cells);
        score = score_cells_per_building_separation(th_cells);
    }

    double RoIScore = score_cells_per_building_separation_RoI(th_cells, target, x, y);

    cout << "Cells per building separation" << endl;
    cout << "Score entire model: " << score << endl;
    cout << "Score RoI: " << RoIScore << endl;
    cout << "Ratio: " << ratio << endl;
    cout << "hcellxwcell: " << meshCFD.hcell << "x" << meshCFD.wcell << endl;
    cout << "\n";

}

void write_blockMeshDict() {

    meshCFD.nx = (meshCFD.domain.xmax()-meshCFD.domain.xmin())/(double)meshCFD.wcell;
    meshCFD.ny = (meshCFD.domain.ymax()-meshCFD.domain.ymin())/(double)meshCFD.wcell;
    meshCFD.nz = (meshCFD.domain.zmax()-meshCFD.domain.zmin())/(double)meshCFD.hcell;

//    cout << "nx " << meshCFD.nx << endl;
//    cout << "ny " << meshCFD.ny << endl;
//    cout << "nz " << meshCFD.nz << endl;

    ifstream input_stream;
    input_stream.open("../data/blockMeshDict");

    ofstream ofile("../output/blockMeshDict");

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
    Point3 p (meshCFD.box1.xmin(), meshCFD.box1.ymin(), meshCFD.box1.zmin());
    return p;
}

bool evaluation_height_II(double eh_target){
    bool valid;
    double eh = meshCFD.hcell*2.5;

    if (eh > eh_target){
        valid = 0;
    } else {
        valid = 1;
    }

    return valid;
}

pair<int, double> evaluation_height_III(double eh_target){

    pair<int, double> extra_refinement;

    double cell0 = meshCFD.hcell;
    double cell1 = meshCFD.hcell/2;
    double cell2 = meshCFD.hcell/4;
    double cell3 = meshCFD.hcell/8;
    double cell4 = meshCFD.hcell/16;

    vector<double> cells {cell0, cell1, cell2, cell3, cell4};

    double eh0=cell0*2.5;
    double eh1=cell1*2.5;
    double eh2=cell2*2.5;
    double eh3=cell3*2.5;
    double eh4=cell4*2.5;

    cout << "eh level 0 " << eh0 << endl;
    cout << "eh level 1 " << eh1 << endl;
    cout << "eh level 2 " << eh2 << endl;
    cout << "eh level 3 " << eh3 << endl;
    cout << "eh level 4 " << eh4 << endl;

    vector<double> ehs {eh0, eh1, eh2, eh3, eh4};

    double di = abs(ehs[0]-eh_target);
    double level = eh0;
    double idx = 0;

    for (int e=0; e<ehs.size(); e++){
        double d = abs(ehs[e]-eh_target);
        if (d < di){
            di = d;
            level = ehs[e];
            idx = e;
        }
    }

    double h0 = cell0*3;
    double h1 = cell1*3;
    double h2 = cell2*3;
    double h3 = cell3*3;
    double h4 = cell4*3;

    vector<double> hs {h0, h1, h2, h3, h4};

    double h = hs[idx];
    cout << "h0 " << h0 << endl;
    cout << "h1 " << h1 << endl;
    cout << "h2 " << h2 << endl;
    cout << "h3 " << h3 << endl;
    cout << "h4 " << h4 << endl;

    cout << "level " << level << endl;
    cout << "idx " << idx << endl;
    cout << "h " << h << endl;

    pair<double, double> zComp = fit_to_cell_dimensions_RB(city.zmin, h, cells[idx], 1);
    h = zComp.second;

    cout << "h " << h << endl;

    extra_refinement = make_pair(idx, h);

    return extra_refinement;
}

void write_snappyHexMeshDict(string filename, double eh_target){

    ifstream input_stream;
    input_stream.open("../data/snappyHexMeshDict");
    ofstream ofile("../output/snappyHexMeshDict");

    int idx = filename.find_last_of(".");
    string name = filename.substr(0, idx);

    vector<CGAL::Bbox_3> domains {meshCFD.domain, meshCFD.box3, meshCFD.box2, meshCFD.box1};

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

            else if (data == "input6") {
                ofile << "        min  (" << meshCFD.box3.xmin() << " " << meshCFD.box3.ymin() << " " << meshCFD.box3.zmin() << ");" << endl;
            }

            else if(data == "input7") {
                ofile << "        max  (" << meshCFD.box3.xmax() << " " << meshCFD.box3.ymax() << " " << meshCFD.box3.zmax() << ");" << endl;
            }


            else if(data == "input8") {
                ofile << "        { file " << "\"" << name << ".eMesh" << "\"; level 1; }" << endl;
            }

            else if(data == "input9") {
                Point3 p = inside_point();
                ofile << "    locationInMesh (" << p.x() << " " << p.y() << " " << p.z()<<");"<< endl;
            }

            else if(data == "input10") {
                if (evaluation_height_II(eh_target) == 0){
                    pair<int, double> ref = evaluation_height_III(eh_target);
                    if (ref.first != 4){
                        CGAL::Bbox_3 ground = domains[ref.first];
                        double h = ref.second;

                        ofile << " " << endl;
                        ofile << "    ground1" << endl;
                        ofile << "    {" << endl;
                        ofile << "        type searchableBox;" << endl;
                        ofile << "        min  (" << domains[0].xmin() << " " << domains[0].ymin() << " " << domains[0].zmin() << ");" << endl;
                        ofile << "        max  (" << ground.xmin() << " " << domains[0].ymax() << " " << h << ");" << endl;
                        ofile << "    }" << endl;
                        ofile << " " << endl;

                        ofile << "    ground2" << endl;
                        ofile << "    {" << endl;
                        ofile << "        type searchableBox;" << endl;
                        ofile << "        min  (" << ground.xmin() << " " << ground.ymax() << " " << domains[0].zmin() << ");" << endl;
                        ofile << "        max  (" << ground.xmax() << " " << domains[0].ymax() << " " << h << ");" << endl;
                        ofile << "    }" << endl;
                        ofile << " " << endl;

                        ofile << "    ground3" << endl;
                        ofile << "    {" << endl;
                        ofile << "        type searchableBox;" << endl;
                        ofile << "        min  (" << ground.xmin() << " " << domains[0].ymin() << " " << domains[0].zmin() << ");" << endl;
                        ofile << "        max  (" << ground.xmax() << " " << ground.ymin() << " " << h << ");" << endl;
                        ofile << "    }" << endl;
                        ofile << " " << endl;

                        ofile << "    ground4" << endl;
                        ofile << "    {" << endl;
                        ofile << "        type searchableBox;" << endl;
                        ofile << "        min  (" << ground.xmax() << " " << domains[0].ymin() << " " << domains[0].zmin() << ");" << endl;
                        ofile << "        max  (" << domains[0].xmax() << " " << domains[0].ymax() << " " << h << ");" << endl;
                        ofile << "    }" << endl;
                        ofile << " " << endl;


                    } else if (ref.first == 4){
                        double h = ref.second;
                        ofile << " " << endl;
                        ofile << "    ground" << endl;
                        ofile << "    {" << endl;
                        ofile << "        type searchableBox;" << endl;
                        ofile << "        min  (" << domains[0].xmin() << " " << domains[0].ymin() << " " << domains[0].zmin() << ");" << endl;
                        ofile << "        max  (" << domains[0].xmax() << " " << domains[0].ymax() << " " << h << ");" << endl;
                        ofile << "    }" << endl;
                        ofile << " " << endl;
                    }

                } else {
                    ofile << " " << endl;
                }

            }

            else if(data == "input11") {
                if (evaluation_height_II(eh_target) == 0){
                    pair<int, double> ref = evaluation_height_III(eh_target);
                    int level = ref.first;

                    if (level != 4){
                        ofile << " " << endl;
                        ofile << "        ground1" << endl;
                        ofile << "        {" << endl;
                        ofile << "            mode inside;" << endl;
                        ofile << "            levels ((1E15 " << level << "));" << endl;
                        ofile << "        }" << endl;
                        ofile << " " << endl;

                        ofile << "        ground2" << endl;
                        ofile << "        {" << endl;
                        ofile << "            mode inside;" << endl;
                        ofile << "            levels ((1E15 " << level << "));" << endl;
                        ofile << "        }" << endl;
                        ofile << " " << endl;

                        ofile << "        ground3" << endl;
                        ofile << "        {" << endl;
                        ofile << "            mode inside;" << endl;
                        ofile << "            levels ((1E15 " << level << "));" << endl;
                        ofile << "        }" << endl;
                        ofile << " " << endl;

                        ofile << "        ground4" << endl;
                        ofile << "        {" << endl;
                        ofile << "            mode inside;" << endl;
                        ofile << "            levels ((1E15 " << level << "));" << endl;
                        ofile << "        }" << endl;
                        ofile << " " << endl;
                    } else if (level == 4){
                        ofile << " " << endl;
                        ofile << "        ground" << endl;
                        ofile << "        {" << endl;
                        ofile << "            mode inside;" << endl;
                        ofile << "            levels ((1E15 4));" << endl;
                        ofile << "        }" << endl;
                        ofile << " " << endl;
                    }

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



void meshing(){
    //User parameters
    meshCFD.nCellsBetweenLevels = 4; //default, add as a user parameter
    meshCFD.z0 = 0.2;
    double th_cells = 700000; //default, add as a user parameter
    double th_separation = 2; //default, add as a user parameter
    double th_ratio = 1.2; //default, add as a user parameter
    int BR = 0; //default, add as a user parameter
    //default, add as a user parameter
    int target = 0;
    double x = 0;
    double y = 0;
    double eh_target = 1.5;
    //vector<double> f1{0, 0, 0, 1.5}; //inlet, outlet, lateral, top
//    vector<double> f1{0.10, 0.10, 0.10, 0.10};
//    vector<double> f2{0.60, 2, 0.60, 0.60};
//    vector<double> f3{1, 3, 1, 1.5};

    vector<double> f1{0, 0, 0, 1.5};
    vector<double> f2{1, 2, 1, 2};
    vector<double> f3{2, 5, 2, 2.5};

    evaluation_height(BR, f1, f2, f3);
    roughness(BR, f1, f2, f3);
    max_number_of_cells(BR, f1, f2, f3, th_cells);
    cells_per_cube_root_building(BR, f1, f2, f3, th_ratio, th_cells, target, x, y);
    if(buildings.size()>1){
        cells_per_building_separation(BR, f1, f2, f3, th_ratio, th_separation, target, x, y);
    }

    cout << "Final cell dimensions:" << endl;
    cout << "hcellxwcell: " << meshCFD.hcell << "x" << meshCFD.wcell << endl;

    double evaluation_height = meshCFD.hcell*2.5;
    cout << "Evaluation height: " << evaluation_height << endl;

    double roughness_check = (meshCFD.hcell/(double)16)/(double)2;
    cout << "Roughness check: z0 is " << meshCFD.z0 << ", hcell level 4 divided by 2 is " << roughness_check << endl;

    double n = n_cells();
    cout << "Number of cells (threshold: " << th_cells << "): " << n << endl;

    double n_score_volume = score_cells_per_cube_root_volume().first*100;
    double v_score_volume = score_cells_per_cube_root_volume().second*100;
    cout << "Number of buildings having 10 cells per cube root volume: " << n_score_volume << "%" << endl;
    cout << "Building volume having 10 cells per cube root volume: " << v_score_volume << "%" << endl;

    double score_separation = score_cells_per_building_separation(th_separation)*100;
    cout << "Number of buildings having 10 cells per building separation: " << score_separation << "%" << endl;




}



double maximum_dimension(Building b){
    vector<double> x_pts;
    vector<double> y_pts;
    vector<double> z_pts;

    for (auto &v: b.pts){
        double x = v.x();
        double y = v.y();
        double z = v.z();
        x_pts.push_back(x);
        y_pts.push_back(y);
        z_pts.push_back(z);
    }

    auto x_minmax = minmax_element(x_pts.begin(), x_pts.end(),[](const double& p1, const double& p2) { return p1 < p2; });
    auto y_minmax = minmax_element(y_pts.begin(), y_pts.end(),[](const double& p1, const double& p2) { return p1 < p2; });
    auto z_minmax = minmax_element(z_pts.begin(), z_pts.end(),[](const double& p1, const double& p2) { return p1 < p2; });

    double length = *x_minmax.second-*x_minmax.first;
    double width = *y_minmax.second-*y_minmax.first;
    double height = *z_minmax.second-*z_minmax.first;
    vector<double> max_dimensions {length, width, height};
    auto dim_max = max_element(max_dimensions.begin(), max_dimensions.end());

    double max = *dim_max;

    return max;
}

SMesh RoI_cylinder(Point2 center, double radius){
    SMesh m;

    double z0 = meshCFD.refinementBox1.zmin();
    double z1 = meshCFD.refinementBox1.zmax();

    int n_circle = radius*0.5;
    double interval = 2.0*M_PI/n_circle;

    vector<Point2> pCircle;

    for (int n = 0; n < n_circle; ++n) {
        double angle = n*interval;
        double x = center.x() + radius * cos(angle);
        double y = center.y() + radius * sin(angle);
        pCircle.emplace_back(Point2(x, y));
    }

    vector<vertex_descriptor> pBottom;
    vector<vertex_descriptor> pTop;

    vertex_descriptor b0 = m.add_vertex(Point3(center.x(), center.y(), z0));
    vertex_descriptor t0 = m.add_vertex(Point3(center.x(), center.y(), z1));

    for (auto p: pCircle){
        pBottom.push_back(m.add_vertex(Point3(p.x(), p.y(), z0)));
        pTop.push_back(m.add_vertex(Point3(p.x(), p.y(), z1)));
    }

    m.add_face(pBottom[pBottom.size()-1], b0, pBottom[0]);
    for (int k = 0; k < pBottom.size()-1; k++){
        m.add_face(pBottom[k], b0, pBottom[k+1]);
    }

    m.add_face(pTop[pTop.size()-1], pTop[0], t0);
    for (int k = 0; k < pTop.size()-1; ++k){
        m.add_face(pTop[k], pTop[k+1], t0);
    }

    int n = pBottom.size();
    m.add_face(pTop[0], pTop[n-1], pBottom[n-1]);
    m.add_face(pBottom[n-1], pBottom[0], pTop[0]);
    for (int i=0; i<n-1; i++){
        m.add_face(pTop[i+1], pTop[i], pBottom[i]);
        m.add_face(pBottom[i], pBottom[i+1], pTop[i+1]);
    }

    return m;
}

void Liu_RoI(){
    double px = (city.xmin+city.xmax)*0.5;
    double py = (city.ymin+city.ymax)*0.5;
    Point2 p(px, py);

    vector<double> max_values;
    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {
        double max = maximum_dimension(b->second);
        max_values.push_back(max);
    }

    auto dim_max = max_element(max_values.begin(), max_values.end());
    double lmax = *dim_max;

    double r = 3*lmax;
    meshCFD.RoICircle = Circle2(p, pow(r, 2));

    SMesh RoI = RoI_cylinder(p, r);
    meshCFD.RoI = RoI;
}

void Liu_RoI(double x, double y){
    double px = x;
    double py = y;
    Point2 p(px, py);

    int id;

    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {

        vector <size_t> a_out = convex_hull(b->second.gs2D);
        vector <Point2> boundary;
        for (size_t ia: a_out) {
            boundary.push_back(b->second.gs2D[ia]);
        }

        switch (CGAL::bounded_side_2(boundary.begin(), boundary.end(), p, K())) {
            case CGAL::ON_BOUNDED_SIDE:
                id = b->first;
                break;
            case CGAL::ON_BOUNDARY:
                id = b->first;
                break;
        }
    }

    double lmax = maximum_dimension(buildings[id]);
    double r = 3*lmax;
    meshCFD.RoICircle = Circle2(p, pow(r, 2));

    SMesh RoI = RoI_cylinder(p, r);
    meshCFD.RoI = RoI;

}

bool incircle(Circle2& circle, Point2& p){
    Point2 center = circle.center();
    FT squared_radius = circle.squared_radius();

    switch(CGAL::compare(CGAL::square(p.x()-center.x())-squared_radius, -CGAL::square(p.y()-center.y())) ){
        case CGAL::LARGER:
            return false;
        case CGAL::SMALLER:
            return true;
        case CGAL::EQUAL:
            return true;
    }
}

vector<int> Buildings_in_RoI(bool target, double x, double y) {

    vector<int> detailed_buildings;

    SMesh RoI;
    if (target == false){
        Liu_RoI();
    } else{
        Liu_RoI(x, y);
    }

    Circle2 circle = meshCFD.RoICircle;

    int count = 0;
    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {
        SMesh building = b->second.mesh;
        SMesh out;

        bool valid_intersection = PMP::do_intersect(meshCFD.RoI, building);
        if (valid_intersection) {
            detailed_buildings.push_back(b->second.id);
        } else {
            bool in = 0;
            if (!b->second.gs2D.empty()){
                for (auto &p: b->second.gs2D){
                    if (!incircle(circle, p)) {
                        in = 1;
                    }
                }

                if (in == 0) {
                    count += 1;
                    detailed_buildings.push_back(b->second.id);
                }
            }
        }
    }

    return detailed_buildings;
}