#include "definitions.h"
#include "input.h"
#include "buildings.h"
#include "meshtools.h"

pair<double, double> fit_domain_to_cell_dimensions(double min, double max, double dcell, int xyz) {

    //ensures that a number of cells (integer) can fit within domain in one direction
    // xyz = 0  -> x-direction
    // xyz = 1 -> y-direction
    // xyz = 2 -> z-direction

    double check = abs(max-min)/(double)dcell; //check if a whole number of cells fits within domain dimension in x, y, or z direction

    if (check - trunc(check) > 0){
        //if not, round min and max values to 5 decimals to avoid mathematical errors
        //the many computations may lead to values with very small decimals (e.g. 0.000000001), which are better to neglect for this purpose
        min = round5(min);
        max = round5(max);

        double L =0; //compute the minimum length in the requested direction L
        if (xyz == 0){
            L = abs(meshCFD.min_domain.xmax()-meshCFD.min_domain.xmin());
        } else if (xyz == 1){
            L = abs(meshCFD.min_domain.ymax()-meshCFD.min_domain.ymin());
        } else if (xyz == 2){
            L = abs(meshCFD.min_domain.zmax()-meshCFD.min_domain.zmin());
        }

        //compute the number of cells that fits within length in the requested direction l
        //round this number to the highest integer and multiply it with the cell size
        //new l is not allowed to be smaller than L, to avoid domain dimensions that are too small
        //if l is initially too small, extra cells are added until it is greater than L
        double l = abs(max-min);
        int n = trunc(l/(double)dcell);
        double k = round5(dcell*n);
        while (k < L){
            k += dcell;
        }

        //adjust max and min values in the required direction
        double c = (min+max)*0.5;
        if (xyz == 2){
            min = min;
            max = min+k;
        } else {
            min = c-k*0.5;
            max = c+k*0.5;
        }
    }

    //return new domain dimensions in the required direction
    return pair(min, max);
}


void get_city_boundaries() {

    //finds max and min values of the urban area (only buildings) in each direction
    //and stores them in "city"
    auto x_minmax = minmax_element(city.vb.begin(), city.vb.end(), [](const Point3 &p1, const Point3 &p2) { return p1.x() < p2.x(); });
    auto y_minmax = minmax_element(city.vb.begin(), city.vb.end(), [](const Point3 &p1, const Point3 &p2) { return p1.y() < p2.y(); });
    auto z_minmax = minmax_element(city.vb.begin(), city.vb.end(), [](const Point3 &p1, const Point3 &p2) { return p1.z() < p2.z(); });

    double xmin = x_minmax.first->x();
    double xmax = x_minmax.second->x();
    double ymin = y_minmax.first->y();
    double ymax = y_minmax.second->y();
    double zmin = z_minmax.first->z();
    double zmax = z_minmax.second->z();

    city.xmin = xmin;
    city.xmax = xmax;
    city.ymin = ymin;
    city.ymax = ymax;
    city.zmin = zmin;
    city.zmax = zmax;
}

void highest_building() {

    //computes height of the highest building
    vector<double> lengths;
    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++){
        double length = building_height(b->second);
        lengths.push_back(length);
    }

    auto max = max_element(lengths.begin(), lengths.end());
    double hmax = *max;

    city.dmax = hmax;
}

void largest_building() {

    //computes the maximum building dimension of the input model

    vector<double> max_dimension;

    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++){

        vector<double> x_pts;
        vector<double> y_pts;
        vector<double> z_pts;

        for (auto &v: b->second.pts){
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

        max_dimension.push_back(max);
    }

    auto max = max_element(max_dimension.begin(), max_dimension.end());
    double dmax = *max;

    city.dmax = dmax;
}

void BR_only_option() {

    //adjusts domain dimensions based on blockage ratio defined by Franke and Baklanov (2007)

    double xmin_i = meshCFD.domain.xmin();
    double xmax_i = meshCFD.domain.xmax();
    double ymin_i = meshCFD.domain.ymin();
    double ymax_i = meshCFD.domain.ymax();
    double zmin_i = meshCFD.domain.zmin();
    double zmax_i = meshCFD.domain.zmax();

    double xmin_final, xmax_final;
    double ymin_final, ymax_final;
    double zmin_final, zmax_final;

    double ymin_BR = ymin_i;
    double ymax_BR = ymax_i;
    double zmax_BR = zmax_i;

    //compute blockage ratio BR
    double A_urban = abs(city.ymax - city.ymin) * abs(city.zmax - city.zmin);
    double A_domain = abs(ymax_BR - ymin_BR) * abs(zmax_BR - zmin_i);
    double BR = A_urban / A_domain;

    //so long as BR is higher than 0.03
    //domain dimensions are increased in the y and z directions
    while (BR > 0.03) {
        ymin_BR = ymin_BR - 0.01;
        ymax_BR = ymax_BR + 0.01;
        zmax_BR = zmax_BR + 0.01;
        A_domain = abs(ymax_BR- ymin_BR) * abs(zmax_BR - zmin_i);
        BR = A_urban / A_domain;

        xmin_final = xmin_i;
        xmax_final = xmax_i;
        ymin_final = ymin_BR;
        ymax_final = ymax_BR;
        zmin_final = zmin_i;
        zmax_final = zmax_BR;
    }

    //ensure that a whole number of cells fits within the new domain
    pair<double, double> xComp = fit_domain_to_cell_dimensions(xmin_final, xmax_final, meshCFD.wcell, 0);
    pair<double, double> yComp = fit_domain_to_cell_dimensions(ymin_final, ymax_final, meshCFD.wcell, 1);
    pair<double, double> zComp = fit_domain_to_cell_dimensions(zmin_final, zmax_final, meshCFD.hcell, 2);

    xmin_final = xComp.first;
    xmax_final = xComp.second;
    ymin_final = yComp.first;
    ymax_final = yComp.second;
    zmin_final = zComp.first;
    zmax_final = zComp.second;

    CGAL::Bbox_3 domain = CGAL::Bbox_3(xmin_final, ymin_final, zmin_final, xmax_final, ymax_final, zmax_final);

    //return the new domain dimensions
    meshCFD.domain = domain;
}

void BRL_BRH_option() {

    //adjusts domain dimensions based on blockage ratios defined by Blocken (2015)

    //compute blockage ratio in the lateral horizontal direction BRL
    double L_urban = abs(city.ymax - city.ymin);
    double L_domain = abs(meshCFD.domain.ymax()-meshCFD.domain.ymin());
    double BR_L = L_urban/L_domain;

    //compute blockage ratio in the lateral vertical direction BRH
    double H_urban = abs(city.zmax-city.zmin);
    double H_domain = abs(meshCFD.domain.zmax()-meshCFD.domain.zmin());
    double BR_H = H_urban/H_domain;

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

    //so long as BRL is higher than 0.17, domain dimensions are increased in the y-direction
    while (BR_L>0.17) {
        ymin_BRh = ymin_BRh-0.01;
        ymax_BRh = ymax_BRh+0.01;
        L_domain = abs(ymax_BRh-ymin_BRh);
        BR_L = (double)L_urban/(double)L_domain;
    }

    double zmax_BRv = zmax_i;

    //so long as BRL is higher than 0.17, domain dimensions are increased in the z-direction
    while (BR_H > 0.17) {
        zmax_BRv = zmax_BRv + 0.01;
        H_domain = abs(zmax_BRv - zmin_i);
        BR_H = H_urban / H_domain;
    }

    xmin_final = xmin_i;
    xmax_final = xmax_i;
    ymin_final = ymin_BRh;
    ymax_final = ymax_BRh;
    zmin_final = zmin_i;
    zmax_final = zmax_BRv;

    //ensure that a whole number of cells fits within the new domain
    pair<double, double> xComp = fit_domain_to_cell_dimensions(xmin_final, xmax_final, meshCFD.wcell, 0);
    pair<double, double> yComp = fit_domain_to_cell_dimensions(ymin_final, ymax_final, meshCFD.wcell, 1);
    pair<double, double> zComp = fit_domain_to_cell_dimensions(zmin_final, zmax_final, meshCFD.hcell, 2);

    xmin_final = xComp.first;
    xmax_final = xComp.second;
    ymin_final = yComp.first;
    ymax_final = yComp.second;
    zmin_final = zComp.first;
    zmax_final = zComp.second;

    CGAL::Bbox_3 domain = CGAL::Bbox_3(xmin_final, ymin_final, zmin_final, xmax_final, ymax_final, zmax_final);

    //return the new domain dimensions
    meshCFD.domain = domain;
}


void computational_domain_I(int BR, int dim){

    //computes initial domain dimensions (used only once at the beginning)
    //final domain cannot be smaller than the domain defined in this function

    get_city_boundaries();

    if (dim == 0){
        highest_building(); //height of the highest is used
    } else {
        largest_building(); //largest building dimension is used
    }

    double hcell = meshCFD.hcell;
    double wcell = meshCFD.wcell;

    //compute domain dimensions
    double dmax = city.dmax;
    double xmin = city.xmin - 5*dmax;
    double xmax = city.xmax + 15*dmax;
    double ymin = city.ymin - 5*dmax;
    double ymax = city.ymax + 5*dmax;
    double zmin;

    if (fTerrain.size() != 0){
        //if there is a terrain, the lowest z-value of the domain is the lowest one of the terrain surfaces
        zmin = city.zmin_terrain;
    } else {
        //otherwise, it is the lowest z-value of the buildings
        zmin = city.zmin;
    }

    double zmax = city.zmax + 5*dmax;

    CGAL::Bbox_3 min_domain(xmin, ymin, zmin, xmax, ymax, zmax);
    meshCFD.min_domain = min_domain;

    //ensure that a whole number of cells fits within the new domain
    pair<double, double> xComp = fit_domain_to_cell_dimensions(xmin, xmax, wcell, 0);
    pair<double, double> yComp = fit_domain_to_cell_dimensions(ymin, ymax, wcell, 1);
    pair<double, double> zComp = fit_domain_to_cell_dimensions(zmin, zmax, hcell, 2);

    xmin = xComp.first;
    xmax = xComp.second;
    ymin = yComp.first;
    ymax = yComp.second;
    zmin = zComp.first;
    zmax = zComp.second;

    CGAL::Bbox_3 domain(xmin, ymin, zmin, xmax, ymax, zmax);
    meshCFD.domain = domain;

    //dependent on the chosen blockage ratios, the domain dimensions are adjusted
    //BR = 0 -> no blockage ratios
    //BR = 1 -> blockage ratio of Franke and Baklanov (2007)
    //BR = 2 -> blockage ratios of Blocken (2015)
    if (BR == 1){
        BR_only_option();
    } else if (BR == 2){
        BRL_BRH_option();
    }
}

void computational_domain_II(int BR){

    //same function as the previous one, except that the minimum domain dimensions are not defined
    //this one is always used after the first domain definition
    double hcell = meshCFD.hcell;
    double wcell = meshCFD.wcell;

    double xmin = meshCFD.domain.xmin();
    double xmax = meshCFD.domain.xmax();
    double ymin = meshCFD.domain.ymin();
    double ymax = meshCFD.domain.ymax();
    double zmin = meshCFD.domain.zmin();
    double zmax = meshCFD.domain.zmax();

    pair<double, double> xComp = fit_domain_to_cell_dimensions(xmin, xmax, wcell, 0);
    pair<double, double> yComp = fit_domain_to_cell_dimensions(ymin, ymax, wcell, 1);
    pair<double, double> zComp = fit_domain_to_cell_dimensions(zmin, zmax, hcell, 2);

    xmin = xComp.first;
    xmax = xComp.second;
    ymin = yComp.first;
    ymax = yComp.second;
    zmin = zComp.first;
    zmax = zComp.second;

    CGAL::Bbox_3 domain(xmin, ymin, zmin, xmax, ymax, zmax);

    meshCFD.domain = domain;

    if (BR == 1){
        BR_only_option();
    } else if (BR == 2){
        BRL_BRH_option();
    }
}
