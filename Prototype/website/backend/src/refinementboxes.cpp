#include "definitions.h"
#include "input.h"
#include "meshtools.h"

pair<double, double> fit_RB_to_cell_dimensions(double min, double max, double dcell, int xyz, CGAL::Bbox_3 min_box) {

    //Same algorithm as fit_domain_to_cell_dimensions in domain.cpp,
    //but then with refinement boxes

    double check = abs(max-min)/(double)dcell;

    if (check - trunc(check) > 0){

        min = round5(min);
        max = round5(max);

        double L=0;
        if (xyz == 0){
            L = abs(min_box.xmax()-min_box.xmin());
        } else if (xyz == 1){
            L = abs(min_box.ymax()-min_box.ymin());
        } else if (xyz == 2){
            L = abs(min_box.zmax()-min_box.zmin());
        }

        double l = abs(max-min);
        int n = trunc(l/(double)dcell);
        double k = round5(dcell*n);
        while (k < L){
            k += dcell;
        }

        double c = (min+max)*0.5;
        if (xyz == 2){
            min = min;
            max = min+k;
        } else {
            min = c-k*0.5;
            max = c+k*0.5;
        }
    }

    return pair(min, max);
}

CGAL::Bbox_3 min_refinement_box(double finlet, double foutlet, double flateral, double ftop) {

    //defines minimum dimensions of a refinement box

    double xmin = city.xmin-finlet*city.dmax;
    double xmax = city.xmax+foutlet*city.dmax;
    double ymin = city.ymin-flateral*city.dmax;
    double ymax = city.ymax+flateral*city.dmax;

    double zmin;
    if (fTerrain.size() != 0){ //similar to the domain, the minimum z-value of terrain faces are used as a minimum z-value for the refinement box
        zmin = city.zmin_terrain;
    } else { //otherwise, the lowest z-value of the urban area is applied
        zmin = city.zmin;
    }

    double zmax = city.zmax+city.dmax*ftop;

    CGAL::Bbox_3 bbox(xmin, ymin, zmin, xmax, ymax, zmax);

    return bbox;
}

CGAL::Bbox_3 refinement_box(double finlet, double foutlet, double flateral, double ftop, CGAL::Bbox_3 min_box) {

    //defines a refinement box

    double hcell = meshCFD.hcell;
    double wcell = meshCFD.wcell;

    double xmin = city.xmin-finlet*city.dmax;
    double xmax = city.xmax+foutlet*city.dmax;
    double ymin = city.ymin-flateral*city.dmax;
    double ymax = city.ymax+flateral*city.dmax;

    double zmin;
    if (fTerrain.size() != 0){
        zmin = city.zmin_terrain;
    } else {
        zmin = city.zmin;
    }

    double zmax = city.zmax+city.dmax*ftop;

    pair<double, double> xComp = fit_RB_to_cell_dimensions(xmin, xmax, wcell, 0, min_box);
    pair<double, double> yComp = fit_RB_to_cell_dimensions(ymin, ymax, wcell, 1, min_box);
    pair<double, double> zComp = fit_RB_to_cell_dimensions(zmin, zmax, hcell, 2, min_box);

    xmin = xComp.first;
    xmax = xComp.second;
    ymin = yComp.first;
    ymax = yComp.second;
    zmin = zComp.first;
    zmax = zComp.second;

    CGAL::Bbox_3 bbox(xmin, ymin, zmin, xmax, ymax, zmax);

    return bbox;
}


void B2_refinement_boxes(vector<double> f1, vector<double> f2){

    //defines two refinement boxes

    CGAL::Bbox_3 min_box1 = min_refinement_box(f1[0], f1[1], f1[2], f1[3]);
    CGAL::Bbox_3 min_box2 = min_refinement_box(f2[0], f2[1], f2[2], f2[3]);

    CGAL::Bbox_3 box1 = refinement_box(f1[0], f1[1], f1[2], f1[3], min_box1);
    CGAL::Bbox_3 box2 = refinement_box(f2[0], f2[1], f2[2], f2[3], min_box2);

    meshCFD.box1 = box1;
    meshCFD.box2 = box2;
}


void B3_refinement_boxes(vector<double> f1, vector<double> f2, vector<double> f3){

    //defines three refinement boxes

    CGAL::Bbox_3 min_box1 = min_refinement_box(f1[0], f1[1], f1[2], f1[3]);
    CGAL::Bbox_3 min_box2 = min_refinement_box(f2[0], f2[1], f2[2], f2[3]);
    CGAL::Bbox_3 min_box3 = min_refinement_box(f3[0], f3[1], f3[2], f3[3]);

    CGAL::Bbox_3 box1 = refinement_box(f1[0], f1[1], f1[2], f1[3], min_box1);
    CGAL::Bbox_3 box2 = refinement_box(f2[0], f2[1], f2[2], f2[3], min_box2);
    CGAL::Bbox_3 box3 = refinement_box(f3[0], f3[1], f3[2], f3[3], min_box3);

    meshCFD.box1 = box1;
    meshCFD.box2 = box2;
    meshCFD.box3 = box3;
}