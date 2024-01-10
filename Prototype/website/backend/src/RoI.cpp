#include "definitions.h"
#include "input.h"
#include "domain.h"
#include "meshtools.h"
#include "buildings.h"

SMesh RoI_cylinder(Point2 center, double radius){

    //creates a cylinder representing the Region of Interest

    SMesh m;

    double z0 = meshCFD.box1.zmin();
    double z1 = meshCFD.box1.zmax();

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

    //defines Region of Interest (RoI) without target building
    //the centre of the input is used
    double px = (city.xmin+city.xmax)*0.5;
    double py = (city.ymin+city.ymax)*0.5;
    Point2 p(px, py);

    //find maximum dimension of buildings
    vector<double> max_values;
    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {
        double max = maximum_dimension(b->second);
        max_values.push_back(max);
    }

    auto dim_max = max_element(max_values.begin(), max_values.end());
    double lmax = *dim_max;

    //create circle and cylinder representing RoI
    double r = 3*lmax;
    meshCFD.RoICircle = Circle2(p, pow(r, 2));

    SMesh RoI = RoI_cylinder(p, r);
    meshCFD.RoI = RoI;
}

void Liu_RoI(double x, double y){

    //defines Region of Interest (RoI) with target building
    //the x and y values defined by users are used as centre of this RoI
    double px = x;
    double py = y;
    Point2 p(px, py);

    //identify the corresponding building of these x and y values
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

    //compute the largest dimension of this building
    double lmax = maximum_dimension(buildings[id]);
    double r = 3*lmax;

    //create circle and cylinder representing RoI
    meshCFD.RoICircle = Circle2(p, pow(r, 2));
    SMesh RoI = RoI_cylinder(p, r);

    meshCFD.RoI = RoI;
}

bool incircle(Circle2& circle, Point2& p){

    //checks if a point is located in circle

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

    //identifies buildings in RoI

    vector<int> detailed_buildings;

    //depending on if there is a target building selected, the RoI is defined
    SMesh RoI;
    if (target == false){
        Liu_RoI();
    } else{
        Liu_RoI(x, y);
    }

    Circle2 circle = meshCFD.RoICircle;

    //identifies which buildings are in this RoI,
    //by using both the cylinder (boundaries area) and circle (interior area)
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

    //return buildings situated in RoI
    return detailed_buildings;
}



