//#include "definitions.h"
//#include "input.h"
//#include "output.h"
//#include "buildings.h"
//
//
//double round5(double d){
//    d = round(d*100000)/100000;
//    return d;
//}
//
////https://www.geeksforgeeks.org/find-number-closest-n-divisible-m/
//pair<double, double> fit2cellsize(double min, double max, double wCell, bool z) {
//
//    double check = (max-min)/wCell;
//
//    if (check - trunc(check) > 0){
//
//        min = round5(min);
//        max = round5(max);
//
//        double l = (max-min);
//        int n = trunc(l/(double)wCell);
//        //cout << "n without trunc " << l/(double)wCell << endl;
//        //cout << "n " << n << endl;
//        double k = round5(wCell*n);
//
//        double c = (min+max)*0.5;
//
//        if (z == 1){
//            min = min;
//            max = min+k;
//        } else {
//            min = c-k*0.5;
//            max = c+k*0.5;
//        }
//    }
//
//    return pair(min, max);
//}
//
//void get_city_boundaries() {
//    //cout << "city size " << city.vb.size() << endl;
//    auto x_minmax = minmax_element(city.vb.begin(), city.vb.end(),[](const Point3& p1, const Point3& p2) { return p1.x() < p2.x(); });
//    auto y_minmax = minmax_element(city.vb.begin(), city.vb.end(),[](const Point3& p1, const Point3& p2) { return p1.y() < p2.y(); });
//    auto z_minmax = minmax_element(city.vb.begin(), city.vb.end(),[](const Point3& p1, const Point3& p2) { return p1.z() < p2.z(); });
//
//    city.xmin = x_minmax.first->x();
//    city.xmax = x_minmax.second->x();
//    city.ymin = y_minmax.first->y();
//    city.ymax = y_minmax.second->y();
////    cout << "city.ymin " << city.ymin << endl;
////    cout << "city.ymax " << city.ymax << endl;
//    city.zmin = z_minmax.first->z();
//    city.zmax = z_minmax.second->z();
//}
//
//double height(Building b){
//    vector<double> z_pts;
//    for (auto &v: b.pts){
//        double z = v.z();
//        z_pts.push_back(z);
//    }
//
//    auto z_minmax = minmax_element(z_pts.begin(), z_pts.end(),[](const double& p1, const double& p2) { return p1 < p2; });
//    double height = *z_minmax.second-*z_minmax.first;
//
//    return height;
//}
//
//void buildings_max_height() {
//
//    vector<double> lengths;
//    for (unordered_map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++){
//        double length = height(b->second);
//        lengths.push_back(length);
//    }
//
//    auto max = max_element(lengths.begin(), lengths.end());
//    double Hmax = *max;
//
//    city.hmax = Hmax;
//}
//
//CGAL::Bbox_3 computational_domain(double wCell){
//    get_city_boundaries();
//    buildings_max_height();
//
//    double Hmax = city.hmax;
//
//    double xmin = city.xmin - 5*Hmax;
//    double xmax = city.xmax + 15*Hmax;
//    double ymin = city.ymin - 5*Hmax;
//    double ymax = city.ymax + 5*Hmax;
//
//    double zmin = city.zmin;
//    double zmax = city.zmax + 5*Hmax;
//
//    pair<double, double> xComp = fit2cellsize(xmin, xmax, wCell, 0);
//    pair<double, double> yComp = fit2cellsize(ymin, ymax, wCell, 0);
//    pair<double, double> zComp = fit2cellsize(zmin, zmax, wCell, 1);
//
//    xmin = xComp.first;
//    xmax = xComp.second;
//    ymin = yComp.first;
//    ymax = yComp.second;
//    zmin = zComp.first;
//    zmax = zComp.second;
//
//    CGAL::Bbox_3 domain(xmin, ymin, zmin, xmax, ymax, zmax);
//
//    meshCFD.domain = domain;
//
//    return domain;
//}
//
//
//
//CGAL::Bbox_3 Refinement_box(double wCell, double finlet, double foutlet, double flateral, double ftop) {
//
////    double dx = city.xmax-city.xmin;
////    double dy = city.ymax-city.ymin;
////
////    double new_dx = dx*XYfactor;
////    double new_dy = dy*XYfactor;
////
////    double xmin = city.xmin-(new_dx-dx)*0.5;
////    double xmax = city.xmax+(new_dx-dx)*0.5;
////    double ymin = city.ymin-(new_dy-dy)*0.5;
////    double ymax = city.ymax+(new_dy-dy)*0.5;
////    double zmin = city.zmin;
////    double zmax = city.zmax+city.hmax*Zfactor;
//
//    double ymin = city.ymin-flateral*city.hmax;
//    double ymax = city.ymax+flateral*city.hmax;
//    double xmin = city.xmin-finlet*city.hmax;
//    double xmax = city.xmax+foutlet*city.hmax;
//
//    double zmin = city.zmin;
//    double zmax = city.zmax+city.hmax*ftop;
//
//    pair<double, double> xComp = fit2cellsize(xmin, xmax, wCell, 0);
//    pair<double, double> yComp = fit2cellsize(ymin, ymax, wCell, 0);
//    pair<double, double> zComp = fit2cellsize(zmin, zmax, wCell, 1);
//
//    xmin = xComp.first;
//    xmax = xComp.second;
//    ymin = yComp.first;
//    ymax = yComp.second;
//    zmin = zComp.first;
//    zmax = zComp.second;
//
//    CGAL::Bbox_3 bbox(xmin, ymin, zmin, xmax, ymax, zmax);
//
//    return bbox;
//}
//
//void RefinementBoxes(double wCell){
////    double finlet, double foutlet, double flateral, double ftop
//
//    CGAL::Bbox_3 box1 = Refinement_box(wCell, 0, 0, 0, 1.5);
//    CGAL::Bbox_3 box2 = Refinement_box(wCell, 2, 6, 2, 3);
//    //CGAL::Bbox_3 box3 = Refinement_box(wCell, 2, 6, 2, 3);
//
//    meshCFD.refinementBox1 = box1;
//    meshCFD.refinementBox2 = box2;
//    //meshCFD.refinementBox3 = box3;
//
////    cout << "REFINEMENT BOX 1" << endl;
////    cout << box1.xmin() << endl;
////    cout << box1.xmax() << endl;
////    cout << box1.ymin() << endl;
////    cout << box1.ymax() << endl;
////    cout << box1.zmin() << endl;
////    cout << box1.zmax() << endl;
////
////    cout << "REFINEMENT BOX 2" << endl;
////    cout << box2.xmin() << endl;
////    cout << box2.xmax() << endl;
////    cout << box2.ymin() << endl;
////    cout << box2.ymax() << endl;
////    cout << box2.zmin() << endl;
////    cout << box2.zmax() << endl;
////
////    cout << "REFINEMENT BOX 3" << endl;
////    cout << box3.xmin() << endl;
////    cout << box3.xmax() << endl;
////    cout << box3.ymin() << endl;
////    cout << box3.ymax() << endl;
////    cout << box3.zmin() << endl;
////    cout << box3.zmax() << endl;
//}
//
//Point3 centroid_building(int id){
//
//    double x=0, y=0, z=0;
//    int n = buildings[id].faces.size()*3;
//    for (auto &f: buildings[id].faces){
//        for (auto &v: f.v){
//            x += vertices[v].p.x();
//            y += vertices[v].p.y();
//            z += vertices[v].p.z();
//        }
//    }
//
//    Point3 centroid((double)x/n, (double)y/n, (double)z/n);
//    return centroid;
//}
//
//bool incircle(Circle2& circle, Point2& p){
//    Point2 center = circle.center();
//    FT squared_radius = circle.squared_radius();
//
//    switch(CGAL::compare(CGAL::square(p.x()-center.x())-squared_radius, -CGAL::square(p.y()-center.y())) ){
//        case CGAL::LARGER:
//            return false;
//        case CGAL::SMALLER:
//            return true;
//        case CGAL::EQUAL:
//            return true;
//    }
//}
//
//SMesh RoI_cylinder(Point2 center, double radius){
//    SMesh m;
//
//    //Refinement grid parameters
//    double z0 = meshCFD.refinementBox1.zmin();
//    double z1 = meshCFD.refinementBox1.zmax();
//
//    //Cylinder parameters
//    int n_disc = 5;
//    int n_circle = radius*0.5;
//    double interval = 2.0*M_PI/n_circle;
//
//    //Distribute points regularly within a circle (pCircle) with center and radius (RoI) in 2D
//    vector<Point2> pCircle;
//
//    for (int n = 0; n < n_circle; ++n) {
//        double angle = n*interval;
//        double x = center.x() + radius * cos(angle);
//        double y = center.y() + radius * sin(angle);
//        pCircle.emplace_back(Point2(x, y));
//    }
//
//    //Create 2 distinct 3D datasets derived from pCircle (2D points), each with different z-values: z0 and z1
//    //Create faces between these two datasets and add to mesh m
//    vector<vertex_descriptor> pBottom;
//    vector<vertex_descriptor> pTop;
//
//    vertex_descriptor b0 = m.add_vertex(Point3(center.x(), center.y(), z0));
//    vertex_descriptor t0 = m.add_vertex(Point3(center.x(), center.y(), z1));
//
//    for (auto p: pCircle){
//        pBottom.push_back(m.add_vertex(Point3(p.x(), p.y(), z0)));
//        pTop.push_back(m.add_vertex(Point3(p.x(), p.y(), z1)));
//    }
//
//    m.add_face(pBottom[pBottom.size()-1], b0, pBottom[0]);
//    for (int k = 0; k < pBottom.size()-1; k++){
//        m.add_face(pBottom[k], b0, pBottom[k+1]);
//    }
//
//    m.add_face(pTop[pTop.size()-1], pTop[0], t0);
//    for (int k = 0; k < pTop.size()-1; ++k){
//        m.add_face(pTop[k], pTop[k+1], t0);
//    }
//
//    int n = pBottom.size();
//    m.add_face(pTop[0], pTop[n-1], pBottom[n-1]);
//    m.add_face(pBottom[n-1], pBottom[0], pTop[0]);
//    for (int i=0; i<n-1; i++){
//        m.add_face(pTop[i+1], pTop[i], pBottom[i]);
//        m.add_face(pBottom[i], pBottom[i+1], pTop[i+1]);
//    }
//
//    return m;
//}
//
//double maximum_dimension(Building b){
//    vector<double> x_pts;
//    vector<double> y_pts;
//    vector<double> z_pts;
//
//    for (auto &v: b.pts){
//        double x = v.x();
//        double y = v.y();
//        double z = v.z();
//        x_pts.push_back(x);
//        y_pts.push_back(y);
//        z_pts.push_back(z);
//    }
//
//    auto x_minmax = minmax_element(x_pts.begin(), x_pts.end(),[](const double& p1, const double& p2) { return p1 < p2; });
//    auto y_minmax = minmax_element(y_pts.begin(), y_pts.end(),[](const double& p1, const double& p2) { return p1 < p2; });
//    auto z_minmax = minmax_element(z_pts.begin(), z_pts.end(),[](const double& p1, const double& p2) { return p1 < p2; });
//
//    double length = *x_minmax.second-*x_minmax.first;
//    double width = *y_minmax.second-*y_minmax.first;
//    double height = *z_minmax.second-*z_minmax.first;
//    vector<double> max_dimensions {length, width, height};
//    auto dim_max = max_element(max_dimensions.begin(), max_dimensions.end());
//
//    double max = *dim_max;
//
////    SMesh sm = b.mesh;
////    std::array<Point3, 8> obb_points;
////    CGAL::oriented_bounding_box(sm, obb_points);
////
////    SMesh obb_sm;
////    CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
////                          obb_points[4], obb_points[5], obb_points[6], obb_points[7], obb_sm);
////    std::ofstream("obb.off") << obb_sm;
////    PMP::triangulate_faces(obb_sm);
////    std::cout << "Volume: " << PMP::volume(obb_sm) << std::endl;
//    return max;
//}
//
//
//SMesh Liu_RoI(){
//    double px = (city.xmin+city.xmax)*0.5;
//    double py = (city.ymin+city.ymax)*0.5;
//    Point2 p(px, py);
//
//    vector<double> max_values;
//    for (unordered_map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {
//        double max = maximum_dimension(b->second);
//        max_values.push_back(max);
//    }
//    auto dim_max = max_element(max_values.begin(), max_values.end());
//    double lmax = *dim_max;
//
//    double r = 3*lmax;
//    meshCFD.RoICircle = Circle2(p, pow(r, 2));
//
//    SMesh RoI = RoI_cylinder(p, r);
//    meshCFD.RoI = RoI;
//
//    return RoI;
//}
//
//vector<size_t> convex_hull(vector<Point2> vertices) {
//    vector<size_t> indices(vertices.size()), out;
//    iota(indices.begin(), indices.end(),0);
//    CGAL::convex_hull_2(indices.begin(), indices.end(), back_inserter(out),Convex_hull_traits_2(CGAL::make_property_map(vertices)));
//    return out;
//}
//
//SMesh Liu_RoI(double x, double y){
//    double px = x;
//    double py = y;
//    Point2 p(px, py);
//
//    int id;
//
//    for (unordered_map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {
//
//        vector <size_t> a_out = convex_hull(b->second.gs2D);
//        vector <Point2> boundary;
//        for (size_t ia: a_out) {
//            boundary.push_back(b->second.gs2D[ia]);
//        }
//
//        switch (CGAL::bounded_side_2(boundary.begin(), boundary.end(), p, K())) {
//            case CGAL::ON_BOUNDED_SIDE:
//                id = b->first;
//                break;
//            case CGAL::ON_BOUNDARY:
//                id = b->first;
//                break;
//        }
//    }
//
//    double lmax = maximum_dimension(buildings[id]);
//    double r = 3*lmax;
//    meshCFD.RoICircle = Circle2(p, pow(r, 2));
//
//    SMesh RoI = RoI_cylinder(p, r);
//    meshCFD.RoI = RoI;
//
//    return RoI;
//}
//
//
//vector<int> Buildings_in_RoI(bool target, double x, double y) {
//
//    vector<int> detailed_buildings;
//
//    SMesh RoI;
//    if (target == false){
//        RoI = Liu_RoI();
//    } else{
//        RoI = Liu_RoI(x, y);
//    }
//
//    Circle2 circle = meshCFD.RoICircle;
//
//    int count = 0;
//
//    for (unordered_map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {
//        SMesh building = b->second.mesh;
//        SMesh out;
//
//        bool valid_intersection = PMP::do_intersect(meshCFD.RoI, building);
//        if (valid_intersection) {
//            detailed_buildings.push_back(b->second.id);
//        } else {
//            bool in = 0;
//            if (!b->second.gs2D.empty()){
//                for (auto &p: b->second.gs2D){
//                    if (!incircle(circle, p)) {
//                        in = 1;
//                    }
//                }
//
//                if (in == 0) {
//                    count += 1;
//                    detailed_buildings.push_back(b->second.id);
//                }
//            }
//        }
//    }
//
//    return detailed_buildings;
//}
//
//void generate_buildings_in_RoI(bool target, double x, double y){
//    vector<int> detailed_buildings = Buildings_in_RoI(target, x, y);
//
//    vector<SMesh> outputs;
//    for (auto db: detailed_buildings){
//        SMesh building = buildings[db].mesh;
//        outputs.push_back(building);
//    }
//
//    ofstream ofile("../output/all_detailed_buildings.obj");
//
//    map<int, vector<int>> f;
//    int fid = 0;
//    int vid = 0;
//
//    for (auto &each: outputs){
//        for (auto fit = each.faces_begin(); fit != each.faces_end(); ++fit){
//            fid += 1;
//            vector<int> pts;
//            CGAL::Vertex_around_face_iterator<SMesh> vbegin, vend;
//            for(boost::tie(vbegin, vend) = vertices_around_face(each.halfedge(*fit), each); vbegin != vend; ++vbegin){
//                Point3 v = each.point(*vbegin);
//                vid += 1;
//                pts.push_back(vid);
//                ofile << std::setprecision(5) << std::fixed << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
//            }
//            f[fid] = pts;
//        }
//    }
//
//    for (auto &ff: f){
//        ofile << "f " << ff.second[0] << " " << ff.second[1] << " " << ff.second[2] << endl;
//    }
//
//    ofile.close();
//}
//
//
//double lowest_volume(vector<int> detailed_buildings){
//
//    vector<double> volumes;
//    for (auto db: detailed_buildings) {
//        SMesh building = buildings[db].mesh;
//        double vol = CGAL::Polygon_mesh_processing::volume(building);
//        if (vol > 5){
//            volumes.push_back(vol);
//        }
//
//    }
//
//    auto vol_min = min_element(volumes.begin(), volumes.end());
//    double min = *vol_min;
//    cout << "min volume " << min << endl;
//
//    return min;
//}
//
//double mean_volume(vector<int> detailed_buildings){
//
//
//    double vdb = 0;
//    for (int db: detailed_buildings) {
//        SMesh building = buildings[db].mesh;
//        double vol = CGAL::Polygon_mesh_processing::volume(building);
//        vdb += vol;
//    }
//
//    double vb = 0;
//    for (unordered_map<int, Building>::iterator db=buildings.begin(); db!=buildings.end(); db++) {
//        SMesh building = db->second.mesh;
//        double vol = CGAL::Polygon_mesh_processing::volume(building);
//        vb += vol;
//    }
//
//    double vdb_mean = vdb/detailed_buildings.size();
//    double vb_mean = vb/buildings.size();
//
//    double v_mean;
//    if (vdb_mean <= vb_mean){
//        v_mean = vdb_mean;
//    } else {
//        v_mean = vb_mean;
//    }
//    //cout << "mean volume " << v_mean << endl;
//    return v_mean;
//}
//
//
//void check_cube_root_volume(){
//
////    double vol = v_mean;
////    double condition = 10*(meshCFD.wCell/16);
////    double cube_root_vol = cbrt(vol);
////
////    int check = 0;
////    if (cube_root_vol<condition && vol > 0){
////        check = 1;
////    }
////
////    if (check == 1  && meshCFD.wCell >= meshCFD.min){
////        if (meshCFD.wCell <= 0.1){
////            meshCFD.wCell = round5(meshCFD.wCell - 0.01);
////        } else  {
////            meshCFD.wCell = round5(meshCFD.wCell - 0.1);
////        }
////        check_cube_root_volume(v_mean);
////    }
//
//    cout << "BuildingVolume: " << endl;
//
//    double n_buildings = buildings.size();
//    double n = 0;
//    double condition = 10*(meshCFD.wCell/8);
//
//    double volume_total = 0;
//    double volume_valid = 0;
//
//    for (unordered_map<int, Building>::iterator db=buildings.begin(); db!=buildings.end(); db++) {
//        SMesh building = db->second.mesh;
//        double vol = CGAL::Polygon_mesh_processing::volume(building);
//        volume_total += vol;
//        double cube_root_vol = cbrt(vol);
//
//        if (cube_root_vol/(meshCFD.wCell/8)>=10){
//            db->second.invalid_volume = 0;
//            volume_valid += vol;
//            n += 1;
//        }
//
//        if (cube_root_vol/(meshCFD.wCell/8)<10){
//            db->second.invalid_volume = 1;
//        }
//    }
//
//    double n_mean = n/n_buildings;
//    double vol_mean = volume_valid/volume_total;
//
//    cout << "n_mean " << n_mean << endl;
//    cout << "vol_mean " << vol_mean << endl;
//
//}
//
//
//void BuildingVolume(bool target, double x, double y){
//
//    check_cube_root_volume();
//    computational_domain(meshCFD.wCell);
//    RefinementBoxes(meshCFD.wCell);
//
//}
//
////DEEL 2
////vector<size_t> convex_hull(vector<Point2> points) {
////
////    vector<size_t> indices(points.size()), out;
////    iota(indices.begin(), indices.end(),0);
////    CGAL::convex_hull_2(indices.begin(), indices.end(), back_inserter(out),Convex_hull_traits_2(CGAL::make_property_map(points)));
////
////    return out;
////}
//
//pair<pair<Point2, Point2>, double> distance_closest_building(vector<Point2> gsBuildingA, vector<Point2> gsBuildingB){
//    double d_min;
//
//    vector<size_t> a_out = convex_hull(gsBuildingA);
//    vector<size_t> b_out = convex_hull(gsBuildingB);
//
//    vector<double> distances;
//    vector<pair<int, int>> pts;
//    for (size_t ia: a_out){
//        for (size_t ib: b_out){
//            double d = sqrt(CGAL::squared_distance(gsBuildingA[ia], gsBuildingB[ib]));
//            pair<int, int> ab = {ia, ib};
//            if (0 <= d && d < 0.001){
//
//                //cout << "d " << d << " point A " << gsBuildingA[ia] << " point B " << gsBuildingB[ib] << endl;
//            }
//
//            distances.push_back(d);
//            pts.push_back(ab);
//        }
//    }
//
//    auto min = min_element(distances.begin(), distances.end());
//    int minIndex = distance(distances.begin(), min);
//    pair<pair<Point2, Point2>, double> r;
//
//    if (min != distances.end()) {
//        d_min = *min;
//        pair<Point2, Point2> ab = {gsBuildingA[pts[minIndex].first], gsBuildingB[pts[minIndex].second]};
//        r = make_pair(ab, d_min); //{ab, d_min};
//    } else {
//        d_min = -1;
//    }
//
//    return r;
//}
//
//
////threshold of 0 is not possible
//double shortest_distance(vector<int> detailed_buildings, double threshold){
//
//    double d = -1;
//
//    if (detailed_buildings.size() > 1) {
//        vector<double> distances;
//
//        for (auto &db: detailed_buildings) {
//            for (auto &dbb: detailed_buildings) {
//                if (db != dbb) {
//                    pair<pair<Point2, Point2>, double> dp = distance_closest_building(buildings[db].gs2D, buildings[dbb].gs2D);
//                    //cout << "check dp " << dp.second << endl;
//                    if (dp.second > threshold) {
//                        distances.push_back(dp.second);
//                    }
//                }
//            }
//        }
//
//        if (!distances.empty()){
//            auto min = min_element(distances.begin(), distances.end());
//            int minIndex = distance(distances.begin(), min);
//            d = distances[minIndex];
//        }
//    }
//
//    cout << "d min " << d << endl;
//
//    return d;
//}
//
//
//
//double mean_distance(vector<int> detailed_buildings, double threshold){
//
//    double d = -1;
//
//
//    if (detailed_buildings.size() > 1) {
//
//        int id_db = 0;
//        double dist_db = 0;
//        for (auto &db: detailed_buildings) {
//            for (auto &dbb: detailed_buildings) {
//                if (db != dbb) {
//                    pair<pair<Point2, Point2>, double> dp = distance_closest_building(buildings[db].gs2D, buildings[dbb].gs2D);
//
//                    //if (dp.second > threshold) {
//                        id_db += 1;
//                        dist_db += dp.second;
//                    //}
//                }
//            }
//        }
//
//        int id_b = 0;
//        double dist_b = 0;
//        for (unordered_map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {
//
//            for (unordered_map<int, Building>::iterator bb=buildings.begin(); bb!=buildings.end(); bb++) {
//                if (b != bb) {
//                    pair<pair<Point2, Point2>, double> dp = distance_closest_building(b->second.gs2D, bb->second.gs2D);
//
//
//                    dist_b += dp.second;
//                   // if (dp.second > threshold) {3
//                        id_b += 1;
//                        dist_b += dp.second;
//                    //}
//                }
//            }
//        }
//
//        if (id_db != 0 && id_b != 0){
//            double ddb_mean = dist_db/id_db;
//            double db_mean = dist_b/id_b;
//            if (ddb_mean <= db_mean){
//                d = ddb_mean;
//            } else {
//                d = db_mean;
//            }
//        }
//    }
//
//    //cout << "mean distance " << d << endl;
//
//    return d;
//}
//
//void shortest_separation(double threshold) {
//
//    double wCell3 = meshCFD.wCell / (double) 8;
//    double wCell2 = meshCFD.wCell / (double) 4;
//    double d_Cell3 = 2 * wCell3 * meshCFD.nCellsBetweenLevels;
//    double n_buildings = buildings.size();
//    double n = 0;
//
//    for (unordered_map<int, Building>::iterator b = buildings.begin(); b != buildings.end(); b++) {
//
//        double volb = CGAL::Polygon_mesh_processing::volume(b->second.mesh);
//        vector<double> distances;
//        vector < pair < pair < Point2, Point2 >, double >> separation_pairs;
//        vector < pair < pair < Point2, Point2 >, double >> invalid_separations;
//
//        for (unordered_map<int, Building>::iterator bb = buildings.begin(); bb != buildings.end(); bb++) {
//            if (b->first != bb->first) {
//                double volbb = CGAL::Polygon_mesh_processing::volume(bb->second.mesh);
//
//                if (b->second.gs2D.size() != 0 && bb->second.gs2D.size() != 0 && volb > 0 && volbb > 0) {
//                    pair<pair<Point2, Point2>, double> dp = distance_closest_building(b->second.gs2D, bb->second.gs2D);
//                    double factor = (dp.second - d_Cell3) / (double) wCell2;
//                    if (dp.second >= 0) {
//                        distances.push_back(dp.second);
//                        separation_pairs.push_back(dp);
//                        if (factor < 2) {
//                            invalid_separations.push_back(dp);
//                        }
//                    }
//                }
//            }
//
//            auto min = min_element(distances.begin(), distances.end());
//            int minIndex = distance(distances.begin(), min);
//            double d_min;
//            pair<pair<Point2, Point2>, double> sp;
//
//            if (min != distances.end()) {
//                d_min = *min;
//                sp = separation_pairs[minIndex];
//            }
//            b->second.d_shortest = sp;
//        }
//
//        b->second.invalid_separations = invalid_separations;
//    }
//}
//
//void check_building_separations(double d_mean, double threshold){
//
////    if (d_mean != -1){
////
////        double wCell3 = meshCFD.wCell/(double)8;
////        double wCell4 = meshCFD.wCell/(double)16;
////        double d_Cell4 = 2*wCell4*meshCFD.nCellsBetweenLevels;
////        double factor = (d_mean-d_Cell4)/(double)wCell3;
////
////        if ( factor < 2 && meshCFD.wCell >= meshCFD.min){
////            if (meshCFD.wCell <= 0.1){
////                meshCFD.wCell = round5(meshCFD.wCell - 0.01);
////            } else {
////                meshCFD.wCell = round5(meshCFD.wCell - 0.1);
////            }
////            check_building_separations(d_mean, threshold);
////        }
////    }
//
//    cout << "BuildingSeparation: " << endl;
//    if (d_mean != -1){
//        shortest_separation(2);
//        double wCell3 = meshCFD.wCell/(double)8;
//        double wCell2 = meshCFD.wCell/(double)4;
//        double d_Cell3 = 2*wCell3*meshCFD.nCellsBetweenLevels;
//        double n_buildings = buildings.size();
//        double n = 0;
//        for (unordered_map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {
//
//            double factor = (b->second.d_shortest.second-d_Cell3)/(double)wCell2;
//            //cout << "b->second.d_shortest.second " << b->second.d_shortest.second << endl;
//            //cout << "factor " << factor << endl;
//                if (factor >= 2){
//                    n += 1;
//                }
//        }
//
//        double n_mean = n/n_buildings;
//
//        cout << "n_mean " << n_mean << endl;
//
//
//
//    }
//}
//
//void BuildingSeparation(bool target, double x, double y, double threshold){
//
//    vector<int> detailed_buildings = Buildings_in_RoI(target, x, y);
//    double d_mean = mean_distance(detailed_buildings, threshold);
//    check_building_separations(d_mean, threshold);
//    computational_domain(meshCFD.wCell);
//    RefinementBoxes(meshCFD.wCell);
//
//}
//
//
////void PedestrianCells(bool target, double x, double y, double threshold, int it){
////    cout << "PedestrianCells " <<  endl;
////    cout << "wCell " << meshCFD.wCell << endl;
////
////    double wCell3 = meshCFD.wCell/(double)8;
////    double z0 = 2*wCell3;
////    double z1 = 4*wCell3;
//////    cout << "z0 " << z0 << endl;
//////    cout << "z1 " << z1 << endl;
////    if (z1 < 1.5 && it < 100){
//////        cout << "first case " << endl;
////        meshCFD.wCell = meshCFD.wCell+0.01;
////        computational_domain(meshCFD.wCell);
////        BuildingVolume(target, x, y);
////        BuildingSeparation(target, x, y, threshold);
////        RefinementBoxes(meshCFD.wCell);
////        it += 1;
////        PedestrianCells(target, x, y, threshold, it);
////    } else if (z0 <= 1.5 && z1 >= 2 && meshCFD.wCell >= 0.01){
//////        cout << "second case " << endl;
////        meshCFD.wCell = meshCFD.wCell;
////        computational_domain(meshCFD.wCell);
////        RefinementBoxes(meshCFD.wCell);
////    } else if (z0>2 && meshCFD.wCell >= 0.01){
//////        cout << "third case " << endl;
////        if (meshCFD.wCell <= 0.1){
////            meshCFD.wCell = round((meshCFD.wCell - 0.01)*10000)/10000;
////        } else {
////            meshCFD.wCell = round((meshCFD.wCell - 0.1)*10000)/10000;
////        }
////        computational_domain(meshCFD.wCell);
////        RefinementBoxes(meshCFD.wCell);
////        PedestrianCells(target, x, y, threshold, it);
////    }
////
////    computational_domain(meshCFD.wCell);
////    RefinementBoxes(meshCFD.wCell);
////
////}
//
//void min_pedestrian(){
//
//    cout << "\nMin pedestrian:" << endl;
//
//    meshCFD.wmin = (1.5*4)/(double)2.5;
//    double wi = (5*4)/(double)2.5;
//    meshCFD.wCell = wi;
//
//    cout << "min pedestrian " << meshCFD.wmin << endl;
//    cout << "max pedestrian " << meshCFD.wCell << "\n" << endl;
//
//}
//
//void min_roughness(){
//
//    cout << "Minimum roughness: " << endl;
//    double wmin_roughness = meshCFD.z0*2*8;
//
//    cout << "wmin " << meshCFD.wmin << endl;
//    cout << "wmin roughness " << wmin_roughness << endl;
//
//    if (wmin_roughness > meshCFD.wmin){
//        meshCFD.wmin = wmin_roughness;
//        cout << "wmin_roughness is higher than wmin" << endl;
//    }
//
//    if (wmin_roughness > meshCFD.wCell){
//        meshCFD.wCell = wmin_roughness;
//        double pedestrian_level = (meshCFD.wCell/(double)4)*2.5;
//        cout << "wmin_roughness is higher than wCell, pedestrian level at 1.5-5m is not respected and is " << pedestrian_level << "m" << endl;
//    }
//
//    cout << "wmin final " << meshCFD.wmin << endl;
//    cout << "wCell final " << meshCFD.wCell << "\n" << endl;
//}
//
//void computational_domain_Blocken() {
//
//    CGAL::Bbox_3 domain;
//
////    cout << "CD dimensions " << meshCFD.domain.xmin() << " " << meshCFD.domain.xmax() << " " << meshCFD.domain.ymin()
////         << " " << meshCFD.domain.ymax() << " " << meshCFD.domain.zmin() << " " << meshCFD.domain.zmax() << endl;
//
//    double L_urban = abs(city.ymax - city.ymin);
//    double L_domain = abs(meshCFD.domain.ymax() - meshCFD.domain.ymin());
//    double BR_horizontal = L_urban / L_domain;
//
//    double H_urban = abs(city.zmax - city.zmin);
//    double H_domain = abs(meshCFD.domain.zmax() - meshCFD.domain.zmin());
//    double BR_vertical = H_urban / H_domain;
//
//    double xmin_i = meshCFD.domain.xmin();
//    double xmax_i = meshCFD.domain.xmax();
//    double ymin_i = meshCFD.domain.ymin();
//    double ymax_i = meshCFD.domain.ymax();
//    double zmin_i = meshCFD.domain.zmin();
//    double zmax_i = meshCFD.domain.zmax();
//
//    double xmin_final, xmax_final;
//    double ymin_final, ymax_final;
//    double zmin_final, zmax_final;
//
//    double ymin_BRh = ymin_i;
//    double ymax_BRh = ymax_i;
//
//    while (BR_horizontal > 0.17) {
//        ymin_BRh = ymin_BRh - 0.1;
//        ymax_BRh = ymax_BRh + 0.1;
//        L_domain = abs(ymax_BRh - ymin_BRh);
//        BR_horizontal = L_urban / L_domain;
//    }
//
//    if (L_domain <= abs(xmax_i - xmin_i)) {
//
//        double zmax_BRv = zmax_i;
//
//        while (BR_vertical > 0.17) {
//            zmax_BRv = zmax_BRv + 0.1;
//            H_domain = abs(zmax_BRv - zmin_i);
//            BR_vertical = H_urban / H_domain;
//        }
//
//        pair<double, double> xComp = fit2cellsize(xmin_i, xmax_i, meshCFD.wCell, 0);
//        pair<double, double> yComp = fit2cellsize(ymin_BRh, ymax_BRh, meshCFD.wCell, 0);
//        pair<double, double> zComp = fit2cellsize(zmin_i, zmax_BRv, meshCFD.wCell, 1);
//
//        xmin_final = xComp.first;
//        xmax_final = xComp.second;
//        ymin_final = yComp.first;
//        ymax_final = yComp.second;
//        zmin_final = zComp.first;
//        zmax_final = zComp.second;
//
//    } else {
////        cout << "L_domain too large" << endl;
//
//        double ymin_BR = ymin_i;
//        double ymax_BR = ymax_i;
//        double zmax_BR = zmax_i;
//
//        double A_urban = abs(city.ymax - city.ymin) * abs(city.zmax - city.zmin);
//        double A_domain = abs(ymax_BR - ymin_BR) * abs(zmax_BR - zmin_i);
//        double BR = A_urban / A_domain;
//
//        while (BR > 0.03) {
//            ymin_BR = ymin_BR - 0.1;
//            ymax_BR = ymax_BR + 0.1;
//            zmax_BR = zmax_BR + 0.1;
//            A_domain = abs(ymax_BR - ymin_BR) * abs(zmax_BR - zmin_i);
//            BR = A_urban / A_domain;
//        }
//
//        pair<double, double> xComp = fit2cellsize(xmin_i, xmax_i, meshCFD.wCell, 0);
//        pair<double, double> yComp = fit2cellsize(ymin_BR, ymax_BR, meshCFD.wCell, 0);
//        pair<double, double> zComp = fit2cellsize(zmin_i, zmax_BR, meshCFD.wCell, 1);
//
//        xmin_final = xComp.first;
//        xmax_final = xComp.second;
//        ymin_final = yComp.first;
//        ymax_final = yComp.second;
//        zmin_final = zComp.first;
//        zmax_final = zComp.second;
//    }
//
//    domain = CGAL::Bbox_3(xmin_final, ymin_final, zmin_final, xmax_final, ymax_final, zmax_final);
//    meshCFD.domain = domain;
//
////    cout << "CD dimensions " << meshCFD.domain.xmin() << " " << meshCFD.domain.xmax() << " " << meshCFD.domain.ymin()
////         << " " << meshCFD.domain.ymax() << " " << meshCFD.domain.zmin() << " " << meshCFD.domain.zmax() << endl;
//}
//
//
//
//double resolution_check() {
//
//    CGAL::Bbox_3 box1 = meshCFD.refinementBox1;
//    CGAL::Bbox_3 box2 = meshCFD.refinementBox2;
//    //CGAL::Bbox_3 box3 = meshCFD.refinementBox3;
//    CGAL::Bbox_3 domain = meshCFD.domain;
//
//    double volB = 0;
//    double volG = 0;
//    double nG = 0;
//    for (unordered_map<int, Building>::iterator b = buildings.begin(); b != buildings.end(); b++) {
//        SMesh building = b->second.mesh;
//        double vol = CGAL::Polygon_mesh_processing::volume(building);
//        volB += vol;
//
//        for (auto &f: b->second.faces){
//            if (f.gs == 0){
//                Triangle3 tri(vertices[f.v[0]].p, vertices[f.v[1]].p, vertices[f.v[2]].p);
//                double area = sqrt(tri.squared_area());
//                //double n = area/(pow(meshCFD.wCell/16, 2));
//                double n = (area/(pow(meshCFD.wCell/8, 2)))*4;
//                double v = n*pow(meshCFD.wCell/8, 3);
//                volG += v;
//                nG += n;
//            }
//        }
//    }
//
//    double vol1 = (abs(box1.xmax()-box1.xmin())*abs(box1.ymax()-box1.ymin())*abs(box1.zmax()-box1.zmin()))-(volB+volG); //*(1+2*meshCFD.nCellsBetweenLevels*meshCFD.wCell);
//    double vol2 = (abs(box2.xmax()-box2.xmin())*abs(box2.ymax()-box2.ymin())*abs(box2.zmax()-box2.zmin()))-(vol1+volB+volG);
//    //double vol3 = (abs(box3.xmax()-box3.xmin())*abs(box3.ymax()-box3.ymin())*abs(box3.zmax()-box3.zmin()))-(vol1+volB+volG+vol2);
//    //double volDomain = (abs(domain.xmax()-domain.xmin())*abs(domain.ymax()-domain.ymin())*abs(domain.zmax()-domain.zmin()))-(vol1+volB+volG+vol2+vol3);
//    double volDomain = (abs(domain.xmax()-domain.xmin())*abs(domain.ymax()-domain.ymin())*abs(domain.zmax()-domain.zmin()))-(vol1+volB+volG+vol2);
//
////    cout << "volB " << volB << endl;
////    cout << "volG " << volG << endl;
////    cout << "vol1 " << vol1 << endl;
////    cout << "vol2 " << vol2 << endl;
////    //cout << "vol3 " << vol3 << endl;
////    cout << "volDomain " << volDomain << endl;
////    cout << "vol computed " << volB+volG+vol1+vol2+volDomain << endl;
////    cout << "vol domain " << (abs(domain.xmax()-domain.xmin())*abs(domain.ymax()-domain.ymin())*abs(domain.zmax()-domain.zmin())) << endl;
//
////    double n1 = vol1/pow((meshCFD.wCell/8),3);
////    double n2 = vol2/pow((meshCFD.wCell/4),3);
////    double n3 = vol3/pow((meshCFD.wCell/2),3);
////    double nDomain = volDomain/pow(meshCFD.wCell,3);
//
//    double n1 = vol1/pow((meshCFD.wCell/4),3);
//    double n2 = vol2/pow((meshCFD.wCell/2),3);
//    double nDomain = volDomain/pow(meshCFD.wCell,3);
//
////    cout << "nG " << nG << endl;
////    cout << "n1 " << n1 << endl;
////    cout << "n2 " << n2 << endl;
////    cout << "n3 " << n3 << endl;
////    cout << "nDomain " << nDomain << endl;
//
//    double N = nG + n1 + n2 + nDomain;
////    cout << "N blockMesh " << (abs(domain.xmax()-domain.xmin())*abs(domain.ymax()-domain.ymin())*abs(domain.zmax()-domain.zmin()))/pow(meshCFD.wCell,3) << endl;
////    cout << "N snappyHexMesh " << N << endl;
//
//    return N;
//
//}
//
//void resolution(double threshold, double it){
//    double N = resolution_check();
//
//    if (N > threshold){
////        cout << "Resolution too high: " << N << endl;
////        cout << "Difference: " << N-threshold << endl;
//    }
//
//    while (N>threshold){
//        meshCFD.wCell = meshCFD.wCell + 0.1;
//        computational_domain(meshCFD.wCell);
//        computational_domain_Blocken();
//        RefinementBoxes(meshCFD.wCell);
//        N = resolution_check();
////        cout << "N " << N << endl;
////        cout << "wcell " << meshCFD.wCell << endl;
//    }
//
//    cout << "Resolution: " << endl;
//    cout << "N blockMesh " << (abs(meshCFD.domain.xmax()-meshCFD.domain.xmin())*abs(meshCFD.domain.ymax()-meshCFD.domain.ymin())*abs(meshCFD.domain.zmax()-meshCFD.domain.zmin()))/pow(meshCFD.wCell,3) << endl;
//    cout << "N snappyHexMesh " << N << "\n" << endl;
//    cout << "wcell " << meshCFD.wCell << endl;
//}
//
//
//
//void Meshing(bool target, double x, double y, double threshold){
//
//    meshCFD.nCellsBetweenLevels = 4;
//    min_pedestrian();
//    meshCFD.z0 = 2;
//    min_roughness();
//
//    computational_domain(meshCFD.wCell);
//    computational_domain_Blocken();
//    RefinementBoxes(meshCFD.wCell);
//
//    double th_resolution = 2400000;
//    double iteration = 0;
//    resolution(th_resolution, iteration);
//
//    BuildingVolume(target, x, y);
//    cout << "\n" << endl;
//    BuildingSeparation(target, x, y, threshold);
//    cout << "\n" << endl;
//
//    cout << "wCell " << meshCFD.wCell << endl;
//            cout << "CD dimensions " << meshCFD.domain.xmin() << " " << meshCFD.domain.xmax() << " " << meshCFD.domain.ymin()
//         << " " << meshCFD.domain.ymax() << " " << meshCFD.domain.zmin() << " " << meshCFD.domain.zmax() << endl;
//    cout << "N blockMesh " << (abs(meshCFD.domain.xmax()-meshCFD.domain.xmin())*abs(meshCFD.domain.ymax()-meshCFD.domain.ymin())*abs(meshCFD.domain.zmax()-meshCFD.domain.zmin()))/pow(meshCFD.wCell,3) << endl;
//
//}
//
//
//
//void min_pedestrian3(){
//
//    cout << "\nMin pedestrian:" << endl;
//
//    meshCFD.wmin = (1.5*8)/(double)2.5;
//    double wi = (5*8)/(double)2.5;
//    meshCFD.wCell = wi;
//
//    cout << "min pedestrian " << meshCFD.wmin << endl;
//    cout << "max pedestrian " << meshCFD.wCell << "\n" << endl;
//
//}
//
//void min_roughness3(){
//
//    cout << "Minimum roughness: " << endl;
//    double wmin_roughness = meshCFD.z0*2*16;
//
//    cout << "wmin " << meshCFD.wmin << endl;
//    cout << "wmin roughness " << wmin_roughness << endl;
//
//    if (wmin_roughness > meshCFD.wmin){
//        meshCFD.wmin = wmin_roughness;
//        cout << "wmin_roughness is higher than wmin" << endl;
//    }
//
//    if (wmin_roughness > meshCFD.wCell){
//        meshCFD.wCell = wmin_roughness;
//        double pedestrian_level = (meshCFD.wCell/(double)8)*2.5;
//        cout << "wmin_roughness is higher than wCell, pedestrian level at 1.5-5m is not respected and is " << pedestrian_level << "m" << endl;
//    }
//
//    cout << "wmin final " << meshCFD.wmin << endl;
//    cout << "wCell final " << meshCFD.wCell << "\n" << endl;
//}
//
//
//void RefinementBoxes3(double wCell){
////    double finlet, double foutlet, double flateral, double ftop
//
//    CGAL::Bbox_3 box1 = Refinement_box(wCell, 0, 0, 0, 1.5);
//    CGAL::Bbox_3 box2 = Refinement_box(wCell, 2, 6, 2, 3);
//    CGAL::Bbox_3 box3 = Refinement_box(wCell, 2, 6, 2, 3);
//
//    meshCFD.refinementBox1 = box1;
//    meshCFD.refinementBox2 = box2;
//    meshCFD.refinementBox3 = box3;
//
////    cout << "REFINEMENT BOX 1" << endl;
////    cout << box1.xmin() << endl;
////    cout << box1.xmax() << endl;
////    cout << box1.ymin() << endl;
////    cout << box1.ymax() << endl;
////    cout << box1.zmin() << endl;
////    cout << box1.zmax() << endl;
////
////    cout << "REFINEMENT BOX 2" << endl;
////    cout << box2.xmin() << endl;
////    cout << box2.xmax() << endl;
////    cout << box2.ymin() << endl;
////    cout << box2.ymax() << endl;
////    cout << box2.zmin() << endl;
////    cout << box2.zmax() << endl;
////
////    cout << "REFINEMENT BOX 3" << endl;
////    cout << box3.xmin() << endl;
////    cout << box3.xmax() << endl;
////    cout << box3.ymin() << endl;
////    cout << box3.ymax() << endl;
////    cout << box3.zmin() << endl;
////    cout << box3.zmax() << endl;
//}
//
//
//double resolution_check3() {
//
//    CGAL::Bbox_3 box1 = meshCFD.refinementBox1;
//    CGAL::Bbox_3 box2 = meshCFD.refinementBox2;
//    CGAL::Bbox_3 box3 = meshCFD.refinementBox3;
//    CGAL::Bbox_3 domain = meshCFD.domain;
//
//    double volB = 0;
//    double volG = 0;
//    double nG = 0;
//    for (unordered_map<int, Building>::iterator b = buildings.begin(); b != buildings.end(); b++) {
//        SMesh building = b->second.mesh;
//        double vol = CGAL::Polygon_mesh_processing::volume(building);
//        volB += vol;
//
//        for (auto &f: b->second.faces){
//            if (f.gs == 0){
//                Triangle3 tri(vertices[f.v[0]].p, vertices[f.v[1]].p, vertices[f.v[2]].p);
//                double area = sqrt(tri.squared_area());
//                double n = (area/(pow(meshCFD.wCell/16, 2)))*4;
//                double v = n*pow(meshCFD.wCell/16, 3);
//                volG += v;
//                nG += n;
//            }
//        }
//    }
//
//    double vol1 = (abs(box1.xmax()-box1.xmin())*abs(box1.ymax()-box1.ymin())*abs(box1.zmax()-box1.zmin()))-(volB+volG); //*(1+2*meshCFD.nCellsBetweenLevels*meshCFD.wCell);
//    double vol2 = (abs(box2.xmax()-box2.xmin())*abs(box2.ymax()-box2.ymin())*abs(box2.zmax()-box2.zmin()))-(vol1+volB+volG);
//    double vol3 = (abs(box3.xmax()-box3.xmin())*abs(box3.ymax()-box3.ymin())*abs(box3.zmax()-box3.zmin()))-(vol1+volB+volG+vol2);
//    double volDomain = (abs(domain.xmax()-domain.xmin())*abs(domain.ymax()-domain.ymin())*abs(domain.zmax()-domain.zmin()))-(vol1+volB+volG+vol2+vol3);
//
////    cout << "volB " << volB << endl;
////    cout << "volG " << volG << endl;
////    cout << "vol1 " << vol1 << endl;
////    cout << "vol2 " << vol2 << endl;
////    //cout << "vol3 " << vol3 << endl;
////    cout << "volDomain " << volDomain << endl;
////    cout << "vol computed " << volB+volG+vol1+vol2+volDomain << endl;
////    cout << "vol domain " << (abs(domain.xmax()-domain.xmin())*abs(domain.ymax()-domain.ymin())*abs(domain.zmax()-domain.zmin())) << endl;
//
//    double n1 = vol1/pow((meshCFD.wCell/8),3);
//    double n2 = vol2/pow((meshCFD.wCell/4),3);
//    double n3 = vol3/pow((meshCFD.wCell/2),3);
//    double nDomain = volDomain/pow(meshCFD.wCell,3);
//
////    cout << "nG " << nG << endl;
////    cout << "n1 " << n1 << endl;
////    cout << "n2 " << n2 << endl;
////    cout << "n3 " << n3 << endl;
////    cout << "nDomain " << nDomain << endl;
//
//    double N = nG + n1 + n2 + n3 + nDomain;
//    N = (abs(meshCFD.domain.xmax()-meshCFD.domain.xmin())*abs(meshCFD.domain.ymax()-meshCFD.domain.ymin())*abs(meshCFD.domain.zmax()-meshCFD.domain.zmin()))/pow(meshCFD.wCell,3);
//    return N;
//}
//
//void resolution3(double threshold, double it){
//    double N = resolution_check3();
//
//    while (N>threshold){
//        meshCFD.wCell = meshCFD.wCell + 0.1;
//        computational_domain(meshCFD.wCell);
//        computational_domain_Blocken();
//        RefinementBoxes3(meshCFD.wCell);
//        N = resolution_check3();
//    }
//
//
//    cout << "Resolution: " << endl;
//    cout << "N blockMesh " << (abs(meshCFD.domain.xmax()-meshCFD.domain.xmin())*abs(meshCFD.domain.ymax()-meshCFD.domain.ymin())*abs(meshCFD.domain.zmax()-meshCFD.domain.zmin()))/pow(meshCFD.wCell,3) << endl;
//    cout << "N snappyHexMesh " << N << "\n" << endl;
//    cout << "wcell " << meshCFD.wCell << endl;
//}
//
//void check_cube_root_volume3(){
//
//    cout << "BuildingVolume: " << endl;
//
//    double n_buildings = buildings.size();
//    double n = 0;
//    double condition = 10*(meshCFD.wCell/16);
//
//    double volume_total = 0;
//    double volume_valid = 0;
//
//    for (unordered_map<int, Building>::iterator db=buildings.begin(); db!=buildings.end(); db++) {
//        SMesh building = db->second.mesh;
//        double vol = CGAL::Polygon_mesh_processing::volume(building);
//        volume_total += vol;
//        double cube_root_vol = cbrt(vol);
//
//        if (cube_root_vol/(meshCFD.wCell/16)>=10){
//            db->second.invalid_volume = 0;
//            volume_valid += vol;
//            n += 1;
//        }
//
//        if (cube_root_vol/(meshCFD.wCell/16)<10){
//            db->second.invalid_volume = 1;
//        }
//    }
//
//    double n_mean = n/n_buildings;
//    double vol_mean = volume_valid/volume_total;
//
//    cout << "n_mean " << n_mean << endl;
//    cout << "vol_mean " << vol_mean << endl;
//
//}
//
//
//void BuildingVolume3(bool target, double x, double y){
//    check_cube_root_volume3();
//    computational_domain(meshCFD.wCell);
//    RefinementBoxes3(meshCFD.wCell);
//}
//
//void check_building_separations3(double d_mean, double threshold){
//
//    cout << "BuildingSeparation: " << endl;
//    if (d_mean != -1){
//        shortest_separation(2);
//        double wCell4 = meshCFD.wCell/(double)16;
//        double wCell3 = meshCFD.wCell/(double)8;
//        double d_Cell4 = 2*wCell4*meshCFD.nCellsBetweenLevels;
//        double n_buildings = buildings.size();
//        double n = 0;
//        for (unordered_map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {
//
//            double factor = (b->second.d_shortest.second-d_Cell4)/(double)wCell3;
//            if (factor >= 2){
//                n += 1;
//            }
//        }
//
//        double n_mean = n/n_buildings;
//
//        cout << "n_mean " << n_mean << endl;
//
//    }
//}
//
//void BuildingSeparation3(bool target, double x, double y, double threshold){
//
//    vector<int> detailed_buildings = Buildings_in_RoI(target, x, y);
//    double d_mean = mean_distance(detailed_buildings, threshold);
//    check_building_separations3(d_mean, threshold);
//    computational_domain(meshCFD.wCell);
//    RefinementBoxes3(meshCFD.wCell);
//
//}
//
//void Meshing3(bool target, double x, double y, double threshold){
//
//    meshCFD.nCellsBetweenLevels = 4;
//    min_pedestrian3();
//    meshCFD.z0 = 0.2;
//    min_roughness3();
//
//    computational_domain(meshCFD.wCell);
//    computational_domain_Blocken();
//    RefinementBoxes3(meshCFD.wCell);
//
//    double th_resolution = 2400000;
//    double iteration = 0;
//    resolution3(th_resolution, iteration);
//
//    BuildingVolume3(target, x, y);
//    cout << "\n" << endl;
//    BuildingSeparation(target, x, y, threshold);
//    cout << "\n" << endl;
//
//    cout << "wCell " << meshCFD.wCell << endl;
//    cout << "CD dimensions " << meshCFD.domain.xmin() << " " << meshCFD.domain.xmax() << " " << meshCFD.domain.ymin()
//         << " " << meshCFD.domain.ymax() << " " << meshCFD.domain.zmin() << " " << meshCFD.domain.zmax() << endl;
//    cout << "N blockMesh " << (abs(meshCFD.domain.xmax()-meshCFD.domain.xmin())*abs(meshCFD.domain.ymax()-meshCFD.domain.ymin())*abs(meshCFD.domain.zmax()-meshCFD.domain.zmin()))/pow(meshCFD.wCell,3) << endl;
//
//    vector<int> detailed_buildings = Buildings_in_RoI(target, x, y);
//    lowest_volume(detailed_buildings);
//    shortest_distance(detailed_buildings, 2);
//
//
//}
//
//
//void write_blockMeshDict() {
//
//    meshCFD.nx = (meshCFD.domain.xmax()-meshCFD.domain.xmin())/(double)meshCFD.wCell;
//    meshCFD.ny = (meshCFD.domain.ymax()-meshCFD.domain.ymin())/(double)meshCFD.wCell;
//    meshCFD.nz = (meshCFD.domain.zmax()-meshCFD.domain.zmin())/(double)meshCFD.wCell;
////    cout << "nx " << meshCFD.nx << endl;
////    cout << "nx " << meshCFD.nx << endl;
////    cout << "ny " << meshCFD.ny << endl;
////    cout << "nz " << meshCFD.nz << endl;
//
//    ifstream input_stream;
//    input_stream.open("../data/blockMeshDict");
//
//    ofstream ofile("../output/blockMeshDict");
//
//    int i = 0;
//    if (input_stream.is_open()) {
//        string line, data;
//
//        while (getline(input_stream, line)) {
//            istringstream tmp_stream(line);
//            tmp_stream >> data;
//
//            if (data == "xMin") {
//                ofile << "    xMin  " << setprecision(5) << fixed << meshCFD.domain.xmin() << ";" << endl;
//            }
//            else if (data == "xMax") {
//                ofile << "    xMax  " << setprecision(5) << fixed << meshCFD.domain.xmax() << ";" << endl;
//            }
//
//            else if(data == "yMin") {
//                ofile << "    yMin  " << setprecision(5) << fixed << meshCFD.domain.ymin() << ";" << endl;
//            }
//            else if(data == "yMax") {
//                ofile << "    yMax  " << setprecision(5) << fixed << meshCFD.domain.ymax() << ";" << endl;
//            }
//
//            else if(data == "zMin") {
//                ofile << "    zMin  " << setprecision(5) << fixed << meshCFD.domain.zmin() << ";" << endl;
//            }
//            else if(data == "zMax") {
//                ofile << "    zMax  " << setprecision(5) << fixed << meshCFD.domain.zmax() << ";" << endl;
//            }
//
//            else if(data == "xCells") {
//                ofile << "    xCells  " << setprecision(0) << fixed << (meshCFD.nx) << ";" << endl;
//            }
//            else if(data == "yCells") {
//                ofile << "    yCells  " << setprecision(0) << fixed << (meshCFD.ny) << ";" << endl;
//            }
//            else if(data == "zCells") {
//                ofile << "    zCells  " << setprecision(0) << fixed << (meshCFD.nz) << ";" << endl;
//            }
//
//            else {
//                ofile << line << endl;
//            }
//        }
//    }
//
//    ofile.close();
//}
//
//
//Point3 inside_point(){
//    Point3 p (meshCFD.refinementBox2.xmin(), meshCFD.refinementBox2.ymin(), meshCFD.refinementBox2.zmin());
//    return p;
//}
//
//void write_snappyHexMeshDict(string filename){
//    ifstream input_stream;
//    input_stream.open("../data/snappyHexMeshDict2");
//    ofstream ofile("../output/snappyHexMeshDict");
//
//    int idx = filename.find_last_of(".");
//    string name = filename.substr(0, idx);
//
//    int i = 0;
//    if (input_stream.is_open()) {
//        string line, data;
//
//        while (getline(input_stream, line)) {
//            istringstream tmp_stream(line);
//            tmp_stream >> data;
//
//            if (data == "input1") {
//                ofile << "        file \"" << filename << "\"" << endl;
//            }
//            else if (data == "input2") {
//                ofile << "        min  (" << setprecision(5) << fixed <<meshCFD.refinementBox1.xmin() << " " << setprecision(5) << fixed <<meshCFD.refinementBox1.ymin() << " " << setprecision(5) << fixed << meshCFD.refinementBox1.zmin() << ");" << endl;
//            }
//
//            else if(data == "input3") {
//                ofile << "        max  (" << setprecision(5) << fixed << meshCFD.refinementBox1.xmax() << " " << setprecision(5) << fixed << meshCFD.refinementBox1.ymax() << " " << setprecision(5) << fixed << meshCFD.refinementBox1.zmax() << ");" << endl;
//            }
//
//            else if (data == "input4") {
//                ofile << "        min  (" << meshCFD.refinementBox2.xmin() << " " << meshCFD.refinementBox2.ymin() << " " << meshCFD.refinementBox2.zmin() << ");" << endl;
//            }
//
//            else if(data == "input5") {
//                ofile << "        max  (" << meshCFD.refinementBox2.xmax() << " " << meshCFD.refinementBox2.ymax() << " " << meshCFD.refinementBox2.zmax() << ");" << endl;
//            }
//
////            else if (data == "input6") {
////                ofile << "        min  (" << meshCFD.refinementBox3.xmin() << " " << meshCFD.refinementBox3.ymin() << " " << meshCFD.refinementBox3.zmin() << ");" << endl;
////            }
//
////            else if(data == "input7") {
////                ofile << "        max  (" << meshCFD.refinementBox3.xmax() << " " << meshCFD.refinementBox3.ymax() << " " << meshCFD.refinementBox3.zmax() << ");" << endl;
////            }
//
//
//            else if(data == "input8") {
//                ofile << "        { file " << "\"" << name << ".eMesh" << "\"; level 1; }" << endl;
//            }
//
//            else if(data == "input9") {
//                Point3 p = inside_point();
//                ofile << "    locationInMesh (" << p.x() << " " << p.y() << " " << p.z()<<");"<< endl;
//            }
//
//            else {
//                ofile << line << endl;
//            }
//        }
//    }
//
//    ofile.close();
//}
//
//
//void write_snappyHexMeshDict3(string filename){
//    ifstream input_stream;
//    input_stream.open("../data/snappyHexMeshDict");
//    ofstream ofile("../output/snappyHexMeshDict");
//
//    int idx = filename.find_last_of(".");
//    string name = filename.substr(0, idx);
//
//    int i = 0;
//    if (input_stream.is_open()) {
//        string line, data;
//
//        while (getline(input_stream, line)) {
//            istringstream tmp_stream(line);
//            tmp_stream >> data;
//
//            if (data == "input1") {
//                ofile << "        file \"" << filename << "\"" << endl;
//            }
//            else if (data == "input2") {
//                ofile << "        min  (" << setprecision(5) << fixed <<meshCFD.refinementBox1.xmin() << " " << setprecision(5) << fixed <<meshCFD.refinementBox1.ymin() << " " << setprecision(5) << fixed << meshCFD.refinementBox1.zmin() << ");" << endl;
//            }
//
//            else if(data == "input3") {
//                ofile << "        max  (" << setprecision(5) << fixed << meshCFD.refinementBox1.xmax() << " " << setprecision(5) << fixed << meshCFD.refinementBox1.ymax() << " " << setprecision(5) << fixed << meshCFD.refinementBox1.zmax() << ");" << endl;
//            }
//
//            else if (data == "input4") {
//                ofile << "        min  (" << meshCFD.refinementBox2.xmin() << " " << meshCFD.refinementBox2.ymin() << " " << meshCFD.refinementBox2.zmin() << ");" << endl;
//            }
//
//            else if(data == "input5") {
//                ofile << "        max  (" << meshCFD.refinementBox2.xmax() << " " << meshCFD.refinementBox2.ymax() << " " << meshCFD.refinementBox2.zmax() << ");" << endl;
//            }
//
//            else if (data == "input6") {
//                ofile << "        min  (" << meshCFD.refinementBox3.xmin() << " " << meshCFD.refinementBox3.ymin() << " " << meshCFD.refinementBox3.zmin() << ");" << endl;
//            }
//
//            else if(data == "input7") {
//                ofile << "        max  (" << meshCFD.refinementBox3.xmax() << " " << meshCFD.refinementBox3.ymax() << " " << meshCFD.refinementBox3.zmax() << ");" << endl;
//            }
//
//
//            else if(data == "input8") {
//                ofile << "        { file " << "\"" << name << ".eMesh" << "\"; level 1; }" << endl;
//            }
//
//            else if(data == "input9") {
//                Point3 p = inside_point();
//                ofile << "    locationInMesh (" << p.x() << " " << p.y() << " " << p.z()<<");"<< endl;
//            }
//
//            else {
//                ofile << line << endl;
//            }
//        }
//    }
//
//    ofile.close();
//}
//
//
