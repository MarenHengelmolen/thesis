#include "definitions.h"
#include "input.h"
#include "val3dity.h"
#include "nlohmann-json/json.hpp"
#include "geometry.h"

int run_val3dity(string input_file, string val3dity_report, double snap_tolerance, double planarity_tolerance, double overlap_tolerance){

    //runs validations from val3dity and creates a json report with the results for separate building and terrain validations
    //from https://github.com/tudelft3d/val3dity
    int valid;
    std::stringstream buffer;
    buffer << std::ifstream(input_file).rdbuf();
    try {
        std::string s = buffer.str();
        bool re = val3dity::is_valid(s, "OBJ", snap_tolerance, planarity_tolerance, 20, overlap_tolerance);

        if (re == true) {
            std::cout << "VALID!" << std::endl;
            valid = 1;
        }
        else {
            std::cout << "INVALID :(" << std::endl;
            valid = 0;
            json report = val3dity::validate(s, "OBJ");
            std::ofstream file(val3dity_report);
            file << report;
            file.close();
        }
    }
    catch (std::exception& ex) {
        std::cerr << ex.what() << std::endl;
    }

    output.valid = valid;

    return valid;
}


vector<pair<pair<Point3, Point3>, double>> identify_short_edges(double threshold){

    vector<pair<pair<Point3, Point3>, double>> short_edges;
    vector<Point3> checkA;
    vector<Point3> checkB;

    for (map<int, Face>::iterator f=fBuildings.begin(); f!= fBuildings.end(); f++) {
            int v1 = f->second.v[0];
            int v2 = f->second.v[1];
            int v3 = f->second.v[2];

            //for each face f, find edges
            pair<Point3, Point3> e1 = make_pair(vertices[v1].p, vertices[v2].p);
            pair<Point3, Point3> e2 = make_pair(vertices[v2].p, vertices[v3].p);
            pair<Point3, Point3> e3 = make_pair(vertices[v3].p, vertices[v1].p);
            vector<pair<Point3, Point3>> edges {e1, e2, e3};

            for (auto &e: edges){
                //for each edge in f, compute distance of each edge e
                double distance = sqrt(CGAL::squared_distance(e.first, e.second));

                //if this distance is lower than a threshold
                //and not equal to 0 (edges with duplicated edges are filtered out)
                //then the edge is stored as short edge
                if (distance <= threshold && distance != 0){

                    //check if edge is already identified as short edge
                    //here an edge is a pair of two vertices e1 and e2
                    //pair e1 and e2 represents the same edge as pair e2 and e1
                    auto itA = find(checkA.begin(), checkA.end(), e.first);
                    auto itB = find(checkA.begin(), checkA.end(), e.second);

                    if (itA != checkA.end())
                    {
                        int index = itA - checkA.begin();
                        if (checkB[index] != e.second){
                            checkA.push_back(e.first);
                            checkB.push_back(e.second);
                            short_edges.push_back(make_pair(make_pair(e.first, e.second), distance));
                        }
                    }

                    else if (itB != checkA.end())
                    {
                        int index = itB - checkA.begin();
                        if (checkB[index] != e.first){
                            checkA.push_back(e.first);
                            checkB.push_back(e.second);
                            short_edges.push_back(make_pair(make_pair(e.first, e.second), distance));
                        }
                    }

                    else {
                        checkA.push_back(e.first);
                        checkB.push_back(e.second);
                        short_edges.push_back(make_pair(make_pair(e.first, e.second), distance));
                    }
                }
            }
    }
    return short_edges;
}

vector<pair<pair<int, int>, double>> identify_sharp_angles(double threshold, string input, string OBJfilename, string TXTfilename){

    vector<pair<pair<int, int>, double>> sharp_angles;

    //generate a mesh from file with buildings and the PMP package
    SMesh mesh;
    if(!PMP::IO::read_polygon_mesh(input, mesh) || !CGAL::is_triangle_mesh(mesh)){
        cerr << "Invalid input." << std::endl;
    }

    //create a new map with faces, using ids defined by the PMP package
    map<int, FaceSA> store_faces;
    for (auto f : mesh.faces()){
        vector<Point3> vertices;
        for (auto v: CGAL::vertices_around_face(mesh.halfedge(f), mesh)) {
            Point3 p = mesh.point(v);
            vertices.push_back(p);
        }
        store_faces[f.idx()].id = f.idx();
        store_faces[f.idx()].v = vertices;
    }

    for (auto f : mesh.faces()){

        //define vertices of face f
        vector<Point3> f_pts;
        vector<int> f_idx;
        map<int, Point3> map_pts;

        for (auto v: CGAL::vertices_around_face(mesh.halfedge(f), mesh)){
            Point3 p = mesh.point(v);
            int idx = v.idx();
            f_pts.push_back(p);
            f_idx.push_back(idx);
            map_pts[idx] = p;
        }

        Point3 f_v1 = f_pts[0];
        Point3 f_v2 = f_pts[1];
        Point3 f_v3 = f_pts[2];

        //find neighbouring faces of f
        vector<int> neighbours;

        for(auto n : CGAL::faces_around_face(mesh.halfedge(f), mesh)) {

            //if id of neighbouring face equals -1, it could indicate that:
            //- there are no neighbouring faces anymore, or
            //- there are more than 3 neighbouring faces
            if (n.idx() == -1){

                //find a possible other neighbour not identified by the PMP package
                //iterate over faces ff in mesh
                for (auto ff: mesh.faces()) {
                    if (f.idx() != ff.idx()){

                        //if the faces f and ff are not the same,
                        //check if they have 2 vertices in common
                        int v_common = 0;
                        for (auto v: CGAL::vertices_around_face(mesh.halfedge(f), mesh)) {
                            for (auto vf: CGAL::vertices_around_face(mesh.halfedge(ff), mesh)) {
                                Point3 p = mesh.point(v);
                                Point3 pf = mesh.point(vf);
                                if (p == pf) {
                                    v_common += 1;
                                }
                            }
                        }


                        if (v_common == 2 && !(find(neighbours.begin(), neighbours.end(), ff.idx()) != neighbours.end())){

                            //if they have 2 vertices in common, and ff is not yet identified by the PMP package
                            //check if they have also 2 uncommon vertices
                            vector<Point3> n_pts;
                            vector<int> pts;
                            vector<Point3> pts_2;
                            pts_2 = f_pts;
                            pts = f_idx;
                            vector<int> unique_values;
                            set<int> duplicates;
                            for (auto v: CGAL::vertices_around_face(mesh.halfedge(ff), mesh)) {
                                Point3 p = mesh.point(v);
                                int idx = v.idx();
                                n_pts.push_back(p);
                                pts_2.push_back(p);
                                pts.push_back(idx);
                                map_pts[idx] = p;
                            }

                            Point3 n_v1 = n_pts[0];
                            Point3 n_v2 = n_pts[1];
                            Point3 n_v3 = n_pts[2];

                            vector<Point3> common_pts;
                            vector<Point3> uncommon_pts;
                            for (auto &each: pts_2) {
                                int count = std::count(pts_2.begin(), pts_2.end(), each);
                                if (count == 2){
                                    auto it = find(common_pts.begin(), common_pts.end(), each);
                                    if (it == common_pts.end()) {
                                        common_pts.push_back(each);
                                    }
                                } else{
                                    uncommon_pts.push_back(each);
                                }
                            }

                            std::sort(common_pts.begin(), common_pts.end());
                            auto unique_end = unique(common_pts.begin(), common_pts.end());
                            pts_2.erase(unique_end, pts_2.end());

                            if (common_pts.size() == 2 && uncommon_pts.size() == 2) {
                                //if they have 2 common and 2 uncommon vertices, compute angles between these two faces
                                if (sqrt(CGAL::squared_area(f_v1, f_v2, f_v3)) > 0 && sqrt(CGAL::squared_area(n_v1, n_v2, n_v3)) > 0 ) {
                                    Point3 P = common_pts[0];
                                    Point3 a0 = common_pts[1];
                                    Point3 a1 = uncommon_pts[0];
                                    Point3 a2 = uncommon_pts[1];

                                    Vector3 b0(a0.x() - P.x(), a0.y() - P.y(), a0.z() - P.z());
                                    Vector3 b1(a1.x() - P.x(), a1.y() - P.y(), a1.z() - P.z());
                                    Vector3 b2(a2.x() - P.x(), a2.y() - P.y(), a2.z() - P.z());

                                    Vector3 A = CGAL::cross_product(b0, b1);
                                    Vector3 B = CGAL::cross_product(b0, b2);

                                    //using their dihedral angle
                                    double phi1 = (acos(CGAL::scalar_product(CGAL::cross_product(b0, b1), CGAL::cross_product(b0, b2)) /
                                                        (pow(A.squared_length(), 0.5) * pow(B.squared_length(), 0.5)))) * (180 / M_PI);
                                    double phi2 = 360 - abs(phi1);


                                    //if one of these angles is lower than threshold and not equal to 0 (to avoid considering identical faces as sharp angles)
                                    //a sharp angle is considered between these two faces
                                    //check if this angle is already identified as sharp angle
                                    //here a sharp angle is a pair of two faces f and ff with their angle
                                    //pair f and ff represents the same sharp angle as pair ff and f
                                    if (abs(phi1) <= threshold && phi1 != 0) {
                                        pair<pair<int, int>, double> pair;
                                        pair = make_pair(make_pair(ff.idx(), f.idx()), phi1);
                                        auto it = find(sharp_angles.begin(), sharp_angles.end(), pair);
                                        if (it == sharp_angles.end()) {
                                            sharp_angles.push_back(make_pair(make_pair(f.idx(), ff.idx()), phi1));
                                        }
                                    }

                                    if (phi2 <= threshold && phi2 != 0) {
                                        pair<pair<int, int>, double> pair;
                                        pair = make_pair(make_pair(ff.idx(), f.idx()), phi2);
                                        auto it = find(sharp_angles.begin(), sharp_angles.end(), pair);
                                        if (!(it != sharp_angles.end())) {
                                            sharp_angles.push_back(make_pair(make_pair(f.idx(), ff.idx()), phi2));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            else{
                //if id of neighbouring face is not equal to -1, perform the same steps
                //except iterating over the faces of the mesh, as the neighbouring face is correctly identified
                neighbours.push_back(n.idx());
                vector <Point3> n_pts;
                vector<int> pts;
                pts = f_idx;
                vector<int> unique_values;
                set<int> duplicates;
                for (auto v: CGAL::vertices_around_face(mesh.halfedge(n), mesh)) {
                    Point3 p = mesh.point(v);
                    int idx = v.idx();
                    n_pts.push_back(p);
                    pts.push_back(idx);
                    map_pts[idx] = p;
                }

                Point3 n_v1 = n_pts[0];
                Point3 n_v2 = n_pts[1];
                Point3 n_v3 = n_pts[2];

                sort(pts.begin(), pts.end());
                set<int> distinct(pts.begin(), pts.end());
                set_difference(pts.begin(), pts.end(), distinct.begin(), distinct.end(), inserter(duplicates, duplicates.end()));
                vector<int> common(duplicates.begin(), duplicates.end());
                auto it = unique(pts.begin(), pts.end());
                pts.erase(it, pts.end());
                unique_values = pts;
                set<int> tmp;
                set_difference(unique_values.begin(), unique_values.end(), common.begin(), common.end(), inserter(tmp, tmp.end()));
                vector<int> uncommon(tmp.begin(), tmp.end());

                if (common.size() == 2 && uncommon.size() == 2) {
                    if (sqrt(CGAL::squared_area(f_v1, f_v2, f_v3)) > 0 && sqrt(CGAL::squared_area(n_v1, n_v2, n_v3)) > 0 ) {
                        Point3 P = map_pts[common[0]];
                        Point3 a0 = map_pts[common[1]];
                        Point3 a1 = map_pts[uncommon[0]];
                        Point3 a2 = map_pts[uncommon[1]];

                        Vector3 b0(a0.x() - P.x(), a0.y() - P.y(), a0.z() - P.z());
                        Vector3 b1(a1.x() - P.x(), a1.y() - P.y(), a1.z() - P.z());
                        Vector3 b2(a2.x() - P.x(), a2.y() - P.y(), a2.z() - P.z());

                        Vector3 A = CGAL::cross_product(b0, b1);
                        Vector3 B = CGAL::cross_product(b0, b2);

                        double phi1 = (acos(CGAL::scalar_product(CGAL::cross_product(b0, b1), CGAL::cross_product(b0, b2)) /
                                            (pow(A.squared_length(), 0.5) * pow(B.squared_length(), 0.5)))) * (180 / M_PI);
                        double phi2 = 360 - abs(phi1);

                        if (abs(phi1) <= threshold && phi1 != 0) {
                            pair<pair<int, int>, double> pair;
                            pair = make_pair(make_pair(n.idx(), f.idx()), phi1);
                            auto it = find(sharp_angles.begin(), sharp_angles.end(), pair);
                            if (!(it != sharp_angles.end())) {
                                sharp_angles.push_back(make_pair(make_pair(f.idx(), n.idx()), phi1));
                            }
                        }

                        if (phi2 <= threshold && phi2 != 0) {
                            pair<pair<int, int>, double> pair;
                            pair = make_pair(make_pair(n.idx(), f.idx()), phi2);
                            auto it = find(sharp_angles.begin(), sharp_angles.end(), pair);
                            if (!(it != sharp_angles.end())) {
                                sharp_angles.push_back(make_pair(make_pair(f.idx(), n.idx()), phi2));
                            }
                        }
                    }
                }
            }
        }
    }

    //create an OBJ and TXT file containing the sharp angles

    //OBJ file
    std::ofstream ofile(OBJfilename);

    map<int, vector<int>> f;
    int fid = 0;
    int vid = 0;

    for (auto &pair: sharp_angles){
        fid += 1;
        vector<int> pts1;
        for (auto &v: store_faces[pair.first.first].v){
            vid += 1;
            pts1.push_back(vid);
            ofile << std::setprecision(5) << std::fixed << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
        }
        f[fid] = pts1;

        fid += 1;
        vector<int> pts2;
        for (auto &v: store_faces[pair.first.second].v){
            vid += 1;
            pts2.push_back(vid);
            ofile << std::setprecision(5) << std::fixed << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
        }
        f[fid] = pts2;
    }

    for (auto &ff: f){
        ofile << "f " << ff.second[0] << " " << ff.second[1] << " " << ff.second[2] << endl;
    }

    ofile.close();

    //TXT file
    ofstream tfile(TXTfilename);

    tfile << fixed << "n " << sharp_angles.size() << endl;

    double SP = 0;
    double n = 0;
    double A = 0;
    double theta = 0;
    for (auto &pair: sharp_angles){

        double angle = pair.second;

        Point3 f1_v1 = store_faces[pair.first.first].v[0];
        Point3 f1_v2 = store_faces[pair.first.first].v[1];
        Point3 f1_v3 = store_faces[pair.first.first].v[2];

        Point3 f2_v1 = store_faces[pair.first.second].v[0];
        Point3 f2_v2 = store_faces[pair.first.second].v[1];
        Point3 f2_v3 = store_faces[pair.first.second].v[2];

        tfile << fixed << "f pair" << endl;
        tfile << fixed << "     f1" << endl;
        tfile << setprecision(5) << fixed << "          v1 " << f1_v1.x() << " " << f1_v1.y() << " " << f1_v1.z()<< endl;
        tfile << setprecision(5) << fixed << "          v2 " << f1_v2.x() << " " << f1_v2.y() << " " << f1_v2.z()<< endl;
        tfile << setprecision(5) << fixed << "          v3 " << f1_v3.x() << " " << f1_v3.y() << " " << f1_v3.z()<< endl;
        tfile << fixed << "     f2" << endl;
        tfile << setprecision(5) << fixed << "          v1 " << f2_v1.x() << " " << f2_v1.y() << " " << f2_v1.z()<< endl;
        tfile << setprecision(5) << fixed << "          v2 " << f2_v2.x() << " " << f2_v2.y() << " " << f2_v2.z()<< endl;
        tfile << setprecision(5) << fixed << "          v3 " << f2_v3.x() << " " << f2_v3.y() << " " << f2_v3.z()<< endl;
        tfile << setprecision(5) << fixed << "     θ    " << angle << endl;
    }

    tfile.close();

    return sharp_angles;
}


vector<pair<int, double>> identify_sliver_triangles(double threshold){

    vector<pair<int, double>> sliver_triangles;
    vector<int> check;

    for (map<int, Face>::iterator f=fBuildings.begin(); f!= fBuildings.end(); f++) {
            //compute the sliver parameter of each face f forming a building
            Point3 p1 = vertices[f->second.v[0]].p;
            Point3 p2 = vertices[f->second.v[1]].p;
            Point3 p3 = vertices[f->second.v[2]].p;

            double area = sqrt(CGAL::squared_area(p1, p2, p3));
            double per1 = sqrt(CGAL::squared_distance(p1, p2));
            double per2 = sqrt(CGAL::squared_distance(p2, p3));
            double per3 = sqrt(CGAL::squared_distance(p3, p1));
            double perimeter = per1 + per2 + per3;

            double sliver_parameter = (2*area)/perimeter;

            //if this sliver parameter is lower than threshold, the vertices forming f are not collinear, and
            //face is not already identified as a sliver,
            //store face in "sliver_triangles"
            if (sliver_parameter <= threshold && !CGAL::collinear(p1, p2, p3)){
                if (find(check.begin(), check.end(), f->second.id) == check.end()){
                    sliver_triangles.push_back(make_pair(f->second.id, sliver_parameter));
                    check.push_back(f->second.id);
                }
            }
    }

    return sliver_triangles;
}


map<int, topo> check_topological_relationship(double threshold, int ground_level, double z_ground){

    map<int, topo> topo_relationships;

    if (!fTerrain.empty()){
        //if there is a terrain, create an AABB tree that store terrain surfaces in 2D
        map<Point2, double> map;

        list<Triangle3> T;
        for (auto &ft: fTerrain){
            Point3 p0 (vi[ft.second.v[0]].p.x(), vi[ft.second.v[0]].p.y(), 0);
            Point3 p1 (vi[ft.second.v[1]].p.x(), vi[ft.second.v[1]].p.y(), 0);
            Point3 p2 (vi[ft.second.v[2]].p.x(), vi[ft.second.v[2]].p.y(), 0);

            double d0 = vi[ft.second.v[0]].p.z();
            double d1 = vi[ft.second.v[1]].p.z();
            double d2 = vi[ft.second.v[2]].p.z();

            map[Point2(p0.x(), p0.y())] = d0;
            map[Point2(p1.x(), p1.y())] = d1;
            map[Point2(p2.x(), p2.y())] = d2;

            Triangle3 tri2D(p0, p1, p2);

            T.push_back(tri2D);
        }
        Tree terrain_tree(T.begin(), T.end());

        int id;

        for (auto &b: buildings) {
            for (auto &gs: buildings[b.first].gs){
                //project each ground surface gs of each building to 2D
                Point3 p0 (vertices[gs.v[0]].p.x(), vertices[gs.v[0]].p.y(), 0);
                Point3 p1 (vertices[gs.v[1]].p.x(), vertices[gs.v[1]].p.y(), 0);
                Point3 p2 (vertices[gs.v[2]].p.x(), vertices[gs.v[2]].p.y(), 0);
                Triangle3 tri(p0, p1, p2);

                //identify intersections between terrain surfaces from the AABB tree and gs
                list<Primitive_id> primitives;
                terrain_tree.all_intersected_primitives(tri, std::back_inserter(primitives));

                for (auto &p: primitives){
                    //iterate over intersecting terrain surfaces and convert them back to 3D
                    double z0 = map[Point2(p->vertex(0).x(), p->vertex(0).y())];
                    double z1 = map[Point2(p->vertex(1).x(), p->vertex(1).y())];
                    double z2 = map[Point2(p->vertex(2).x(), p->vertex(2).y())];

                    Point3 pp0(p->vertex(0).x(), p->vertex(0).y(), z0);
                    Point3 pp1(p->vertex(1).x(), p->vertex(1).y(), z1);
                    Point3 pp2(p->vertex(2).x(), p->vertex(2).y(), z2);

                    Triangle3 tri3D(pp0, pp1, pp2);

                    //position of each vertex of gs relative to the intersecting terrain surface (in 2D) is evaluated
                    Plane3 PlaneTri3D(pp0, pp1, pp2);
                    CGAL::Oriented_side or0 = PlaneTri3D.oriented_side(p0);
                    CGAL::Oriented_side or1 = PlaneTri3D.oriented_side(p1);
                    CGAL::Oriented_side or2 = PlaneTri3D.oriented_side(p2);

                    Point3 p03D(p0.x(), p0.y(), vertices[gs.v[0]].p.z());
                    Point3 p13D(p1.x(), p1.y(), vertices[gs.v[1]].p.z());
                    Point3 p23D(p2.x(), p2.y(), vertices[gs.v[2]].p.z());

                    //distance between each vertex of gs and plane defined by intersecting terrain surface (in 2D) is computed
                    double d0 = sqrt(CGAL::squared_distance (p03D, PlaneTri3D));
                    double d1 = sqrt(CGAL::squared_distance (p13D, PlaneTri3D));
                    double d2 = sqrt(CGAL::squared_distance (p23D, PlaneTri3D));

                    //for each vertex of gs, if it is located above the intersecting terrain surface (in 2D),
                    //or the distance between them is greater than threshold
                    //a topological relationship error is considered
                    //vertex of gs and intersecting terrain surface (in 2D) are stored
                    if (d0 >= threshold || or0 == CGAL::ON_POSITIVE_SIDE){
                        id += 1;
                        topo_relationships[id].id_building = b.first;
                        topo_relationships[id].pBuilding = p03D;
                        topo_relationships[id].triTerrain = tri3D;

                        if (or0 == CGAL::ON_POSITIVE_SIDE){
                            topo_relationships[id].error = 1;
                        } else {
                            topo_relationships[id].warning = 1;
                        }

                    }

                    if (d1 >= threshold || or1 == CGAL::ON_POSITIVE_SIDE){
                        id += 1;
                        topo_relationships[id].id_building = b.first;
                        topo_relationships[id].pBuilding = p13D;
                        topo_relationships[id].triTerrain = tri3D;
                        topo_relationships[id].d = sqrt(CGAL::squared_distance(tri3D, p23D));

                        if (or1 == CGAL::ON_POSITIVE_SIDE){
                            topo_relationships[id].error = 1;
                        } else {
                            topo_relationships[id].warning = 1;
                        }

                    }

                    if (d2 >= threshold || or2 == CGAL::ON_POSITIVE_SIDE){
                        id += 1;
                        topo_relationships[id].id_building = b.first;
                        topo_relationships[id].pBuilding = p23D;
                        topo_relationships[id].triTerrain = tri3D;
                        topo_relationships[id].d = sqrt(CGAL::squared_distance(tri3D, p23D));

                        if (or2 == CGAL::ON_POSITIVE_SIDE){
                            topo_relationships[id].error = 1;
                        } else {
                            topo_relationships[id].warning = 1;
                        }
                    }
                }
            }
        }
    } else {

        //if there is no terrain
        double zmin;
        if (ground_level == 0){ //and no ground level is given by users, the minimum z-value in the model is selected as zmin
            auto z_minmax = minmax_element(city.vb.begin(), city.vb.end(), [](const Point3 &p1, const Point3 &p2) { return p1.z() < p2.z(); });
            zmin = z_minmax.first->z();
        } else { //otherwise, the ground level inserted by users is used as zmin
            zmin = z_ground;
        }

        int id = 0;
        for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {
            for (auto &gs: b->second.gs){
                //for each vertex of each ground surface of each building,
                //its position is evaluated with regarded to zmin
                //if this vertex is higher or have a distance a greater than threshold,
                //a topological relationships error is considered
                //this vertex and the same vertex but then with a z-value of zmin are stored
                for (auto &v: gs.v){
                    if (vertices[v].p.z() < zmin-threshold || vertices[v].p.z() > zmin){
                        id += 1;
                        topo_relationships[id].id_building = b->first;
                        topo_relationships[id].pBuilding = vertices[v].p;
                        topo_relationships[id].pz = Point3(vertices[v].p.x(), vertices[v].p.y(), zmin);

                        if (vertices[v].p.z() > zmin){
                            topo_relationships[id].error = 1;
                        }
                        else{
                            topo_relationships[id].warning = 1;
                        }
                    }
                }
            }
        }
    }

    return topo_relationships;
}


vector<pair<int, int>> identify_overlapping_buildings(){

    vector<pair<int,int>> overlapped;

    int i = 0;

    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {
        SMesh b_mesh = b->second.mesh;
        //SMesh b_gs_mesh = b->second.gs_mesh;
        for (map<int, Building>::iterator bb=buildings.begin(); bb!=buildings.end(); bb++) {
            if (b->first != bb->first){
                SMesh bb_mesh = bb->second.mesh;

                //compute intersection in 3D between building b and building bb
                //for this method, two conditions are required:
                // - buildings must bound a volume,
                // - buildings must not be self-intersecting
                if (CGAL::Polygon_mesh_processing::volume(b_mesh) > 0 && CGAL::Polygon_mesh_processing::volume(bb_mesh) > 0) {

                    if (PMP::do_intersect(b_mesh, bb_mesh)) {
                        if (!CGAL::Polygon_mesh_processing::does_self_intersect(b_mesh) && !CGAL::Polygon_mesh_processing::does_self_intersect(bb_mesh)) {
                            SMesh out;
                            PMP::corefine_and_compute_intersection(b_mesh, bb_mesh, out);

                            double volume = CGAL::Polygon_mesh_processing::volume(out);

                            //if the volume of this intersection is higher zero and buildings b and bb are not already stored as overlapping buildings,
                            //they are stored as overlapping buildings.
                            //here overlapping buildings are stored as a pair
                            //pair b and bb is the same as pair bb and b
                            if (volume > 0) {
                                i += 1;
                                pair<int, int> check;
                                check = make_pair(bb->first, b->first);
                                auto it = find(overlapped.begin(), overlapped.end(), check);
                                if (!(it != overlapped.end())) {
                                    overlapped.push_back(make_pair(b->first, bb->first));
                                }
                            }
                        } else{
                            //if the two conditions mentioned previously are not met,
                            //the shared areas in 2D between the ground surfaces of buildings b and bb are computed
                            //if this area is higher than 0, buildings b and bb are considered overlapping
                            double area = 0;
                            for (auto &f: b->second.gs){

                                Point3 p1 (vertices[f.v[0]].p.x(),vertices[f.v[0]].p.y(), 0);
                                Point3 p2 (vertices[f.v[1]].p.x(),vertices[f.v[1]].p.y(), 0);
                                Point3 p3 (vertices[f.v[2]].p.x(),vertices[f.v[2]].p.y(), 0);

                                Triangle3 tri1(p1, p2, p3);
                                for (auto &ff: bb->second.gs) {
                                    Point3 p4 (vertices[ff.v[0]].p.x(),vertices[ff.v[0]].p.y(), 0);
                                    Point3 p5 (vertices[ff.v[1]].p.x(),vertices[ff.v[1]].p.y(), 0);
                                    Point3 p6 (vertices[ff.v[2]].p.x(),vertices[ff.v[2]].p.y(), 0);

                                    Triangle3 tri2(p4, p5, p6);

                                    auto tri_in = CGAL::intersection(tri1, tri2);

                                    if (Triangle3* s = boost::get<Triangle3>(&*tri_in)){
                                        Triangle3 t = *s;
                                        Point3 v1 = t.vertex(0);
                                        Point3 v2 = t.vertex(1);
                                        Point3 v3 = t.vertex(2);

                                        area += sqrt(CGAL::squared_area(v1, v2, v3));
                                    }
                                }
                            }

                            if (area > 0){
                                i += 1;
                                pair<int, int> check;
                                check = make_pair(bb->first, b->first);
                                auto it = find(overlapped.begin(), overlapped.end(), check);
                                if (!(it != overlapped.end())) {
                                    overlapped.push_back(make_pair(b->first, bb->first));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return overlapped;
}
