#include "definitions.h"
#include "input.h"
#include "output.h"

void define_buildings(string input){

    save2obj_iBuildings(input);

    SMesh mesh;

    if(!PMP::IO::read_polygon_mesh(input, mesh) || !CGAL::is_triangle_mesh(mesh))
    {
        cerr << "Invalid input." << std::endl;
    }

    std::vector<SMesh> ccmeshes;

    PMP::split_connected_components(mesh, ccmeshes);

    int id = 0;
    int t = 1;
    int k = 1;

    for (auto &mesh: ccmeshes){

        id += 1;
        buildings[id].id = id;
        buildings[id].mesh = mesh;

        vector<int> pts;

        for (auto fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit){

                int fidx = fit->idx();
                int fid = fidx+t;

                fBuildings[fid].id = fid;
                fBuildings[fid].id_building = id;

                CGAL::Vertex_around_face_iterator<SMesh> vbegin, vend;
                for(boost::tie(vbegin, vend) = vertices_around_face(mesh.halfedge(*fit), mesh); vbegin != vend; ++vbegin){
                    Point3 v = mesh.point(*vbegin);

                    int vidx = vbegin->idx();
                    int vid = vidx+k;

                    vertices[vid].id = vid;
                    vertices[vid].p = v;
                    fBuildings[fid].v.push_back(vid);

                    pts.push_back(vid);
                }

                buildings[id].faces.push_back(fBuildings[fid]);
        }

        sort(pts.begin(), pts.end());
        auto it = unique(pts.begin(), pts.end());
        pts.erase(it, pts.end());

        for (auto &p: pts){
            buildings[id].pts.push_back(vertices[p].p);
        }

        t += mesh.num_faces();
        k += mesh.num_vertices();
    }
}

Vector3 normal_vector(Point3 a, Point3 b, Point3 c){

    Vector3 v1(b.x()-a.x(), b.y()-a.y(), b.z()-a.z());
    Vector3 v2(c.x()-a.x(), c.y()-a.y(), c.z()-a.z());

    Vector3 normal = CGAL::cross_product(v1, v2);
    normal = normal/sqrt(normal.squared_length());
    return normal;
}

void define_ground_surfaces(double z, double theta){

    for (map<int, Building>::iterator b=buildings.begin(); b!= buildings.end(); b++) {
        vector<Point3> lowest_faces;
        for (auto &f: b->second.faces) {

            if (vertices[f.v[0]].p.z() <= z && vertices[f.v[1]].p.z() <= z && vertices[f.v[2]].p.z() <= z){
                    Vector3 normal = normal_vector(vertices[f.v[0]].p, vertices[f.v[1]].p, vertices[f.v[2]].p);

                    double thetax = abs(atan(normal.x()/normal.z())*(180/M_PI));
                    double thetay = abs(atan(normal.y()/normal.z())*(180/M_PI));

                    if (thetax <= theta && thetay <= theta && normal.z() < 0 ){
                        f.gs = 1;
                        b->second.gs.push_back(f);
                        for (auto &p: f.v){
                            b->second.gs2D.push_back(Point2(vertices[p].p.x(), vertices[p].p.y()));
                        }
                    }

            } else {
                f.gs = 0;
            }
        }
    }
}

vector<size_t> convex_hull(vector<Point2> vertices) {
    vector<size_t> indices(vertices.size()), out;
    iota(indices.begin(), indices.end(),0);
    CGAL::convex_hull_2(indices.begin(), indices.end(), back_inserter(out),Convex_hull_traits_2(CGAL::make_property_map(vertices)));
    return out;
}

Point2 centroid(int b){

    double x = 0;
    double y = 0;
    for (auto &p: buildings[b].pts){
        x += p.x();
        y += p.y();
    }
    x = x/(double)buildings[b].pts.size();
    y = y/(double)buildings[b].pts.size();

    Point2 c(x, y);

    return c;
}


void identify_neighbours_building(int b) {

    vector<int> neighbours;
    vector <Point_2> vertices;
    vector <Point2> centroids;

    if (buildings[b].gs2D.size() != 0) {

        map<Point2, int> map_centroids;

        for (auto &bb: buildings) {
            Point2 c = centroid(bb.first);
            double d = sqrt(CGAL::squared_distance(c, centroid(b)));
            if (b != bb.first && d <= 1000) {
                if (bb.second.gs2D.size() != 0) {
                    Point_2 cbb = centroid(bb.first);
                    map_centroids[cbb] = bb.first;
                    centroids.push_back(cbb);
                }
            }
        }

        if (centroids.size() != 0) {
            Point2 cb = centroid(b);
            centroids.push_back(cb);

            VD vd(centroids.begin(), centroids.end());

            //cout << "check validity " << vd.is_valid() << endl;

            assert(vd.is_valid());

            VD::Locate_result lr = vd.locate(cb);

            if (Face_handle *f = boost::get<Face_handle>(&lr)) {
                Ccb_halfedge_circulator ec_start = (*f)->ccb();
                Ccb_halfedge_circulator ec = ec_start;
                do {
                    VD::Face_handle neighborFace = ec->twin()->face();
                    neighbours.push_back(map_centroids[neighborFace->dual()->point()]);
                    ++ec;
                } while (ec != ec_start);
            }
        }
    }
    sort(neighbours.begin(), neighbours.end());
    neighbours.erase(unique(neighbours.begin(), neighbours.end()), neighbours.end());

    buildings[b].neighbours = neighbours;
}


void identify_neighbours_buildings(){
    if (buildings.size() > 1){
        for (map<int, Building>::iterator b = buildings.begin(); b != buildings.end(); b++) {
            identify_neighbours_building(b->first);
        }
    }
}

void distances_for_histogram(string file_with_distances){

    identify_neighbours_buildings();

    ofstream tfile(file_with_distances);

    for (map<int, Building>::iterator b=buildings.begin(); b!= buildings.end(); b++) {

        for (auto &bb: b->second.neighbours) {
            vector<size_t> a_out = convex_hull(b->second.gs2D);
            vector<size_t> b_out = convex_hull(buildings[bb].gs2D);

            for (size_t ia: a_out){
                for (size_t ib: b_out){
                    double d = sqrt(CGAL::squared_distance(b->second.gs2D[ia], buildings[bb].gs2D[ib]));

                    tfile << setprecision(5) << fixed << d << endl;


                }
            }
        }
    }
    tfile.close();


//    ofstream ofile("../output/target_building.obj");
//
//    map<int, vector<int>> f;
//    int fid = 0;
//    int vid = 0;
//
//    for (auto &ff: buildings[1].faces){
//        fid += 1;
//        vector<int> pts;
//
//        for (auto &v: ff.v){
//            Point3 p = vertices[v].p;
//            vid += 1;
//            pts.push_back(vid);
//            ofile << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
//        }
//
//        f[fid] = pts;
//
//    }
//
//    for (auto &ff: f){
//        ofile << "f " << ff.second[0] << " " << ff.second[1] << " " << ff.second[2] << endl;
//    }
//
//    ofile.close();
//
//    string name = "../output/neigh.obj";
//    ofstream o2file(name);
//
//
//    int fid2 = 0;
//    int vid2 = 0;
//
//    for (auto &b: buildings[1].neighbours){
//
//        map<int, vector<int>> f;
//
//        for (auto &ff: buildings[b].faces){
//            fid2 += 1;
//            vector<int> pts;
//
//            for (auto &v: ff.v){
//                Point3 p = vertices[v].p;
//                vid2 += 1;
//                pts.push_back(vid2);
//                o2file << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
//            }
//
//            f[fid2] = pts;
//
//        }
//
//        for (auto &ff: f){
//            o2file << "f " << ff.second[0] << " " << ff.second[1] << " " << ff.second[2] << endl;
//        }
//
//
//    }
//    o2file.close();
}
