#include "definitions.h"
#include "input.h"
#include "output.h"

void define_buildings(string input){

    save2obj_iBuildings(input); //store buildings in separate OBJ file

    SMesh mesh;
    if(!PMP::IO::read_polygon_mesh(input, mesh) || !CGAL::is_triangle_mesh(mesh)) { //generate a mesh from file previously created
        cerr << "Invalid input." << endl;
    }
    vector<SMesh> ccmeshes;
    PMP::split_connected_components(mesh, ccmeshes); //identifies and splits connected components in mesh, and stores them in ccmeshes

    int id = 0;
    int t = 1;
    int k = 1;
    for (auto &mesh: ccmeshes){ //store each mesh as a building

        id += 1;
        buildings[id].id = id; //store id
        buildings[id].mesh = mesh; //store mesh

        vector<int> pts; //store vertex ids

        for (auto fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit){ //store faces in fBuildings

                int fidx = fit->idx();
                int fid = fidx+t;

                fBuildings[fid].id = fid;
                fBuildings[fid].id_building = id;

                CGAL::Vertex_around_face_iterator<SMesh> vbegin, vend;
                for(boost::tie(vbegin, vend) = vertices_around_face(mesh.halfedge(*fit), mesh); vbegin != vend; ++vbegin){ //store vertices in vertices
                    Point3 v = mesh.point(*vbegin);

                    int vidx = vbegin->idx();
                    int vid = vidx+k;

                    vertices[vid].id = vid;
                    vertices[vid].p = v;
                    fBuildings[fid].v.push_back(vid); //store vertices per face in fBuildings

                    pts.push_back(vid); //store vertex ids in pts
                }

                buildings[id].faces.push_back(fBuildings[fid]); //store faces per building
            }

            sort(pts.begin(), pts.end()); //remove duplicates in pts
            auto it = unique(pts.begin(), pts.end());
            pts.erase(it, pts.end());

            for (auto &p: pts){ //store vertices per building using pts
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

    //surfaces with all three vertices lower or equal than z and angle smaller than θ are selected

    for (map<int, Building>::iterator b=buildings.begin(); b!= buildings.end(); b++) {

        vector<Point3> lowest_faces;

        for (auto &f: b->second.faces) {

            if (vertices[f.v[0]].p.z() <= z && vertices[f.v[1]].p.z() <= z && vertices[f.v[2]].p.z() <= z){

                    Vector3 normal = normal_vector(vertices[f.v[0]].p, vertices[f.v[1]].p, vertices[f.v[2]].p);
                    double thetax = abs(atan(normal.x()/normal.z())*(180/M_PI));
                    double thetay = abs(atan(normal.y()/normal.z())*(180/M_PI));

                    if (thetax <= theta && thetay <= theta && normal.z() < 0 ){
                        f.gs = 1; //label face as a ground surface
                        b->second.gs.push_back(f); //store ground surfaces per building
                        for (auto &p: f.v){ //store 2D vertices of ground surfaces per building
                            b->second.gs2D.push_back(Point2(vertices[p].p.x(), vertices[p].p.y()));
                        }

                    }

            } else {
                f.gs = 0; //label face as not being a ground surface
            }
        }
    }
}

double building_height(Building b){

    //compute building height

    vector<double> z_pts; //store all z-values
    for (auto &v: b.pts){
        double z = v.z();
        z_pts.push_back(z);
    }

    auto z_minmax = minmax_element(z_pts.begin(), z_pts.end(),[](const double& p1, const double& p2) { return p1 < p2; });
    double height = *z_minmax.second-*z_minmax.first; //compute difference highest and lowest z-values

    return height;
}

vector<size_t> convex_hull(vector<Point2> vertices) {

    //generate convex hull from a set of 2D points

    vector<size_t> indices(vertices.size()), out;
    iota(indices.begin(), indices.end(),0);
    CGAL::convex_hull_2(indices.begin(), indices.end(), back_inserter(out),Convex_hull_traits_2(CGAL::make_property_map(vertices)));

    return out;
}

void distances_between_buildings(double threshold){

    //computes and stores distances between buildings
    //these distances are stored in separations for each building

    for (map<int, Building>::iterator b=buildings.begin(); b!= buildings.end(); b++) { //building A

        vector<pair<pair<Point2, Point2>, double>> separations;

        for (map<int, Building>::iterator bb=buildings.begin(); bb!= buildings.end(); bb++) { //building B

            //generate convex hulls by using the ground surface vertices of building A and B
            vector<size_t> a_out = convex_hull(b->second.gs2D);
            vector<size_t> b_out = convex_hull(bb->second.gs2D);

            for (size_t ia: a_out){

                for (size_t ib: b_out){

                    ////compute distance between a vertex from convex hull building A and one from convex hull building B
                    double d = sqrt(CGAL::squared_distance(b->second.gs2D[ia], bb->second.gs2D[ib]));

                    //create a pair with ids from building A and B
                    pair<Point2, Point2> ab = {b->second.gs2D[ia], bb->second.gs2D[ib]};

                    if (d >= threshold){
                        //if this distance is higher than threshold, the actual distance value is stored with pair formed previously
                        separations.push_back(make_pair(ab, d));
                    } else if (d < threshold) {
                        //otherwise, -1 is stored with this pair to indicate that this distance does not need to be taken into account
                        separations.push_back(make_pair(ab, -1));
                    }
                }
            }
        }
        b->second.separations = separations; //store separations per building
    }
}


double perimeter_building(vector<Point2> gsA) {

    //approximate building perimeters based on their ground surfaces in 2D

    list<Point> points; //create a list from points ground surfaces
    for (auto &p: gsA){
        points.push_back(p);
    }

    //create an alpha-shape with the list created previously
    Alpha_shape_2 A(points.begin(), points.end(), FT(10000),Alpha_shape_2::GENERAL);
    vector<Segment2> segments;
    *A.find_optimal_alpha(1);

    double perimeter = 0;

    //compute length of each edges forming the alpha shape and add to perimeter
    Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin(), end = A.alpha_shape_edges_end();
    for( ; it!=end; ++it) {
        perimeter += sqrt(A.segment(*it).squared_length());
    }

    perimeter = perimeter;

    return perimeter;
}

double maximum_dimension(Building b){

    //compute maximum dimension of building b

    //store each x, y, z values in separate vectors
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

    //compute max and min x, y, and z values
    auto x_minmax = minmax_element(x_pts.begin(), x_pts.end(),[](const double& p1, const double& p2) { return p1 < p2; });
    auto y_minmax = minmax_element(y_pts.begin(), y_pts.end(),[](const double& p1, const double& p2) { return p1 < p2; });
    auto z_minmax = minmax_element(z_pts.begin(), z_pts.end(),[](const double& p1, const double& p2) { return p1 < p2; });

    //Use these values to compute length, width, and height
    double length = *x_minmax.second-*x_minmax.first;
    double width = *y_minmax.second-*y_minmax.first;
    double height = *z_minmax.second-*z_minmax.first;

    //Select the highest value between length, width, and height
    vector<double> max_dimensions {length, width, height};
    auto dim_max = max_element(max_dimensions.begin(), max_dimensions.end());
    double max = *dim_max;

    return max;
}