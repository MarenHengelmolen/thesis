#include "definitions.h"
#include "input.h"
#include "buildings.h"
#include <map>

////SHORT EDGE IDENTIFICATION
vector<pair<pair<Point3, Point3>, double>> identify_short_edges(double threshold){

    vector<pair<pair<Point3, Point3>, double>> short_edges;
    vector<Point3> checkA;
    vector<Point3> checkB;

    for (map<int, Face>::iterator f=fBuildings.begin(); f!= fBuildings.end(); f++) {
            int v1 = f->second.v[0];
            int v2 = f->second.v[1];
            int v3 = f->second.v[2];

            pair<Point3, Point3> e1 = make_pair(vertices[v1].p, vertices[v2].p);
            pair<Point3, Point3> e2 = make_pair(vertices[v2].p, vertices[v3].p);
            pair<Point3, Point3> e3 = make_pair(vertices[v3].p, vertices[v1].p);

            vector<pair<Point3, Point3>> edges{e1, e2, e3};

            for (auto &e: edges){

                double distance = sqrt(CGAL::squared_distance(e.first, e.second));
                if (distance <= threshold){

                    auto itA = find(checkA.begin(), checkA.end(), e.first);
                    auto itB = find(checkA.begin(), checkA.end(), e.second);

                    if (itA != checkA.end())
                    {
                        int index = itA - checkA.begin();
                        if (checkB[index] == e.second){
                            //pass
                        }
                        else {
                            checkA.push_back(e.first);
                            checkB.push_back(e.second);
                            short_edges.push_back(make_pair(make_pair(e.first, e.second), distance));
                        }
                    }

                    else if (itB != checkA.end())
                    {
                        int index = itB - checkA.begin();
                        if (checkB[index] == e.first){
                            //pass
                        }
                        else {
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

//SHARP ANGLE IDENTIFICATION
vector<pair<pair<int, int>, double>> identify_sharp_angles(double threshold){
    vector<pair<pair<int, int>, double>> sharp_angles;

    for (map<int, Face>::iterator f=fBuildings.begin(); f!= fBuildings.end(); f++) {
        for (auto &neigh: f->second.neigh){

            Point3 P = vertices[f->second.v_common[neigh][0]].p;
            Point3 a0 = vertices[f->second.v_common[neigh][1]].p;
            Point3 a1 = vertices[f->second.v_uncommon[neigh][0]].p;
            Point3 a2 = vertices[f->second.v_uncommon[neigh][1]].p;

            Vector3 b0 (a0.x()-P.x(), a0.y()-P.y(), a0.z()-P.z());
            Vector3 b1 (a1.x()-P.x(), a1.y()-P.y(), a1.z()-P.z());
            Vector3 b2 (a2.x()-P.x(), a2.y()-P.y(), a2.z()-P.z());

            //Must be replaced by sliver triangles
            Vector3 A = CGAL::cross_product(b0, b1);
            Vector3 B = CGAL::cross_product(b0, b2);
            Vector3 Anorm = A/pow(A.squared_length(), 0.5);
            Vector3 Bnorm = B/pow(B.squared_length(), 0.5);

            if ( A.squared_length() != 0 && B.squared_length() != 0 && Anorm != Bnorm){
                if (A.x() == B.x() && A.y() == B.y() && A.z() == B.z()){
                    //Pass
                } else {
                    double phi = (acos(CGAL::scalar_product(CGAL::cross_product(b0, b1), CGAL::cross_product(b0, b2)) /
                                           (pow(A.squared_length(), 0.5) * pow(B.squared_length(), 0.5)))) * (180 / M_PI);
                    if (0.01 < phi && phi <= threshold) {
                        pair<pair<int, int>, double> pair;
                        pair = make_pair(make_pair(neigh, f->second.id), phi);

                        auto it = find(sharp_angles.begin(), sharp_angles.end(), pair);
                        if (it != sharp_angles.end()){
                            //Pair is found
                        }

                        else {
                            sharp_angles.push_back(make_pair(make_pair(f->second.id, neigh), phi));
                        }
                    }
                }
            }
        }
    }

    return sharp_angles;
}

//SLIVER TRIANGLE IDENTIFICATION
vector<pair<int, double>> identify_sliver_triangles(double threshold){

    vector<pair<int, double>> sliver_triangles;
    vector<int> check;

    for (map<int, Face>::iterator f=fBuildings.begin(); f!= fBuildings.end(); f++) {

            Point3 p1 = vertices[f->second.v[0]].p;
            Point3 p2 = vertices[f->second.v[1]].p;
            Point3 p3 = vertices[f->second.v[2]].p;

            double area = sqrt(CGAL::squared_area(p1, p2, p3));
            double per1 = sqrt(CGAL::squared_distance(p1, p2));
            double per2 = sqrt(CGAL::squared_distance(p2, p3));
            double per3 = sqrt(CGAL::squared_distance(p3, p1));
            double perimeter = per1 + per2 + per3;

            double sliver_parameter = (2*area)/perimeter;

            if (sliver_parameter < 0.1){
                if (find(check.begin(), check.end(), f->second.id) != check.end()){
                    //Pass
                } else{
                    sliver_triangles.push_back(make_pair(f->second.id, sliver_parameter));
                    check.push_back(f->second.id);
                }
            }
    }

    return sliver_triangles;
}


map<int, topo> check_topological_relationship(double threshold){

    map<int, topo> topo_relationships;

    define_ground_surfaces();

    if (!fTerrain.empty()){
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
                Point3 p0 (vertices[gs.v[0]].p.x(), vertices[gs.v[0]].p.y(), 0);
                Point3 p1 (vertices[gs.v[1]].p.x(), vertices[gs.v[1]].p.y(), 0);
                Point3 p2 (vertices[gs.v[2]].p.x(), vertices[gs.v[2]].p.y(), 0);
                Triangle3 tri(p0, p1, p2);

                list<Primitive_id> primitives;
                terrain_tree.all_intersected_primitives(tri, std::back_inserter(primitives));

                for (auto &p: primitives){
                    double z0 = map[Point2(p->vertex(0).x(), p->vertex(0).y())];
                    double z1 = map[Point2(p->vertex(1).x(), p->vertex(1).y())];
                    double z2 = map[Point2(p->vertex(2).x(), p->vertex(2).y())];

                    Point3 pp0(p->vertex(0).x(), p->vertex(0).y(), z0);
                    Point3 pp1(p->vertex(1).x(), p->vertex(1).y(), z1);
                    Point3 pp2(p->vertex(2).x(), p->vertex(2).y(), z2);

                    Triangle3 tri3D(pp0, pp1, pp2);

                    double d0 = sqrt(CGAL::squared_distance (p0, tri3D));
                    double d1 = sqrt(CGAL::squared_distance (p1, tri3D));
                    double d2 = sqrt(CGAL::squared_distance (p2, tri3D));

                    if (d0 >= 0 && d0 < threshold){
                        id += 1;
                        topo_relationships[id].id_building = id;
                        topo_relationships[id].pBuilding = p0;
                        topo_relationships[id].triTerrain = tri3D;

                    } else if (d1 >= 0 && d1 < threshold){
                        id += 1;
                        topo_relationships[id].id_building = id;
                        topo_relationships[id].pBuilding = p1;
                        topo_relationships[id].triTerrain = tri3D;

                    } else if (d2 >= 0 && d2 < threshold){
                        id += 1;
                        topo_relationships[id].id_building = id;
                        topo_relationships[id].pBuilding = p2;
                        topo_relationships[id].triTerrain = tri3D;
                    }
                }
            }
        }
    }


    return topo_relationships;
}



