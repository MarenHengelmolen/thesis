#include "definitions.h"
#include "input.h"
#include "output.h"
#include <omp.h>

void define_buildings(){

    const string input = "../output/buildings.obj";
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
                fBuildings[fid].type = 1;

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

        for (auto fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit){

            CGAL::Face_around_face_iterator<SMesh> fbegin, fend;
            for(boost::tie(fbegin, fend) = faces_around_face(mesh.halfedge(*fit), mesh); fbegin != fend; ++fbegin){
                vector<int> pts;
                vector<int> unique_values;
                set<int> duplicates;

                int fid = fit->idx()+t;
                int neigh_idx = fbegin->idx()+t;


                pts = fBuildings[fid].v;
                for (auto &v: fBuildings[neigh_idx].v){
                    pts.push_back(v);
                }

                sort(pts.begin(), pts.end());

                set<int> distinct(pts.begin(), pts.end());
                set_difference(pts.begin(), pts.end(), distinct.begin(), distinct.end(), inserter(duplicates, duplicates.end()));

                vector<int> common(duplicates.begin(), duplicates.end());

                auto it = unique(pts.begin(), pts.end());

                pts.erase(it, pts.end());

                unique_values = pts;
                set<int> tmp;
                set_difference(unique_values.begin(), unique_values.end(), common.begin(), common.end(),inserter(tmp, tmp.end()));

                vector<int> uncommon(tmp.begin(),tmp.end());

                if (common.size() == 2 && uncommon.size() == 2){
                    fBuildings[fid].neigh.push_back(neigh_idx);
                    fBuildings[fid].v_common[neigh_idx] = common;
                    fBuildings[fid].v_uncommon[neigh_idx] = uncommon;
                }
            }
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

void define_ground_surfaces(){

    for (map<int, Building>::iterator b=buildings.begin(); b!= buildings.end(); b++) {
        vector<Point3> lowest_faces;
        for (auto &f: b->second.faces) {
           //add 5 as a threshold
            if (vertices[f.v[0]].p.z() < 5 && vertices[f.v[1]].p.z() < 5 && vertices[f.v[2]].p.z() < 5){
                    Vector3 normal = normal_vector(vertices[f.v[0]].p, vertices[f.v[1]].p, vertices[f.v[2]].p);
                    if (pow(normal.z(), 2)/(pow(normal.x(), 2)*pow(normal.y(), 2)) >= 1 && normal.z() < 0 ){
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

