#include "definitions.h"
#include "input.h"
#include "cfd.h"
#include "geometry.h"

void save2obj_terrain(string output_file){

    ofstream ofile(output_file);

    map<int, vector<int>> f;
    int fid = 0;
    int vid = 0;

    for (auto &ff: fTerrain){
            fid += 1;
            vector<int> pts;
            for (auto &v: ff.second.v){
                Point3 p = vi[v].p;
                vid += 1;
                pts.push_back(vid);
                ofile << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
            }
            f[fid] = pts;
    }

    for (auto &ff: f){
        ofile << "f " << ff.second[0] << " " << ff.second[1] << " " << ff.second[2] << endl;
    }

    ofile.close();
}

void save2obj_iBuildings(string output_file){

    ofstream ofile(output_file);

    map<int, vector<int>> f;
    int fid = 0;
    int vid = 0;

    for (auto &ff: fiBuildings){
            fid += 1;
            vector<int> pts;

            for (auto &v: ff.second.v){
                Point3 p = vi[v].p;
                vid += 1;
                pts.push_back(vid);
                ofile << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
            }

            f[fid] = pts;

    }

    for (auto &ff: f){
        ofile << "f " << ff.second[0] << " " << ff.second[1] << " " << ff.second[2] << endl;
    }

    ofile.close();
}

void save2obj_building(string output_file, int b){

    ofstream ofile(output_file);

    map<int, vector<int>> f;
    int fid = 0;
    int vid = 0;

    for (auto &ff: buildings[b].faces){

        fid += 1;
        vector<int> pts1;
        for (auto &v: ff.v){
            Point3 p = vertices[v].p;
            vid += 1;
            pts1.push_back(vid);

            ofile << setprecision(5) << fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
        }
        f[fid] = pts1;
    }

    for (auto &ff: f){
        ofile << "f " << ff.second[0] << " " << ff.second[1] << " " << ff.second[2] << endl;
    }

    ofile.close();
}

void computational_domain(string OBJfilename, string TXTfilename, double wCell){

    double xmin = meshCFD.domain.xmin();
    double xmax = meshCFD.domain.xmax();
    double ymin = meshCFD.domain.ymin();
    double ymax = meshCFD.domain.ymax();
    double zmin = meshCFD.domain.zmin();
    double zmax = meshCFD.domain.zmax();

    Point3 p1 (xmin, ymin, zmin);
    Point3 p2 (xmin, ymin, zmax);
    Point3 p3 (xmin, ymax, zmin);
    Point3 p4 (xmin, ymax, zmax);

    Point3 p5 (xmax, ymin, zmin);
    Point3 p6 (xmax, ymin, zmax);
    Point3 p7 (xmax, ymax, zmin);
    Point3 p8 (xmax, ymax, zmax);

    vector<Point3> points = {p1, p2, p3, p4, p5, p6, p7, p8};

    //save to OBJ
    ofstream ofile(OBJfilename);

    for (auto &p: points){
        ofile << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    }

    ofile << "f " << 1 << " " << 7 << " " << 5 << endl;
    ofile << "f " << 1 << " " << 3 << " " << 7 << endl;
    ofile << "f " << 1 << " " << 4 << " " << 3 << endl;
    ofile << "f " << 1 << " " << 2 << " " << 4 << endl;
    ofile << "f " << 3 << " " << 8 << " " << 7 << endl;
    ofile << "f " << 3 << " " << 4 << " " << 8 << endl;
    ofile << "f " << 5 << " " << 7 << " " << 8 << endl;
    ofile << "f " << 5 << " " << 8 << " " << 6 << endl;
    ofile << "f " << 1 << " " << 5 << " " << 6 << endl;
    ofile << "f " << 1 << " " << 6 << " " << 2 << endl;
    ofile << "f " << 2 << " " << 6 << " " << 8 << endl;
    ofile << "f " << 2 << " " << 8 << " " << 4 << endl;

    ofile.close();

    //save to TXT
    ofstream tfile(TXTfilename);

    tfile << setprecision(5) << fixed << "xmin  " << xmin << endl;
    tfile << setprecision(5) << fixed << "xmax  " << xmax << endl;
    tfile << setprecision(5) << fixed << "ymin  " << ymin << endl;
    tfile << setprecision(5) << fixed << "ymax  " << ymax << endl;
    tfile << setprecision(5) << fixed << "zmin  " << zmin << endl;
    tfile << setprecision(5) << fixed << "zmax  " << zmax << endl;

    for (auto &p: points) {
        tfile << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    }

    tfile.close();
}

void shortedges(string OBJfilename, string TXTfilename, double threshold) {
    vector<pair<pair<Point3, Point3>, double>> short_edges = identify_short_edges(threshold);

    //Save to OBJ
    std::ofstream ofile(OBJfilename);

    map<int, vector<int>> l;
    int lid = 0;
    int vid = 0;

    for (auto &pair: short_edges) {
        lid += 1;
        vector<int> pts;

        Point3 p1 = pair.first.first;
        vid += 1;
        pts.push_back(vid);
        ofile << std::setprecision(5) << std::fixed << "v " << p1.x() << " " << p1.y() << " " << p1.z() << std::endl;

        Point3 p2 = pair.first.second;
        vid += 1;
        pts.push_back(vid);
        ofile << std::setprecision(5) << std::fixed << "v " << p2.x() << " " << p2.y() << " " << p2.z() << std::endl;

        l[lid] = pts;
    }

    for (auto &ll: l){
        ofile << "l " << ll.second[0] << " " << ll.second[1] << endl;
    }

    ofile.close();

    //Save to TXT
    std::ofstream tfile(TXTfilename);

    tfile << fixed << "n " << short_edges.size() << endl;

    for (auto &pair: short_edges){

        Point3 p1 = pair.first.first;
        Point3 p2 = pair.first.second;
        double distance = pair.second;

        tfile << setprecision(5) << fixed << "v1 " << p1.x() << " " << p1.y() << " " << p1.z() << endl;
        tfile << setprecision(5) << fixed << "v2 " << p2.x() << " " << p2.y() << " " << p2.z() << endl;
        tfile << setprecision(5) << fixed << "d " << distance << endl;

    }

    tfile.close();
}

void sharp_angles(string OBJfilename, string TXTfilename,  double threshold) {
    vector<pair<pair<int, int>, double>> sharp_angles = identify_sharp_angles(threshold);

    //Save to OBJ
    std::ofstream ofile(OBJfilename);

    map<int, vector<int>> f;
    int fid = 0;
    int vid = 0;

    for (auto &pair: sharp_angles){
        fid += 1;
        vector<int> pts1;
        for (auto &v: fBuildings[pair.first.first].v){
            Point3 p = vertices[v].p;
            vid += 1;
            pts1.push_back(vid);
            ofile << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
        }
        f[fid] = pts1;

        fid += 1;
        vector<int> pts2;
        for (auto &v: fBuildings[pair.first.second].v){
            Point3 p = vertices[v].p;
            vid += 1;
            pts2.push_back(vid);
            ofile << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
        }
        f[fid] = pts2;
    }

    for (auto &ff: f){
        ofile << "f " << ff.second[0] << " " << ff.second[1] << " " << ff.second[2] << endl;
    }

    ofile.close();

    //save to TXT
    ofstream tfile(TXTfilename);

    tfile << fixed << "n " << sharp_angles.size() << endl;

    for (auto &pair: sharp_angles){

        Face f1 = fBuildings[pair.first.first];
        Face f2 = fBuildings[pair.first.second];
        double angle = pair.second;

        tfile << fixed << "f pair" << endl;
        tfile << fixed << "     f1" << endl;
        tfile << setprecision(5) << fixed << "          v1 " << vertices[f1.v[0]].p.x() << " " << vertices[f1.v[0]].p.y() << " " << vertices[f1.v[0]].p.z()<< endl;
        tfile << setprecision(5) << fixed << "          v2 " << vertices[f1.v[1]].p.x() << " " << vertices[f1.v[1]].p.y() << " " << vertices[f1.v[1]].p.z()<< endl;
        tfile << setprecision(5) << fixed << "          v3 " << vertices[f1.v[2]].p.x() << " " << vertices[f1.v[2]].p.y() << " " << vertices[f1.v[2]].p.z()<< endl;
        tfile << fixed << "     f2" << endl;
        tfile << setprecision(5) << fixed << "          v1 " << vertices[f2.v[0]].p.x() << " " << vertices[f2.v[0]].p.y() << " " << vertices[f2.v[0]].p.z()<< endl;
        tfile << setprecision(5) << fixed << "          v2 " << vertices[f2.v[1]].p.x() << " " << vertices[f2.v[1]].p.y() << " " << vertices[f2.v[1]].p.z()<< endl;
        tfile << setprecision(5) << fixed << "          v3 " << vertices[f2.v[2]].p.x() << " " << vertices[f2.v[2]].p.y() << " " << vertices[f2.v[2]].p.z()<< endl;
        tfile << setprecision(5) << fixed << "     θ    " << angle << endl;
    }

    tfile.close();
}

void sliver_triangles(string OBJfilename, string TXTfilename,  double threshold) {
    vector<pair<int, double>> sliver_triangles = identify_sliver_triangles(threshold);

    //Save to OBJ
    std::ofstream ofile(OBJfilename);

    map<int, vector<int>> f;
    int fid = 0;
    int vid = 0;

    for (auto &pair: sliver_triangles){

        fid += 1;
        vector<int> pts;
        for (auto &v: fBuildings[pair.first].v){
            Point3 p = vertices[v].p;
            vid += 1;
            pts.push_back(vid);
            ofile << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
        }
        f[fid] = pts;

    }

    for (auto &ff: f){
        ofile << "f " << ff.second[0] << " " << ff.second[1] << " " << ff.second[2] << endl;
    }

    ofile.close();

    //save to TXT
    ofstream tfile(TXTfilename);

    tfile << fixed << "n " << sliver_triangles.size() << endl;

    for (auto &pair: sliver_triangles){

        Face f = fBuildings[pair.first];
        double sliver_parameter = pair.second;

        tfile << fixed << "f" << endl;
        tfile << setprecision(5) << fixed << "      v1 " << vertices[f.v[0]].p.x() << " " << vertices[f.v[0]].p.y() << " " << vertices[f.v[0]].p.z()<< endl;
        tfile << setprecision(5) << fixed << "      v2 " << vertices[f.v[1]].p.x() << " " << vertices[f.v[1]].p.y() << " " << vertices[f.v[1]].p.z()<< endl;
        tfile << setprecision(5) << fixed << "      v3 " << vertices[f.v[2]].p.x() << " " << vertices[f.v[2]].p.y() << " " << vertices[f.v[2]].p.z()<< endl;
        tfile << setprecision(5) << fixed << "      SP " << sliver_parameter << endl;
    }

    tfile.close();
}

//void distances_between_buildings(string OBJfilename, string TXTfilename, double threshold){
//
//    vector<pair<pair<Point2, Point2>, double>> small_distances = identify_small_distances(threshold);
//
//    //save to OBJ
//    ofstream ofile(OBJfilename);
//
//    map<int, vector<int>> l;
//    int lid = 0;
//    int vid = 0;
//
//    for (auto &pair: small_distances) {
//
//        lid += 1;
//        vector<int> pts;
//
//        Point2 p1 = pair.first.first;
//        vid += 1;
//        pts.push_back(vid);
//        ofile << setprecision(5) << fixed << "v " << p1.x() << " " << p1.y() << " " << 0 << endl;
//
//        Point2 p2 = pair.first.second;
//        vid += 1;
//        pts.push_back(vid);
//        ofile << setprecision(5) << fixed << "v " << p2.x() << " " << p2.y() << " " << 0 << endl;
//
//        l[lid] = pts;
//    }
//
//    for (auto &ll: l){
//        ofile << "l " << ll.second[0] << " " << ll.second[1] << endl;
//    }
//
//    ofile.close();
//
//    //save to TXT
//    ofstream tfile(TXTfilename);
//
//    tfile << fixed << "n " << small_distances.size() << endl;
//
//    for (auto &pair: small_distances){
//
//        Point2 p1 = pair.first.first;
//        Point2 p2 = pair.first.second;
//        double distance = pair.second;
//
//        tfile << setprecision(5) << fixed << "v1 " << p1.x() << " " << p1.y() << endl;
//        tfile << setprecision(5) << fixed << "v2 " << p2.x() << " " << p2.y() << endl;
//        tfile << setprecision(5) << fixed << "d " << distance << endl;
//
//    }
//
//    tfile.close();
//
//}
//
//

void save2obj_gs(string output_file){

    ofstream ofile(output_file);

    map<int, vector<int>> f;
    int fid = 0;
    int vid = 0;

    for (map<int, Building>::iterator b=buildings.begin(); b!= buildings.end(); b++) {
        for (auto &ff: b->second.gs) {
            fid += 1;
            vector<int> pts;
            for (auto &v: ff.v) {
                Point3 p = vertices[v].p;
                vid += 1;
                pts.push_back(vid);
                ofile << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
            }
            f[fid] = pts;
        }
    }

    for (auto &ff: f){
        ofile << "f " << ff.second[0] << " " << ff.second[1] << " " << ff.second[2] << endl;
    }

    ofile.close();
}



void topological_relationships(string OBJfilename, string TXTfilename,  double threshold) {
    map<int, topo> relationships = check_topological_relationship(threshold);

    if (!relationships.empty()){
        std::ofstream ofile(OBJfilename);

        map<int, vector<int>> f;
        int fid = 0;
        int vid = 0;

        for (map<int, topo>::iterator t=relationships.begin(); t!= relationships.end(); t++) {
            Point3 pb = t->second.pBuilding;
            ofile << std::setprecision(5) << std::fixed << "v " << pb.x() << " " << pb.y() << " " << pb.z() << std::endl;
            vid += 1;

            Triangle3 tri = t->second.triTerrain;
            fid += 1;
            vector<int> pts;
            for (int i=0; i<3; i++) {
                Point3 p = tri.vertex(i);
                vid += 1;
                pts.push_back(vid);
                ofile << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
            }
            f[fid] = pts;
        }

        for (auto &ff: f){
            ofile << "f " << ff.second[0] << " " << ff.second[1] << " " << ff.second[2] << endl;
        }

        ofile.close();
    }
    //Save to OBJ


//    //save to TXT
//    ofstream tfile(TXTfilename);
//
//    tfile << fixed << "n " << sharp_angles.size() << endl;
//
//    for (auto &pair: sharp_angles){
//
//        Face f1 = fBuildings[pair.first.first];
//        Face f2 = fBuildings[pair.first.second];
//        double angle = pair.second;
//
//        tfile << fixed << "f pair" << endl;
//        tfile << fixed << "     f1" << endl;
//        tfile << setprecision(5) << fixed << "          v1 " << vertices[f1.v[0]].p.x() << " " << vertices[f1.v[0]].p.y() << " " << vertices[f1.v[0]].p.z()<< endl;
//        tfile << setprecision(5) << fixed << "          v2 " << vertices[f1.v[1]].p.x() << " " << vertices[f1.v[1]].p.y() << " " << vertices[f1.v[1]].p.z()<< endl;
//        tfile << setprecision(5) << fixed << "          v3 " << vertices[f1.v[2]].p.x() << " " << vertices[f1.v[2]].p.y() << " " << vertices[f1.v[2]].p.z()<< endl;
//        tfile << fixed << "     f2" << endl;
//        tfile << setprecision(5) << fixed << "          v1 " << vertices[f2.v[0]].p.x() << " " << vertices[f2.v[0]].p.y() << " " << vertices[f2.v[0]].p.z()<< endl;
//        tfile << setprecision(5) << fixed << "          v2 " << vertices[f2.v[1]].p.x() << " " << vertices[f2.v[1]].p.y() << " " << vertices[f2.v[1]].p.z()<< endl;
//        tfile << setprecision(5) << fixed << "          v3 " << vertices[f2.v[2]].p.x() << " " << vertices[f2.v[2]].p.y() << " " << vertices[f2.v[2]].p.z()<< endl;
//        tfile << setprecision(5) << fixed << "     θ    " << angle << endl;
//    }
//
//    tfile.close();
}

//void region_of_interest(string OBJfilename, string TXTfilename, bool target, double x, double y){
//    vector<int> detailed_buildings = Buildings_in_RoI(target, x, y);
//
//    //Save to OBJ
//    std::ofstream ofile(OBJfilename);
//
//    map<int, vector<int>> f;
//    int fid = 0;
//    int vid = 0;
//
//    for (auto &db: detailed_buildings) {
//        for (auto &ff: buildings[db].faces){
//            fid += 1;
//            vector<int> pts;
//            for (auto &v: ff.v) {
//                Point3 p = vertices[v].p;
//                vid += 1;
//                pts.push_back(vid);
//                ofile << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
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
//
//    //save to TXT
//    ofstream tfile(TXTfilename);
//
//    tfile << fixed << "n " << detailed_buildings.size() << endl;
//
//    for (auto &db: detailed_buildings){
//        Point3 c = centroid_building(db);
//        tfile << setprecision(5) << fixed << "b " << c.x() << " " << c.y() << endl;
//    }
//
//    tfile.close();
//}


void buildingseparations(string OBJfilename1, string OBJfilename2) {


    std::ofstream o1file(OBJfilename1);

    map<int, vector<int>> f;
    int fid = 0;
    int vid1 = 0;

    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {
        if (b->second.invalid_separations.size() != 0) {

            for (auto &ff: b->second.faces){

                fid += 1;
                vector<int> pts1;
                for (auto &v: ff.v){
                    Point3 p = vertices[v].p;
                    vid1 += 1;
                    pts1.push_back(vid1);

                    o1file << setprecision(5) << fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
                }
                f[fid] = pts1;
            }
        }
    }

    for (auto &ff: f){
        o1file << "f " << ff.second[0] << " " << ff.second[1] << " " << ff.second[2] << endl;
    }

    o1file.close();

    //Save to OBJ
    std::ofstream o2file(OBJfilename2);

    map<int, vector<int>> l;
    int lid = 0;
    int vid2 = 0;

    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {
        if (b->second.invalid_separations.size() != 0) {

            for (auto &sp: b->second.invalid_separations) {
                lid += 1;
                vector<int> pts;

                Point2 p1 = sp.first.first;
                vid2 += 1;
                pts.push_back(vid2);
                o2file << std::setprecision(5) << std::fixed << "v " << p1.x() << " " << p1.y() << " " << 0 << std::endl;

                Point2 p2 = sp.first.second;
                vid2 += 1;
                pts.push_back(vid2);
                o2file << std::setprecision(5) << std::fixed << "v " << p2.x() << " " << p2.y() << " " << 0 << std::endl;

                l[lid] = pts;
            }
        }
    }


    for (auto &ll: l){
        o2file << "l " << ll.second[0] << " " << ll.second[1] << endl;
    }

    o2file.close();

}

void building_volume(string output_file){

    ofstream ofile(output_file);

    map<int, vector<int>> f;
    int fid = 0;
    int vid = 0;

    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {

        if (b->second.invalid_volume == 1) {
            for (auto &ff: b->second.faces) {
                fid += 1;
                vector<int> pts1;
                for (auto &v: ff.v) {
                    Point3 p = vertices[v].p;
                    vid += 1;
                    pts1.push_back(vid);

                    ofile << setprecision(5) << fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
                }
                f[fid] = pts1;
            }
        }
    }

    for (auto &ff: f){
        ofile << "f " << ff.second[0] << " " << ff.second[1] << " " << ff.second[2] << endl;
    }

    ofile.close();
}
