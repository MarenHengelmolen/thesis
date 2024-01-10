#include "definitions.h"
#include "input.h"
#include "geometry.h"
#include "RoI.h"

void store_rotated_model(string output_file){

    //creates an OBJ file of the rotated input model

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

    ofile << "g Buildings" << endl;

    for (auto &ff: f){
        ofile << "f " << ff.second[0] << " " << ff.second[1] << " " << ff.second[2] << endl;
    }


    map<int, vector<int>> ft;
    int ftid = 0;
    int vtid = 0;

    for (auto &ff: fTerrain){
        ftid += 1;
        vector<int> pts;
        for (auto &v: ff.second.v){
            Point3 p = vi[v].p;
            vtid += 1;
            pts.push_back(vtid);
            ofile << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
        }
        ft[ftid] = pts;
    }

    ofile << "g Terrain" << endl;

    for (auto &ff: ft){
        ofile << "f " << ff.second[0] << " " << ff.second[1] << " " << ff.second[2] << endl;
    }

    ofile.close();

}

void save2obj_iBuildings(string output_file){

    //creates OBJ file including buildings from input file

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



void computational_domain(string OBJfilename){

    //creates an OBJ file containing the computational domain

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
}

void shortedges(string OBJfilename, string TXTfilename, double threshold) {

    //creates an OBJ and TXT file containing short edges

    vector<pair<pair<Point3, Point3>, double>> short_edges = identify_short_edges(threshold);
    output.n_short_edges = short_edges.size();

    cout << "Number of short edges: " << short_edges.size() << endl;

    //OBJ
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

    //TXT
    std::ofstream tfile(TXTfilename);

    tfile << fixed << "n " << short_edges.size() << endl;
    int n = 0;
    double dist = 0;
    for (auto &pair: short_edges){

        Point3 p1 = pair.first.first;
        Point3 p2 = pair.first.second;
        double distance = pair.second;
        n+= 1;
        dist += distance;

        tfile << setprecision(5) << fixed << "      v1 " << p1.x() << " " << p1.y() << " " << p1.z() << endl;
        tfile << setprecision(5) << fixed << "      v2 " << p2.x() << " " << p2.y() << " " << p2.z() << endl;
        tfile << setprecision(5) << fixed << "      l " << distance << endl;
        tfile << "\n" << endl;

    }

    tfile.close();
}

void sharp_angles(string OBJfilename, string TXTfilename,  double threshold, string input) {

    //creates an OBJ and TXT file containing sharp angles

    vector<pair<pair<int, int>, double>> sharp_angles = identify_sharp_angles(threshold, input, OBJfilename, TXTfilename);

    cout << "Number of sharp angles: " << sharp_angles.size() << endl;

    output.n_sharp_angles = sharp_angles.size();

    double d = 0;
    for (auto &angle: sharp_angles){
        d += angle.second;
    }
}

void sliver_triangles(string OBJfilename, string TXTfilename,  double threshold) {

    //creates an OBJ and TXT file containing sliver triangles

    vector<pair<int, double>> sliver_triangles = identify_sliver_triangles(threshold);

    cout << "Number of sliver triangles: " << sliver_triangles.size() << endl;

    output.n_slivers = sliver_triangles.size();

    //Save to OBJ
    std::ofstream ofile(OBJfilename);

    map<int, vector<int>> f;
    int fid = 0;
    int vid = 0;

    double SP = 0;
    double n = 0;
    double A = 0;

    for (auto &pair: sliver_triangles){

        Point3 p1 = vertices[fBuildings[pair.first].v[0]].p;
        Point3 p2 = vertices[fBuildings[pair.first].v[1]].p;
        Point3 p3 = vertices[fBuildings[pair.first].v[2]].p;

        double area = sqrt(CGAL::squared_area(p1, p2, p3));
        double per1 = sqrt(CGAL::squared_distance(p1, p2));
        double per2 = sqrt(CGAL::squared_distance(p2, p3));
        double per3 = sqrt(CGAL::squared_distance(p3, p1));
        double perimeter = per1 + per2 + per3;

        double sliver_parameter = (2*area)/perimeter;

        n += 1;
        SP += sliver_parameter;
        A += area;

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
        tfile << "\n" << endl;
    }

    tfile.close();
}

void save2obj_gs(string output_file){

    //creates an OBJ including ground surfaces

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



void topological_relationships(string OBJfilename, string TXTfilename,  double threshold, int ground_level, double z_ground) {

    //creates an OBJ and TXT file containing topological errors

    map<int, topo> relationships = check_topological_relationship(threshold, ground_level, z_ground);

    cout << "Number of topological errors: " << relationships.size() << endl;

    std::ofstream ofile(OBJfilename);

    if (!relationships.empty() && !fTerrain.empty()){

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


    } else if (!relationships.empty() && fTerrain.empty()){

        map<int, vector<int>> l;
        int lid = 0;
        int vid = 0;

        for (map<int, topo>::iterator t=relationships.begin(); t!= relationships.end(); t++) {
            lid += 1;
            vector<int> pts;

            Point3 p1 = t->second.pBuilding;
            vid += 1;
            pts.push_back(vid);
            ofile << std::setprecision(5) << std::fixed << "v " << p1.x() << " " << p1.y() << " " << p1.z() << std::endl;

            Point3 p2 = t->second.pz;
            vid += 1;
            pts.push_back(vid);
            ofile << std::setprecision(5) << std::fixed << "v " << p2.x() << " " << p2.y() << " " << p2.z() << std::endl;

            l[lid] = pts;
        }

        for (auto &ll: l){
            ofile << "l " << ll.second[0] << " " << ll.second[1] << endl;
        }


    }

    ofile.close();

    int errors = 0;
    int warnings = 0;

    for (auto &topo: relationships){
        if (topo.second.error == 1){
            errors += 1;
        } else {
            warnings += 1;
        }
    }

    output.n_topo = relationships.size();
    output.n_topo_errors = errors;
    output.n_topo_warnings = warnings;

    //TXT
    ofstream tfile(TXTfilename);
    tfile << fixed << "n " << relationships.size() << endl;
    tfile << fixed << "\n        centroid (2D)" << endl;
    if (fTerrain.size() != 0){
        for (map<int, topo>::iterator t=relationships.begin(); t!= relationships.end(); t++) {
            tfile << fixed << "     vb " <<  t->second.pBuilding << endl;
            Triangle3 tri = t->second.triTerrain;
            tfile << fixed << "      t " <<  CGAL::centroid(tri) << endl;
            tfile << fixed << "      d " <<  t->second.d << endl;
            tfile << "\n" << endl;
        }
    } else {
        for (map<int, topo>::iterator t=relationships.begin(); t!= relationships.end(); t++) {
            tfile << fixed << "     vb " <<  t->second.pBuilding << endl;
            tfile << fixed << "     vt " << t->second.pz << endl;
            tfile << fixed << "      d " <<  t->second.d << endl;
            tfile << "\n" << endl;
        }
    }

    tfile.close();
}

Point2 centroid(int target_building){

    //computes centroid of a building

    vector<Face> faces = buildings[target_building].faces;
    double x=0, y=0, z=0;
    int n = buildings[target_building].faces.size()*3; //faces are triangles
    for (auto &f: faces){
        for (auto &v: f.v){
            x += vertices[v].p.x();
            y += vertices[v].p.y();
            z += vertices[v].p.z();
        }
    }

    Point3 centroid3d ((double)x/n, (double)y/n, (double)z/n);
    Point2 centroid2d (centroid3d.x(), centroid3d.y());
    return centroid2d;
}

void region_of_interest(string OBJfilename, string TXTfilename, bool target, double x, double y){

    //creates an OBJ and TXT file containing buildings in the Region of Interest

    vector<int> detailed_buildings = Buildings_in_RoI(target, x, y);
    output.n_RoI = detailed_buildings.size();
    output.per_RoI = (detailed_buildings.size()/(double)buildings.size())*100;

    //OBJ
    std::ofstream ofile(OBJfilename);

    map<int, vector<int>> f;
    int fid = 0;
    int vid = 0;

    for (auto &db: detailed_buildings) {
        for (auto &ff: buildings[db].faces){
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

    //TXT
    ofstream tfile(TXTfilename);

    tfile << fixed << "n " << detailed_buildings.size() << endl;
    tfile << fixed << "        b centroid (2D)" << endl;
    for (auto &db: detailed_buildings){
        Point2 c = centroid(db);

        tfile << setprecision(5) << fixed << "b " << c << endl;
    }

    tfile.close();
}


void region_of_interest_cylinder(string OBJfilename){

    //creates an OBJ and TXT file containing the Region of Interest

    ofstream ofile(OBJfilename);

    SMesh RoI = meshCFD.RoI;

    map<int, vector<int>> f;
    int fid = 0;
    int vid = 0;

    for (auto fit = RoI.faces_begin(); fit != RoI.faces_end(); ++fit){
        fid += 1;
        vector<int> pts;

        CGAL::Vertex_around_face_iterator<SMesh> vbegin, vend;
        for(boost::tie(vbegin, vend) = vertices_around_face(RoI.halfedge(*fit), RoI); vbegin != vend; ++vbegin) {
            Point3 v = RoI.point(*vbegin);
            vid += 1;
            pts.push_back(vid);
            ofile << std::setprecision(5) << std::fixed << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
        }

        f[fid] = pts;
    }

    for (auto &ff: f){
        ofile << "f " << ff.second[0] << " " << ff.second[1] << " " << ff.second[2] << endl;
    }

    ofile.close();
}


void buildingseparations(string OBJinvalid, string OBJseparations, string TXTfile) {

    //creates an OBJ and TXT file containing invalid building separations

    std::ofstream o1file(OBJinvalid);

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
    std::ofstream o2file(OBJseparations);

    map<int, vector<int>> l;
    int lid = 0;
    int vid2 = 0;
    int i = 0;

    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {
        if (b->second.invalid_separations.size() != 0) {
            i += 1;
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

    //TXT
    ofstream tfile(TXTfile);
    tfile << fixed << "Separations" << endl;
    tfile << fixed << "n " << i << endl;
    int n = 0;
    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {
        if (b->second.invalid_separations.size() != 0) {
            n += 1;
            for (auto &sep: b->second.invalid_separations){
                tfile << fixed << "       v1 " << sep.first.first << endl;
                tfile << fixed << "       v2 " << sep.first.second << endl;
                tfile << fixed << "       d  " << sep.second << endl;
                tfile << fixed << "\n" << sep.second << endl;
            }
        }
    }

    tfile.close();

}

void building_volume(string OBJfile, string TXTfile){

    //creates an OBJ and TXT file containing invalid building volumes

    //OBJ and TXT
    ofstream ofile(OBJfile);

    map<int, vector<int>> f;
    int fid = 0;
    int vid = 0;
    int i = 0;

    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {

        if (b->second.invalid_volume == 1) {
            i += 1;
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

    //TXT
    ofstream tfile(TXTfile);
    tfile << fixed << "n " << i << endl;
    tfile << fixed << "\n      b centroid (2D)" << endl;
    for (map<int, Building>::iterator b=buildings.begin(); b!=buildings.end(); b++) {
        if (b->second.invalid_volume == 1) {
            tfile << fixed << "      b " << centroid(b->first) << endl;
        }
    }
    tfile.close();
}

void overlapping_buildings(string OBJfile, string TXTfile) {

    //creates an OBJ and TXT file containing overlapping buildings

    vector<pair<int, int>> overlapped = identify_overlapping_buildings();
    vector<int> selected;

    std::ofstream ofile(OBJfile);

    map<int, vector<int>> f;
    int fid = 0;
    int vid = 0;

    for (auto &p: overlapped){
        if (!count(selected.begin(), selected.end(), p.first)) {
            selected.push_back(p.first);
        }
        if (!count(selected.begin(), selected.end(), p.second)) {
            selected.push_back(p.second);
        }
    }
    cout << "Number overlapping buildings " << selected.size() << endl;

    output.overlap = selected.size();

    for (auto &db: selected) {
        for (auto &ff: buildings[db].faces){
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

    ofstream tfile(TXTfile);

    tfile << fixed << "n " << overlapped.size() << endl;
    tfile << fixed << "      centroid (2D)" << endl;

    for (auto &p: overlapped){
        tfile << fixed << "      b1 " << centroid(p.first) << endl;
        tfile << fixed << "      b2 " << centroid(p.second) << endl;
        tfile << fixed << "\n" << endl;
    }

    tfile.close();

}


void output_for_UI(string TXTfilename){

    //creates a TXT file used to return results to users in UI

    ofstream tfile(TXTfilename);

    tfile << setprecision(5) << fixed << "valid                     " << output.valid << endl;
    tfile << setprecision(5) << fixed << "n_topo                    " << output.n_topo << endl;
    tfile << setprecision(5) << fixed << "n_topo_errors             " << output.n_topo_errors << endl;
    tfile << setprecision(5) << fixed << "n_topo_warnings           " << output.n_topo_warnings << endl;
    tfile << setprecision(5) << fixed << "n_sharp_angles            " << output.n_sharp_angles << endl;
    tfile << setprecision(5) << fixed << "n_short_edges             " << output.n_short_edges << endl;
    tfile << setprecision(5) << fixed << "n_slivers                 " << output.n_slivers << endl;
    tfile << setprecision(5) << fixed << "domain                    " << output.domain.xmin() << " " << output.domain.xmax() << " " << output.domain.ymin() << " " << output.domain.ymax() << " " << output.domain.zmin() << " " << output.domain.zmax() << endl;
    tfile << setprecision(5) << fixed << "domain_GR_noNmax          " << output.domain_GR_noNmax.xmin() << " " << output.domain_GR_noNmax.xmax() << " " << output.domain_GR_noNmax.ymin() << " " << output.domain_GR_noNmax.ymax() << " " << output.domain_GR_noNmax.zmin() << " " << output.domain_GR_noNmax.zmax() << endl;
    tfile << setprecision(5) << fixed << "domain_GR_Nmax            " << output.domain_GR_Nmax.xmin() << " " << output.domain_GR_Nmax.xmax() << " " << output.domain_GR_Nmax.ymin() << " " << output.domain_GR_Nmax.ymax() << " " << output.domain_GR_Nmax.zmin() << " " << output.domain_GR_Nmax.zmax() << endl;
    tfile << setprecision(5) << fixed << "cell                      " << output.cell.first << " " << output.cell.second << endl;
    tfile << setprecision(5) << fixed << "cell_GR_noNmax            " << output.cell_GR_noNmax.first << " " << output.cell_GR_noNmax.second << endl;
    tfile << setprecision(5) << fixed << "cell_GR_Nmax              " << output.cell_GR_Nmax.first << " " << output.cell_GR_Nmax.second << endl;
    tfile << setprecision(5) << fixed << "N                         " << output.N << endl;
    tfile << setprecision(5) << fixed << "N_GR_noNmax               " << output.N_GR_noNmax << endl;
    tfile << setprecision(5) << fixed << "N_GR_Nmax                 " << output.N_GR_Nmax << endl;
    tfile << setprecision(5) << fixed << "eh                        " << output.eh.first << " " << output.eh.second << endl;
    tfile << setprecision(5) << fixed << "eh_GR_noNmax              " << output.eh_GR_noNmax << endl;
    tfile << setprecision(5) << fixed << "eh_GR_Nmax                " << output.eh_GR_Nmax << endl;
    tfile << setprecision(5) << fixed << "n_RoI                     " << output.n_RoI << endl;
    tfile << setprecision(5) << fixed << "per_RoI                   " << output.per_RoI << endl;

    tfile << setprecision(5) << fixed << "BV                        " << output.BV << endl;
    tfile << setprecision(5) << fixed << "BV_valid                  " << output.BV_valid << endl;
    tfile << setprecision(5) << fixed << "BV_valid_GR_noNmax        " << output.BV_valid_GR_noNmax << endl;
    tfile << setprecision(5) << fixed << "BV_valid_GR_Nmax          " << output.BV_valid_GR_Nmax << endl;

    tfile << setprecision(5) << fixed << "S_valid                   " << output.S_valid << endl;
    tfile << setprecision(5) << fixed << "S_valid_GR_noNmax         " << output.S_valid_GR_noNmax << endl;
    tfile << setprecision(5) << fixed << "S_valid_GR_Nmax           " << output.S_valid_GR_Nmax << endl;

    tfile << setprecision(5) << fixed << "BR                        " << output.BR << endl;
    tfile << setprecision(5) << fixed << "BRL                       " << output.BRL << endl;
    tfile << setprecision(5) << fixed << "BRH                       " << output.BRH << endl;
    tfile << setprecision(5) << fixed << "hx                        " << output.hx<< endl;
    tfile << setprecision(5) << fixed << "hy                        " << output.hy << endl;
    tfile << setprecision(5) << fixed << "hz                        " << output.hz << endl;

    tfile << setprecision(5) << fixed << "BR_noNmax                 " << output.BR_noNmax << endl;
    tfile << setprecision(5) << fixed << "BRL_noNmax                " << output.BRL_noNmax << endl;
    tfile << setprecision(5) << fixed << "BRH_noNmax                " << output.BRH_noNmax << endl;
    tfile << setprecision(5) << fixed << "hx_noNmax                 " << output.hx_noNmax<< endl;
    tfile << setprecision(5) << fixed << "hy_noNmax                 " << output.hy_noNmax << endl;
    tfile << setprecision(5) << fixed << "hz_noNmax                 " << output.hz_noNmax << endl;

    tfile << setprecision(5) << fixed << "BR_Nmax                   " << output.BR_Nmax << endl;
    tfile << setprecision(5) << fixed << "BRL_Nmax                  " << output.BRL_Nmax << endl;
    tfile << setprecision(5) << fixed << "BRH_Nmax                  " << output.BRH_Nmax << endl;
    tfile << setprecision(5) << fixed << "hx_Nmax                   " << output.hx_Nmax<< endl;
    tfile << setprecision(5) << fixed << "hy_Nmax                   " << output.hy_Nmax << endl;
    tfile << setprecision(5) << fixed << "hz_Nmax                   " << output.hz_Nmax << endl;

    tfile << setprecision(5) << fixed << "dmax                      " << output.dmax << endl;

    tfile << setprecision(5) << fixed << "overlap                      " << output.overlap << endl;

    tfile.close();
};

