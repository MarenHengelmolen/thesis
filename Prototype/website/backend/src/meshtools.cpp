#include "definitions.h"
#include "input.h"
#include "buildings.h"

double round5(double d){
    d = round(d*100000)/100000;
    return d;
}

double n_cells_block(){

    //computes number of cells in background mesh (blockMesh)

    CGAL::Bbox_3 domain = meshCFD.domain;
    double V_domain = abs(domain.xmax()-domain.xmin())*abs(domain.ymax()-domain.ymin())*abs(domain.zmax()-domain.zmin());
    double V_cell = meshCFD.wcell*meshCFD.wcell*meshCFD.hcell;
    double N = (double)V_domain/(double)V_cell;

    return N;
}

map<int, pair<double, double>>  B2_ground_volume(){

    //computes volume under ground surfaces

    double N_refined = 0; //number of cells with the highest refinement level (needed to check whether there are enough cells for 10 cells per cube root building volume)

    //cell areas with different refinement levels
    double A_cell3 = pow(meshCFD.wcell/8, 2);
    double A_cell2 = pow(meshCFD.wcell/4, 2);
    double A_cell1 = pow(meshCFD.wcell/2, 2);
    double A_cell0 = pow(meshCFD.wcell, 2);

    //cell volumes with different refinement levels
    double V_cell3 = A_cell3*(meshCFD.hcell/8);
    double V_cell2 = A_cell2*(meshCFD.hcell/4);
    double V_cell1 = A_cell1*(meshCFD.hcell/2);
    double V_cell0 = A_cell0*meshCFD.hcell;

    double volume = 0;
    double v_box1 = 0;
    double v_box2 = 0;
    double v_domain = 0;
    double v_domain_first = 0;
    double zmin = city.zmin_terrain;

    //split areas of refinement boxes and domain into triangles
    //add these triangles to AABB tree

    list<Triangle3> T1;
    list<Triangle3> T2;
    list<Triangle3> Td;

    CGAL::Bbox_3 box1 = meshCFD.box1;
    CGAL::Bbox_3 box2 = meshCFD.box2;
    CGAL::Bbox_3 domain = meshCFD.domain;

    Point3 a1 (box1.xmin(), box1.ymin(), 0);
    Point3 b1 (box1.xmin(), box1.ymax(), 0);
    Point3 c1 (box1.xmax(), box1.ymax(), 0);
    Point3 d1 (box1.xmax(), box1.ymin(), 0);

    Point3 a2 (box2.xmin(), box2.ymin(), 0);
    Point3 b2 (box2.xmin(), box2.ymax(), 0);
    Point3 c2 (box2.xmax(), box2.ymax(), 0);
    Point3 d2 (box2.xmax(), box2.ymin(), 0);

    Point3 ad (domain.xmin(), domain.ymin(), 0);
    Point3 bd (domain.xmin(), domain.ymax(), 0);
    Point3 cd (domain.xmax(), domain.ymax(), 0);
    Point3 dd (domain.xmax(), domain.ymin(), 0);

    T1.push_back(Triangle3(a1, b1, c1));
    T1.push_back(Triangle3(c1, d1, a1));
    T2.push_back(Triangle3(a2, b2, c2));
    T2.push_back(Triangle3(c2, d2, a2));
    Td.push_back(Triangle3(ad, bd, cd));
    Td.push_back(Triangle3(cd, dd, ad));

    Tree tree1(T1.begin(), T1.end());
    Tree tree2(T2.begin(), T2.end());
    Tree treed(Td.begin(), Td.end());

    //create plane passing through the bottom boundary of the domain
    Vector3 VPlane(0, 0, 1);
    Plane3 Plane(Point3(domain.xmin(), domain.ymin(), zmin), VPlane);

    double n1=0, n2=0, nD=0;
    double v1=0, v2=0, vD=0;

    for (auto &&ff: fTerrain) {

        //for each terrain face,
        //compute its centroid and use this point to approximate its distance with the plane previously defined,
        //multiply its area with this distance to approximate the volume under this face

        Point3 p1 = vi[ff.second.v[0]].p;
        Point3 p2 = vi[ff.second.v[1]].p;
        Point3 p3 = vi[ff.second.v[2]].p;
        Triangle3 tri2D(Point3(p1.x(), p1.y(), 0), Point3(p2.x(), p2.y(), 0), Point3(p3.x(), p3.y(), 0));

        Triangle3 t(p1, p2, p3);
        Point3 c = centroid(t);
        double d = sqrt(CGAL::squared_distance(Plane, c));
        double v = sqrt(t.squared_area()) * d;
        volume += v;

        double area = sqrt(t.squared_area());

        //identify with which area the face is intersecting in 2D: one of the refinement boxes or the remaining area of the domain

        //check if all vertices from the terrain surface are intersecting with the domain
        //compute number of cells per face, which depends on the area in which they are located
        //more cell layers with different refinement levels are needed in refinement box 2 than refinement box 1
        if (treed.do_intersect(tri2D.vertex(0)) && treed.do_intersect(tri2D.vertex(1)) && treed.do_intersect(tri2D.vertex(2))){

            if (treed.do_intersect(tri2D)) {
                v_domain_first += v;

                double na = (area/(double)A_cell3)*4;
                N_refined += na;
                double nb = (area/(double)A_cell2)*4;
                double va = na*V_cell3+nb*V_cell2;

                double nc = (area/(double)A_cell1)*4;
                double vc = nc*V_cell1;

                if (tree2.do_intersect(tri2D)){
                    ff.second.domain = 0;

                    if (tree1.do_intersect(tri2D)){
                        ff.second.box1 = 1;
                        ff.second.box2 = 0;
                        v_box1 += v;
                        n1 += (na+nb);
                        v1 += va;
                    } else {
                        ff.second.box2 = 1;
                        ff.second.box1 = 0;
                        v_box2 += v;
                        n2 += (na+nb+nc);
                        v2 += (va+vc);
                    }
                }else{
                    ff.second.box2 = 0;
                    ff.second.domain = 1;
                    v_domain += v;
                    double nd = (area/A_cell0)*4;
                    double vd = nd*V_cell0;
                    nD += (na+nb+nc+nd);
                    vD += va+vc+nd*vd;
                }
            }
        }
    }

    map<int, pair<double, double>> nv;

    nv[1] = make_pair(n1, v1); //number and total volume of cells in refinement box 1
    nv[2] = make_pair(n2, v2); // " in refinement box 2
    nv[3] = make_pair(nD, vD); // " in refinement box 3

    //store total volume under terrain surfaces per refinement box and remaining area of the domain
    meshCFD.v_box1 = v_box1;
    meshCFD.v_box2 = v_box2;
    meshCFD.v_domain = v_domain;

    //store number of cells with the highest refinement level
    meshCFD.N_refined += N_refined;

    return nv;

}


map<int, pair<double, double>>  B3_ground_volume(){

    //same algorithm as the previous one, but then with 3 refinement boxes

    double N_refined = 0;
    double A_cell4 = pow(meshCFD.wcell/16, 2);
    double A_cell3 = pow(meshCFD.wcell/8, 2);
    double A_cell2 = pow(meshCFD.wcell/4, 2);
    double A_cell1 = pow(meshCFD.wcell/2, 2);
    double A_cell0 = pow(meshCFD.wcell, 2);

    double V_cell4 = A_cell3*(meshCFD.hcell/16);
    double V_cell3 = A_cell3*(meshCFD.hcell/8);
    double V_cell2 = A_cell2*(meshCFD.hcell/4);
    double V_cell1 = A_cell1*(meshCFD.hcell/2);
    double V_cell0 = A_cell0*meshCFD.hcell;

    double volume = 0;
    double v_box1 = 0;
    double v_box2 = 0;
    double v_box3 = 0;
    double v_domain = 0;
    double v_domain_first = 0;
    double zmin = city.zmin_terrain;

    list<Triangle3> T1;
    list<Triangle3> T2;
    list<Triangle3> T3;
    list<Triangle3> Td;

    CGAL::Bbox_3 box1 = meshCFD.box1;
    CGAL::Bbox_3 box2 = meshCFD.box2;
    CGAL::Bbox_3 box3 = meshCFD.box3;
    CGAL::Bbox_3 domain = meshCFD.domain;

    Point3 a1 (box1.xmin(), box1.ymin(), 0);
    Point3 b1 (box1.xmin(), box1.ymax(), 0);
    Point3 c1 (box1.xmax(), box1.ymax(), 0);
    Point3 d1 (box1.xmax(), box1.ymin(), 0);

    Point3 a2 (box2.xmin(), box2.ymin(), 0);
    Point3 b2 (box2.xmin(), box2.ymax(), 0);
    Point3 c2 (box2.xmax(), box2.ymax(), 0);
    Point3 d2 (box2.xmax(), box2.ymin(), 0);

    Point3 a3 (box3.xmin(), box3.ymin(), 0);
    Point3 b3 (box3.xmin(), box3.ymax(), 0);
    Point3 c3 (box3.xmax(), box3.ymax(), 0);
    Point3 d3 (box3.xmax(), box3.ymin(), 0);

    Point3 ad (domain.xmin(), domain.ymin(), 0);
    Point3 bd (domain.xmin(), domain.ymax(), 0);
    Point3 cd (domain.xmax(), domain.ymax(), 0);
    Point3 dd (domain.xmax(), domain.ymin(), 0);

    T1.push_back(Triangle3(a1, b1, c1));
    T1.push_back(Triangle3(c1, d1, a1));
    T2.push_back(Triangle3(a2, b2, c2));
    T2.push_back(Triangle3(c2, d2, a2));
    T3.push_back(Triangle3(a3, b3, c3));
    T3.push_back(Triangle3(c3, d3, a3));
    Td.push_back(Triangle3(ad, bd, cd));
    Td.push_back(Triangle3(cd, dd, ad));

    Tree tree1(T1.begin(), T1.end());
    Tree tree2(T2.begin(), T2.end());
    Tree tree3(T3.begin(), T3.end());
    Tree treed(Td.begin(), Td.end());

    Vector3 VPlane(0, 0, 1);
    Plane3 Plane(Point3(domain.xmin(), domain.ymin(), zmin), VPlane);

    double n1=0, n2=0, n3=0, nD=0;
    double v1=0, v2=0, v3=0, vD=0;

    for (auto &&ff: fTerrain) {
        Point3 p1 = vi[ff.second.v[0]].p;
        Point3 p2 = vi[ff.second.v[1]].p;
        Point3 p3 = vi[ff.second.v[2]].p;
        Triangle3 tri2D(Point3(p1.x(), p1.y(), 0), Point3(p2.x(), p2.y(), 0), Point3(p3.x(), p3.y(), 0));

        Triangle3 t(p1, p2, p3);
        Point3 c = centroid(t);
        double d = sqrt(CGAL::squared_distance(Plane, c));
        double v = sqrt(t.squared_area()) * d;
        volume += v;

        double area = sqrt(t.squared_area());

        if (treed.do_intersect(tri2D.vertex(0)) && treed.do_intersect(tri2D.vertex(1)) && treed.do_intersect(tri2D.vertex(2))){// p1.x() <= domain.xmax() && p1.x() >= domain.xmin() && p1.y() <= domain.ymax() && p1.y() >= domain.ymin() && p2.x() <= domain.xmax() && p2.x() >= domain.xmin() && p2.y() <= domain.ymax() && p2.y() >= domain.ymin() && p3.x() <= domain.xmax() && p3.x() >= domain.xmin() && p3.y() <= domain.ymax() && p3.y() >= domain.ymin()){

            if (treed.do_intersect(tri2D)) {
                v_domain_first += v;

                double na = (area/(double)A_cell4)*4;
                N_refined += na;
                double nb = (area/(double)A_cell3)*4;
                double va = na*V_cell4+nb*V_cell4;

                double nc = (area/(double)A_cell2)*4;
                double vc = nc*V_cell2;

                double nz = (area/(double)A_cell1)*4;
                double vz = nz*V_cell1;

                if (tree3.do_intersect(tri2D)){
                    ff.second.domain = 0;

                    if (tree2.do_intersect(tri2D)){
                        ff.second.box3 = 0;

                        if (tree1.do_intersect(tri2D)){
                            ff.second.box1 = 1;
                            ff.second.box2 = 0;
                            v_box1 += v;
                            n1 += (na+nb);
                            v1 += va;
                        }else{
                            ff.second.box1 = 0;
                            ff.second.box2 = 1;
                            v_box2 += v;
                            n2 += (na+nb+nc);
                            v2 += (va+vc);
                        }
                    }else{
                        ff.second.box3 = 1;
                        ff.second.box2 = 0;
                        ff.second.box1 = 0;

                        v_box3 += v;
                        n3 += (na+nb+nc+nz);
                        v3 += (va+vc+vz);
                    }
                }else{
                    ff.second.domain = 1;
                    ff.second.box3 = 0;
                    ff.second.box2 = 0;
                    ff.second.box1 = 0;

                    v_domain += v;
                    double nd = (area/A_cell0)*4;
                    double vd = nd*V_cell0;
                    nD += (na+nb+nc+nz+nd);
                    vD += va+vc+vz+nd*vd;
                }
            }
        }
    }

    map<int, pair<double, double>> nv;

    nv[1] = make_pair(n1, v1);
    nv[2] = make_pair(n2, v2);
    nv[3] = make_pair(n3, v3);
    nv[3] = make_pair(nD, vD);

    meshCFD.v_box1 = v_box1;
    meshCFD.v_box2 = v_box2;
    meshCFD.v_box3 = v_box3;
    meshCFD.v_domain = v_domain;

    meshCFD.N_refined += N_refined;

    return nv;

}

double B2_n_cells_snappy(int gr, double height) {

    //approximates the number of cells in mesh with two refinement boxes

    meshCFD.N_refined = 0;

    CGAL::Bbox_3 box1 = meshCFD.box1;
    CGAL::Bbox_3 box2 = meshCFD.box2;
    CGAL::Bbox_3 domain = meshCFD.domain;

    double N_refined = 0; //number of cells with the highest refinement level (needed to check whether there are enough cells for 10 cells per cube root building volume)

    double volB = 0; //total volume of buildings
    double volG = 0; //total volume of cells around buildings
    double nG = 0; //number of cells around buildings

    //compute number of cells around buildings
    for (map<int, Building>::iterator b = buildings.begin(); b != buildings.end(); b++) {
        SMesh building = b->second.mesh;
        double vol = CGAL::Polygon_mesh_processing::volume(building);
        volB += vol;

        for (auto &f: b->second.faces){
            if (f.gs == 0){
                Triangle3 tri(vertices[f.v[0]].p, vertices[f.v[1]].p, vertices[f.v[2]].p);
                double area = sqrt(tri.squared_area());
                double n = (area/(meshCFD.wcell/8*meshCFD.hcell/8))*4; //number of cells per face with refinement level 3
                N_refined += n;
                double n_extra = (area/(meshCFD.wcell/4*meshCFD.hcell/4))*4; //number of cells per face with refinement level 2
                double v = n*pow(meshCFD.wcell/8, 2)*meshCFD.hcell/8; //volume of cells per face with refinement level 3
                double v_extra = n_extra*pow(meshCFD.wcell/4, 2)*meshCFD.hcell/4; //volume of cells per face with refinement level 2
                volG += v ;
                volG += v_extra;
                nG += n;
                nG += n_extra;
            }
        }
    }

    //cells area with different refinement levels
    double A_cell3 = pow(meshCFD.wcell/8, 2);
    double A_cell2 = pow(meshCFD.wcell/4, 2);
    double A_cell1 = pow(meshCFD.wcell/2, 2);
    double A_cell0 = pow(meshCFD.wcell, 2);

    //cells volumes with different refinement levels
    double V_cell3 = A_cell3*(meshCFD.hcell/8);
    double V_cell2 = A_cell2*(meshCFD.hcell/4);
    double V_cell1 = A_cell1*(meshCFD.hcell/2);
    double V_cell0 = A_cell0*meshCFD.hcell;

    //volumes per refinement boxes and remaining area of the domain
    double vol1 = ((box1.xmax()-box1.xmin())*(box1.ymax()-box1.ymin())*(box1.zmax()-box1.zmin()))-(volB+volG);
    double vol2 = ((box2.xmax()-box2.xmin())*(box2.ymax()-box2.ymin())*(box2.zmax()-box2.zmin()))-(vol1+volB+volG);
    double volD = ((domain.xmax()-domain.xmin())*(domain.ymax()-domain.ymin())*(domain.zmax()-domain.zmin()))-(vol1+volB+volG+vol2);

    double N = 0;

    if (gr == 0){ //without ground refinement

        double n1;
        double n2;
        double nDomain;

        if (fTerrain.size() == 0){ //without terrain

            //compute number of cells per refinement box and remaining area of the domain
            n1 = vol1/(pow(meshCFD.wcell/4, 2)*meshCFD.hcell/4);
            n2 = vol2/(pow(meshCFD.wcell/2, 2)*meshCFD.hcell/2);
            nDomain = volD/(pow(meshCFD.wcell, 2)*meshCFD.hcell);

            N = nG + n1 + n2 + nDomain;

        } else { //with terrain

            //volume of cells around terrain surfaces per refinement box and remaining area of the domain
            map<int, pair<double, double>> nv = B2_ground_volume();

            //volumes under terrain surfaces per refinement box and remaining area of the domain
            double v_box1 = meshCFD.v_box1;
            double v_box2 = meshCFD.v_box2;
            double v_domain = meshCFD.v_domain;

            //substract these two volumes from the corresponding area of the domain
            vol1 = vol1 - (v_box1+nv[1].second);
            vol2 = vol2 - (v_box2+nv[2].second);
            volD = volD - (v_domain+nv[3].second);

            //compute number of cells in each area of the domain
            n1 = vol1/(pow(meshCFD.wcell/4, 2)*meshCFD.hcell/4);
            n2 = vol2/(pow(meshCFD.wcell/2, 2)*meshCFD.hcell/2);
            nDomain = volD/(pow(meshCFD.wcell, 2)*meshCFD.hcell);

            double nT = nv[1].first+nv[2].first+nv[3].first;

            //sum up all number of cells to find the total number of cells
            N = nG + n1 + n2 + nDomain + nT;

        }

    } else { //with ground refinement

        double AB = 0;
        double perB = 0;
        //compute and sum up perimeters of each building based on alpha shapes
        //calculate total area of buildings
        for (map<int, Building>::iterator b = buildings.begin(); b != buildings.end(); b++) {
            if (!b->second.gs2D.empty()){
                for (auto &f: b->second.faces){
                   if (f.gs == 1){
                        Triangle3 tri(vertices[f.v[0]].p, vertices[f.v[1]].p, vertices[f.v[2]].p);
                        double area = sqrt(tri.squared_area());
                        AB += area;
                    }
                }
                perB += perimeter_building(b->second.gs2D);
            }
        }

        //Box 1
        //compute area refinement box 1 minus the area of the cells around the buildings and building areas
        //the latter is done by using the building perimeters and cell areas with refinement level 3
        double n1_geometry = (perB/(meshCFD.wcell/8))*4;
        N_refined += n1_geometry;
        double A_geometry = n1_geometry*A_cell3;
        double A_box1 = (box1.xmax()-box1.xmin())*(box1.ymax()-box1.ymin());
        double A_ground1 = A_box1-(AB+A_geometry);

        //compute number of cells n_ground1 from ground refinement in refinement box 1
        double n1_ground1 = (A_ground1/A_cell3)*height;
        n1_ground1 += N_refined;
        double n2_ground1 = (A_ground1/A_cell2)*4;
        double n_ground1 = n1_ground1 + n2_ground1;

        //compute the volume of cells taken by n_ground1 and subtract this from volume refinement box 1
        //to be able to compute number of cells within this box
        double V_ground1 = n1_ground1*V_cell3 + n2_ground1*V_cell2;
        vol1 = vol1 - V_ground1;
        double n1 = vol1/V_cell2;

        //Box 2
        //follow the same steps as with box 1
        //except that an extra layer of cells is needed for a smooth transition between refinement level 3 to 1
        double A_box2 = ((box2.xmax()-box2.xmin())*(box2.ymax()-box2.ymin()));
        double A_ground2 = A_box2 -A_box1;
        double n1_ground2 = (A_ground2/A_cell3)*height;
        N_refined += n1_ground2;
        double n2_ground2 = (A_ground2/A_cell2)*4;
        double n3_ground2 = (A_ground2/A_cell1)*4;
        double n_ground2 = n1_ground2 + n2_ground2 + n3_ground2;

        double V_ground2 = n1_ground2*V_cell3 + n2_ground2*V_cell2 + n3_ground2*V_cell1;
        vol2 = vol2 - V_ground2;
        double n2 = vol2/V_cell1;

        //Domain
        //follow the same steps as with box 1
        //except that two extra layers of cells are needed for a smooth transition between refinement level 3 to 0
        double A_domain = (domain.xmax()-domain.xmin())*(domain.ymax()-domain.ymin());
        double A_groundD = A_domain-A_box2;
        double n1_groundD = (A_groundD/A_cell3)*height;
        N_refined += n1_groundD;
        double n2_groundD = (A_groundD/A_cell2)*4;
        double n3_groundD = (A_groundD/A_cell1)*4;
        double n4_groundD = (A_groundD/A_cell0)*4;
        double n_groundD = n1_groundD + n2_groundD + n3_groundD + n4_groundD;

        double V_groundD = n1_groundD*V_cell3 + n2_groundD*V_cell2 + n3_groundD*V_cell1 + n4_groundD*V_cell0;
        volD = volD - V_groundD;
        double nD = volD/V_cell0;

        //sum up all the computed number of cells
        N = nG + n1 + n2 + nD + n_ground1 + n_ground2 + n_groundD;
    }

    meshCFD.N_refined = N_refined;

    return N;
}


double B3_n_cells_snappy(int gr, double height) {

    //same algorithm as the previous one, but then with 3 refinement boxes

    meshCFD.N_refined = 0;

    CGAL::Bbox_3 box1 = meshCFD.box1;
    CGAL::Bbox_3 box2 = meshCFD.box2;
    CGAL::Bbox_3 box3 = meshCFD.box3;
    CGAL::Bbox_3 domain = meshCFD.domain;

    double N_refined = 0;

    double volB = 0;
    double volG = 0;
    double nG = 0;
    for (map<int, Building>::iterator b = buildings.begin(); b != buildings.end(); b++) {
        SMesh building = b->second.mesh;
        double vol = CGAL::Polygon_mesh_processing::volume(building);
        volB += vol;

        for (auto &f: b->second.faces){
            if (f.gs == 0){
                Triangle3 tri(vertices[f.v[0]].p, vertices[f.v[1]].p, vertices[f.v[2]].p);
                double area = sqrt(tri.squared_area());
                double n = (area/(meshCFD.wcell/16*meshCFD.hcell/16))*4;
                N_refined += n;
                double n_extra = (area/(meshCFD.wcell/8*meshCFD.hcell/8))*4;
                double v = n*pow(meshCFD.wcell/16, 2)*meshCFD.hcell/16;
                double v_extra = n_extra*pow(meshCFD.wcell/8, 2)*meshCFD.hcell/8;
                volG += v ;
                volG += v_extra;
                nG += n;
                nG += n_extra;
            }
        }
    }

    double A_cell4 = pow(meshCFD.wcell/16, 2);
    double A_cell3 = pow(meshCFD.wcell/8, 2);
    double A_cell2 = pow(meshCFD.wcell/4, 2);
    double A_cell1 = pow(meshCFD.wcell/2, 2);
    double A_cell0 = pow(meshCFD.wcell, 2);

    double V_cell4 = A_cell4*(meshCFD.hcell/16);
    double V_cell3 = A_cell3*(meshCFD.hcell/8);
    double V_cell2 = A_cell2*(meshCFD.hcell/4);
    double V_cell1 = A_cell1*(meshCFD.hcell/2);
    double V_cell0 = A_cell0*meshCFD.hcell;

    double vol1 = ((box1.xmax()-box1.xmin())*(box1.ymax()-box1.ymin())*(box1.zmax()-box1.zmin()))-(volB+volG);
    double vol2 = ((box2.xmax()-box2.xmin())*(box2.ymax()-box2.ymin())*(box2.zmax()-box2.zmin()))-(vol1+volB+volG);
    double vol3 = ((box3.xmax()-box3.xmin())*(box3.ymax()-box3.ymin())*(box3.zmax()-box3.zmin()))-(vol1+vol2+volB+volG);
    double volD = ((domain.xmax()-domain.xmin())*(domain.ymax()-domain.ymin())*(domain.zmax()-domain.zmin()))-(vol1+volB+volG+vol2+vol3);

    double N;

    if (gr == 0){

        double n1;
        double n2;
        double n3;
        double nDomain;

        if (fTerrain.size() == 0){

            n1 = vol1/(pow(meshCFD.wcell/8, 2)*meshCFD.hcell/8);
            n2 = vol2/(pow(meshCFD.wcell/4, 2)*meshCFD.hcell/4);
            n3 = vol3/(pow(meshCFD.wcell/2, 2)*meshCFD.hcell/2);
            nDomain = volD/(pow(meshCFD.wcell, 2)*meshCFD.hcell);

            N = nG + n1 + n2 + n3 + nDomain;

        } else {

            map<int, pair<double, double>> nv = B3_ground_volume();

            double v_box1 = meshCFD.v_box1;
            double v_box2 = meshCFD.v_box2;
            double v_box3 = meshCFD.v_box3;
            double v_domain = meshCFD.v_domain;

            vol1 = vol1 - (v_box1+nv[1].second);
            vol2 = vol2 - (v_box2+nv[2].second);
            vol3 = vol3 - (v_box3+nv[3].second);
            volD = volD - (v_domain+nv[4].second);

            n1 = vol1/(pow(meshCFD.wcell/8, 2)*meshCFD.hcell/8);
            n2 = vol2/(pow(meshCFD.wcell/4, 2)*meshCFD.hcell/4);
            n3 = vol3/(pow(meshCFD.wcell/2, 2)*meshCFD.hcell/2);
            nDomain = volD/(pow(meshCFD.wcell, 2)*meshCFD.hcell);

            double nT = nv[1].first+nv[2].first+nv[3].first+nv[4].first;

            N = nG + n1 + n2 + n3 + nDomain + nT;

        }

    } else {
        double AB = 0;
        double perB = 0;
        for (map<int, Building>::iterator b = buildings.begin(); b != buildings.end(); b++) {
            if (!b->second.gs2D.empty()){
                for (auto &f: b->second.faces){
                    if (f.gs == 1){
                        Triangle3 tri(vertices[f.v[0]].p, vertices[f.v[1]].p, vertices[f.v[2]].p);
                        double area = sqrt(tri.squared_area());
                        AB += area;
                    }
                }
                perB += perimeter_building(b->second.gs2D);
            }
        }

        double n1_geometry = (perB/(meshCFD.wcell/16))*4;
        N_refined += n1_geometry;
        double n_geometry = n1_geometry;

        double A_geometry = n1_geometry*A_cell4;
        double A_box1 = (box1.xmax()-box1.xmin())*(box1.ymax()-box1.ymin());
        double A_ground1 = A_box1-(AB+A_geometry);

        double n1_ground1 = (A_ground1/A_cell4)*height;
        n1_ground1 += N_refined;
        double n2_ground1 = (A_ground1/A_cell3)*4;
        double n_ground1 = n1_ground1 + n2_ground1;

        double V_ground1 = n1_ground1*V_cell4 + n2_ground1*V_cell3;
        vol1 = vol1 - V_ground1;
        double n1 = vol1/V_cell3;

        //Box 2
        double A_box2 = ((box2.xmax()-box2.xmin())*(box2.ymax()-box2.ymin()));
        double A_ground2 = A_box2 -A_box1;
        double n1_ground2 = (A_ground2/A_cell4)*height;
        N_refined += n1_ground2;
        double n2_ground2 = (A_ground2/A_cell3)*4;
        double n3_ground2 = (A_ground2/A_cell2)*4;
        double n_ground2 = n1_ground2 + n2_ground2 + n3_ground2;

        double V_ground2 = n1_ground2*V_cell4 + n2_ground2*V_cell3 + n3_ground2*V_cell2;
        vol2 = vol2 - V_ground2;
        double n2 = vol2/V_cell2;

        //Box 3
        double A_box3 = ((box3.xmax()-box3.xmin())*(box3.ymax()-box3.ymin()));
        double A_ground3 = A_box3 -A_box2;
        double n1_ground3 = (A_ground3/A_cell4)*height;
        N_refined += n1_ground3;
        double n2_ground3 = (A_ground3/A_cell3)*4;
        double n3_ground3 = (A_ground3/A_cell2)*4;
        double n4_ground3 = (A_ground3/A_cell1)*4;
        double n_ground3 = n1_ground3 + n2_ground3 + n3_ground3 + n4_ground3;

        double V_ground3 = n1_ground3*V_cell4 + n2_ground3*V_cell3 + n3_ground3*V_cell2 + n4_ground3*V_cell1;
        vol3 = vol3 - V_ground3;
        double n3 = vol3/V_cell1;

        //Domain
        double A_domain = (domain.xmax()-domain.xmin())*(domain.ymax()-domain.ymin());
        double A_groundD = A_domain-A_box3;
        double n1_groundD = (A_groundD/A_cell4)*height;
        N_refined += n1_groundD;
        double n2_groundD = (A_groundD/A_cell3)*4;
        double n3_groundD = (A_groundD/A_cell2)*4;
        double n4_groundD = (A_groundD/A_cell1)*4;
        double n5_groundD = (A_groundD/A_cell0)*4;
        double n_groundD = n1_groundD + n2_groundD + n3_groundD + n4_groundD + n5_groundD;

        double V_groundD = n1_groundD*V_cell4 + n2_groundD*V_cell3 + n3_groundD*V_cell2 + n4_groundD*V_cell1 + n5_groundD*V_cell0;
        volD = volD - V_groundD;
        double nD = volD/V_cell0;

        N = nG + n1 + n2 + nD + n_ground1 + n_ground2 + n_ground3 + n_groundD;
    }

    meshCFD.N_refined = N_refined;

    return N;
}






