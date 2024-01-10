#include "definitions.h"
#include "input.h"

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
