#include "definitions.h"
#include "input.h"

Citymodel city;

map<int, Face> fiBuildings;
map<int, Vertex> vi;

map<int, Building> buildings;
map<int, Face> fBuildings;
map<int, Face> fTerrain;
map<int, Vertex> vertices;

MeshCFD meshCFD;


void OBJ_data(string input_file, double wd) {

    ifstream input_stream;
    input_stream.open(input_file);

    int id_face = 0, id_vertex = 0;

    if (input_stream.is_open()) {
        string line, data;

        vector<Point3> pts;
        int objectType;

        while (getline(input_stream, line)) {
            istringstream tmp_stream(line);
            tmp_stream >> data;

            if (data == "v"){
                id_vertex += 1;
                double x, y, z;
                tmp_stream >> x >> y >> z;
                double wind = -wd*(M_PI/180);
                double xwind = x*cos(wind)-y*sin(wind);
                double ywind = x*sin(wind)+y*cos(wind);
                Point3 p(xwind, ywind, z);
                vi[id_vertex].id = id_vertex;
                vi[id_vertex].p = p;
            }

            else if (data == "g"){
                string type;
                tmp_stream >> type;

                //objectType 0:terrain, 1:building, 2:others
                if (type == "Terrain" || type == "Water") {
                    objectType = 0;
                }

                if (type == "Buildings"){
                    objectType = 1;
                }

                else if (type != "Terrain" && type != "Buildings" && type != "Water" ){
                    objectType = 2;
                }
            }


            else if (data == "o"){
                objectType = 1;
            }

            else if (data == "f"){
                id_face += 1;
                int v1, v2, v3;
                tmp_stream >> v1 >> v2 >> v3;

                if (objectType == 0){ //type 0:terrain
                    fTerrain[id_face].v.push_back(v1);
                    fTerrain[id_face].v.push_back(v2);
                    fTerrain[id_face].v.push_back(v3);
                }

                if (objectType == 1){ //type 1:building
                    fiBuildings[id_face].v.push_back(v1);
                    fiBuildings[id_face].v.push_back(v2);
                    fiBuildings[id_face].v.push_back(v3);
                    city.vb.push_back(vi[v1].p);
                    city.vb.push_back(vi[v2].p);
                    city.vb.push_back(vi[v3].p);
                }
            }
        }
    }
};


void STL_data(string input_file, double wd) {
    ifstream input_stream;
    input_stream.open(input_file);

    if (input_stream.is_open()){
        string line;
        string data;

        int id_face = 0;
        int id_vertex = 0;
        while (getline(input_stream, line)){
            std::istringstream tmp_stream(line);
            tmp_stream >> data;

            if (data == "facet"){
                id_face += 1;
                fiBuildings[id_face].id = id_face;
            }

            else if (data == "vertex"){
                double x, y, z;
                tmp_stream >> x >> y >> z;
                double wind = -wd*(M_PI/180);
                double xwind = x*cos(wind)-y*sin(wind);
                double ywind = x*sin(wind)+y*cos(wind);
                Point3 p(xwind, ywind, z);
                id_vertex += 1;
                vi[id_vertex].id = id_vertex;
                vi[id_vertex].p = p;
                fiBuildings[id_face].v.push_back(id_vertex);
                city.vb.push_back(p);
            }
        }
    }
};