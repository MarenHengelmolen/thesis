#include "definitions.h"
#include "input.h"

Citymodel city; //stores data of the urban area

map<int, Face> fiBuildings; //stores building faces (before building definition)
map<int, Vertex> vi; //stores vertices (before building definition)

map<int, Building> buildings; //stores buildings
map<int, Face> fBuildings; //stores building faces (after building definition)
map<int, Face> fTerrain; //stores terrain faces
map<int, Vertex> vertices; //stores vertices (after building definition)

MeshCFD meshCFD; //stores mesh parameters
Output output; //stores output data and paths to return output files

void OBJ_data(string input_file, double fd, double unit) {

    //reads and stores data from OBJ input file

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

            if (data == "v"){ //stores rotated vertices based on flow direction fd
                id_vertex += 1;
                double x, y, z;
                tmp_stream >> x >> y >> z;
                x = x*unit;
                y = y*unit;
                z = z*unit;

                double wind = -(90-fd)*(M_PI/180);
                double xwind = x*cos(wind)-y*sin(wind);
                double ywind = x*sin(wind)+y*cos(wind);
                Point3 p(xwind, ywind, z);
                vi[id_vertex].id = id_vertex;
                vi[id_vertex].p = p;
            }

            else if (data == "g"){ //stores faces
                string type;
                tmp_stream >> type;

                // if face is stored under:
                //- terrain, objectType -> 0
                //- building, objectType -> 1
                //- others, objectType -> 2
                if (type != "Terrain" || type == "Water" || type == "Vegetation") {
                    objectType = 0;
                }

                if (type == "Buildings"){
                    objectType = 1;
                }

                else if (type != "Terrain" && type != "Buildings" && type != "Water" && type != "Vegetation"){
                    objectType = 2;
                }
            }

            else if (data == "o"){
                objectType = 1;
            }

            else if (data == "f"){ //stores vertices per face in fTerrain or fiBuildings, dependent if they are terrain or building features
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
                    city.vb.push_back(vi[v1].p); //stores vertices in city
                    city.vb.push_back(vi[v2].p);
                    city.vb.push_back(vi[v3].p);
                }
            }
        }
    }

};


void STL_data(string input_file, double fd, double unit) {
    ifstream input_stream;
    input_stream.open(input_file);

    //reads and stores data from STL input file

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
                x = x*unit;
                y = y*unit;
                z = z*unit;
                double wind = -(90-fd)*(M_PI/180);
                double xwind = x*cos(wind)-y*sin(wind);
                double ywind = x*sin(wind)+y*cos(wind);
                Point3 p(xwind, ywind, z);
                id_vertex += 1;
                vi[id_vertex].id = id_vertex; //stores rotated vertices based on flow direction fd
                vi[id_vertex].p = p;
                fiBuildings[id_face].v.push_back(id_vertex); //stores vertices per face
                city.vb.push_back(p); //stores vertices in city
            }
        }
    }
};