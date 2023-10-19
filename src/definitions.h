#ifndef THESIS_V3_DEFINITIONS_H
#define THESIS_V3_DEFINITIONS_H

using namespace std;

#include <set>

//CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

#include <CGAL/Triangle_2.h>
#include <CGAL/Plane_3.h>

typedef K::Point_3                  Point3;
typedef K::Vector_3                 Vector3;
typedef K::Point_2                  Point2;
typedef K::Triangle_2               Triangle2;
typedef K::Circle_2                 Circle2;
typedef K::Sphere_3                 Sphere3;
typedef K::Triangle_3               Triangle3;
typedef K::Plane_3                  Plane3;


//CONVEX HULL
#include <CGAL/convex_hull_2.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
typedef CGAL::Convex_hull_traits_adapter_2<K, CGAL::Pointer_property_map<Point2>::type > Convex_hull_traits_2;

#include <CGAL/Plane_3.h> //ground surface to_2d

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>

//Alpha shape
typedef K::FT                                                FT;
typedef K::Segment_2                                         Segment2;

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/algorithm.h>

typedef CGAL::Alpha_shape_vertex_base_2<K>                   Vb;
typedef CGAL::Alpha_shape_face_base_2<K>                     Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>          Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds>                Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2>                 Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator            Alpha_shape_edges_iterator;

namespace PMP = CGAL::Polygon_mesh_processing;
typedef CGAL::Surface_mesh<Point3>                      SMesh;
typedef SMesh::Vertex_index                             vertex_descriptor;
typedef SMesh::Face_index                               face_descriptor;
typedef SMesh::Vertex_iterator                          Vertex_iterator;
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
typedef boost::graph_traits<SMesh>::vertex_descriptor           vertex_descriptor;
typedef K::Vector_3                                             Vector;
typedef boost::property_map<SMesh, CGAL::vertex_point_t>::type  VPMap;
typedef SMesh::template Property_map<vertex_descriptor, Vector> VNMap;
typedef CGAL::Triple<int, int, int> Triangle_int;




#include <CGAL/compute_average_spacing.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_points.h>

typedef std::pair<Point3, Vector3> PointVectorPair;
typedef CGAL::Parallel_if_available_tag Concurrency_tag;
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/poisson_surface_reconstruction.h>
typedef std::pair<Point3, Vector3> Pwn;
typedef CGAL::Polyhedron_3<K> Polyhedron;
#include <CGAL/IO/read_points.h>

#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

typedef std::list<Triangle3>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
typedef Tree::Primitive_id Primitive_id;


#include <CGAL/point_generators_3.h>
typedef CGAL::Random_points_in_triangle_3<Point3> Point_generator;

#include <CGAL/point_generators_2.h>

//#include <CGAL/Triangulation_2.h>
//typedef CGAL::Triangulation_2<K>                    Triangulation;
//typedef Triangulation::Vertex_handle                Vertex_handle;
//typedef Triangulation::Point                        Point;
//typedef Triangulation::Finite_vertex_handles        Finite_vertex_handles;
//typedef Triangulation::Finite_vertices_iterator     Finite_vertices_iterator;
//typedef Finite_vertex_handles::iterator             Finite_vertex_handles_iterator;
//typedef Triangulation::Face_handle Face_handle;

//typedef Triangulation::Face_iterator Face_iterator;


#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
typedef CGAL::Triangulation_vertex_base_2<K> VertexBase;
typedef CGAL::Constrained_triangulation_face_base_2<K> FaceBase;
struct FaceInfo {
    int nesting_level;
    bool in_domain(){
        return nesting_level%2 == 1;
    }
};

typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, K, FaceBase> FaceBaseWithInfo;
typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBaseWithInfo> TriangulationDataStructure;
typedef CGAL::Exact_predicates_tag Tag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TriangulationDataStructure, Tag> Triangulation;

typedef Triangulation::Face_handle Face_handle;




#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <boost/format.hpp>
namespace params = CGAL::parameters;
typedef boost::graph_traits<SMesh>::edge_descriptor            edge_descriptor;

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <CGAL/Polygon_2.h>
typedef CGAL::Polygon_2<K> Polygon2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator               Alpha_shape_edges_iterator;

#include <CGAL/optimal_bounding_box.h>



#endif
