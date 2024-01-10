#ifndef BACKEND_DEFINITIONS_H
#define BACKEND_DEFINITIONS_H

using namespace std;

//CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Plane_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
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

//PMP
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

namespace PMP = CGAL::Polygon_mesh_processing;
typedef CGAL::Surface_mesh<Point3>                              SMesh;
typedef SMesh::Vertex_index                                     vertex_descriptor;
typedef boost::graph_traits<SMesh>::vertex_descriptor           vertex_descriptor;

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
typedef CGAL::Delaunay_triangulation_2<K,Tds>                Triangulation2;
typedef CGAL::Alpha_shape_2<Triangulation2>                  Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator            Alpha_shape_edges_iterator;

//AABB Tree
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
typedef std::list<Triangle3>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
typedef Tree::Primitive_id Primitive_id;

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

#include <boost/format.hpp>
namespace params = CGAL::parameters;

typedef Triangulation::Point          Point;

#include <CGAL/intersections.h>


#endif
