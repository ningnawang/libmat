#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/point_generators_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel CGAL_K;
typedef CGAL_K::Point_3 CGAL_Point;
typedef CGAL::Polyhedron_3<CGAL_K> CGAL_Mesh;