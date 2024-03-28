#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <geogram/delaunay/periodic_delaunay_3d.h>

#include "medial_sphere.h"
#include "params.h"

//////////
// CGAL

// Regular Vertex info (dual to sphere)
class RVI {
 public:
  // matching MedialSphere::id
  // remains -1 if 8 bbox points
  int all_id = -1;
};

// Regular Tetrahedron info
class RTI {
 public:
  int id = -1;  // assigned while iterating RT, not unique for each run
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kf;
typedef Kf::FT Weight;
typedef Kf::Point_3 Point;
typedef Kf::Weighted_point_3 Weighted_point;

typedef CGAL::Regular_triangulation_vertex_base_3<Kf> Vb0_rt;
typedef CGAL::Triangulation_vertex_base_with_info_3<RVI, Kf, Vb0_rt> Vb_rt;
typedef CGAL::Regular_triangulation_cell_base_3<Kf> Cb0_rt;
typedef CGAL::Triangulation_cell_base_with_info_3<RTI, Kf, Cb0_rt> Cb_rt;
typedef CGAL::Triangulation_data_structure_3<Vb_rt, Cb_rt> Tds_rt;
typedef CGAL::Regular_triangulation_3<Kf, Tds_rt> Rt;
typedef CGAL::Triangulation_3<Kf, Tds_rt> Tr;

typedef Kf::Point_3 Point_rt;
typedef Rt::Vertex_iterator Vertex_iterator_rt;
typedef Rt::Vertex_handle Vertex_handle_rt;
typedef Rt::Cell_iterator Cell_iterator_rt;
typedef Rt::Cell_handle Cell_handle_rt;
typedef Rt::Cell_circulator Cell_circulator_rt;
typedef Rt::Facet_circulator Facet_circulator_rt;

typedef Rt::Finite_cells_iterator Finite_cells_iterator_rt;
typedef Rt::Finite_facets_iterator Finite_facets_iterator_rt;
typedef Rt::Finite_edges_iterator Finite_edges_iterator_rt;
typedef Rt::Finite_vertices_iterator Finite_vertices_iterator_rt;
typedef Rt::Tetrahedron Tetrahedron_rt;

///////////////
class RegularTriangulationNN : public Rt, public GEO::Counted {
 public:
  ~RegularTriangulationNN(){};

 public:
  inline void clean() {
    this->clear();
    nb_vertices = 0;
  }

 protected:
  // no use for now
  // nb_vertices >= number_of_vertices()
  // nb_vertices = all_medial_spheres.size()
  int nb_vertices;
};

void generate_RT_CGAL_given_spheres(const Parameter& params,
                                    const std::vector<float>& spheres,
                                    RegularTriangulationNN& rt, bool is_debug);

void generate_RT_CGAL_and_mark_valid_spheres(
    const Parameter& params, std::vector<MedialSphere>& all_medial_spheres,
    RegularTriangulationNN& rt, std::set<int>& valid_sphere_ids);

void generate_RT_CGAL_and_purge_spheres(
    const Parameter& params, std::vector<MedialSphere>& all_medial_spheres,
    RegularTriangulationNN& rt);

// -----------------------------------------------------------------------------
int get_RT_spheres_and_neighbors(const int num_spheres,
                                 const RegularTriangulationNN& rt,
                                 std::vector<int>& site_knn,
                                 std::vector<int>& is_sphere_valid,
                                 bool is_debug);

int get_RT_partial_spheres_and_neighbors(
    const std::set<int>& sphere_ids, const RegularTriangulationNN& rt,
    std::vector<MedialSphere> all_medial_spheres,
    std::vector<int>& map_site2msphere, std::set<int>& spheres_and_1rings,
    std::vector<int>& site_knn, bool is_debug = false);

int get_RT_vertex_neighbors(const RegularTriangulationNN& rt, const int& n_site,
                            std::vector<int>& site_knn);

// update MedialSphere::rt_neigh_ids_prev
int update_spheres_RT_neighbors(const RegularTriangulationNN& rt,
                                std::vector<MedialSphere>& all_medial_spheres);

void get_PC_vertices(const RegularTriangulationNN& rt,
                     std::vector<Vector3>& pc_vertices);