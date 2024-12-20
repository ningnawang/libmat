#ifndef H_UPDATING_H
#define H_UPDATING_H

#include <assert.h>

#include <cmath>
#include <cstdlib>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "common_geogram.h"
#include "medial_mesh.h"
#include "medial_sphere.h"
#include "shrinking.h"
// Default parameters:
// alpha1 = 0.01;  // energy of distance to tangent point
// alpha2 = 1;     // energy of distance to tangent plane
// alpha3 = 1;     // energy of distance to concave line

struct iter_params {
  bool is_check_new_tan_plane;
  double alpha1;
  double alpha2;
  double alpha3;
  double break_threshold;
  int itr_limit;
  bool is_break_use_energy_over_sq_radius;

  iter_params() {
    is_check_new_tan_plane = true;
    alpha1 = 0.01;
    alpha2 = 1.;
    alpha3 = 1.;
    break_threshold = SCALAR_ZERO_3;
    itr_limit = 30;
    is_break_use_energy_over_sq_radius = true;
  };
};

// Two steps:
// 1. update medial sphere center and radius
// 2. update tangent planes and concave lines
bool iterate_sphere(
    const SurfaceMesh& sf_mesh, const AABBWrapper& aabb_wrapper,
    const std::set<aint2>& fe_sf_fs_pairs,
    const std::vector<FeatureEdge>& feature_edges, MedialSphere& mat_p,
    bool is_debug,
    // const bool is_check_new_tan_plane = true,
    // double alpha1 = 0.01, double alpha2 = 1, double alpha3 = 1.,
    // const double break_threshold = SCALAR_ZERO_3,
    // const int itr_limit = 30
    iter_params params = iter_params());

// Two steps (reversed):
// 1. update medial sphere center and radius
// 2. update tangent planes and concave lines
bool iterate_sphere_reversed(
    const SurfaceMesh& sf_mesh, const AABBWrapper& aabb_wrapper,
    const std::set<aint2>& fe_sf_fs_pairs,
    const std::vector<FeatureEdge>& feature_edges, MedialSphere& mat_p,
    bool is_debug,
    // bool is_check_new_tan_plane = true, double alpha1 = 0.01,
    // double alpha2 = 1, double alpha3 = 1.,
    // const double break_threshold = SCALAR_ZERO_3,
    // const bool is_break_use_energy_over_sq_radius = true,
    // const int itr_limit = 30
    iter_params params = iter_params());

// void init_shrink_and_update(const SurfaceMesh& sf_mesh, const TetMesh&
// tet_mesh,
//                             std::vector<MedialSphere>& all_medial_spheres,
//                             int num_init_spheres, bool is_debug = false);

// relax spheres to centroids of powercell then iterate spheres
// type_to_handle:
// 0 - all spheres except for external features
// 1 - SphereType::T_2 only
// 2 - internal features only
void relax_and_iterate_spheres(const SurfaceMesh& sf_mesh,
                               const std::vector<FeatureEdge>& feature_edges,
                               std::vector<MedialSphere>& all_medial_spheres,
                               const bool site_is_transposed,
                               std::vector<float>& site_updated,
                               int type_to_handle, bool is_debug);
// Optimal Delaunay Triangulation
void relax_and_iterate_spheres_ODT(
    const SurfaceMesh& sf_mesh, const std::vector<FeatureEdge>& feature_edges,
    const MedialMesh& mmesh, std::vector<MedialSphere>& all_medial_spheres,
    int type_to_handle, bool is_debug);

// relax spheres to centroids of filtered neighbors (Laplacian smoothing)
void relax_and_iterate_spheres_Laplacian(
    const SurfaceMesh& sf_mesh, const std::vector<FeatureEdge>& feature_edges,
    const MedialMesh& mmesh, std::vector<MedialSphere>& all_medial_spheres,
    bool is_debug);

void relax_and_iterate_spheres_both(
    const SurfaceMesh& sf_mesh, const std::vector<FeatureEdge>& feature_edges,
    const MedialMesh& mmesh, std::vector<MedialSphere>& all_medial_spheres,
    const bool site_is_transposed, std::vector<float>& site_updated,
    bool is_debug);
#endif  // __UPDATING_H__