#pragma once
// Internal helpers shared by the functional API (not installed).
//
// The legacy modules assume that the input lives in [0, scale_max]^3
// (scale_max = 1000) and that the surface is stored in a SurfaceMesh with
// AABB trees and (optionally) detected sharp/concave features. The helpers
// here build that state from plain Eigen matrices, following the same order
// as the MATFP/MATTopo pipelines.

#include <Eigen/Core>
#include <memory>

#include "libmat/inputs/input_types.h"
#include "libmat/inputs/params.h"

namespace libmat {
namespace internal {

// Must be called before any geogram call. Safe to call many times.
void ensure_geogram_initialized();

// Affine map from user coordinates to the normalized [0, scale_max]^3 box:
// x_normalized = (x - offset) * factor
struct Normalization {
  Eigen::RowVector3d offset = Eigen::RowVector3d::Zero();
  double factor = 1.;

  Eigen::RowVector3d to_normalized(const Eigen::RowVector3d& p) const {
    return (p - offset) * factor;
  }
  Eigen::RowVector3d from_normalized(const Eigen::RowVector3d& p) const {
    return p / factor + offset;
  }
};

// Computes the normalization of V and fills the bbox fields of params
// (bbox before/after normalization and the 8 bbox corners).
Normalization normalize_and_fill_params(const Eigen::MatrixXd& V,
                                        Parameter& params,
                                        Eigen::MatrixXd& V_normalized);

// Fills the bbox fields of params from already normalized vertices.
void fill_params_bbox(const Eigen::MatrixXd& V_normalized, Parameter& params);

// Builds sf_mesh from (V, F) the same way load_sf_mesh_from_internal() does
// (vertex attribute "tet_vid" = sf2tet[v]), and additionally stores the input
// face index in the facet attribute "orig_fid". The mesh is Morton-reordered.
void build_sf_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                   const std::vector<int>& sf2tet, SurfaceMesh& sf_mesh);

// Surface-only context: the "tet mesh" holds the surface vertices only (no
// tets) and SurfaceMesh::sf2tet_vs_mapping is the identity, so that feature
// edges are reported in input vertex indices.
struct SurfaceContext {
  Parameter params;
  Normalization normalization;
  std::unique_ptr<TetMesh> tet_mesh;
  std::unique_ptr<SurfaceMesh> sf_mesh;
};

// V, F: input triangle mesh (closed, manifold) in user coordinates.
void build_surface_context(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                           bool detect_features, double thres_concave,
                           double thres_convex, SurfaceContext& ctx);

// Validation helpers, throw std::invalid_argument
void check_cols(const Eigen::MatrixXd& M, int cols, const char* name);
void check_cols(const Eigen::MatrixXi& M, int cols, const char* name);
void check_indices(const Eigen::MatrixXi& M, int num_vertices,
                   const char* name);

}  // namespace internal
}  // namespace libmat
