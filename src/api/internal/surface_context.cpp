#include "api/internal/surface_context.h"

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/common.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh_reorder.h>

#include <bitset>
#include <mutex>
#include <stdexcept>
#include <string>

#include "libmat/inputs/sharp_feature_detection.h"
#include "libmat/matfun/shrinking.h"

namespace libmat {
namespace internal {

void ensure_geogram_initialized() {
  static std::once_flag flag;
  std::call_once(flag, []() {
    // no signal handlers: they would hijack the python interpreter
    GEO::initialize(GEO::GEOGRAM_NO_HANDLER);
    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("algo");
    GEO::Logger::instance()->set_quiet(true);
  });
}

void fill_params_bbox(const Eigen::MatrixXd& VN, Parameter& params) {
  const Eigen::RowVector3d pmin = VN.colwise().minCoeff();
  const Eigen::RowVector3d pmax = VN.colwise().maxCoeff();
  params.xmin = pmin[0], params.ymin = pmin[1], params.zmin = pmin[2];
  params.xmax = pmax[0], params.ymax = pmax[1], params.zmax = pmax[2];
  params.bbox_diag_l = (pmax - pmin).norm();
  // 8 bbox points, same order as load_tet()
  params.bb_points.clear();
  for (int i = 0; i < 8; i++) {
    std::bitset<sizeof(int) * 8> flag(i);
    for (int j = 0; j < 3; j++)
      params.bb_points.push_back(flag.test(j) ? pmax[j] : pmin[j]);
  }
}

Normalization normalize_and_fill_params(const Eigen::MatrixXd& V,
                                        Parameter& params,
                                        Eigen::MatrixXd& VN) {
  if (V.rows() == 0) throw std::invalid_argument("V is empty");
  const Eigen::RowVector3d pmin = V.colwise().minCoeff();
  const Eigen::RowVector3d pmax = V.colwise().maxCoeff();
  const double maxside = (pmax - pmin).maxCoeff();
  if (!(maxside > 0)) throw std::invalid_argument("V has a degenerate bbox");

  Normalization n;
  n.offset = pmin;
  n.factor = params.scale_max / maxside;
  VN = (V.rowwise() - n.offset) * n.factor;

  params.xmin_orig = pmin[0], params.ymin_orig = pmin[1];
  params.zmin_orig = pmin[2];
  params.xmax_orig = pmax[0], params.ymax_orig = pmax[1];
  params.zmax_orig = pmax[2];
  params.bbox_diag_l_prev = (pmax - pmin).norm();
  params.scale_maxside_orig = maxside;
  fill_params_bbox(VN, params);
  return n;
}

void build_sf_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                   const std::vector<int>& sf2tet, SurfaceMesh& sf_mesh) {
  sf_mesh.clear(false, false);
  GEO::Attribute<int> tet_vid_attr(sf_mesh.vertices.attributes(), "tet_vid");
  GEO::Attribute<int> orig_fid_attr(sf_mesh.facets.attributes(), "orig_fid");

  sf_mesh.vertices.create_vertices(V.rows());
  for (int i = 0; i < V.rows(); ++i) {
    GEO::vec3& p = sf_mesh.vertices.point(i);
    p[0] = V(i, 0), p[1] = V(i, 1), p[2] = V(i, 2);
    tet_vid_attr[i] = sf2tet.empty() ? i : sf2tet[i];
  }
  sf_mesh.facets.create_triangles(F.rows());
  for (int f = 0; f < F.rows(); ++f) {
    for (int lv = 0; lv < 3; ++lv) sf_mesh.facets.set_vertex(f, lv, F(f, lv));
    orig_fid_attr[f] = f;
  }

  std::vector<aint2> edges;
  edges.reserve(3 * F.rows());
  for (int f = 0; f < F.rows(); ++f) {
    for (int j = 0; j < 3; j++) {
      aint2 e = {{F(f, j), F(f, (j + 1) % 3)}};
      if (e[0] > e[1]) std::swap(e[0], e[1]);
      edges.push_back(e);
    }
  }
  vector_unique(edges);
  sf_mesh.edges.create_edges(edges.size());
  for (uint e = 0; e < edges.size(); e++) {
    sf_mesh.edges.set_vertex(e, 0, edges[e][0]);
    sf_mesh.edges.set_vertex(e, 1, edges[e][1]);
  }

  GEO::mesh_reorder(sf_mesh, GEO::MESH_ORDER_MORTON);
  sf_mesh.facets.connect();
}

void build_surface_context(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                           bool detect_features, double thres_concave,
                           double thres_convex, SurfaceContext& ctx) {
  ensure_geogram_initialized();
  check_cols(V, 3, "V");
  check_cols(F, 3, "F");
  check_indices(F, V.rows(), "F");
  if (F.rows() == 0) throw std::invalid_argument("F is empty");

  ctx.params = Parameter();
  ctx.params.thres_concave = thres_concave;
  ctx.params.thres_convex = thres_convex;
  Eigen::MatrixXd VN;
  ctx.normalization = normalize_and_fill_params(V, ctx.params, VN);

  // surface-only "tet mesh": vertices only, identity sf2tet mapping
  ctx.tet_mesh = std::make_unique<TetMesh>("");
  ctx.tet_mesh->tet_vertices.resize(3 * VN.rows());
  for (int i = 0; i < VN.rows(); ++i)
    for (int j = 0; j < 3; ++j) ctx.tet_mesh->tet_vertices[3 * i + j] = VN(i, j);

  ctx.sf_mesh = std::make_unique<SurfaceMesh>();
  SurfaceMesh& sf_mesh = *ctx.sf_mesh;
  build_sf_mesh(VN, F, {}, sf_mesh);
  // the AABB tree reorders sf_mesh again, so it must be built before any
  // facet index (e.g. feature edges) is cached
  pre_and_init_aabb(sf_mesh, sf_mesh.aabb_wrapper);
  sf_mesh.reload_sf2tet_vs_mapping();
  std::map<int, int> tet2sf_vs;
  std::map<int, std::set<int>> sf_vs2fids;
  load_sf_tet_mapping(sf_mesh, tet2sf_vs, sf_vs2fids,
                      ctx.tet_mesh->tet_vs2sf_fids);

  if (detect_features) {
    detect_mark_sharp_features(ctx.params, sf_mesh, *ctx.tet_mesh);
    pre_and_init_feature_aabb(*ctx.tet_mesh, sf_mesh.aabb_wrapper);
  }
}

void check_cols(const Eigen::MatrixXd& M, int cols, const char* name) {
  if (M.cols() != cols)
    throw std::invalid_argument(std::string(name) + " must have " +
                                std::to_string(cols) + " columns, got " +
                                std::to_string(M.cols()));
}
void check_cols(const Eigen::MatrixXi& M, int cols, const char* name) {
  if (M.cols() != cols)
    throw std::invalid_argument(std::string(name) + " must have " +
                                std::to_string(cols) + " columns, got " +
                                std::to_string(M.cols()));
}
void check_indices(const Eigen::MatrixXi& M, int num_vertices,
                   const char* name) {
  if (M.size() == 0) return;
  if (M.minCoeff() < 0 || M.maxCoeff() >= num_vertices)
    throw std::invalid_argument(std::string(name) +
                                " contains out-of-range vertex indices");
}

}  // namespace internal
}  // namespace libmat
