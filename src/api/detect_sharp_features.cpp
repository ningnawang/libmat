#include "libmat/detect_sharp_features.h"

#include "api/internal/surface_context.h"

namespace libmat {

void detect_sharp_features(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                           double thres_concave, double thres_convex,
                           Eigen::MatrixXi& E, Eigen::VectorXi& E_type,
                           Eigen::VectorXi& E_line, Eigen::VectorXi& corners,
                           Eigen::VectorXi& corner_type) {
  internal::SurfaceContext ctx;
  internal::build_surface_context(V, F, true /*detect_features*/,
                                  thres_concave, thres_convex, ctx);
  const TetMesh& tet_mesh = *ctx.tet_mesh;

  // FeatureLine ids of sharp and concave lines are numbered independently,
  // shift concave line ids so that they are unique
  const int num_se_lines = tet_mesh.se_lines.size();
  const auto& fes = tet_mesh.feature_edges;
  E.resize(fes.size(), 2);
  E_type.resize(fes.size());
  E_line.resize(fes.size());
  for (int i = 0; i < (int)fes.size(); ++i) {
    E(i, 0) = fes[i].t2vs_group[0];
    E(i, 1) = fes[i].t2vs_group[1];
    E_type(i) = fes[i].type == EdgeType::CE ? 2 : 1;
    E_line(i) = fes[i].get_fl_id() + (fes[i].type == EdgeType::CE ? num_se_lines : 0);
  }

  std::vector<int> c, ct;
  for (const int v : tet_mesh.corners_se_tet) c.push_back(v), ct.push_back(1);
  for (const int v : tet_mesh.corners_ce_tet) {
    if (tet_mesh.corners_se_tet.count(v)) continue;
    c.push_back(v), ct.push_back(2);
  }
  corners = Eigen::Map<Eigen::VectorXi>(c.data(), c.size());
  corner_type = Eigen::Map<Eigen::VectorXi>(ct.data(), ct.size());
}

}  // namespace libmat
