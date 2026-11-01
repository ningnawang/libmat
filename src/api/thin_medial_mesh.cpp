#include "libmat/thin_medial_mesh.h"

#include <set>
#include <stdexcept>
#include <vector>

#include "api/internal/surface_context.h"
#include "libmat/matbase/medial_mesh.h"
#include "libmat/matfun/thinning.h"

namespace libmat {

void thin_medial_mesh(const Eigen::MatrixXd& C, const Eigen::VectorXd& R,
                      const Eigen::MatrixXi& E, const Eigen::MatrixXi& F,
                      const Eigen::VectorXd& importance, double threshold,
                      Eigen::MatrixXi& E_out, Eigen::MatrixXi& F_out) {
  internal::check_cols(C, 3, "C");
  if (R.size() != C.rows())
    throw std::invalid_argument("R must have one entry per row of C");
  if (E.size() > 0) internal::check_cols(E, 2, "E");
  internal::check_cols(F, 3, "F");
  internal::check_indices(E, C.rows(), "E");
  internal::check_indices(F, C.rows(), "F");
  if (importance.size() != F.rows())
    throw std::invalid_argument("importance must have one entry per face");

  // MedialMesh owns (and deletes) its vertices
  auto* spheres = new std::vector<MedialSphere>();
  spheres->reserve(C.rows());
  for (int i = 0; i < C.rows(); ++i)
    spheres->emplace_back(i, Vector3(C(i, 0), C(i, 1), C(i, 2)), R(i),
                          SphereType::T_2);
  MedialMesh mat(*spheres);

  // edges and faces, skipping duplicates
  for (int e = 0; e < E.rows(); ++e) {
    int eid = -1;
    if (E(e, 0) == E(e, 1) || mat.Edge(E(e, 0), E(e, 1), eid)) continue;
    mat.create_edge(E(e, 0), E(e, 1), 1 /*dup_cnt*/);
  }
  std::set<aint3> face_set;
  for (int f = 0; f < F.rows(); ++f) {
    aint3 vs = {{F(f, 0), F(f, 1), F(f, 2)}};
    std::sort(vs.begin(), vs.end());
    if (vs[0] == vs[1] || vs[1] == vs[2] || face_set.count(vs)) continue;
    face_set.insert(vs);
    const int fid = mat.create_face(vs, 1 /*dup_cnt*/, -1,
                                    true /*is_auto_create_edge*/);
    mat.faces[fid].importance = importance(f);
  }

  // medial tets: four vertices whose four triangles are all medial faces
  std::vector<std::set<int>> neighbors(C.rows());
  for (const auto& edge : mat.edges) {
    neighbors[edge.vertices_[0]].insert(edge.vertices_[1]);
    neighbors[edge.vertices_[1]].insert(edge.vertices_[0]);
  }
  std::set<aint4> tets;
  for (const auto& vs : face_set) {
    for (const int d : neighbors[vs[0]]) {
      if (d <= vs[2]) continue;  // each tet is found from its smallest face
      if (!neighbors[vs[1]].count(d) || !neighbors[vs[2]].count(d)) continue;
      if (face_set.count({{vs[0], vs[1], d}}) &&
          face_set.count({{vs[0], vs[2], d}}) &&
          face_set.count({{vs[1], vs[2], d}}))
        tets.insert({{vs[0], vs[1], vs[2], d}});
    }
  }
  for (const auto& t : tets)
    mat.create_tet(t, 1 /*dup_cnt*/, -1, false /*is_auto_create_face*/);

  if (!mat.tets.empty())
    Thinning::prune(*spheres, mat, threshold, false /*is_dup_cnt*/,
                    true /*is_import_given*/);

  std::vector<aint2> edges_out;
  for (const auto& edge : mat.edges)
    if (!edge.is_deleted) edges_out.push_back(edge.vertices_);
  std::vector<aint3> faces_out;
  for (const auto& face : mat.faces)
    if (!face.is_deleted) faces_out.push_back(face.vertices_);
  E_out.resize(edges_out.size(), 2);
  for (int i = 0; i < (int)edges_out.size(); ++i)
    E_out(i, 0) = edges_out[i][0], E_out(i, 1) = edges_out[i][1];
  F_out.resize(faces_out.size(), 3);
  for (int i = 0; i < (int)faces_out.size(); ++i)
    for (int j = 0; j < 3; ++j) F_out(i, j) = faces_out[i][j];
}

}  // namespace libmat
