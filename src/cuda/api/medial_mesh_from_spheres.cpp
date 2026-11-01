// CUDA implementation of libmat::medial_mesh_from_spheres(): restricted power
// diagram on the GPU (RPD3D_GPU), then medial mesh extraction as in
// MATFP/MATTopo.

#include "libmat/medial_mesh_from_spheres.h"

#include <memory>
#include <stdexcept>

#include "api/internal/surface_context.h"
#include "libmat/cuda/rpd3d/rpd_api.h"
#include "libmat/inputs/sharp_feature_detection.h"
#include "libmat/io/io.h"
#include "libmat/matbase/medial_mesh.h"
#include "libmat/matfun/shrinking.h"
#include "libmat/matfun/thinning.h"

namespace libmat {

void medial_mesh_from_spheres(const Eigen::MatrixXd& V,
                              const Eigen::MatrixXi& T,
                              const Eigen::MatrixXd& C,
                              const Eigen::VectorXd& R, Eigen::MatrixXi& E,
                              Eigen::MatrixXi& F,
                              Eigen::VectorXd& importance) {
  internal::ensure_geogram_initialized();
  internal::check_cols(V, 3, "V");
  internal::check_cols(T, 4, "T");
  internal::check_indices(T, V.rows(), "T");
  internal::check_cols(C, 3, "C");
  if (R.size() != C.rows())
    throw std::invalid_argument("R must have one entry per row of C");

  // RPD3D_GPU deletes tet_mesh, sf_mesh, params and spheres it points to,
  // so all of them are heap allocated and owned by rpd
  auto params = std::make_unique<Parameter>();
  Eigen::MatrixXd VN;
  const internal::Normalization nrm =
      internal::normalize_and_fill_params(V, *params, VN);

  auto tet_mesh = std::make_unique<TetMesh>("");
  tet_mesh->tet_vertices.resize(3 * VN.rows());
  for (int i = 0; i < VN.rows(); ++i)
    for (int j = 0; j < 3; ++j) tet_mesh->tet_vertices[3 * i + j] = VN(i, j);
  tet_mesh->tet_indices.resize(4 * T.rows());
  for (int i = 0; i < T.rows(); ++i)
    for (int j = 0; j < 4; ++j) tet_mesh->tet_indices[4 * i + j] = T(i, j);
  load_v2tets(tet_mesh->tet_vertices, tet_mesh->tet_indices, tet_mesh->v2tets);

  // boundary surface, with sf2tet stored in the "tet_vid" attribute
  std::vector<std::array<float, 3>> sf_vertices;
  std::vector<std::array<int, 3>> sf_faces;
  std::vector<int> sf2tet;
  get_surface_from_tet(tet_mesh->tet_vertices, tet_mesh->tet_indices,
                       sf_vertices, sf_faces, sf2tet);
  Eigen::MatrixXd SV(sf2tet.size(), 3);
  for (int i = 0; i < (int)sf2tet.size(); ++i) SV.row(i) = VN.row(sf2tet[i]);
  Eigen::MatrixXi SF(sf_faces.size(), 3);
  for (int i = 0; i < (int)sf_faces.size(); ++i)
    for (int j = 0; j < 3; ++j) SF(i, j) = sf_faces[i][j];

  auto sf_mesh = std::make_unique<SurfaceMesh>();
  internal::build_sf_mesh(SV, SF, sf2tet, *sf_mesh);
  pre_and_init_aabb(*sf_mesh, sf_mesh->aabb_wrapper, true /*is_reorder*/);
  sf_mesh->reload_sf2tet_vs_mapping();
  std::map<int, int> tet2sf_vs;
  std::map<int, std::set<int>> sf_vs2fids;
  load_sf_tet_mapping(*sf_mesh, tet2sf_vs, sf_vs2fids,
                      tet_mesh->tet_vs2sf_fids);
  detect_mark_sharp_features(*params, *sf_mesh, *tet_mesh,
                             false /*is_sf_file_exist*/);
  pre_and_init_feature_aabb(*tet_mesh, sf_mesh->aabb_wrapper);
  load_tet_adj_info(tet_mesh->v2tets, tet_mesh->tet_indices,
                    tet_mesh->tet_vs2sf_fids, sf_mesh->facets.nb(),
                    tet_mesh->v_adjs, tet_mesh->e_adjs, tet_mesh->f_adjs,
                    tet_mesh->f_ids);

  auto spheres = std::make_unique<std::vector<MedialSphere>>();
  spheres->reserve(C.rows());
  for (int i = 0; i < C.rows(); ++i) {
    const Eigen::RowVector3d c = nrm.to_normalized(C.row(i));
    spheres->emplace_back(i, Vector3(c[0], c[1], c[2]), R(i) * nrm.factor,
                          SphereType::T_2);
  }

  RPD3D_GPU rpd;
  std::vector<MedialSphere>& all_spheres = *spheres;
  rpd.init(tet_mesh.release(), sf_mesh.release(), params.release());
  rpd.set_spheres(spheres.release());
  rpd.calculate();
  rpd.update_spheres_power_cells(true /*is_compute_se_sfids*/);

  // MedialMesh also deletes its vertices: detach them before it goes out of
  // scope, rpd owns the spheres
  struct DetachGuard {
    MedialMesh& m;
    ~DetachGuard() { m.vertices = nullptr; }
  };
  MedialMesh mmesh;
  DetachGuard guard{mmesh};
  mmesh.genearte_medial_spheres(all_spheres);
  mmesh.generate_medial_faces();
  mmesh.generate_medial_edges(*rpd.sf_mesh, false);
  Thinning::load_all_mat_face_importance_globally(mmesh);

  std::vector<aint2> edges;
  for (const auto& e : mmesh.edges)
    if (!e.is_deleted) edges.push_back(e.vertices_);
  E.resize(edges.size(), 2);
  for (int i = 0; i < (int)edges.size(); ++i)
    E(i, 0) = edges[i][0], E(i, 1) = edges[i][1];

  std::vector<int> fids;
  for (int f = 0; f < (int)mmesh.faces.size(); ++f)
    if (!mmesh.faces[f].is_deleted) fids.push_back(f);
  F.resize(fids.size(), 3);
  importance.resize(fids.size());
  for (int i = 0; i < (int)fids.size(); ++i) {
    const auto& face = mmesh.faces[fids[i]];
    for (int j = 0; j < 3; ++j) F(i, j) = face.vertices_[j];
    importance(i) = face.importance;
  }
}

}  // namespace libmat
