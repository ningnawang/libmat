#pragma once

#include "input_types.h"
#include "medial_sphere.h"
#include "triangulation.h"

class RPD3D_GPU {
 public:
  RPD3D_GPU(){};
  ~RPD3D_GPU();
  void init(const TetMesh* _tet_mesh, const SurfaceMesh* _sf_mesh,
            const Parameter* _params);
  void set_spheres(std::vector<MedialSphere>* _all_medial_spheres);

 public:
  void calculate();
  void calculate_partial(int& num_itr_global, int& num_sphere_added,
                         const bool is_given_all_tets);

 public:
  void update_spheres_power_cells(bool is_compute_se_sfids = true);
  void load_partial_spheres_to_sites(const std::vector<int>& map_site2msphere,
                                     const std::set<int>& spheres_and_1rings);
  void update_convex_cells_voro_and_tet_ids(
      const std::vector<int>& map_site2msphere,
      const std::vector<int>& map_tet_new2orig,
      std::vector<ConvexCellHost>& cells_partials, bool is_debug);

  std::vector<ConvexCellHost> merge_convex_cells(
      const std::set<int>& valid_sphere_ids,
      const std::set<int>& spheres_and_1rings,
      const std::vector<ConvexCellHost>& convex_cells_prev,
      const std::vector<ConvexCellHost>& convex_cells_new, bool is_debug);

  void load_partial_tet_given_spheres(
      const std::vector<int>& old_tet_indices,
      const std::vector<int>& old_tet_fids,
      const std::vector<int>& old_tet_f_adjs,
      const std::set<int>& given_spheres,
      const std::vector<MedialSphere>& all_medial_spheres,
      std::vector<int>& partial_tet_indices, std::vector<int>& map_tet_new2orig,
      std::vector<int>& partial_tet_fids, std::vector<int>& partial_tet_f_adjs,
      bool is_debug);

 public:
  std::vector<ConvexCellHost> powercells;
  bool is_debug = false;

 public:
  const TetMesh* tet_mesh;
  const SurfaceMesh* sf_mesh;
  const Parameter* params;
  std::vector<MedialSphere>* all_medial_spheres;

  bool site_is_transposed;
  bool is_given_all_tets;

  RegularTriangulationNN rt;
  std::vector<float> site;
  std::vector<float> site_weights;
  std::vector<uint> site_flags;
  std::vector<int> site_knn;
  std::vector<float> site_cell_vol;

  int num_itr_rpd;
  int n_site;  // number of sites
  int site_k;  // maximum number of neighboring spheres
};
