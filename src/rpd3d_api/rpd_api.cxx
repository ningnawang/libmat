#include "rpd_api.h"

#include <assert.h>

#include "voronoi.h"

RPD3D_GPU::~RPD3D_GPU() {
  delete tet_mesh;
  delete sf_mesh;
  delete params;
  delete all_medial_spheres;
}

void RPD3D_GPU::init(const TetMesh* _tet_mesh, const SurfaceMesh* _sf_mesh,
                     const Parameter* _params) {
  this->tet_mesh = _tet_mesh;
  this->sf_mesh = _sf_mesh;
  this->params = _params;
  num_itr_rpd = 0;
  n_site = -1;
  site_k = -1;
}

void RPD3D_GPU::set_spheres(std::vector<MedialSphere>* _all_medial_spheres) {
  this->all_medial_spheres = _all_medial_spheres;
}

// ------------------------------------------------------------------------------------
// public functions
// ------------------------------------------------------------------------------------
void RPD3D_GPU::calculate() {
  // compute RT and mark valid spheres
  std::set<int> valid_sphere_ids;
  generate_RT_CGAL_and_mark_valid_spheres(
      *this->params, *this->all_medial_spheres, this->rt, valid_sphere_ids);
  if (is_debug)
    printf(
        "generate_RT_CGAL_and_mark_valid_spheres done, valid_sphere_ids: %zu\n",
        valid_sphere_ids.size());

  // //////////////////////
  // // load spheres_N in spheres_changed
  // std::set<int> spheres_changed;
  // load_changed_spheres(num_itr_rpd, all_medial_spheres, spheres_changed);
  // if (spheres_changed.empty()) {
  //   printf("[Partial RPD] no sphere changed\n");
  //   return;
  // }
  // if (is_debug)
  //   printf("[Partial RPD] spheres_changed: %ld \n", spheres_changed.size());

  //////////////////////
  // Note: We only care about (spheres_N + 1ring) stored in spheres_and_1rings.
  //       To cut pcells of 1ring we need 2rings as well,
  //       map_site2msphere stores (spheres_N + 1ring + 2ring)
  //       but we don't care about 2rings RPDs.
  //
  //       During RPD, we will save (sphere_N + 1ring) in site_flags,
  //       and save (sphere_N + 1ring + 2ring) in site and n_site
  //
  // id matches MedialSphere::id
  //
  // TODO: update this func
  std::vector<int> map_site2msphere;  // (spheres_N + 1ring + 2ring)
  std::set<int> spheres_and_1rings;   // (spheres_N + 1ring)
  this->site_k = get_RT_partial_spheres_and_neighbors(
      valid_sphere_ids, rt, *this->all_medial_spheres, map_site2msphere,
      spheres_and_1rings, site_knn, is_debug);
  if (map_site2msphere.empty() || site_k == -1) {
    printf("[Partial RPD] empty site \n");
    return;
  }

  if (is_debug)
    printf(
        "[Partial RPD] spheres_and_1rings: %ld, map_site2msphere (2ring): "
        "%ld\n",
        spheres_and_1rings.size(), map_site2msphere.size());

  // printf("[Partial RPD] spheres_changed: (");
  // for (const auto& csid : spheres_changed) printf("%d, ", csid);
  // printf(")\n");
  // printf("[Partial RPD] spheres_and_1rings: (");
  // for (const auto& csid : spheres_and_1rings) printf("%d, ", csid);
  // printf(")\n");

  //////////////////////
  // For debug: load all_medial_spheres
  //
  // for (const auto& msphere : all_medial_spheres) {
  //   map_site2msphere.push_back(msphere.id);
  //   spheres_and_1rings.insert(msphere.id);
  // }
  // n_site = map_site2msphere.size();
  // // compute site KNN
  // params.site_k = get_RT_vertex_neighbors(rt, n_site, site_knn);

  //////////////////////
  // load to sites:       (sphere_N + 1ring + 2ring) -> map_site2msphere
  // load to site_flags:  (sphere_N + 1ring) -> spheres_and_1rings
  this->load_partial_spheres_to_sites(map_site2msphere, spheres_and_1rings);
  printf("load_partial_spheres_to_sites done, n_site: %d\n", n_site);

  //////////////////////
  // load all tets, let map_tet_new2orig empty
  const std::vector<int> map_tet_new2orig;

  //////////////////////
  // compute parial RPD using CUDA
  // make sure load_tet_adj_info() was called ahead
  //
  // site_flags: store (sphere_N + 1ring)
  // site:       store (sphere_N + 1ring + 2ring)
  assert(!this->tet_mesh->v_adjs.empty() && !this->tet_mesh->e_adjs.empty() &&
         !this->tet_mesh->f_adjs.empty() && !this->tet_mesh->f_ids.empty());
  std::vector<ConvexCellHost> cells_partials = compute_clipped_voro_diagram_GPU(
      this->num_itr_rpd, this->tet_mesh->tet_vertices,
      this->tet_mesh->tet_indices, this->tet_mesh->v2tets,
      this->tet_mesh->v_adjs, this->tet_mesh->e_adjs, this->tet_mesh->f_adjs,
      this->tet_mesh->f_ids, this->site, this->n_site, this->site_weights,
      this->site_flags, this->site_knn, this->site_k, site_cell_vol,
      site_is_transposed, 1 /*nb_iter*/, all_medial_spheres->size());
  if (is_debug) printf("compute compute_clipped_voro_diagram_GPU done \n");

  // update cells_partials:
  // 1. ConvexCellHost::voro_id
  // 2. ConvexCellHost::tet_id
  // 3. ConvexCellHost::clip_id2_data_trans
  // to match MedialSphere::id and TetMesh::tet_indices
  // (sphere + 1ring + 2ring)
  update_convex_cells_voro_and_tet_ids(map_site2msphere, map_tet_new2orig,
                                       cells_partials, is_debug);
  if (is_debug) printf("#cells_partials: %zu \n", cells_partials.size());

  // re-assign ConvexCellHost::id
  // merge old convex_cells and current convex_cells
  // (keep only sphere + 1ring, NO 2ring)
  std::vector<ConvexCellHost> merged_cells =
      merge_convex_cells(valid_sphere_ids, spheres_and_1rings, powercells,
                         cells_partials, is_debug);
  powercells.clear();
  powercells = merged_cells;
  this->num_itr_rpd++;
}

void RPD3D_GPU::calculate_partial(int& num_itr_global, int& num_sphere_added,
                                  const bool is_given_all_tets) {
  if (num_sphere_added == 0) {
    printf("No sphere added, not run RPD\n");
    return;
  }
  printf("-------- Cal RPD Partial -- num_itr_global %d, num_sphere_added:%d\n",
         num_itr_global, num_sphere_added);

  // // remove duplicated medial spheres
  // remove_duplicated_medial_spheres(all_medial_spheres);
  // compute RT and mark valid spheres
  // also will update
  // 1. MedialSphere::rt_neigh_ids_prev
  // 2. MedialSphere::rt_neigh_ids
  std::set<int> valid_sphere_ids;
  generate_RT_CGAL_and_mark_valid_spheres(
      *this->params, *this->all_medial_spheres, this->rt, valid_sphere_ids);
  if (is_debug)
    printf(
        "generate_RT_CGAL_and_mark_valid_spheres done, valid_sphere_ids: %zu\n",
        valid_sphere_ids.size());

  //////////////////////
  // load spheres_N in spheres_changed
  // spheres_changed does not contains sphere in spheres_invalid
  // but contains neighbors of spheres in spheres_invalid
  std::set<int> spheres_changed;
  // deleted spheres
  // this will make sure load_partial_tet_given_spheres() load
  // partial tets for these deleted spheres
  std::set<int> spheres_invalid;
  load_changed_spheres(num_itr_global, *this->all_medial_spheres,
                       spheres_changed, spheres_invalid, is_debug);
  if (spheres_changed.empty()) {
    printf("[Partial RPD] no sphere changed\n");
    return;
  }
  // if (is_debug)
  printf("[Partial RPD] spheres_changed: %ld, spheres_invalid: %ld \n",
         spheres_changed.size(), spheres_invalid.size());

  //////////////////////
  // Note: We only care about (spheres_N + 1ring) stored in spheres_and_1rings.
  //       To cut pcells of 1ring we need 2rings as well,
  //       map_site2msphere stores (spheres_N + 1ring + 2ring)
  //       but we don't care about 2rings RPDs.
  //
  //       During RPD, we will save (sphere_N + 1ring) in site_flags,
  //       and save (sphere_N + 1ring + 2ring) in site and n_site
  //
  // id matches MedialSphere::id
  std::vector<int> map_site2msphere;  // (spheres_N + 1ring + 2ring)
  std::set<int> spheres_and_1rings;   // (spheres_N + 1ring)
  this->site_k = get_RT_partial_spheres_and_neighbors(
      spheres_changed, this->rt, *this->all_medial_spheres, map_site2msphere,
      spheres_and_1rings, this->site_knn, this->is_debug);
  if (map_site2msphere.empty() || this->site_k == -1) {
    printf("[Partial RPD] empty site \n");
    return;
  }

  if (is_debug)
    printf(
        "[Partial RPD] spheres_and_1rings: %ld, map_site2msphere (2ring): "
        "%ld\n",
        spheres_and_1rings.size(), map_site2msphere.size());

  // printf("[Partial RPD] spheres_changed: (");
  // for (const auto& csid : spheres_changed) printf("%d, ", csid);
  // printf(")\n");
  // printf("[Partial RPD] spheres_and_1rings: (");
  // for (const auto& csid : spheres_and_1rings) printf("%d, ", csid);
  // printf(")\n");

  //////////////////////
  // For debug: load all_medial_spheres
  //
  // for (const auto& msphere : all_medial_spheres) {
  //   map_site2msphere.push_back(msphere.id);
  //   spheres_and_1rings.insert(msphere.id);
  // }
  // n_site = map_site2msphere.size();
  // // compute site KNN
  // params.site_k = get_RT_vertex_neighbors(rt, n_site, site_knn);

  //////////////////////
  // load to sites:       (sphere_N + 1ring + 2ring) -> map_site2msphere
  // load to site_flags:  (sphere_N + 1ring) -> spheres_and_1rings
  this->load_partial_spheres_to_sites(map_site2msphere, spheres_and_1rings);
  printf("load_partial_spheres_to_sites done, n_site: %d\n", this->n_site);

  //////////////////////
  // select tets relates to map_site2msphere
  // (sphere + 1ring), no need 2ring
  std::vector<int> partial_tet_indices;
  std::vector<int> partial_tet_fids, partial_tet_f_adjs; /*no need map*/
  std::vector<int> map_tet_new2orig;
  // for initial run
  if (is_given_all_tets) {
    partial_tet_indices = this->tet_mesh->tet_indices;
    partial_tet_fids = this->tet_mesh->f_ids;
    partial_tet_f_adjs = this->tet_mesh->f_adjs;
  } else {
    // remember to load tets from deleted spheres!
    std::set<int> target_spheres = spheres_and_1rings;
    target_spheres.insert(spheres_invalid.begin(), spheres_invalid.end());
    load_partial_tet_given_spheres(
        this->tet_mesh->tet_indices, this->tet_mesh->f_ids,
        this->tet_mesh->f_adjs, target_spheres, *this->all_medial_spheres,
        partial_tet_indices, map_tet_new2orig, partial_tet_fids,
        partial_tet_f_adjs, is_debug);
  }

  if (is_debug)
    printf("partial #tet: %ld, orignal #tet: %ld\n",
           partial_tet_indices.size() / 4,
           this->tet_mesh->tet_indices.size() / 4);

  //////////////////////
  // compute parial RPD using CUDA
  // make sure load_tet_adj_info() was called ahead
  //
  // site_flags: store (sphere_N + 1ring)
  // site:       store (sphere_N + 1ring + 2ring)
  assert(!this->tet_mesh->v_adjs.empty() && !this->tet_mesh->e_adjs.empty() &&
         !this->tet_mesh->f_adjs.empty() && !this->tet_mesh->f_ids.empty());
  assert(partial_tet_indices.size() == partial_tet_fids.size());
  std::vector<ConvexCellHost> cells_partials = compute_clipped_voro_diagram_GPU(
      num_itr_global, this->tet_mesh->tet_vertices, partial_tet_indices,
      this->tet_mesh->v2tets, this->tet_mesh->v_adjs, this->tet_mesh->e_adjs,
      partial_tet_f_adjs, partial_tet_fids, this->site, this->n_site,
      this->site_weights, this->site_flags, this->site_knn, this->site_k,
      site_cell_vol, site_is_transposed, 1 /*nb_iter*/,
      all_medial_spheres->size());

  if (is_debug) printf("compute compute_clipped_voro_diagram_GPU done \n");
  if (is_debug)
    printf("#map_site2msphere: %zu, #map_tet_new2orig: %zu\n",
           map_site2msphere.size(), map_tet_new2orig.size());

  // update cells_partials:
  // 1. ConvexCellHost::voro_id
  // 2. ConvexCellHost::tet_id
  // 3. ConvexCellHost::clip_id2_data_trans
  // to match MedialSphere::id and TetMesh::tet_indices
  // (sphere + 1ring + 2ring)
  update_convex_cells_voro_and_tet_ids(map_site2msphere, map_tet_new2orig,
                                       cells_partials, is_debug);

  if (is_debug)
    printf("#powercells: %zu, #cells_partials: %zu \n", powercells.size(),
           cells_partials.size());

  // re-assign ConvexCellHost::id
  // merge old convex_cells and current convex_cells
  // (keep only sphere + 1ring, NO 2ring)
  std::vector<ConvexCellHost> merged_cells =
      merge_convex_cells(valid_sphere_ids, spheres_and_1rings, powercells,
                         cells_partials, is_debug);
  powercells.clear();
  powercells = merged_cells;

  this->update_spheres_power_cells();
  num_itr_global++;
  num_sphere_added = 0;  // reset
}

// ------------------------------------------------------------------------------------
// private functions
// ------------------------------------------------------------------------------------
void RPD3D_GPU::update_spheres_power_cells(bool is_compute_se_sfids) {
  if (is_compute_se_sfids) {
    // for SE spheres and corners
    // also update covered sf_mesh fids
    //
    // call before update_power_cells()
    // to make sure TangentPlane::sf_fids_covered will be updated
    //
    // TODO: move this function to a better place!!!!
    update_se_tangent_planes(*this->sf_mesh, this->tet_mesh->feature_edges,
                             *this->all_medial_spheres,
                             true /*is_clear_existing*/);
  }

  // update powercells in multi-threads
  update_power_cells(*this->sf_mesh, this->powercells,
                     *this->all_medial_spheres, this->tet_mesh->tet_es2fe_map,
                     this->is_debug);
}

// will update
// 1. this->site
// 2. this->site_weights
// 3. this->site_flags
// 4. this->n_site
void RPD3D_GPU::load_partial_spheres_to_sites(
    const std::vector<int>& map_site2msphere,
    const std::set<int>& spheres_and_1rings) {
  // (map_site2msphere: site id to MedialSphere::id)
  if (map_site2msphere.empty() || spheres_and_1rings.empty()) return;
  this->site.clear();
  this->site_weights.clear();
  this->site_flags.clear();
  this->n_site = map_site2msphere.size();
  int dim = 4;
  std::cout << "n_site: " << this->n_site << ", dim: " << dim << std::endl;

  this->site.resize(this->n_site * 3);
  this->site_weights.resize(this->n_site);
  this->site_flags.resize(this->n_site);
  for (int i = 0; i < this->n_site; ++i) {
    const int id = map_site2msphere[i];
    const auto& msphere = this->all_medial_spheres->at(id);
    assert(id == msphere.id);
    assert(!msphere.is_deleted);
    this->site[i] = msphere.center[0];
    this->site[i + this->n_site] = msphere.center[1];
    this->site[i + (this->n_site << 1)] = msphere.center[2];

    // add weight (sq_radius) info
    this->site_weights[i] = std::pow(float(msphere.radius), 2);

    // add flag for sphere
    this->site_flags[i] = SiteFlag::no_flag;
    if (spheres_and_1rings.find(msphere.id) != spheres_and_1rings.end())
      this->site_flags[i] = SiteFlag::is_selected;

    // printf("site %d has sphere: (%lf %lf %lf %lf) with flag %d \n", i,
    //        msphere.center[0], msphere.center[1], msphere.center[2],
    //        msphere.radius, site_flags[i]);
  }
}

// update convex_cells_new:
// 1. ConvexCellHost::voro_id
// 2. ConvexCellHost::tet_id
// 3. ConvexCellHost::clip_id2_data_trans
// to match MedialSphere::id and TetMesh::tet_indices
void RPD3D_GPU::update_convex_cells_voro_and_tet_ids(
    const std::vector<int>& map_site2msphere,
    const std::vector<int>& map_tet_new2orig,
    std::vector<ConvexCellHost>& convex_cells_new, bool is_debug) {
  if (map_site2msphere.empty() && map_tet_new2orig.empty()) return;

  for (auto& convex_cell : convex_cells_new) {
    if (!map_site2msphere.empty()) {
      // printf(
      //     "updating convex_cell.voro_id %d, cell_id %d, tet_id %d, status %d,
      //     " "thread %d\n", convex_cell.voro_id, convex_cell.id,
      //     convex_cell.tet_id, convex_cell.status, convex_cell.thread_id);
      /////////////
      // 1. update ConvexCellHost::voro_id
      assert(convex_cell.voro_id > -1 &&
             convex_cell.voro_id < map_site2msphere.size());
      convex_cell.voro_id = map_site2msphere.at(convex_cell.voro_id);
      /////////////
      // 3. update ConvexCellHost::clip_id2_data_trans
      // update halfplane of two seeds
      // where we store [site_id_min, site_id_max]
      FOR(i, convex_cell.nb_p) {
        cint2& clip_id2 = convex_cell.clip_id2_data_trans[i];
        if (clip_id2.y == -1)
          continue;  // plane inherit from tets, don't update
        assert(clip_id2.x > -1 && clip_id2.x < map_site2msphere.size());
        assert(clip_id2.y > -1 && clip_id2.y < map_site2msphere.size());
        clip_id2.x = map_site2msphere.at(clip_id2.x);
        clip_id2.y = map_site2msphere.at(clip_id2.y);
      }
    }

    /////////////
    // 2. update ConvexCellHost::tet_id
    if (!map_tet_new2orig.empty()) {
      convex_cell.tet_id = map_tet_new2orig.at(convex_cell.tet_id);
    }

    // if (convex_cell.voro_id == 377)
    //   printf(
    //       " ---------- [UpdateConvexCells] convex_cell_new %d has voro_id %d,
    //       " "tet_id %d \n", convex_cell.id, convex_cell.voro_id,
    //       convex_cell.tet_id);
  }
}

std::vector<ConvexCellHost> RPD3D_GPU::merge_convex_cells(
    const std::set<int>& valid_sphere_ids,
    const std::set<int>& spheres_and_1rings,
    const std::vector<ConvexCellHost>& convex_cells_prev,
    const std::vector<ConvexCellHost>& convex_cells_new, bool is_debug) {
  std::vector<ConvexCellHost> merged_convex_cells;
  for (const auto& cell_new : convex_cells_new) {
    // filter by spheres + 1ring, No 2ring
    if (spheres_and_1rings.find(cell_new.voro_id) == spheres_and_1rings.end())
      continue;

    merged_convex_cells.push_back(cell_new);
    merged_convex_cells.back().id = merged_convex_cells.size() - 1;

    // if (cell_new.voro_id == 377)
    //   printf("[Merge ConvexCells] cell_new %d has voro_id %d \n",
    //          merged_convex_cells.back().id,
    //          merged_convex_cells.back().voro_id);
  }

  for (const auto& cell_prev : convex_cells_prev) {
    // if prev voro_id became invalid, then skip
    if (valid_sphere_ids.find(cell_prev.voro_id) == valid_sphere_ids.end())
      continue;

    // do not store cells of (spheres + 1ring)
    if (spheres_and_1rings.find(cell_prev.voro_id) != spheres_and_1rings.end())
      continue;
    // store prev
    merged_convex_cells.push_back(cell_prev);
    merged_convex_cells.back().id = merged_convex_cells.size() - 1;

    // if (cell_prev.voro_id == 377) {
    //   printf("[Merge ConvexCells] cell_prev %d has voro_id %d, tet_id: %d\n",
    //          merged_convex_cells.back().id,
    //          merged_convex_cells.back().voro_id,
    //          merged_convex_cells.back().tet_id);
    //   if (merged_convex_cells.back().tet_id == 900) {
    //     merged_convex_cells.back().print_info();
    //   }
    // }
  }

  printf("[Merge ConvexCells] #prev: %zu, #new: %zu, #merged: %zu\n",
         convex_cells_prev.size(), convex_cells_new.size(),
         merged_convex_cells.size());
  return merged_convex_cells;
}

// see load_tet_adj_info() for old_tet_indices and old_tet_fids
void RPD3D_GPU::load_partial_tet_given_spheres(
    const std::vector<int>& old_tet_indices,
    const std::vector<int>& old_tet_fids,
    const std::vector<int>& old_tet_f_adjs, const std::set<int>& given_spheres,
    const std::vector<MedialSphere>& all_medial_spheres,
    std::vector<int>& partial_tet_indices, std::vector<int>& map_tet_new2orig,
    std::vector<int>& partial_tet_fids, std::vector<int>& partial_tet_f_adjs,
    bool is_debug) {
  if (is_debug) {
    printf("[Partial Tet] given_spheres: (");
    for (const auto& csid : given_spheres) printf("%d, ", csid);
    printf(")\n");
  }

  // map_tet_new2orig stores new tet_ids to orignal tet_ids
  // see load_tet_adj_info()
  assert(old_tet_indices.size() == old_tet_fids.size());
  assert(old_tet_f_adjs.size() == old_tet_fids.size());
  partial_tet_indices.clear();
  map_tet_new2orig.clear();
  partial_tet_fids.clear();
  partial_tet_f_adjs.clear();
  std::set<int> tet_old_visited;
  for (const auto& all_id : given_spheres) {
    const auto& msphere = all_medial_spheres.at(all_id);
    for (const auto& tet_id_old : msphere.pcell.tet_ids) {
      if (tet_old_visited.find(tet_id_old) != tet_old_visited.end()) continue;
      tet_old_visited.insert(tet_id_old);
      map_tet_new2orig.push_back(tet_id_old);
      for (int j = 0; j < 4; j++) {
        partial_tet_indices.push_back(old_tet_indices.at(tet_id_old * 4 + j));
        partial_tet_fids.push_back(old_tet_fids.at(tet_id_old * 4 + j));
        partial_tet_f_adjs.push_back(old_tet_f_adjs.at(tet_id_old * 4 + j));
      }

      // int tet_id_new = map_tet_new2orig.size() - 1;
      // if (tet_id_old == 1471 || tet_id_new == 716) {
      // printf(
      //     "[Partial Tet] orig_tet %d (%d,%d,%d,%d) updated to partial_tet "
      //     "%d (%d,%d,%d,%d) \n",
      //     tet_id_old, old_tet_indices.at(tet_id_old * 4),
      //     old_tet_indices.at(tet_id_old * 4 + 1),
      //     old_tet_indices.at(tet_id_old * 4 + 2),
      //     old_tet_indices.at(tet_id_old * 4 + 3), tet_id_new,
      //     partial_tet_indices.at(tet_id_new * 4),
      //     partial_tet_indices.at(tet_id_new * 4 + 1),
      //     partial_tet_indices.at(tet_id_new * 4 + 2),
      //     partial_tet_indices.at(tet_id_new * 4 + 3));
      // }
    }
  }
  assert(map_tet_new2orig.size() == partial_tet_indices.size() / 4);
  assert(partial_tet_indices.size() == partial_tet_fids.size());
}
