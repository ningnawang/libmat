#include "fix_topo.h"

#include "queue"
/////////////////////////////////////////////
// Helper Functions
/////////////////////////////////////////////

/**
 * max_fid = sf_mesh.facets.nb() -1
 * Defines the maximum
 */
int get_random_fid_in_range(const std::set<int> bfids, const int max_fid = -1) {
  while (true) {
    int tmp_idx = RANDOM_INT(0, bfids.size());
    int rand_fidx = *std::next(bfids.begin(), tmp_idx);
    if (max_fid != -1 && rand_fidx > max_fid) continue;
    return rand_fidx;
  }
}

v2int get_random_v2fid(const std::vector<v2int> surf_v2fids) {
  int tmp_idx = RANDOM_INT(0, surf_v2fids.size());
  v2int rand_v2fid = *std::next(surf_v2fids.begin(), tmp_idx);
  return rand_v2fid;
}

/**
 * max_fid = sf_mesh.facets.nb() -1
 */
int get_fid_max_to_point(const GEO::Mesh &sf_mesh, const std::set<int> bfids,
                         const Vector3 &scenter, const int max_fid = -1) {
  int fid_max_dist = -1;
  double max_dist = -1;
  for (auto &fidx : bfids) {
    if (max_fid != -1 && fidx > max_fid) continue;
    Vector3 fcentroid = GEO::Geom::mesh_facet_center(sf_mesh, fidx);
    double dist = GEO::Geom::distance(scenter, fcentroid);
    if (dist > max_dist) {
      fid_max_dist = fidx;
      max_dist = dist;
    }
  }
  return fid_max_dist;
}

/**
 * @brief Get the CC boundary fids object
 *
 * @param cell_to_tfids
 * @param cc_cells
 * @param cc_boundary_fids return boundary fids of each CC, first
 * #sf_mesh.facets.nb()-1 matches GEO::Mesh, later are unique fids from orignal
 * tet
 */
void get_CC_boundary_fids(const std::map<int, std::set<int>> &cell_to_tfids,
                          const std::vector<std::set<int>> &cc_cells,
                          std::vector<std::set<int>> &cc_boundary_fids) {
  cc_boundary_fids.clear();
  // save boundary fids for each CC
  // Note: some fids match sf_mesh, but some are just assigned uniquely in
  // function load_num_adjacent_cells_and_ids() and match nothing
  std::set<int> b_fids;
  for (const auto &one_cc_cells : cc_cells) {
    b_fids.clear();
    for (int cell_id : one_cc_cells) {
      const auto &fids = cell_to_tfids.at(cell_id);
      b_fids.insert(*fids.begin());
    }
    cc_boundary_fids.push_back(b_fids);
  }
}

/**
 * @brief Check # of connected components for each power cell (of each sphere).

 * @param msphere
 * @param is_debug
 * @return true if #CC>1
 * @return false
 */
bool is_to_fix_voro_cc(MedialSphere &msphere, bool is_debug) {
  // printf("calling is_to_fix_voro_cc ... with #spheres: %zu \n",
  //        all_medial_spheres.size());

  // assert(!msphere.is_deleted);
  // // calculate number of CC
  // int num_cc = get_CC_given_neighbors<int>(msphere.pcell.cell_ids,
  //                                          msphere.pcell.cell_neighbors,
  //                                          msphere.pcell.cc_cells);

  // // TODO: remove this
  // get_CC_boundary_fids(msphere.pcell.cell_to_tfids, msphere.pcell.cc_cells,
  //                      msphere.pcell.cc_bfids);
  // get_CC_surf_vfids(msphere.pcell.cell_to_surfv2fid, msphere.pcell.cc_cells,
  //                   msphere.pcell.cc_surf_v2fids);

  assert(!msphere.pcell.cc_cells.empty());
  assert(!msphere.pcell.cc_surf_v2fids.empty());
  int num_cc = msphere.pcell.cc_cells.size();
  if (is_debug)
    printf("[CellCC] msphere: %d has num_cc: %d \n", msphere.id, num_cc);
  if (num_cc > 1) {
    msphere.pcell.topo_status = Topo_Status::high_cell_cc;
    return true;
  }
  assert(num_cc > 0);
  return false;
}

/**
 * @brief Check # of Euler Characteristics for each power cell (of each sphere).
 *
 * @param convex_cells_host
 * @param msphere
 * @param is_debug
 * @return true if #euler<1
 * @return false
 */
bool is_to_fix_voro_euler(std::vector<ConvexCellHost> &convex_cells_host,
                          MedialSphere &msphere, bool is_debug) {
  const std::set<int> &all_cell_ids = msphere.pcell.cell_ids;
  assert(all_cell_ids.size() > 0);
  // collect all convex cells
  for (uint cell_id : all_cell_ids) {
    auto &convex_cell = convex_cells_host.at(cell_id);
    assert(is_convex_cell_valid(convex_cell));
    // update convex_cell.euler
    convex_cell.cal_cell_euler();
    msphere.euler_sum += convex_cell.euler;
    msphere.num_cells += 1;
  }
  msphere.euler = msphere.euler_sum - msphere.num_cells;
  if (is_debug)
    printf(
        "[CellEuler] msphere: %d has euler: %f, euler_without cells: %f, "
        "num_cells: %d \n",
        msphere.id, msphere.euler, msphere.euler_sum, msphere.num_cells);
  // if (msphere.euler - 1 < -1.0 * SCALAR_ZERO_1) {
  if (msphere.euler < 0.5) {
    msphere.pcell.topo_status = Topo_Status::low_cell_euler;
    return true;
  }
  return false;
}

/**
 * @brief Check # of connected components for each halfplane (of each sphere)
 *
 * @param msphere
 * @param is_debug
 * @return true if #facet_cc > 1
 * @return false
 */
bool is_to_fix_facet_cc(const std::vector<MedialSphere> &all_medial_spheres,
                        MedialSphere &msphere, bool is_debug) {
  // calculate number of CC for each [msphere.id, neigh_id] pair
  int num_neigh_to_fix = 0;
  // must call update_pc_facet_cc_info() ahead
  assert(!msphere.pcell.facet_cc_cells.empty());
  assert(!msphere.pcell.facet_cc_surf_v2fids.empty());

  for (auto &pair : msphere.pcell.facet_cc_cells) {
    bool is_neigh_to_fix = false;  // check for each [msphere.id, neigh_id] pair

    // halfplane defined by [msphere.id, neigh_id]
    const int neigh_id = pair.first;
    int num_facet_cc = pair.second.size();
    if (num_facet_cc > 1 + SCALAR_ZERO_1) {
      is_neigh_to_fix = true;
      num_neigh_to_fix++;
    }
    if (!is_neigh_to_fix) continue;

    if (is_debug) {
      printf("[FacetCC] msphere %d with neigh_id: %d has num_facet_cc: %d \n",
             msphere.id, neigh_id, num_facet_cc);
      for (int i = 0; i < num_facet_cc; i++) {
        printf("[FacetCC] mspehre %d FacetCC %d has cell_ids: {", msphere.id,
               i);
        for (const auto &cell_id : pair.second.at(i)) printf("%d, ", cell_id);
        printf("]\n");
      }
    }

    // // Skip Conditions: for CAD models
    // // Checking neighbor sphere before changing topo_status
    // // 1. topo error is Topo_Status::high_facet_cc
    // // 2. two pcells cover the same concave line
    // // We do not add any new sphere. The concave part may split the
    // // halfplane facet into two connected components, but we do not want to
    // // handle this special case.
    // const auto &msphere_neigh = all_medial_spheres.at(neigh_id);
    // if (is_neigh_to_fix && !msphere.pcell.ce_covered_lvids.empty() &&
    //     !msphere_neigh.pcell.ce_covered_lvids.empty()) {
    //   if (is_debug)
    //     printf("[FacetCC] two mspheres (%d,%d) cover concave edges!!!! \n",
    //            msphere.id, neigh_id);
    //   std::unordered_set<int> ce_lines, ce_lines_neigh;
    //   // printf("[FacetCC] sphere %d cover concave line: (", msphere.id);
    //   for (const auto &code5 : msphere.pcell.ce_covered_lvids) {
    //     ce_lines.insert(code5[3]);  // ce_line_id
    //     // printf("%d, ", code5[3]);
    //   }
    //   // printf(") \n");
    //   // printf("[FacetCC] sphere %d cover concave line:
    //   (",msphere_neigh.id); for (const auto &code5_neigh :
    //   msphere_neigh.pcell.ce_covered_lvids) {
    //     ce_lines_neigh.insert(code5_neigh[3]);  // ce_line_id
    //     // printf("%d, ", code5_neigh[3]);
    //   }
    //   // printf(") \n");

    //   // sphere and neigh_sphere has common concave line ids
    //   std::vector<int> common_ce_lines;
    //   set_intersection<int>(ce_lines, ce_lines_neigh, common_ce_lines);
    //   if (!common_ce_lines.empty()) {
    //     if (is_debug) {
    //       printf("[FacetCC] two mspheres (%d,%d) cover common concave lines:
    //       (",
    //              msphere.id, neigh_id);
    //       for (const int id : common_ce_lines) {
    //         printf("%d, ", id);
    //       }
    //       printf(") \n");
    //     }
    //     // no need to fix this facet
    //     is_neigh_to_fix = false;
    //     num_neigh_to_fix--;
    //     continue;
    //   }
    // }  // if is_neigh_to_fix

    // need to fix this facet defined by [msphere.id, neigh_id]
    msphere.pcell.facet_neigh_is_fixed[neigh_id] = false;  // not fixed yet
    msphere.pcell.topo_status = Topo_Status::high_facet_cc;
  }  // for msphere.pcell.facet_cc_cells

  return num_neigh_to_fix > 0;
}

bool is_to_fix_facet_euler(std::vector<ConvexCellHost> &convex_cells_host,
                           MedialSphere &msphere, bool is_debug) {
  bool is_to_fix = false;
  for (auto &pair : msphere.pcell.facet_neigh_to_cells) {
    // one halfplane defined by [msphere.id, neigh_id]
    const int neigh_id = pair.first;
    const std::set<int> &halfplane_cells = pair.second;
    double facet_euler_sum = 0.0;
    // collect all convex cells
    for (uint cell_id : halfplane_cells) {
      auto &convex_cell = convex_cells_host.at(cell_id);
      assert(is_convex_cell_valid(convex_cell));
      // update convex_cell.euler
      float facet_euler = convex_cell.cal_halfplane_facet_euler(neigh_id);
      facet_euler_sum += facet_euler;
    }
    if (is_debug)
      printf(
          "[FacetEuler] msphere: %d with neigh_id: %d, facet_euler_sum: %f \n",
          msphere.id, neigh_id, facet_euler_sum);
    if (facet_euler_sum - 1 < -0.1) {
      printf(
          "[FacetEuler FIX] msphere: %d with neigh_id: %d, facet_euler_sum: "
          "%f, facet_euler_sum-1: %f\n",
          msphere.id, neigh_id, facet_euler_sum, facet_euler_sum - 1);
      msphere.pcell.facet_neigh_is_fixed[neigh_id] = false;  // not fixed yet
      msphere.pcell.topo_status = Topo_Status::low_facet_euler;
      is_to_fix = true;
    }
  }  // for each halfplane

  return is_to_fix;
}

/**
 * @brief Check # of connected components for each shared edge (2 halfplanes)
 *
 * @param msphere
 * @param is_debug
 * @return true if #edge_cc > 1
 * @return false
 */
// ninwang: NNNNNNNNNNNNNNNO this function is not tested!!!!
bool is_to_fix_edge_cc(MedialSphere &msphere, bool is_debug) {
  // calculate number of CC for aggregated facets
  bool is_to_fix = false;
  for (auto &pair : msphere.pcell.edge_cc_cells) {
    // two halfplane defined by:
    // 1. [msphere.id, neigh_id_min]
    // 2. [msphere.id, neigh_id_max]
    const aint2 neigh_min_max = pair.first;
    int num_edge_cc = pair.second.size();  // NNNNNo test
    printf(
        "[EdgeCC] msphere %d with two neigh_ids: [%d,%d] has num_edge_cc:%d\n",
        msphere.id, neigh_min_max[0], neigh_min_max[1], num_edge_cc);
    if (num_edge_cc > 1) {
      msphere.pcell.e_is_fixed[neigh_min_max] = false;  // not fixed yet
      msphere.pcell.topo_status = Topo_Status::high_edge_cc;
      is_to_fix = true;
      printf("[EdgeCC] CHECK THIS MODEL!");
      assert(false);
    }
  }  // for each shared edge
  return is_to_fix;
}

// If 'is_merge_to_ce = true', then the new sphere's pin is pushed to concave
// line, this sphere is created as a concave sphere.
// (will set 'is_del_near_ce = false')
//
// If 'is_merge_to_ce = false', then 'is_del_near_ce = true', will delete sphere
// if too close to concave lines.
bool add_new_sphere_given_v2fid(const int num_itr_global,
                                const SurfaceMesh &sf_mesh,
                                const TetMesh &tet_mesh,
                                const v2int v2fid_chosen,
                                std::vector<MedialSphere> &all_medial_spheres,
                                const bool is_merge_to_ce, bool is_debug) {
  // did not find a proper fid
  if (v2fid_chosen.second < 0) return false;
  bool is_del_near_ce = true;  // will be false if is_merge_to_ce = true
  bool is_del_near_se = false;
  Vector3 p = v2fid_chosen.first;

  // check if pin point close to any concave lines
  // if so, then create a concave sphere
  if (is_merge_to_ce) {
    // do not delete sphere if too close to concave line
    is_del_near_ce = false;
    Vector3 p_nearest;
    double p_sq_dist_ce;
    int p_feid = sf_mesh.aabb_wrapper.get_nearest_point_on_ce(p, p_nearest,
                                                              p_sq_dist_ce);
    if (is_debug)
      printf(
          "[FixAdd CCSphere] pin (%f,%f,%f) dist to concave edge %d "
          "p_sq_dist_ce %f\n",
          p[0], p[1], p[2], p_feid, p_sq_dist_ce);
    // since scale is [0,1000]^3
    if (p_feid != UNK_FACE && p_sq_dist_ce <= SCALAR_1) {
      // int p_fid = MedialSphere::convert_ss(p_feid);  // let it be negative
      if (is_debug)
        printf(
            "[FixAdd CCSphere] pin (%f,%f,%f) too close to concave edge %d "
            "with p_sq_dist_ce %f, create a new SphereType::T_X_c sphere \n",
            p[0], p[1], p[2], p_feid, p_sq_dist_ce);
      // create a concave sphere
      int new_sphere_id = create_new_concave_sphere_given_pin_wrapper(
          sf_mesh, tet_mesh.feature_edges, p_nearest, p_feid,
          all_medial_spheres, 1 /*SphereType::T_X_c*/, false /*is_debug*/);
      if (new_sphere_id < 0) {
        printf("[FixAdd CCSphere] failed \n");
        return false;
      } else {
        if (is_debug)
          printf("[FixAdd CCSphere] newly added sphere %d \n", new_sphere_id);
        return true;
      }
    }
  }  // if is_merge_to_ce

  // normal sphere shrinking
  is_debug = false;
  Vector3 p_normal = get_mesh_facet_normal(sf_mesh, v2fid_chosen.second);
  MedialSphere new_msphere(all_medial_spheres.size(), p, p_normal,
                           SphereType::T_2, num_itr_global);
  new_msphere.ss.p_fid = v2fid_chosen.second;
  if (new_msphere.id == 262 || new_msphere.id == 79 || new_msphere.id == 155)
    is_debug = true;
  else
    is_debug = false;
  if (!shrink_sphere(sf_mesh, sf_mesh.aabb_wrapper, sf_mesh.fe_sf_fs_pairs,
                     tet_mesh.feature_edges, new_msphere, -1 /*itr_limit*/,
                     is_del_near_ce, is_del_near_se, is_debug /*is_debug*/)) {
    // printf("[FixAdd] failed to add sphere\n");
    return false;
  }
  new_msphere.pcell.topo_status = Topo_Status::unkown;
  new_msphere.type = SphereType::T_2;
  // all_medial_spheres.push_back(new_msphere);
  bool is_good = add_new_sphere_validate(all_medial_spheres, new_msphere,
                                         true /*is_small_threshold*/);
  if (!is_good) printf("[FixAdd] failed to validate sphere\n");

  if (is_debug)
    printf(
        "[Fix Add] add new_msphere %d, is_good %d, select v2fid_chosen %d as "
        "pin\n ",
        all_medial_spheres.back().id, is_good, v2fid_chosen.second);
  return is_good;
}

/////////////////////////////////////////////
// Main Functions (new)
/////////////////////////////////////////////
void check_cc_and_euler(
    std::vector<ConvexCellHost> &convex_cells_host,
    std::vector<MedialSphere> &all_medial_spheres,
    std::set<int> &spheres_to_fix,
    bool is_debug) {  // seed_id -> # of connected components

  std::vector<bool> spheres_to_fix_tmp(all_medial_spheres.size(), false);
  auto run_thread_check = [&](int sphere_id, int check_id) {
    auto &msphere = all_medial_spheres[sphere_id];
    if (msphere.is_deleted) return;
    // skip pin invisible spheres for concave edges
    if (msphere.is_on_ce_pin()) return;
    // skip those newly added spheres to fix through
    // fix_topo_by_adding_new_sphere()
    if (msphere.pcell.topo_status == Topo_Status::unkown) return;
    // already need to fix
    if (spheres_to_fix_tmp[sphere_id]) return;
    // not enough powercell to check
    const std::set<int> &all_cell_ids = msphere.pcell.cell_ids;
    if (all_cell_ids.size() <= 1) {
      return;
    }

    // if (msphere.id == 573)
    //   is_debug = true;
    // else
    //   is_debug = false;

    bool is_to_fix = false;
    // check Cells
    if (check_id == 0)
      is_to_fix = is_to_fix_voro_cc(msphere, is_debug);
    else if (check_id == 1) {
      // Note: must call is_to_fix_voro_cc() ahead to load
      // msphere.pcell.cc_bfids in get_CC_boundary_fids()
      // printf("calling is_to_fix_voro_euler ...\n");
      is_to_fix = is_to_fix_voro_euler(convex_cells_host, msphere, is_debug);
    }
    // check Facets
    else if (check_id == 2) {
      // printf("calling is_to_fix_facet_cc ...\n");
      is_to_fix = is_to_fix_facet_cc(all_medial_spheres, msphere, is_debug);
    } else if (check_id == 3) {
      // Note: must call is_to_fix_facet_cc() ahead to load
      // msphere.pcell.facet_cc_boundary_fids in get_facet_CC_boundary_fids()
      // printf("calling is_to_fix_facet_euler...\n");
      is_to_fix = is_to_fix_facet_euler(convex_cells_host, msphere, is_debug);
    }
    // check Edges
    else if (check_id == 4) {
      // printf("calling is_to_fix_edge_cc ...\n");
      is_to_fix = is_to_fix_edge_cc(msphere, is_debug);
    }

    // if (is_to_fix) spheres_to_fix.insert(msphere.id);
    spheres_to_fix_tmp[sphere_id] = is_to_fix;
  };  // run thread

  printf("---------------- checking #CC and Euler ...\n");
  GEO::parallel_for(0, all_medial_spheres.size(),
                    [&](int sphere_id) { run_thread_check(sphere_id, 0); });
  GEO::parallel_for(0, all_medial_spheres.size(),
                    [&](int sphere_id) { run_thread_check(sphere_id, 2); });
  GEO::parallel_for(0, all_medial_spheres.size(),
                    [&](int sphere_id) { run_thread_check(sphere_id, 1); });

  spheres_to_fix.clear();
  for (int i = 0; i < spheres_to_fix_tmp.size(); i++) {
    if (spheres_to_fix_tmp[i]) spheres_to_fix.insert(i);
  }
}

bool fix_topo_is_delete_or_skip(std::vector<MedialSphere> &all_medial_spheres,
                                const int sphere_id, int &num_sphere_change) {
  MedialSphere &msphere = all_medial_spheres[sphere_id];
  if (msphere.is_deleted) return true;
  // -----------------------------------------------------------------------
  // Topo_Status::high_cell_cc too many times
  // -----------------------------------------------------------------------
  if (msphere.itr_topo_fix > 3 &&
      msphere.pcell.topo_status == Topo_Status::high_cell_cc) {
    // still no fix, then delete
    msphere.is_deleted = true;
    num_sphere_change++;
    printf(
        "[Fix] msphere %d topo cannot be fixed, type %d, itr_topo_fix: %d, "
        "delete \n",
        msphere.id, msphere.type, msphere.itr_topo_fix);
    return true;
  }
  // -----------------------------------------------------------------------
  // Topo_Status::low_cell_euler do not fix, just delete,
  // pcell is non-manifold if #euler < 0 and #CC == 1
  // -----------------------------------------------------------------------
  if (msphere.pcell.topo_status == Topo_Status::low_cell_euler) {
    msphere.is_deleted = true;
    num_sphere_change++;
    printf(
        "[Fix] msphere %d topo skip fixing, topo_status %d, itr_topo_fix: %d, "
        "delete \n",
        msphere.id, msphere.pcell.topo_status, msphere.itr_topo_fix);
    return true;
  }

  // -----------------------------------------------------------------------
  // Powercell FaceCC
  // make sure the number of FacetCC is the same for msphere.id and neigh_id
  // [msphere.id, neigh_id] -> { set of cell_ids in one facet CC }
  // -----------------------------------------------------------------------
  // done by function delete_degenerated_medial_spheres()

  // -----------------------------------------------------------------------
  // Powercell EdgeCC
  // make sure the number of EdgeCC is the same for sphere1, sphere2 and sphere3
  // [sphere0, sphere1, sphere2] -> {set of cell_ids in one edge CC}
  // -----------------------------------------------------------------------
  // done by function delete_degenerated_medial_spheres()

  // -----------------------------------------------------------------------
  // Topo_Status::high_facet_cc too many times, then skip fixing
  // This might occurs around the concave lines
  // -----------------------------------------------------------------------
  if (msphere.itr_topo_fix > 4 &&
      msphere.pcell.topo_status == Topo_Status::high_facet_cc) {
    // mark as fixed
    for (const auto &neigh_is_fixed : msphere.pcell.facet_neigh_is_fixed) {
      int neigh_id = neigh_is_fixed.first;
      all_medial_spheres.at(neigh_id).fcc_fixed(msphere.id);
    }
    printf(
        "[Fix] msphere %d topo skip fixing, topo_status %d, itr_topo_fix: %d, "
        "skip\n",
        msphere.id, msphere.pcell.topo_status, msphere.itr_topo_fix);
    return true;
  }

  return false;
}

void fix_topo_facet_cc_euler(
    const int num_itr_global,
    const std::vector<ConvexCellHost> &convex_cells_host,
    const SurfaceMesh &sf_mesh, const TetMesh &tet_mesh, const int sphere_id,
    std::vector<MedialSphere> &all_medial_spheres, bool is_debug) {
  MedialSphere &msphere = all_medial_spheres[sphere_id];
  // only for these two cases
  assert(msphere.pcell.topo_status == Topo_Status::high_facet_cc ||
         msphere.pcell.topo_status == Topo_Status::low_facet_euler);
  const auto &facet_cc_cells = msphere.pcell.facet_cc_cells;
  const auto &facet_cc_surf_v2fids = msphere.pcell.facet_cc_surf_v2fids;
  // if current sphere is a concave sphere,
  // then we will check if newly added sphere close to concave line
  // this check will be large since we have bbox 1000^3
  bool is_merge_to_ce = false;
  if (msphere.pcell.topo_status == Topo_Status::high_facet_cc &&
      !msphere.pcell.ce_covered_lvids.empty()) {
    is_merge_to_ce = true;
  }
  // fix each halfplane faceCC [sphere_id, neigh_id]
  for (const auto &neigh_is_fixed : msphere.pcell.facet_neigh_is_fixed) {
    int neigh_id = neigh_is_fixed.first;
    const auto &msphere_neigh = all_medial_spheres.at(neigh_id);
    // Many skipping conditions
    //
    // Condition1:
    // if neigh_id has fixed the problem
    if (neigh_is_fixed.second == true) {
      printf("[FixFacet] halfplane [%d, %d] has been fixed \n", sphere_id,
             neigh_id);
      continue;
    }
    // Condition2:
    // fix Topo_Status::high_cell_cc first
    if (msphere_neigh.pcell.topo_status == Topo_Status::high_cell_cc) {
      printf("[FixFacet] neigh_id %d has high CC, fix it first \n", neigh_id);
      continue;
    }
    // Condition3:
    // We only fix FacetCC (Topo_Status::high_facet_cc) if both spheres
    // has the same problem, otherwise denegerate happens， delete
    // sphere
    //
    // should be already handled by fix_topo_is_delete_or_skip()
    if (!msphere_neigh.fcc_is_to_fix(sphere_id)) {
      printf("[FixFacet] neigh_id %d not contain faceCC of sphere_id %d\n",
             neigh_id, sphere_id);
      printf("[FixFacet] [%d,%d] has FacetCC [%d,%d] \n", sphere_id, neigh_id,
             msphere.pcell.facet_cc_cells.at(neigh_id).size(),
             msphere_neigh.pcell.facet_cc_cells.at(sphere_id).size());
      assert(false);
    }
    assert(msphere_neigh.fcc_is_to_fix(sphere_id));

    // cell_ids grouped by each facet CC
    const auto &fcc_cells_in_group = facet_cc_cells.at(neigh_id);
    const auto &fcc_v2fids_in_group = facet_cc_surf_v2fids.at(neigh_id);
    // we want to add ${#facet_cc -1} new spheres if #facet_cc is wrong
    // or add ${#facet_cc} spheres if #facet_euler is wrong
    int num_new_spheres = fcc_cells_in_group.size() - 1;
    if (msphere.pcell.topo_status == Topo_Status::low_facet_euler) {
      num_new_spheres = fcc_cells_in_group.size();
    }
    if (is_debug)
      printf(
          "[FixFacet] msphere %d adding new sphere for neigh_id: %d, "
          "num_new_spheres: %d \n",
          sphere_id, neigh_id, num_new_spheres);

    // add new spheres to fix facet_cc or facet_euler
    for (int i = 0; i < num_new_spheres; i++) {
      // select a random point on the halfplane of a random cell
      // that contains the facetCC
      int tmp_rand_id = RANDOM_INT(0, fcc_cells_in_group[i].size());
      int tmp_cell_id = *std::next(fcc_cells_in_group[i].begin(), tmp_rand_id);
      const auto &cc_trans = convex_cells_host.at(tmp_cell_id);
      // select a random point on halfplane face
      bool is_point_found = false;
      Vector3 point(0.f, 0.f, 0.f);
      FOR(lvid, cc_trans.nb_v) {
        cuchar4 v = cc_trans.ver_trans_const(lvid);
        int surf_fid = -1;
        FOR(i, 3) {
          uchar lf = cc_trans.ith_plane_const(lvid, i);
          assert(cc_trans.active_clipping_planes[lf] > 0);
          cint2 hp = cc_trans.clip_id2_const(lf);
          if (hp.y == -1) continue;  // not halfplane
          int cc_neigh_id = hp.x == cc_trans.voro_id ? hp.y : hp.x;
          if (cc_neigh_id != neigh_id) continue;  // not target halfplane
          is_point_found = true;
          point = to_vec(cc_trans.compute_vertex_coordinates(cmake_uchar3(v)));
          // printf(
          //     "[FixFacet] msphere %d found cell %d with lv %d on lf %d,
          //     point (%f,%f,%f)\n", sphere_id, cc_trans.id, lvid, lf,
          //     point[0], point[1], point[2]);
          break;
        }
        if (is_point_found) break;
      }
      assert(is_point_found);
      // get surface v2fid closet to the given point
      v2int v2fid_chosen =
          get_v2fid_min_to_point(sf_mesh, fcc_v2fids_in_group[i], point);
      bool is_good = add_new_sphere_given_v2fid(
          num_itr_global, sf_mesh, tet_mesh, v2fid_chosen, all_medial_spheres,
          is_merge_to_ce, is_debug);
      // Note: to make sure msphere is still valid
      // add_new_sphere_given_v2fid() add new sphere to
      // all_medial_spheres which make msphere access not stable
      auto &msphere = all_medial_spheres.at(sphere_id);
      printf("[FixFacet] msphere %d newly added sphere %d is_good: %d\n",
             sphere_id, all_medial_spheres.back().id, is_good);
      if (is_debug) {
        printf(
            "[FixFacet] msphere %d found v2fid_chosen %d, pin (%f,%f,%f) "
            "sf_mesh facets size: %d, fcc_cells_in_group size: %ld, "
            "fcc_v2fids_in_group size: %ld, random_cell_id: %d \n",
            sphere_id, v2fid_chosen.second, v2fid_chosen.first[0],
            v2fid_chosen.first[1], v2fid_chosen.first[2], sf_mesh.facets.nb(),
            fcc_cells_in_group.size(), fcc_v2fids_in_group.size(), cc_trans.id);
      }
    }  // done for num_new_spheres

    // Note: to make sure msphere is still valid
    auto &msphere = all_medial_spheres.at(sphere_id);
    // to makes sure sphere neigh_id won't repeatly add a new sphere
    all_medial_spheres.at(neigh_id).fcc_fixed(sphere_id);
  }  // for each halfplane faceCC [sphere_id, neigh_id]
}

void fix_topo_cell_cc_euler(
    const int num_itr_global,
    const std::vector<ConvexCellHost> &convex_cells_host,
    const SurfaceMesh &sf_mesh, const TetMesh &tet_mesh, const int sphere_id,
    std::vector<MedialSphere> &all_medial_spheres, bool is_debug) {
  MedialSphere &msphere = all_medial_spheres[sphere_id];
  assert(msphere.pcell.topo_status == Topo_Status::high_cell_cc ||
         msphere.pcell.topo_status == Topo_Status::low_cell_euler);

  std::vector<v2int> v2fid_chosen_vec;
  // if sphere is on sharp edge, then add sphere wisely
  if (msphere.is_on_se() &&
      msphere.pcell.topo_status == Topo_Status::high_cell_cc) {
    // collect surface info
    // CC that contains SE will be skipped
    std::vector<v2int> surf_v2fids;
    for (const auto &one_cc_cells : msphere.pcell.cc_cells) {
      surf_v2fids.clear();
      for (int cell_id : one_cc_cells) {
        // // if one cell in CC is to be skipped, then skip the whole CC
        // if (cells_to_skip.find(cell_id) != cells_to_skip.end()) {
        //   surf_v2fids.clear();
        //   if (is_debug)
        //     printf("[Fix] msphere %d skip cell %d \n", msphere.id,
        //     cell_id);
        //   break;  // check next CC
        // }
        if (msphere.pcell.cell_to_surfv2fid.find(cell_id) ==
            msphere.pcell.cell_to_surfv2fid.end())
          continue;
        const auto &v2fid_pairs = msphere.pcell.cell_to_surfv2fid.at(cell_id);
        surf_v2fids.insert(surf_v2fids.end(), v2fid_pairs.begin(),
                           v2fid_pairs.end());
      }  // for cells in one CC

      if (surf_v2fids.empty()) continue;
      v2int v2fid_chosen =
          get_v2fid_max_to_point(sf_mesh, surf_v2fids, msphere.center);
      v2fid_chosen_vec.push_back(v2fid_chosen);
    }

  } else {
    // Note: some bfid match sf_mesh, but some are just assigned uniquely
    // in function load_num_adjacent_cells_and_ids() and match nothing
    for (const auto &surf_v2fids : msphere.pcell.cc_surf_v2fids) {
      if (surf_v2fids.empty()) continue;
      v2int v2fid_chosen =
          get_v2fid_max_to_point(sf_mesh, surf_v2fids, msphere.center);
      v2fid_chosen_vec.push_back(v2fid_chosen);
      // // add a new sphere for each CC using sphere shrinking
      // if (msphere.pcell.topo_status == Topo_Status::high_cell_cc) {
      //   v2fid_chosen = get_random_v2fid(surf_v2fids);
      // } else if (msphere.pcell.topo_status ==
      // Topo_Status::low_cell_euler) { v2fid_chosen =
      //     get_v2fid_max_to_point(sf_mesh, surf_v2fids, msphere.center);
      // }
    }
  }

  // To add new spheres
  // 1. select suface pins
  int num_new_spheres = msphere.pcell.cc_surf_v2fids.size() - 1;
  // if low euler, might be a hole, then must add at least 1
  if (msphere.pcell.topo_status == Topo_Status::low_cell_euler)
    num_new_spheres = msphere.pcell.cc_surf_v2fids.size();
  if (num_new_spheres > v2fid_chosen_vec.size())
    num_new_spheres = v2fid_chosen_vec.size();
  if (is_debug)
    printf(
        "[FixCell] msphere %d has msphere.pcell.cc_surf_v2fids: %ld, "
        "num_new_spheres: %d\n",
        msphere.id, msphere.pcell.cc_surf_v2fids.size(), num_new_spheres);
  // sort v2fid_chosen_vec based on max distance to sphere center
  std::sort(v2fid_chosen_vec.begin(), v2fid_chosen_vec.end(),
            [&](v2int lhs, v2int rhs) {
              return (lhs.first - msphere.center).length2() >
                     (rhs.first - msphere.center).length2();
            });
  if (is_debug) {
    printf("[FixCell] found v2fid_sorted: [");
    for (const auto &v2fid_sorted : v2fid_chosen_vec) {
      printf("%d, ", v2fid_sorted.second);
    }
    printf("]\n");
  }
  // 2. add new spheres
  for (int i = 0; i < num_new_spheres; i++) {
    v2int &v2fid_chosen = v2fid_chosen_vec[i];
    if (is_debug)
      printf(
          "[FixCell] adding %d/%d, msphere %d found v2fid_chosen %d, sf_mesh "
          "facets size:%d\n",
          i + 1, num_new_spheres, msphere.id, v2fid_chosen.second,
          sf_mesh.facets.nb());
    bool is_good = add_new_sphere_given_v2fid(
        num_itr_global, sf_mesh, tet_mesh, v2fid_chosen, all_medial_spheres,
        true /*is_merge_to_ce*/, is_debug);
    printf("[FixCell] msphere %d newly added sphere %d is_good: %d\n",
           msphere.id, all_medial_spheres.back().id, is_good);
  }
}

/**
 * @brief Adding new spheres to fix either CC or Euler.
 * @param sf_mesh
 * @param spheres_to_fix
 * @param all_medial_spheres
 * @param is_debug
 */
int fix_topo_by_adding_new_sphere(
    const int num_itr_global, const int num_topo_itr,
    const std::vector<ConvexCellHost> &convex_cells_host,
    const SurfaceMesh &sf_mesh, const TetMesh &tet_mesh,
    const std::set<int> &spheres_to_fix,
    std::vector<MedialSphere> &all_medial_spheres, bool is_debug) {
  printf("---------------- Fixing #CC and Euler ...\n");
  printf("[Fix] found %ld spheres to fix \n", spheres_to_fix.size());
  int num_sphere_change = 0;
  int old_num_spheres = all_medial_spheres.size();

  // fix by adding new spheres
  for (const int mid : spheres_to_fix) {
    auto &msphere = all_medial_spheres.at(mid);  // yes, coppy
    if (msphere.is_deleted) continue;            // dont care
    assert(msphere.pcell.topo_status != Topo_Status::unkown &&
           msphere.pcell.topo_status != Topo_Status::ok);
    printf("[Fix] msphere %d has topo_status: %d, itr_topo_fix: %d\n",
           msphere.id, msphere.pcell.topo_status, msphere.itr_topo_fix);

    msphere.itr_topo_fix++;
    if (fix_topo_is_delete_or_skip(all_medial_spheres, mid, num_sphere_change))
      continue;

    // -----------------------------------------------------------------------
    // Topo_Status::high_facet_cc and Topo_Status::low_facet_euler
    // -----------------------------------------------------------------------
    if (msphere.pcell.topo_status == Topo_Status::high_facet_cc ||
        msphere.pcell.topo_status == Topo_Status::low_facet_euler) {
      fix_topo_facet_cc_euler(num_itr_global, convex_cells_host, sf_mesh,
                              tet_mesh, mid, all_medial_spheres, is_debug);
    }
    // -----------------------------------------------------------------------
    // Topo_Status::high_cell_cc and Topo_Status::low_cell_euler
    // we only add num_new_spheres
    // -----------------------------------------------------------------------
    else if (msphere.pcell.topo_status == Topo_Status::high_cell_cc ||
             msphere.pcell.topo_status == Topo_Status::low_cell_euler) {
      fix_topo_cell_cc_euler(num_itr_global, convex_cells_host, sf_mesh,
                             tet_mesh, mid, all_medial_spheres, is_debug);
    } else {
      // should not happen
      printf("[Fix] msphere %d topo_status %d is not handled \n", msphere.id,
             msphere.pcell.topo_status);
      assert(false);
    }
  }  // for spheres_to_fix

  /////////////////////////////////////////////////////////////
  // sanity check
  // remove spheres if still Topo_Status::high_cell_cc
  // and no new sphere can be added
  int updated_total =
      all_medial_spheres.size() - old_num_spheres + num_sphere_change;
  if (updated_total == 0) {
    printf("[Fix] no new spheres added, check high_cell_cc spheres \n");
    for (int i = 0; i < all_medial_spheres.size(); i++) {
      auto &msphere = all_medial_spheres[i];
      if (msphere.is_deleted) continue;
      if (msphere.pcell.topo_status == Topo_Status::high_cell_cc) {
        printf(
            "[Fix] msphere %d cannot add sphere, topo_status %d, itr_topo_fix: "
            "%d, delete \n",
            msphere.id, msphere.pcell.topo_status, msphere.itr_topo_fix);
        msphere.is_deleted = true;
        num_sphere_change++;
        updated_total++;
      }
    }
  }  // if updated_total

  printf("[Fix] added %d->%d new spheres, changed spheres: %d\n",
         old_num_spheres, all_medial_spheres.size(), num_sphere_change);
  return updated_total;
}