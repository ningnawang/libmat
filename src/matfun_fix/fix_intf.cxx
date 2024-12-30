#include "fix_intf.h"

#include "common_geogram.h"
#include "fix_topo.h"
// helper functions
bool add_new_T2_sphere_from_common(
    const int num_itr_global, const SurfaceMesh& sf_mesh,
    const TetMesh& tet_mesh, MedialEdge& medge_cand,
    std::vector<MedialSphere>& all_medial_spheres, bool is_debug) {
  const MedialSphere& msphere1 = all_medial_spheres.at(medge_cand.vertices_[0]);
  const MedialSphere& msphere2 = all_medial_spheres.at(medge_cand.vertices_[1]);

  if (medge_cand.common_tan_pls.empty()) {
    if (is_debug) {
      printf("[FIX_INTF_T2] medge (%d,%d) has 0 common_tan_pls\n", msphere1.id,
             msphere2.id);
    }
    return false;
  }

  bool is_good = false;
  for (const auto& tan_pl : medge_cand.common_tan_pls) {
    v2int v2fid_chosen;
    v2fid_chosen.first = tan_pl.tan_point;
    v2fid_chosen.second = tan_pl.fid;
    is_good = insert_new_sphere_given_v2fid(
        num_itr_global, sf_mesh, tet_mesh, v2fid_chosen, all_medial_spheres,
        true /*is_merge_to_ce*/, false /*is_debug*/);
    if (is_good) return true;
  }
  return is_good;
}

bool check_and_fix_if_both_on_intf(
    const int num_itr_global, const SurfaceMesh& sf_mesh,
    const TetMesh& tet_mesh, MedialEdge& medge_cand,
    std::vector<MedialSphere>& all_medial_spheres, bool is_debug) {
  const MedialSphere& msphere1 = all_medial_spheres.at(medge_cand.vertices_[0]);
  const MedialSphere& msphere2 = all_medial_spheres.at(medge_cand.vertices_[1]);
  assert(msphere1.is_on_intf() && msphere2.is_on_intf());
  // if (medge_cand.v1_non_common_tan_pls.size() == 0 &&
  //     medge_cand.v2_non_common_tan_pls.size() == 0 &&
  //     medge_cand.common_tan_pls.size() > 2) {
  if (medge_cand.common_tan_pls.size() >= 2) {
    if (is_debug)
      printf(
          "[FIX_INTF_BOTH] medge (%d,%d) both on intf and common_tan_pls %zu, "
          "skip checking\n",
          msphere1.id, msphere2.id, medge_cand.common_tan_pls.size());
    return true;
  }

  if (is_debug)
    printf("[FIX_INTF_BOTH] medge (%d,%d) needs to add T2 spheres\n",
           msphere1.id, msphere2.id);
  return add_new_T2_sphere_from_common(num_itr_global, sf_mesh, tet_mesh,
                                       medge_cand, all_medial_spheres,
                                       false /*is_debug*/);
}

bool check_and_fix_intf_by_adding_new_spheres(
    const int num_itr_global, const SurfaceMesh& sf_mesh,
    const TetMesh& tet_mesh, MedialEdge& medge_cand,
    std::vector<MedialSphere>& all_medial_spheres, bool is_debug) {
  MedialSphere* msphere1 = &all_medial_spheres.at(medge_cand.vertices_[0]);
  MedialSphere* msphere2 = &all_medial_spheres.at(medge_cand.vertices_[1]);
  medge_cand.is_vs_swapped = false;
  // make sure msphere1 is on internal/external feature if any
  if (msphere2->is_on_intf() && !msphere1->is_on_intf()) {
    msphere1 = &all_medial_spheres.at(medge_cand.vertices_[1]);
    msphere2 = &all_medial_spheres.at(medge_cand.vertices_[0]);
    medge_cand.is_vs_swapped = true;
  } else if (msphere2->is_on_ce() && !msphere1->is_on_ce()) {
    // put T_2^c in front of T_1_2
    msphere1 = &all_medial_spheres.at(medge_cand.vertices_[1]);
    msphere2 = &all_medial_spheres.at(medge_cand.vertices_[0]);
    medge_cand.is_vs_swapped = true;
  } else if (msphere2->tan_cc_lines.size() > msphere1->tan_cc_lines.size()) {
    msphere1 = &all_medial_spheres.at(medge_cand.vertices_[1]);
    msphere2 = &all_medial_spheres.at(medge_cand.vertices_[0]);
    medge_cand.is_vs_swapped = true;
  }

  // if (msphere1->id == 125 || msphere2->id == 125)
  //   is_debug = true;
  // else
  // is_debug = false;

  if (is_debug)
    printf("[FIX_INTF] checking internal sphere for medge (%d,%d)\n",
           msphere1->id, msphere2->id);

  ////////////////////////////
  // Pre-step: checking common/diff tangent info
  std::vector<TangentPlane>* common_tan_pls = &medge_cand.common_tan_pls;
  std::vector<TangentPlane>* A_non_common_tan_pls =
      &medge_cand.v1_non_common_tan_pls;
  std::vector<TangentPlane>* B_non_common_tan_pls =
      &medge_cand.v2_non_common_tan_pls;
  if (medge_cand.is_vs_swapped) {
    A_non_common_tan_pls = &medge_cand.v2_non_common_tan_pls;
    B_non_common_tan_pls = &medge_cand.v1_non_common_tan_pls;
  }
  // assert(!common_tan_pls->empty() || !A_non_common_tan_pls->empty() ||
  //        !B_non_common_tan_pls->empty());
  if (common_tan_pls->empty() && A_non_common_tan_pls->empty() &&
      B_non_common_tan_pls->empty())
    return true;

  // no need to add a new sphere case 1
  if (B_non_common_tan_pls->empty() || A_non_common_tan_pls->empty()) {
    if (is_debug)
      printf(
          "[CREATE_more] medial edge (%d,%d) no need to add new sphere, "
          "empty\n",
          msphere1->id, msphere2->id);
    return true;
  }

  // no need to add a new sphere case 2
  int total_tan = common_tan_pls->size() + A_non_common_tan_pls->size() +
                  B_non_common_tan_pls->size();
  if (total_tan < 3) {
    if (is_debug)
      printf(
          "[CREATE_more] medial edge (%d,%d) no need to add new sphere, , "
          "common_tan_pls size: %zu, A_non_common_tan_pls: %zu, "
          "B_non_common_tan_pls: %zu \n",
          msphere1->id, msphere2->id, common_tan_pls->size(),
          A_non_common_tan_pls->size(), B_non_common_tan_pls->size());
    return true;
  }

  // create new sphere
  bool is_good = false;
  Vector3 mid_center = (msphere1->center + msphere2->center) / 2;
  MedialSphere new_sphere(all_medial_spheres.size(), mid_center, INIT_RADIUS,
                          SphereType::T_3_MORE, num_itr_global);
  if (is_debug)
    printf("[FIX_INTF] trying to add new_sphere %d\n", new_sphere.id);

  ////////////////////////////
  // Step1: 1. both not on internal feature
  //        2. both on external feature, but not the same SL (checked
  //        ahead) aggregate tangent info then try to add T_N/T_N_c sphere
  if (!msphere1->is_on_intf() && !msphere2->is_on_intf() ||
      msphere1->is_on_extf() && msphere2->is_on_extf()) {
    if (is_debug) printf("[FIX_INTF] step 1: aggregate all tangent info \n");

    // add tangent planes, including adjacent planes of tan_cc_line
    if (is_debug) {
      printf("common_tan_pls: ----------------------------\n");
      for (const auto& common : *common_tan_pls) {
        new_sphere.tan_planes.push_back(common);
        common.print_info();
      }
    }
    for (const auto& a_diff : *A_non_common_tan_pls)
      new_sphere.tan_planes.push_back(a_diff);
    for (const auto& b_diff : *B_non_common_tan_pls)
      new_sphere.tan_planes.push_back(b_diff);
    // add tangent concave lines, since iterate_sphere() may not detect
    for (const auto& a_tan_cc_line : msphere1->tan_cc_lines)
      new_sphere.new_cc_line_no_dup(a_tan_cc_line);
    for (const auto& b_tan_cc_line : msphere2->tan_cc_lines)
      new_sphere.new_cc_line_no_dup(b_tan_cc_line);
    // purge using tan_cc_lines
    new_sphere.purge_and_delete_tan_planes();

    // try to add T_N/T_N_c sphere after aggregation
    is_good = iterate_sphere(sf_mesh, sf_mesh.aabb_wrapper,
                             sf_mesh.fe_sf_fs_pairs, tet_mesh.feature_edges,
                             new_sphere, is_debug /*is_debug*/);

    if (is_debug)
      printf("[FIX_INTF] step 1 is_good: %d, with tan_planes: %zu\n", is_good,
             new_sphere.tan_planes.size());
  }

  ////////////////////////////
  // Step 2: otherwise
  //         1. step 1 failed
  //         2. both on external feature
  //         2. only one sphere not on internal feature
  //         try add a seam sphere (eg. T3) than junction sphere (eg. T4)
  // Note:  if msphere2 is not on internal feature, we need to
  //        check if two spheres are on the same sheet, if so, then not add
  //        (check if msphere2 is already covered by msphere1)
  if (!is_good &&
      (!is_good || (msphere1->is_on_intf() || msphere1->is_on_ce()) &&
                       !msphere2->is_on_intf())) {
    if (is_debug)
      printf("[FIX_INTF] step 2: aggregate partial tangent info \n");

    // will add only one A_non-common
    for (const auto& A_non_common : *A_non_common_tan_pls) {
      // reset
      new_sphere.tan_cc_lines.clear();
      new_sphere.tan_planes.clear();
      // add common
      for (const auto& common : *common_tan_pls) {
        new_sphere.tan_planes.push_back(common);
      }
      // add B_non_common
      for (const auto& B_non_common : *B_non_common_tan_pls) {
        new_sphere.tan_planes.push_back(B_non_common);
      }
      // add only one A_non-common
      new_sphere.tan_planes.push_back(A_non_common);
      if (is_debug)
        printf("[FIX_INTF] try to add A_non_common %d \n", A_non_common.fid);
      // add one A_non-common
      if (iterate_sphere(sf_mesh, sf_mesh.aabb_wrapper, sf_mesh.fe_sf_fs_pairs,
                         tet_mesh.feature_edges, new_sphere,
                         is_debug /*is_debug*/)) {
        if (is_debug)
          printf("[FIX_INTF] step 2: add T_lower sphere %d for medge (%d,%d)\n",
                 new_sphere.id, msphere1->id, msphere2->id);
        is_good = true;
        break;
      } else {  // else: remove new_tan_pl, for next for loop
        if (is_debug) printf("[FIX_INTF] step 2: re-try next A_non_common\n");
      }
    }  // for A_non_common_tan_pls
    if (is_debug) printf("[FIX_INTF] step 2 is_good: %d\n", is_good);
  }

  if (!is_good) {
    if (is_debug)
      printf(
          "[FIX_INTF] failed to add internal feature sphere for medge "
          "(%d,%d)\n",
          msphere1->id, msphere2->id);
    return false;
  }

  // iterate_sphere may merge tangent planes as well
  if (new_sphere.tan_planes.size() + new_sphere.tan_cc_lines.size() < 3) {
    if (is_debug)
      printf(
          "[FIX_INTF] medial edge (%d,%d) no need to add sphere after "
          "iteration, tan_plane size: %ld, tan_cc_line size: %ld\n",
          msphere1->id, msphere2->id, new_sphere.tan_planes.size(),
          new_sphere.tan_cc_lines.size());
    return true;
  }

  // internal feature condition 3
  // not add sphere if newly added sphere is already covered by msphere1
  std::vector<TangentPlane> tmp_common_tan_pls, tmp_A_non_common_tan_pls,
      tmp_B_non_common_tan_pls;
  A_B_spheres_common_diff_tangent_info_from_surface(
      sf_mesh, *msphere1, new_sphere, tmp_common_tan_pls,
      tmp_A_non_common_tan_pls, tmp_B_non_common_tan_pls);
  if (tmp_B_non_common_tan_pls.empty()) {
    if (is_debug)
      printf(
          "[FIX_INTF] medial edge (%d,%d) no need to add new sphere %d, "
          " same tangent info as old spheres \n",
          msphere1->id, msphere2->id, new_sphere.id);
    return true;
  }

  // //////////
  // // step 2: if failed, then check if two spheres belongs to the same
  // sheet,
  // //          if not then add a T_2 sphere
  // if (!is_good) {
  //   // return false;
  //   is_good = check_mat_sheet_and_fix(sf_mesh, tet_mesh, msphere1,
  //   msphere2,
  //                                     new_sphere, is_debug);
  // if (!is_good) return false;
  // }

  if (!add_new_sphere_validate(all_medial_spheres, new_sphere)) return false;
  return true;
}

bool random_check_edges_and_fix_internal_feature(
    const int num_itr_global, const SurfaceMesh& sf_mesh,
    const TetMesh& tet_mesh, const MedialMesh& mmesh,
    std::vector<MedialSphere>& all_medial_spheres,
    std::set<aint2>& checked_mmesh_edges, bool is_debug) {
  // is_debug = true;

  // randomly select k unchecked edges to check
  // int random_k = mmesh.edges.size() / 20;
  int random_k = mmesh.edges.size() / 100;
  std::vector<MedialEdge> unchecked_edges;
  int num_old_spheres = all_medial_spheres.size();
  int num_new_spheres = 0;
  std::vector<MedialEdge> random_medges;
  while ((num_new_spheres = all_medial_spheres.size() - num_old_spheres) <
             random_k &&
         random_k > 0) {
    // update unchecked_edges
    unchecked_edges.clear();
    for (const auto& medge : mmesh.edges) {
      aint2 mvs = medge.vertices_;
      std::sort(mvs.begin(), mvs.end());
      if (checked_mmesh_edges.find(mvs) != checked_mmesh_edges.end()) continue;
      unchecked_edges.push_back(medge);
    }
    if (is_debug)
      printf("[FIX_INTF] unchecked_edges: %zu, checked_mmesh_edges: %zu \n",
             unchecked_edges.size(), checked_mmesh_edges.size());

    // update random_edges
    num_old_spheres = all_medial_spheres.size();
    random_k = random_k - num_new_spheres;
    if (random_k > unchecked_edges.size()) random_k = unchecked_edges.size();
    random_medges.clear();
    for (int i = 0; i < random_k; i++) {
      MedialEdge& rand_e =
          unchecked_edges.at(RANDOM_INT(0, unchecked_edges.size()));
      while (checked_mmesh_edges.find(rand_e.vertices_) !=
             checked_mmesh_edges.end()) {
        rand_e = unchecked_edges.at(RANDOM_INT(0, unchecked_edges.size()));
      }
      random_medges.push_back(rand_e);
    }
    if (is_debug)
      printf("[FIX_INTF] select %zu random_medges to fix, random_k: %d \n",
             random_medges.size(), random_k);
    assert(random_medges.size() == random_k);

    // loop random_edges
    for (auto& medge_cand : random_medges) {
      if (is_debug)
        printf("[FIX_INTF] checking random edge (%d,%d) \n",
               medge_cand.vertices_[0], medge_cand.vertices_[1]);
      MedialSphere& msphere1 = all_medial_spheres.at(medge_cand.vertices_[0]);
      MedialSphere& msphere2 = all_medial_spheres.at(medge_cand.vertices_[1]);

      if (msphere1.id == 1319 && msphere2.id == 1319)
        is_debug = true;
      else
        is_debug = false;

      // skip checked edges
      if (checked_mmesh_edges.find(medge_cand.vertices_) !=
          checked_mmesh_edges.end()) {
        if (is_debug)
          printf("[FIX_INTF] medial edge (%d,%d) has checked, pass! \n",
                 medge_cand.vertices_[0], medge_cand.vertices_[1]);
        continue;
      }
      checked_mmesh_edges.insert(medge_cand.vertices_);

      bool is_good_or_add = false;
      // not handling corners so far
      if (msphere1.is_on_corner() || msphere2.is_on_corner()) continue;

      // if both on intf, then skip
      if (msphere1.is_on_intf() && msphere2.is_on_intf()) {
        check_and_fix_if_both_on_intf(num_itr_global, sf_mesh, tet_mesh,
                                      medge_cand, all_medial_spheres,
                                      false /*is_debug*/);
        continue;
      }

      // case 1: both spheres on same sharp line, not corners
      if (msphere1.is_on_se() && msphere2.is_on_se() &&
          is_two_mspheres_on_same_se(msphere1, msphere2)) {
        continue;
      }

      // case 2: general case
      is_good_or_add = check_and_fix_intf_by_adding_new_spheres(
          num_itr_global, sf_mesh, tet_mesh, medge_cand, all_medial_spheres,
          is_debug);

      // // [no] this would add too many spheres
      // // case 3: if failed, then add a T_2 sphere
      // if (!is_good_or_add && (msphere1.is_on_intf() ||
      // msphere2.is_on_intf())) {
      //   if (is_debug)
      //     printf(
      //         "[FIX_INTF] medial edge (%d,%d) failed to add intf sphere, add
      //         " "T2 sphere instead\n ", medge_cand.vertices_[0],
      //         medge_cand.vertices_[1]);
      //   is_good_or_add = add_new_T2_sphere_from_common(
      //       num_itr_global, sf_mesh, tet_mesh, medge_cand,
      //       all_medial_spheres, false /*is_debug*/);
      // }

      if (is_debug)
        printf("[FIX_INTF] is_good_or_add: %d for edge (%d,%d) \n",
               is_good_or_add, msphere1.id, msphere2.id);
    }  // for random_edges
  }

  // printf("[FIX_INTF] all %d medial edges good, INTF is good!\n",
  //        mmesh.edges.size());
  return false;
}

// [no use]
int delete_T2_if_TN(const int num_itr_global, const SurfaceMesh& sf_mesh,
                    const TetMesh& tet_mesh,
                    std::vector<MedialSphere>& all_medial_spheres,
                    bool is_debug) {
  std::vector<MedialSphere> new_medial_spheres(all_medial_spheres.size());

  // for parallization
  auto run_thread = [&](int mid) {
    auto& msphere = all_medial_spheres.at(mid);
    if (!msphere.is_on_sheet()) {
      new_medial_spheres[mid] = msphere;
      return;
    }
    assert(!msphere.covered_sf_fids_in_group.empty());
    if (msphere.get_sf_covered_group_size() < 3) {
      new_medial_spheres[mid] = msphere;
      return;
    }

    bool is_to_update = false;
    // update T2 -> TN
    if (msphere.type == SphereType::T_2 &&
        msphere.get_sf_covered_group_size() > 2) {
      is_to_update = true;
      if (msphere.pcell.ce_covered_lvids.empty()) {
        // only update when not T_2_c
        msphere.type = SphereType::T_3_MORE;
      }
    }
    // update T_2_c -> T_N_c
    if (msphere.type == SphereType::T_2_c &&
        msphere.get_sf_covered_group_size() > 2) {
      is_to_update = true;
      msphere.type = SphereType::T_N_c;
    }

    // just update type, no real update
    new_medial_spheres[mid] = msphere;
    if (is_to_update) {
      new_medial_spheres[mid].is_deleted = true;
      printf("[REPLACE] delete sphere %d\n", msphere.id);
    }
  };

  // for (int i = 0; i < all_medial_spheres.size(); i++) {
  // run_thread(i);
  // }
  GEO::parallel_for(0, all_medial_spheres.size(),
                    [&](int mid) { run_thread(mid); });

  // replace all_medial_spheres by new_medial_spheres
  assert(all_medial_spheres.size() == new_medial_spheres.size());
  all_medial_spheres.clear();
  all_medial_spheres = new_medial_spheres;

  int num_replace = 0, num_deleted = 0;
  for (const auto& msphere : all_medial_spheres) {
    if (msphere.itr_cnt == num_itr_global) num_replace++;
    if (msphere.is_deleted) num_deleted++;
  }

  printf(
      "[REPLACE] replaced %d spheres, deleted %d to internal feature spheres\n",
      num_replace, num_deleted);
  return num_deleted;
}