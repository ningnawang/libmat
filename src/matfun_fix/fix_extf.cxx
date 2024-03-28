#include "fix_extf.h"

#include "sharp_feature_detection.h"

// we add new sphere centers on one sharp edges group
// (belongs to one sharp line)
// using given fe_abs_len to roughly control the density
// (adding 1 feature sphere at a time is too time consuming for iterations)
void add_new_centers_given_one_se_group(
    const std::vector<aint4>& se_one_group,
    const std::map<int, Vector3>& vs_pos_map, const float fe_abs_len,
    bool is_add_one_empty, std::vector<v2int2>& new_centers_fl, bool is_debug) {
  // NOTE: do not clear new_centers_fl!!!!
  float len_tmp = fe_abs_len;
  int cur_fl_id = se_one_group[0][2];  // index 2
  for (const auto& one_se : se_one_group) {
    assert(vs_pos_map.find(one_se[0]) != vs_pos_map.end());
    assert(vs_pos_map.find(one_se[1]) != vs_pos_map.end());
    assert(one_se[2] == cur_fl_id);  // belongs to the same fl
    // save the FeatureEdge::id as well
    aint2 se_info = {{cur_fl_id, one_se[3]}};
    Vector3 e0 = vs_pos_map.at(one_se[0]);
    Vector3 e1 = vs_pos_map.at(one_se[1]);
    double edge_len = (e1 - e0).length();
    Vector3 edge_dir = GEO::normalize(e1 - e0);
    if (is_debug)
      printf("[FIX_EXTF] edge (%d,%d) has edge_len %f \n", one_se[0], one_se[1],
             edge_len);
    if (edge_len > fe_abs_len) {
      int step_max = std::floor(edge_len / fe_abs_len);
      if (step_max == 1) {
        Vector3 newc = 1. / 2. * (e0 + e1);
        new_centers_fl.push_back(std::make_pair(newc, se_info));
        if (is_debug)
          printf(
              "[FIX_EXTF] edge (%d,%d) cur_fl_id %d, create new center with "
              "middle (%f,%f,%f) \n",
              one_se[0], one_se[1], cur_fl_id, newc[0], newc[1], newc[2]);
        continue;
      }
      for (int step = 1; step < step_max; step++) {
        Vector3 newc = e0 + step * fe_abs_len * edge_dir;
        new_centers_fl.push_back(std::make_pair(newc, se_info));
        if (is_debug)
          printf(
              "[FIX_EXTF] edge (%d,%d) cur_fl_id %d, create new center with "
              "step %d, (%f,%f,%f) \n",
              one_se[0], one_se[1], cur_fl_id, step, newc[0], newc[1], newc[2]);
      }
    } else if (edge_len < len_tmp) {
      len_tmp = len_tmp - edge_len;
    } else {
      // Vector3 newc = e0 + len_tmp * edge_dir;
      // use middle to avoid reverse edge_dir
      Vector3 newc = 1. / 2. * (e0 + e1);
      new_centers_fl.push_back(std::make_pair(newc, se_info));
      len_tmp = fe_abs_len;
      if (is_debug)
        printf(
            "[FIX_EXTF] edge (%d,%d) cur_fl_id %d, create new center with newc "
            "(%f.%f,%f)\n",
            one_se[0], one_se[1], cur_fl_id, newc[0], newc[1], newc[2]);
    }

    // add random middle if empty
    // this is for adding at least one zero-radius sphere than nothing
    if (is_add_one_empty && new_centers_fl.empty()) {
      aint4 mid_e = se_one_group[0];
      Vector3 newc =
          1. / 2. * (vs_pos_map.at(mid_e[0]) + vs_pos_map.at(mid_e[1]));
      new_centers_fl.push_back(std::make_pair(newc, se_info));
      if (is_debug)
        printf(
            "[FIX_EXTF] edge (%d,%d) cur_fl_id %d, create random newc "
            "(%f.%f,%f)\n",
            one_se[0], one_se[1], cur_fl_id, newc[0], newc[1], newc[2]);
    }
  }  // for se_one_group
}

void fix_corner_sphere(float fe_abs_len, MedialSphere& msphere,
                       const std::vector<std::vector<aint4>>& se_tvs_grouped,
                       const std::map<int, Vector3>& tvs_pos_map,
                       std::vector<v2int2>& new_centers_fl, bool is_debug) {
  if (is_debug)
    printf("[FIX_EXTF corner] calling fix_corner for msphere %d\n", msphere.id);
  assert(msphere.is_on_corner());
  assert(!msphere.corner_fls.empty());  // call init_corner_spheres() ahead
  new_centers_fl.clear();
  int num_covered_se_group = se_tvs_grouped.size();
  std::set<int> visited_se_lines;
  for (int i = 0; i < num_covered_se_group; i++) {
    const auto& se_one_group = se_tvs_grouped[i];
    int cur_se_line_id = se_one_group[0][2];  // inherited from tet_mesh
    // skip if current feature line is handled
    if (visited_se_lines.find(cur_se_line_id) != visited_se_lines.end())
      continue;
    visited_se_lines.insert(cur_se_line_id);

    // if SE not adjacent to corner
    // then add at least 1 new center
    bool is_add_one_empty = true;
    // if SE is adjacent to corner
    // then only add new center if too long, add 0 otherwise
    // (set is_add_one_empty = false)
    if (msphere.corner_fls.find(cur_se_line_id) != msphere.corner_fls.end())
      is_add_one_empty = false;

    // add new centers
    add_new_centers_given_one_se_group(se_one_group, tvs_pos_map, fe_abs_len,
                                       is_add_one_empty, new_centers_fl,
                                       is_debug /*is_debug*/);
  }
}

// only add new
void fix_se_sphere(const float fe_abs_len, MedialSphere& msphere,
                   const std::vector<MedialSphere>& all_medial_spheres,
                   const std::vector<std::vector<aint4>>& se_tvs_grouped,
                   const std::map<int, Vector3>& tvs_pos_map,
                   const std::map<int, std::set<int>>& fl2corner_sphere,
                   std::vector<v2int2>& new_centers_fl, bool is_debug) {
  if (is_debug)
    printf("[FIX_EXTF SE] calling fix_se_sphere for msphere %d\n", msphere.id);
  assert(msphere.is_on_se());
  new_centers_fl.clear();
  int num_covered_se_group = se_tvs_grouped.size();
  // stays there, to avoid new insertions
  if (num_covered_se_group == 1) return;
  // if other covered_se_line contains the same corner
  // and msphere is close to that corner
  // then we stop inserting
  std::set<int> prev_corner_ids, cur_corner_ids;
  if (fl2corner_sphere.find(msphere.se_line_id) != fl2corner_sphere.end())
    prev_corner_ids = fl2corner_sphere.at(msphere.se_line_id);
  std::set<int> visited_se_lines;
  for (int i = 0; i < num_covered_se_group; i++) {
    const auto& se_one_group = se_tvs_grouped[i];
    int cur_se_line_id = se_one_group[0][2];  // inherited from tet_mesh
    // skip if current feature line is handled
    if (visited_se_lines.find(cur_se_line_id) != visited_se_lines.end())
      continue;
    visited_se_lines.insert(cur_se_line_id);

    // if SE not adjacent to sphere
    // then add at least 1 new center
    bool is_add_one_empty = true;
    // if SE is adjacent to sphere
    // then only add new center if too long, add 0 otherwise
    // (set is_add_one_empty = false)
    if (cur_se_line_id == msphere.se_line_id) is_add_one_empty = false;

    // if multiple SE joint on a corner
    // and msphere is close to this corner, then stop adding new SE spheres
    bool is_add_new_centers = true;
    cur_corner_ids.clear();
    if (cur_se_line_id != msphere.se_line_id && !prev_corner_ids.empty()) {
      if (fl2corner_sphere.find(cur_se_line_id) != fl2corner_sphere.end())
        cur_corner_ids = fl2corner_sphere.at(cur_se_line_id);
      std::vector<int> corners_common;
      set_intersection(prev_corner_ids, cur_corner_ids, corners_common);
      for (const auto& cid : corners_common) {
        double dist = get_distance_between_two_vectors(
            all_medial_spheres.at(cid).center, msphere.center);
        if (dist < SCALAR_1) {
          if (is_debug)
            printf(
                "[FIX_EXTF SE] msphere %d close to corner %d with dist %f, not "
                "add for cur_se_line_id %d, msphere.se_line_id %d \n",
                msphere.id, cid, dist, cur_se_line_id, msphere.se_line_id);
          is_add_new_centers = false;
          break;
        }
      }  // for corners_common
    }

    if (is_add_new_centers) {
      // add new centers
      add_new_centers_given_one_se_group(se_one_group, tvs_pos_map, fe_abs_len,
                                         is_add_one_empty, new_centers_fl,
                                         is_debug /*is_debug*/);
    }
  }
}

// may 1. update current, 2. add new
// return number of updated spheres
int fix_non_extf_sphere(const int num_itr_global, float fe_abs_len,
                        MedialSphere& msphere,
                        const std::vector<std::vector<aint4>>& se_tvs_grouped,
                        const std::map<int, Vector3>& tvs_pos_map,
                        const std::vector<MedialSphere>& all_medial_spheres,
                        const std::map<int, std::set<int>>& fl2corner_sphere,
                        std::vector<v2int2>& new_centers_fl, bool is_debug) {
  if (is_debug)
    printf("[FIX_EXTF NonExtf] calling fix_non_extf_sphere for msphere %d\n",
           msphere.id);
  assert(!msphere.is_on_extf());
  new_centers_fl.clear();
  int num_updated = 0;
  int num_covered_se_group = se_tvs_grouped.size();
  int se_group_start_id = 0;
  ////////////////
  // step 1: if sphere's radius is too small
  //         push to middle of sharp line of se_tvs_grouped[0]
  if (msphere.radius < SCALAR_SE_MERGE_RADIUS) {
    const auto& se_one_group = se_tvs_grouped[se_group_start_id++];
    int cur_se_line_id = se_one_group[0][2];
    int cur_se_id = se_one_group[0][3];
    float radius = SCALAR_FEATURE_RADIUS;
    int mid_e_idx = std::floor(se_one_group.size() / 2);
    aint4 mid_e = *std::next(se_one_group.begin(), mid_e_idx);
    Vector3 mid_center =
        1. / 2. * (tvs_pos_map.at(mid_e[0]) + tvs_pos_map.at(mid_e[1]));

    // if close to its corners, then delete than update
    if (fl2corner_sphere.find(cur_se_line_id) != fl2corner_sphere.end()) {
      for (const int corner_id : fl2corner_sphere.at(cur_se_line_id)) {
        const auto corner_sphere = all_medial_spheres.at(corner_id);
        double dist =
            get_distance_between_two_vectors(corner_sphere.center, mid_center);
        if (dist < SCALAR_1) {
          printf(
              "[FIX_EXTF NonExtf] non-feature msphere %d samll radius, too "
              "close to corner %d, delete\n ",
              msphere.id, corner_sphere.id);
          msphere.is_deleted = true;
        }
      }
    }

    // shift sphere on current se group
    // this would update msphere.type as well
    msphere.center = mid_center;
    msphere.radius = radius;
    msphere.type = SphereType::T_1_2;
    msphere.se_line_id = cur_se_line_id;
    msphere.se_edge_id = cur_se_id;
    msphere.itr_cnt = num_itr_global;  // update
    apply_perturb(msphere.center);
    num_updated++;
    if (is_debug && !msphere.is_deleted)
      printf(
          "[FIX_EXTF NonExtf] non-feature msphere %d samll radius, pusehd to "
          "sharp edge: (%f,%f,%f)\n",
          msphere.id, mid_center[0], mid_center[1], mid_center[2]);
  }

  ////////////////
  // step 2: if sphere contains multiple sharp lines
  //         then add at least 1 new center
  // insert sphere for each grouped sharp edges (per sharp line)
  //
  // Note: some tets may share the same feature edge,
  //       we handle one feature line only one time.
  std::set<int> visited_se_lines;
  for (int i = se_group_start_id; i < num_covered_se_group; i++) {
    const auto& se_one_group = se_tvs_grouped[i];
    int cur_se_line_id = se_one_group[0][2];  // inherited from tet_mesh
    // skip if current feature line is handled
    if (visited_se_lines.find(cur_se_line_id) != visited_se_lines.end())
      continue;
    visited_se_lines.insert(cur_se_line_id);
    if (is_debug) {
      printf(
          "[FIX_EXTF NonExtf] non-feature msphere %d, type: %d, "
          "msphere.se_line_id: %d, has %d-th se_one_group size %ld of "
          "cur_se_line_id %d: ",
          msphere.id, msphere.type, msphere.se_line_id, i, se_one_group.size(),
          cur_se_line_id);
      for (auto& se : se_one_group) {
        printf("(%d,%d,%d), ", se[0], se[1], se[2]);
      }
      printf("\n");
    }

    // add at least 1 new center (is_add_one_empty = true)
    add_new_centers_given_one_se_group(se_one_group, tvs_pos_map, fe_abs_len,
                                       true /*is_add_one_empty*/,
                                       new_centers_fl, is_debug /*is_debug*/);
    if (is_debug) {
      printf(
          "[FIX_EXTF NonExtf] msphere %d, fe_abs_len: %f, create new_centers: "
          "%ld\n",
          msphere.id, fe_abs_len, new_centers_fl.size());
    }
  }  // for se_tvs_grouped

  return num_updated;
}

int check_and_fix_external_feature(
    const int num_itr_global, const Parameter& param,
    const std::vector<float>& tet_vertices,
    const std::vector<ConvexCellHost>& convex_cells_host,
    const GEO::Mesh& sf_mesh, const std::map<aint4, int>& tet_vs_lfs2tvs_map,
    const std::map<int, std::set<int>>& fl2corner_sphere,
    std::vector<MedialSphere>& all_medial_spheres, bool is_debug) {
  auto add_new_se_sphere =
      [](const int num_itr_global, const int start_row, const Vector3& center,
         const double radius, const aint2 cur_se_info,
         std::vector<MedialSphere>& new_medial_spheres, bool is_debug) {
        uint mid = start_row + new_medial_spheres.size();
        MedialSphere new_msphere(mid, center, radius, SphereType::T_1_2 /*SE*/,
                                 num_itr_global);
        new_msphere.se_line_id = cur_se_info[0];
        new_msphere.se_edge_id = cur_se_info[1];
        apply_perturb(new_msphere.center);
        // TODO: add tangent planes
        new_medial_spheres.push_back(new_msphere);
        if (is_debug)
          printf(
              "[FIX_EXTF] add new new_msphere %d on sharp line %d, sharp edge "
              "%d\n",
              new_msphere.id, new_msphere.se_line_id, new_msphere.se_edge_id);
      };

  float fe_abs_len = param.bbox_diag_l * param.fe_rel_len;
  float radius = SCALAR_FEATURE_RADIUS;
  int num_changed = 0 /*update+add*/, num_added = 0;
  std::vector<MedialSphere> new_medial_spheres;
  for (auto& msphere : all_medial_spheres) {
    // skip if not cover any SE
    const auto& se_covered_lvids = msphere.pcell.se_covered_lvids;
    if (se_covered_lvids.empty()) continue;

    // convert sharp edge representation from
    // convexcell [cell_id, lvid1, lvid2, se_line_id, se_id]
    // to tet [tvid1, tvid2, se_line_id, se_id]
    // 1. if cell_lvid is from tet, then keep the tvid
    // 2. if cell_lvid is newly cutted, then assign negative id (vs_new_ids)
    //
    // for multiple convexcells, multple id may share the same vertex
    // but who cares
    //
    // Note: some tets may share the same feature edge!! (should be rare for
    // sharp edges)
    std::set<aint4> se_tvs_all;
    std::map<int, Vector3> tvs_pos_map;
    convert_fes_from_ccell_to_tet(msphere, tet_vertices, convex_cells_host,
                                  tet_vs_lfs2tvs_map, se_covered_lvids,
                                  se_tvs_all, tvs_pos_map, false /*is_debug*/);

    // group se_covered_lvids into groups (based on CC)
    std::vector<std::vector<aint4>> se_tvs_grouped;
    get_grouped_fes_in_order(se_tvs_all, se_tvs_grouped);

    std::vector<v2int2> new_centers_fl;  // <center, se_fl_id>
    if (msphere.is_on_corner()) {
      fix_corner_sphere(fe_abs_len, msphere, se_tvs_grouped, tvs_pos_map,
                        new_centers_fl, is_debug);
    } else if (msphere.is_on_se()) {
      fix_se_sphere(fe_abs_len, msphere, all_medial_spheres, se_tvs_grouped,
                    tvs_pos_map, fl2corner_sphere, new_centers_fl, is_debug);
    } else {
      num_changed += fix_non_extf_sphere(
          num_itr_global, fe_abs_len, msphere, se_tvs_grouped, tvs_pos_map,
          all_medial_spheres, fl2corner_sphere, new_centers_fl, is_debug);
    }

    // adding new feature spheres
    for (int cidx = 0; cidx < new_centers_fl.size(); cidx++) {
      v2int2 newc_fl_info = new_centers_fl[cidx];
      add_new_se_sphere(num_itr_global, all_medial_spheres.size(),
                        newc_fl_info.first, radius, newc_fl_info.second,
                        new_medial_spheres, is_debug);
    }  // for new_centers

  }  // for all_medial_spheres

  // add new spheres
  for (auto& new_msphere : new_medial_spheres) {
    if (add_new_sphere_validate(all_medial_spheres, new_msphere, is_debug)) {
      num_changed++;
      num_added++;
    }
  }
  printf(
      "[FIX_EXTF] changed %d including %d added new external feature spheres\n",
      num_changed, num_added);

  return num_changed;
}