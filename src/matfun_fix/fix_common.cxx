#include "fix_common.h"

/**
 * @brief convert sharp edge representation from convexcell [cell_id, lvid1,
 * lvid2, fe_line_id, fe_id] to tet [tvid1, tvid2, fe_line_id, fe_id]
 *
 * 1. if cell_lvid is from tet, then keep the tvid
 * 2. if cell_lvid is newly cutted, then assign negative id (vs_new_ids)
 *
 * for multiple convexcells, multple id may share the same vs_new_ids,
 * but who cares.
 *
 * @param msphere no use, only for debug
 * @param tet_vertices
 * @param convex_cells_host
 * @param se_covered_lvids se_lvids format: [cell_id, lvid1, lvid2,
 * fe_line_id, fe_id]
 * @param fe_tvs_all return, format [tvid1/vs_new_id1, tvid2/vs_new_id2,
 * fe_line_id, fe_id]
 * @param tvs_pos_map return, tvid/vs_new_id -> Vector3
 */
void convert_fes_from_ccell_to_tet(
    const MedialSphere& msphere, const std::vector<float>& tet_vertices,
    const std::vector<ConvexCellHost>& convex_cells_host,
    const std::map<aint4, int>& tet_vs_lfs2tvs_map,
    const std::set<aint5>& se_covered_lvids, std::set<aint4>& fe_tvs_all,
    std::map<int, Vector3>& tvs_pos_map, bool is_debug) {
  // [tvid1/vs_new_id1, tvid2/vs_new_id2, fe_line_id, fe_id]
  fe_tvs_all.clear();
  tvs_pos_map.clear();  // tvid/vs_new_ids -> Vector3
  int vs_new_ids = -2;  // using negative, avoid UNK_INT(-1)
  // se_lvids format: [cell_id, lvid1, lvid2, fe_line_id, fe_id]
  for (const aint5& se_lvids : se_covered_lvids) {
    int fe_line_id = se_lvids[3];  // index 3
    int fe_id = se_lvids[4];       // index 4
    aint2 se_vs = {{UNK_INT, UNK_INT}};
    const auto& cc_trans = convex_cells_host.at(se_lvids[0]);  // index 0
    for (int i = 0; i < 2; i++) {
      int lvid = se_lvids[i + 1];  // index 1 or 2
      cuchar4 v = cc_trans.ver_trans_const(lvid);
      aint3 lfids = {{v.x, v.y, v.z}};
      lfids = get_sorted(lfids);
      aint4 v_key = {{cc_trans.tet_id, lfids[0], lfids[1], lfids[2]}};
      if (is_debug)
        printf(
            "[Convert] msphere %d has cell %d lvid %d in fe_line_id "
            "%d has v_key (%d,%d,%d,%d) \n",
            msphere.id, se_lvids[0], lvid, fe_line_id, v_key[0], v_key[1],
            v_key[2], v_key[3]);

      if (tet_vs_lfs2tvs_map.find(v_key) != tet_vs_lfs2tvs_map.end()) {
        // from tet mesh, has tvid
        se_vs[i] = tet_vs_lfs2tvs_map.at(v_key);
        tvs_pos_map[se_vs[i]] = Vector3(tet_vertices.at(se_vs[i] * 3),
                                        tet_vertices.at(se_vs[i] * 3 + 1),
                                        tet_vertices.at(se_vs[i] * 3 + 2));
      } else {
        // newly created vertex, assign new index
        se_vs[i] = vs_new_ids--;
        tvs_pos_map[se_vs[i]] = to_vec(cc_trans.compute_vertex_coordinates(
            cmake_uchar3(cc_trans.ver_trans_const(lvid))));
      }
    }  // for 2 endpoints
    assert(se_vs[0] != UNK_INT && se_vs[1] != UNK_INT);
    std::sort(se_vs.begin(), se_vs.end());
    fe_tvs_all.insert({{se_vs[0], se_vs[1], fe_line_id, fe_id}});

    if (is_debug)
      printf(
          "[Convert] msphere %d has cell %d covers feature edge in lvids: "
          "(%d,%d), tvs:(%d,%d), fe_line_id: %d, fe_id: %d\n",
          msphere.id, se_lvids[0], se_lvids[1], se_lvids[2], se_vs[0], se_vs[1],
          fe_line_id, fe_id);
  }  // for one sharp edge in a cell
}

/**
 * @brief Get the grouped feature edges, grouped in geometric sorted order
 *
 * @param fe_to_visit input, a set of aint3, format of <tvs_min, tvs_max,
 * fe_line_id, fe_id>. ** No use of fe_id for now **.
 * @param fe_groups return
 * @return int
 */
int get_grouped_fes_in_order(const std::set<aint4>& fe_to_visit,
                             std::vector<std::vector<aint4>>& fe_groups) {
  fe_groups.clear();

  // split into groups based on fe_line_id
  // <tvid1, tvid2> -> <fe_line_id, fe_id>
  std::map<aint2, aint2> fe2group_map;
  // fe_line_id -> set of <tvid1, tvid2>
  std::map<int, std::set<aint2>> group2fe_map;
  for (const aint4 fe : fe_to_visit) {
    fe2group_map[{{fe[0], fe[1]}}] = {{fe[2], fe[3]}};
    group2fe_map[fe[2]].insert({{fe[0], fe[1]}});
  }

  // note: the 2nd element (fe_line_id) of aint3 in fe_groups needs to be
  // re-asigned
  std::vector<std::vector<aint2>> fe_groups_tmp, fe_group_dup;
  std::set<int> _;
  for (const auto& pair : group2fe_map) {
    // int fe_line_id = pair.first;
    auto& fe_to_visit_tmp = pair.second;
    get_grouped_feature_edges(fe_to_visit_tmp, _, fe_groups_tmp);
    fe_group_dup.insert(fe_group_dup.end(), fe_groups_tmp.begin(),
                        fe_groups_tmp.end());
  }
  // re-assign original fe_line_id in fe_groups
  // convert a set of aint2 -> aint4
  std::vector<aint4> one_group;
  for (auto& one_fe_group_dup : fe_group_dup) {
    one_group.clear();
    for (auto& one_fe_dup : one_fe_group_dup) {
      aint2& orig_ids = fe2group_map.at(one_fe_dup);
      one_group.push_back(
          {{one_fe_dup[0], one_fe_dup[1], orig_ids[0], orig_ids[1]}});
    }
    fe_groups.push_back(one_group);
  }
  return fe_groups.size();
}