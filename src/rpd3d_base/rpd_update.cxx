#include "rpd_update.h"

/////////////////////////////////////////////////////////////////////////////////////
// Update PowerCells
/////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief For each powercell, we collect its pc_face centroid that on
 * surface/boundary matches GEO::Mesh. Note that pc_face centroid is not
 * triangle centroid.
 *
 * @param max_surf_fid defines the maximum fid on GEO::Mesh
 * @param cc_trans
 * @param lfa powercell's local ACTIVE face id
 * @param sf_fid corresponding surface fid matches GEO::Mesh
 * @param result
 * @return true
 * @return false
 */
v2int get_cell_v2surffid(const uint max_surf_fid,
                         const ConvexCellHost& cc_trans, const uchar lfa,
                         int sf_fid) {
  assert(cc_trans.is_pc_explicit);
  assert(cc_trans.is_active);
  assert(lfa >= 0 && lfa <= cc_trans.nb_p);  // lfa <= lfid
  assert(sf_fid <= max_surf_fid);

  const auto& pc_points = cc_trans.pc_points;
  const auto& pc_local_active_faces = cc_trans.pc_local_active_faces;
  assert(lfa <= pc_local_active_faces.size());

  // find pc_face centroid
  cfloat3 centroid = {0.0f, 0.0f, 0.0f};  // centroid of facet
  const auto& pc_lf_vertices = pc_local_active_faces.at(lfa);
  for (int lv : pc_lf_vertices) {
    const auto& pc_vertex = pc_points.at(lv);
    centroid = cplus3(centroid, pc_vertex);
  }
  centroid = cdivide3(centroid, pc_lf_vertices.size());
  Vector3 p(centroid.x, centroid.y, centroid.z);
  return std::pair<Vector3, int>(p, sf_fid);
}

aint2 convert_e_lfids_to_lvids(
    const std::vector<std::vector<int>>& pc_local_active_faces,
    const std::vector<int>& pc_lf2active_map, const aint2& e_lfids,
    bool is_debug) {
  assert(e_lfids[0] < pc_lf2active_map.size() &&
         e_lfids[1] < pc_lf2active_map.size());
  int lafid1 = pc_lf2active_map.at(e_lfids[0]);
  int lafid2 = pc_lf2active_map.at(e_lfids[1]);
  assert(lafid1 != -1 && lafid2 != -1);
  const auto& lf1_vs = to_set(pc_local_active_faces.at(lafid1));
  const auto& lf2_vs = to_set(pc_local_active_faces.at(lafid2));
  if (is_debug) {
    printf("[UpdateVoro] lf1 %d, laf1: %d has lvs: (", e_lfids[0], lafid1);
    for (const auto& lvs1 : lf1_vs) printf("%d, ", lvs1);
    printf(")\n");
    printf("[UpdateVoro] lf2 %d laf2: %d, has lvs: (", e_lfids[1], lafid2);
    for (const auto& lvs2 : lf2_vs) printf("%d, ", lvs2);
    printf(")\n");
  }

  std::vector<int> common_vs;
  set_intersection<int>(lf1_vs, lf2_vs, common_vs);
  assert(common_vs.size() == 2);
  aint2 result = {{common_vs[0], common_vs[1]}};
  return get_sorted(result);
}

/**
 * @brief Given all convex cells belong to a power cell, we want to update the
 * info about their adjacencies.
 *
 * @param voro_convex_cells
 * @param powercell
 * @param tet_es2fe_map from TetMesh::tet_es2fe_map, format <tid, lfid_min,
 * lfid_max> -> <fe_type, fe_id, fe_line_id>.
 */
void get_all_voro_info(const std::vector<ConvexCellHost>& voro_convex_cells,
                       PowerCell& powercell, const uint max_surf_fid,
                       const std::map<aint3, aint3>& tet_es2fe_map) {
  // cell_id -> {neighboring cell ids}
  auto& cell_neighbors = powercell.cell_neighbors;
  // original tet fid -> 1 or 2 cell ids
  // auto& tfid_to_cells = powercell.tfid_to_cells;
  std::map<int, std::set<int>> tfid_to_cells;
  // cell id -> orignal tet fids
  auto& cell_to_tfids = powercell.cell_to_tfids;
  // halfplane seed_neigh_id -> list of cell ids
  auto& facet_neigh_to_cells = powercell.facet_neigh_to_cells;
  // halfplane [seed_neigh_min, seed_neigh_max] -> list of cell ids
  auto& e_to_cells = powercell.e_to_cells;
  // [neigh_id_min, neigh_id_max] ->
  // {a set of convex edges [aint2, aint2] all powercell edges}
  // aint2 is [cell_id, lvid (from cc_trans.nb_v)]
  auto& edge_2endvertices = powercell.edge_2endvertices;
  // <pc_face_centroid, surf_fid> pair (<Vector3, int>), all matching GEO::Mesh
  auto& cell_to_surfv2fid = powercell.cell_to_surfv2fid;
  // all touched sharp edges, store [cell_id, lvid1, lvid2, fe_line_id, fe_id]
  auto& se_covered_lvids = powercell.se_covered_lvids;
  // each sharp line endpoint has a unique neigh_id
  // [cell_id, lvid, neigh_id, se_line_id] -> pos
  auto& se_line_endpos = powercell.se_line_endpos;
  // all touched concave edges, store [cell_id, lvid1, lvid2, fe_line_id, fe_id]
  auto& ce_covered_lvids = powercell.ce_covered_lvids;
  // [cell_id, lvid (from cc_trans.nb_v)] -> [neigh_id1, neigh_id2, neigh_id3]
  auto& vertex_2id = powercell.vertex_2id;
  // [cell_id, lvid (from cc_trans.nb_v)] -> <pos, sf_fid>
  auto& vertex_2pos = powercell.vertex_2pos;

  // collect info of vertex/facet/edge <-> cell_ids
  for (uint i = 0; i < voro_convex_cells.size(); i++) {
    auto& cc_trans = voro_convex_cells[i];
    int cell_id = cc_trans.id;
    assert(cc_trans.is_active_updated);
    const auto& active_clipping_planes = cc_trans.active_clipping_planes;

    // powercell facet <-> cell
    int lfa = 0;  // local active fid (matches pc_local_active_faces)
    FOR(plane, cc_trans.nb_p) {
      // not active
      if (active_clipping_planes[plane] <= 0) continue;
      cint2 hp = cc_trans.clip_id2_const(plane);
      if (hp.y != -1) {  // store halfplanes
        int neigh_seed_id = hp.x == cc_trans.voro_id ? hp.y : hp.x;
        facet_neigh_to_cells[neigh_seed_id].insert(cell_id);
      } else {
        // centroid-facet pairs on surface (matching GEO::Mesh)
        if (hp.x <= max_surf_fid) {
          v2int v2fid_pair =
              get_cell_v2surffid(max_surf_fid, cc_trans, lfa, hp.x);
          cell_to_surfv2fid[cell_id].push_back(v2fid_pair);
        }
        // not halfplane, facet from orignal tet
        tfid_to_cells[hp.x].insert(cell_id);
        cell_to_tfids[cell_id].insert(hp.x);
      }
      lfa++;
    }

    // powercell vertex -> pos
    // [neigh_seed1, neigh_seed2, neigh_seed3],
    // where neigh_seed1 can be -1, then v is on surface
    //
    // the vertex key is [cell_id, lvid (from cc_trans.nb_v)]
    aint3 v_unique_id = {{-1, -1, -1}};
    std::set<int> v_neigh_seeds;
    FOR(vid, cc_trans.nb_v) {
      v_neigh_seeds.clear();
      cuchar4 v = cc_trans.ver_trans_const(vid);
      int surf_fid = -1;
      // store neighboring seeds & surface fid
      FOR(i, 3) {
        uchar lf = cc_trans.ith_plane_const(vid, i);
        assert(active_clipping_planes[lf] > 0);
        cint2 hp = cc_trans.clip_id2_const(lf);
        if (hp.y != -1) {  // halfplane
          int neigh_seed_id = hp.x == cc_trans.voro_id ? hp.y : hp.x;
          v_neigh_seeds.insert(neigh_seed_id);
        } else {  // store surface face id
          if (hp.x <= max_surf_fid) surf_fid = hp.x;
        }
      }
      // Note:
      // 1. if v_neigh_seeds.size() == 3,
      //    then vid is inside the shape dual to a medial tet
      // 2. if v_neigh_seeds.size() == 2,
      //    then vid is on the surface
      // 3. if v_neigh_seeds.size() < 2,
      //    not a powercell vertex, just a convex cell vertex
      //    or on sharp edges, do not store
      if (v_neigh_seeds.size() < 2) continue;
      if (v_neigh_seeds.size() == 2) v_neigh_seeds.insert(-1);
      int i = 0;
      for (const auto& n : v_neigh_seeds) v_unique_id[i++] = n;
      std::sort(v_unique_id.begin(), v_unique_id.end());
      Vector3 point = to_vec(cc_trans.compute_vertex_coordinates(
          cmake_uchar3(cc_trans.ver_trans_const(vid))));
      aint2 key = {{cc_trans.id, vid}};  // [cell_id, lvid (from cc_trans.nb_v)]
      vertex_2id[key] = v_unique_id;
      // surf_fid can be -1 if not on surface
      vertex_2pos[key] = std::pair<Vector3, int>(point, surf_fid);
      // if (cell_id == 16394) {
      //   printf(
      //       "[MVertex] msphere %d has powercell vertex v_unique_id:"
      //       "(%d,%d,%d), in convex cell %d and lvid: %d \n ",
      //       cc_trans.voro_id, v_unique_id[0], v_unique_id[1], v_unique_id[2],
      //       cc_trans.id, vid);
      // }
    }  // for cc_trans.nb_v

    // powercell edge <-> cell
    aint2 e_unique_id = {{-1, -1}};  // [neigh_seed_min, neigh_seed_max]
    FOR(eid, cc_trans.nb_e) {
      if (cc_trans.active_edges[eid] <= 0) continue;
      // find [seed_neigh_min, seed_neigh_max] as edge's unique id
      cuchar2 e_fid2 = cc_trans.edge_id2_const(eid);
      cint2 hp1 = cc_trans.clip_id2_const(e_fid2.x);  // first halfplane
      cint2 hp2 = cc_trans.clip_id2_const(e_fid2.y);  // second halfplane

      // if (cell_id == 19967) {
      //   printf("cell %d has edge (%d,%d) \n", cell_id, hp1.x, hp2.x);
      // }

      // store covered feature edges (sharp/concave)
      // bool is_debug = cc_trans.voro_id == 1707 ? true : false;
      bool is_debug = false;
      if (hp1.y == -1 && hp2.y == -1) {  // edge on orignal tet
        // store local fid, not tet fid
        aint2 edge_lfids = {{e_fid2.x, e_fid2.y}};
        std::sort(edge_lfids.begin(), edge_lfids.end());
        aint3 tlfs = {{cc_trans.tet_id, edge_lfids[0], edge_lfids[1]}};
        // if (cell_id == 147) {
        //   printf("[UpdateVoro] cell %d has tet %d edge_lfids (%d,%d) \n",
        //          cell_id, cc_trans.tet_id, edge_lfids[0], edge_lfids[1]);
        // }
        if (tet_es2fe_map.find(tlfs) == tet_es2fe_map.end()) continue;
        // found a feature edge, get lvids of the edge
        aint2 edge_lvids = convert_e_lfids_to_lvids(
            cc_trans.pc_local_active_faces, cc_trans.pc_lf2active_map,
            edge_lfids, is_debug);
        aint3 fe_info = tet_es2fe_map.at(tlfs);  // <fe_type, fe_id, fe_line_id>
        // format [cell_id, lvid1, lvid2, fe_line_id, fe_id]
        aint5 key_tmp = {
            {cell_id, edge_lvids[0], edge_lvids[1], fe_info[2], fe_info[1]}};

        if (fe_info[0] == EdgeType::SE)
          // found sharp edge, store
          se_covered_lvids.insert(key_tmp);
        else if (fe_info[0] == EdgeType::CE)
          // found concave edge, store
          ce_covered_lvids.insert(key_tmp);
        //////
        // store the endpoint of sharp line if exists only 1 halfplane
        FOR(lv, 2) {
          int lvid = edge_lvids[lv];
          FOR(i, 3) {  // check 3 local faces
            cint2 hp_i =
                cc_trans.clip_id2_const(cc_trans.ith_plane_const(lvid, i));
            if (hp_i.y == -1) continue;  // not a halfplane
            int neigh_seed_id = hp_i.x == cc_trans.voro_id ? hp_i.y : hp_i.x;
            Vector3 point = to_vec(cc_trans.compute_vertex_coordinates(
                cmake_uchar3(cc_trans.ver_trans_const(lvid))));
            se_line_endpos[{{cell_id, lvid, neigh_seed_id, fe_info[2]}}] =
                point;
            if (is_debug)
              printf("se_line_endpos has position: (%f,%f,%f) \n", point[0],
                     point[1], point[2]);
            break;
          }
        }
        if (is_debug)
          printf(
              "[UpdateVoro] cell %d has tet %d sharp edge in lfids (%d,%d), in "
              "lvids: (%d,%d) \n",
              cell_id, cc_trans.tet_id, edge_lfids[0], edge_lfids[1],
              edge_lvids[0], edge_lvids[1]);
      }  // done edge on surface boundary

      // only cares about halfplanes
      if (hp1.y == -1 || hp2.y == -1) continue;
      e_unique_id[0] = hp1.x == cc_trans.voro_id ? hp1.y : hp1.x;
      e_unique_id[1] = hp2.x == cc_trans.voro_id ? hp2.y : hp2.x;
      std::sort(e_unique_id.begin(), e_unique_id.end());
      e_to_cells[e_unique_id].insert(cell_id);

      // update edge_2endvertices
      // the vertex key is [cell_id, lvid (from cc_trans.nb_v)]
      std::vector<int> end_vertices;
      int e_fid1_active = cc_trans.pc_lf2active_map.at(e_fid2.x);
      int e_fid2_active = cc_trans.pc_lf2active_map.at(e_fid2.y);
      assert(e_fid1_active != -1 && e_fid2_active != -1);
      const auto& e_f1_vs =
          to_set(cc_trans.pc_local_active_faces.at(e_fid1_active));
      const auto& e_f2_vs =
          to_set(cc_trans.pc_local_active_faces.at(e_fid2_active));
      set_intersection(e_f1_vs, e_f2_vs, end_vertices);
      assert(end_vertices.size() == 2);
      // make sure two end vertices are on powercell edge, and exists in
      // vertex_2id, may be convex cell vertices
      aint2 end_v1_key = {{cc_trans.id, end_vertices[0]}};
      aint2 end_v2_key = {{cc_trans.id, end_vertices[1]}};
      if (vertex_2id.find(end_v1_key) == vertex_2id.end()) {
        printf(
            "[UpdateVoro] cell %d, voro_id %d, tet %d, end_v1_key: (%d,%d), "
            "end_v2_key: (%d,%d), vertex_2id.size %ld \n",
            cc_trans.id, cc_trans.voro_id, cc_trans.tet_id, end_v1_key[0],
            end_v1_key[1], end_v2_key[0], end_v2_key[1], vertex_2id.size());
        cc_trans.print_info();
        assert(false);
      }
      assert(vertex_2id.find(end_v2_key) != vertex_2id.end());
      edge_2endvertices[e_unique_id].push_back({{end_v1_key, end_v2_key}});
    }  // for cc_trans.nb_e

    // if (cc_trans.id == 4039 ||
    //     (cc_trans.voro_id == 1483 && cc_trans.tet_id == 1107)) {
    //   cc_trans.print_info();
    // }
  }  // for voro_convex_cells

  // update cell_neighbors
  // cells neighboring inherited from orignal tet face is what we care
  // adjacency through halfplanes is for different seeds/powercells
  for (const auto& pair : tfid_to_cells) {
    auto& cells_set = pair.second;
    // printf("face %d has cells_set size %zu \n", pair.first,
    // cells_set.size());
    assert(cells_set.size() > 0 && cells_set.size() < 3);
    if (cells_set.size() == 1) continue;  // boundary, no neighbor cell
    int cell1 = *cells_set.begin();
    int cell2 = *cells_set.rbegin();
    cell_neighbors[cell1].insert(cell2);
    cell_neighbors[cell2].insert(cell1);
  }  // for tfid_to_cells

  // // print cell_neighbors
  // if (powercell.voro_id == 1478 || powercell.voro_id == 1483) {
  //   printf("--------------- powercell %d\n", powercell.voro_id);

  //   // print facet_neigh_to_cells
  //   for (const auto& pair : facet_neigh_to_cells) {
  //     printf("neigh_seed %d has cells: (", pair.first);
  //     for (const auto& c : pair.second) printf("%d, ", c);
  //     printf(")\n");
  //   }

  //   printf("----\n");
  //   for (const auto& pair : cell_neighbors) {
  //     // if (pair.first != 4039) continue;
  //     printf("cell %d has neighbors: (", pair.first);
  //     for (const auto& n : pair.second) printf("%d, ", n);
  //     printf(")\n");
  //   }
  // }

  // printf("voro_convex_cells: %ld, cell_neighbors: %ld \n",
  //        voro_convex_cells.size(), cell_neighbors.size());
}

// Mainly for updating msphere.pcell.surf_v2fid_in_groups
void update_sphere_surf_v2fids(const SurfaceMesh& sf_mesh,
                               MedialSphere& msphere, bool is_debug) {
  // already did in MedailSphere::topo_clear()
  // msphere.pcell.surf_v2fid_in_groups.clear();
  const std::map<int, std::vector<v2int>>& cell_to_surfv2fid =
      msphere.pcell.cell_to_surfv2fid;
  if (cell_to_surfv2fid.empty()) return;
  std::vector<std::vector<v2int>>& surf_v2fid_in_groups =
      msphere.pcell.surf_v2fid_in_groups;

  const auto& fe_sf_pairs_not_cross = sf_mesh.fe_sf_fs_pairs;
  auto is_skip_se_neighbor = [&](const int f, const int nf) {
    aint2 ref_fs_pair = {{f, nf}};
    std::sort(ref_fs_pair.begin(), ref_fs_pair.end());
    if (fe_sf_pairs_not_cross.find(ref_fs_pair) !=
        fe_sf_pairs_not_cross.end()) {
      if (is_debug)
        printf("[UpdateVoro] face %d skip nf %d since sharing a sharp edge \n",
               f, nf);
      return true;  // skip checking its neighbor
    }
    return false;
  };

  // store surface fids in powercell of msphere
  std::set<int> sf_fids_to_visit;
  std::map<int, v2int> sf_fid2v_map;
  for (const auto& pair : cell_to_surfv2fid) {
    const auto& one_group_v2fids = pair.second;
    for (const auto& one_v2fids : one_group_v2fids) {
      int fid = one_v2fids.second;
      sf_fids_to_visit.insert(fid);
      sf_fid2v_map[fid] = one_v2fids;
    }  // one_group_v2fids
  }

  std::map<int, std::set<int>> sf_fid_neighbors;
  for (const auto& fid : sf_fids_to_visit) {
    for (GEO::index_t le = 0; le < sf_mesh.facets.nb_vertices(fid); le++) {
      GEO::index_t nfid = sf_mesh.facets.adjacent(fid, le);
      if (nfid == GEO::NO_FACET) continue;
      if (is_skip_se_neighbor(fid, nfid)) continue;
      // if not in sf_fids_to_visit, then skip
      if (sf_fids_to_visit.find(nfid) == sf_fids_to_visit.end()) continue;
      // if (is_debug)
      //   printf("[UpdateVoro] fid %d has neighbor face nfid %d\n", fid, nfid);
      sf_fid_neighbors[fid].insert(nfid);
    }  // for facets.nb_vertices
  }  // for sf_fids_to_visit

  // find fids in groups
  std::vector<std::set<int>> covered_sf_fids_in_groups;
  get_CC_given_neighbors(sf_fids_to_visit, sf_fid_neighbors,
                         covered_sf_fids_in_groups);
  // store back to surf_v2fid_in_groups
  std::vector<v2int> one_v2fid_group;
  for (const auto& one_fid_group : covered_sf_fids_in_groups) {
    one_v2fid_group.clear();
    for (const int& fid : one_fid_group) {
      one_v2fid_group.push_back(sf_fid2v_map.at(fid));
    }
    surf_v2fid_in_groups.push_back(one_v2fid_group);
  }
}

// helper function for update_pc_facet_cc_info()
void get_facet_CC_surf_v2fids(
    const std::map<int, std::vector<v2int>>& cell_to_surfv2fid,
    const std::map<int, std::vector<std::set<int>>>& facet_cc_cells,
    std::map<int, std::vector<std::vector<v2int>>>& facet_cc_surf_v2fids) {
  facet_cc_surf_v2fids.clear();
  // save surface centroid_fid pairs for each facet CC, groupd in one
  // neigh_id/halfplane
  //
  // Note: all fid in vertex_fid pairs matches GEO::Mesh
  std::vector<v2int> surf_v2fids;
  for (const auto& pair : facet_cc_cells) {
    int neigh_id = pair.first;
    auto& one_halfplane_cells = pair.second;
    for (const auto& one_cc_cells : one_halfplane_cells) {
      surf_v2fids.clear();
      for (int cell_id : one_cc_cells) {
        if (cell_to_surfv2fid.find(cell_id) == cell_to_surfv2fid.end())
          continue;
        const auto& v2fids = cell_to_surfv2fid.at(cell_id);
        surf_v2fids.insert(surf_v2fids.end(), v2fids.begin(), v2fids.end());
      }
      facet_cc_surf_v2fids[neigh_id].push_back(surf_v2fids);
    }
  }
}

// update
// 1. PowerCell::facet_cc_cells
// 2. PowerCell::facet_cc_surf_v2fids
// call after get_all_voro_info()
void update_pc_facet_cc_info(PowerCell& pcell, const int msphere_id) {
  assert(!pcell.facet_neigh_to_cells.empty());
  for (auto& pair : pcell.facet_neigh_to_cells) {
    // halfplane defined by [msphere.id, neigh_id]
    const int neigh_id = pair.first;
    const std::set<int>& halfplane_cells = pair.second;

    int num_facet_cc = get_CC_given_neighbors<int>(
        halfplane_cells, pcell.cell_neighbors, pcell.facet_cc_cells[neigh_id]);
    // sort pcell.facet_cc_cells[neigh_id] based on the size of group (in set)
    auto& new_name = pcell.facet_cc_cells[neigh_id];
    std::sort(new_name.begin(), new_name.end(),
              [](const std::set<int>& a, const std::set<int>& b) {
                return a.size() < b.size();
              });
    get_facet_CC_surf_v2fids(pcell.cell_to_surfv2fid, pcell.facet_cc_cells,
                             pcell.facet_cc_surf_v2fids);

    if (msphere_id == 377 && neigh_id == 649 ||
        msphere_id == 649 && neigh_id == 377) {
      printf("[UpdatePCFaceCC] sid %d, neigh_id %d, has cells: [", msphere_id,
             neigh_id);
      for (int cid : halfplane_cells) {
        printf("%d, ", cid);
      }
      printf("], num_facet_cc: %d\n", num_facet_cc);
    }
  }
}

// helper function for update_pc_cc_info()
void get_CC_surf_vfids(
    const std::map<int, std::vector<v2int>>& cell_to_surfv2fid,
    const std::vector<std::set<int>>& cc_cells,
    std::vector<std::vector<v2int>>& cc_surf_v2fids) {
  cc_surf_v2fids.clear();
  // printf("cell_to_surfv2fid size: %ld\n", cell_to_surfv2fid.size());
  // save surface vertex2fid mapping for each CC
  // Note: we only care about pc's face centroid on surface, so its fids match
  // GEO::Mesh Note: each vertex only map to one surface fid (good enough)
  std::vector<v2int> surf_v2fids;
  for (const auto& one_cc_cells : cc_cells) {
    surf_v2fids.clear();
    for (int cell_id : one_cc_cells) {
      if (cell_to_surfv2fid.find(cell_id) == cell_to_surfv2fid.end()) continue;
      const auto& v2fid_pairs = cell_to_surfv2fid.at(cell_id);
      surf_v2fids.insert(surf_v2fids.end(), v2fid_pairs.begin(),
                         v2fid_pairs.end());
    }
    cc_surf_v2fids.push_back(surf_v2fids);
  }
}

// update:
// 1. PowerCell::cc_cells
// 2. PowerCell::cc_surf_v2fids
// call after get_all_voro_info()
void update_pc_cc_info(PowerCell& pcell) {
  // calculate number of CC
  int num_cc = get_CC_given_neighbors<int>(pcell.cell_ids, pcell.cell_neighbors,
                                           pcell.cc_cells);
  get_CC_surf_vfids(pcell.cell_to_surfv2fid, pcell.cc_cells,
                    pcell.cc_surf_v2fids);
  assert(num_cc > 0);
}

// update
// PowerCell::edge_cc_cells
// call after get_all_voro_info()
void update_pc_edge_cc_info(PowerCell& pcell) {
  assert(!pcell.facet_neigh_to_cells.empty());
  for (auto& pair : pcell.e_to_cells) {
    // two halfplane defined by:
    // 1. [msphere.id, neigh_id_min]
    // 2. [msphere.id, neigh_id_max]
    const aint2 neigh_min_max = pair.first;
    const std::set<int>& halfplane_cells = pair.second;
    int num_edge_cc =
        get_CC_given_neighbors<int>(halfplane_cells, pcell.cell_neighbors,
                                    pcell.edge_cc_cells[neigh_min_max]);
    assert(num_edge_cc > 0);
  }
}

// [no use]
void update_sphere_type_based_on_pcells(const SurfaceMesh& input_mesh,
                                        MedialSphere& msphere, bool is_debug) {
  if (msphere.type == SphereType::T_3_MORE) {
    // TODO: use msphere.pcell.surf_v2fid_in_groups instead!!
    //
    // group pcells' surface fids in CC not cross sharp edges
    for (const auto& one_surf_cc : msphere.pcell.cc_surf_v2fids) {
      std::set<int> facets_to_group;
      for (const v2int& tmp : one_surf_cc) {
        facets_to_group.insert(tmp.second);  // CC contains surface fids
      }
      std::vector<std::set<int>> cc_facets;
      get_CC_given_neighbors<int>(facets_to_group,
                                  input_mesh.sf_fid_neighs_no_cross, cc_facets);
      int num_surf_cc = cc_facets.size();
      int old_type = msphere.type;
      if (num_surf_cc < 3) {
        msphere.type = SphereType::T_2;
        if (is_debug) {
          printf(
              "[Pcell CC] update msphere %d type %d->%d with num_surf_cc: %d, "
              "contains surface fid in CC: \n",
              msphere.id, old_type, msphere.type, num_surf_cc);
          for (const auto& fset : cc_facets) {
            printf("[");
            for (const auto& fid : fset) {
              printf("%d, ", fid);
            }
            printf("]\n");
          }
        }
      }  // num_surf_cc < 3
    }  // for msphere.pcell.cc_surf_v2fids
  }

  // TODO: update T_N_c
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// map_site2msphere: site id -> MedialSphere::id
// map_tet_new2orig: tet id partial/new -> old
// TetMesh::tet_es2fe_map:
// format  <tid, lfid_min, lfid_max> -> <fe_type, fe_id, fe_line_id>.
void update_power_cells(const SurfaceMesh& sf_mesh,
                        std::vector<ConvexCellHost>& convex_cells_host,
                        std::vector<MedialSphere>& all_medial_spheres,
                        const std::map<aint3, aint3>& tet_es2fe_map,
                        bool is_debug) {
  // for parallization
  auto run_thread = [&](int mid, const uint max_surf_fid) {
    MedialSphere& msphere = all_medial_spheres.at(mid);
    if (msphere.is_deleted) return;
    const std::set<int>& all_cell_ids = msphere.pcell.cell_ids;
    if (all_cell_ids.empty()) {
      if (is_debug)
        printf("[UpdateVoro] sphere %d has zero convex cells!!\n", msphere.id);
      msphere.is_deleted = true;  // delete?
      return;
    }

    // collect all convex cells and update their adjacencies
    std::vector<ConvexCellHost> all_ccells;
    for (uint cell_id : all_cell_ids) {
      auto& convex_cell = convex_cells_host.at(cell_id);
      all_ccells.push_back(convex_cell);
    }

    msphere.pcell.voro_id = msphere.id;
    get_all_voro_info(all_ccells, msphere.pcell, max_surf_fid, tet_es2fe_map);
    update_sphere_surf_v2fids(sf_mesh, msphere, is_debug);
    update_pc_cc_info(msphere.pcell);
    update_pc_facet_cc_info(msphere.pcell, msphere.id);
    update_pc_edge_cc_info(msphere.pcell);
    // // TODO: use msphere.pcell.surf_v2fid_in_groups instead!!
    // update_sphere_type_based_on_pcells(sf_mesh, msphere, is_debug);
    msphere.update_sphere_covered_sf_fids(sf_mesh, is_debug);
  };

  printf("updating power cells ... using all #convex_cells: %zu \n",
         convex_cells_host.size());

  const uint max_surf_fid = sf_mesh.facets.nb() - 1;
  // clear existing power cells
  for (auto& msphere : all_medial_spheres) {
    msphere.topo_clear();
  }

  // seed_id -> {id of convex_cells_host}
  for (uint i = 0; i < convex_cells_host.size(); i++) {
    auto& convex_cell = convex_cells_host.at(i);
    int cell_id = convex_cell.id;
    assert(is_convex_cell_valid(convex_cell));
    auto& msphere = all_medial_spheres.at(convex_cell.voro_id);
    msphere.pcell.cell_ids.insert(cell_id);
    msphere.pcell.tet_ids.insert(convex_cell.tet_id);

    // if (convex_cell.voro_id == 925 && convex_cell.tet_id == 4226) {
    //   printf("++++++ voro_id %d and tet_id %d has cell %d\n",
    //          convex_cell.voro_id, convex_cell.tet_id, convex_cell.id);
    // }

    // some clipping planes may not exist in tri but
    // we still store it, here is to filter those planes
    if (!convex_cell.is_active_updated) convex_cell.reload_active();
    if (!convex_cell.is_pc_explicit) convex_cell.reload_pc_explicit();
  }

  // update
  // seed_id -> vector of convex cells
  GEO::parallel_for(0, all_medial_spheres.size(),
                    [&](int i) { run_thread(i, max_surf_fid); });
  // for (int mid = 0; mid < all_medial_spheres.size(); mid++) {
  //   run_thread(mid, max_surf_fid);
  // }
}