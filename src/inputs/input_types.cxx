#include "input_types.h"

void ConcaveCorner::print_info() const {
  printf("---- ConcaveCorner corner_tvid: %d, adj_fe_ids: [", tvid);
  for (const auto adj_fe : adj_fe_ids) {
    printf("%d, ", adj_fe);
  }
  printf("]\n");
}

void FeatureLine::print_info() const {
  printf("------ FeatureLine: %d, length: %f\n", id, length);
  // fe_id
  printf("feature edges: (");
  for (const auto& fe_id : fe_ids) {
    printf("%d, ", fe_id);
  }
  printf(")\n");
  // tvs_neighbors
  printf("tvs_neighbors: \n");
  for (const auto& pair : tvs_neighbors) {
    printf("%d -> [%d, %d]\n", pair.first, pair.second[0], pair.second[1]);
  }
  printf("-------\n");
}

void FL_Sample::print_info() const {
  printf("+++ FL_Sample: fe_id: %d, point (%f,%f,%f)\n", fe_id, point[0],
         point[1], point[2]);
}

v2int FeatureLine::get_fe_endpoint_given_dir(const FeatureEdge& fe,
                                             int dir) const {
  assert(dir == 0 || dir == 1);
  int tvid_left = fe.t2vs_group[0];
  int tvid_right = fe.t2vs_group[1];
  Vector3 pos_left = fe.t2vs_pos[0];
  Vector3 pos_right = fe.t2vs_pos[1];
  const aint2& min_neighbors = tvs_neighbors.at(tvid_left);
  // tvid_left and tvid_right must be neighbors
  assert(min_neighbors[0] == tvid_right || min_neighbors[1] == tvid_right);
  if (min_neighbors[0] == tvid_right) {  // swap
    swap_in_place<int>(tvid_left, tvid_right);
    swap_in_place<Vector3>(pos_left, pos_right);
  }
  if (dir == 0) return std::pair<Vector3, int>(pos_left, tvid_left);
  if (dir == 1) return std::pair<Vector3, int>(pos_right, tvid_right);
  assert(false);
}

int FeatureLine::get_next_fe_given_dir(int fe_cur, int dir) {
  int id_found = binary_search(fe_ids, fe_cur);
  if (dir == 0)
    id_found--;
  else if (dir == 1)
    id_found++;
  else
    assert(false);
  // check bound
  if (id_found < 0 || id_found >= fe_ids.size()) return -1;
  return fe_ids.at(id_found);
}

void SurfaceMesh::reload_sf2tet_vs_mapping() {
  // load sf2tet_vs_mapping
  const GEO::Attribute<int> tet_vid_attr(this->vertices.attributes(),
                                         "tet_vid");
  sf2tet_vs_mapping.clear();
  sf2tet_vs_mapping.resize(this->vertices.nb());
  for (int v = 0; v < this->vertices.nb(); v++) {
    sf2tet_vs_mapping[v] = tet_vid_attr[v];
  }
}

void SurfaceMesh::cache_sf_fid_neighs() {
  // load input_mesh adjacent info
  sf_fid_neighs.clear();
  // std::map<int, std::set<int>> conn_tris;
  for (int i = 0; i < this->facets.nb(); i++) {
    // printf("Input mesh facet %d has neighbors: [", i);
    FOR(le, 3) {
      int nfid = this->facets.adjacent(i, le);
      sf_fid_neighs[i].insert(nfid);
      // printf("%d, ", nfid);
    }
    // printf("]\n");
  }
}

// must call after detect_mark_sharp_features()
void SurfaceMesh::cache_sf_fid_neighs_no_cross() {
  if (fe_sf_fs_pairs.empty()) return;
  // load input_mesh adjacent info, not cross sharp features
  sf_fid_neighs_no_cross.clear();
  for (int i = 0; i < this->facets.nb(); i++) {
    FOR(le, 3) {
      int nfid = this->facets.adjacent(i, le);
      aint2 fpair = {{i, nfid}};
      std::sort(fpair.begin(), fpair.end());
      if (fe_sf_fs_pairs.find(fpair) != fe_sf_fs_pairs.end()) continue;
      sf_fid_neighs_no_cross[i].insert(nfid);
    }
  }
}

// must call after detect_mark_sharp_features()
void SurfaceMesh::cache_sf_fid_krings_no_cross_se_only(const int k) {
  if (fe_sf_fs_pairs_se_only.empty()) return;
  // load input_mesh adjacent info, not cross sharp features
  sf_fid_krings_no_cross_se_only.clear();
  std::set<int> kring_neighbors;
  for (int fid = 0; fid < this->facets.nb(); fid++) {
    collect_kring_neighbors_given_fid_se_only(k, fid, kring_neighbors);
    sf_fid_krings_no_cross_se_only[fid] = kring_neighbors;
  }
}

bool SurfaceMesh::get_sf_fid_krings_no_cross_se_only(
    const int fid_given, std::set<int>& kring_neighbors) const {
  kring_neighbors.clear();
  if (sf_fid_krings_no_cross_se_only.empty()) return false;
  if (sf_fid_krings_no_cross_se_only.find(fid_given) ==
      sf_fid_krings_no_cross_se_only.end())
    return false;
  kring_neighbors = sf_fid_krings_no_cross_se_only.at(fid_given);
  return true;
}

bool SurfaceMesh::collect_kring_neighbors_given_fid(
    const int k, int tan_fid, std::set<int>& kring_neighbors) const {
  assert(tan_fid >= 0);
  kring_neighbors.clear();
  get_k_ring_neighbors_no_cross(*this, this->fe_sf_fs_pairs, tan_fid, k,
                                kring_neighbors, true /*is_clear_cur*/,
                                false /*is_debug*/);
  return kring_neighbors.empty() ? false : true;
}

bool SurfaceMesh::collect_kring_neighbors_given_fid_se_only(
    const int k, int tan_fid, std::set<int>& kring_neighbors) const {
  assert(tan_fid >= 0);
  kring_neighbors.clear();
  get_k_ring_neighbors_no_cross(*this, this->fe_sf_fs_pairs_se_only, tan_fid, k,
                                kring_neighbors, true /*is_clear_cur*/,
                                false /*is_debug*/);
  return kring_neighbors.empty() ? false : true;
}

void SurfaceMesh::update_fe_sf_fs_pairs_to_ce_id(
    const std::vector<FeatureEdge>& feature_edges) {
  this->fe_sf_fs_pairs_to_ce_id.clear();
  for (const auto& fe : feature_edges) {
    if (fe.type == EdgeType::CE) {
      aint2 fs_pair = {{fe.adj_sf_fs_pair[0], fe.adj_sf_fs_pair[1]}};
      std::sort(fs_pair.begin(), fs_pair.end());
      this->fe_sf_fs_pairs_to_ce_id[fs_pair] = fe.id;
    }
  }
}

void SurfaceMesh::collect_fid_centroids(
    const std::set<int>& given_fids, std::vector<v2int>& one_group_fids) const {
  one_group_fids.clear();
  for (const int fid : given_fids) {
    one_group_fids.push_back(
        std::pair<Vector3, int>(get_mesh_facet_centroid(*this, fid), fid));
  }
}

bool AABBWrapper::init_mesh_from_edges(const std::vector<float>& input_vertices,
                                       const std::vector<aint2>& edges,
                                       const std::vector<int>& fe_ids,
                                       GEO::Mesh& mesh) {
  assert(edges.size() == fe_ids.size());
  GEO::Attribute<int> attr_fe_ids(mesh.facets.attributes(), "fe_id");
  if (edges.empty()) {
    mesh.vertices.clear();
    mesh.vertices.create_vertices(1);
    mesh.vertices.point(0) = GEO::vec3(0, 0, 0);
    mesh.facets.clear();
    mesh.facets.create_triangles(1);
    mesh.facets.set_vertex(0, 0, 0);
    mesh.facets.set_vertex(0, 1, 0);
    mesh.facets.set_vertex(0, 2, 0);
    return false;
  } else {
    mesh.vertices.clear();
    mesh.vertices.create_vertices((int)edges.size() * 2);
    int cnt = 0;
    for (auto& e : edges) {
      for (int j = 0; j < 2; j++) {
        GEO::vec3& p = mesh.vertices.point(cnt++);
        p[0] = input_vertices[e[j] * 3 + 0];
        p[1] = input_vertices[e[j] * 3 + 1];
        p[2] = input_vertices[e[j] * 3 + 2];
      }
    }
    mesh.facets.clear();
    mesh.facets.create_triangles(
        (int)edges.size());  // degenerated facets -> edge
    for (int i = 0; i < edges.size(); i++) {
      mesh.facets.set_vertex(i, 0, i * 2);
      mesh.facets.set_vertex(i, 1, i * 2);
      mesh.facets.set_vertex(i, 2, i * 2 + 1);
      attr_fe_ids[i] = fe_ids.at(i);
    }
  }
  // mesh_reorder(mesh, GEO::MESH_ORDER_MORTON);
  return true;
}

// is_reorder == true => will reorder sf_mesh
void AABBWrapper::init_sf_mesh_and_tree(GEO::Mesh& _sf_mesh, bool is_reorder) {
  sf_tree = std::make_shared<GEO::MeshFacetsAABB>(_sf_mesh, is_reorder);
}

void AABBWrapper::init_feature_meshes_and_trees(
    const std::vector<float>& input_vertices,
    const std::vector<FeatureEdge>& feature_edges) {
  std::vector<aint2> se_edges, ce_edges;
  std::vector<int> se_feids, ce_feids;
  for (const auto& fe : feature_edges) {
    if (fe.type == EdgeType::SE) {
      se_edges.push_back({{fe.t2vs_group[0], fe.t2vs_group[1]}});
      se_feids.push_back(fe.id);
    } else if (fe.type == EdgeType::CE) {
      ce_edges.push_back({{fe.t2vs_group[0], fe.t2vs_group[1]}});
      ce_feids.push_back(fe.id);
    }
  }
  printf("se_edges size %ld, ce_edges: %ld \n", se_edges.size(),
         ce_edges.size());
  // for sharp edges
  se_mesh.clear(false, false);
  is_se_mesh_exist =
      init_mesh_from_edges(input_vertices, se_edges, se_feids, se_mesh);
  if (!se_edges.empty()) assert(is_se_mesh_exist);
  se_tree = std::make_shared<GEO::MeshFacetsAABB>(se_mesh, true /*reorder*/);

  // for concave edges
  ce_mesh.clear(false, false);
  is_ce_mesh_exist =
      init_mesh_from_edges(input_vertices, ce_edges, ce_feids, ce_mesh);
  if (!ce_edges.empty()) assert(is_ce_mesh_exist);
  ce_tree = std::make_shared<GEO::MeshFacetsAABB>(ce_mesh, true /*reorder*/);
}

/**
 * tet2sf_vs_mapping:   id mapping from tet_vertices to sf_mesh
 * sf_vs2fids:          sf_mesh vertices to adjacent facets
 * tet_vs2sf_fids:      tet vertex to sf_mesh fids
 */
void load_sf_tet_mapping(const GEO::Mesh& sf_mesh,
                         std::map<int, int>& tet2sf_vs_mapping,
                         std::map<int, std::set<int>>& sf_vs2fids,
                         std::map<int, std::set<int>>& tet_vs2sf_fids) {
  tet2sf_vs_mapping.clear();
  sf_vs2fids.clear();
  const GEO::Attribute<int> tet_vid_attr(sf_mesh.vertices.attributes(),
                                         "tet_vid");
  for (int v = 0; v < sf_mesh.vertices.nb(); v++) {
    tet2sf_vs_mapping[tet_vid_attr[v]] = v;
  }

  for (int f = 0; f < sf_mesh.facets.nb(); f++) {
    int f_nb_v = sf_mesh.facets.nb_vertices(f);
    assert(f_nb_v == 3);
    for (int lv = 0; lv < f_nb_v; lv++) {
      int vid = sf_mesh.facets.vertex(f, lv);
      int tvid = tet_vid_attr[vid];
      sf_vs2fids[vid].insert(f);
      tet_vs2sf_fids[tvid].insert(f);
    }
  }
}

// Given fid_given, we fetch the k_ring neighboring faces
// not crossing sharp edges
void get_k_ring_neighbors_no_cross(const GEO::Mesh& sf_mesh,
                                   const std::set<aint2>& fe_sf_pairs_not_cross,
                                   const int fid_given, const int k,
                                   std::set<int>& k_ring_fids,
                                   bool is_clear_cur, bool is_debug) {
  auto is_skip_se_neighbor = [&](const int f, const int nf) {
    aint2 ref_fs_pair = {{f, nf}};
    std::sort(ref_fs_pair.begin(), ref_fs_pair.end());
    if (fe_sf_pairs_not_cross.find(ref_fs_pair) !=
        fe_sf_pairs_not_cross.end()) {
      if (is_debug)
        printf("[K_RING] face %d skip nf %d since sharing a sharp edge \n", f,
               nf);
      return true;  // skip checking its neighbor
    }
    return false;
  };
  if (is_debug)
    printf("calling get_k_ring_neighbor_no_se for fid_given: %d \n", fid_given);
  assert(fid_given >= 0);

  if (is_clear_cur) k_ring_fids.clear();
  k_ring_fids.insert(fid_given);
  std::set<int> new_added_fids, new_added_fids_copy;
  new_added_fids.insert(fid_given);

  int i = 0;
  while (i < k) {
    new_added_fids_copy = new_added_fids;
    new_added_fids.clear();
    for (const auto& fid : new_added_fids_copy) {
      for (GEO::index_t le = 0; le < sf_mesh.facets.nb_vertices(fid); le++) {
        GEO::index_t nfid = sf_mesh.facets.adjacent(fid, le);
        if (nfid == GEO::NO_FACET) continue;
        if (is_skip_se_neighbor(fid, nfid)) continue;
        if (is_debug) printf("fid %d has neighbor face nfid %d\n", fid, nfid);
        k_ring_fids.insert(nfid);
        new_added_fids.insert(nfid);
      }  // for facets.nb_vertices
    }  // for k_ring_fids
    i++;
  }  // while (i < k)
  if (is_debug) {
    printf("[K_RING] fid %d found k_ring %d neighbors: %ld\n", fid_given, k,
           k_ring_fids.size());
  }
}

int store_feature_edge(const TetMesh& tet_mesh, const SurfaceMesh& sf_mesh,
                       const EdgeType& fe_type, const aint3& t2vs_group,
                       std::vector<FeatureEdge>& feature_edges) {
  int tv1 = t2vs_group[0];
  int tv2 = t2vs_group[1];
  assert(tet_mesh.tet_vs2sf_fids.find(tv1) != tet_mesh.tet_vs2sf_fids.end());
  assert(tet_mesh.tet_vs2sf_fids.find(tv2) != tet_mesh.tet_vs2sf_fids.end());

  Vector3 tv1_pos(tet_mesh.tet_vertices.at(tv1 * 3),
                  tet_mesh.tet_vertices.at(tv1 * 3 + 1),
                  tet_mesh.tet_vertices.at(tv1 * 3 + 2));
  Vector3 tv2_pos(tet_mesh.tet_vertices.at(tv2 * 3),
                  tet_mesh.tet_vertices.at(tv2 * 3 + 1),
                  tet_mesh.tet_vertices.at(tv2 * 3 + 2));

  std::vector<int> adj_sf_fids;
  set_intersection<int>(tet_mesh.tet_vs2sf_fids.at(tv1),
                        tet_mesh.tet_vs2sf_fids.at(tv2), adj_sf_fids);
  assert(adj_sf_fids.size() == 2);
  std::sort(adj_sf_fids.begin(), adj_sf_fids.end());

  FeatureEdge fe(feature_edges.size(), fe_type, t2vs_group);
  fe.t2vs_pos = {{tv1_pos, tv2_pos}};
  fe.length = GEO::length(tv2_pos - tv1_pos);
  fe.adj_sf_fs_pair = {{adj_sf_fids[0], adj_sf_fids[1]}};
  fe.adj_normals = {{get_mesh_facet_normal(sf_mesh, adj_sf_fids[0]),
                     get_mesh_facet_normal(sf_mesh, adj_sf_fids[1])}};
  fe.adj_tan_points = {{get_mesh_facet_centroid(sf_mesh, adj_sf_fids[0]),
                        get_mesh_facet_centroid(sf_mesh, adj_sf_fids[1])}};
  feature_edges.push_back(fe);
  return fe.id;
}

void print_one_fe(const int i, const std::vector<aint2>& one_fe_group) {
  printf("%d: [", i);
  for (const aint2& se : one_fe_group) {
    printf("(%d,%d), ", se[0], se[1]);
  }
  printf("]\n");
}

void print_fe_group(const std::vector<std::vector<aint2>>& fe_groups) {
  printf("[DetectFeatures] fe_groups: \n");
  FOR(i, fe_groups.size()) {
    const auto& one_fe_group = fe_groups.at(i);
    print_one_fe(i, one_fe_group);
  }
}

void store_feature_line(const TetMesh& tet_mesh, const SurfaceMesh& sf_mesh,
                        const EdgeType& fe_type,
                        const std::vector<std::vector<aint2>> fe_tet_groups,
                        std::vector<FeatureEdge>& feature_edges,
                        std::vector<FeatureLine>& feature_lines,
                        std::set<aint4>& fe_tet_info,
                        const std::set<int>& corners_se_tet,
                        std::map<int, std::set<int>>& corner2fl,
                        std::map<int, std::set<int>>& corner2fe,
                        bool is_debug) {
  auto get_common_vs = [](const aint2& e1, const aint2& e2) {
    FOR(j, 2) {
      if (e1[j] == e2[0] || e1[j] == e2[1]) return e1[j];
    }
    return -1;
  };
  if (is_debug) print_fe_group(fe_tet_groups);
  feature_lines.clear();
  for (int fl_id = 0; fl_id < fe_tet_groups.size(); fl_id++) {
    // Note: one_fe_group is pre-sorted in an order
    const std::vector<aint2>& one_fe_group = fe_tet_groups[fl_id];
    if (is_debug) print_one_fe(fl_id, one_fe_group);
    // create new feature line
    FeatureLine new_fl(fl_id, fe_type);
    new_fl.length = 0.f;
    for (int i = 0; i < one_fe_group.size(); i++) {
      // create new feature edge
      const aint2& one_fe = one_fe_group.at(i);
      if (is_debug)
        printf("processing fe_line %d: fe_edge: (%d,%d) \n", fl_id, one_fe[0],
               one_fe[1]);

      aint3 t2vs_group = {{one_fe[0], one_fe[1], fl_id}};
      int fe_id = store_feature_edge(tet_mesh, sf_mesh, fe_type, t2vs_group,
                                     feature_edges);
      new_fl.fe_ids.push_back(fe_id);
      new_fl.length += feature_edges.at(fe_id).length;

      // for corners, we store the mapping
      if (corners_se_tet.find(one_fe[0]) != corners_se_tet.end()) {
        corner2fl[one_fe[0]].insert(fl_id);
        corner2fe[one_fe[0]].insert(fe_id);
      }
      if (corners_se_tet.find(one_fe[1]) != corners_se_tet.end()) {
        corner2fl[one_fe[1]].insert(fl_id);
        corner2fe[one_fe[1]].insert(fe_id);
      }

      // update FeatureLine::t2vs_to_feid
      aint2 t2vs = {{one_fe[0], one_fe[1]}};
      std::sort(t2vs.begin(), t2vs.end());
      new_fl.t2vs_to_feid[t2vs] = fe_id;

      // update fe_tet_info
      fe_tet_info.insert({{t2vs[0], t2vs[1], fe_id, fl_id}});

      // update FeatureLine::tvs_neighbors (part1)
      FOR(j, 2) {
        int evid = one_fe[j];
        int evid_another = one_fe[(j + 1) % 2];
        if (new_fl.tvs_neighbors.find(evid) == new_fl.tvs_neighbors.end()) {
          // create new neighbors <neighbor, -1>
          new_fl.tvs_neighbors[evid] = {{evid_another, -1}};
        } else {
          // check and update order <tvid_left, tvid_right>
          aint2& neighbors = new_fl.tvs_neighbors[evid];
          if (neighbors[1] != -1 || neighbors[0] == -1) {
            printf("evid %d has neighbors [%d, %d], wrong!!\n", evid,
                   neighbors[0], neighbors[1]);
            assert(false);
          }
          neighbors = {{neighbors[0], evid_another}};
        }
      }
    }  // for one_fe_group

    // update FeatureLine::tvs_neighbors (part2)
    // swap the first order
    // update the first element of update FeatureLine::tvs_neighbors
    int first = 0, last = one_fe_group.size() - 1;
    const aint2& fe_first = one_fe_group.at(first);
    const aint2& fe_last = one_fe_group.at(last);
    int to_swap = -1;
    if (first == last || first == last - 1) {  // only 1/2 fe in fl
      // first element of first fe
      to_swap = one_fe_group.at(first)[0];
    } else {
      // find common tvid
      to_swap = get_common_vs(fe_first, fe_last);
    }
    if (is_debug) printf("found to_swap tvid: %d \n", to_swap);
    if (to_swap == -1) {
      // fl is not a loop
      bool is_swap = false;
      // swap 'first' with neighbors[1] == -1
      FOR(j, 2) {
        int evid = fe_first.at(j);
        if (new_fl.tvs_neighbors.find(evid) == new_fl.tvs_neighbors.end())
          assert(false);
        auto& neighbors = new_fl.tvs_neighbors.at(evid);
        if (neighbors[1] == -1) {
          neighbors = {{neighbors[1], neighbors[0]}};
          is_swap = true;
          break;
        }
      }
      assert(is_swap);
    } else {
      // fl is a loop
      new_fl.is_loop = true;
      if (new_fl.tvs_neighbors.find(to_swap) == new_fl.tvs_neighbors.end())
        assert(false);
      auto& neighbors = new_fl.tvs_neighbors.at(to_swap);
      neighbors = {{neighbors[1], neighbors[0]}};
    }

    // done
    if (is_debug) new_fl.print_info();
    feature_lines.push_back(new_fl);
  }
}

// helper function for sample_points_on_feature_line
// recursive sample point given len and dir
// dir: 0 or 1, matching FeatureLine::tvs_neighbors::aint2
void sample_recursive_helper(const std::vector<FeatureEdge>& feature_edges,
                             FeatureLine& fl_given, FL_Sample& sample,
                             const int dir, double len, bool is_debug) {
  // stop adding
  if (len <= 0.f) return;

  const FeatureEdge& fe_cur = feature_edges.at(sample.fe_id);
  v2int fe_endpoint_dir = fl_given.get_fe_endpoint_given_dir(fe_cur, dir);
  Vector3& fe_end_pos = fe_endpoint_dir.first;

  double len_max = GEO::length(fe_end_pos - sample.point);
  if (is_debug) {
    printf("[recursive] fl_given %d, dir %d, len %f, len_max: %f\n",
           fl_given.id, dir, len, len_max);
    printf("[recursive] fe_cur: %d (%d, %d), fe_endpoint_dir: %d\n", fe_cur.id,
           fe_cur.t2vs_group[0], fe_cur.t2vs_group[1], fe_endpoint_dir.second);
  }

  // move sample towards dir, traverse grouped feature edges
  if (len_max <= len) {
    sample.point = fe_end_pos;
    sample.fe_id = fl_given.get_next_fe_given_dir(sample.fe_id, dir);
    if (sample.fe_id == -1)  // no next fe
      return;
    if (is_debug)
      printf("[recursive] loop again, next fe_id: %d\n", sample.fe_id);
    sample_recursive_helper(feature_edges, fl_given, sample, dir, len - len_max,
                            is_debug);
  } else {
    // store next sample
    Vector3 dir_vec = GEO::normalize(fe_end_pos - sample.point);
    if (is_debug)
      printf("[recursive] sample.point: (%f,%f,%f), dir_vec: (%f,%f,%f) \n",
             sample.point[0], sample.point[1], sample.point[2], dir_vec[0],
             dir_vec[1], dir_vec[2]);
    sample.point = sample.point + len * dir_vec;
    fl_given.samples.push_back(sample);
    if (is_debug) sample.print_info();
  }
}

// Given a FeatureLine as a set of FeatureEdges,
// we start from a random point on line, then expand into two directions.
void sample_points_on_feature_line(
    const std::vector<FeatureEdge>& feature_edges, FeatureLine& fl_given,
    double esp_len, bool is_debug) {
  assert(fl_given.length > 0.f);
  // get a random point as init
  fl_given.samples.clear();
  double len_random = fl_given.length * RANDOM_01();
  double num_pins = std::floor(fl_given.length / esp_len);
  if (num_pins < 3) {
    double new_eps_len = esp_len / 3;
    if (is_debug)
      printf("[FL sample] FL %d update esp_len %f->%f \n", fl_given.id, esp_len,
             new_eps_len);
    esp_len = new_eps_len;
  }

  if (is_debug) {
    fl_given.print_info();
    printf("[FL Sample] FL %d length: %f, esp_len: %f, len_random: %f \n",
           fl_given.id, fl_given.length, esp_len, len_random);
  }

  if (is_debug) printf("---- selecting random point \n");
  auto& fe_first = feature_edges.at(fl_given.fe_ids.at(0));
  FL_Sample sample(fe_first.t2vs_pos[0], fe_first.id);  // start from first fe
  sample_recursive_helper(feature_edges, fl_given, sample, 1 /*dir*/,
                          len_random, is_debug);
  if (is_debug)
    printf("fl_given.samples 2 size: %ld\n", fl_given.samples.size());
  if (fl_given.samples.empty()) return;
  assert(fl_given.samples.size() == 1);

  // expand samples from init
  // ninwang: no need to worry about loop, since loop still stores as a vector
  if (is_debug) printf("---- expanding samples \n");
  auto sample1 = fl_given.samples.front();  // copy
  auto sample0 = sample1;                   // copy
  // dir = 1, ->
  while (sample1.fe_id != -1)
    sample_recursive_helper(feature_edges, fl_given, sample1, 1 /*dir*/,
                            esp_len, is_debug);
  if (is_debug) printf("done dir 1\n");
  // dir = 0, <-
  while (sample0.fe_id != -1)
    sample_recursive_helper(feature_edges, fl_given, sample0, 0 /*dir*/,
                            esp_len, is_debug);
  if (is_debug) printf("done dir 0\n");
  if (is_debug)
    printf("[FL sample] FL %d add %d sample points\n", fl_given.id,
           fl_given.samples.size());
}