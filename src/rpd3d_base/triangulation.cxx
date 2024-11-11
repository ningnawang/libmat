#include "triangulation.h"

// no use
void generate_RT_geogram(const int& n_site, const std::vector<float>& site,
                         const std::vector<float>& site_weights) {
  assert(n_site == site.size() / 3);
  std::vector<double> site_double, site_weights_double;
  for (const auto& s : site) site_double.push_back(s);
  for (const auto& w : site_weights) site_weights_double.push_back(w);
  GEO::SmartPointer<GEO::PeriodicDelaunay3d> rt =
      new GEO::PeriodicDelaunay3d(false /*periodic*/);
  rt->set_vertices(n_site, site_double.data());
  rt->set_weights(site_weights_double.data());
  rt->compute();

  // must use GEO::vector for neighbors
  std::vector<GEO::vector<GEO::index_t>> site_neighbors(n_site);
  int max_num_neigh = -1;
  for (uint i = 0; i < n_site; i++) {
    auto& neighbors = site_neighbors[i];
    rt->get_neighbors(i, site_neighbors[i]);
    if (neighbors.size() > max_num_neigh) {
      max_num_neigh = neighbors.size();
    }
  }
}

void generate_RT_CGAL_given_spheres(
    const Parameter& params, const std::vector<float>& spheres /*4D (x,y,z,r)*/,
    const std::vector<bool>& is_sphere_deleted, RegularTriangulationNN& rt,
    bool is_debug) {
  int num_spheres = spheres.size() / 4;
  if (is_debug)
    printf("[RT spheres] generate RT for %d spheres\n", num_spheres);
  rt.clean();
  // add all medial spheres
  for (int mid = 0; mid < num_spheres; mid++) {
    if (is_sphere_deleted.at(mid)) continue;
    Point_rt p(spheres.at(mid * 4 + 0), spheres.at(mid * 4 + 1),
               spheres.at(mid * 4 + 2));
    Weight weight = std::pow(spheres.at(mid * 4 + 3), 2);
    Vertex_handle_rt vh;
    vh = rt.insert(Weighted_point(p, weight));
    if (vh == nullptr) continue;  // THIS IS IMPORTANT
    // update RT
    vh->info().all_id = mid;
  }

  // add 8 bbox, vh->info().all_id = -1
  assert(params.bb_points.size() / 3 == 8);
  for (int i = 0; i < 8; i++) {
    Point_rt p(params.bb_points[i * 3], params.bb_points[i * 3 + 1],
               params.bb_points[i * 3 + 2]);
    Weight weight = SCALAR_FEATURE_RADIUS;
    Vertex_handle_rt vh;
    vh = rt.insert(Weighted_point(p, weight));
    if (vh == nullptr) continue;  // THIS IS IMPORTANT
  }
  if (is_debug)
    printf("[RT] number_of_vertices - 8: %ld, number_of_finite_edges: %ld\n",
           rt.number_of_vertices() - 8, rt.number_of_finite_edges());
}

void generate_RT_CGAL_and_mark_valid_spheres(
    const Parameter& params, std::vector<MedialSphere>& all_medial_spheres,
    RegularTriangulationNN& rt, std::set<int>& valid_sphere_ids) {
  valid_sphere_ids.clear();
  int num_spheres = all_medial_spheres.size();
  printf("[RT] generate RT for %d spheres\n", num_spheres);
  rt.clean();
  // add all medial spheres
  for (int mid = 0; mid < num_spheres; mid++) {
    MedialSphere& msphere = all_medial_spheres.at(mid);
    // reset site_id to default
    msphere.site_id = -1;
    msphere.is_rt_valid = false;
    // do not add deleted sphere in RT
    // so later we can purge them after RT
    if (msphere.is_deleted) continue;
    // skip unkown type of sphere
    if (msphere.type == SphereType::T_UNK) continue;
    Point_rt p(msphere.center[0], msphere.center[1], msphere.center[2]);
    Weight weight = std::pow(msphere.radius, 2);
    Vertex_handle_rt vh;
    vh = rt.insert(Weighted_point(p, weight));
    if (vh == nullptr) continue;  // THIS IS IMPORTANT
    // update RT
    vh->info().all_id = msphere.id;
  }

  // add 8 bbox, vh->info().all_id = -1
  assert(params.bb_points.size() / 3 == 8);
  for (int i = 0; i < 8; i++) {
    Point_rt p(params.bb_points[i * 3], params.bb_points[i * 3 + 1],
               params.bb_points[i * 3 + 2]);
    Weight weight = SCALAR_FEATURE_RADIUS;
    Vertex_handle_rt vh;
    vh = rt.insert(Weighted_point(p, weight));
    if (vh == nullptr) continue;  // THIS IS IMPORTANT
  }
  printf("[RT] number_of_vertices - 8: %ld, number_of_finite_edges: %ld\n",
         rt.number_of_vertices() - 8, rt.number_of_finite_edges());

  // mark valid spheres
  int num_valid = 0;
  for (Finite_vertices_iterator_rt vit = rt.finite_vertices_begin();
       vit != rt.finite_vertices_end(); vit++) {
    // skip 8 bbox points
    if (vit->info().all_id == -1) continue;
    // assign valid_id
    Vertex_handle_rt vh = vit;
    int all_id = vh->info().all_id;
    assert(all_id != -1);
    MedialSphere& msphere = all_medial_spheres.at(all_id);  // copy
    msphere.is_rt_valid = true;
    valid_sphere_ids.insert(msphere.id);
    num_valid++;
  }

  // update RT validness
  // NOTE: update all, not just changed
  // if (num_valid != all_medial_spheres.size()) {
  for (auto& msphere : all_medial_spheres) {
    if (!msphere.is_rt_prev_valid && msphere.is_rt_valid) {
      // some spheres became valid (from invalid)
      msphere.rt_change_status = RtValidStates::Invalid2Valid;
    } else if (msphere.is_rt_prev_valid && !msphere.is_rt_valid) {
      // some spheres became invalid (from valid)
      msphere.rt_change_status = RtValidStates::Valid2Invalid;
    } else
      msphere.rt_change_status = RtValidStates::NoChange;
    // update msphere.is_rt_prev_valid
    // msphere.is_rt_valid will be updated while loading to RT
    msphere.is_rt_prev_valid = msphere.is_rt_valid;
  }
  // }

  printf("[RT] marked valid spheres %d/%d, rt.number_of_vertices: %ld\n",
         num_valid, num_spheres, rt.number_of_vertices());
  assert(num_valid == rt.number_of_vertices() - 8);

  // update
  // 1. MedialSphere::rt_neigh_ids_prev
  // 2. MedialSphere::rt_neigh_ids [no use]
  update_spheres_RT_neighbors(rt, all_medial_spheres);
}

/**
 * @brief Given sphere + 8 bbox, we generate RT. Since some spheres may not
 * exist in RT (bcs of weights), we purge the all_medial_spheres to those valid
 * in RT. Please call remove_duplicated_medial_spheres() ahead.
 *
 * @param params provides 8 bbox info
 * @param all_medial_spheres will be updated after RT
 * @param rt
 */
void generate_RT_CGAL_and_purge_spheres(
    const Parameter& params, std::vector<MedialSphere>& all_medial_spheres,
    RegularTriangulationNN& rt) {
  int num_spheres = all_medial_spheres.size();
  printf("[RT] generate RT for %d spheres\n", num_spheres);
  rt.clean();
  // add all medial spheres
  for (int mid = 0; mid < num_spheres; mid++) {
    const MedialSphere& msphere = all_medial_spheres.at(mid);
    // do not add deleted sphere in RT
    // so later we can purge them after RT
    if (msphere.is_deleted) continue;
    // skip unkown type of sphere
    if (msphere.type == SphereType::T_UNK) continue;
    Point_rt p(msphere.center[0], msphere.center[1], msphere.center[2]);
    Weight weight = std::pow(msphere.radius, 2);
    Vertex_handle_rt vh;
    vh = rt.insert(Weighted_point(p, weight));
    if (vh == nullptr) continue;  // THIS IS IMPORTANT
    // update RT
    vh->info().all_id = msphere.id;
  }

  // add 8 bbox, vh->info().all_id = -1
  assert(params.bb_points.size() / 3 == 8);
  for (int i = 0; i < 8; i++) {
    Point_rt p(params.bb_points[i * 3], params.bb_points[i * 3 + 1],
               params.bb_points[i * 3 + 2]);
    Weight weight = SCALAR_FEATURE_RADIUS;
    Vertex_handle_rt vh;
    vh = rt.insert(Weighted_point(p, weight));
    if (vh == nullptr) continue;  // THIS IS IMPORTANT
  }
  printf("[RT] number_of_vertices - 8: %ld, number_of_finite_edges: %ld\n",
         rt.number_of_vertices() - 8, rt.number_of_finite_edges());
  // no need to purge
  if (num_spheres == rt.number_of_vertices() - 8) return;

  // purge non-exist RT vertices (spheres)
  std::vector<MedialSphere> valid_medial_spheres;
  for (Finite_vertices_iterator_rt vit = rt.finite_vertices_begin();
       vit != rt.finite_vertices_end(); vit++) {
    // skip 8 bbox points
    if (vit->info().all_id == -1) continue;
    // start purging
    int valid_id = valid_medial_spheres.size();
    Vertex_handle_rt vh = vit;
    int all_id = vh->info().all_id;
    assert(all_id != -1);
    MedialSphere new_msphere = all_medial_spheres.at(all_id);  // copy
    vh->info().all_id = valid_id;
    new_msphere.id = valid_id;
    new_msphere.is_rt_valid = true;
    valid_medial_spheres.push_back(new_msphere);
  }
  all_medial_spheres.clear();
  all_medial_spheres = valid_medial_spheres;
  printf("[RT] purged spheres %d->%ld, rt.number_of_vertices: %ld\n",
         num_spheres, valid_medial_spheres.size(), rt.number_of_vertices());
  assert(all_medial_spheres.size() == rt.number_of_vertices() - 8);

  // update RT validness
  // NOTE: update all, not just changed
  // if (num_valid != all_medial_spheres.size()) {
  for (auto& msphere : all_medial_spheres) {
    if (!msphere.is_rt_prev_valid && msphere.is_rt_valid) {
      // some spheres became valid (from invalid)
      msphere.rt_change_status = RtValidStates::Invalid2Valid;
    } else if (msphere.is_rt_prev_valid && !msphere.is_rt_valid) {
      // some spheres became invalid (from valid)
      msphere.rt_change_status = RtValidStates::Valid2Invalid;
    } else
      msphere.rt_change_status = RtValidStates::NoChange;
    // update msphere.is_rt_prev_valid
    // msphere.is_rt_valid will be updated while loading to RT
    msphere.is_rt_prev_valid = msphere.is_rt_valid;
  }
  // }
}

// -------------------------------------------------------------------------------------------
//  helper function
void convert_to_flat_site_knn(
    const int& n_site, const int& num_neigh_max,
    const std::map<int, std::set<int>>& sphere_neighbors,
    std::vector<int>& site_knn) {
  // convert to flat 2D matrix
  // 1. dim: (num_neigh_max+1) x n_site
  // 2. each column j store all neighbors of sphere all_medial_spheres.at(j)
  // 3. init as -1
  site_knn.clear();
  site_knn.resize((num_neigh_max + 1) * n_site, -1);
  for (int site_id = 0; site_id < n_site; site_id++) {
    if (sphere_neighbors.find(site_id) == sphere_neighbors.end()) {
      // not find, non-valie spheres
      // keep -1
      continue;
    }
    const std::set<int>& neighbors = sphere_neighbors.at(site_id);
    for (int lneigh_id = 0; lneigh_id < neighbors.size(); lneigh_id++) {
      site_knn[lneigh_id * n_site + site_id] =
          *std::next(neighbors.begin(), lneigh_id);
    }
  }

  // printf("site_knn matrix: \n\t");
  // for (uint sid = 0; sid < n_site; sid++) {     // column
  //   for (uint i = 0; i < num_neigh_max; i++) {  // row
  //     printf("%d ", site_knn[i * n_site + sid]);
  //   }
  //   printf("\n\t ");
  // }
  // printf("\n");
}

// return num_neigh_max
// site_knn size:         (num_neigh_max+1) x n_site
// is_sphere_valid size:  num_spheres
int get_RT_spheres_and_neighbors(const int num_spheres,
                                 const RegularTriangulationNN& rt,
                                 std::vector<int>& site_knn,
                                 std::vector<int>& is_sphere_valid,
                                 bool is_debug) {
  is_sphere_valid.clear();
  is_sphere_valid.resize(num_spheres, 0);
  int num_neigh_max = 0;
  std::map<int, std::set<int>> sphere_neighbors;
  std::vector<Vertex_handle_rt> one_neighs;
  if (is_debug) printf("[RT Neighbors] loading num_spheres %d \n", num_spheres);
  for (Finite_vertices_iterator_rt vit = rt.finite_vertices_begin();
       vit != rt.finite_vertices_end(); ++vit) {
    int all_id = vit->info().all_id;
    // skip 8 bbox points
    if (all_id == -1) continue;
    // mark valid spheres
    is_sphere_valid[all_id] = 1;
    // get neighbors
    one_neighs.clear();
    rt.finite_adjacent_vertices(vit, std::back_inserter(one_neighs));
    if (one_neighs.size() > num_neigh_max) {
      num_neigh_max = one_neighs.size();
    }

    if (is_debug)
      printf("[RT Neighbors] given sphere %d has neighbors: [", all_id);
    // save neighbor all_ids
    for (auto& neigh_handle : one_neighs) {
      int neigh_all_id = neigh_handle->info().all_id;
      // skip 8 bbox points
      if (neigh_all_id == -1) continue;
      sphere_neighbors[all_id].insert(neigh_all_id);
      sphere_neighbors[neigh_all_id].insert(all_id);  // both sides
      if (is_debug) printf("%d, ", neigh_all_id);
    }
    if (is_debug) printf("]\n");
  }

  // convert to flat 2D matrix
  // update site_knn
  convert_to_flat_site_knn(num_spheres, num_neigh_max, sphere_neighbors,
                           site_knn);

  if (is_debug)
    printf("[RT Neigh] maximum %d neighbors per sphere\n", num_neigh_max);
  return num_neigh_max;
}

// -------------------------------------------------------------------------------------------
// helper function
void convert_to_flat_site_knn(
    const int& n_site, const int& num_neigh_max,
    const std::vector<std::set<int>>& sphere_neighbors,
    std::vector<int>& site_knn) {
  // convert to flat 2D matrix
  // 1. dim: (num_neigh_max+1) x n_site
  // 2. each column j store all neighbors of sphere all_medial_spheres.at(j)
  // 3. init as -1
  site_knn.clear();
  site_knn.resize((num_neigh_max + 1) * n_site, -1);
  for (int site_id = 0; site_id < sphere_neighbors.size(); site_id++) {
    const auto& neighbors = sphere_neighbors.at(site_id);
    for (int neigh_id = 0; neigh_id < neighbors.size(); neigh_id++) {
      // not true!!!
      // assert(neigh_id <= num_neigh_max);
      site_knn[neigh_id * n_site + site_id] =
          *std::next(neighbors.begin(), neigh_id);
    }
  }

  // printf("site_knn matrix: \n\t");
  // for (uint sid = 0; sid < n_site; sid++) {     // column
  //   for (uint i = 0; i < num_neigh_max; i++) {  // row
  //     printf("%d ", site_knn[i * n_site + sid]);
  //   }
  //   printf("\n\t ");
  // }
  // printf("\n");
}

// helper function of get_RT_neighbors_given_spheres()
int update_spheres_RT_neighbors(const RegularTriangulationNN& rt,
                                std::vector<MedialSphere>& all_medial_spheres) {
  // get sphere neighbors
  // sphere -> a set of neighbors
  std::map<int, std::set<int>> sphere_neighbors;
  std::vector<Vertex_handle_rt> one_neighs;
  int num_neigh_max = 0;
  for (Finite_vertices_iterator_rt vit = rt.finite_vertices_begin();
       vit != rt.finite_vertices_end(); ++vit) {
    // skip 8 bbox points
    if (vit->info().all_id == -1) continue;
    one_neighs.clear();
    rt.finite_adjacent_vertices(vit, std::back_inserter(one_neighs));
    if (one_neighs.size() > num_neigh_max) {
      num_neigh_max = one_neighs.size();
    }

    int all_id = vit->info().all_id;
    for (auto& neigh_handle : one_neighs) {
      // skip 8 bbox points
      if (neigh_handle->info().all_id == -1) continue;
      sphere_neighbors[all_id].insert(neigh_handle->info().all_id);
    }
  }

  for (auto& msphere : all_medial_spheres) {
    msphere.rt_neigh_ids_prev = msphere.rt_neigh_ids;
    msphere.rt_neigh_ids.clear();
    if (sphere_neighbors.find(msphere.id) == sphere_neighbors.end()) continue;
    msphere.rt_neigh_ids = sphere_neighbors.at(msphere.id);
  }

  printf("maximum %d neighbors per sphere\n", num_neigh_max);
  return num_neigh_max;
}

// All sphere ids mathicng MedialSphere::id
// given spheres in given_sphere_ids
// return:
// 1. spheres_and_1rings: sphere + 1ring neighbors
// 2. sphere_neighbors: keys of spheres_N + 1ring + 2ring neighbors
int get_RT_neighbors_given_spheres(
    const RegularTriangulationNN& rt, const std::set<int>& given_sphere_ids,
    std::set<int>& spheres_and_1rings,
    std::map<int, std::set<int>>& sphere_neighbors, bool is_debug) {
  int num_neigh_max = 0;
  // sphere_neighbors.clear(); // ninwang: on purpose
  std::vector<Vertex_handle_rt> one_neighs;
  for (Finite_vertices_iterator_rt vit = rt.finite_vertices_begin();
       vit != rt.finite_vertices_end(); ++vit) {
    int all_id = vit->info().all_id;
    // skip 8 bbox points
    if (all_id == -1) continue;
    // ninwang: we store both side, so check them all
    // // skip if already stored
    // if (sphere_neighbors.find(all_id) != sphere_neighbors.end()) continue;
    // skip if not given spheres
    if (given_sphere_ids.find(all_id) == given_sphere_ids.end()) continue;
    spheres_and_1rings.insert(all_id);

    // get neighbors
    one_neighs.clear();
    rt.finite_adjacent_vertices(vit, std::back_inserter(one_neighs));
    if (one_neighs.size() > num_neigh_max) {
      num_neigh_max = one_neighs.size();
    }

    if (is_debug) printf("given sphere %d has neighbors: [", all_id);

    // save neighbor all_ids
    for (auto& neigh_handle : one_neighs) {
      int neigh_all_id = neigh_handle->info().all_id;
      // skip 8 bbox points
      if (neigh_all_id == -1) continue;
      sphere_neighbors[all_id].insert(neigh_all_id);
      sphere_neighbors[neigh_all_id].insert(all_id);  // both sides
      spheres_and_1rings.insert(neigh_all_id);
      if (is_debug) printf("%d, ", neigh_all_id);
    }

    if (is_debug) printf("]\n");
  }
  return num_neigh_max;
}

// will update map_site2msphere and MedialSphere::site_id
// return num_neigh_max, can be -1 if no new sphere found!!!
int get_RT_partial_spheres_and_neighbors(
    const std::set<int>& sphere_ids, const RegularTriangulationNN& rt,
    std::vector<MedialSphere> all_medial_spheres,
    std::vector<int>& map_site2msphere, std::set<int>& spheres_and_1rings,
    std::vector<int>& site_knn, bool is_debug) {
  auto add_to_site_and_update = [&](MedialSphere& msphere,
                                    bool is_debug = false) {
    // will reset in generate_RT_CGAL_and_mark_valid_spheres()
    if (msphere.site_id != -1) return msphere.site_id;
    // add to map (site id to MedialSphere::id))
    msphere.site_id = map_site2msphere.size();
    map_site2msphere.push_back(msphere.id);
    // if (is_debug || msphere.id == 377 || msphere.id == 649 ||
    //     msphere.id == 706) {
    //   printf("[RT Neigh] msphere %d changed to site_id: %d \n", msphere.id,
    //          msphere.site_id);
    // }
    return msphere.site_id;
  };

  if (is_debug)
    printf("[Partial RT Neighbors] sphere_ids %ld to get parital neighbors\n",
           sphere_ids.size());
  if (sphere_ids.empty()) return -1;

  // update map_site2msphere:
  // (site id to MedialSphere::id))
  // stores set of spheres_N + 1ring + 2ring neighbors,
  // use for calculating RPD later
  map_site2msphere.clear();
  spheres_and_1rings.clear();

  /////////////////////////////////
  // Get partial spheres_N + 1ring + 2ring neighbors
  std::set<int> _;
  // store spheres_N + 1ring + 2ring, both side neighboring
  std::map<int, std::set<int>> map_spheres_and_12rings;

  get_RT_neighbors_given_spheres(rt, sphere_ids, spheres_and_1rings,
                                 map_spheres_and_12rings, false /*is_debug*/);
  if (is_debug) {
    for (const auto& key_val : map_spheres_and_12rings) {
      printf("sphere %d has all_id 1-ring neighbors: [", key_val.first);
      for (const auto& v : key_val.second) {
        printf("%d, ", v);
      }
      printf("]\n");
    }
    printf("----------------------------------------------\n");
  }

  /////////////////////////////////
  // Part 2: get neighbors of map_spheres_and_12rings (2-ring)
  get_RT_neighbors_given_spheres(rt, spheres_and_1rings, _ /*no need*/,
                                 map_spheres_and_12rings, false /*is_debug*/);
  if (is_debug) {
    for (const auto& key_val : map_spheres_and_12rings) {
      printf("sphere %d has all_id 2-ring neighbors: [", key_val.first);
      for (const auto& v : key_val.second) {
        printf("%d, ", v);
      }
      printf("]\n");
    }
    printf("----------------------------------------------\n");
  }

  /////////////////////////////////
  // 1. store map_site2msphere, for spheres_N + 1ring + 2ring neighbors
  // 2. store site_neighbors (for site_knn later),
  // 3. get num_neigh_max
  int num_neigh_max = 0;
  std::vector<std::set<int>> site_neighbors(map_spheres_and_12rings.size());
  for (const auto& kv : map_spheres_and_12rings) {
    int all_id = kv.first;
    int site_id = add_to_site_and_update(all_medial_spheres.at(all_id));
    if (num_neigh_max < kv.second.size()) num_neigh_max = kv.second.size();
    for (const auto& n_all_id : kv.second) {
      int n_site_id = add_to_site_and_update(all_medial_spheres.at(n_all_id));
      // get_RT_neighbors_given_spheres() take cares of two-sided neighbors
      site_neighbors[site_id].insert(n_site_id);
      site_neighbors[n_site_id].insert(site_id);
    }
  }
  printf(
      "[RT Neigh] site_neighbors size %ld, map_site2msphere: %ld, "
      "map_spheres_and_12rings: %ld \n",
      site_neighbors.size(), map_site2msphere.size(),
      map_spheres_and_12rings.size());
  // all stores sphere + 1&2-ring neighbors
  assert(map_site2msphere.size() == map_spheres_and_12rings.size());

  if (is_debug) {
    for (int i = 0; i < site_neighbors.size(); i++) {
      printf("site %d has site_id neighbors: [", i);
      for (const auto& v : site_neighbors.at(i)) {
        printf("%d, ", v);
      }
      printf("]\n");
    }
  }

  // convert to flat 2D matrix
  int n_site = map_site2msphere.size();
  convert_to_flat_site_knn(n_site, num_neigh_max, site_neighbors, site_knn);

  printf("[RT Neigh] maximum %d neighbors per sphere\n", num_neigh_max);
  return num_neigh_max;
}

/**
 * @brief
 *
 * @param rt
 * @param n_site
 * @param site_knn
 * @return int return the site_k, maximum size of sphere neighbors
 */
int get_RT_vertex_neighbors(const RegularTriangulationNN& rt, const int& n_site,
                            std::vector<int>& site_knn) {
  // n_site == all_medial_spheres.size();
  // get sphere neighbors
  // sphere -> a set of neighbors
  std::vector<std::set<int>> sphere_neighbors(n_site);
  std::vector<Vertex_handle_rt> one_neighs;
  int num_neigh_max = 0;
  for (Finite_vertices_iterator_rt vit = rt.finite_vertices_begin();
       vit != rt.finite_vertices_end(); ++vit) {
    // skip 8 bbox points
    if (vit->info().all_id == -1) continue;
    one_neighs.clear();
    rt.finite_adjacent_vertices(vit, std::back_inserter(one_neighs));
    if (one_neighs.size() > num_neigh_max) {
      num_neigh_max = one_neighs.size();
    }

    int all_id = vit->info().all_id;
    for (auto& neigh_handle : one_neighs) {
      // skip 8 bbox points
      if (neigh_handle->info().all_id == -1) continue;
      sphere_neighbors[all_id].insert(neigh_handle->info().all_id);
    }
  }

  // convert to flat 2D matrix
  // 1. dim: (num_neigh_max+1) x n_site
  // 2. each column j store all neighbors of sphere all_medial_spheres.at(j)
  // 3. init as -1
  convert_to_flat_site_knn(n_site, num_neigh_max, sphere_neighbors, site_knn);

  // // for debug
  // printf("site_knn matrix: \n\t ");
  // for (uint i = 0; i < (num_neigh_max + 1); i++) {
  //   printf("\n\t ");
  //   // for (uint seed = 0; seed < n_site; seed++) {
  //   for (uint seed = 740; seed < 741; seed++) {
  //     printf("%d ", site_knn[i * n_site + seed]);
  //   }
  // }
  // printf("\n");

  printf("maximum %d neighbors per sphere\n", num_neigh_max);
  return num_neigh_max;
}

// Each (non-restricted) powercell is dual to a RT tet
void get_PC_vertices(const RegularTriangulationNN& rt,
                     std::vector<Vector3>& pc_vertices) {
  pc_vertices.clear();
  for (Finite_cells_iterator_rt fci = rt.finite_cells_begin();
       fci != rt.finite_cells_end(); fci++) {
    RegularTriangulationNN::Bare_point bp =
        rt.dual(fci);  // dunno why we do not need *fci here
    pc_vertices.push_back(Vector3(bp[0], bp[1], bp[2]));
  }
  printf("found pc_vertices: %ld\n", pc_vertices.size());
}