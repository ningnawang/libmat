#include "fix_geo.h"

#include "fix_geo_error.h"

/////////////////////////////////////////////
// Main Functions
/////////////////////////////////////////////

// // No use
// // may update the order of MedialFace::dual_edge_endpoints
// bool is_slab_geometry_ok(const MedialFace mface, const double dist_threshold)
// {
//   assert(mface.is_valid_st);
//   assert(mface.is_sorted_de);

//   // call sort_dual_edges() and sort_dual_edges_all() ahead
//   // so that always size == 2
//   assert(mface.dual_edge_endpoints.size() == 2);

//   // check each simple_triangle (st) & dual_edge_endpoints (de) pair
//   double avg_dist = 0.f;
//   int num_surf_endpoints = 0;  // can be either 1 or 2
//   for (int i = 0; i < mface.dual_edge_endpoints.size(); i++) {
//     // this happens when only 1 endpoint on surface
//     // but we sort endpoints based on SimpleTriangle so that size always 2
//     if (mface.dual_edge_endpoints[i].second == -1) continue;
//     Vector3 st_point = get_triangle_centroid(mface.st[i].v[0],
//     mface.st[i].v[1],
//                                              mface.st[i].v[2]);
//     const Vector3 de_point = mface.dual_edge_endpoints[i].first;
//     double dist = GEO::distance(st_point, de_point);
//     avg_dist += dist;
//     num_surf_endpoints++;
//   }
//   avg_dist /= num_surf_endpoints;

//   std::vector<int> fvids;
//   for (const auto& v : mface.vertices_) fvids.push_back(v);
//   printf(
//       "[SlabCheck] checking mface %d with mspheres (%d,%d,%d) has dist %f to
//       " "surface, and %d endpoints of dual_edge on surface, dist_threshold:
//       %f "
//       "\n",
//       mface.fid, fvids[0], fvids[1], fvids[2], avg_dist, num_surf_endpoints,
//       dist_threshold);

//   if (avg_dist > dist_threshold) {
//     return false;
//   }
//   return true;
// }

// /**
//  * @brief [No use] For each medial face, check the surface from the centroid
//  of
//  * a simple triangle to its corresponding dual powercell edge endpoint on
//  * surface.
//  *
//  * @param sf_mesh
//  * @param params
//  * @param mmesh
//  * @param all_medial_spheres
//  * @param slabs_to_fix
//  * @param is_debug
//  */
// void check_mm_geo(const SurfaceMesh& sf_mesh, const Parameter& params,
//                   MedialMesh& mmesh,
//                   std::vector<MedialSphere>& all_medial_spheres,
//                   std::set<int>& slabs_to_fix, bool is_debug) {
//   // check medial mesh
//   assert(mmesh.vertices != nullptr);
//   if (mmesh.faces.empty()) return;
//   //   for (auto& mface : mmesh.faces) {
//   //     printf("mface %d has vertices: (%d,%d,%d) end_points size: %ld\n",
//   //            mface.fid, mface.vertices_[0], mface.vertices_[1],
//   //            mface.vertices_[2], mface.dual_edge_endpoints.size());
//   //   }

//   mmesh.sort_dual_edges_all(sf_mesh);
//   // slabs_to_fix stores MedialFace::fid
//   slabs_to_fix.clear();

//   assert(params.bbox_diag_l != 0.f);
//   double dist_threshold = params.bbox_diag_l * params.hd_rel_slab;
//   for (auto& mface : mmesh.faces) {
//     if (!mface.is_valid_st) continue;
//     if (!mface.is_sorted_de) continue;
//     // std::vector<int> fvids;
//     // for (const auto& v : mface.vertices_) fvids.push_back(v);
//     // printf("[MM_Check] checking mface %d with mspheres: (%d,%d,%d)\n",
//     //        mface.fid, fvids[0], fvids[1], fvids[2]);
//     if (!is_slab_geometry_ok(mface, dist_threshold)) {
//       slabs_to_fix.insert(mface.fid);
//     }
//   }
// }

// // No use
// void fix_mm_by_adding_new_sphere(const SurfaceMesh& sf_mesh,
//                                  const std::set<int>& slabs_to_fix,
//                                  const MedialMesh& mmesh,
//                                  std::vector<MedialSphere>&
//                                  all_medial_spheres, bool is_debug) {
//   printf("---------------- Fixing MedialMesh ...\n");
//   printf("[FixMM] found %ld slab to fix \n", slabs_to_fix.size());

//   // fix by adding new spheres
//   for (int mm_fid : slabs_to_fix) {
//     assert(mm_fid >= 0 && mm_fid < mmesh.faces.size());
//     auto& mface = mmesh.faces.at(mm_fid);
//     assert(mface.dual_edge_endpoints.size() > 0);
//     // select any dual_edge endpoint as pin
//     v2int v2fid_chosen;
//     v2fid_chosen.second = -1;
//     FOR(i, 2) {
//       if (mface.dual_edge_endpoints[i].second == -1) continue;
//       v2fid_chosen = mface.dual_edge_endpoints[i];
//       break;
//     }
//     assert(v2fid_chosen.second != -1);
//     printf(
//         "[FixMM] select endpoint v2fid_chosen %d as pin, sf_mesh #facets:
//         %d\n", v2fid_chosen.second, sf_mesh.facets.nb());

//     insert_new_sphere_given_v2fid(sf_mesh, v2fid_chosen, all_medial_spheres,
//                                is_debug);
//   }
// }

////////////////////////////////////////////////////////////////////////////////////////////////////
// ninwang: do not use the joint edge of 3 powercells to measure
// there might be some distance even slab envelop already on surface

/**
 * @brief For each medial face, check the surface from the centroid of
 * a simple triangle to its nearest point on surface.
 *
 * @param aabb_wrapper
 * @param mface
 * @param dist_threshold
 * @param v2fid_chosen return nearest v2fid if false
 * @return true
 * @return false
 */
bool is_slab_geometry_ok_aabb(const AABBWrapper& aabb_wrapper,
                              const MedialFace& mface,
                              const double dist_threshold,
                              v2int& v2fid_chosen) {
  assert(mface.is_valid_st);

  // find the minimal, not average
  double min_dist = DBL_MAX;
  // double avg_dist = 0.f;
  for (int i = 0; i < 2; i++) {
    // Vector3 st_point = get_triangle_centroid(mface.st[i].v[0],
    // mface.st[i].v[1],
    //                                          mface.st[i].v[2]);
    // Vector3 nearest_p;

    // int fid = -1;
    // bool is_int = aabb_wrapper.get_ray_nearest_intersection(
    //     st_point, mface.st[i].normal, nearest_p, fid);
    // if (!is_int || fid == -1) {
    //   printf(
    //       "[SlabCheck] mface %d with mspheres (%d,%d,%d) has st %d has no "
    //       "intersection to surface??? check this \n",
    //       mface.fid, mface.vertices_[0], mface.vertices_[1],
    //       mface.vertices_[2], i);
    //   continue;
    // }
    // double dist = (st_point - nearest_p).length();

    // double sq_dist;
    // int fid =
    //     aabb_wrapper.get_nearest_point_on_sf(st_point, nearest_p, sq_dist);
    // // avg_dist += std::sqrt(sq_dist);
    // double dist = std::sqrt(sq_dist);

    // if (dist < min_dist) min_dist = dist;
    // v2fid_chosen.first = nearest_p;
    // v2fid_chosen.second = fid;
    if (mface.st[i].dist_to_sf < min_dist) min_dist = mface.st[i].dist_to_sf;
    v2fid_chosen.first = mface.st[i].nearest_point;
    v2fid_chosen.second = mface.st[i].nearest_sf_fid;
  }
  // avg_dist /= 2;
  // if (avg_dist > dist_threshold) {
  if (min_dist > dist_threshold) {
    printf(
        "[SlabCheck] mface %d with mspheres (%d,%d,%d) has min_dist %f to "
        "surface > dist_threshold: %f, choose fid %d to insert new T2 sphere\n",
        mface.fid, mface.vertices_[0], mface.vertices_[1], mface.vertices_[2],
        min_dist, dist_threshold, v2fid_chosen.second);
    assert(v2fid_chosen.second != -1);
    return false;
  }

  return true;
}

// check the distance between centers of:
// 1. middle sphere of medge
// 2. the new sphere inserted pinned at the joint of two RPCs
void check_and_insert_cone_geometry(
    const SurfaceMesh& sf_mesh, const TetMesh& tet_mesh,
    const MedialEdge& medge, const double dist_threshold,
    std::vector<MedialSphere>& all_medial_spheres) {
  // step 1: find the middle sphere of medge
  const auto& msphere1 = all_medial_spheres.at(medge.vertices_[0]);
  const auto& msphere2 = all_medial_spheres.at(medge.vertices_[1]);
  Vector3 mid_center = 1. / 2. * (msphere1.center + msphere2.center);
  double mid_radius = 1. / 2. * (msphere1.radius + msphere2.radius);

  // step 2: find the new sphere inserted pinned at the joint of two RPCs
  std::vector<std::vector<v2int>> all_v2fid_chosen;  // copy
  if (msphere1.pcell.facet_cc_surf_v2fids.empty() &&
      msphere2.pcell.facet_cc_surf_v2fids.empty())
    return;
  auto itr = msphere1.pcell.facet_cc_surf_v2fids.find(msphere2.id);
  if (itr == msphere1.pcell.facet_cc_surf_v2fids.end()) {
    itr = msphere2.pcell.facet_cc_surf_v2fids.find(msphere1.id);
    if (itr == msphere2.pcell.facet_cc_surf_v2fids.end()) return;
    all_v2fid_chosen = msphere2.pcell.facet_cc_surf_v2fids.at(msphere1.id);
  } else {
    all_v2fid_chosen = msphere1.pcell.facet_cc_surf_v2fids.at(msphere2.id);
  }
  if (all_v2fid_chosen.empty() || all_v2fid_chosen[0].empty()) return;
  // find the fid closest to any sphere
  v2int v2fid_chosen =
      get_v2fid_max_to_point(sf_mesh, all_v2fid_chosen[0], msphere1.center);
  if (v2fid_chosen.second < 0) return;
  // printf("[ConeCheck] medge (%d,%d) choose fid %d to check new T2 sphere\n",
  //        msphere1.id, msphere2.id, v2fid_chosen.second);
  // sphere shrinking
  Vector3 p = v2fid_chosen.first;
  Vector3 p_normal = get_mesh_facet_normal(sf_mesh, v2fid_chosen.second);
  MedialSphere new_msphere(all_medial_spheres.size(), p, p_normal,
                           v2fid_chosen.second /*pin_fid*/);
  if (!shrink_sphere(sf_mesh, sf_mesh.aabb_wrapper, sf_mesh.fe_sf_fs_pairs,
                     tet_mesh.feature_edges, new_msphere, -1 /*itr_limit*/,
                     true /*is_del_near_ce*/, false /*is_del_near_se*/,
                     false /*is_debug*/)) {
    // printf(
    //     "[ConeCheck] medge (%d,%d) failed to check new T2 sphere with fid: "
    //     "%d\n",
    //     msphere1.id, msphere2.id, v2fid_chosen.second);
    return;
  }
  new_msphere.pcell.topo_status = Topo_Status::unkown;
  new_msphere.type == SphereType::T_2;

  // step 3: check distance
  double dist = GEO::distance(mid_center, new_msphere.center);
  if (dist > dist_threshold) {
    printf(
        "[ConeCheck] simple medge (%d,%d) has dist %f > dist_threshold: %f, "
        "add new T2 sphere with fid %d\n",
        msphere1.id, msphere2.id, dist, dist_threshold, v2fid_chosen.second);
    // step 4: insert new sphere if distance is large
    add_new_sphere_validate(all_medial_spheres, new_msphere);
  }
}

// check the distance between centers of:
// 1. middle external feature sphere of medge
// 2. middle of two RPCs
void check_and_insert_cone_extf_geometry(
    const MedialEdge& medge, const double dist_threshold,
    std::vector<MedialSphere>& all_medial_spheres) {
  // step 1: find the middle external feature sphere of medge
  const auto& msphere1 = all_medial_spheres.at(medge.vertices_[0]);
  const auto& msphere2 = all_medial_spheres.at(medge.vertices_[1]);
  if (!msphere1.is_on_extf() || !msphere2.is_on_extf()) return;
  Vector3 mid_center = 1. / 2. * (msphere1.center + msphere2.center);

  // step 2: middle of two RPCs
  if (msphere1.pcell.se_line_endpos.empty() &&
      msphere2.pcell.se_line_endpos.empty())
    return;
  for (const auto& pair : msphere1.pcell.se_line_endpos) {
    int neigh_id = pair.first[2];
    if (neigh_id != msphere2.id) continue;
    Vector3 pos = pair.second;
    Vector3 pert_dir = GEO::normalize(msphere2.center - msphere1.center);
    double pert_max = GEO::length(msphere2.center - msphere1.center) / 4.;
    MedialSphere new_msphere(all_medial_spheres.size(), pos,
                             SCALAR_FEATURE_RADIUS, SphereType::T_1_2 /*SE*/);
    new_msphere.pcell.topo_status = Topo_Status::unkown;
    new_msphere.se_line_id = pair.first[3];
    apply_perturb_given_dir_max(new_msphere.center, pert_dir, pert_max);

    // step 3: check distance
    double dist = GEO::distance(mid_center, new_msphere.center);
    if (dist > dist_threshold) {
      printf(
          "[ConeCheck] extf medge (%d,%d) has dist %f > dist_threshold: %f, "
          "add new T_1_2 sphere %d\n",
          msphere1.id, msphere2.id, dist, dist_threshold, new_msphere.id);
      // step 4: insert new sphere if distance is large
      add_new_sphere_validate(all_medial_spheres, new_msphere);
    }
  }
}

// [no use]
// check if the volume of medial sphere's RPC is too big
bool check_and_insert_sphere_geometry(
    const int num_itr_global, const SurfaceMesh& sf_mesh,
    const TetMesh& tet_mesh, const MedialSphere& msphere,
    const double vol_threshold, std::vector<MedialSphere>& all_medial_spheres,
    bool is_debug) {
  if (msphere.id == 4) is_debug = true;
  if (is_debug)
    printf("[SphereCheck] checking the volume of sphere %d\n", msphere.id);
  if (msphere.pcell.pcell_vol == -1) return true;
  if (msphere.pcell.cc_surf_v2fids.empty()) return true;

  double sphere_vol =
      4.0 / 3.0 * PI * msphere.radius * msphere.radius * msphere.radius;
  double vol_ratio = msphere.pcell.pcell_vol / sphere_vol;

  if (is_debug)
    printf(
        "[SphereCheck] msphere %d found vol_ratio %f, vol_threshold: %f, "
        "sphere_vol %f, pcell_vol: %f\n",
        msphere.id, vol_ratio, vol_threshold, sphere_vol,
        msphere.pcell.pcell_vol);

  if (vol_ratio <= vol_threshold) return true;
  for (const auto& surf_v2fids : msphere.pcell.cc_surf_v2fids) {
    if (surf_v2fids.empty()) continue;
    v2int v2fid_chosen;
    v2fid_chosen = get_v2fid_max_to_point(sf_mesh, surf_v2fids, msphere.center);
    if (is_debug)
      printf("[SphereCheck] msphere %d found fid %d, surf_v2fids size: %ld\n",
             msphere.id, v2fid_chosen.second, surf_v2fids.size());

    insert_new_sphere_given_v2fid(num_itr_global, sf_mesh, tet_mesh,
                                  v2fid_chosen, all_medial_spheres,
                                  false /*is_merge_to_ce*/, is_debug);
    break;
  }
  return false;
}

//----------------------------------------------------------------------------
// Main Function
void check_and_fix_mm_geo(const int num_itr_global, const SurfaceMesh& sf_mesh,
                          const TetMesh& tet_mesh, const Parameter& params,
                          const MedialMesh& mmesh,
                          std::vector<MedialSphere>& all_medial_spheres,
                          std::vector<double>& samples,
                          std::vector<float>& samples_dist2mat,
                          std::vector<aint3>& samples_clostprim,
                          bool is_debug) {
  // check medial mesh
  assert(mmesh.vertices != nullptr);
  // if (mmesh.faces.empty()) return;

  assert(params.bbox_diag_l != 0.f);
  double dist_threshold = params.bbox_diag_l * params.geo_error_bb_per / 100;
  const double sampling_step = params.bbox_diag_l * params.geo_rel;

  std::vector<int> sample_fids;
  sample2prims_and_compute_dist2mat(
      sf_mesh, sampling_step, all_medial_spheres, mmesh, samples, sample_fids,
      samples_dist2mat, samples_clostprim, true /*is_filter_se*/);

  int num_samples = samples.size() / 3;
  assert(samples_dist2mat.size() == num_samples);
  assert(samples_clostprim.size() == num_samples);

  v2int v2fid_chosen;
  v2fid_chosen.second = -1;
  printf("[fix_geo] dist_threshold: %f, bbox_diag_l: %f\n", dist_threshold,
         params.bbox_diag_l);
  for (int i = 0; i < num_samples; i++) {
    float dist2mat = samples_dist2mat.at(i);
    if (std::abs(dist2mat) > dist_threshold) {
      // we do not care negative dist for now
      // seems like mostly happens around concave lines
      // but will be cleaned after thinning
      // if (dist2mat > dist_threshold) {
      if (is_debug)
        printf(
            "[fix_geo] sample %d has dist2mat %f > dist_threshold %f, try to "
            "add new sphere...\n",
            i, dist2mat, dist_threshold);

      Vector3 s(samples.at(i * 3), samples.at(i * 3 + 1),
                samples.at(i * 3 + 2));
      v2fid_chosen.first = s;
      v2fid_chosen.second = sample_fids.at(i);
      insert_new_sphere_given_v2fid(
          num_itr_global, sf_mesh, tet_mesh, v2fid_chosen, all_medial_spheres,
          true /*is_merge_to_ce*/, false /*is_debug*/);
    }
  }
}