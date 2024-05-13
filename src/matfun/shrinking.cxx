/**
 * @file shrinking.cxx
 * @author Ningna Wang (ningna.wang@utdallas.com)
 * @brief  Calculate a medial ball for a given oriented point using the
 * shrinking ball algorithm (https://3d.bk.tudelft.nl/rypeters/pdfs/16candg.pdf
 * section 3.2) Original code:
 * https://github.com/tudelft3d/masbcpp/blob/master/src/compute_ma_processing.cpp
 * @version 0.1
 * @date 2022-12-24
 *
 * @copyright Copyright (c) 2022
 *
 */
#include "shrinking.h"

#include <stdlib.h>

float compute_radius(const Vector3& p, const Vector3& n, const Vector3& q,
                     bool is_debug = false) {
  // Compute radius of the ball that touches points p and q and whose center
  // falls on the normal n from p
  // double d = GEO::Geom::distance(p, q);
  double d = (p - q).length();
  double cos_theta = GEO::dot(n, p - q) / d;

  if (is_debug)
    printf(
        "p: (%f,%f,%f), q: (%f,%f,%f), n: (%f,%f,%f), d: %f, cos_theta: %f \n",
        p[0], p[1], p[2], q[0], q[1], q[2], n[0], n[1], n[2], d, cos_theta);

  return float(d / (2. * cos_theta));
}

bool shrink_post_process(const SurfaceMesh& sf_mesh,
                         const std::vector<FeatureEdge>& feature_edges,
                         MedialSphere& msphere, bool is_del_near_ce,
                         bool is_del_near_se, bool is_debug) {
  bool is_good = true;
  // Work 1:
  // update tangent planes / tangent cc_lines
  //
  // update tangent plane of p
  if (is_debug)
    printf("msphere %d updating p, is_p_on_ce %d, p_fid %d \n", msphere.id,
           msphere.ss.is_p_on_ce(), msphere.ss.get_p_fid());
  if (msphere.is_on_ce_pre_or_fix()) {
    // at least one pin on concave line
    assert(msphere.tan_cc_lines.size() == 1);
    assert(msphere.ss.is_p_on_ce());
  } else {
    msphere.update_tan_planes_from_ss_params(sf_mesh, true /*is_update_p*/,
                                             false /*is_update_q*/);
  }
  // update tangent plane of q
  if (is_debug)
    printf("msphere %d updating q, is_q_on_ce %d, q_fid %d \n", msphere.id,
           msphere.ss.is_q_on_ce(), msphere.ss.get_q_fid());
  if (msphere.ss.is_q_on_ce()) {
    msphere.update_tan_cc_lines_from_ss_params(
        sf_mesh, feature_edges, false /*is_update_p*/, true /*is_update_q*/);
  } else {
    msphere.update_tan_planes_from_ss_params(sf_mesh, false /*is_update_p*/,
                                             true /*is_update_q*/);
  }
  // sanity check
  for (const auto& tan_pl : msphere.tan_planes) {
    if (tan_pl.fid < 0) {  // sanity check
      msphere.print_info();
      assert(false);
    }
  }
  // update MedialSphere::covered_sf_fids_in_group
  msphere.update_sphere_covered_sf_fids(sf_mesh, false /*is_debug*/);

  // Work 2:
  // if any tangent point is too close to concave lines, then delete
  // this is because we have inserted enough spheres around concave lines
  // by setting pin point p on concave lines and shrink
  if (is_del_near_ce) {
    for (const auto& tan_pl : msphere.tan_planes) {
      msphere.min_sq_dist_2cc =
          sf_mesh.aabb_wrapper.get_sq_dist_to_ce(tan_pl.tan_point);
      if (msphere.min_sq_dist_2cc <= SQ_DIST_TO_CC) {
        printf("[shrink] msphere %d too close to concave lines %f, bad \n",
               msphere.id, msphere.min_sq_dist_2cc);
        is_good = false;
        break;
      }
    }  // for msphere.tan_planes
    return is_good;
  }

  // // Work 3:
  // // if any tangent point is too close to sharp edges, then delete
  // if (is_del_near_se) {
  //   double sq_dist_2se =
  //   sf_mesh.aabb_wrapper.get_sq_dist_to_se(msphere.center); if (sq_dist_2se
  //   <= SCALAR_1) {
  //     printf("[shrink] msphere %d too close to sharp edge %f, bad \n",
  //            msphere.id, sq_dist_2se);
  //     return false;
  //   }
  // }

  // ninwang: no need, shrink will handle concave line
  //
  // // Check 2:
  // // Update sphere of type SphereType::T_2_c
  // // check if q is also around concave line
  // // if so, then remove from tangent plane and save as concave line
  // // (sphere tangent to 2 concave lines)
  // if (msphere.type == SphereType::T_2_c) {
  //   // at least one pin on concave line
  //   assert(msphere.tan_cc_lines.size() == 1);
  //   double q_sq_dist;
  //   Vector3 q_tan_point = msphere.ss.q;
  //   int feid = aabb_wrapper.project_to_ce_get_feid(q_tan_point, q_sq_dist);
  //   const FeatureEdge& one_fe = feature_edges.at(feid);
  //   if (q_sq_dist <= SQ_DIST_TO_CC) {
  //     // clear tangent plane
  //     msphere.tan_planes.clear();
  //     TangentConcaveLine new_cc_line(msphere.tan_cc_lines.size(), one_fe);
  //     new_cc_line.tan_point = q_tan_point;
  //     new_cc_line.normal = GEO::normalize(q_tan_point - msphere.center);
  //     new_cc_line.energy = q_sq_dist;
  //     new_cc_line.energy_over_sq_radius =
  //         q_sq_dist / std::pow(msphere.radius, 2);
  //     msphere.new_cc_line_no_dup(new_cc_line);
  //     if (is_debug) {
  //       printf("[shrink] update q on concave line %d, id_fe: %d\n",
  //              new_cc_line.id, new_cc_line.id_fe);
  //       msphere.print_tan_cc_lines();
  //     }
  //     // msphere.dilate_sphere_radius();
  //   }
  // }  // SphereType::T_2_c
  return is_good;
}

// p_fid / q_fid:
// positive: fid
// negative: concave eid, index matching FeatureEdge::id
bool shrink_sphere(const SurfaceMesh& sf_mesh, const AABBWrapper& aabb_wrapper,
                   const std::set<aint2>& fe_sf_fs_pairs,
                   const std::vector<FeatureEdge>& feature_edges,
                   MedialSphere& msphere, int itr_limit, bool is_del_near_ce,
                   bool is_del_near_se, bool is_debug) {
  uint iteration_limit = 30;
  bool itr_given = false;
  if (itr_limit != -1 && itr_limit > 0) {
    iteration_limit = itr_limit;
    itr_given = true;
  }
  bool is_p_on_ce = msphere.ss.is_p_on_ce();
  Vector3 p = msphere.ss.p, n = msphere.ss.p_normal;
  int p_fid = msphere.ss.p_fid;  // might be negative!!! let it be
  assert(p_fid != -1);
  Vector3 q, q_n;
  float r = msphere.radius;
  double sq_dist = -1., sq_dist_ce = -1.;
  Vector3 c_next = msphere.center;
  int q_fid_prev = -1, q_fid;
  bool is_good = true;
  unsigned int num_itr = 1;
  int k = 5;                      // k_ring sf_mesh face neighbors
  if (sf_mesh.is_non_cad) k = 2;  // use smaller k for non-cad model
  std::set<int> k_ring_fids;      // can be EMPTY if p is on ce
  if (is_p_on_ce) {
    aint2 fe_adj_sf_fs_pair =
        feature_edges.at(msphere.ss.get_p_fid()).adj_sf_fs_pair;
    get_k_ring_neighbors_no_cross(sf_mesh, fe_sf_fs_pairs, fe_adj_sf_fs_pair[0],
                                  k, k_ring_fids, false /*is_clear_cur*/,
                                  false /*is_debug*/);
    get_k_ring_neighbors_no_cross(sf_mesh, fe_sf_fs_pairs, fe_adj_sf_fs_pair[1],
                                  k, k_ring_fids, false /*is_clear_cur*/,
                                  false /*is_debug*/);
  } else {
    get_k_ring_neighbors_no_cross(sf_mesh, fe_sf_fs_pairs, p_fid, k,
                                  k_ring_fids, true /*is_clear_cur*/,
                                  false /*is_debug*/);
  }
  if (is_debug) {
    printf("-------p's p_fid %d, k_ring_fids:\n", p_fid);
    print_set<int>(k_ring_fids);
  }

  if (!all_finite(c_next)) return false;
  while (true) {
    // Stop iteration if this looks like an infinite loop:
    if (num_itr > iteration_limit) {
      if (itr_given) {
        if (is_debug)
          printf(
              "[shrink] good break, msphere %d reaches given iteration limits "
              "%d\n",
              msphere.id, iteration_limit);
        is_good = true;
      } else {
        if (is_debug)
          printf("[shrink] bad break, msphere %d reaches iteration limits %d\n",
                 msphere.id, iteration_limit);
        is_good = false;
      }
      break;
    }

    // q_fid_prev = msphere.ss.q_fid;
    bool is_q_on_ce = false;
    q_fid_prev = msphere.ss.get_q_fid();
    q_fid = aabb_wrapper.get_nearest_point_on_sf(c_next, q, sq_dist);
    q_n = get_mesh_facet_normal(sf_mesh, q_fid);
    if (is_debug)
      printf(
          "[shrink] msphere %d has is_p_on_ce: %d, p_fid: %d, found "
          "q_fid/q_feid: %d->%d\n",
          msphere.id, is_p_on_ce, p_fid, msphere.ss.q_fid, q_fid);
    // check if q closer to any concave edge
    Vector3 q_copy = q;
    int feid = aabb_wrapper.get_nearest_point_on_ce(c_next, q_copy, sq_dist_ce);
    if (is_debug)
      printf("sq_dist_ce %f to feid %d, sq_dist %f\n", sq_dist_ce, sq_dist,
             feid);
    if (feid != UNK_FACE && sq_dist_ce <= sq_dist + SCALAR_ZERO_6) {
      is_q_on_ce = true;
      sq_dist = sq_dist_ce;
      q = q_copy;
      q_fid = MedialSphere::convert_ss(feid);  // let it be negative
      q_n = GEO::normalize(q - c_next);
      if (is_debug)
        printf(
            "[shrink] msphere %d found is_q_on_ce: %d, q_feid: %d, "
            "sq_dist_ce: %f <= sq_dist: %f\n",
            msphere.id, is_q_on_ce, q_fid, sq_dist_ce, sq_dist);
    }

    // This should handle all (special) cases where we want to break the loop
    // - normal case when ball no longer shrinks
    // - the case where q==p
    // - any duplicate point cases
    // Detail see:
    // https://github.com/tudelft3d/masbcpp/blob/master/src/compute_ma_processing.cpp#L88
    double dist_pq = GEO::Geom::distance(p, q);
    // minus a small value for radius
    double sq_radius_eps = (r - SCALAR_ZERO_3) * (r - SCALAR_ZERO_3);
    // use smaller threashold for non-cad model
    if (sf_mesh.is_non_cad)
      sq_radius_eps = (r - SCALAR_ZERO_6) * (r - SCALAR_ZERO_6);
    // model scaled [0, 1000]^3, use larger threashold
    if (sq_radius_eps < SCALAR_ZERO_1) {
      if (is_debug)
        printf("[shrink] bad break, msphere %d shrink too small\n", msphere.id);
      is_good = false;
      break;
    }
    if (sq_dist >= sq_radius_eps || dist_pq < SCALAR_1 || p_fid == q_fid) {
      if (is_debug) {
        printf(
            "[shrink] good break, msphere %d found its breaking condition, "
            "sq_dist: %f, sq_radius_eps: %f, dist_pq: %f\n",
            msphere.id, sq_dist, sq_radius_eps, dist_pq);
        msphere.print_ss_info();
      }
      is_good = true;
      break;
    }
    // if both p and q on concave lines, check number of cc_line group
    if (is_p_on_ce && is_q_on_ce) {
      const auto& p_ce = feature_edges.at(msphere.ss.get_p_fid());
      const auto& q_ce = feature_edges.at(
          MedialSphere::convert_ss(q_fid));  // q_fid is negative
      if (p_ce == q_ce) {  // same as TangentConcaveLine::operator==()
        if (is_debug) {
          printf(
              "[shrink] good break, msphere %d found its breaking condition, "
              "p_ce_group: %d == q_ce_group: %d\n",
              msphere.id, p_ce.get_fl_id(), q_ce.get_fl_id());
          msphere.print_ss_info();
        }
        is_good = true;
        break;
      }
    }
    // check if q_fid within p's k_ring_fids
    // break if within k-ring and normals are not too deviated
    // k_ring_fids can be empty of p is on ce
    double n_angle = angle_between_two_vectors_in_degrees(n, q_n);
    aint2 q_fid_adj = {{-1, q_fid}};
    double angle_eps = EPS_DEGREE_30;
    if (sf_mesh.is_non_cad) angle_eps = EPS_DEGREE_10;
    if (is_q_on_ce)
      q_fid_adj =
          feature_edges.at(MedialSphere::convert_ss(q_fid)).adj_sf_fs_pair;
    if (is_debug)
      printf("q_fid %d, q_fid_adj (%d,%d)\n", q_fid, q_fid_adj[0],
             q_fid_adj[1]);
    if ((k_ring_fids.find(q_fid_adj[0]) != k_ring_fids.end() ||
         k_ring_fids.find(q_fid_adj[1]) != k_ring_fids.end()) &&
        n_angle <= angle_eps) {
      if (is_debug) {
        printf(
            "[shrink] good break, msphere %d found is_q_on_ce %d, q_fid %d, "
            "q_fid_adj (%d,%d), within %d-ring of p_fid %d, and normal angle: "
            "%f\n",
            msphere.id, is_q_on_ce, q_fid, q_fid_adj[0], q_fid_adj[1], k, p_fid,
            n_angle);
        msphere.print_ss_info();
      }
      is_good = true;
      break;
    }

    // Compute next ball center
    r = compute_radius(p, n, q, is_debug);
    c_next = p - n * r;

    if (!all_finite(c_next)) {
      if (is_debug) {
        printf(
            "[shrink] bad break, msphere %d has q_fid: %d->%d, c_next "
            "(%f,%f,%f) with radius: %f\n",
            msphere.id, msphere.ss.q_fid, q_fid, c_next[0], c_next[1],
            c_next[2], r);
        msphere.print_ss_info();
      }
      is_good = false;
      break;
    }

    if (is_debug) {
      double dist_to_p = GEO::distance(p, c_next);
      double dist_to_q = GEO::distance(q, c_next);
      printf(
          "[shrink] sphere %d updating, pos (%f,%f,%f) -> (%f,%f,%f), radius "
          "%f->%f. is_q_on_ce: %d, q_fid: %d->%d, sq_dist: %f, sq_radius_eps: "
          "%f, dist_pq: %f, dist_to_p: %f, dist_to_q: %f\n",
          msphere.id, msphere.center[0], msphere.center[1], msphere.center[2],
          c_next[0], c_next[1], c_next[2], msphere.radius, r, is_q_on_ce,
          msphere.ss.get_q_fid(), q_fid, sq_dist, sq_radius_eps, dist_pq,
          dist_to_p, dist_to_q);

      assert(std::abs(dist_to_p - r) < SCALAR_ZERO_2 &&
             std::abs(dist_to_q - r) < SCALAR_ZERO_2);
    }

    msphere.center = c_next;
    msphere.radius = r;
    msphere.ss.q = q;
    msphere.ss.q_fid = q_fid;  // q_fid could be negative
    msphere.ss.q_normal = q_n;
    num_itr++;
  }  // while true

  if (is_debug)
    printf("[shrink] Done. msphere: %d has num_itr: %d \n", msphere.id,
           num_itr);

  // post-process
  if (is_good)
    is_good = shrink_post_process(sf_mesh, feature_edges, msphere,
                                  is_del_near_ce, is_del_near_se, is_debug);
  if (is_debug) {
    printf("[shrink] shrink_post_process done\n");
    msphere.print_ss_info();
  }

  if (is_good == false || num_itr == 0) return false;
  // update type only when not given
  if (msphere.type == SphereType::T_UNK) msphere.type = SphereType::T_2;
  return true;
}

void shrink_spheres(const SurfaceMesh& sf_mesh, const AABBWrapper& aabb_wrapper,
                    const std::set<aint2>& fe_sf_fs_pairs,
                    const std::vector<FeatureEdge>& feature_edges,
                    std::vector<MedialSphere>& all_medial_spheres,
                    int itr_limit, bool is_del_near_ce, bool is_del_near_se,
                    bool is_debug) {
  if (is_debug) printf("calling shrink_spheres ... \n");
  int num_active = 0;
  for (unsigned int i = 0; i < all_medial_spheres.size(); i++) {
    MedialSphere& msphere = all_medial_spheres[i];
    bool is_success = shrink_sphere(sf_mesh, aabb_wrapper, fe_sf_fs_pairs,
                                    feature_edges, msphere, itr_limit,
                                    is_del_near_ce, is_del_near_se, is_debug);
    if (!is_success) {
      msphere.is_deleted = true;
      continue;
    }
    num_active++;
  }
  if (is_debug)
    printf("[shrink] medial sphere active: %d/%ld \n", num_active,
           all_medial_spheres.size());
}

void pre_and_init_aabb(GEO::Mesh& sf_mesh, AABBWrapper& aabb_wrapper) {
  GEO::compute_normals(sf_mesh);
  aabb_wrapper.init_sf_mesh_and_tree(sf_mesh, true /*is_reorder*/);
}

void pre_and_init_feature_aabb(const TetMesh& tet_mesh,
                               AABBWrapper& aabb_wrapper) {
  const std::vector<FeatureEdge>& feature_edges = tet_mesh.feature_edges;
  aabb_wrapper.init_feature_meshes_and_trees(tet_mesh.tet_vertices,
                                             feature_edges);
}

void init_and_shrink(const SurfaceMesh& sf_mesh, const TetMesh& tet_mesh,
                     std::vector<MedialSphere>& all_medial_spheres,
                     int num_init_spheres, int itr_limit, bool is_debug) {
  if (is_debug)
    printf("[init_and_shrink] calling init_and_shrink with num %d \n",
           num_init_spheres);
  all_medial_spheres.clear();
  std::map<int, std::set<int>> v2fs;  // vid -> fids
  for (uint f = 0; f < sf_mesh.facets.nb(); f++) {
    for (uint lv = 0; lv < 3; ++lv) {
      uint v = sf_mesh.facets.vertex(f, lv);
      v2fs[v].insert(f);
    }
  }

  // select non_fe_adj random surface facets
  std::set<int> rand_fids;
  int h = sf_mesh.facets.nb();
  for (int i = 0; i < num_init_spheres + 10; i++) {
    int rand_fidx = RANDOM_INT(0, h);
    // if (fe_adj_fs.find(rand_fidx) != fe_adj_fs.end()) continue;
    rand_fids.insert(rand_fidx);
    if (rand_fids.size() == num_init_spheres) break;
  }
  if (rand_fids.size() != num_init_spheres) {
    printf("[init_and_shrink] init %ld/%d spheres as we can\n",
           rand_fids.size(), num_init_spheres);
  }

  // shrink spheres
  for (auto& rand_fidx : rand_fids) {
    // here we use random point in a triangle as pin point (instead of centroid)
    // to avoid the symmetric cases (easily happens in model Cube)
    // https://stackoverflow.com/a/21722167
    Vector3 p = get_random_point_given_facet(
        sf_mesh.vertices.point(sf_mesh.facets.vertex(rand_fidx, 0)),
        sf_mesh.vertices.point(sf_mesh.facets.vertex(rand_fidx, 1)),
        sf_mesh.vertices.point(sf_mesh.facets.vertex(rand_fidx, 2)));
    Vector3 p_normal = get_mesh_facet_normal(sf_mesh, rand_fidx);
    MedialSphere msphere(all_medial_spheres.size(), p, p_normal);
    msphere.ss.p_fid = rand_fidx;
    msphere.itr_cnt = 0;  // init spheres
    add_new_sphere_validate(all_medial_spheres, msphere);
    if (is_debug)
      printf("[init_and_shrink] choosing rand_fidx: %d, p: (%f, %f, %f)\n",
             rand_fidx, p[0], p[1], p[2]);
  }
  if (is_debug)
    printf("[init_and_shrink] init %ld medial spheres \n",
           all_medial_spheres.size());
  shrink_spheres(sf_mesh, sf_mesh.aabb_wrapper, sf_mesh.fe_sf_fs_pairs,
                 tet_mesh.feature_edges, all_medial_spheres, itr_limit,
                 true /*is_del_near_ce*/, true /*is_del_near_se*/, is_debug);
}

// will update TetMesh::fl2corner_sphere
int init_corner_spheres(const int num_itr_global, TetMesh& tet_mesh,
                        std::vector<MedialSphere>& all_medial_spheres) {
  const auto& corners_se_tet = tet_mesh.corners_se_tet;
  const auto& corner2fl = tet_mesh.corner2fl;
  const auto& corner2fe = tet_mesh.corner2fe;
  auto& fl2corner_sphere = tet_mesh.fl2corner_sphere;
  fl2corner_sphere.clear();
  if (corners_se_tet.empty()) return 0;
  for (const auto& c_tid : corners_se_tet) {
    assert(corner2fl.find(c_tid) != corner2fl.end());
    assert(corner2fe.find(c_tid) != corner2fe.end());
    Vector3 center(tet_mesh.tet_vertices[c_tid * 3],
                   tet_mesh.tet_vertices[c_tid * 3 + 1],
                   tet_mesh.tet_vertices[c_tid * 3 + 2]);
    MedialSphere corner_sphere(all_medial_spheres.size(), center,
                               SCALAR_FEATURE_RADIUS, SphereType::T_1_N,
                               num_itr_global);
    corner_sphere.corner_fls = corner2fl.at(c_tid);
    corner_sphere.corner_fes = corner2fe.at(c_tid);
    for (const auto& fl_id : corner_sphere.corner_fls)
      fl2corner_sphere[fl_id].insert(corner_sphere.id);
    apply_perturb(corner_sphere.center);
    corner_sphere.dilate_sphere_radius();
    // assert(corner_sphere_validate(all_medial_spheres, corner_sphere));
    all_medial_spheres.push_back(corner_sphere);
  }
  printf("[init_corner] init %ld corner spheres \n", corners_se_tet.size());
  return corners_se_tet.size();
}

int create_new_spheres_for_concave_corner(
    const SurfaceMesh& sf_mesh, const std::vector<FeatureEdge>& feature_edges,
    const ConcaveCorner& cc_corner, const double cc_len_eps,
    const double cc_normal_eps, std::vector<MedialSphere>& all_medial_spheres,
    bool is_debug) {
  int corner_tvid = cc_corner.tvid;
  if (is_debug)
    printf("[CCorner Sphere] concave corner %d creating new T_2 spheres ...\n",
           corner_tvid);
  // get position of conave corner
  int num_adj_fes = cc_corner.adj_fe_ids.size();
  Vector3 pin_corner =
      feature_edges.at(cc_corner.adj_fe_ids.at(0)).t2vs_group[0] == corner_tvid
          ? feature_edges.at(cc_corner.adj_fe_ids.at(0)).t2vs_pos[0]
          : feature_edges.at(cc_corner.adj_fe_ids.at(0)).t2vs_pos[1];
  double angle_avg = 0;
  for (const auto& one_fe_id : cc_corner.adj_fe_ids) {
    const auto one_fe = feature_edges.at(one_fe_id);
    // get average angle
    const std::array<Vector3, 2>& adj_normals = one_fe.adj_normals;
    const double angle =
        angle_between_two_vectors_in_degrees(adj_normals[0], adj_normals[1]);
    angle_avg += angle;
  }
  angle_avg /= num_adj_fes;
  int num_new_spheres = std::ceil(angle_avg / cc_normal_eps * 2) + 1;

  // create new T_2_c spheres
  std::vector<Vector3> new_normals;  // size will be num_new_spheres
  std::vector<int> new_sphere_ids;
  sample_k_vectors_given_N_vectors(cc_corner.adj_fe_dir, num_new_spheres,
                                   new_normals);

  for (int i = 0; i < num_new_spheres; i++) {
    // get random adjacent fe
    int rand_idx = RANDOM_INT(0, num_adj_fes);
    const auto& one_adj_fe =
        feature_edges.at(cc_corner.adj_fe_ids.at(rand_idx));
    Vector3 pert_dir = cc_corner.adj_fe_dir.at(rand_idx);
    double pert_max = std::min(cc_corner.adj_fe_len.at(rand_idx), cc_len_eps) /
                      2.;  // smaller region
    if (is_debug)
      printf(
          "[CCorner Sphere] concave corner tvid %d choose rand_idx: %d, "
          "pert_dir: (%f,%f,%f), pert_max: %f \n",
          corner_tvid, rand_idx, pert_dir[0], pert_dir[1], pert_dir[2],
          pert_max);

    // apply perturbation
    Vector3 pin_per =
        apply_perturb_given_dir_max(pin_corner, pert_dir, pert_max);
    MedialSphere new_sphere(all_medial_spheres.size(), pin_per, new_normals[i],
                            SphereType::T_2_c);
    // if (new_sphere.id == 448)
    //   is_debug = true;
    // else
    // is_debug = false;
    new_sphere.ss.set_p_fid(one_adj_fe.id, true /*is_on_ce*/);
    new_sphere.update_tan_cc_lines_from_ss_params(
        sf_mesh, feature_edges, true /*is_update_p*/, false /*is_update_q*/);
    assert(!new_sphere.tan_cc_lines.empty());
    if (!shrink_sphere(sf_mesh, sf_mesh.aabb_wrapper, sf_mesh.fe_sf_fs_pairs,
                       feature_edges, new_sphere, -1 /*itr_limit*/,
                       false /*is_del_near_ce*/, false /*is_del_near_se*/,
                       false /*is_debug*/))
      continue;

    // check if new_sphere.ss.q is:
    // 1. one adj_fe neighbor of corner
    // 2. too close to cc_corner(new_sphere.ss.p)
    if (new_sphere.ss.is_q_on_ce()) {
      int q_eid = new_sphere.ss.get_q_fid();
      const auto& adj_fe_ids_set = to_set(cc_corner.adj_fe_ids);
      // check if q_eid on corner neighbors
      if (adj_fe_ids_set.find(q_eid) != adj_fe_ids_set.end()) {
        double dist_pq = GEO::Geom::distance(new_sphere.ss.p, new_sphere.ss.q);
        if (dist_pq < pert_max) {
          if (is_debug)
            printf(
                "[CCorner Sphere] new_sphere %d has p with p_eid %d too close "
                "to q q_eid %d, dist_pq: %f > pert_max %f \n",
                new_sphere.id, new_sphere.ss.get_p_fid(),
                new_sphere.ss.get_q_fid(), dist_pq, pert_max);
          new_sphere.is_deleted = true;
        }
      }
    }  // check new_sphere.ss.q

    if (add_new_sphere_validate(all_medial_spheres, new_sphere)) {
      new_sphere_ids.push_back(new_sphere.id);
      if (is_debug)
        printf(
            "[CCorner Sphere] created new non-pin sphere %d, pin_per: "
            "(%f,%f,%f)\n",
            new_sphere.id, pin_per[0], pin_per[1], pin_per[2]);
    }
  }

  printf(
      "[CCorner Sphere] concave corner %d created new CCorner spheres: %ld \n",
      corner_tvid, new_sphere_ids.size());
  return new_sphere_ids.size();
}

// Wrapper function for create_new_concave_sphere_given_pin
int create_new_concave_sphere_given_pin_wrapper(
    const SurfaceMesh& sf_mesh, const std::vector<FeatureEdge>& feature_edges,
    const Vector3& pin_point, const int fe_id,
    std::vector<MedialSphere>& all_medial_spheres, int sphere_type,
    bool is_debug) {
  int try_cnt = 10;
  int new_sphere_id = -1;
  while (try_cnt > 0) {
    try_cnt--;
    new_sphere_id = create_new_concave_sphere_given_pin(
        sf_mesh, feature_edges, pin_point, fe_id, all_medial_spheres,
        sphere_type, is_debug);
    if (new_sphere_id != -1) break;
  }
  return new_sphere_id;
}

// Create one new T_2_c/T_X_c spheres given a pin sample.
// For each given pin point, only sample 1 random normal.
//
// sphere_type:
//    0: SphereType::T_2_c
//    1: SphereType::T_X_c, for newly added concave spheres
int create_new_concave_sphere_given_pin(
    const SurfaceMesh& sf_mesh, const std::vector<FeatureEdge>& feature_edges,
    const Vector3& pin_point, const int fe_id,
    std::vector<MedialSphere>& all_medial_spheres, int sphere_type,
    bool is_debug) {
  // apply perturbation
  const FeatureEdge& one_fe = feature_edges.at(fe_id);
  const std::array<Vector3, 2>& adj_normals = one_fe.adj_normals;
  Vector3 new_normal =
      sample_random_vector_given_two_vectors(adj_normals[0], adj_normals[1]);
  MedialSphere new_sphere(
      all_medial_spheres.size(), pin_point, new_normal,
      sphere_type == 0 ? SphereType::T_2_c : SphereType::T_X_c, 0 /*itr_cnt*/);
  new_sphere.ss.set_p_fid(fe_id, true /*is_on_ce*/);
  new_sphere.update_tan_cc_lines_from_ss_params(
      sf_mesh, feature_edges, true /*is_update_p*/, false /*is_update_q*/);
  assert(!new_sphere.tan_cc_lines.empty());
  if (!shrink_sphere(sf_mesh, sf_mesh.aabb_wrapper, sf_mesh.fe_sf_fs_pairs,
                     feature_edges, new_sphere, -1 /*itr_limit*/,
                     false /*is_del_near_ce*/, false /*is_del_near_se*/,
                     is_debug /*is_debug*/))
    return -1;
  if (add_new_sphere_validate(all_medial_spheres, new_sphere)) {
    if (is_debug)
      printf(
          "[CC Sphere] created new non-pin sphere %d, type %d, pin_point: "
          "(%f,%f,%f)\n",
          new_sphere.id, new_sphere.type, pin_point[0], pin_point[1],
          pin_point[2]);
    return new_sphere.id;
  }
  return -1;
}

// For each given pin point, only sample 1 random normal.
int create_new_spheres_given_pin_sample(
    const SurfaceMesh& sf_mesh, const std::vector<FeatureEdge>& feature_edges,
    const FL_Sample& pin_sample, const double cc_len_eps,
    std::vector<MedialSphere>& all_medial_spheres, bool is_debug) {
  const Vector3& pin = pin_sample.point;
  const FeatureEdge& one_fe = feature_edges.at(pin_sample.fe_id);
  const std::array<Vector3, 2>& adj_normals = one_fe.adj_normals;
  const double angle =
      angle_between_two_vectors_in_degrees(adj_normals[0], adj_normals[1]);
  // int num_new_spheres = std::ceil(angle / cc_normal_eps) + 1;
  // num_new_spheres = 1;

  // Step 1: create new T_2_c spheres
  // std::vector<Vector3> new_normals;  // size will be num_new_spheres
  std::vector<int> new_sphere_ids;
  // for (int i = 0; i < num_new_spheres; i++) {
  int new_sphere_id = create_new_concave_sphere_given_pin_wrapper(
      sf_mesh, feature_edges, pin, one_fe.id, all_medial_spheres,
      0 /*SphereType::T_2_c*/, is_debug);
  if (new_sphere_id != -1) new_sphere_ids.push_back(new_sphere_id);

  if (is_debug)
    printf("[CC Sphere] cc_line %d created new cc sphere: %ld \n", one_fe.id,
           new_sphere_ids.size());
  return new_sphere_ids.size();
}

// Insert new spheres that pin at concave edges.
// We specify pin points on concave line with length cc_len_eps (user defined),
// and for each pin point, we specify the number of newly sampled normals with
// cc_normal_eps (user defined)
//
// All FeatureEdge are grouped into FeatureLine.
//
// For each FeatureLine, we start from a random position (to avoid degeneration
// with another FeatureLine), then expand in two directions with step length
// cc_len_eps to sample pin point.
//
// For each given pin point, only sample 1 random normal.
void insert_spheres_for_concave_lines_new(
    const SurfaceMesh& sf_mesh, const std::vector<ConcaveCorner>& cc_corners,
    const std::vector<FeatureEdge>& feature_edges,
    std::vector<FeatureLine>& ce_lines,
    std::vector<MedialSphere>& all_medial_spheres,
    const double cc_len_eps /*=length, scaled in [0, Parameter::scale_max]*/,
    bool is_debug) {
  // store if two end points of given one_cc_line is visited as pin
  printf("[CC Sphere] cc_len_eps: %f\n", cc_len_eps);

  // 1. for concave lines
  for (auto& one_ce_line : ce_lines) {
    sample_points_on_feature_line(feature_edges, one_ce_line, cc_len_eps,
                                  is_debug /*is_debug*/);
    // create new spheres given samples on CE
    for (const auto& pin_sample : one_ce_line.samples) {
      create_new_spheres_given_pin_sample(sf_mesh, feature_edges, pin_sample,
                                          cc_len_eps, all_medial_spheres,
                                          is_debug);
    }
    // break;
  }

  // // step 2: for concave corners
  // for (const auto& cc_corner : cc_corners) {
  //   create_new_spheres_for_concave_corner(
  //       sf_mesh, feature_edges, cc_corner, cc_len_eps, cc_normal_eps,
  //       all_medial_spheres, true /*is_debug*/);
  // }
}