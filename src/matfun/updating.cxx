#include "updating.h"

/////////////////////////////////////////////////////////////////////////////////////
// Helper Functions
/////////////////////////////////////////////////////////////////////////////////////
bool try_add_tangent_plane(const int num_itr, MedialSphere& mat_p,
                           TangentPlane& tan_pl_q, bool is_debug) {
  // 1. check if covered by existing concave lines
  //    if so then do not add
  for (const auto& one_cc_line : mat_p.tan_cc_lines) {
    one_cc_line.purge_one_tan_plane(tan_pl_q);
    if (tan_pl_q.is_deleted) {
      if (is_debug)
        printf(
            "[Add_TanPL] mat_p %d tangent plane fid %d is covered by "
            "tan_cc_line with id_fe %d, id_fl %d, NOT add \n",
            mat_p.id, tan_pl_q.fid, one_cc_line.id_fe, one_cc_line.id_fl);
      return false;
    }
  }
  // 2. check if covered by existing tangent planes
  //    if so then delete existing one
  if (!tan_pl_q.is_deleted) {
    std::vector<TangentPlane> new_tan_planes;
    new_tan_planes.push_back(tan_pl_q);
    for (auto& tan_pl : mat_p.tan_planes) {
      if (tan_pl.is_same_normal(tan_pl_q.normal)) {
        tan_pl.is_deleted = true;
        if (is_debug) {
          printf("[Add_TanPL] mat_p %d num_itr: %d remove tangent plane: \n",
                 mat_p.id, num_itr);
          tan_pl.print_info();
        }
      } else {
        new_tan_planes.push_back(tan_pl);
      }
    }
    if (is_debug) {
      printf("[Add_TanPL] mat_p %d num_itr %d added new tangent plane: \n",
             mat_p.id, num_itr);
      tan_pl_q.print_info();
    }
    mat_p.tan_planes = new_tan_planes;
  }
  return true;
}

bool try_add_tangent_cc_line(const int num_itr, MedialSphere& mat_p,
                             TangentConcaveLine& tan_cc_line, bool is_debug) {
  // 1. check if already exits
  for (auto& exit_cc_line : mat_p.tan_cc_lines) {
    if (exit_cc_line == tan_cc_line) {
      // NOTE: do not update!!!
      //       in case the sphere is sliding along the whole concave_line
      // exit_cc_line = tan_cc_line;  // may update TangentConcaveLine::id_fl
      // etc
      if (is_debug)
        printf(
            "[Add_TanCCLine] mat_p %d tan_cc_line id_fl %d already exists, NOT "
            "add but update \n",
            mat_p.id, exit_cc_line.id_fl);
      return false;
    }
  }
  // 2. add and purge existing tangent planes
  tan_cc_line.purge_tan_planes(mat_p.tan_planes, is_debug /*is_debug*/);
  mat_p.remove_deleted_tangents(false /*is_run_cc*/);
  mat_p.tan_cc_lines.push_back(tan_cc_line);

  return true;
}

// return:
// true - break the loop
// false - keep looping
// is_good will tell us if this update is valid
//
// load tangent planes (pN, nN)
// load concave lines (c_pN, c_nN)
bool is_break_iteration(const int iteration_limit, const int num_itr,
                        const double alpha1, const double alpha2,
                        const double alpha3, const double break_threshold,
                        const bool is_break_use_energy_over_sq_radius,
                        const std::vector<FeatureEdge>& feature_edges,
                        const SurfaceMesh& sf_mesh,
                        const AABBWrapper& aabb_wrapper, MedialSphere& mat_p,
                        double& sum_energy_over_sq_radius,
                        std::vector<Vector3>& pN, std::vector<Vector3>& nN,
                        std::vector<Vector3>& c_pN, std::vector<Vector3>& c_nN,
                        bool& is_good, bool is_check_new_tan_plane,
                        bool is_debug) {
  pN.clear();
  nN.clear();
  c_pN.clear();
  c_nN.clear();
  //////////////////////////
  // Break conditions: 2 parts
  // check breaking condition
  if (is_check_new_tan_plane) {
    // fetch the nearest point q and qnormal
    bool is_q_on_ce = false;
    Vector3 q;
    double sq_dist = -1., sq_dist_ce = -1.;
    int q_fid = aabb_wrapper.get_nearest_point_on_sf(mat_p.center, q, sq_dist);
    Vector3 qnormal = get_mesh_facet_normal(sf_mesh, q_fid);
    if (is_debug) {
      printf(
          "[Iterate Break] mat_p %d num_itr %d has q_fid: %d with q: "
          "(%f,%f,%f), qnormal: (%f,%f,%f) \n",
          mat_p.id, num_itr, q_fid, q[0], q[1], q[2], qnormal[0], qnormal[1],
          qnormal[2]);
    }

    // check if q closer to any concave edge
    Vector3 q_copy = q;
    int feid =
        aabb_wrapper.get_nearest_point_on_ce(mat_p.center, q_copy, sq_dist_ce);
    if (is_debug)
      printf("[Iterate Break] sq_dist_ce %f, sq_dist %f\n", sq_dist_ce,
             sq_dist);
    if (feid != UNK_FACE && sq_dist_ce <= sq_dist + SCALAR_ZERO_6) {
      is_q_on_ce = true;
      sq_dist = sq_dist_ce;
      q = q_copy;
      q_fid = feid;  // let it be positive
      qnormal = GEO::normalize(q - mat_p.center);
      if (is_debug)
        printf(
            "[Iterate Break] msphere %d found is_q_on_ce: %d, q_feid: %d, "
            "sq_dist_ce: %f <= sq_dist: %f\n",
            mat_p.id, is_q_on_ce, q_fid, sq_dist_ce, sq_dist);
    }

    // add new tangent info
    if (is_q_on_ce) {
      assert(q_fid > -1 && q_fid < feature_edges.size());
      // create a tangent cc_line
      TangentConcaveLine tan_cc_line(sf_mesh, mat_p.tan_cc_lines.size(),
                                     feature_edges.at(q_fid));
      bool is_add = try_add_tangent_cc_line(num_itr, mat_p, tan_cc_line,
                                            is_debug /*is_debug*/);
      if (is_debug) {
        printf("[Iterate Break] is_add %d new tan_cc_line: \n", is_add);
        if (is_add) tan_cc_line.print_info();
      }
    } else {
      // create a new tangent plane
      TangentPlane tan_pl_q(sf_mesh, qnormal, q, q_fid);
      bool is_add = try_add_tangent_plane(num_itr, mat_p, tan_pl_q,
                                          is_debug /*is_debug*/);
      if (is_debug) {
        printf("[Iterate Break] is_add %d new tan_pl: \n", is_add);
        if (is_add) tan_pl_q.print_info();
      }
    }
  }  // if is_check_new_tan_plane

  //////////////////////////
  // Break condition part 1
  // iterate ends when the energy reaches minima
  double prev_energy_over_sq_radius = sum_energy_over_sq_radius;
  double sum_energy =
      mat_p.update_all_energy_values(alpha1, alpha2, alpha3, is_debug);
  sum_energy_over_sq_radius = sum_energy / std::pow(mat_p.radius, 2);

  if (is_debug) {
    printf(
        "[Iterate Break] mat_p %d num_itr: %d tangent infos before checking "
        "break conditions: \n",
        mat_p.id, num_itr);
    mat_p.print_tan_planes();
    mat_p.print_tan_cc_lines();
  }

  if (is_debug) {
    printf(
        "[Iterate Break] mat_p %d num_itr: %d, has sum_energy %f, "
        "sum_energy_over_sq_radius: %f \n",
        mat_p.id, num_itr, sum_energy, sum_energy_over_sq_radius);
  }
  // this minima is small, good
  if (is_break_use_energy_over_sq_radius &&
      (sum_energy_over_sq_radius < break_threshold)) {
    if (is_debug)
      printf(
          "[Iterate Break] sum_energy_over_sq_radius %f < break_threshold "
          "%f, minima is really small, good \n",
          sum_energy_over_sq_radius, break_threshold);
    is_good = true;
    return true;  // break
  } else if (sum_energy < break_threshold) {
    if (is_debug)
      printf(
          "[Iterate Break] sum_energy %f < break_threshold %f, minima is "
          "really small, good\n",
          sum_energy, break_threshold);
    is_good = true;
    return true;  // break
  }

  // // change is too small, not good
  // if (prev_energy_over_sq_radius != -1 &&
  //     std::abs(sum_energy_over_sq_radius - prev_energy_over_sq_radius) <=
  //         SCALAR_ZERO_5) {
  //   if (is_debug)
  //     printf(
  //         "[Iterate Break] sum_energy_over_sq_radius %f == "
  //         "prev_energy_over_sq_radius %f, break_threshold: %f,"
  //         "change too small, not good \n",
  //         sum_energy_over_sq_radius, prev_energy_over_sq_radius,
  //         break_threshold);
  //   is_good = false;
  //   return true;  // break
  // }

  //////////////////////////
  // Break condition part 2
  if (num_itr > iteration_limit || mat_p.radius < SCALAR_ZERO_3) {
    if (is_debug)
      printf(
          "[Iterate Break] break, mat_p %d reaches iteration limits %d or "
          "radius too small with sum energy: %f, mat_p.radius %f \n",
          mat_p.id, num_itr, sum_energy, mat_p.radius);
    is_good = false;
    return true;  // break
  }

  // For all iterations
  // prepare for next iteration
  for (const auto& tan_pl : mat_p.tan_planes) {
    pN.push_back(tan_pl.tan_point);
    nN.push_back(tan_pl.normal);
    if (is_debug)
      printf(
          "[Iterate Break] mat_p %d num_itr: %d, has tan_pl fid %d energy %f, "
          "reiterate \n",
          mat_p.id, num_itr, tan_pl.fid, tan_pl.energy);
  }  // for mat_p.tan_planes
  for (const auto& one_cc_line : mat_p.tan_cc_lines) {
    c_pN.push_back(one_cc_line.tan_point);
    c_nN.push_back(one_cc_line.normal);
    if (is_debug) {
      printf(
          "[Iterate Break] mat_p %d num_itr: %d, has cc_line with energy %f, "
          "reiterate \n",
          mat_p.id, num_itr, one_cc_line.energy);
    }
  }  // mat_p.tan_cc_lines

  is_good = false;
  return false;  // not break
}

// This energy function has 3 parts with weights:
// 1. distance to tangent point -> alpha 1 (pN, nN)
// 2. distance to tangent plane -> alpha 2 (pN, nN)
// 3. distance to concave line -> alpha 3 (c_pN, c_nN)
bool update_sphere_given_plN_opt(
    const double alpha1, const double alpha2, const std::vector<Vector3>& pN,
    const std::vector<Vector3>& nN, const double alpha3,
    const std::vector<Vector3>& c_pN, const std::vector<Vector3>& c_nN,
    Vector3& new_center, double& new_radius, bool is_debug) {
  if (pN.size() != nN.size() || c_pN.size() != c_nN.size()) {
    fprintf(stderr,
            "[TAN_PLN] unmatched plane parameters, pN: %ld, nN: %ld, c_pN: "
            "%ld, c_nN: %ld\n ",
            pN.size(), nN.size(), c_pN.size(), c_nN.size());
    return false;
  }

  if (is_debug) {
    printf("calling update_sphere_given_plN_opt ...\n");
    printf("[TAN_PLN] with N: %ld, alpha1: %f, alpha2: %f, alpha3: %f\n",
           pN.size(), alpha1, alpha2, alpha3);
  }
  Matrix3 A, I;
  Vector3 B, C, E;  // zeros
  A.load_zero();
  I.load_identity();
  double D = 0., F = 0.;
  // if any tangent plane
  for (int i = 0; i < pN.size(); i++) {
    const Vector3& p = pN[i];
    const Vector3& n = nN[i];
    Matrix3 common = I * alpha1 + vec_vec_trans(n, n) * alpha2;
    A += common;
    B += (alpha1 + alpha2) * n;
    C += (alpha1 + alpha2) * n;
    D += alpha1 + alpha2;
    E += GEO::mult(common, p);
    F += (alpha1 + alpha2) * GEO::dot(n, p);
  }
  // if any concave line
  for (int i = 0; i < c_pN.size(); i++) {
    const auto& p = c_pN[i];
    const auto& n = c_nN[i];
    A += I * alpha3;
    B += alpha3 * n;
    C += alpha3 * n;
    D += alpha3;
    E += alpha3 * p;
    F += alpha3 * GEO::dot(n, p);
  }

  std::array<double, 16> N_data = {{A(0, 0), A(0, 1), A(0, 2), B[0], A(1, 0),
                                    A(1, 1), A(1, 2), B[1], A(2, 0), A(2, 1),
                                    A(2, 2), B[2], C[0], C[1], C[2], D}};
  Matrix4 N(N_data.data());
  Vector4 b(E[0], E[1], E[2], F);

  Matrix4 N_inv;
  bool invertible = N.compute_inverse(N_inv);
  if (!invertible) {
    printf("[TAN_PLN] Matrix 4x4 not invertible \n");
    // printf(
    //     "[TAN_PLN] Matrix 4x4: \n"
    //     "({},{},{},{})\n({},{},{},{})\n({},{},{},{})\n({},{},{},{})\n",
    //     N(0, 0), N(0, 1), N(0, 2), N(0, 3), N(1, 0), N(1, 1), N(1, 2),
    // N(1,
    //     3), N(2, 0), N(2, 1), N(2, 2), N(2, 3), N(3, 0), N(3, 1), N(3, 2),
    //     N(3, 3));
    return false;
  }

  // update mat_p
  Vector4 c_r = GEO::mult(N_inv, b);
  if (is_debug) {
    printf("[TAN_PLN] mat new sphere: center (%f,%f,%f), radius: %f \n", c_r[0],
           c_r[1], c_r[2], c_r[3]);
  }
  if (c_r[3] < SCALAR_ZERO_3) return false;
  new_center = Vector3(c_r[0], c_r[1], c_r[2]);
  new_radius = c_r[3];
  return true;
}

void get_faces_params(const GEO::Mesh& sf_mesh, std::set<int>& k_ring_fids,
                      std::map<int, std::array<Vector3, 4>>& faces_params) {
  // face_pts_normal stores a vector of triangle
  // each triangle we store fid -> three points (pA, pB, pC) and normal nN
  faces_params.clear();
  for (const int& fid : k_ring_fids) {
    std::array<Vector3, 4> params;
    for (int lv = 0; lv < 3; lv++) {
      params[lv] = sf_mesh.vertices.point(sf_mesh.facets.vertex(fid, lv));
    }
    params[3] = get_mesh_facet_normal(sf_mesh, fid);
    faces_params[fid] = params;
  }
}

// Given sphere (theta, r), find the closest tangent point pN
// on triangle (pA, pB, pC) with normal nN
// by mininumizing energy function:
//
// This energy function has 2 parts with weights:
// 1. distance to tangent point -> alpha 1
// 2. distance to tangent plane -> alpha 2
double get_closest_tangent_point_opt(const double alpha1, const double alpha2,
                                     const Vector3& pA, const Vector3& pB,
                                     const Vector3& pC, const Vector3& nN,
                                     const Vector3& theta, const double& radius,
                                     Vector3& pN, bool is_debug) {
  Vector3 pAB = pB - pA;
  Vector3 pAC = pC - pA;
  //   Matrix3 N = nN * nN.transpose();  // symmetric matrix
  Matrix3 N = vec_vec_trans(nN, nN);
  Vector3 M1 = theta + radius * nN;
  Vector3 M2 = GEO::mult(N, theta) + radius * nN;

  // NOTE:
  // here N is a symmetric matrix
  // so "pAB.transpose() * N * pAB" (not compiling) can be written
  // as "(N * pAB).dot(pAB)", otherwise its wrong
  Matrix2 K;
  K(0, 0) =
      alpha1 * (GEO::dot(pAB, pAB)) + alpha2 * GEO::dot(GEO::mult(N, pAB), pAB);
  K(0, 1) =
      alpha1 * (GEO::dot(pAB, pAC)) + alpha2 * GEO::dot(GEO::mult(N, pAB), pAC);
  K(1, 0) =
      alpha1 * (GEO::dot(pAB, pAC)) + alpha2 * GEO::dot(GEO::mult(N, pAB), pAC);
  K(1, 1) =
      alpha1 * (GEO::dot(pAC, pAC)) + alpha2 * GEO::dot(GEO::mult(N, pAC), pAC);
  Vector2 b;
  b[0] = GEO::dot(alpha1 * M1 + alpha2 * M2, pAB) - alpha1 * GEO::dot(pA, pAB) -
         alpha2 * GEO::dot(GEO::mult(N, pA), pAB);
  b[1] = GEO::dot(alpha1 * M1 + alpha2 * M2, pAC) - alpha1 * GEO::dot(pA, pAC) -
         alpha2 * GEO::dot(GEO::mult(N, pA), pAC);

  Matrix2 K_inv;
  bool invertible = K.compute_inverse(K_inv);
  if (!invertible) {
    if (is_debug) {
      fprintf(stderr, "[CLOSEST_PNT] Matrix 2x2 not invertible \n");
      printf("%f, %f \n %f, %f \n", K(0, 0), K(0, 1), K(1, 0), K(1, 1));
    }
    return DBL_MAX;
  }
  // update mat_p
  Vector2 t1_t2 = GEO::mult(K_inv, b);
  double sum = t1_t2[0] + t1_t2[1];
  if (is_debug)
    printf("[CLOSEST_PNT] t1: %f, t2: %f, sum: %f", t1_t2[0], t1_t2[1], sum);
  if (sum <= 1 && sum >= 0 && t1_t2[0] <= 1 && t1_t2[0] >= 0 && t1_t2[1] <= 1 &&
      t1_t2[1] >= 0) {
    pN = pA + t1_t2[0] * pAB + t1_t2[1] * pAC;
    double E =
        TangentPlane::get_energy_value(theta, radius, alpha1, alpha2, pN, nN);
    if (is_debug)
      printf("[CLOSEST_PNT] return pN (%f,%f,%f) with energy E %f\n", pN[0],
             pN[1], pN[2], E);
    return E;
  }

  // closest point not in triangle, check 3 vertices
  double EA =
      TangentPlane::get_energy_value(theta, radius, alpha1, alpha2, pA, nN);
  double EB =
      TangentPlane::get_energy_value(theta, radius, alpha1, alpha2, pB, nN);
  double EC =
      TangentPlane::get_energy_value(theta, radius, alpha1, alpha2, pC, nN);
  if (EA <= EB && EA <= EC) {
    if (is_debug)
      printf("[CLOSEST_PNT] return pA (%f,%f,%f) with energy EA %f\n", pA[0],
             pA[1], pA[2], EA);
    pN = pA;
    return EA;
  }
  if (EB <= EA & EB <= EC) {
    if (is_debug)
      printf("[CLOSEST_PNT] return pB (%f,%f,%f) with energy EB %f\n", pB[0],
             pB[1], pB[2], EB);
    pN = pB;
    return EB;
  }

  if (is_debug)
    printf("[CLOSEST_PNT] return pC (%f,%f,%f) with energy EC %f\n", pC[0],
           pC[1], pC[2], EC);
  pN = pC;
  return EC;
}

void update_tangent_points_on_tan_pls(const SurfaceMesh& sf_mesh,
                                      const std::set<aint2>& fe_sf_fs_pairs,
                                      MedialSphere& mat_p, const double alpha1,
                                      const double alpha2, bool is_debug) {
  if (is_debug) {
    printf("calling update_tangent_points_on_tan_pls for mat_p %d \n",
           mat_p.id);
  }

  const int k = 2;  // k_ring sf_mesh face neighbors
  std::set<int> k_ring_fids;
  std::map<int, std::array<Vector3, 4>> faces_params;

  // store fids that are already stored
  // if another tangent plane found the same fids to update
  // then we delete duplicated tangent plane
  std::set<int> tangent_fids_now;
  std::vector<TangentPlane> tan_planes_new;

  for (auto& tan_pl : mat_p.tan_planes) {
    if (is_debug) {
      printf("[Update TAN_POINT] before update: \n");
      tan_pl.print_info();
    }
    get_k_ring_neighbors_no_cross(sf_mesh, fe_sf_fs_pairs, tan_pl.fid, k,
                                  k_ring_fids, true /*is_clear_cur*/,
                                  false /*is_debug*/);
    get_faces_params(sf_mesh, k_ring_fids, faces_params);
    // find tangent point that minimize the energy
    for (const auto& f_param_pair : faces_params) {
      int fid_tmp = f_param_pair.first;
      auto& f_param = f_param_pair.second;
      Vector3 pN_tmp;
      // Note that here we use alpha1=1, alpha2=0 to update tangent points
      double e_tmp = get_closest_tangent_point_opt(
          alpha1, alpha2, f_param[0], f_param[1], f_param[2], f_param[3],
          mat_p.center, mat_p.radius, pN_tmp, false /*is_debug*/);
      if (e_tmp < tan_pl.energy) {
        // update tangent point
        tan_pl.energy = e_tmp;
        tan_pl.fid = fid_tmp;
        tan_pl.tan_point = pN_tmp;
        tan_pl.normal = f_param[3];
        if (is_debug)
          printf(
              "[Update TAN_POINT] mat_p %d found lower energy %f with fid %d "
              "\n",
              mat_p.id, e_tmp, fid_tmp);
      }
    }  // for faces_params

    if (is_debug) {
      printf("[Update TAN_POINT] after update: \n");
    }

    ////////////////
    // after update
    // 1. remove duplicated tangent planes
    if (tangent_fids_now.find(tan_pl.fid) != tangent_fids_now.end()) {
      if (is_debug) {
        printf(
            "[Update TAN_POINT] mat_p {} with fid {} is stored, delete tan_pl",
            mat_p.id, tan_pl.fid);
      }
      tan_pl.is_deleted = true;
      continue;
    }
    // 2. remove similar tangent planes
    for (const auto& tan_pl_new : tan_planes_new) {
      if (tan_pl_new.is_same_tan_pl(tan_pl)) {
        if (is_debug) {
          printf(
              "[Update TAN_POINT] mat_p {} whose tangent plane already stored "
              "similar planes, skip storing",
              mat_p.id);
        }
        tan_pl.is_deleted = true;
        break;
      }
    }
    if (!tan_pl.is_deleted) {
      // update covered sf_mesh fids
      tan_pl.update_covered_sf_fids(sf_mesh);
      tan_planes_new.push_back(tan_pl);
      tangent_fids_now.insert(tan_pl.fid);
    }

    if (is_debug) tan_pl.print_info();
  }  // for mat_p.tan_planes

  mat_p.tan_planes = tan_planes_new;
}

// Given sphere (theta, r), find the closest tangent point pN
// and normal nN on concave segement (m1, m2) with adjacent normals (n1, n2)
// by solving an equation
double get_closest_concave_point(const double alpha3, const Vector3& m1,
                                 const Vector3& m2, const Vector3& n1,
                                 const Vector3& n2, const Vector3& theta,
                                 const double& radius, Vector3& pN, Vector3& nN,
                                 bool is_debug) {
  Vector3 m12 = m2 - m1;
  Vector3 vM = GEO::normalize(m12);
  double t = (GEO::dot(theta, vM) - GEO::dot(m1, vM)) / GEO::dot(m12, vM);
  double E = DBL_MAX;
  if (t >= 0 && t <= 1) {
    pN = m1 + t * m12;
    nN = GEO::normalize(pN - theta);
    // nN must in range of (n1, n2)
    if (is_vector_within_range_of_two_vectors(n1, n2, nN)) {
      E = TangentConcaveLine::get_energy_value(theta, radius, alpha3, pN, nN);
      if (is_debug)
        printf(
            "[CONCAVE_PNT] return pN (%f,%f,%f) with energy E %f, t: %f, nN: "
            "(%f,%f,%f)\n",
            pN[0], pN[1], pN[2], E, t, nN[0], nN[1], nN[2]);
      return E;
    }
  }
  // check boundary points/normals on segment
  // find the best combination with smallest energy
  // (pN, nN) pairs
  bool is_replaced = false;
  E = DBL_MAX;
  int pair_idx = -1;
  std::vector<std::array<Vector3, 2>> pN_nN_pairs{
      /*{{m1, nN}}, {{m2, nN}},*/
      {{pN, n1}}, {{pN, n2}}, {{m1, n1}}, {{m1, n2}}, {{m2, n1}}, {{m2, n2}}};
  for (int i = 0; i < pN_nN_pairs.size(); i++) {
    const auto& pair = pN_nN_pairs[i];
    double E_tmp = TangentConcaveLine::get_energy_value(theta, radius, alpha3,
                                                        pair[0], pair[1]);
    if (E_tmp < E) {
      pair_idx = i;
      E = E_tmp;
      pN = pair[0];
      nN = pair[1];
      // nN = (pN - theta).normalized();
      is_replaced = true;
    }
  }

  // till now, is_vector_within_range_of_two_vectors(n1, n2, nN) must be false
  if (!is_replaced) {
    // keep the pN, force nN to be n1
    // otherwise nN may not be in range of [n1, n2]
    printf("[CONCAVE_PNT] force normal to be n1 (%f,%f,%f)\n", n1[0], n1[1],
           n1[2]);
    nN = n1;
  }

  if (is_debug)
    printf(
        "[CONCAVE_PNT] return pair_idx: %d, pN (%f,%f,%f) nN: (%f,%f,%f) with "
        "energy E %f, t: %f, is_replaced: %d\n",
        pair_idx, pN[0], pN[1], pN[2], nN[0], nN[1], nN[2], E, t, is_replaced);

  return E;
}

void update_tangent_points_on_cc_lines(MedialSphere& mat_p, const double alpha3,
                                       bool is_debug) {
  if (is_debug && mat_p.tan_cc_lines.size() > 0)
    printf("calling update_tangent_points_on_cc_lines for mat_p %d...\n",
           mat_p.id);

  for (auto& one_cc_line : mat_p.tan_cc_lines) {
    get_closest_concave_point(
        alpha3, one_cc_line.t2vs_pos[0], one_cc_line.t2vs_pos[1],
        one_cc_line.adj_normals[0], one_cc_line.adj_normals[1], mat_p.center,
        mat_p.radius, one_cc_line.tan_point, one_cc_line.normal, is_debug);
  }
}

/////////////////////////////////////////////////////////////////////////////////////
// Main Functions
/////////////////////////////////////////////////////////////////////////////////////
// Default parameters:
// alpha1 = 0.01;  // energy of distance to tangent point
// alpha2 = 1;     // energy of distance to tangent plane
// alpha3 = 1;     // energy of distance to concave line
//
// Two steps:
// 1. update medial sphere center and radius
// 2. update tangent planes and concave lines
bool iterate_sphere(const SurfaceMesh& sf_mesh, const AABBWrapper& aabb_wrapper,
                    const std::set<aint2>& fe_sf_fs_pairs,
                    const std::vector<FeatureEdge>& feature_edges,
                    MedialSphere& mat_p, bool is_debug,
                    iter_params iter_param) {
  if (is_debug) {
    printf("[Iterate] calling iterate_sphere for mat_p: %d, is_debug: %d \n",
           mat_p.id, is_debug);
  }

  int num_itr = 0;
  double sum_energy_over_sq_radius = -1;
  std::vector<Vector3> pN, nN;      // for tangent planes
  std::vector<Vector3> c_pN, c_nN;  // for tangent concave lines

  if (is_debug) {
    printf("[Iterate] before iteration, mat_p: %d has tangent planes: %ld \n",
           mat_p.id, mat_p.tan_planes.size());
    mat_p.print_info();
  }

  bool is_good = false;
  while (true) {
    // update is_good
    // also update (pN, nN)
    bool is_break = is_break_iteration(
        iter_param.itr_limit, num_itr, iter_param.alpha1, iter_param.alpha2,
        iter_param.alpha3, iter_param.break_threshold,
        iter_param.is_break_use_energy_over_sq_radius, feature_edges, sf_mesh,
        aabb_wrapper, mat_p, sum_energy_over_sq_radius, pN, nN, c_pN, c_nN,
        is_good, iter_param.is_check_new_tan_plane, is_debug);
    if (is_break) break;
    if (is_debug)
      printf("[Iterate] Update sphere then update tangent points \n");
    // Step 1: update sphere given tangent points
    Vector3 new_center(0., 0., 0.);
    double new_radius = 0.;
    bool is_success = update_sphere_given_plN_opt(
        iter_param.alpha1, iter_param.alpha2, pN, nN, iter_param.alpha3, c_pN,
        c_nN, new_center, new_radius, is_debug);
    if (is_debug)
      printf("[Iterate] update_sphere_given_plN_opt is_success: %d\n",
             is_success);
    if (!is_success || new_radius < SCALAR_ZERO_3) {
      //////////////////////////
      // Break condition part 4
      is_good = false;
      break;
    }
    mat_p.save_old_center_radius();
    mat_p.center = new_center;
    mat_p.radius = new_radius;

    // Step 2: update tangent points given new sphere
    update_tangent_points_on_tan_pls(sf_mesh, fe_sf_fs_pairs, mat_p,
                                     iter_param.alpha1, iter_param.alpha2,
                                     is_debug /*is_debug*/);
    update_tangent_points_on_cc_lines(mat_p, iter_param.alpha3, is_debug);
    if (is_debug) {
      printf("[Iterate] mat_p %d has tangent elements after update: \n",
             mat_p.id);
      mat_p.print_tan_planes();
      mat_p.print_tan_cc_lines();
    }

    num_itr++;
    // mat_p.store_sq_energy(aabb_wrapper);
  }  // while true

  // exist?
  if (is_good && mat_p.radius >= INIT_RADIUS - SCALAR_ZERO_3) {
    is_good = false;
    is_debug = true;
  }

  // purge tan_planes by tan_cc_lines
  mat_p.purge_and_delete_tan_planes();
  // update MedialSphere::covered_sf_fids_in_group
  mat_p.update_sphere_covered_sf_fids(sf_mesh, false /*is_debug*/);

  // // Post checking spheres
  // // if any tangent plane is too close to concave lines, then delete
  // // this is because we have inserted enough spheres around concave lines
  // // by setting pin point p on concave lines and shrink
  // //
  // // Note: this will not impact tan_cc_lines
  // for (const auto& tan_pl : mat_p.tan_planes) {
  //   mat_p.min_sq_dist_2cc = aabb_wrapper.get_sq_dist_to_ce(tan_pl.tan_point);
  //   if (mat_p.min_sq_dist_2cc <= SQ_DIST_TO_CC) {
  //     printf("[Iterate] mat_p %d too close to concave lines %f, not good \n",
  //            mat_p.id, mat_p.min_sq_dist_2cc);
  //     is_good = false;
  //     break;
  //   }
  // }  // for mat_p.tan_planes

  // Only place for adding SphereType::T_N_c spheres
  if (mat_p.tan_cc_lines.size() > 0 &&
      mat_p.tan_cc_lines.size() + mat_p.tan_planes.size() > 2) {
    mat_p.type = SphereType::T_N_c;
  } else if (mat_p.tan_planes.size() > 2) {
    mat_p.type = SphereType::T_3_MORE;
  }

  if (is_debug)
    printf("[Iterate] mat_p: %d iterated %d time, is_good: %d \n", mat_p.id,
           num_itr, is_good);

  // if (is_good == false || num_itr == 0) {
  // ninwang: num_itr==0 seems fine during some relaxation?
  if (is_good == false) {
    return false;
  }
  return true;
}

// Default parameters:
// alpha1 = 0.01;  // energy of distance to tangent point
// alpha2 = 1;     // energy of distance to tangent plane
// alpha3 = 1;     // energy of distance to concave line
//
// Two steps (reversed):
// 1. update tangent planes and concave lines
// 2. update medial sphere center and radius
bool iterate_sphere_reversed(const SurfaceMesh& sf_mesh,
                             const AABBWrapper& aabb_wrapper,
                             const std::set<aint2>& fe_sf_fs_pairs,
                             const std::vector<FeatureEdge>& feature_edges,
                             MedialSphere& mat_p, bool is_debug,
                             iter_params params) {
  if (is_debug) {
    printf(
        "[Iterate Reversed] calling iterate_sphere for mat_p: %d, is_debug: "
        "%d\n",
        mat_p.id, is_debug);
    printf("[Iterate Reversed] Update tangent points then update sphere \n");
  }

  int num_itr = 0;
  double sum_energy_over_sq_radius = -1;
  std::vector<Vector3> pN, nN;      // for tangent planes
  std::vector<Vector3> c_pN, c_nN;  // for tangent concave lines

  if (is_debug) {
    printf(
        "[Iterate Reversed] before iteration, mat_p: %d has tangent planes: "
        "%ld \n",
        mat_p.id, mat_p.tan_planes.size());
    mat_p.print_info();
  }

  bool is_good = false;
  while (true) {
    // update is_good
    // also update (pN, nN)
    bool is_break = is_break_iteration(
        params.itr_limit, num_itr, params.alpha1, params.alpha2, params.alpha3,
        params.break_threshold, params.is_break_use_energy_over_sq_radius,
        feature_edges, sf_mesh, aabb_wrapper, mat_p, sum_energy_over_sq_radius,
        pN, nN, c_pN, c_nN, is_good, params.is_check_new_tan_plane, is_debug);
    if (is_break) break;

    // Step 2: update tangent points given new sphere
    update_tangent_points_on_tan_pls(sf_mesh, fe_sf_fs_pairs, mat_p,
                                     params.alpha1, params.alpha2,
                                     is_debug /*is_debug*/);
    update_tangent_points_on_cc_lines(mat_p, params.alpha3, is_debug);
    if (is_debug) {
      printf(
          "[Iterate Reversed] mat_p %d has tangent elements after update: \n",
          mat_p.id);
      mat_p.print_tan_planes();
      mat_p.print_tan_cc_lines();
    }

    // Step 1: update sphere given tangent points
    Vector3 new_center(0., 0., 0.);
    double new_radius = 0.;
    bool is_success = update_sphere_given_plN_opt(
        params.alpha1, params.alpha2, pN, nN, params.alpha3, c_pN, c_nN,
        new_center, new_radius, is_debug);
    if (is_debug)
      printf("[Iterate Reversed] update_sphere_given_plN_opt is_success: %d\n",
             is_success);
    if (!is_success || new_radius < SCALAR_ZERO_3) {
      //////////////////////////
      // Break condition part 4
      is_good = false;
      break;
    }
    mat_p.save_old_center_radius();
    mat_p.center = new_center;
    mat_p.radius = new_radius;

    num_itr++;
    // mat_p.store_sq_energy(aabb_wrapper);
  }  // while true

  // Post checking spheres
  // if any tangent plane is too close to concave lines, then delete
  // this is because we have inserted enough spheres around concave lines
  // by setting pin point p on concave lines and shrink
  //
  // Note: this will not impact tan_cc_lines
  for (const auto& tan_pl : mat_p.tan_planes) {
    mat_p.min_sq_dist_2cc = aabb_wrapper.get_sq_dist_to_ce(tan_pl.tan_point);
    if (mat_p.min_sq_dist_2cc <= SQ_DIST_TO_CC) {
      if (is_debug)
        printf(
            "[Iterate Reversed] mat_p %d too close to concave lines %f, not "
            "good \n",
            mat_p.id, mat_p.min_sq_dist_2cc);
      is_good = false;
      break;
    }
  }  // for mat_p.tan_planes

  // Only place for adding SphereType::T_N_c spheres
  if (mat_p.tan_cc_lines.size() > 0 &&
      mat_p.tan_cc_lines.size() + mat_p.tan_planes.size() >= 3) {
    mat_p.type = SphereType::T_N_c;
    // mat_p.dilate_sphere_radius();
  }

  // update MedialSphere::covered_sf_fids_in_group
  mat_p.update_sphere_covered_sf_fids(sf_mesh, false /*is_debug*/);

  if (is_debug)
    printf("[Iterate Reversed] mat_p: %d iterated %d time, is_good: %d \n",
           mat_p.id, num_itr, is_good);

  // if (is_good == false || num_itr == 0) {
  // ninwang: num_itr==0 seems fine during some relaxation?
  if (is_good == false) {
    return false;
  }
  return true;
}

// void init_shrink_and_update(const SurfaceMesh& sf_mesh, const TetMesh&
// tet_mesh,
//                             std::vector<MedialSphere>& all_medial_spheres,
//                             int num_init_spheres, bool is_debug) {
//   // shrink once
//   init_and_shrink(sf_mesh, tet_mesh, all_medial_spheres, num_init_spheres,
//                   2 /*itr_limit*/, false);

//   double alpha1 = 0.01;  // energy of distance to tangent point
//   double alpha2 = 1;     // energy of distance to tangent plane
//   double alpha3 = 1;     // energy of distance to tangent plane
//   double break_threshold = SCALAR_ZERO_3;
//   int iteration_limit = 30;
//   uint num_active = 0;
//   for (auto& mat_p : all_medial_spheres) {
//     // mat_p.update_tan_planes_from_ss_params();
//     assert(!mat_p.tan_planes.empty());
//     iterate_sphere(sf_mesh, sf_mesh.aabb_wrapper, sf_mesh.fe_sf_fs_pairs,
//     mat_p,
//                    false /*is_debug*/, alpha1, alpha2, alpha3,
//                    break_threshold, iteration_limit);
//     if (!mat_p.is_deleted) num_active++;
//   }
//   printf("[Init] init %u spheres\n", num_active);
// }

// type_to_handle:
// 0 - all spheres except for external features
// 1 - SphereType::T_2 only
// 2 - internal features only
void relax_and_iterate_spheres(const SurfaceMesh& sf_mesh,
                               const std::vector<FeatureEdge>& feature_edges,
                               std::vector<MedialSphere>& all_medial_spheres,
                               const bool site_is_transposed,
                               std::vector<float>& site_updated,
                               int type_to_handle, bool is_debug) {
  int n_site = site_updated.size() / 3;
  assert(n_site == all_medial_spheres.size());
  assert(type_to_handle == 0 || type_to_handle == 1 || type_to_handle == 2);
  uint num_active = 0;

  for (auto& msphere : all_medial_spheres) {
    if (msphere.is_deleted) continue;
    if (msphere.is_on_corner() || msphere.is_on_se()) continue;
    if (msphere.is_on_ce_pin() || msphere.is_on_ce_pre()) continue;
    if (type_to_handle == 1 && msphere.type != SphereType::T_2) continue;
    if (type_to_handle == 2 && !msphere.is_on_intf()) continue;
    // do not check new tangent planes for SphereType::T_N_c spheres
    bool is_check_new_tan_plane =
        msphere.type == SphereType::T_N_c ? false : true;
    // if (msphere.id == 88)
    //   is_debug = true;
    // else
    is_debug = false;

    int i = msphere.id;
    MedialSphere new_msphere = msphere;  // copy
    // Note: keep the old radius to iterate,
    // otherwise threshold energy/radius will break
    new_msphere.save_old_center_radius(false /*is_clear*/);

    // update new center to the centroid of the power cell
    if (site_is_transposed) {
      new_msphere.center[0] = site_updated[i];
      new_msphere.center[1] = site_updated[i + n_site];
      new_msphere.center[2] = site_updated[i + (n_site << 1)];
    } else {
      new_msphere.center[0] = site_updated[3 * i];
      new_msphere.center[1] = site_updated[3 * i + 1];
      new_msphere.center[2] = site_updated[3 * i + 2];
    }

    if (is_debug)
      printf(
          "[Relax] before relax, msphere %d change center (%f,%f,%f) to "
          "(%f,%f,%f)\n",
          new_msphere.id, new_msphere.old_center[0], new_msphere.old_center[1],
          new_msphere.old_center[2], new_msphere.center[0],
          new_msphere.center[1], new_msphere.center[2]);

    // then iterate spheres
    // update tangent planes (step2) then update sphere (step1)
    if (!iterate_sphere_reversed(sf_mesh, sf_mesh.aabb_wrapper,
                                 sf_mesh.fe_sf_fs_pairs, feature_edges,
                                 new_msphere, is_debug /*is_debug*/)) {
      // msphere.is_deleted = true;
      // printf("[Relax] after relax, msphere %d is deleted \n", msphere.id);

      // revert back to old sphere
      new_msphere = msphere;
      printf(
          "[Relax] after relax, msphere %d is deleted, reverted to old "
          "sphere\n",
          msphere.id);
    } else {
      // update the current msphere
      msphere = new_msphere;
    }
    if (msphere.is_deleted) continue;
    num_active++;

    // TODO: update type of sphere after relaxation?
    if (msphere.type == SphereType::T_2 && msphere.tan_planes.size() >= 3) {
      msphere.type = SphereType::T_3_MORE;
    }

    if (is_debug)
      printf("[Relax] after relax, msphere %d change center to (%f,%f,%f)\n",
             msphere.id, msphere.center[0], msphere.center[1],
             msphere.center[2]);

  }  // for all_medial_spheres

  printf("[Relax] relax %u spheres\n", num_active);
}

// Optimized Delaunay Triangulation
// type_to_handle:
// 0 - all spheres except for external features
// 1 - SphereType::T_2 only
// 2 - internal features only
void relax_and_iterate_spheres_ODT(
    const SurfaceMesh& sf_mesh, const std::vector<FeatureEdge>& feature_edges,
    const MedialMesh& mmesh, std::vector<MedialSphere>& all_medial_spheres,
    int type_to_handle, bool is_debug) {
  assert(type_to_handle == 0 || type_to_handle == 1 || type_to_handle == 2);
  uint num_active = 0;

  for (auto& msphere : all_medial_spheres) {
    if (msphere.is_deleted) continue;
    if (msphere.is_on_corner() || msphere.is_on_se()) continue;
    if (msphere.is_on_ce_pin() || msphere.is_on_ce_pre()) continue;
    if (type_to_handle == 1 && msphere.type != SphereType::T_2) continue;
    if (type_to_handle == 2 && !msphere.is_on_intf()) continue;
    // do not check new tangent planes for SphereType::T_N_c spheres
    bool is_check_new_tan_plane =
        msphere.type == SphereType::T_N_c ? false : true;
    // if (msphere.id == 88)
    //   is_debug = true;
    // else
    is_debug = false;

    int i = msphere.id;
    MedialSphere new_msphere = msphere;  // copy
    // Note: keep the old radius to iterate,
    // otherwise threshold energy/radius will break
    new_msphere.save_old_center_radius(false /*is_clear*/);
    new_msphere.center = Vector3(0.0, 0.0, 0.0);

    // update new center to the centroid of the power cell
    double area_sum = 0;
    for (const auto& mfid : new_msphere.faces_) {
      const MedialFace& mface = mmesh.faces.at(mfid);
      if (mface.is_on_boundary)
        new_msphere.center += mface.area * mface.centroid;
      else
        new_msphere.center += mface.area * mface.circumcenter;
      area_sum += mface.area;
    }
    new_msphere.center /= area_sum;

    if (is_debug)
      printf(
          "[Relax ODT] before relax, msphere %d change center (%f,%f,%f) to "
          "(%f,%f,%f)\n",
          new_msphere.id, new_msphere.old_center[0], new_msphere.old_center[1],
          new_msphere.old_center[2], new_msphere.center[0],
          new_msphere.center[1], new_msphere.center[2]);

    // then iterate spheres
    // update tangent planes (step2) then update sphere (step1)
    if (!iterate_sphere_reversed(sf_mesh, sf_mesh.aabb_wrapper,
                                 sf_mesh.fe_sf_fs_pairs, feature_edges,
                                 new_msphere, is_debug /*is_debug*/)) {
      // msphere.is_deleted = true;
      // printf("[Relax ODT] after relax, msphere %d is deleted \n",
      // msphere.id);

      // revert back to old sphere
      new_msphere = msphere;
      printf(
          "[Relax ODT] after relax, msphere %d is deleted, reverted to old "
          "sphere\n",
          msphere.id);
    } else {
      // update the current msphere
      msphere = new_msphere;
    }
    if (msphere.is_deleted) continue;
    num_active++;

    // TODO: update type of sphere after relaxation?
    if (msphere.type == SphereType::T_2 && msphere.tan_planes.size() >= 3) {
      msphere.type = SphereType::T_3_MORE;
    }

    if (is_debug)
      printf(
          "[Relax ODT] after relax, msphere %d change center to (%f, %f, "
          "%f)\n ",
          msphere.id, msphere.center[0], msphere.center[1], msphere.center[2]);

  }  // for all_medial_spheres

  printf("[Relax ODT] relax %u spheres\n", num_active);
}

// relax spheres to centroids of filtered neighbors (Laplacian smoothing)
void relax_and_iterate_spheres_Laplacian(
    const SurfaceMesh& sf_mesh, const std::vector<FeatureEdge>& feature_edges,
    const MedialMesh& mmesh, std::vector<MedialSphere>& all_medial_spheres,
    bool is_debug) {
  uint num_active = 0;
  std::map<int, std::set<int>> sphere_neighbors;
  for (const auto& me : mmesh.edges) {
    sphere_neighbors[me.vertices_[0]].insert(me.vertices_[1]);
    sphere_neighbors[me.vertices_[1]].insert(me.vertices_[0]);
  }

  for (auto& msphere : all_medial_spheres) {
    if (msphere.is_deleted) continue;
    if (msphere.is_on_corner() || msphere.is_on_se()) continue;
    if (msphere.is_on_ce_pin()) continue;
    if (sphere_neighbors.at(msphere.id).size() < 1) continue;
    if (msphere.id == 6)
      is_debug = true;
    else
      is_debug = false;

    // Note: keep the old radius to iterate,
    // otherwise threshold energy/radius will break
    msphere.save_old_center_radius(false /*is_clear*/);

    // update new center to the centroid of filtered neighbors
    // 1. T2 consider all neighboring spheres
    // 2. T3 only consider TN/corners
    // 3. T4/TN no move
    Vector3 new_center(0., 0., 0.);
    const auto& all_neighbors = sphere_neighbors.at(msphere.id);
    if (msphere.type == SphereType::T_3_MORE &&
        msphere.tan_planes.size() > 3) {  // T4/TN
      continue;                           // no move
    }
    std::set<int> neighbors_filtered;
    if (msphere.type == SphereType::T_3_MORE) {  // T3
      for (const int& nid : all_neighbors) {
        const auto& neigh_sphere = all_medial_spheres.at(nid);
        if (neigh_sphere.type == SphereType::T_3_MORE ||
            neigh_sphere.is_on_corner()) {
          new_center += neigh_sphere.center;
          neighbors_filtered.insert(nid);
        }
      }
    } else if (msphere.type == SphereType::T_2) {
      for (const int& nid : all_neighbors)
        new_center += all_medial_spheres.at(nid).center;
      neighbors_filtered = all_neighbors;
    } else {
      // dunno how
      if (is_debug)
        printf("[Relax] msphere %d of type %d, no relax\n", msphere.id,
               msphere.type);
      continue;
    }
    // update new center
    msphere.center = new_center / neighbors_filtered.size();

    // if (is_debug)
    //   printf(
    //       "[Relax] before relax, msphere %d change center (%f,%f,%f) to "
    //       "(%f,%f,%f), found %ld neighbors_filtered \n",
    //       msphere.id, msphere.old_center[0], msphere.old_center[1],
    //       msphere.old_center[2], msphere.center[0], msphere.center[1],
    //       msphere.center[2], neighbors_filtered.size());

    // then iterate spheres
    // update tangent planes (step2) then update sphere (step1)
    if (!iterate_sphere_reversed(sf_mesh, sf_mesh.aabb_wrapper,
                                 sf_mesh.fe_sf_fs_pairs, feature_edges, msphere,
                                 is_debug /*is_debug*/)) {
      msphere.is_deleted = true;
      printf("[Relax] after relax, msphere %d is deleted \n", msphere.id);
    }
    if (msphere.is_deleted) continue;
    num_active++;

    // TODO: update type of sphere after relaxation?
    if (msphere.type == SphereType::T_2 && msphere.tan_planes.size() >= 3) {
      msphere.type = SphereType::T_3_MORE;
    }

    if (is_debug)
      printf("[Relax] after relax, msphere %d change center to (%f,%f,%f)\n",
             msphere.id, msphere.center[0], msphere.center[1],
             msphere.center[2]);

  }  // for all_medial_spheres

  printf("[Relax] relax %u spheres\n", num_active);
}

// update center based on
// 1. centroid of powercell
// 2. centroid of neighbors (Laplacian smoothing)
void relax_and_iterate_spheres_both(
    const SurfaceMesh& sf_mesh, const std::vector<FeatureEdge>& feature_edges,
    const MedialMesh& mmesh, std::vector<MedialSphere>& all_medial_spheres,
    const bool site_is_transposed, std::vector<float>& site_updated,
    bool is_debug) {
  int n_site = site_updated.size() / 3;
  assert(n_site == all_medial_spheres.size());
  uint num_active = 0;

  std::map<int, std::set<int>> sphere_neighbors;
  for (const auto& me : mmesh.edges) {
    sphere_neighbors[me.vertices_[0]].insert(me.vertices_[1]);
    sphere_neighbors[me.vertices_[1]].insert(me.vertices_[0]);
  }

  for (auto& msphere : all_medial_spheres) {
    if (msphere.is_deleted) continue;
    if (msphere.is_on_corner() || msphere.is_on_se()) continue;
    // if (msphere.is_on_ce_pin() || msphere.is_on_ce()) continue;
    if (msphere.is_on_ce_pin()) continue;
    // if (msphere.id == 155)
    // is_debug = true;
    // else
    // is_debug = false;

    int i = msphere.id;
    // Note: keep the old radius to iterate,
    // otherwise threshold energy/radius will break
    msphere.save_old_center_radius(false /*is_clear*/);

    ///////////////////////
    // step 1
    // update new center to the centroid of the power cell
    Vector3 new_center1(0., 0., 0.);
    if (site_is_transposed) {
      new_center1[0] = site_updated[i];
      new_center1[1] = site_updated[i + n_site];
      new_center1[2] = site_updated[i + (n_site << 1)];
    } else {
      new_center1[0] = site_updated[3 * i];
      new_center1[1] = site_updated[3 * i + 1];
      new_center1[2] = site_updated[3 * i + 2];
    }

    ///////////////////////
    // step 2
    // update new center to the centroid of filtered neighbors
    // 1. T2 consider all neighboring spheres
    // 2. T3 only consider TN/corners
    // 3. T4/TN no move
    Vector3 new_center2(0., 0., 0.);
    const auto& all_neighbors = sphere_neighbors.at(msphere.id);
    if (msphere.type == SphereType::T_3_MORE &&
        msphere.tan_planes.size() > 3) {  // T4/TN
      continue;                           // no move
    }
    std::set<int> neighbors_filtered;
    if (msphere.type == SphereType::T_3_MORE) {  // T3
      for (const int& nid : all_neighbors) {
        const auto& neigh_sphere = all_medial_spheres.at(nid);
        if (neigh_sphere.type == SphereType::T_3_MORE ||
            neigh_sphere.is_on_corner()) {
          new_center2 += neigh_sphere.center;
          neighbors_filtered.insert(nid);
        }
      }
    } else if (msphere.type == SphereType::T_2) {
      for (const int& nid : all_neighbors)
        new_center2 += all_medial_spheres.at(nid).center;
      neighbors_filtered = all_neighbors;
    } else {
      // dunno how
      if (is_debug)
        printf("[Relax] msphere %d of type %d, no relax\n", msphere.id,
               msphere.type);
      continue;
    }
    new_center2 /= neighbors_filtered.size();

    ///////////////////////
    // step 3
    // update new center with weights
    // 1. (alpha1) power will push spheres from bigger power to lower power
    // (center to sharp edges)
    // 2. (alpha2) Laplacian will shrink spheres
    // double alpha1 = 1. / 3., alpha2 = 2. / 3.;
    // double alpha1 = 2. / 3., alpha2 = 1. / 3.;
    double alpha1 = 1., alpha2 = 0;
    if (msphere.type == SphereType::T_3_MORE)
      msphere.center = new_center2;
    else
      msphere.center = alpha1 * new_center1 + alpha2 * new_center2;

    if (is_debug)
      printf(
          "[Relax] before relax, msphere %d change center (%f,%f,%f) to "
          "(%f,%f,%f), new_center1: (%f,%f,%f), new_center2: (%f,%f,%f)\n",
          msphere.id, msphere.old_center[0], msphere.old_center[1],
          msphere.old_center[2], msphere.center[0], msphere.center[1],
          msphere.center[2], new_center1[0], new_center1[1], new_center1[2],
          new_center2[0], new_center2[1], new_center2[2]);

    // then iterate spheres
    // update tangent planes (step2) then update sphere (step1)
    if (!iterate_sphere_reversed(sf_mesh, sf_mesh.aabb_wrapper,
                                 sf_mesh.fe_sf_fs_pairs, feature_edges, msphere,
                                 is_debug /*is_debug*/)) {
      msphere.is_deleted = true;
      printf("[Relax] after relax, msphere %d is deleted \n", msphere.id);
    }
    if (msphere.is_deleted) continue;
    num_active++;

    // TODO: update type of sphere after relaxation?
    if (msphere.type == SphereType::T_2 && msphere.tan_planes.size() >= 3) {
      msphere.type = SphereType::T_3_MORE;
    }

    if (is_debug)
      printf("[Relax] after relax, msphere %d change center to (%f,%f,%f)\n",
             msphere.id, msphere.center[0], msphere.center[1],
             msphere.center[2]);

  }  // for all_medial_spheres

  printf("[Relax] relax %u spheres\n", num_active);
}