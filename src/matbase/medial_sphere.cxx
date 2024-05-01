#include "medial_sphere.h"

#include <assert.h>
#include <geogram/basic/process.h>

#include "common_cxx.h"

/////////////////////////////////////////////////////////////////////////////////////
// TangentPlane
/////////////////////////////////////////////////////////////////////////////////////
TangentPlane::TangentPlane(const SurfaceMesh& sf_mesh, const Vector3& _normal,
                           const Vector3& _point, const int _fid) {
  normal = _normal;
  tan_point = _point;
  fid = _fid;
  energy = DBL_MAX;
  energy_over_sq_radius = DBL_MAX;
  is_deleted = false;
  this->update_covered_sf_fids(sf_mesh);
}

void TangentPlane::print_info() const {
  printf(
      "TangentPlane: is_deleted %d, fid: %d energy: %.10e, "
      "energy_over_sq_radius: %.10e, tan_point (%f,%f,%f), normal: (%f,%f,%f), "
      "sf_fids_covered size: %zu\n",
      is_deleted, fid, energy, energy_over_sq_radius, tan_point[0],
      tan_point[1], tan_point[2], normal[0], normal[1], normal[2],
      sf_fids_covered.size());
  printf("sf_fids_covered: [");
  for (const int fid : sf_fids_covered) {
    printf("%d, ", fid);
  }
  printf("]\n");
}

bool TangentPlane::is_same_normal(const Vector3& bnormal,
                                  const double eps_degree) const {
  return is_vector_same_direction(normal, bnormal, eps_degree);
}
bool TangentPlane::is_same_normal(const Vector3& anormal,
                                  const Vector3& bnormal,
                                  const double eps_degree) {
  return is_vector_same_direction(anormal, bnormal, eps_degree);
}

// 1. normals are simialr
// 2. tangent point is close (if !bool is_only_normal)
bool TangentPlane::is_same_tan_pl(const TangentPlane& tan_pl2,
                                  bool is_only_normal,
                                  double eps_p_dist) const {
  if (!is_same_normal(tan_pl2.normal)) {
    // if (is_debug) printf("NOT same normal!!\n");
    return false;
  }
  // if (is_debug) printf("same normal!!\n");
  if (is_only_normal) return true;
  double dist = get_distance_between_two_vectors(tan_point, tan_pl2.tan_point);
  if (dist > eps_p_dist) return false;
  return true;
}

// static function
double TangentPlane::get_energy_value(const Vector3& theta,
                                      const double& radius, const double alpha1,
                                      const double alpha2,
                                      const Vector3& tan_point,
                                      const Vector3& normal) {
  Vector3 T = theta + radius * normal - tan_point;
  double K = GEO::dot(tan_point - theta, normal) - radius;
  return alpha1 * GEO::dot(T, T) + alpha2 * std::pow(K, 2);
}
// sphere (theta, radiu)
double TangentPlane::update_energy_value(const Vector3& theta,
                                         const double& radius,
                                         const double alpha1,
                                         const double alpha2) {
  energy = get_energy_value(theta, radius, alpha1, alpha2, tan_point, normal);
  energy_over_sq_radius = energy / std::pow(radius, 2);
  return energy;
}

void TangentPlane::update_tan_point(const Vector3& _new_tan_point) {
  tan_point = (tan_point + _new_tan_point) / 2.;
}

// mainly for updating fids
// but also update points and normal
// make sure they are projected on the surface
void TangentPlane::update_by_sf_mesh(const GEO::Mesh& sf_mesh,
                                     const AABBWrapper& aabb_wrapper) {
  Vector3 near_p;
  int p_fid = aabb_wrapper.project_to_sf_get_nearest_face(tan_point);
  const Vector3 pnormal = get_mesh_facet_normal(sf_mesh, p_fid);
  fid = p_fid;
  normal = pnormal;
}

bool TangentPlane::update_covered_sf_fids(const SurfaceMesh& sf_mesh, int k) {
  return sf_mesh.collect_kring_neighbors_given_fid(k, this->fid,
                                                   this->sf_fids_covered);
}

/////////////////////////////////////////////////////////////////////////////////////
// TangentConcaveLine
/////////////////////////////////////////////////////////////////////////////////////
TangentConcaveLine::TangentConcaveLine(const SurfaceMesh& sf_mesh, int _id,
                                       const FeatureEdge& fe) {
  id = _id;
  id_fe = fe.id;
  id_fl = fe.t2vs_group[2];
  t2vs_ids = {{fe.t2vs_group[0], fe.t2vs_group[1]}};
  t2vs_pos = fe.t2vs_pos;
  adj_sf_fs_pair = fe.adj_sf_fs_pair;
  adj_normals = fe.adj_normals;
  direction = GEO::normalize(t2vs_pos[1] - t2vs_pos[0]);
  tan_point = (t2vs_pos[0] + t2vs_pos[1]) / 2.;
  double t = RANDOM_01();
  normal = GEO::normalize(t * adj_normals[0] + (1. - t) * adj_normals[1]);
  energy = DBL_MAX;
  energy_over_sq_radius = DBL_MAX;
  is_deleted = false;
  this->update_covered_sf_fids(sf_mesh);
}

void TangentConcaveLine::print_info() const {
  printf(
      "---- TangentConcaveLine: id: %d, id_fe: %d, id_fl: %d, is_deleted, "
      "energy: %f, energy_over_sq_radius: %f, tan_point (%f,%f,%f), normal: "
      "(%f,%f,%f) ",
      id, id_fe, id_fl, is_deleted, energy, energy_over_sq_radius, tan_point[0],
      tan_point[1], tan_point[2], normal[0], normal[1], normal[2]);
  // printf("adj_normals (%f,%f,%f) - (%f,%f,%f), ", adj_normals[0][0],
  //        adj_normals[0][1], adj_normals[0][2], adj_normals[1][0],
  //        adj_normals[1][1], adj_normals[1][2]);
  printf("adj_sf_fs_pair: (%d,%d), ", adj_sf_fs_pair[0], adj_sf_fs_pair[1]);
  printf("sf_fids_covered_two size: (");
  for (const auto& fids_one : sf_fids_covered_two) {
    // printf("[");
    // for (const int fid : fids_one) {
    //   printf("%d, ", fid);
    // }
    // printf("], ");
    printf("%zu, ", fids_one.size());
  }
  printf(")\n");
}

bool TangentConcaveLine::operator==(TangentConcaveLine const& b) const {
  // belongs to the same grouped concave line
  if (id_fe == b.id_fe) return true;
  if (id_fl == b.id_fl) return true;
  return false;

  // Vector3 vec = (b.tan_point - tan_point).normalized();
  // if (adj_sf_fs_pair == b.adj_sf_fs_pair) return true;

  // // more likely to be different
  // // adj_sf_fs_pair might be different, but represent the same line
  // if (!is_vector_same_direction(direction, b.direction, esp_degree_5) &&
  //     !is_vector_oppo_direction(direction, b.direction, esp_degree_5))
  //   return false;

  // if (!is_vector_same_direction(vec, direction, esp_degree_5) &&
  //     !is_vector_oppo_direction(vec, direction, esp_degree_5))
  //   return false;

  // if (!is_vector_same_direction(vec, b.direction, esp_degree_5) &&
  //     !is_vector_oppo_direction(vec, b.direction, esp_degree_5))
  //   return false;

  // // vec should parallel to both direction and b.direction
  // return true;
}

// void TangentConcaveLine::print_info() const {
//   printf(
//       "TangentConcaveLine: energy: %f, energy_over_sq_radius: %f, "
//       "adj_sf_fs_pair: (%d,%d), tan_point (%f,%f,%f), normal: "
//       "(%f,%f,%f) \n",
//       energy, energy_over_sq_radius, adj_sf_fs_pair[0], adj_sf_fs_pair[1],
//       tan_point[0], tan_point[1], tan_point[2], normal[0], normal[1],
//       normal[2]);
//   printf("adj_normals (%f,%f,%f) - (%f,%f,%f)", adj_normals[0][0],
//          adj_normals[0][1], adj_normals[0][2], adj_normals[1][0],
//          adj_normals[1][1], adj_normals[1][2]);
// }

// if concave line is a curve
// then each line should cover more adjacent faces
bool TangentConcaveLine::is_normal_covered_by_adj_fs(
    const Vector3& n, double esp_degree_given) const {
  if (is_vector_same_direction(adj_normals[0], n, esp_degree_given))
    return true;
  if (is_vector_same_direction(adj_normals[1], n, esp_degree_given))
    return true;
  return false;
}

// // given two points of a segment
// void TangentConcaveLine::update_tangent_point(const GEO::vec3& v1_p,
//                                               const GEO::vec3& v2_p) {
//   GEO::vec3 v_mid = 1. / 2. * (v1_p + v2_p);
//   Vector3 p2(v_mid[0], v_mid[1], v_mid[2]);
//   if (!is_tan_point_updated)
//     tan_point = p2;
//   else
//     tan_point = 1. / 2. * (p2 + tan_point);
// }

bool TangentConcaveLine::purge_one_tan_plane(
    TangentPlane& one_tan_plane) const {
  // tan_pl has common sf_mesh fids as tan_cc_line
  for (const auto& sf_fids_covered : sf_fids_covered_two) {
    if (has_intersection(sf_fids_covered, one_tan_plane.sf_fids_covered)) {
      one_tan_plane.is_deleted = true;
      return true;
    }
  }
  return false;
}

void TangentConcaveLine::purge_tan_planes(std::vector<TangentPlane>& tan_planes,
                                          bool is_debug) const {
  int cnt = 0;
  for (auto& tan_pl : tan_planes) {
    if (purge_one_tan_plane(tan_pl)) cnt++;
  }
  if (is_debug) printf("purged %d tangent planes\n", cnt);
}

// static function
double TangentConcaveLine::get_energy_value(const Vector3& theta,
                                            const double& radius,
                                            const double alpha3,
                                            const Vector3& tan_point,
                                            const Vector3& normal) {
  Vector3 T = theta + radius * normal - tan_point;
  return alpha3 * GEO::dot(T, T);
}
// sphere (theta, radiu)
double TangentConcaveLine::update_energy_value(const Vector3& theta,
                                               const double& radius,
                                               const double alpha3) {
  // sanity_check();
  energy = get_energy_value(theta, radius, alpha3, tan_point, normal);
  energy_over_sq_radius = energy / std::pow(radius, 2);
  return energy;
}

bool TangentConcaveLine::update_covered_sf_fids(const SurfaceMesh& sf_mesh,
                                                int k) {
  this->sf_fids_covered_two.clear();
  this->sf_fids_covered_two.resize(2);
  FOR(i, 2) {  // two adjacent fids
    int adj_fid = this->adj_sf_fs_pair[i];
    sf_mesh.collect_kring_neighbors_given_fid(k, adj_fid,
                                              this->sf_fids_covered_two[i]);
  }
}

/////////////////////////////////////////////////////////////////////////////////////
// MedialSphere
/////////////////////////////////////////////////////////////////////////////////////
MedialSphere::MedialSphere(int _id, Vector3 _pin, Vector3 _pin_normal,
                           SphereType _type, int _itr) {
  id = _id;
  is_rt_valid = false;
  is_rt_prev_valid = false;
  rt_change_status = RtValidStates::NoChange;
  site_id = -1;
  ss.p = _pin;
  ss.p_normal = _pin_normal;
  radius = INIT_RADIUS;
  type = _type ? _type : SphereType::T_UNK;
  center = ss.p - ss.p_normal * radius;
  pcell.topo_status = Topo_Status::unkown;
  itr_cnt = _itr;
}
MedialSphere::MedialSphere(int _id, Vector3 _center, double _radius,
                           SphereType _type, int _itr) {
  id = _id;
  is_rt_valid = false;
  is_rt_prev_valid = false;
  rt_change_status = RtValidStates::NoChange;
  site_id = -1;
  center = _center;
  radius = _radius;
  type = _type;
  pcell.topo_status = Topo_Status::unkown;
  itr_cnt = _itr;
}

void MedialSphere::clear_mm() {
  this->edges_.clear();
  this->faces_.clear();
}

void MedialSphere::print_ss_info() const {
  printf("ss_params: p: (%f,%f,%f), q: (%f,%f,%f), q_fid: %d, p_fid: %d \n",
         ss.p[0], ss.p[1], ss.p[2], ss.q[0], ss.q[1], ss.q[2], ss.q_fid,
         ss.p_fid);
}

void MedialSphere::print_info() const {
  printf("------ MedialSphere Info ------\n");
  printf(
      "id: %d, is_deleted: %d, dup_cnt: %d, is_on_se: %d, is_on_corner: %d, "
      "is_on_ce_pre: %d, is_on_ce_pre_or_fix: %d, se_line_id: %d\n",
      id, is_deleted, dup_cnt, is_on_se(), is_on_corner(), is_on_ce_pre(),
      is_on_ce_pre_or_fix(), se_line_id);
  printf(
      "center: (%f,%f,%f), radius: %f, type: %d, is_radius_dilated: %d, "
      "min_sq_dist_2cc: %f\n",
      center[0], center[1], center[2], radius, type, is_radius_dilated,
      min_sq_dist_2cc);
  print_ss_info();
  printf("pcell: topo_status: %d, covered_sf_fids_in_group: %zu\n",
         pcell.topo_status, covered_sf_fids_in_group.size());
  printf("RT is_rt_valid: %d, rt_change_status: %d\n", is_rt_valid,
         rt_change_status);
  printf("rt_neigh_ids_prev: [");
  for (const auto rt_neigh_id_prev : rt_neigh_ids_prev)
    printf("%d, ", rt_neigh_id_prev);
  printf("]\n");
  printf("rt_neigh_ids: [");
  for (const auto rt_neigh_id : rt_neigh_ids) printf("%d, ", rt_neigh_id);
  printf("]\n");
  print_tan_planes();
  print_tan_cc_lines();
  print_covered_sf_fids_in_group();
}

void MedialSphere::print_covered_sf_fids_in_group() const {
  printf("msphere %d has covered_sf_fids_in_group: \n", this->id);
  for (const auto& one_fid_group : covered_sf_fids_in_group) {
    printf("(");
    for (const auto& fid : one_fid_group) {
      printf("%d ", fid.second);
    }
    printf(")\n");
  }
}

double MedialSphere::update_all_energy_values(const double alpha1,
                                              const double alpha2,
                                              const double alpha3,
                                              bool is_debug) {
  double sum = 0.;
  for (auto& tan_pl : tan_planes) {
    if (tan_pl.is_deleted) continue;
    tan_pl.update_energy_value(center, radius, alpha1, alpha2);
    sum += tan_pl.energy;
  }

  for (auto& tan_cc_line : tan_cc_lines) {
    if (tan_cc_line.is_deleted) continue;
    tan_cc_line.update_energy_value(center, radius, alpha3);
    sum += tan_cc_line.energy;
    // if (is_debug) {
    //   printf(
    //       "In update_all_energy_values, tan_cc_line %d, fe %d, fl %d, has "
    //       "energy: %f \n",
    //       tan_cc_line.id, tan_cc_line.id_fe, tan_cc_line.id_fl,
    //       tan_cc_line.energy);
    // }
  }
  return sum;
}

void MedialSphere::print_tan_planes() const {
  printf("------- MedialSphere %d Tangent Planes: %ld\n", id,
         tan_planes.size());
  for (const auto& tan_pl : tan_planes) {
    tan_pl.print_info();
  }
}

void MedialSphere::print_tan_cc_lines() const {
  printf("------- MedialSphere %d Tangent ConcaveLines: %ld\n", id,
         tan_cc_lines.size());
  for (const auto& tan_cc_line : tan_cc_lines) {
    tan_cc_line.print_info();
  }
}

int MedialSphere::get_tan_element_size() const {
  return tan_planes.size() + tan_cc_lines.size();
}

int MedialSphere::get_sf_covered_group_size() const {
  // each tan_cc_line covers 2 sf_mesh group
  // we count these 2 as 1
  return covered_sf_fids_in_group.size() - tan_cc_lines.size();
}

void MedialSphere::save_old_center_radius(bool is_clear) {
  old_center = center;
  old_radius = radius;
  if (is_clear) {
    center = Vector3(0.0, 0.0, 0.0);
    radius = 0.;
  }
}

void MedialSphere::update_tan_planes_from_ss_params(const SurfaceMesh& sf_mesh,
                                                    bool is_update_p,
                                                    bool is_update_q) {
  assert(is_update_p || is_update_q);
  if (is_update_p) {
    TangentPlane tan_pl1(sf_mesh, ss.p_normal, ss.p, ss.p_fid);
    tan_planes.push_back(tan_pl1);
  }
  if (is_update_q) {
    TangentPlane tan_pl2(sf_mesh, ss.q_normal, ss.q, ss.q_fid);
    tan_planes.push_back(tan_pl2);
  }
}

void MedialSphere::update_tan_cc_lines_from_ss_params(
    const SurfaceMesh& sf_mesh, const std::vector<FeatureEdge>& fe,
    bool is_update_p, bool is_update_q) {
  assert(is_update_p || is_update_q);
  if (is_update_p) {
    TangentConcaveLine new_cc_line(sf_mesh, tan_cc_lines.size(),
                                   fe.at(ss.get_p_fid()));
    new_cc_line.tan_point = ss.p;
    new_cc_line.normal = ss.p_normal;
    new_cc_line.energy = 0.0;
    new_cc_line.energy_over_sq_radius = 0.0;
    tan_cc_lines.push_back(new_cc_line);
  }

  if (is_update_q) {
    TangentConcaveLine new_cc_line(sf_mesh, tan_cc_lines.size(),
                                   fe.at(ss.get_q_fid()));
    new_cc_line.tan_point = ss.q;
    new_cc_line.normal = ss.q_normal;
    new_cc_line.energy = 0.0;
    new_cc_line.energy_over_sq_radius = 0.0;
    // new_cc_line_no_dup(new_cc_line);
    tan_cc_lines.push_back(new_cc_line);
  }
}

// also update TangentPlane::sf_fids_covered
void MedialSphere::new_tan_plane_no_dup(const SurfaceMesh& sf_mesh,
                                        const Vector3& _normal,
                                        const Vector3& _point, const int _fid,
                                        bool is_only_normal) {
  bool is_debug = false;
  // bool is_debug = id == 237 ? true : false;
  // if (is_debug) print_tan_planes();
  TangentPlane tan_pl(sf_mesh, _normal, _point, _fid);
  // if (is_debug) {
  //   printf("+++");
  //   tan_pl.print_info();
  // }
  bool is_add = true;
  for (auto& curr_tan_pl : tan_planes) {
    if (curr_tan_pl.is_same_tan_pl(tan_pl, is_only_normal)) {
      // if (is_debug) printf("duplicated!! not add \n");
      curr_tan_pl.update_tan_point(tan_pl.tan_point);
      is_add = false;
      break;
    }
  }
  // if (is_debug) printf("is_add: %d \n", is_add);
  if (!is_add) return;
  // NOTE: do this!!
  tan_pl.update_covered_sf_fids(sf_mesh);
  tan_planes.push_back(tan_pl);
}

// Note: this function will maximize new tangent concave line's energy
bool MedialSphere::new_cc_line_no_dup(const TangentConcaveLine& one_cc_line) {
  bool is_added = true;
  for (const auto& curr_cc_line : tan_cc_lines) {
    if (curr_cc_line.id_fl == one_cc_line.id_fl ||
        curr_cc_line.id_fe == one_cc_line.id_fe) {
      // printf(
      //     "Not add: curr_cc_line id_fl: %d, id_fe: %d, one_cc_line "
      //     "id_fl: %d, id_fe: %d \n",
      //     curr_cc_line.id_fl, curr_cc_line.id_fe,
      //     one_cc_line.id_fl, one_cc_line.id_fe);
      is_added = false;
      break;
    }
  }
  if (is_added) {
    tan_cc_lines.push_back(one_cc_line);
    tan_cc_lines.back().id = tan_cc_lines.size() - 1;
    // !!!! prepare for new iterations
    tan_cc_lines.back().energy = DBL_MAX;
    tan_cc_lines.back().energy_over_sq_radius = DBL_MAX;
  }
  return is_added;
}

// tangent pair: (point, normal)
void MedialSphere::get_sphere_all_tangent_pairs_includes_cc_lines(
    std::vector<std::array<Vector3, 2>>& tan_pairs) const {
  tan_pairs.clear();
  if (tan_cc_lines.size() > 0) {
    for (const auto& cc_line : tan_cc_lines) {
      tan_pairs.push_back({{cc_line.tan_point, cc_line.adj_normals[0]}});
      tan_pairs.push_back({{cc_line.tan_point, cc_line.adj_normals[1]}});
    }
  }
  for (const auto& tan_pl : tan_planes) {
    tan_pairs.push_back({{tan_pl.tan_point, tan_pl.normal}});
  }
}

// tangent pair: (point, normal, fid)
void MedialSphere::get_sphere_all_tangent_pairs_includes_cc_lines(
    std::vector<avec2int>& tan_pairs) const {
  tan_pairs.clear();
  if (tan_cc_lines.size() > 0) {
    for (const auto& cc_line : tan_cc_lines) {
      avec2int tmp1;
      tmp1.p = cc_line.tan_point;
      tmp1.n = cc_line.adj_normals[0];
      tmp1.f = cc_line.adj_sf_fs_pair[0];
      tan_pairs.push_back(tmp1);

      avec2int tmp2;
      tmp2.p = cc_line.tan_point;
      tmp2.n = cc_line.adj_normals[1];
      tmp2.f = cc_line.adj_sf_fs_pair[1];
      tan_pairs.push_back(tmp2);
    }
  }
  for (const auto& tan_pl : tan_planes) {
    avec2int tmp;
    tmp.p = tan_pl.tan_point;
    tmp.n = tan_pl.normal;
    tmp.f = tan_pl.fid;
    tan_pairs.push_back(tmp);
  }
}

// mainly for updating fids
void MedialSphere::update_tan_planes_by_sf_mesh(
    const GEO::Mesh& sf_mesh, const AABBWrapper& aabb_wrapper) {
  for (auto& tan_pl : tan_planes) {
    tan_pl.update_by_sf_mesh(sf_mesh, aabb_wrapper);
  }
}

// [no use]
void MedialSphere::update_tangent_covered_fids_by_sf_mesh(
    const SurfaceMesh& sf_mesh, int k) {
  for (auto& tan_pl : tan_planes) tan_pl.update_covered_sf_fids(sf_mesh, k);
  for (auto& tan_cc_line : tan_cc_lines)
    tan_cc_line.update_covered_sf_fids(sf_mesh, k);
}

void MedialSphere::purge_and_delete_tan_planes() {
  for (const auto& tan_cc_line : tan_cc_lines) {
    tan_cc_line.purge_tan_planes(tan_planes);
  }
  remove_deleted_tangents(false /*is_run_cc*/);
}

void MedialSphere::remove_deleted_tangents(bool is_run_cc) {
  std::vector<TangentPlane> new_tan_planes;
  for (const auto& tan_pl : tan_planes)
    if (!tan_pl.is_deleted) new_tan_planes.push_back(tan_pl);
  this->tan_planes = new_tan_planes;

  if (is_run_cc) {
    std::vector<TangentConcaveLine> new_cc_lines;
    for (const auto& tan_cc : tan_cc_lines)
      if (!tan_cc.is_deleted) new_cc_lines.push_back(tan_cc);
    this->tan_cc_lines = new_cc_lines;
  }
}

// Dilate the sphere radius to make it more robust
void MedialSphere::dilate_sphere_radius() {
  // make it more important
  this->radius *= 1.02;
  is_radius_dilated = true;
}

bool MedialSphere::is_same_tangent_info(const MedialSphere& msphereB) {
  if (this->tan_planes.size() != msphereB.tan_planes.size()) return false;
  if (this->tan_cc_lines.size() != msphereB.tan_cc_lines.size()) return false;

  // check tan_cc_line
  int cnt_same = 0;
  for (const auto& tan_cc_line_A : tan_cc_lines) {
    for (const auto& tan_cc_line_B : msphereB.tan_cc_lines) {
      if (tan_cc_line_A == tan_cc_line_B) {
        cnt_same++;
        break;
      }
    }
  }
  if (cnt_same != msphereB.tan_cc_lines.size()) return false;
  cnt_same = 0;
  for (const auto& tan_pl_A : tan_planes) {
    for (const auto& tan_pl_B : msphereB.tan_planes) {
      if (tan_pl_A.is_same_tan_pl(tan_pl_B, true /*is_only_normal*/)) {
        cnt_same++;
        break;
      }
    }
  }
  if (cnt_same != msphereB.tan_planes.size()) return false;
  return true;
}

// Update MedialSphere::covered_sf_fids_in_group
void MedialSphere::update_sphere_covered_sf_fids(const SurfaceMesh& sf_mesh,
                                                 bool is_debug) {
  auto merge_two_v2fid_vec = [&](std::vector<v2int>& v2fids1,
                                 const std::vector<v2int>& v2fids2) {
    // merge v2fids2 to v2fids1, without duplication
    std::set<int> visited_fids;
    for (const auto& v2fid1 : v2fids1) visited_fids.insert(v2fid1.second);
    for (const auto& v2fid2 : v2fids2) {
      if (visited_fids.find(v2fid2.second) != visited_fids.end()) continue;
      // save v2fid2
      v2fids1.push_back(v2fid2);
    }
  };
  auto check_and_merge_two_v2fid =
      [&](const std::set<int>& fids_to_check,
          std::vector<v2int>& one_group_fids,
          std::vector<std::vector<v2int>>& covered_sf_fids_in_group,
          std::map<int, int>& all_saved_fids) {
        assert(!fids_to_check.empty());
        std::vector<v2int> tmp_v2fid_group;
        sf_mesh.collect_fid_centroids(fids_to_check, tmp_v2fid_group);
        // if already saved, then get existing group
        int exist_group_id = -1;
        for (const auto& fid : fids_to_check) {
          if (all_saved_fids.find(fid) == all_saved_fids.end()) continue;
          exist_group_id = all_saved_fids.at(fid);
          break;
        }
        if (exist_group_id != -1) {
          // merge with current group
          one_group_fids = covered_sf_fids_in_group.at(exist_group_id);
          merge_two_v2fid_vec(one_group_fids, tmp_v2fid_group);
          // replace existing
          covered_sf_fids_in_group.at(exist_group_id) = one_group_fids;
        } else {
          // directly fetch covered sf_fids from tangent plane
          one_group_fids = tmp_v2fid_group;
          covered_sf_fids_in_group.push_back(one_group_fids);
          // assign new group_id
          exist_group_id = covered_sf_fids_in_group.size() - 1;
        }
        assert(!one_group_fids.empty());
        assert(exist_group_id != -1);

        // add as saved, only care about newly added
        for (const auto& fid : fids_to_check)
          all_saved_fids[fid] = exist_group_id;
      };

  // init
  this->covered_sf_fids_in_group.clear();
  std::vector<v2int> one_group_fids, _;
  // fid -> local_group_id (matching this->covered_sf_fids_in_group)
  std::map<int, int> all_saved_fids;

  // if (is_debug)
  //   printf("[Update_fids] fetching msphere %d's pcell covered_sf_fids\n",
  //          this->id);
  // // step1: collect from pcell's surf_v2fid_in_groups
  // //        store <pc_sf_centroid, sf_fid>
  // for (const auto& surf_v2fid : this->pcell.surf_v2fid_in_groups) {
  //   one_group_fids.clear();
  //   one_group_fids.insert(one_group_fids.end(), surf_v2fid.begin(),
  //                         surf_v2fid.end());
  //   this->covered_sf_fids_in_group.push_back(one_group_fids);
  //   for (const auto& one_v2int : one_group_fids)
  //     all_saved_fids[one_v2int.second] =
  //         this->covered_sf_fids_in_group.size() - 1;
  // }
  // // collect from pcell's 1ring
  // std::set<int> kring_neighbors, __;
  // for (auto& one_group_fids : this->covered_sf_fids_in_group) {
  //   kring_neighbors.clear();
  //   // get 1ring of one_group_fids
  //   for (const auto& v2fid : one_group_fids) {
  //     __.clear();
  //     sf_mesh.collect_kring_neighbors_given_fid(1 /*k*/, v2fid.second, __);
  //     kring_neighbors.insert(__.begin(), __.end());
  //   }
  //   sf_mesh.collect_fid_centroids(kring_neighbors, _);
  //   merge_two_v2fid_vec(one_group_fids, _);
  // }

  if (is_debug)
    printf(
        "[Update_fids] fetching msphere %d's tangent element's "
        "covered_sf_fids\n",
        this->id);

  // step2: merge with this->covered_sf_fids_in_group
  //        [collect k-ring of tangent fids (not crossing feature edges) after
  //        iterate/shrinking]
  for (const auto& tan_pl : this->tan_planes) {
    // should already updated
    assert(!tan_pl.sf_fids_covered.empty());
    check_and_merge_two_v2fid(tan_pl.sf_fids_covered, one_group_fids,
                              this->covered_sf_fids_in_group, all_saved_fids);
  }
  for (const auto& tan_cc_line : this->tan_cc_lines) {
    FOR(i, 2) {  // two adjacent fids
      int adj_fid = tan_cc_line.adj_sf_fs_pair[i];
      const auto& covered_fids = tan_cc_line.sf_fids_covered_two[i];
      // should already updated
      assert(!covered_fids.empty());
      check_and_merge_two_v2fid(covered_fids, one_group_fids,
                                this->covered_sf_fids_in_group, all_saved_fids);
    }  // for two adjacent fids
  }

  if (is_debug)
    printf("[Surf_fids] msphere %d collected covered_sf_fids_in_group: %zu\n",
           this->id, this->covered_sf_fids_in_group.size());
}

void MedialSphere::pcell_insert(int cell_id) { pcell.cell_ids.insert(cell_id); }

void MedialSphere::topo_clear() {
  // cells
  pcell.cell_ids.clear();
  pcell.tet_ids.clear();
  pcell.cell_neighbors.clear();
  pcell.cell_to_surfv2fid.clear();
  pcell.cc_cells.clear();
  // pcell.cc_bfids.clear();
  pcell.cc_surf_v2fids.clear();
  // facets
  // pcell.tfid_to_cells.clear();
  pcell.facet_neigh_to_cells.clear();
  pcell.cell_to_tfids.clear();
  pcell.facet_cc_cells.clear();
  pcell.facet_cc_surf_v2fids.clear();
  pcell.facet_neigh_is_fixed.clear();
  // powercell edges
  pcell.e_to_cells.clear();
  pcell.e_is_fixed.clear();
  pcell.edge_cc_cells.clear();
  pcell.edge_2endvertices.clear();
  // powercell vertices
  pcell.vertex_2id.clear();
  pcell.vertex_2pos.clear();
  // sharp edges
  pcell.se_covered_lvids.clear();
  pcell.se_line_endpos.clear();
  // concave edges
  pcell.ce_covered_lvids.clear();
  // surface faces
  pcell.surf_v2fid_in_groups.clear();

  pcell.topo_status = Topo_Status::ok;
  euler = 0;
  euler_sum = 0;
  num_cells = 0;
}

// We only consider Topo_Status::high_facet_cc as true
// iff both spheres have Topo_Status::high_facet_cc
bool MedialSphere::fcc_is_to_fix(int neigh_id) const {
  return pcell.facet_neigh_is_fixed.find(neigh_id) !=
         pcell.facet_neigh_is_fixed.end();
}

void MedialSphere::fcc_fixed(int neigh_id) {
  // only update when the topo_type is the same
  // if (pcell.topo_status == Topo_Status::high_facet_cc)
  // assert(pcell.facet_neigh_is_fixed.find(neigh_id) !=
  //        pcell.facet_neigh_is_fixed.end());
  pcell.facet_neigh_is_fixed[neigh_id] = true;
}

/**
 * Compares two 4D points with respect to the lexicographic order.
 * s1, s2 are the coordinates of the two 4D points.
 * true if s1 is strictly before s2 in the lexicographic order.
 * false otherwise.
 */
// SCALAR_ZERO seems to small, use SCALAR_ZERO_N3 instead
bool MedialSphere::operator<(const MedialSphere& m2) const {
  double threshold = SCALAR_ZERO_2;
  if (std::abs(center[0] - m2.center[0]) > threshold) {
    return center[0] < m2.center[0];
  }
  if (std::abs(center[1] - m2.center[1]) > threshold) {
    return center[1] < m2.center[1];
  }
  if (std::abs(center[2] - m2.center[2]) > threshold) {
    return center[2] < m2.center[2];
  }
  if (std::abs(radius - m2.radius) > threshold) {
    return radius < m2.radius;
  }

  /*
   * Read this if you see heap overflow!!
   * Compare requirement (https://en.cppreference.com/w/cpp/named_req/Compare)
   * For all a, comp(a,a)==false
   * If comp(a,b)==true then comp(b,a)==false
   * otherwise the std::sort would return heap-buffer-overflow
   */

  if (is_deleted && !m2.is_deleted) {
    return true;  // // cur < m2, not deleted sphere larger
  }

  // // check internal features sizes
  // // keep the one with more internal feature neighbors
  // if (intf_adj_intf.size() != m2.intf_adj_intf.size()) {
  //   return intf_adj_intf.size() < m2.intf_adj_intf.size();
  // }

  // keep feature spheres large
  if (m2.type == SphereType::T_1_N) {
    if (type != SphereType::T_1_N) {
      return true;  // cur < m2
    }
  }
  if (type == SphereType::T_1_N) {
    if (m2.type != SphereType::T_1_N) {
      return false;  // m2 < cur
    }
  }

  // keep sphere crossing concave lines large
  if (m2.type == SphereType::T_2_c) {
    if (type != SphereType::T_2_c) return true;  // cur < m2
  }
  if (m2.type == SphereType::T_X_c) {
    if (type != SphereType::T_X_c) return true;  // cur < m2
  }

  // also add type info, we want to keep the one with higher type
  if (type != SphereType::T_UNK && m2.type != SphereType::T_UNK) {
    return std::abs(type) < std::abs(m2.type);
  }
  if (type > 0 && type > 0) {
    return type < m2.type;
  }

  // when a == b, then always return false to avoid heap overflow
  return false;
}

// for removing medial spheres that are too close
// same as operator== but given threshold
bool MedialSphere::is_sphere_too_close(const MedialSphere& m2,
                                       double threshold) const {
  // our mesh has been normalized through
  // MeshIO::normalize_mesh()
  double diff = (center - m2.center).length();
  if (diff > threshold) return false;
  diff = std::abs(radius - m2.radius);
  if (diff > threshold) return false;
  return true;
}

// for removing medial spheres that are too close
// used by remove_duplicated_medial_spheres() ONLY
bool MedialSphere::operator==(const MedialSphere& m2) const {
  // our mesh has been normalized through
  // MeshIO::normalize_mesh()
  double threshold = SCALAR_ZERO_1;
  double diff = (center - m2.center).length();
  if (diff > threshold) return false;
  diff = std::abs(radius - m2.radius);
  if (diff > threshold) return false;
  return true;
}

bool MedialSphere::operator!=(const MedialSphere& m2) const {
  // if (id != m2.id) return true;
  double diff = (center - m2.center).length();
  if (diff > SCALAR_ZERO_1) return true;
  diff = std::abs(radius - m2.radius);
  if (diff > SCALAR_ZERO_1) return true;
  return false;
}

bool MedialSphere::is_on_se() const {
  if (type == SphereType::T_1_2) return true;
  return false;
}
bool MedialSphere::is_on_corner() const {
  if (type == SphereType::T_1_N) return true;
  return false;
}
bool MedialSphere::is_on_ce_pin() const {  // fake spheres
  if (type == SphereType::T_c) return true;
  return false;
}
bool MedialSphere::is_on_ce_pre() const {
  if (type == SphereType::T_2_c) return true;
  return false;
}
bool MedialSphere::is_on_ce_pre_or_fix() const {
  if (type == SphereType::T_2_c || type == SphereType::T_X_c) return true;
  return false;
}
bool MedialSphere::is_on_ce() const {
  if (type == SphereType::T_2_c || type == SphereType::T_X_c ||
      type == SphereType::T_N_c)
    return true;
  return false;
}
bool MedialSphere::is_on_intf() const {
  if (type == SphereType::T_N || type == SphereType::T_N_c) return true;
  return false;
}
bool MedialSphere::is_on_extf() const {
  if (is_on_se()) return true;
  if (is_on_corner()) return true;
  return false;
}

bool MedialSphere::is_on_sheet() const {
  if (is_on_se() || is_on_corner()) return false;
  if (is_on_intf()) return false;
  return true;
}

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
  if (msphere.type == SphereType::T_N) {
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

// setup different threshold
// we only use smaller threshold for topo_fix (is_small_threshold = true)
// others (extf_fix, intf_fix, geo_fix) use larger thrshold
// to avoid adding too many sphere at one place
bool validate_new_sphere(const std::vector<MedialSphere>& all_medial_spheres,
                         const MedialSphere& new_sphere,
                         bool is_small_threshold, bool is_debug) {
  double threshold = SCALAR_1;
  if (is_small_threshold) threshold = SCALAR_ZERO_1;

  // make sure new_sphere is not duplicated with any sphere
  for (const auto& msphere : all_medial_spheres) {
    if (msphere.is_sphere_too_close(new_sphere, threshold)) {
      if (is_debug)
        printf(
            "[NewSphereAdd Failed] new_sphere (%f,%f,%f,%f) too close to sphere"
            "%d (%f,%f,%f,%f) using threshold %f (is_small_threshold %d), not "
            "add\n",
            new_sphere.center[0], new_sphere.center[1], new_sphere.center[2],
            new_sphere.radius, msphere.id, msphere.center[0], msphere.center[1],
            msphere.center[2], msphere.radius, threshold, is_small_threshold);
      return false;
    }
  }
  return true;
}

bool add_new_sphere_validate(std::vector<MedialSphere>& all_medial_spheres,
                             MedialSphere& new_sphere, bool is_small_threshold,
                             bool is_debug) {
  if (new_sphere.is_deleted) return false;
  if (!validate_new_sphere(all_medial_spheres, new_sphere, is_small_threshold,
                           is_debug))
    return false;
  // save to add
  int old_sphere_id = new_sphere.id;
  new_sphere.id = all_medial_spheres.size();
  // update type:
  // internal feature condition
  // see check_and_fix_intf_by_adding_new_spheres()
  if (new_sphere.tan_cc_lines.size() > 0 &&
      new_sphere.tan_cc_lines.size() + new_sphere.tan_planes.size() > 2)
    new_sphere.type = SphereType::T_N_c;
  else if (new_sphere.tan_planes.size() > 2)
    new_sphere.type = SphereType::T_N;
  if (is_debug)
    printf("[NewSphereAdd Success] new_sphere added %d, old_id: %d, type %d\n",
           new_sphere.id, old_sphere_id, new_sphere.type);
  all_medial_spheres.push_back(new_sphere);
  return true;
}

// This function will update MedialSphere::id
void purge_deleted_medial_spheres(
    std::vector<MedialSphere>& all_medial_spheres) {
  std::vector<MedialSphere> purged;
  int new_id = 0, num_deleted = 0;
  for (const auto& msphere : all_medial_spheres) {
    if (msphere.is_deleted) {
      num_deleted++;
      continue;
    }
    purged.push_back(msphere);
    purged.back().id = new_id++;
  }

  if (num_deleted == 0) return;
  all_medial_spheres.clear();
  all_medial_spheres = purged;
  printf("[Purge] purged %d/%ld deleted spheres \n", num_deleted,
         all_medial_spheres.size());
}

// mark MedialSphere::is_delete as true
// then generate_RT_CGAL_and_purge_spheres() function will purge
// call ahead!
int remove_duplicated_medial_spheres(
    std::vector<MedialSphere>& all_medial_spheres) {
  // Remove duplicates
  // sort spheres and delete duplicated spheres
  printf("[Duplicate] start sorting given medial spheres ... \n");
  std::vector<MedialSphere> sorted_partial_medial_spheres = all_medial_spheres;
  std::sort(sorted_partial_medial_spheres.begin(),
            sorted_partial_medial_spheres.end());
  int curr_idx = 0;  // init
  std::set<int> deleted_set;
  for (int n_idx = 1; n_idx < sorted_partial_medial_spheres.size(); n_idx++) {
    MedialSphere& mat_curr = sorted_partial_medial_spheres[curr_idx];
    MedialSphere& mat_next = sorted_partial_medial_spheres[n_idx];
    // sanity check
    if (mat_curr != all_medial_spheres[mat_curr.id]) {
      mat_curr.print_info();
      all_medial_spheres[mat_curr.id].print_info();
      printf("[Duplicate] sorted_mat != all_medial_spheres \n");
      assert(false);
    }
    // if any already deleted
    if (mat_curr.is_deleted) {
      curr_idx = n_idx;
      continue;
    }
    if (mat_next.is_deleted) {
      continue;
    }
    // mat_next would have higher type than mat_curr,
    // so we delete mat_curr
    if (mat_curr == mat_next) {
      // if (mat_curr.intf_adj_intf.size() <= 2) {
      // this is not a internal junction, can be deleted
      // mark delete
      mat_curr.is_deleted = true;  // this is just a copy
      all_medial_spheres[mat_curr.id].is_deleted =
          true;  // this is the real deletion
      deleted_set.insert(mat_curr.id);
      // logger().debug("[Delete] delete mat_curr {}, mat_next {}",
      // mat_curr.id, mat_next.id);

      // // update internal feature adjacency
      // update_intf_adjacency(mat_curr.id, all_medial_spheres);
      // } else {
      //   logger().debug(
      //       "[Delete] mat_curr {} is an internal junction, cannot be deletd,
      //       " "intf_adj_intf: {}, mat_next {}", mat_curr.id,
      //       mat_curr.intf_adj_intf, mat_next.id);
      // }
    }
    curr_idx = n_idx;
  }

  printf("[Duplicate] Removed %ld duplicated MAT spheres \n",
         deleted_set.size());
  return deleted_set.size();
}

// check if two medial spheres on the same sharp edge
bool is_two_mspheres_on_same_se(const MedialSphere& msphere1,
                                const MedialSphere& msphere2) {
  if (!msphere1.is_on_se() || !msphere2.is_on_se()) return false;
  // std::set<int> s1_se_group, s2_se_group;
  // // format [cell_id, lvid1, lvid2, fe_line_id, fe_id]
  // for (const aint5& se_lvds : msphere1.pcell.se_covered_lvids)
  //   s1_se_group.insert(se_lvds[3]);  // index 3
  // for (const aint5& se_lvds : msphere2.pcell.se_covered_lvids)
  //   s2_se_group.insert(se_lvds[3]);  // index 3

  // printf("[FIX_INTF] medial edge (%d,%d) checking covered se_groups: \n",
  //        msphere1.id, msphere2.id);
  // printf("s1_se_group: ");
  // for (const int group : s1_se_group) printf("%d, ", group);
  // printf("\n");
  // printf("s2_se_group: ");
  // for (const int group : s2_se_group) printf("%d, ", group);
  // printf("\n");
  // if (s1_se_group != s2_se_group) {
  //   return false;
  // }

  if (msphere1.se_line_id != msphere2.se_line_id) return false;
  return true;
}

// must call before update_power_cells()
void update_se_tangent_planes(const SurfaceMesh& sf_mesh,
                              const std::vector<FeatureEdge>& feature_edges,
                              std::vector<MedialSphere>& all_medial_spheres,
                              bool is_clear_existing) {
  int cnt = 0;
  for (auto& msphere : all_medial_spheres) {
    if (!msphere.is_on_extf()) continue;

    // clear the old tangent info
    if (is_clear_existing) msphere.tan_planes.clear();

    // For SE spheres
    if (msphere.is_on_se()) {
      // update based on MedialSphere::se_edge_id
      if (msphere.se_edge_id < 0 ||
          msphere.se_edge_id >= feature_edges.size()) {
        continue;
      }
      const auto fe = feature_edges.at(msphere.se_edge_id);
      FOR(i, 2) {
        // also update TangentPlane::sf_fids_covered
        msphere.new_tan_plane_no_dup(sf_mesh, fe.adj_normals[i], msphere.center,
                                     fe.adj_sf_fs_pair[i],
                                     true /*is_only_normal*/);
      }
      continue;
    }
    // For corners
    if (msphere.is_on_corner()) {
      // update based on MedialSphere::corner_fes
      for (const auto& fe_id : msphere.corner_fes) {
        assert(fe_id > -1 && fe_id < feature_edges.size());
        const auto fe = feature_edges.at(fe_id);
        FOR(i, 2) {
          // also update TangentPlane::sf_fids_covered
          msphere.new_tan_plane_no_dup(sf_mesh, fe.adj_normals[i],
                                       msphere.center, fe.adj_sf_fs_pair[i],
                                       true /*is_only_normal*/);
        }
      }
      continue;
    }  // on corner
    cnt++;
  }

  printf("[UpdateTangentPlanes] updated tangent planes of %d spheres\n", cnt);
}

void copy_powercell_volume(const std::vector<float>& all_cell_vol,
                           std::vector<MedialSphere>& all_medial_spheres) {
  assert(all_cell_vol.size() == all_medial_spheres.size());
  for (int i = 0; i < all_medial_spheres.size(); i++) {
    auto& msphere = all_medial_spheres.at(i);
    msphere.pcell.pcell_vol = all_cell_vol.at(i);
  }
}

// load sphere that RT status has changed
// 1. newly added with num_itr_global
// 2. RT from invalid to valid: add to changed_spheres
// 3. RT from valid to invalid: add neighbors to changed_spheres
//
// changed_spheres does not contains sphere in invalid_spheres
// but contains neighbors of spheres in invalid_spheres
void load_changed_spheres(const int num_itr_global,
                          const std::vector<MedialSphere> all_medial_spheres,
                          std::set<int>& changed_spheres,
                          std::set<int>& invalid_spheres, bool is_debug) {
  changed_spheres.clear();
  invalid_spheres.clear();
  for (const auto& msphere : all_medial_spheres) {
    if (msphere.is_rt_valid && msphere.itr_cnt == num_itr_global) {
      changed_spheres.insert(msphere.id);
      // if contains old neighbors, then add old neighbors as well
      // this this for the case when msphere whose pos/radius has been updated
      // but msphere is not newly added
      for (const auto& prev_neigh_id : msphere.rt_neigh_ids_prev)
        changed_spheres.insert(prev_neigh_id);
    } else if (msphere.rt_change_status == RtValidStates::Invalid2Valid)
      changed_spheres.insert(msphere.id);
    else if (msphere.rt_change_status == RtValidStates::Valid2Invalid)
      invalid_spheres.insert(msphere.id);
    else {
      // no change do nothing
    }
  }

  if (is_debug)
    printf("[ChangeStatus] changed_spheres: %ld, invalid_spheres: %ld \n",
           changed_spheres.size(), invalid_spheres.size());

  // Add invalid_spheres' neighbors to changed_spheres
  //
  // Note: invalid spheres NOT EXIST in RT anymore,
  //       so we cannot use RT to trace their neighboring spheres,
  //       here we use MedialSpheres::pcell::facet_neigh_to_cells instead
  for (const auto& invalid_id : invalid_spheres) {
    const auto& inval_sphere = all_medial_spheres.at(invalid_id);
    if (inval_sphere.pcell.facet_neigh_to_cells.empty()) continue;
    // neigh_id -> { set of cell_ids in one facet CC }
    for (const auto& pair : inval_sphere.pcell.facet_neigh_to_cells) {
      // already mapped to MedialSphere::id during previous call of
      // function update_convex_cells_voro_and_tet_ids()
      const int& seed_neigh_id = pair.first;
      assert(seed_neigh_id >= 0 && seed_neigh_id < all_medial_spheres.size());
      // add invalid_spheres' neighbors to changed_spheres
      if (all_medial_spheres.at(seed_neigh_id).is_rt_valid)
        changed_spheres.insert(seed_neigh_id);
    }
  }

  // not true
  // // sanity check
  // for (int i : invalid_spheres) {
  //   assert(changed_spheres.find(i) == changed_spheres.end());
  // }

  if (is_debug)
    printf(
        "[ChangeStatus] changed_spheres: %ld after adding invalid_spheres's "
        "neighbors \n",
        changed_spheres.size());
}

std::vector<ConvexCellHost> merge_convex_cells(
    const std::vector<MedialSphere>& all_medial_spheres,
    const std::set<int>& spheres_and_1rings,
    const std::vector<ConvexCellHost>& convex_cells_prev,
    const std::vector<ConvexCellHost>& convex_cells_new, bool is_debug) {
  std::vector<ConvexCellHost> merged_convex_cells;
  for (const auto& cell_new : convex_cells_new) {
    // filter by spheres + 1ring, No 2ring
    if (spheres_and_1rings.find(cell_new.voro_id) == spheres_and_1rings.end())
      continue;

    merged_convex_cells.push_back(cell_new);
    merged_convex_cells.back().id = merged_convex_cells.size() - 1;

    // if (cell_new.voro_id == 377)
    //   printf("[Merge ConvexCells] cell_new %d has voro_id %d \n",
    //          merged_convex_cells.back().id,
    //          merged_convex_cells.back().voro_id);
  }

  for (const auto& cell_prev : convex_cells_prev) {
    // if prev voro_id became invalid, then skip

    // do not store cells of (spheres + 1ring)
    if (spheres_and_1rings.find(cell_prev.voro_id) != spheres_and_1rings.end())
      continue;
    // store prev
    merged_convex_cells.push_back(cell_prev);
    merged_convex_cells.back().id = merged_convex_cells.size() - 1;

    // if (cell_prev.voro_id == 377) {
    //   printf("[Merge ConvexCells] cell_prev %d has voro_id %d, tet_id: %d\n",
    //          merged_convex_cells.back().id,
    //          merged_convex_cells.back().voro_id,
    //          merged_convex_cells.back().tet_id);
    //   if (merged_convex_cells.back().tet_id == 900) {
    //     merged_convex_cells.back().print_info();
    //   }
    // }
  }

  printf("[Merge ConvexCells] #prev: %zu, #new: %zu, #merged: %zu\n",
         convex_cells_prev.size(), convex_cells_new.size(),
         merged_convex_cells.size());
  return merged_convex_cells;
}

// Note:
// 1. And common tangent plane will be stored in common_tan_pls. For tangent
// point of common tangent plane, we project the middle of two spheres
// centers A&B to common surface fids.
// [This projection is important !!!]
// To make sure internal feature sphere added will be around middle of A&B.
//
// 2. Any non-common tangent plane of A is diff from B store them in
// A_non_common_tan_pls, vice versa.
void A_B_spheres_common_diff_tangent_info_from_surface(
    const SurfaceMesh& sf_mesh, const MedialSphere& mat_A,
    const MedialSphere& mat_B, std::vector<TangentPlane>& common_tan_pls,
    std::vector<TangentPlane>& A_non_common_tan_pls,
    std::vector<TangentPlane>& B_non_common_tan_pls, bool is_debug) {
  auto get_fids_in_vec = [&](const std::vector<v2int>& one_v2int_group,
                             std::set<int>& fids_one_group) {
    fids_one_group.clear();
    for (const auto& one_v2int : one_v2int_group)
      fids_one_group.insert(one_v2int.second);
  };

  auto find_centroid_given_fid = [&](const std::vector<v2int>& one_v2int_group,
                                     int fid_given) {
    for (const auto& one_v2int : one_v2int_group) {
      if (one_v2int.second == fid_given) return one_v2int.first;
    }
  };

  auto filter_v2fid_given_fids = [&](const std::vector<v2int>& one_v2int_group,
                                     std::set<int>& fids_given) {
    std::vector<v2int> filtered_v2fid;
    for (const auto& one_v2int : one_v2int_group) {
      if (fids_given.find(one_v2int.second) != fids_given.end())
        filtered_v2fid.push_back(one_v2int);
    }
    return filtered_v2fid;
  };

  // project point p on fids, and create a new tangent plane
  auto project_and_store_tan_pl = [&](const std::set<int>& fids,
                                      const Vector3& p,
                                      bool is_update_covered_fids = false) {
    v2int v2fid_chosen;
    v2fid_chosen.second =
        project_point_onto_triangles(sf_mesh, fids, p, v2fid_chosen.first);
    TangentPlane common_tan_pl(
        sf_mesh, get_mesh_facet_normal(sf_mesh, v2fid_chosen.second),
        v2fid_chosen.first, v2fid_chosen.second);
    // only update if asked, otherwise use fids
    if (is_update_covered_fids)
      common_tan_pl.update_covered_sf_fids(sf_mesh);
    else
      common_tan_pl.sf_fids_covered = fids;
    return common_tan_pl;
  };

  if (mat_A.id == 678 && mat_B.id == 125 || mat_A.id == 125 && mat_B.id == 678)
    is_debug = true;
  else
    is_debug = false;

  // fetch all tangent info for A/B
  common_tan_pls.clear();
  A_non_common_tan_pls.clear();
  B_non_common_tan_pls.clear();
  Vector3 AB_mid_point = (mat_A.center + mat_B.center) / 2;
  std::set<int> A_one_fids_group, B_one_fids_group, B_common_indices;

  // // update MedialSphere::covered_sf_fids_in_group
  // if (mat_A.covered_sf_fids_in_group.empty())
  //   mat_A.update_sphere_covered_sf_fids(sf_mesh, false /*is_debug*/);
  // if (mat_B.covered_sf_fids_in_group.empty())
  //   mat_B.update_sphere_covered_sf_fids(sf_mesh, false /*is_debug*/);
  if (mat_A.covered_sf_fids_in_group.empty() ||
      mat_B.covered_sf_fids_in_group.empty()) {
    mat_A.print_covered_sf_fids_in_group();
    mat_B.print_covered_sf_fids_in_group();
  }

  // if (is_debug)
  if (mat_A.covered_sf_fids_in_group.empty() ||
      mat_B.covered_sf_fids_in_group.empty()) {
    printf(
        "[AB_common_diff] sphere %d has sf_fid_groups: %zu, sphere %d has "
        "sf_fid_groups: %zu\n",
        mat_A.id, mat_A.covered_sf_fids_in_group.size(), mat_B.id,
        mat_B.covered_sf_fids_in_group.size());
    mat_A.print_info();
    mat_B.print_info();
  }

  assert(!mat_A.covered_sf_fids_in_group.empty() &&
         !mat_B.covered_sf_fids_in_group.empty());

  // get A's diff from B
  for (int a = 0; a < mat_A.covered_sf_fids_in_group.size(); a++) {
    const auto& A_one_v2int_group = mat_A.covered_sf_fids_in_group.at(a);
    bool is_store_a = true;
    get_fids_in_vec(A_one_v2int_group, A_one_fids_group);

    // loop B's covered fids
    for (int b = 0; b < mat_B.covered_sf_fids_in_group.size(); b++) {
      const auto& B_one_v2int_group = mat_B.covered_sf_fids_in_group.at(b);
      get_fids_in_vec(B_one_v2int_group, B_one_fids_group);
      // get A&B common
      std::set<int> common_fid_one_group;
      set_intersection(A_one_fids_group, B_one_fids_group,
                       common_fid_one_group);
      if (!common_fid_one_group.empty()) {
        // save index for mat_B.covered_sf_fids_in_group
        B_common_indices.insert(b);
        ////////
        // save common tangent plane
        // project middle of two spheres A&B to common fids
        TangentPlane common_tan_pl =
            project_and_store_tan_pl(common_fid_one_group, AB_mid_point);
        common_tan_pls.push_back(common_tan_pl);
        is_store_a = false;
        if (is_debug) {
          printf(
              "[AB_common_diff] A %d has same sf_group as B %d with common "
              "size %d, a chosen fid %d, p (%f,%f,%f)\n",
              mat_A.id, mat_B.id, common_fid_one_group.size(),
              common_tan_pl.fid, common_tan_pl.tan_point[0],
              common_tan_pl.tan_point[1], common_tan_pl.tan_point[2]);
        }
        break;
      }
    }  // for mat_B.covered_sf_fids_in_group

    if (is_store_a) {
      get_fids_in_vec(A_one_v2int_group, A_one_fids_group);
      // get project of mat_A.center onto A_one_fids_group
      TangentPlane tan_pl_A =
          project_and_store_tan_pl(A_one_fids_group, mat_A.center);
      A_non_common_tan_pls.push_back(tan_pl_A);
    }
  }  // for matA

  // get B's diff from A
  for (int b = 0; b < mat_B.covered_sf_fids_in_group.size(); b++) {
    // remove common
    if (B_common_indices.find(b) != B_common_indices.end()) continue;
    // left non-common
    const auto& B_one_v2int_group = mat_B.covered_sf_fids_in_group.at(b);
    get_fids_in_vec(B_one_v2int_group, B_one_fids_group);
    // get project of mat_B.center onto B_one_fids_group
    TangentPlane tan_pl_B =
        project_and_store_tan_pl(B_one_fids_group, mat_B.center);
    B_non_common_tan_pls.push_back(tan_pl_B);
  }  // for mat_B

  if (is_debug)
    printf(
        "[B Diff] mat_A %d has diff %zu, mat_B %d has diff %d, common %zu "
        "tangent planes\n",
        mat_A.id, A_non_common_tan_pls.size(), mat_B.id,
        B_non_common_tan_pls.size(), common_tan_pls.size());
}

v2int get_v2fid_max_to_point(const GEO::Mesh& sf_mesh,
                             const std::vector<v2int> surf_v2fids,
                             const Vector3& point) {
  v2int v2fid_max_dist;
  double max_dist = -1;
  for (auto& v2fid : surf_v2fids) {
    double dist = GEO::Geom::distance(point, v2fid.first);
    if (dist > max_dist) {
      v2fid_max_dist = v2fid;
      max_dist = dist;
    }
  }
  return v2fid_max_dist;
}

v2int get_v2fid_min_to_point(const GEO::Mesh& sf_mesh,
                             const std::vector<v2int> surf_v2fids,
                             const Vector3& point) {
  v2int v2fid_min_dist;
  double min_dist = DBL_MAX;
  for (auto& v2fid : surf_v2fids) {
    double dist = GEO::Geom::distance(point, v2fid.first);
    if (dist < min_dist) {
      v2fid_min_dist = v2fid;
      min_dist = dist;
    }
  }
  return v2fid_min_dist;
}

int delete_degenerated_medial_spheres(
    std::vector<MedialSphere>& all_medial_spheres, bool is_debug) {
  // for tracking deleted spheres
  std::vector<bool> is_sphere_deleted_this_turn(all_medial_spheres.size(),
                                                false);
  auto run_thread_deletion = [&](int sphere_id) {
    MedialSphere& msphere = all_medial_spheres.at(sphere_id);
    // do not update is_sphere_deleted_this_turn
    if (msphere.is_deleted) return;
    if (msphere.pcell.cell_ids.empty()) {
      printf("[Degenerate] sphere %d has zero convex cells, return\n",
             msphere.id);
      return;
    }

    // -----------------------------------------------------------------------
    // Powercell FaceCC
    // make sure the number of FacetCC is the same for msphere.id and neigh_id
    // [msphere.id, neigh_id] -> { set of cell_ids in one facet CC }
    // -----------------------------------------------------------------------
    for (const auto& cur_facet_cc : msphere.pcell.facet_cc_cells) {
      const int& neigh_id = cur_facet_cc.first;
      if (all_medial_spheres.at(neigh_id).is_deleted) continue;
      const auto& neigh_facet_cc_cells =
          all_medial_spheres.at(neigh_id).pcell.facet_cc_cells;

      // neigh_id must also have this halfplane
      // if not, then delete msphere
      if (neigh_facet_cc_cells.find(msphere.id) == neigh_facet_cc_cells.end()) {
        // delete msphere.id
        printf(
            "[Degenerate] halfplane [%d,%d] not exist for sphere %d, delete "
            "sphere %d  \n",
            msphere.id, neigh_id, neigh_id, msphere.id);
        // msphere.fcc_fixed(neigh_id);  // no need to fix
        msphere.is_deleted = true;
        is_sphere_deleted_this_turn[msphere.id] = true;
        return;
      }

      // if the count of FacetCC for [msphere.id, neigh_id] is different,
      // then delete one of them
      if (cur_facet_cc.second.size() !=
          neigh_facet_cc_cells.at(msphere.id).size()) {
        int sphere_id_to_delete = msphere.id;
        if (cur_facet_cc.second.size() >
            neigh_facet_cc_cells.at(msphere.id).size()) {
          // delete msphere.id
          // msphere.fcc_fixed(neigh_id);  // no need to fix
          msphere.is_deleted = true;
        } else {
          // delete neigh_id
          sphere_id_to_delete = neigh_id;
          // all_medial_spheres.at(neigh_id).fcc_fixed(
          //     msphere.id);  // no need to fix
          all_medial_spheres.at(neigh_id).is_deleted = true;
        }
        printf(
            "[Degenerate] halfplane [%d,%d] has FacetCC: (%d,%d), degenerate "
            "happens delete %d \n",
            msphere.id, neigh_id, cur_facet_cc.second.size(),
            neigh_facet_cc_cells.at(msphere.id).size(), sphere_id_to_delete);
        is_sphere_deleted_this_turn[sphere_id_to_delete] = true;
        return;
      }
    }  // powercel faceCC

    // -----------------------------------------------------------------------
    // Powercell EdgeCC
    // make sure the number of EdgeCC is the same for sphere1, sphere2 and
    // sphere3 [sphere0, sphere1, sphere2] -> {set of cell_ids in one edge CC}
    // -----------------------------------------------------------------------
    for (const auto& cur_edge_cc : msphere.pcell.edge_cc_cells) {
      const aint2& neigh_min_max = cur_edge_cc.first;
      if (all_medial_spheres.at(neigh_min_max[0]).is_deleted ||
          all_medial_spheres.at(neigh_min_max[1]).is_deleted)
        continue;

      // two neigh sphere1 and sphere2
      int num_edge_cc12[2] = {-1, -1};
      for (int i = 0; i < 2; i++) {
        int neigh1_id = neigh_min_max[i];
        int neigh2_id = neigh_min_max[(i + 1) % 2];
        const auto& neigh_edge_cc_cells =
            all_medial_spheres.at(neigh1_id).pcell.edge_cc_cells;
        aint2 neigh_neigh_min_max = {{msphere.id, neigh2_id}};
        std::sort(neigh_neigh_min_max.begin(), neigh_neigh_min_max.end());

        // cannot find other two
        // delete neigh1_id
        if (neigh_edge_cc_cells.find(neigh_neigh_min_max) ==
            neigh_edge_cc_cells.end()) {
          printf(
              "[Degenerate] poweredge [%d,%d,%d], sphere %d cannot find other "
              "two [%d,%d], delete %d\n ",
              msphere.id, neigh_min_max[0], neigh_min_max[1], neigh1_id,
              neigh_neigh_min_max[0], neigh_neigh_min_max[1], neigh1_id);
          all_medial_spheres.at(neigh1_id).is_deleted = true;
          is_sphere_deleted_this_turn[neigh1_id] = true;
          return;
        }
        // else: update num_edge_cc12
        num_edge_cc12[i] = neigh_edge_cc_cells.at(neigh_neigh_min_max).size();
      }  // for two neigh spheres
      assert(num_edge_cc12[0] != -1 && num_edge_cc12[1] != -1);

      // if the count of EdgeCC for [sphere0, sphere1, sphere2] is different,
      int num_edge_cc0 = cur_edge_cc.second.size();
      if (num_edge_cc0 != num_edge_cc12[0] ||
          num_edge_cc0 != num_edge_cc12[1]) {
        // int my_num_edge_cc[3] = {num_edge_cc0, num_edge_cc12[0],
        // num_edge_cc12[1]}; int my_spheres[3] = {msphere.id,
        // neigh_min_max[0], neigh_min_max[1]}; int sphere_to_delete =
        //     my_spheres[*std::max_element(my_num_edge_cc, my_num_edge_cc +
        //     3)];
        // all_medial_spheres.at(sphere_to_delete).is_deleted = true;
        // num_sphere_change++;
        // printf(
        //     "[Degenerate] poweredge [%d,%d,%d] has FacetCC: (%d,%d,%d),
        //     degenerate, " "delete sphere %d\n", msphere.id, neigh_min_max[0],
        //     neigh_min_max[1], num_edge_cc0, num_edge_cc12[0],
        //     num_edge_cc12[1], sphere_to_delete);

        // delete them all
        msphere.is_deleted = true;
        all_medial_spheres.at(neigh_min_max[0]).is_deleted = true;
        all_medial_spheres.at(neigh_min_max[1]).is_deleted = true;
        is_sphere_deleted_this_turn[msphere.id] = true;
        is_sphere_deleted_this_turn[neigh_min_max[0]] = true;
        is_sphere_deleted_this_turn[neigh_min_max[1]] = true;
        printf(
            "[Degenerate] poweredge [%d,%d,%d] has FacetCC: (%d,%d,%d), "
            "degenerate, delete all three spheres\n",
            msphere.id, neigh_min_max[0], neigh_min_max[1], num_edge_cc0,
            num_edge_cc12[0], num_edge_cc12[1]);
        return;
      }
    }  // powercell edgeCC
  };

  // for (int i = 0; i < all_medial_spheres.size(); i++) {
  //   run_thread_deletion(i);
  // }
  GEO::parallel_for(0, all_medial_spheres.size(),
                    [&](int sphere_id) { run_thread_deletion(sphere_id); });

  // count the sphere deleted this turn
  int num_sphere_change = 0;
  for (int i = 0; i < is_sphere_deleted_this_turn.size(); i++) {
    if (is_sphere_deleted_this_turn[i]) num_sphere_change++;
  }
  printf("[Degenerate] deleted %d spheres this turn\n", num_sphere_change);
  return num_sphere_change;
}