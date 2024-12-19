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
  if (tan_pl2.fid == fid) return true;
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
  assert(this->fid >= 0);
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
      id, id_fe, id_fl, is_deleted, energy == DBL_MAX ? 1e10 : energy,
      energy_over_sq_radius == DBL_MAX ? 1e10 : energy_over_sq_radius,
      tan_point[0], tan_point[1], tan_point[2], normal[0], normal[1],
      normal[2]);
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

// NOTE: keep it same as FeatureEdge::operator==()
bool TangentConcaveLine::operator==(TangentConcaveLine const& b) const {
  if (id_fe == b.id_fe) return true;
  // belongs to the same grouped concave line
  // this line can be a curve!!!
  if (id_fl == b.id_fl) {
    // we need to check the segment direction
    // if the direction deviate too much, then we consider them as different
    if (!is_vector_same_direction(direction, b.direction, EPS_DEGREE_30) &&
        !is_vector_oppo_direction(direction, b.direction, EPS_DEGREE_30))
      return false;
    return true;
  }
  return false;
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
    // 1. surface covered has intersection
    // 2. tan_pl normal are similar to at least 1 tan_cc_line adjacent normals
    if (has_intersection(sf_fids_covered, one_tan_plane.sf_fids_covered) &&
        (angle_between_two_vectors_in_degrees(
             one_tan_plane.normal, this->adj_normals[0]) < EPS_DEGREE_30 ||
         angle_between_two_vectors_in_degrees(
             one_tan_plane.normal, this->adj_normals[1]) < EPS_DEGREE_30)) {
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
    assert(adj_fid >= 0);
    sf_mesh.collect_kring_neighbors_given_fid(k, adj_fid,
                                              this->sf_fids_covered_two[i]);
  }
}

/////////////////////////////////////////////////////////////////////////////////////
// MedialSphere
/////////////////////////////////////////////////////////////////////////////////////
MedialSphere::MedialSphere(int _id, Vector3 _pin, Vector3 _pin_normal,
                           int _pin_fid, SphereType _type, int _itr) {
  id = _id;
  is_rt_valid = false;
  is_rt_prev_valid = false;
  rt_change_status = RtValidStates::NoChange;
  site_id = -1;
  ss.p = _pin;
  ss.p_normal = _pin_normal;
  ss.p_fid = _pin_fid;
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
  if (is_on_se()) printf("se_line_id: %d\n", se_line_id);
  if (is_on_corner()) print_set<int>(corner_fls, "corner_fls");
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

  // printf("rt_neigh_ids_prev: [");
  // for (const auto rt_neigh_id_prev : rt_neigh_ids_prev)
  //   printf("%d, ", rt_neigh_id_prev);
  // printf("]\n");
  // printf("rt_neigh_ids: [");
  // for (const auto rt_neigh_id : rt_neigh_ids) printf("%d, ", rt_neigh_id);
  // printf("]\n");
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

void MedialSphere::copy(const MedialSphere& b) {
  this->center = b.center;
  this->radius = b.radius;
  this->type = b.type;
  this->tan_planes = b.tan_planes;
  this->tan_cc_lines = b.tan_cc_lines;
  this->se_edge_id = b.se_edge_id;
  this->se_line_id = b.se_line_id;
  this->corner_fls = b.corner_fls;
  this->corner_fes = b.corner_fes;
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
  if (type == SphereType::T_N || type == SphereType::T_N_c ||
      type == SphereType::T_N_JUNC)
    return true;
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
bool MedialSphere::is_on_junction() const {
  if (type == SphereType::T_N_JUNC) return true;
  return false;
}

bool MedialSphere::is_on_seam_endpoint() const {
  if (is_on_corner()) return true;
  if (is_on_junction()) return true;
  return false;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
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

bool is_two_mspheres_on_same_sl_including_corners(
    const MedialSphere& msphere1, const MedialSphere& msphere2) {
  if (!msphere1.is_on_extf() || !msphere2.is_on_extf()) return false;
  if (msphere1.is_on_se() && msphere2.is_on_se()) {
    if (msphere1.se_line_id != msphere2.se_line_id) return false;
    return true;
  }
  if (msphere1.is_on_corner() && msphere2.is_on_corner()) {
    if (has_intersection(msphere1.corner_fes, msphere2.corner_fes)) return true;
    return false;
  }
  if (msphere1.is_on_se() && msphere2.is_on_corner()) {
    if (msphere2.corner_fls.find(msphere1.se_line_id) !=
        msphere2.corner_fls.end())
      return true;
    return false;
  }
  if (msphere1.is_on_corner() && msphere2.is_on_se()) {
    if (msphere1.corner_fls.find(msphere2.se_line_id) !=
        msphere1.corner_fls.end())
      return true;
    return false;
  }
  return false;
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

// std::vector<ConvexCellHost> merge_convex_cells(
//     const std::vector<MedialSphere>& all_medial_spheres,
//     const std::set<int>& spheres_and_1rings,
//     const std::vector<ConvexCellHost>& convex_cells_prev,
//     const std::vector<ConvexCellHost>& convex_cells_new, bool is_debug) {
//   std::vector<ConvexCellHost> merged_convex_cells;
//   for (const auto& cell_new : convex_cells_new) {
//     // filter by spheres + 1ring, No 2ring
//     if (spheres_and_1rings.find(cell_new.voro_id) ==
//     spheres_and_1rings.end())
//       continue;

//     merged_convex_cells.push_back(cell_new);
//     merged_convex_cells.back().id = merged_convex_cells.size() - 1;

//     // if (cell_new.voro_id == 377)
//     //   printf("[Merge ConvexCells] cell_new %d has voro_id %d \n",
//     //          merged_convex_cells.back().id,
//     //          merged_convex_cells.back().voro_id);
//   }

//   for (const auto& cell_prev : convex_cells_prev) {
//     // if prev voro_id became invalid, then skip

//     // do not store cells of (spheres + 1ring)
//     if (spheres_and_1rings.find(cell_prev.voro_id) !=
//     spheres_and_1rings.end())
//       continue;
//     // store prev
//     merged_convex_cells.push_back(cell_prev);
//     merged_convex_cells.back().id = merged_convex_cells.size() - 1;

//     // if (cell_prev.voro_id == 377) {
//     //   printf("[Merge ConvexCells] cell_prev %d has voro_id %d, tet_id:
//     %d\n",
//     //          merged_convex_cells.back().id,
//     //          merged_convex_cells.back().voro_id,
//     //          merged_convex_cells.back().tet_id);
//     //   if (merged_convex_cells.back().tet_id == 900) {
//     //     merged_convex_cells.back().print_info();
//     //   }
//     // }
//   }

//   printf("[Merge ConvexCells] #prev: %zu, #new: %zu, #merged: %zu\n",
//          convex_cells_prev.size(), convex_cells_new.size(),
//          merged_convex_cells.size());
//   return merged_convex_cells;
// }

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