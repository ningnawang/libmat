#pragma once
#include <assert.h>
#include <geogram/basic/geometry.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/basic/matrix.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>

#include <array>
#include <ctime>
#include <fstream>
#include <queue>
#include <string>
#include <unordered_set>

#include "common_cxx.h"

typedef GEO::vec2 Vector2;
typedef GEO::vec3 Vector3;  // all coordinates are initialized to 0 (zero).
typedef GEO::vec4 Vector4;
typedef GEO::vec3i Vector3i;
typedef GEO::mat2 Matrix2;
typedef GEO::mat3 Matrix3;  //  initializes the matrix to the identity matrix
typedef GEO::mat4 Matrix4;
typedef std::pair<Vector3, int>
    v2int;  // mapping from vertex to one sf_mesh fid
typedef std::pair<Vector3, aint2> v2int2;
typedef std::array<Vector3, 2> avec2;
struct avec2int {
  Vector3 p;
  Vector3 n;
  int f;
};

inline Vector3 to_vec(const cfloat4& v) { return Vector3(v.x, v.y, v.z); }

inline Vector3 get_triangle_centroid(const Vector3& v1, const Vector3& v2,
                                     const Vector3& v3) {
  return (v1 + v2 + v3) / 3.;
};

// Return unit normal given 3 points
inline Vector3 get_normal(const Vector3& a, const Vector3& b,
                          const Vector3& c) {
  // facet abc is clockwise/counter-clockwise
  // (a-c) x (b-c) is clockwise/counter-clockwise
  return GEO::normalize(GEO::cross(a - c, b - c));
}

// unit vector from a to b
inline Vector3 get_direction(const Vector3& a, const Vector3& b) {
  return GEO::normalize(b - a);
}

inline void get_random_vector_between_two_vectors(const Vector3& nA,
                                                  const Vector3& nB,
                                                  Vector3& nX) {
  // double rand_k = ((double)std::rand() / (RAND_MAX));
  Scalar rand_k = RANDOM_01();
  nX = GEO::normalize(rand_k * nA + (1. - rand_k) * nB);
  // logger().debug("rand_k: {}, nA ({},{},{}), nB: ({},{},{}), nX: ({},{},{})",
  // 	rand_k, nA[0], nA[1], nA[2], nB[0], nB[1], nB[2], nX[0], nX[1], nX[2]);
}

inline double angle_between_two_vectors_in_degrees(const Vector3& a,
                                                   const Vector3& b) {
  // a and b should be normalized, so a.dot(b) in range [-1, 1]
  // However, with floating-point math, this is not necessarily true
  // so we need to clamp into [-1, 1]
  Vector3 a_normalized = GEO::normalize(a);
  Vector3 b_normalized = GEO::normalize(b);
  double ab_dot = GEO::dot(a_normalized, b_normalized);
  if (ab_dot <= -1.0) {
    return 180.;
  } else if (ab_dot >= 1.0) {
    return 0.;
  }
  // ab_dot in (-1, 1)
  double diff_angle = std::acos(ab_dot);
  // logger().debug("calculate angle ({},{},{}), ({},{},{}), diff_angle {}",
  //     a[0], a[1], a[2], b[0], b[1], b[2], diff_angle
  // );
  diff_angle *= (180. / PI);
  return diff_angle;
}

inline Vector3 sample_random_vector_given_two_vectors(const Vector3& nA,
                                                      const Vector3& nB) {
  double t_rand = RANDOM_01();
  while (t_rand <= 0 || t_rand >= 1) t_rand = RANDOM_01();
  return GEO::normalize((1. - t_rand) * nA + t_rand * nB);
}

// total k
// skip 2 boundary normals: nA and nB
inline void sample_k_vectors_given_two_vectors_perturb(
    const Vector3& nA, const Vector3& nB, const int k,
    std::vector<Vector3>& nXs) {
  nXs.clear();
  if (k < 1) return;
  // skip nA and nB
  double step = 1. / std::max(k + 1, 2);  // in (0, 1/2]
  // t in the range of [step, 1)
  // somehow not working fine when < 1.
  for (double t = step; t < 0.99999999999; t += step) {
    // // apply perturbation to t
    // double t_tmp = t + RANDOM_01() * (step / 5.);
    // while (t_tmp >= 1.) t_tmp = t + RANDOM_01() * (step / 5.);
    // while (t_tmp > 1.) {
    //   t_tmp = t_tmp - 1;
    // }
    // // printf("t: %f, t_tmp: %f \n", t, t_tmp);
    double t_tmp = t;  // no perturbation
    Vector3 nX = (1. - t_tmp) * nA + t_tmp * nB;
    // perturb nX to any direction, not only nA and nB
    int rand_idx = RANDOM_INT(0, 3);
    nX[rand_idx] = nX[rand_idx] + RANDOM_01() / 5;
    nX = GEO::normalize(nX);
    // nX[0] = nX[0] + RANDOM_01() / 10;
    // nX[1] = nX[1] + RANDOM_01() / 10;
    // nX[2] = nX[2] + RANDOM_01() / 10;
    nXs.push_back(nX);
  }
  if (nXs.size() == k - 1) {
    double t = RANDOM_01();
    Vector3 nX = GEO::normalize((1. - t) * nA + t * nB);
    nXs.push_back(nX);
  }
  if (k != nXs.size()) {
    printf("[ERROR] did not created k: %d normlas with step %f, created %ld \n",
           k, step, nXs.size());
    assert(false);
  }
}

// total k
inline void sample_k_vectors_given_N_vectors(const std::vector<Vector3>& n_vec,
                                             const int k,
                                             std::vector<Vector3>& nXs) {
  nXs.clear();
  if (k < 1) return;
  std::vector<double> n_t(n_vec.size() - 1, -1);
  double t_sum = 0.f;
  FOR(i, k) {
    do {
      t_sum = 0.f;
      FOR(j, n_t.size()) {
        n_t[j] = RANDOM_01();
        t_sum += n_t[j];
      }
    } while (t_sum <= 0. || t_sum >= 0.9999999999);  // t_sum = 1 is fine
    Vector3 nX(0.f, 0.f, 0.f);
    FOR(j, n_t.size()) nX += n_t[j] * n_vec[j];
    nX += (1 - t_sum) * n_vec[n_vec.size() - 1];
    nXs.push_back(GEO::normalize(nX));
  }
  if (k != nXs.size()) {
    printf("[ERROR] did not created k: %d normals, created %ld \n", k,
           nXs.size());
    assert(false);
  }
}

// exclude 2 boundary points: nA and nB
inline void sample_new_points_given_len(const Vector3& start_pos,
                                        const Vector3& end_pos,
                                        const double len_ideal,
                                        std::vector<Vector3>& new_points) {
  new_points.clear();
  if (len_ideal <= 0.) return;
  double len = GEO::distance(start_pos, end_pos);
  int step = std::floor(len / len_ideal);
  // add new points given len_ideal
  Vector3 dir = get_direction(start_pos, end_pos);
  for (int j = 1; j < step; j++) {
    Vector3 new_pos = start_pos + j * len_ideal * dir;
    new_points.push_back(new_pos);
  }
}

// Project p onto line defined by [v0, v1]
inline void project_point_onto_line(const Vector3& p, const Vector3& v0,
                                    const Vector3& v1, Vector3& p_proj,
                                    double& dist) {
  Vector3 v0p(p - v0), v0v1(v1 - v0);
  p_proj = GEO::dot(v0 + v0p, v0v1) / GEO::dot(v0v1, v0v1) * v0v1;
  dist = (p - p_proj).length();

  // TODO: using GEO::Geom::point_segment_squared_distance()
  // double lambda1, lambda2;
  // double sq_dist = GEO::Geom::point_segment_squared_distance(p, v0, v1,
  // p_proj, lambda1, lambda2);
  // dist = std::sqrt(sq_dist);
}

// Project p onto triangle defined by [v1, v2, v3]
inline void project_point_onto_triangle(const Vector3& p, const Vector3& v1,
                                        const Vector3& v2, const Vector3& v3,
                                        Vector3& p_proj, double& sq_dist) {
  double lambda1, lambda2, lambda3;
  sq_dist = GEO::Geom::point_triangle_squared_distance(
      p, v1, v2, v3, p_proj, lambda1, lambda2, lambda3);
}

// Project p onto a triangle of sf_mesh
inline void project_point_onto_triangle(const GEO::Mesh& sf_mesh, const int fid,
                                        const Vector3& p, Vector3& p_proj,
                                        double& sq_dist) {
  const Vector3 v0 = sf_mesh.vertices.point(sf_mesh.facets.vertex(fid, 0));
  const Vector3 v1 = sf_mesh.vertices.point(sf_mesh.facets.vertex(fid, 1));
  const Vector3 v2 = sf_mesh.vertices.point(sf_mesh.facets.vertex(fid, 2));
  double lambda1, lambda2, lambda3;
  sq_dist = GEO::Geom::point_triangle_squared_distance(
      p, v0, v1, v2, p_proj, lambda1, lambda2, lambda3);
}

// Project p onto a set of triangles of sf_mesh
// return the fid that projected on
inline int project_point_onto_triangles(const GEO::Mesh& sf_mesh,
                                        const std::set<int>& fids,
                                        const Vector3& p, Vector3& p_proj) {
  double min_sq_dist = DBL_MAX, sq_dist = DBL_MAX;
  Vector3 p_proj_tmp;
  int fid_return = -1;
  for (const int fid : fids) {
    project_point_onto_triangle(sf_mesh, fid, p, p_proj_tmp, sq_dist);
    if (sq_dist < min_sq_dist) {
      min_sq_dist = sq_dist;
      p_proj = p_proj_tmp;
      fid_return = fid;
    }
  }
  return fid_return;
}

// Return normalized vector
inline Vector3 get_mesh_facet_normal(const GEO::Mesh& mesh, const int fidx) {
  Vector3 fn = GEO::normalize(GEO::Geom::mesh_facet_normal(mesh, fidx));
  return fn;
}

inline Vector3 get_mesh_facet_centroid(const GEO::Mesh& mesh, const int fid) {
  assert(fid >= 0 || fid < mesh.facets.nb());
  return get_triangle_centroid(mesh.vertices.point(mesh.facets.vertex(fid, 0)),
                               mesh.vertices.point(mesh.facets.vertex(fid, 1)),
                               mesh.vertices.point(mesh.facets.vertex(fid, 2)));
}

// Return normalized vector
inline Vector3 get_mesh_vertex_normal(const GEO::Mesh& mesh, const int vidx) {
  Vector3 vn = GEO::normalize(GEO::Geom::mesh_vertex_normal(mesh, vidx));
  return vn;
}

inline Vector3 get_random_point_given_facet(const Vector3& p1,
                                            const Vector3& p2,
                                            const Vector3& p3) {
  // https://stackoverflow.com/a/21722167
  Scalar r1 = RANDOM_01();
  Scalar r2 = RANDOM_01();
  if (r1 + r2 > 1.f) {
    r1 = 1.f - r1;
    r2 = 1.f - r2;
  }
  std::array<Scalar, 3> w = {{1 - r1 - r2, r1, r2}};
  return w[0] * p1 + w[1] * p2 + w[2] * p3;
}

inline Vector3 get_random_point_given_edge(const Vector3& p1,
                                           const Vector3& p2) {
  // https://stackoverflow.com/a/21722167
  Scalar r = RANDOM_01();
  std::array<Scalar, 2> w = {{1 - r, r}};
  return w[0] * p1 + w[1] * p2;
}

inline bool all_finite(const Vector3& p) {
  for (uint i = 0; i < 3; i++) {
    if (std::isnan(p[i]) || !std::isfinite(p[i])) return false;
  }
  return true;
}

// make sure f1_normal and f2_normal are normalized!!
inline bool is_share_a_sharp_edge(const Vector3& f1_normal,
                                  const Vector3& f2_normal) {
  double cos = GEO::dot(f1_normal, f2_normal);
  // sharp edge: theta in [90, 180), cos in (-1, 0]
  if (cos > -1 && cos < SCALAR_ZERO_3) {
    return true;  // skip checking its neighbor
  }
  return false;
}

inline bool is_vector_same_direction(const Vector3& a, const Vector3& b,
                                     const double degree) {
  // angle betwen a and b is in [0, degree]
  // const double halfC = M_PI / 180;HALF_PI
  return GEO::dot(a, b) / (a.length() * b.length()) >=
         std::cos(degree * HALF_PI);
}

// check if x is within range of (a, b) vectors
inline bool is_vector_within_range_of_two_vectors(const Vector3& a,
                                                  const Vector3& b,
                                                  const Vector3& x) {
  // https://stackoverflow.com/a/43384516
  Vector3 ab = GEO::cross(a, b);  // right system, angle from a to b
  Vector3 ax = GEO::cross(a, x);
  Vector3 xb = GEO::cross(x, b);
  if (is_vector_same_direction(ab, ax, 90) &&
      is_vector_same_direction(ab, xb, 90))
    return true;
  return false;
}

inline Vector3 get_centroid(const std::vector<Vector3>& vs) {
  Vector3 c(0., 0., 0.);
  if (vs.size() == 0) return c;
  for (int i = 0; i < vs.size(); i++) c += vs[i];
  c /= vs.size();
  return c;
};

inline double get_distance_between_two_vectors(const Vector3& a,
                                               const Vector3& b) {
  return (a - b).length();
}

// perturbation of the coordinates
// to avoid boundary error (degenerations) while calculating RPD
inline void apply_perturb(Vector3& pos) {
  // int i = RANDOM_INT(0, 3);
  for (int i = 0; i < 3; i++) pos[i] += RANDOM_01() / 10.;
}
inline Vector3 apply_perturb(const Vector3& pos) {
  Vector3 pos_per = pos;
  // int i = RANDOM_INT(0, 3);
  for (int i = 0; i < 3; i++) pos_per[i] += RANDOM_01() / 10.;
  return pos_per;
}
// perturb pos towards pert_dir with len_perturb
inline Vector3 apply_perturb_given_dir_max(const Vector3& pos,
                                           const Vector3& pert_dir,
                                           const double pert_max) {
  Vector3 pos_per = pos;
  double len_perturb = RANDOM_01() * pert_max;
  // printf("len_perturb: %f \n", len_perturb);
  pos_per += len_perturb * pert_dir;
  return pos_per;
}

/**
 * @brief Compute A * transpose(B)
 *
 * @param A 3x3 vector
 * @param B 3x3 vector
 * @return Matrix3
 */
inline Matrix3 vec_vec_trans(Vector3 A, Vector3 B) {
  std::array<double, 9> result = {{A[0] * B[0], A[0] * B[1], A[0] * B[2],
                                   A[1] * B[0], A[1] * B[1], A[1] * B[2],
                                   A[2] * B[0], A[2] * B[1], A[2] * B[2]}};
  return Matrix3(result.data());
}

inline Vector3 convert_std2geo(const adouble3& p) {
  return Vector3(p[0], p[1], p[2]);
}