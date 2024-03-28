#include <math.h>

#include <algorithm>

#include "convex_cell.h"
// #include "stopwatch.h"

// ###################  ConvexCell   ######################
// niwnang: no use
__device__ ConvexCell::ConvexCell(const int p_seed, const float* p_pts,
                                  const size_t pitch, Status* p_status)
    : pts_pitch(pitch), pts(p_pts) {
  float eps = .1f;
  float xmin = -eps;
  float ymin = -eps;
  float zmin = -eps;
  float xmax = 1000 + eps;
  float ymax = 1000 + eps;
  float zmax = 1000 + eps;
  first_boundary_ = END_OF_LIST;
  FOR(i, _MAX_P_) boundary_next(i) = END_OF_LIST;
  voro_id = p_seed;
  if (pts_pitch)
    voro_seed = {pts[voro_id], pts[voro_id + pts_pitch],
                 pts[voro_id + (pts_pitch << 1)], 1};
  else
    voro_seed = {pts[3 * voro_id], pts[3 * voro_id + 1], pts[3 * voro_id + 2],
                 1};
  status = p_status;
  *status = success;

  clip(0) = make_float7(1.0, 0.0, 0.0, -xmin, F_TET_ADJ_DEFAULT, UNK_FACE_ID,
                        UNK_FACE_ID);
  clip(1) = make_float7(-1.0, 0.0, 0.0, xmax, F_TET_ADJ_DEFAULT, UNK_FACE_ID,
                        UNK_FACE_ID);
  clip(2) = make_float7(0.0, 1.0, 0.0, -ymin, F_TET_ADJ_DEFAULT, UNK_FACE_ID,
                        UNK_FACE_ID);
  clip(3) = make_float7(0.0, -1.0, 0.0, ymax, F_TET_ADJ_DEFAULT, UNK_FACE_ID,
                        UNK_FACE_ID);
  clip(4) = make_float7(0.0, 0.0, 1.0, -zmin, F_TET_ADJ_DEFAULT, UNK_FACE_ID,
                        UNK_FACE_ID);
  clip(5) = make_float7(0.0, 0.0, -1.0, zmax, F_TET_ADJ_DEFAULT, UNK_FACE_ID,
                        UNK_FACE_ID);
  nb_p = 6;

  ver(0) = make_uchar4(2, 5, 0, UNK_UCHAR);
  ver(1) = make_uchar4(5, 3, 0, UNK_UCHAR);
  ver(2) = make_uchar4(1, 5, 2, UNK_UCHAR);
  ver(3) = make_uchar4(5, 1, 3, UNK_UCHAR);
  ver(4) = make_uchar4(4, 2, 0, UNK_UCHAR);
  ver(5) = make_uchar4(4, 0, 3, UNK_UCHAR);
  ver(6) = make_uchar4(2, 4, 1, UNK_UCHAR);
  ver(7) = make_uchar4(4, 3, 1, UNK_UCHAR);
  nb_v = 8;
}

// ninwang: no use
__device__ ConvexCell::ConvexCell(
    const int p_seed, const float* p_pts, const size_t pitch, Status* p_status,
    const int tid, const float* vert /*tet vertices*/, const size_t vert_pitch,
    const int* idx /*tet indices*/, const size_t idx_pitch)
    : pts_pitch(pitch), pts(p_pts) {
  first_boundary_ = END_OF_LIST;
  FOR(i, _MAX_P_) boundary_next(i) = END_OF_LIST;
  voro_id = p_seed;
  if (pts_pitch)
    voro_seed = {pts[voro_id], pts[voro_id + pts_pitch],
                 pts[voro_id + (pts_pitch << 1)], 1};
  else
    voro_seed = {pts[3 * voro_id], pts[3 * voro_id + 1], pts[3 * voro_id + 2],
                 1};
  status = p_status;
  *status = success;

  FOR(i, 4) {  // tet vertices
    float3 vertices[3];
    FOR(j, 3) {
      vertices[j] = {
          vert[idx[tid + (tet_faces_lvid[i][j] * idx_pitch)]],
          vert[idx[tid + (tet_faces_lvid[i][j] * idx_pitch)] + vert_pitch],
          vert[idx[tid + (tet_faces_lvid[i][j] * idx_pitch)] +
               (vert_pitch << 1)]};
    }
    clip(i) = make_float7(tri2plane(vertices), F_TET_ADJ_DEFAULT, UNK_FACE_ID,
                          UNK_FACE_ID);
  }
  nb_p = 4;

  ver(0) = make_uchar4(1, 3, 2, UNK_UCHAR);
  ver(1) = make_uchar4(0, 1, 2, UNK_UCHAR);  // wrong ?
  ver(2) = make_uchar4(0, 3, 1, UNK_UCHAR);
  ver(3) = make_uchar4(0, 2, 3, UNK_UCHAR);  // wrong ?
  nb_v = 4;
}

// ninwang: no use
__device__ ConvexCell::ConvexCell(const VoronoiCell& vc, Status* p_status) {
  first_boundary_ = END_OF_LIST;
  FOR(i, _MAX_P_) boundary_next(i) = END_OF_LIST;
  voro_id = -1;
  voro_seed = vc.voro_seed;
  nb_p = vc.nb_p;
  nb_v = vc.nb_v;
  FOR(i, nb_p)
  clip(i) =
      make_float7(vc.clip[i], F_TET_ADJ_DEFAULT, UNK_FACE_ID, UNK_FACE_ID);
  FOR(i, nb_v)
  ver(i) = make_uchar4(vc.ver[i], UNK_UCHAR);
  status = p_status;
  *status = vc.status;
}

// ###################  Classes   ######################
// ninwang: load #adjacent cells for cell's vertices and edges,
//          each cell is suppose to be a tet
__device__ ConvexCell::ConvexCell(
    const int p_seed, const float* p_pts, const size_t pitch,
    const float* p_weights, const uint* p_flags, Status* p_status,
    const int tid, const float* vert /*tet vertices*/, const int n_vert,
    const size_t vert_pitch, const int* idx /*tet indices*/,
    const size_t idx_pitch, const int* v_adjs, const int* e_adjs,
    const int* f_adjs, const int* f_ids)
    : pts_pitch(pitch), pts(p_pts), pts_weights(p_weights) {
  first_boundary_ = END_OF_LIST;
  FOR(i, _MAX_P_) boundary_next(i) = END_OF_LIST;
  voro_id = p_seed;
  if (pts_pitch)
    voro_seed = {pts[voro_id], pts[voro_id + pts_pitch],
                 pts[voro_id + (pts_pitch << 1)], 1};
  else
    voro_seed = {pts[3 * voro_id], pts[3 * voro_id + 1], pts[3 * voro_id + 2],
                 1};
  voro_seed.w = p_weights[voro_id];  // sq_radius
  voro_flag = p_flags[voro_id];
  tet_id = tid;
  status = p_status;
  *status = success;

  // load dual vertices
  FOR(i, 4) {  // tet vertices
    float3 vertices[3];
    FOR(j, 3) {
      vertices[j] = {
          vert[idx[tid + (tet_faces_lvid[i][j] * idx_pitch)]],
          vert[idx[tid + (tet_faces_lvid[i][j] * idx_pitch)] + vert_pitch],
          vert[idx[tid + (tet_faces_lvid[i][j] * idx_pitch)] +
               (vert_pitch << 1)]};
    }
    // some tet faces are on boundary so only #adj=1
    // store another array f_adj.size() = idx.size() -> 4 faces/vertices per tet
    // matching tet_faces_lvid, for storing #adj for each face
    //
    // For face's unique id inherit from tet (and sf_mesh), we store {id, -1}.
    // For halfplane/clippping plane, we store {seed_id_min, seed_id_max}
    clip(i) = make_float7(tri2plane(vertices), f_adjs[tid * 4 + i],
                          f_ids[tid * 4 + i], UNK_FACE_ID);
    // int vids[3];
    // FOR(j, 3) { vids[j] = idx[tid + (tet_faces_lvid[i][j] * idx_pitch)]; }
    // if (voro_id == 1 && tet_id == 2849)
    //   printf(
    //       "after creation, clip %d, fid: %d, vids: (%d, %d, %d), #adj: %f, "
    //       "(%f, %f, %f, %f) \n",
    //       i, tid * 4 + i, vids[0], vids[1], vids[2], clip(i).h, clip(i).x,
    //       clip(i).y, clip(i).z, clip(i).w);
  }
  nb_p = 4;

  // load dual triangle (clip1, clip2, clip3, #adjacent_cells)
  int tr_adjs[4] = {UNK_INT, UNK_INT, UNK_INT, UNK_INT};
  FOR(lvid, 4) {
    int vid = idx[tid + (lvid * idx_pitch)];
    assert(vid < n_vert && vid >= 0);
    if (tr_adjs[lvid] != UNK_INT) continue;
    tr_adjs[lvid] = v_adjs[vid];

    // if (voro_id == 0 && tid == 1)
    //   printf("after creation, tr %d, vid: %d, #adj: %d \n", lvid, vid,
    //          v_adjs[vid]);
  }

  // store 4 vertices as dual triangles
  // for each vertex, the triangle order matters (pointing outwards)
  // but the order of vertices does not matter
  // so (0,1,2) can either be ver(1) or ver(3)
  FOR(i, 4) assert(tr_adjs[i] != UNK_INT);
  ver(0) = make_uchar4(1, 3, 2, tr_adjs[0]);
  ver(1) = make_uchar4(0, 2, 3, tr_adjs[1]);  // yes!
  ver(2) = make_uchar4(0, 3, 1, tr_adjs[2]);
  ver(3) = make_uchar4(0, 1, 2, tr_adjs[3]);  // yes!
  nb_v = 4;

  // load edges, (clip1, clip2, #adjacent_cells)
  nb_e = 0;
  for (uint lv1 = 0; lv1 < 4; lv1++) {
    for (uint lv2 = lv1 + 1; lv2 < 4; lv2++) {
      int v1 = idx[tid + (lv1 * idx_pitch)];
      int v2 = idx[tid + (lv2 * idx_pitch)];
      assert(v1 < n_vert && v2 < n_vert);
      int e_adj = get_e_adj(e_adjs, n_vert, v1, v2);
      assert(e_adj != UNK_INT);
      // e.g. vertices (0,1) share faces/clip (1,2)
      if (lv1 == 0 && lv2 == 1) edge(nb_e) = make_uchar3(2, 3, e_adj);
      if (lv1 == 0 && lv2 == 2) edge(nb_e) = make_uchar3(1, 3, e_adj);
      if (lv1 == 0 && lv2 == 3) edge(nb_e) = make_uchar3(1, 2, e_adj);
      if (lv1 == 1 && lv2 == 2) edge(nb_e) = make_uchar3(0, 3, e_adj);
      if (lv1 == 1 && lv2 == 3) edge(nb_e) = make_uchar3(0, 2, e_adj);
      if (lv1 == 2 && lv2 == 3) edge(nb_e) = make_uchar3(0, 1, e_adj);
      nb_e++;
    }
  }
  assert(nb_e == 6);

  euler = UNK_FLOAT;
}

__device__ bool ConvexCell::is_security_radius_reached_legacy(float4 last_neig,
                                                              bool is_debug) {
  // finds furthest voro vertex distance2
  float v_dist = 0;
  FOR(i, nb_v) {
    float4 pc = compute_vertex_coordinates(ver3(i));
    float4 diff = minus4(pc, voro_seed);
    float d2 = dot3(diff, diff);  // TODO safe to put dot4 here, diff.w = 0
    v_dist = max(d2, v_dist);
  }
  // compare to new neighbors distance2
  float4 diff = minus4(
      last_neig,
      voro_seed);  // TODO it really should take index of the neighbor instead
                   // of the float4, then would be safe to put dot4
  float d2 = dot3(diff, diff);
  //
  // change the following line if you want to compute power diagram to
  return (d2 > 6 * v_dist);
  //
  // return (d2 > 4 * v_dist);
}

// taking weights into account
__device__ bool ConvexCell::is_security_radius_reached(float4 last_neig,
                                                       bool is_debug) {
  // finds furthest voro vertex distance2
  float v_dist = 0;
  FOR(i, nb_v) {
    float4 pc = compute_vertex_coordinates(ver3(i));
    float4 diff = minus4(pc, voro_seed);
    float d2 = dot3(diff, diff);  // TODO safe to put dot4 here, diff.w = 0
    v_dist = max(d2, v_dist);
  }
  // compare to 'halfplane point' to new neighbors distance2.
  // 'halfplane point' is the joint point of halfplane and
  // line (voro_seed, last_neig)
  float4 diff = minus4(voro_seed, last_neig);
  float r2_diff = voro_seed.w - last_neig.w;  // weight is already r^2
  float w = (dot3(diff, diff) - r2_diff) / (2 * dot3(diff, diff));
  float4 p_halfplane = plus4(mul4(w, diff), last_neig);  // don't care about w
  float4 voro_p_diff = minus4(p_halfplane, voro_seed);   // don't care about w
  float d2 = dot3(voro_p_diff, voro_p_diff);

  if (is_debug)
    printf(
        "[SecurityRadius] seed %d (%f,%f,%f,%f) and last_neig (%f,%f,%f,%f) "
        "has v_dist %f and d2 %f \n",
        voro_id, voro_seed.x, voro_seed.y, voro_seed.z, voro_seed.w,
        last_neig.x, last_neig.y, last_neig.z, last_neig.w, v_dist, d2);
  return (d2 > 4 * v_dist);
  // return (d2 > 6 * v_dist);
}

__device__ inline uchar& ConvexCell::ith_plane(uchar v, int i) {
  return reinterpret_cast<uchar*>(&(ver(v)))[i];
}

__device__ bool ConvexCell::is_vertex_perturb(uchar3 v) {
  float4 pi1 = clip4(v.x);
  float4 pi2 = clip4(v.y);
  float4 pi3 = clip4(v.z);
  float4 result;
  result.x =
      -det3x3(pi1.w, pi1.y, pi1.z, pi2.w, pi2.y, pi2.z, pi3.w, pi3.y, pi3.z);
  result.y =
      -det3x3(pi1.x, pi1.w, pi1.z, pi2.x, pi2.w, pi2.z, pi3.x, pi3.w, pi3.z);
  result.z =
      -det3x3(pi1.x, pi1.y, pi1.w, pi2.x, pi2.y, pi2.w, pi3.x, pi3.y, pi3.w);
  result.w =
      det3x3(pi1.x, pi1.y, pi1.z, pi2.x, pi2.y, pi2.z, pi3.x, pi3.y, pi3.z);

  float4 result_persp = make_float4(result.x / result.w, result.y / result.w,
                                    result.z / result.w, 1);

  // if (voro_id == 58 && tet_id == 1143 && v.x == 6 && v.y == 5 && v.z == 4) {
  //   printf(
  //       "[ERROR] tet_id %d, voro_id %d, vertex (%d,%d,%d) w: %f, xyz "
  //       "(%f,%f,%f)\n",
  //       tet_id, voro_id, v.x, v.y, v.z, result.w, result.x, result.y,
  //       result.z);
  //   printf("p1: (%f,%f,%f,%f)\n", pi1.x, pi1.y, pi1.z, pi1.w);
  //   printf("p2: (%f,%f,%f,%f)\n", pi2.x, pi2.y, pi2.z, pi2.w);
  //   printf("p3: (%f,%f,%f,%f)\n", pi3.x, pi3.y, pi3.z, pi3.w);
  //   printf("result_persp: (%f,%f,%f,%f)\n", result_persp.x, result_persp.y,
  //          result_persp.z, result_persp.w);
  // }

  // ninwang: this must be rare
  // if (w < SCALAR_ZERO_6 && w > -SCALAR_ZERO_6) {
  if (result.w == 0.f) {
    printf(
        "[ERROR] tet_id %d, voro_id %d, vertex (%d,%d,%d) w: %f, xyz "
        "(%f,%f,%f) needs_perturb!\n",
        tet_id, voro_id, v.x, v.y, v.z, result.w, result.x, result.y, result.z);
    *status = needs_perturb;
    return true;
  }

  return false;
}

// ninwang: must in sync with ConvexCellHost::compute_vertex_coordinates()
__device__ float4
ConvexCell::compute_vertex_coordinates(uchar3 v, bool persp_divide) const {
  float4 pi1 = clip4(v.x);
  float4 pi2 = clip4(v.y);
  float4 pi3 = clip4(v.z);
  float4 result;
  result.x =
      -det3x3(pi1.w, pi1.y, pi1.z, pi2.w, pi2.y, pi2.z, pi3.w, pi3.y, pi3.z);
  result.y =
      -det3x3(pi1.x, pi1.w, pi1.z, pi2.x, pi2.w, pi2.z, pi3.x, pi3.w, pi3.z);
  result.z =
      -det3x3(pi1.x, pi1.y, pi1.w, pi2.x, pi2.y, pi2.w, pi3.x, pi3.y, pi3.w);
  result.w =
      det3x3(pi1.x, pi1.y, pi1.z, pi2.x, pi2.y, pi2.z, pi3.x, pi3.y, pi3.z);

  // if (tet_id == 166 && voro_id == 5 && v.x == 5) {
  //   printf(
  //       "tet_id %d, voro_id %d, vertex (%d,%d,%d) result: (%f, %f, %f, %f), "
  //       "persp_divide: %d \n",
  //       tet_id, voro_id, v.x, v.y, v.z, result.x, result.y, result.z,
  //       result.w, persp_divide);
  // }

  // ninwang: this must be rare
  if (result.w == 0.f) {
    *status = needs_perturb;
  }

  if (persp_divide)
    return make_float4(result.x / result.w, result.y / result.w,
                       result.z / result.w, 1);
  return result;
}

__device__ bool ConvexCell::cc_vertex_is_in_conflict_float(
    uchar3 v, float4 eqn, bool is_debug) const {
  float4 pi1 = clip4(v.x);
  float4 pi2 = clip4(v.y);
  float4 pi3 = clip4(v.z);

  // check if new plane (=eqn) shares with pi1/pi2/pi3
  if ((eqn.x == pi1.x && eqn.y == pi1.y && eqn.z == pi1.z && eqn.w == pi1.w) ||
      (eqn.x == pi2.x && eqn.y == pi2.y && eqn.z == pi2.z && eqn.w == pi2.w) ||
      (eqn.x == pi3.x && eqn.y == pi3.y && eqn.z == pi3.z && eqn.w == pi3.w)) {
    // printf("[conflict] voro_id: %d, v: (%d, %d, %d), duplicated, no
    // conflict\n", voro_id, v.x, v.y, v.z);
    return false;  // no conflict
  }

  float det = det4x4(pi1.x, pi2.x, pi3.x, eqn.x, pi1.y, pi2.y, pi3.y, eqn.y,
                     pi1.z, pi2.z, pi3.z, eqn.z, pi1.w, pi2.w, pi3.w, eqn.w);

  // check if v too close to halfplane
  float4 voro_vertex = compute_vertex_coordinates(v);
  float3 dir = make_float3(eqn.x, eqn.y, eqn.z);  // plane normal
  float dirNorm = sqrt(dot3(dir, dir));
  float dist = dot4(voro_vertex, eqn);
  float dist_norm = dist / dirNorm;
  float dist_abs = fabsf(dist_norm);

  if (is_debug) {
    if (dist_norm > 0 && det > 0 || dist_norm < 0 && det < 0) {
      printf("++++++++++ ERROR: sign(dist_norm) %f == sign(det) %f \n",
             dist_norm, det);
    }
    printf(
        "[conflict] v:(%d,%d,%d), pos:(%f,%f,%f,%f), eqn:(%f,%f,%f,%f), "
        "dist: %f, dist_norm: %f, distNorm: %f, dist_abs: %lf\n",
        v.x, v.y, v.z, voro_vertex.x, voro_vertex.y, voro_vertex.z,
        voro_vertex.w, eqn.x, eqn.y, eqn.z, eqn.w, dist, dist_norm, dirNorm,
        dist_abs);
    if (dist_abs <= 0.01)
      printf(
          "[conflict OUT1] voro_id: %d, v: (%d, %d, %d) dist_abs: %lf, det: "
          "%lf \n",
          voro_id, v.x, v.y, v.z, dist_abs, det);
    else if (det > 0.0f)
      printf("[conflict IN2] voro_id: %d, v: (%d, %d, %d) det: %lf \n", voro_id,
             v.x, v.y, v.z, det);
    else
      printf("[conflict OUT2] voro_id: %d, v: (%d, %d, %d) det: %lf \n",
             voro_id, v.x, v.y, v.z, det);
    printf("[conflict] pi1: (%f, %f, %f, %f)\n", pi1.x, pi1.y, pi1.z, pi1.w);
    printf("[conflict] pi2: (%f, %f, %f, %f)\n", pi2.x, pi2.y, pi2.z, pi2.w);
    printf("[conflict] pi3: (%f, %f, %f, %f)\n", pi3.x, pi3.y, pi3.z, pi3.w);
    printf("[conflict] eqn: (%f, %f, %f, %f)\n", eqn.x, eqn.y, eqn.z, eqn.w);
  }

  // if (dist_norm > 0 && det > 0 || dist_norm < 0 && det < 0) {
  //   return false;
  // }
  // if (dist_abs < 0.01) {
  //   return false;
  // }

#ifdef USE_ARITHMETIC_FILTER
  float maxx = max4(fabsf(pi1.x), fabsf(pi2.x), fabsf(pi3.x), fabsf(eqn.x));
  float maxy = max4(fabsf(pi1.y), fabsf(pi2.y), fabsf(pi3.y), fabsf(eqn.y));
  float maxz = max4(fabsf(pi1.z), fabsf(pi2.z), fabsf(pi3.z), fabsf(eqn.z));

  // The constant is computed by the program
  // in predicate_generator/
  float eps = 6.6876506e-05 * maxx * maxy * maxz;

  float min_max;
  float max_max;
  get_minmax3(min_max, max_max, maxx, maxy, maxz);

  eps *= (max_max * max_max);

  if (fabsf(det) < eps) {
    *status = needs_exact_predicates;
  }
#endif

  return (det > 0.0f);
}

__device__ bool ConvexCell::cc_vertex_is_in_conflict_double(
    uchar3 v, float4 eqn_f, bool is_debug) const {
  float4 pi1_f = clip4(v.x);
  float4 pi2_f = clip4(v.y);
  float4 pi3_f = clip4(v.z);

  // check if new plane (=eqn) shares with pi1/pi2/pi3
  if ((eqn_f.x == pi1_f.x && eqn_f.y == pi1_f.y && eqn_f.z == pi1_f.z &&
       eqn_f.w == pi1_f.w) ||
      (eqn_f.x == pi2_f.x && eqn_f.y == pi2_f.y && eqn_f.z == pi2_f.z &&
       eqn_f.w == pi2_f.w) ||
      (eqn_f.x == pi3_f.x && eqn_f.y == pi3_f.y && eqn_f.z == pi3_f.z &&
       eqn_f.w == pi3_f.w)) {
    return false;  // no conflict
  }

  double4 eqn = make_double4(eqn_f.x, eqn_f.y, eqn_f.z, eqn_f.w);
  double4 pi1 = make_double4(pi1_f.x, pi1_f.y, pi1_f.z, pi1_f.w);
  double4 pi2 = make_double4(pi2_f.x, pi2_f.y, pi2_f.z, pi2_f.w);
  double4 pi3 = make_double4(pi3_f.x, pi3_f.y, pi3_f.z, pi3_f.w);

  double det = det4x4(pi1.x, pi2.x, pi3.x, eqn.x, pi1.y, pi2.y, pi3.y, eqn.y,
                      pi1.z, pi2.z, pi3.z, eqn.z, pi1.w, pi2.w, pi3.w, eqn.w);

  // if (is_debug) {
  //   float4 voro_vertex = compute_vertex_coordinates(v);
  //   printf("[conflict] v: (%d, %d, %d), pos (%f, %f, %f)\n", v.x, v.y, v.z,
  //          voro_vertex.x, voro_vertex.y, voro_vertex.z);
  //   if (det > 0.0f)
  //     printf("[conflict IN] voro_id: %d, v: (%d, %d, %d) det: %lf \n",
  //     voro_id,
  //            v.x, v.y, v.z, det);
  //   else
  //     printf("[conflict OUT] voro_id: %d, v: (%d, %d, %d) det: %lf \n",
  //     voro_id,
  //            v.x, v.y, v.z, det);
  //   printf("[conflict] pi1: (%f, %f, %f, %f)\n", pi1.x, pi1.y, pi1.z, pi1.w);
  //   printf("[conflict] pi2: (%f, %f, %f, %f)\n", pi2.x, pi2.y, pi2.z, pi2.w);
  //   printf("[conflict] pi3: (%f, %f, %f, %f)\n", pi3.x, pi3.y, pi3.z, pi3.w);
  //   printf("[conflict] eqn: (%f, %f, %f, %f)\n", eqn.x, eqn.y, eqn.z, eqn.w);
  // }

#ifdef USE_ARITHMETIC_FILTER
  double maxx = max4(fabs(pi1.x), fabs(pi2.x), fabs(pi3.x), fabs(eqn.x));
  double maxy = max4(fabs(pi1.y), fabs(pi2.y), fabs(pi3.y), fabs(eqn.y));
  double maxz = max4(fabs(pi1.z), fabs(pi2.z), fabs(pi3.z), fabs(eqn.z));

  // The constant is computed by the program
  // in predicate_generator/
  double eps = 1.2466136531027298e-13 * maxx * maxy * maxz;

  double min_max;
  double max_max;
  get_minmax3(min_max, max_max, maxx, maxy, maxz);

  eps *= (max_max * max_max);

  if (fabs(det) < eps) {
    *status = needs_exact_predicates;
  }
#endif

  return (det > 0.0f);
}

// pts_weights stores sq_radii (defaults 1)
__device__ float4 ConvexCell::point_from_index(int idx) {
  if (pts_pitch)
    return {pts[idx], pts[idx + pts_pitch], pts[idx + (pts_pitch << 1)],
            pts_weights[idx]};
  else
    return {pts[3 * idx], pts[3 * idx + 1], pts[3 * idx + 2], pts_weights[idx]};
}

__device__ int ConvexCell::find_edge_id(uchar clip1, uchar clip2) {
  uchar cmin = min(clip1, clip2);
  uchar cmax = max(clip1, clip2);
  for (uint eid = 0; eid < nb_e; eid++) {
    // yes, we loop them all, slow?
    const auto& e = edge(eid);
    if (e.x == cmin && e.y == cmax) {
      return eid;
    }
  }
  return -1;
}

__device__ void ConvexCell::new_vertex(uchar i, uchar j, uchar k) {
  if (nb_v + 1 >= _MAX_T_) {
    *status = triangle_overflow;
    return;
  }
  // find #adjacent cells
  int eid1 = find_edge_id(i, j);
  int eid2 = find_edge_id(i, k);
  int eid3 = find_edge_id(j, k);
  if (eid1 < 0 || eid2 < 0 || eid3 < 0)
    printf(
        "seed %d, tet_id %d, has clip (%d, %d, %d) has eid (%d, %d, %d), "
        "nb_e: "
        "%d \n",
        voro_id, tet_id, i, j, k, eid1, eid2, eid3, nb_e);
  assert(eid1 >= 0 && eid2 >= 0 && eid3 >= 0);
  uchar nb_adjs = max(max(edge(eid1).z, edge(eid2).z), edge(eid3).z);
  assert(nb_adjs > 0);
  // assert(nb_adjs == 4 || nb_adjs == 1 || nb_adjs == 2);

  ver(nb_v) = make_uchar4(i, j, k, nb_adjs);
  // if (voro_id == 2) {
  //   float4 voro_vertex =
  //   compute_vertex_coordinates(make_uchar3(ver(nb_v)));
  //   printf("[new_vertex] new v: %d, (%d, %d, %d), (%f, %f, %f) \n", nb_v,
  //   i, j,
  //          k, voro_vertex.x, voro_vertex.y, voro_vertex.z);
  // }

  nb_v++;
}

// power distance, take weight into account
// halfplane: ax+by+cz+d>0, we store (a,b,c,d)
// halfplane normal point to voro_seed
//
// Note: normal vector dir=(a,b,c) NO need to normalize!!
__device__ int ConvexCell::new_plane(int seed_id) {
  if (nb_p >= _MAX_P_) {
    *status = vertex_overflow;
    return -1;
  }

  // here we do not use plane_from_point_and_normal()
  // because it is hard to compute point on halfplane
  float4 B = point_from_index(seed_id);
  float4 dir = minus4(voro_seed, B);  // normal point to voro_seed
  float4 ave2 = plus4(voro_seed, B);
  // we add weights for power diagram
  float dot = dot3(ave2, dir) + (B.w - voro_seed.w);  // ninwang: should be + !!
  // if (voro_id == 5 && (seed_id == 19 || seed_id == 29) && tet_id == 166) {
  //   float d = -dot / 2.f;
  //   printf(
  //       "voro_id: %d (%f,%f,%f), seed_id: %d (%f,%f,%f), voro_seed.w: %f,
  //       B.w: "
  //       "%f, dot: %f, a: %f, b: %f, c: %f, d: %f\n ",
  //       voro_id, voro_seed.x, voro_seed.y, voro_seed.z, seed_id, B.x, B.y,
  //       B.z, voro_seed.w, B.w, dot, dir.x, dir.y, dir.z, d);
  // }
  int f_id_min = voro_id, f_id_max = seed_id;
  if (seed_id < voro_id) {
    f_id_min = seed_id;
    f_id_max = voro_id;
  }
  clip(nb_p) = make_float7(dir.x, dir.y, dir.z, -dot / 2.f, F_CELL_ADJ_DEFAULT,
                           f_id_min, f_id_max);
  nb_p++;
  return nb_p - 1;
}

__device__ int ConvexCell::new_plane(float4 eqn) {
  if (nb_p >= _MAX_P_) {
    *status = vertex_overflow;
    return -1;
  }
  clip(nb_p) = make_float7(eqn, F_CELL_ADJ_DEFAULT, UNK_FACE_ID, UNK_FACE_ID);
  nb_p++;
  return nb_p - 1;
}

__device__ void ConvexCell::new_edge(uchar clip1, uchar clip2) {
  if (nb_e >= _MAX_E_) {
    *status = edge_overflow;
    return;
  }

  uchar cmin = min(clip1, clip2);
  uchar cmax = max(clip1, clip2);
  float max_adjs = max(clip(clip1).h, clip(clip2).h);
  assert(max_adjs == 1 || max_adjs == 2);
  edge(nb_e) = make_uchar3(cmin, cmax, static_cast<uchar>(max_adjs));
  nb_e++;
}

__device__ void ConvexCell::compute_boundary() {
  // clean circular list of the boundary
  FOR(i, _MAX_P_) boundary_next(i) = END_OF_LIST;
  first_boundary_ = END_OF_LIST;

  int nb_iter = 0;
  uchar t = nb_v;
  while (nb_r > 0) {
    if (nb_iter++ > 65535) {
      *status = inconsistent_boundary;
      return;
    }
    bool is_in_border[3];
    bool next_is_opp[3];
    FOR(e, 3)
    is_in_border[e] = (boundary_next(ith_plane(t, e)) != END_OF_LIST);
    FOR(e, 3)
    next_is_opp[e] =
        (boundary_next(ith_plane(t, (e + 1) % 3)) == ith_plane(t, e));

    bool new_border_is_simple = true;
    // check for non manifoldness
    FOR(e, 3)
    if (!next_is_opp[e] && !next_is_opp[(e + 1) % 3] &&
        is_in_border[(e + 1) % 3])
      new_border_is_simple = false;

    // check for more than one boundary ... or first triangle
    if (!next_is_opp[0] && !next_is_opp[1] && !next_is_opp[2]) {
      if (first_boundary_ == END_OF_LIST) {
        FOR(e, 3)
        boundary_next(ith_plane(t, e)) = ith_plane(t, (e + 1) % 3);
        first_boundary_ = ver(t).x;
      } else
        new_border_is_simple = false;
    }

    if (!new_border_is_simple) {
      t++;
      if (t == nb_v + nb_r) t = nb_v;
      continue;
    }

    // link next
    FOR(e, 3)
    if (!next_is_opp[e])
      boundary_next(ith_plane(t, e)) = ith_plane(t, (e + 1) % 3);

    // destroy link from removed vertices
    FOR(e, 3) if (next_is_opp[e] && next_is_opp[(e + 1) % 3]) {
      if (first_boundary_ == ith_plane(t, (e + 1) % 3))
        first_boundary_ = boundary_next(ith_plane(t, (e + 1) % 3));
      boundary_next(ith_plane(t, (e + 1) % 3)) = END_OF_LIST;
    }

    // remove triangle from R, and restart iterating on R
    swap_cuda(ver(t), ver(nb_v + nb_r - 1));
    t = nb_v;
    nb_r--;
  }
}

__device__ void ConvexCell::clip_by_plane(int neigh_seed_id) {
  int cur_p = new_plane(neigh_seed_id);  // add new plane equation
  if (*status == vertex_overflow) return;
  float4 eqn = clip4(cur_p);
  nb_r = 0;

  bool is_debug = false;
  if (tet_id == 3175 && voro_id == 63 && neigh_seed_id == 864) is_debug = true;
  // if (is_debug) {  // for debug only
  //   printf(
  //       "[clip_by_plane2] voro_id: %d, neigh_idx: %d, eqn: "
  //       "(%f, %f, %f, %f)\n ",
  //       voro_id, neigh_seed_id, eqn.x, eqn.y, eqn.z, eqn.w);

  //   int i = 0;
  //   while (i < nb_v) {  // for all vertices of the cell
  //     cc_vertex_is_in_conflict(ver3(i), eqn, is_debug);
  //     printf(
  //         "[clip_by_plane2] voro_id: %d, vertex %d: (%d,%d,%d) checked its "
  //         "conflict\n",
  //         voro_id, i, ver(i).x, ver(i).y, ver(i).z);
  //     i++;
  //   }
  //   return;
  // }

  int i = 0;
  while (i < nb_v) {  // for all vertices of the cell
    if (cc_vertex_is_in_conflict(ver3(i), eqn, is_debug)) {
      if (is_debug) {
        // in conflict = closer to neigh_seed than current voro_seed
        printf(
            "[clip_by_plane1] voro_id: %d, vertex %d: (%d, %d, %d) in "
            "conflict\n",
            voro_id, i, ver(i).x, ver(i).y, ver(i).z);
      }
      nb_v--;
      swap_cuda(ver(i), ver(nb_v));
      nb_r++;
    } else
      i++;
  }

  // if (is_debug) {
  //   printf(
  //       "[clip_by_plane1] voro_id %d, tet_id: %d, neigh_seed_id: %d, eqn: "
  //       "(%f, %f, %f, %f), nb_r: %d \n",
  //       voro_id, tet_id, neigh_seed_id, eqn.x, eqn.y, eqn.z, eqn.w, nb_r);
  //   FOR(i, nb_v + nb_r) {
  //     float4 voro_vertex = compute_vertex_coordinates(make_uchar3(ver(i)));
  //     printf("[clip_by_plane1] v %d: (%d, %d, %d), #adj %d, (%f, %f,
  //     %f)\n", i,
  //            ver(i).x, ver(i).y, ver(i).z, ver(i).w, voro_vertex.x,
  //            voro_vertex.y, voro_vertex.z);
  //   }
  // }

  if (*status == needs_exact_predicates) {
    return;
  }

  if (nb_r == 0) {  // if no clips, then remove the plane equation
    nb_p--;
    return;
  }

  if (nb_v == 0) {
    *status = no_intersection;
    return;
  }

  // Step 2: compute cavity boundary
  compute_boundary();
  if (*status != success) return;
  if (first_boundary_ == END_OF_LIST) return;

  // Step 3 (ninwang): create edges, before new_vertex
  uchar cir = first_boundary_;
  do {
    new_edge(cur_p, cir);
    if (*status != success) return;  // have this!
    cir = boundary_next(cir);
  } while (cir != first_boundary_);

  // Step 4: Triangulate cavity, loop again
  cir = first_boundary_;
  do {
    uchar3 new_v = make_uchar3(cur_p, cir, boundary_next(cir));
    new_vertex(new_v.x, new_v.y, new_v.z);
    is_vertex_perturb(new_v);
    if (*status != success) return;

    cir = boundary_next(cir);
  } while (cir != first_boundary_);
}

__device__ void ConvexCell::clip_by_plane(float4 eqn) {
  int cur_p = new_plane(eqn);
  if (*status == vertex_overflow) return;
  nb_r = 0;

  int i = 0;
  while (i < nb_v) {  // for all vertices of the cell
    if (cc_vertex_is_in_conflict(ver3(i), eqn, false)) {
      nb_v--;
      swap_cuda(ver(i), ver(nb_v));
      nb_r++;
    } else
      i++;
  }

  if (*status == needs_exact_predicates) {
    return;
  }

  if (nb_r == 0) {  // if no clips, then remove the plane equation
    nb_p--;
    return;
  }

  if (nb_v == 0) {
    *status = no_intersection;
    return;
  }
  // Step 2: compute cavity boundary
  compute_boundary();
  if (*status != success && *status != security_radius_not_reached) return;
  if (first_boundary_ == END_OF_LIST) return;

  // Step 3 (ninwang): create edges, before new_vertex
  uchar cir = first_boundary_;
  do {
    new_edge(cur_p, cir);
    // if (*status != success) return;
    cir = boundary_next(cir);
  } while (cir != first_boundary_);

  // Step 4: Triangulate cavity, loop again
  cir = first_boundary_;
  do {
    new_vertex(cur_p, cir, boundary_next(cir));
    if (*status != success && *status != security_radius_not_reached) return;

    cir = boundary_next(cir);
  } while (cir != first_boundary_);
}

__device__ void ConvexCell::print_info() const {
  // TODO: euler always -1
  printf(
      "Cell of voro_id: %d, tet_id: %d, sq_radius: %f, nb_v: %d, nb_p: %d, "
      "nb_e: %d, status: %d, euler: %f\n",
      voro_id, tet_id, voro_seed.w, nb_v, nb_p, nb_e, *status, euler);
}

__device__ void ConvexCell::print_info_detail() const {
  // TODO: euler always -1
  printf(
      "Cell of voro_id: %d, tet_id: %d, sq_radius: %f, nb_v: %d, nb_p: %d, "
      "nb_e: %d, status: %d, euler: %f\n",
      voro_id, tet_id, voro_seed.w, nb_v, nb_p, nb_e, *status, euler);

  FOR(i, nb_p) {
    // printf("clip %d, #adj %f, seed1/fid1: %f, seed2/fid2: %f\n", i,
    // clip(i).h,
    //        clip(i).g, clip(i).k);
    printf("clip %d, a: %f, b: %f c: %f, d: %f\n", i, clip(i).x, clip(i).y,
           clip(i).z, clip(i).w);
  }

  FOR(i, nb_v) {
    printf("vertex %d: (%d, %d, %d), #adj %d\n", i, ver(i).x, ver(i).y,
           ver(i).z, ver(i).w);
  }

  FOR(i, nb_e) {
    printf("edge %d: (%d, %d), #adj %d\n", i, edge(i).x, edge(i).y, edge(i).z);
  }
}

__device__ void ConvexCell::print_vs() const {
  FOR(i, nb_v) {
    printf("vertex %d: (%d, %d, %d), #adj %d\n", i, ver(i).x, ver(i).y,
           ver(i).z, ver(i).w);
  }

  // FOR(i, nb_e) {
  //   printf("edge %d: (%d, %d), #adj %d\n", i, edge(i).x, edge(i).y,
  //   edge(i).z);
  // }
}

// TODO: no use, the active_edges need to
//       be collected in CPU
__device__ float ConvexCell::get_cell_euler() {
  float sum_v = 0.f, sum_f = 0.f, sum_e = 0.f;

  // vertices
  FOR(vid, nb_v) {
    int v_adj = static_cast<int>(ver(vid).w);
    assert(v_adj > 0);
    sum_v += 1.f / v_adj;
  }

  // some clipping planes may not exist in tri but
  // we still store it, here is to filter those planes
  int active_clipping_planes[_MAX_P_] = {0};
  FOR(vid, nb_v) {
    active_clipping_planes[ver(vid).x]++;
    active_clipping_planes[ver(vid).y]++;
    active_clipping_planes[ver(vid).z]++;
  }

  // faces
  FOR(fid, nb_p) {
    assert(fid < _MAX_P_);
    if (active_clipping_planes[fid] <= 0) continue;
    assert(clip(fid).h == F_CELL_ADJ_DEFAULT ||
           clip(fid).h == F_TET_ADJ_DEFAULT);
    sum_f += 1. / clip(fid).h;
  }

  // edges
  // here we do not filter those inactive edges caused
  // by two non-incident planes (used to be incident)
  // and keep the filter in the
  FOR(eid, nb_e) {
    uchar3 e = edge(eid);
    assert(e.x < _MAX_P_ && e.y < _MAX_P_);
    if (active_clipping_planes[e.x] <= 0 || active_clipping_planes[e.y] <= 0)
      continue;
    sum_e += 1. / e.z;
  }

  // TODO: not correct, inactive edges
  //       are not filterd
  euler = sum_v - sum_e + sum_f;

  // if (voro_id == 1 && sum_v > 1.90 && sum_v < 1.92)
  //   printf("seed %d has euler %f: sum_v %f, sum_e %f, sum_f %f \n",
  //   voro_id,
  //          euler, sum_v, sum_e, sum_f);
}

// ###################  ConvexCellTransfer   ######################
__host__ bool is_convex_cell_valid(const ConvexCellTransfer& cell) {
  if (cell.status != Status::success &&
      cell.status != Status::security_radius_not_reached)
    return false;
  return true;
}

// ###################  KERNEL   ######################
__device__ void copy(const ConvexCell& cc, ConvexCellTransfer& cc_trans,
                     const int threadid) {
  cc_trans.is_active = true;
  cc_trans.status = *cc.status;
  cc_trans.thread_id = threadid;
  cc_trans.voro_id = cc.voro_id;
  cc_trans.tet_id = cc.tet_id;
  cc_trans.euler = cc.euler;
  cc_trans.weight = cc.voro_seed.w;
  cc_trans.nb_v = cc.nb_v;
  cc_trans.nb_p = cc.nb_p;
  cc_trans.nb_e = cc.nb_e;
  FOR(i, cc.nb_v) cc_trans.ver_data_trans[i] = ver(i);
  FOR(i, cc.nb_p) cc_trans.clip_data_trans[i] = clip5(i);
  FOR(i, cc.nb_p) cc_trans.clip_id2_data_trans[i] = clip_id2(i);
  FOR(i, cc.nb_e) cc_trans.edge_data[i] = edge(i);
}

__device__ void compute_voro_cell(const float* site, const size_t site_pitch,
                                  const int* site_knn,
                                  const size_t site_knn_pitch, int seed,
                                  VoronoiCell& vc) {
  // create BBox
  ConvexCell cc(seed, site, site_pitch, &vc.status);

  for (uint v = 1 - (site_pitch == 0); v <= _K_ - (site_pitch == 0); ++v) {
    int idx = site_pitch ? seed + v * site_knn_pitch : _K_ * seed + v;
    int neigh_idx = site_knn[idx];

    cc.clip_by_plane(neigh_idx);
    if (cc.is_security_radius_reached(cc.point_from_index(neigh_idx))) {
      break;
    }
    if (vc.status != success) {
      return;
    }
  }
  // check security radius
  int last_idx =
      site_pitch ? seed + _K_ * site_knn_pitch : _K_ * (seed + 1) - 1;
  if (!cc.is_security_radius_reached(cc.point_from_index(site_knn[last_idx]))) {
    vc.status = security_radius_not_reached;
  }

  vc.nb_v = cc.nb_v;
  vc.nb_p = cc.nb_p;
  FOR(i, cc.nb_v)
  vc.ver[i] = ver3(i);
  FOR(i, cc.nb_p)
  vc.clip[i] = clip4(i);
  vc.voro_seed = cc.voro_seed;
}

__device__ void get_tet_decomposition_of_vertex(ConvexCell& cc, int t,
                                                float4* P) {
  float4 C = cc.voro_seed;
  float4 A = cc.compute_vertex_coordinates(ver3(t));
  FOR(i, 3) {
    // if (cc.voro_id == 1 && cc.tet_id == 2849) {
    //   float4 plane = clip4(cc.ith_plane(t, i));
    //   float4 n = make_float4(plane.x, plane.y, plane.z, 0);
    //   float n_2 = dot4(n, n);
    //   printf(
    //       "seed %d tet_id %d has clip %d: (%f,%f,%f,%f) with n:
    //       (%f,%f,%f,%f) " "and n_2: %f \n", cc.voro_id, cc.tet_id, i,
    //       plane.x, plane.y, plane.z, plane.w, n.x, n.y, n.z, n.w, n_2);
    // }
    P[2 * i] = project_on_plane(C, clip4(cc.ith_plane(t, i)));
  }
  FOR(i, 3)
  P[2 * i + 1] = project_on_plane(
      A, plane_from_point_and_normal(
             C, cross3(minus4(P[2 * i], C), minus4(P[(2 * (i + 1)) % 6], C))));
}

__device__ void atomic_add_bary_and_volume(ConvexCell& cc, int seed,
                                           float* cell_bary_sum,
                                           const size_t cell_bary_sum_pitch,
                                           float* cell_vol) {
  float4 tet_bary;
  float tet_vol;
  float4 bary_sum = {0.0f, 0.0f, 0.0f, 0.0f};
  float cur_cell_vol = 0;
  float4 P[6];
  float4 C = cc.voro_seed;

  FOR(t, cc.nb_v) {
    float4 A = cc.compute_vertex_coordinates(ver3(t));
    get_tet_decomposition_of_vertex(cc, t, P);
    FOR(i, 6) {
      // if (seed == 588 && cc.tet_id == 805)
      //   printf(
      //       "seed %d has tet_decom (%f,%f,%f,%f), (%f,%f,%f,%f), "
      //       "(%f,%f,%f,%f), (%f,%f,%f,%f) \n",
      //       seed, P[i].x, P[i].y, P[i].z, P[i].w, P[(i + 1) % 6].x,
      //       P[(i + 1) % 6].y, P[(i + 1) % 6].z, P[(i + 1) % 6].w, C.x, C.y,
      //       C.z, C.w, A.x, A.y, A.z, A.w);
      get_tet_volume_and_barycenter(tet_bary, tet_vol, P[i], P[(i + 1) % 6], C,
                                    A);
      bary_sum = plus4(bary_sum, mul3(tet_vol, tet_bary));
      cur_cell_vol += tet_vol;
      // if (seed == 588 && cc.tet_id == 805)
      //   printf("seed %d has tet_vol %f, tet_bary: (%f,%f,%f,%f)\n", seed,
      //          tet_vol, tet_bary.x, tet_bary.y, tet_bary.z, tet_bary.w);
    }
  }

  if (abs(cur_cell_vol) < 0.1) {
    *cc.status = no_intersection;
    return;
  }

  if (cell_bary_sum_pitch) {
    atomicAdd(cell_bary_sum + seed, bary_sum.x);
    atomicAdd(cell_bary_sum + seed + cell_bary_sum_pitch, bary_sum.y);
    atomicAdd(cell_bary_sum + seed + (cell_bary_sum_pitch << 1), bary_sum.z);
  } else {
    atomicAdd(cell_bary_sum + 3 * seed, bary_sum.x);
    atomicAdd(cell_bary_sum + 3 * seed + 1, bary_sum.y);
    atomicAdd(cell_bary_sum + 3 * seed + 2, bary_sum.z);
  }
  atomicAdd(cell_vol + seed, cur_cell_vol);

  // if (seed == 588) {
  //   float3 cell_bary;
  //   cell_bary.x = bary_sum.x / cur_cell_vol;
  //   cell_bary.y = bary_sum.y / cur_cell_vol;
  //   cell_bary.z = bary_sum.z / cur_cell_vol;
  //   printf(
  //       "seed %d has tet_id: %d, bary_sum (%f,%f,%f,%f), cur_cell_vol %f, "
  //       "cell_bary: (%f,%f,%f), cell_bary_sum: (%f,%f,%f), cell_vol: %f\n",
  //       seed, cc.tet_id, bary_sum.x, bary_sum.y, bary_sum.z, bary_sum.w,
  //       cur_cell_vol, cell_bary.x, cell_bary.y, cell_bary.z,
  //       cell_bary_sum[seed], cell_bary_sum[seed + cell_bary_sum_pitch],
  //       cell_bary_sum[seed + (cell_bary_sum_pitch << 1)], cell_vol[seed]);
  // }
}

__device__ void choose_random_bary_and_volume(ConvexCell& cc, int seed,
                                              float* cell_bary_sum,
                                              const size_t cell_bary_sum_pitch,
                                              float* cell_vol) {
  float4 tet_bary;
  float tet_vol;
  float4 bary_sum = {0.0f, 0.0f, 0.0f, 0.0f};
  float cur_cell_vol = 0;
  float4 P[6];
  float4 C = cc.voro_seed;

  FOR(t, cc.nb_v) {
    float4 A = cc.compute_vertex_coordinates(ver3(t));
    get_tet_decomposition_of_vertex(cc, t, P);
    FOR(i, 6) {
      // if (seed == 588 && cc.tet_id == 805)
      //   printf(
      //       "seed %d has tet_decom (%f,%f,%f,%f), (%f,%f,%f,%f), "
      //       "(%f,%f,%f,%f), (%f,%f,%f,%f) \n",
      //       seed, P[i].x, P[i].y, P[i].z, P[i].w, P[(i + 1) % 6].x,
      //       P[(i + 1) % 6].y, P[(i + 1) % 6].z, P[(i + 1) % 6].w, C.x, C.y,
      //       C.z, C.w, A.x, A.y, A.z, A.w);
      get_tet_volume_and_barycenter(tet_bary, tet_vol, P[i], P[(i + 1) % 6], C,
                                    A);
      bary_sum = plus4(bary_sum, mul3(tet_vol, tet_bary));
      cur_cell_vol += tet_vol;
      // if (seed == 588 && cc.tet_id == 805)
      //   printf("seed %d has tet_vol %f, tet_bary: (%f,%f,%f,%f)\n", seed,
      //          tet_vol, tet_bary.x, tet_bary.y, tet_bary.z, tet_bary.w);
    }
  }

  // if (abs(cur_cell_vol) < 0.1) {
  if (abs(cur_cell_vol) < 5) {
    *cc.status = no_intersection;
    return;
  }

  // TODO: use cuda curand_init()
  cell_bary_sum[seed] = bary_sum.x;
  cell_bary_sum[seed + cell_bary_sum_pitch] = bary_sum.y;
  cell_bary_sum[seed + (cell_bary_sum_pitch << 1)] = bary_sum.z;
  cell_vol[seed] = cur_cell_vol;

  // if (seed == 588) {
  //   float3 cell_bary;
  //   cell_bary.x = bary_sum.x / cur_cell_vol;
  //   cell_bary.y = bary_sum.y / cur_cell_vol;
  //   cell_bary.z = bary_sum.z / cur_cell_vol;
  //   printf(
  //       "seed %d has tet_id: %d, bary_sum (%f,%f,%f,%f), cur_cell_vol %f, "
  //       "cell_bary: (%f,%f,%f), cell_bary_sum: (%f,%f,%f), cell_vol: %f\n",
  //       seed, cc.tet_id, bary_sum.x, bary_sum.y, bary_sum.z, bary_sum.w,
  //       cur_cell_vol, cell_bary.x, cell_bary.y, cell_bary.z,
  //       cell_bary_sum[seed], cell_bary_sum[seed + cell_bary_sum_pitch],
  //       cell_bary_sum[seed + (cell_bary_sum_pitch << 1)], cell_vol[seed]);
  // }
}

// Warning: this method does not assign unique id to new plane
__device__ void clip_voro_cell_by_tet(ConvexCell& cc, int tid, float* vert,
                                      size_t vert_pitch, int* idx,
                                      size_t idx_pitch) {
  FOR(i, 4) {
    float3 vertices[3];
    FOR(j, 3) {
      vertices[j] = {
          vert[idx[tid + (tet_faces_lvid[i][j] * idx_pitch)]],
          vert[idx[tid + (tet_faces_lvid[i][j] * idx_pitch)] + vert_pitch],
          vert[idx[tid + (tet_faces_lvid[i][j] * idx_pitch)] +
               (vert_pitch << 1)]};
    }

    cc.clip_by_plane(tri2plane(vertices));

    if (*cc.status == no_intersection) return;
  }
}

//----------------------------------KERNEL
__global__ void voro_cell_test_GPU_param(const float* site, const int n_site,
                                         const size_t site_pitch,
                                         const int* site_knn,
                                         const size_t site_knn_pitch,
                                         VoronoiCell* voronoi_cells) {
  int seed = blockIdx.x * blockDim.x + threadIdx.x;
  if (seed >= n_site) return;
  compute_voro_cell(site, site_pitch, site_knn, site_knn_pitch, seed,
                    voronoi_cells[seed]);
}

// each thread will clip a ConvexCell
// and each Voronoi cell of one site may contains multiple ConvexCells
//
// e_features: UNK_INT(-1) / SHARP_EDGE(1) / CONCAVE_EDGE(2)
__global__ void clipped_voro_cell_test_GPU_param_tet(
    const float* site, const int n_site, const size_t site_pitch /*always>0*/,
    const float* site_weights, const uint* site_flags, const int* site_knn,
    const size_t site_knn_pitch, const int site_k, const float* vert,
    const int n_vert, const size_t vert_pitch, const int* idx, const int n_tet,
    const size_t idx_pitch, const int* v_adjs, const int* e_adjs,
    const int* f_adjs, const int* f_ids, const int* tet_knn,
    const size_t tet_knn_pitch, const int tet_k, Status* gpu_stat,
    VoronoiCell* voronoi_cells, ConvexCellTransfer* convex_cells_dev,
    float* cell_bary_sum, const size_t cell_bary_sum_pitch, float* cell_vol) {
  bool is_debug = false;
  FOR(i, n_vert) { assert(v_adjs[i] > 0); }

  int thread = blockIdx.x * blockDim.x + threadIdx.x;  // thread id
  int tid = thread / tet_k;                            // tet index

  // we define n_grids = n_tet * tet_k / VORO_BLOCK_SIZE + 1;
  // so tid may == n_tet, this would cause problem later
  if (tid >= n_tet) {
    if (is_debug)
      printf("[clipped] tid %d out of bound [0, %d) \n", tid, n_tet);
    // to avoid random value assigned to thread
    convex_cells_dev[thread].status = Status::early_return;
    return;
  }
  int seed = tet_knn[(thread % tet_k) * tet_knn_pitch + tid];
  // happens when #(tet_related_spheres) < tet_k
  // and we init as -1
  if (seed < 0 || seed >= n_site) {
    if (is_debug)
      printf("[clipped] seed %d out of bound [0, %d) \n", seed, n_site);
    // to avoid random value assigned to thread
    convex_cells_dev[thread].status = Status::early_return;
    return;
  }

  if (site_flags[seed] == SiteFlag::no_flag) {
    printf("[clipped] seed %d not partial spheres, no clip\n", seed);
    // to avoid random value assigned to thread
    convex_cells_dev[thread].status = Status::early_return;
    return;
  }

  // if (tid == 4 && seed == 4) is_debug = true;
  ConvexCell cc(seed, site, site_pitch, site_weights, site_flags,
                &(gpu_stat[thread]), tid, vert, n_vert, vert_pitch, idx,
                idx_pitch, v_adjs, e_adjs, f_adjs, f_ids);

  if (is_debug) {
    printf("[clipped] processing tid: %d, seed: %d\n", tid, seed);
    cc.print_info();
  }

  // ninwang:
  // The j-th row of column 'seed' in 'site_k x n_site' matrix is the j-th
  // random neighbor of sphere 'seed' from RT (not sorted).
  //
  // Since #neighbors <= site_k, filled with -1 if #neighbors < site_k,
  // so the real last_row is the last non-(-1) row
  int last_row = site_k - 1;  // will be updated later
  // ninwang: start with first random neighbor
  uint v_begin = 0;
  uint v_end = last_row;

  // // step 1: clip CE pin neighboring spheres first
  // //         also update last_row
  // for (uint v = v_begin; v <= v_end; ++v) {
  //   int idx = seed + v * site_knn_pitch;
  //   int neigh_idx = site_knn[idx];
  //   // ninwang: for RT neighbors, reached the last
  //   if (neigh_idx == -1) {
  //     last_row = v == 0 ? 0 : v - 1;
  //     break;
  //   }
  //   if (neigh_idx[site_flags] != SiteFlag::is_ce_pin) continue;
  //   cc.clip_by_plane(neigh_idx);
  //   if (is_debug) {
  //     printf(
  //         "[clipped1] v: %d, tid: %d, seed: %d, neigh_idx: %d, n_site: %d
  //         \n", v, tid, seed, neigh_idx, n_site);
  //     cc.print_info();
  //   }
  //   if (*cc.status != success) {
  //     return;
  //   }
  // }

  // step 2: clip other spheres
  //         last_row is updated
  for (uint v = v_begin; v <= last_row; ++v) {
    int idx = seed + v * site_knn_pitch;
    int neigh_idx = site_knn[idx];
    // ninwang: for RT neighbors, reached the last
    if (neigh_idx == -1) break;
    cc.clip_by_plane(neigh_idx);
    if (is_debug) {
      printf(
          "[clipped2 done] v: %d, tid: %d, seed: %d, neigh_idx: %d, n_site: "
          "%d, status %d \n",
          v, tid, seed, neigh_idx, n_site, *cc.status);

      // if (neigh_idx == 181) {
      //   cc.print_info_detail();
      //   break;
      // } else {
      //   cc.print_info();
      // }
    }

    if (*cc.status == needs_exact_predicates) {
      cc.print_info();
    }

    if (*cc.status != success) {
      // to avoid random value assigned to thread
      convex_cells_dev[thread].status = Status::early_return;
      return;
    }

    // ninwang: no need to check security radius
    // all site neighbors are precomputed by RT
    // if (cc.is_security_radius_reached(cc.point_from_index(neigh_idx),
    //                                   is_debug)) {
    //   if (is_debug) {
    //     printf(
    //         "[clipped2] v: %d, tid: %d, seed: %d, neigh_idx: %d, "
    //         "security_radius_reached \n",
    //         v, tid, seed, neigh_idx);
    //   }
    //   break;
    // }
  }

  if (is_debug) {
    printf("seed %d has cell_vol %f \n", seed, cell_vol[seed]);
    cc.print_info();
  }

  // ninwang: no need to check security radius
  // all site neighbors are precomputed by RT
  //
  // // check security radius
  // int last_idx =
  //     site_pitch ? seed + last_row * site_knn_pitch : last_row * (seed + 1)
  //     - 1;
  // // if cell is no_intersection, then do not update status
  // if (*cc.status != no_intersection &&
  //     !cc.is_security_radius_reached(cc.point_from_index(site_knn[last_idx])))
  //     {
  //   *cc.status = security_radius_not_reached;
  // }

  if (*cc.status != no_intersection) {
    // ninwang: debug
    // if (is_debug) cc.print_info();

    // ninwang: save cc
    copy(cc, convex_cells_dev[thread], thread);

    // ninwang:
    // the sum/vol may push the center far far away
    // from the power cell
    atomic_add_bary_and_volume(cc, seed, cell_bary_sum, cell_bary_sum_pitch,
                               cell_vol);

    // use random convex cell barycenter instead
    // choose_random_bary_and_volume(cc, seed, cell_bary_sum,
    // cell_bary_sum_pitch,
    //                               cell_vol);
    // cc.get_cell_euler();
  }
}

__global__ void clipped_voro_cell_test_GPU_param(
    float* vert, int n_vert, size_t vert_pitch, int* idx, int n_tet,
    size_t idx_pitch, int* tet_knn, size_t tet_knn_pitch, int tet_k,
    Status* gpu_stat, VoronoiCell* voronoi_cells, float* cell_bary_sum,
    const size_t cell_bary_sum_pitch, float* cell_vol) {
  int thread = blockIdx.x * blockDim.x + threadIdx.x;
  int tid = thread / tet_k;  // tet index
  if (tid >= n_tet) return;

  int seed = tet_knn[(thread % tet_k) * tet_knn_pitch + tid];

  ConvexCell cc(voronoi_cells[seed], &(gpu_stat[thread]));

  clip_voro_cell_by_tet(cc, tid, vert, vert_pitch, idx, idx_pitch);

  if (*cc.status != no_intersection)
    atomic_add_bary_and_volume(cc, seed, cell_bary_sum, cell_bary_sum_pitch,
                               cell_vol);
}