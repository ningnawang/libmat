#ifndef H_CONVEX_CELL_H
#define H_CONVEX_CELL_H

#include <assert.h>
#include <cuda_runtime.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "common_cuda.h"
#include "voronoi_common.h"

// [no use]
struct VoronoiCell {
  Status status;
  uchar nb_v, nb_p;
  uchar3 ver[_MAX_T_];
  float4 clip[_MAX_P_];
  float4 voro_seed;
};

// 4 faces of tet abcd: cbd acd bad abc (in vids)
__constant__ int tet_faces_lvid[4][3] = {
    {2, 1, 3}, {0, 2, 3}, {1, 0, 3}, {0, 1, 2}};
// 6 edges of tet (in vids)
__constant__ int tet_edges_lvid[6][2] = {{2, 3}, {1, 3}, {1, 2},
                                         {0, 3}, {0, 2}, {0, 1}};

// For orignal tet mesh
//
// Edge Eulers are stored as a diagonal adjacency matrix with size n
// given two sorted vertex indices (vmin, vmax), we can get
// the index for this edge
//
// Edge e -> (vid_min, vid_max)
// n -> #vertices
// idx(e) = idx(vid_min, vid_max)
//             = n + (n-1) + ... + (n-vid_min) - (n - vid_max)
inline __host__ __device__ int get_edge_idx(int v1, int v2, int n) {
  int vmin = v1;
  int vmax = v2;
  if (v1 > v2) {
    vmin = v2;
    vmax = v1;
  }
  int idx = 0;
  for (uint j = 0; j <= vmin; j++) {
    idx += (n - j);
  }
  idx -= (n - vmax);
  return idx;
}
inline __host__ __device__ int get_e_adj(const int* e_adjs, const int v_size,
                                         const int v1, const int v2) {
  int eid = get_edge_idx(v1, v2, v_size);
  int edge_max = v_size * (1 + v_size) / 2 + 1;
  assert(eid >= 0 && eid < edge_max);
  return e_adjs[eid];
}

// NOTE: #adjacent_cells are int!!
//
// memory pool for chained lists of triangles
// cell's vertex -> dual triangle
// [clip1, clip2, clip3, #adjacent_cells]
__shared__ uchar4 ver_data[VORO_BLOCK_SIZE * _MAX_T_];
// cell's clipping plane -> dual vertex
// clipping planes, some might not exist anymore
// (not in any dual triangle), but we still store it
// check 'active_clipping_planes' in get_convex_cell()
//
// we also store a unique id for clipping plane
// 1. plane inherit from tets, we store [fid, -1]
// 2. halfplane of two seeds, we store [seed_id_min, seed_id_max]
//
// format: [a, b, c, d, #adjacent_cells, fid1, fid2]
//         [x, y, z, w, h,               g,    k   ]
__shared__ float7 clip_data[VORO_BLOCK_SIZE * _MAX_P_];
// cell edges in [clip1, clip2, #adjacent_cells, feature_id]
// feature: UNK_INT(-1) / EdgeType::SE or EdgeType::CE
// some may not exist anymore
__shared__ uchar3 edge_data[VORO_BLOCK_SIZE * _MAX_E_];
// dual vertex IDs point directly to the right locaction
// in the circular lsit
__shared__ uchar boundary_next_data[VORO_BLOCK_SIZE * _MAX_P_];

inline __device__ uchar4& ver(int v) {
  assert(v < _MAX_T_);
  return ver_data[threadIdx.x * _MAX_T_ + v];
}

inline __device__ uchar& boundary_next(int v) {
  return boundary_next_data[threadIdx.x * _MAX_P_ + v];
}
inline __device__ float7& clip(int p) {
  return clip_data[threadIdx.x * _MAX_P_ + p];
}
inline __device__ uchar3& edge(int e) {
  assert(e >= 0);
  return edge_data[threadIdx.x * _MAX_E_ + e];
}
// return a copy, only used in const functions!
inline __device__ uchar3 ver3(int v) { return make_uchar3(ver(v)); }
inline __device__ float4 clip4(int p) { return make_float4(clip(p)); }
inline __device__ float5 clip5(int p) { return make_float5(clip(p)); }
inline __device__ int2 clip_id2(int p) {
  return make_int2(__float2int_rn(clip(p).g), __float2int_rn(clip(p).k));
}

struct ConvexCell {
  // no use
  __device__ ConvexCell(const int p_seed, const float* p_pts,
                        const size_t pitch, Status* p_status);
  __device__ ConvexCell(const int p_seed, const float* p_pts,
                        const size_t pitch, Status* p_status, const int tid,
                        const float* vert, const size_t vert_pitch,
                        const int* idx, const size_t idx_pitch);
  __device__ ConvexCell(const VoronoiCell& vc, Status* p_status);
  // in use
  __device__ ConvexCell(const int p_seed, const float* p_pts,
                        const size_t pitch, const float* p_weights,
                        const uint* p_flags, Status* p_status, const int tid,
                        const float* vert, const int n_vert,
                        const size_t vert_pitch, const int* idx,
                        const size_t idx_pitch, const int* v_adjs,
                        const int* e_adjs, const int* f_adjs, const int* f_ids);
  __device__ void clip_by_plane(int neigh_seed_id);
  __device__ void clip_by_plane(float4 eqn);
  __device__ bool is_vertex_perturb(uchar3 v);
  __device__ float4 compute_vertex_coordinates(uchar3 t,
                                               bool persp_divide = true) const;
  __device__ inline uchar& ith_plane(uchar t, int i);
  // halfplane, defined by neigh_seed_id and voro_id
  __device__ int new_plane(int neigh_seed_id);
  __device__ int new_plane(float4 eqn);
  __device__ int find_edge_id(uchar clip1, uchar clip2);
  __device__ void new_vertex(uchar i, uchar j, uchar k);  // dual triangle
  __device__ void new_edge(uchar clip1, uchar clip2);
  __device__ void compute_boundary();
  __device__ bool is_security_radius_reached_legacy(float4 last_neig,
                                                    bool is_debug = false);
  __device__ bool is_security_radius_reached(float4 last_neig,
                                             bool is_debug = false);
  __device__ float4 point_from_index(int idx);
  __device__ void print_info() const;
  __device__ void print_info_detail() const;
  __device__ void print_vs() const;
  __device__ float get_cell_euler();  // TODO: no use

  __device__ bool cc_vertex_is_in_conflict(uchar3 v, float4 eqn,
                                           bool is_debug) const {
    // return cc_vertex_is_in_conflict_float(v, eqn, is_debug);
    return cc_vertex_is_in_conflict_double(v, eqn, is_debug);
  }
  __device__ bool cc_vertex_is_in_conflict_float(uchar3 v, float4 eqn,
                                                 bool is_debug) const;
  __device__ bool cc_vertex_is_in_conflict_double(uchar3 v, float4 eqn,
                                                  bool is_debug) const;

  Status* status;
  uchar nb_v;  // number of vertices
  uchar nb_r;
  const float* pts = nullptr;  // for storing all sites
  const size_t pts_pitch = 0;
  const float* pts_weights = nullptr;  // for storing site weights (sq_radius)
  int voro_id;       // may not matching MedialSphere::id in all_medial_spheres,
                     // if given partial spheres
  int tet_id;        // index of tet, to filter multiple seed&tet pairs
  float4 voro_seed;  // (x,y,z,sq_radius)
  uint voro_flag;
  uchar nb_p;  // number of clipping planes
  uchar first_boundary_;
  uchar nb_e;
  float euler = -1;

  // // for debug, asigned by host
  // int tid;
};

//----------------------------------Host
//-------------------------TODO: use ConvexCell instead
struct ConvexCellTransfer {
  __host__ ConvexCellTransfer(){};
  Status status;
  int thread_id;
  int voro_id;  //  matching all_medial_spheres
  int tet_id;
  float weight = -1;  // sq_radius
  bool is_active = false;
  // [Option1]: let CPU compute explicit
  // (so that we can remove duplicated voro points)
  uchar nb_v, nb_p, nb_e;
  // format: [clip1, clip2, clip3, #adjacent_cells]
  uchar4 ver_data_trans[_MAX_T_];
  // format: [a, b, c, d, #adjacent_cells]
  //         [x, y, z, w, h              ]
  float5 clip_data_trans[_MAX_P_];
  // format: [fid1, fid2]
  // we store a unique id for clipping plane
  // 1. plane inherit from tets, we store [fid, -1]
  // 2. halfplane of two seeds, we store [seed_id_min, seed_id_max]
  int2 clip_id2_data_trans[_MAX_P_];
  // format: [clip1, clip2, #adjacent_cells]
  uchar3 edge_data[_MAX_E_];
  float euler = -1;
  float cell_vol = -1;

  // for check CC and Euler, debug and GUI
  int id = UNK_INT;  // matching convex_cells_host_non_dup
};
__host__ bool is_convex_cell_valid(const ConvexCellTransfer& cell);

//---------------------------------- Functions
__global__ void clipped_voro_cell_test_GPU_param_tet(
    const float* site, const int n_site,
    const size_t site_pitch /*always true*/, const float* site_weights,
    const uint* site_flags, const int* site_knn, const size_t site_knn_pitch,
    const int site_k, const float* vert, const int n_vert,
    const size_t vert_pitch, const int* idx, const int n_tet,
    const size_t idx_pitch, const int* v_adjs, const int* e_adjs,
    const int* f_adjs, const int* f_ids, const int* tet_knn,
    const size_t tet_knn_pitch, const int tet_k, Status* gpu_stat,
    VoronoiCell* voronoi_cells, ConvexCellTransfer* convex_cells_dev,
    float* cell_bary_sum, const size_t cell_bary_sum_pitch, float* cell_vol);

__global__ void clipped_voro_cell_test_GPU_param(
    float* vert, int n_vert, size_t vert_pitch, int* idx, int n_tet,
    size_t idx_pitch, int* tet_knn, size_t tet_knn_pitch, int tet_k,
    Status* gpu_stat, VoronoiCell* voronoi_cells, float* cell_bary_sum,
    const size_t cell_bary_sum_pitch, float* cell_vol);

#endif  // __CONVEX_CELL_H__
