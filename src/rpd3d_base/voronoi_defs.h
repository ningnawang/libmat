#pragma once

#include <assert.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "common_cxx.h"  // no cuda, no geogram
#include "voronoi_common.h"

//----------------------------------Host
struct ConvexCellHost {
  ConvexCellHost(){};
  cfloat4 compute_vertex_coordinates(cuchar3 t, bool persp_divide = true) const;
  void print_info() const;

  // return references
  inline cuchar4& ver_trans(int t) { return ver_data_trans[t]; }
  inline cfloat5& clip_trans(int p) { return clip_data_trans[p]; }
  inline cuchar3& edge(int e) { return edge_data[e]; }
  inline uchar& ith_plane(uchar t, int i) {
    return reinterpret_cast<uchar*>(&(ver_trans(t)))[i];
  }

  // return constants
  inline cuchar4 ver_trans_const(int t) const { return ver_data_trans[t]; }
  inline cfloat4 clip_trans_const(int f) const {
    return cmake_float4(clip_data_trans[f].x, clip_data_trans[f].y,
                        clip_data_trans[f].z, clip_data_trans[f].w);
  }
  inline cint2 clip_id2_const(int p) const { return clip_id2_data_trans[p]; }
  inline cuchar2 edge_id2_const(int e) const {
    return cmake_uchar2(edge_data[e].x, edge_data[e].y);
  }
  inline uchar ith_plane_const(uchar t, int i) const {
    cuchar4 v = ver_trans_const(t);
    return reinterpret_cast<uchar*>(&(v))[i];
  }

  // topology
  bool is_active_updated = false;
  std::vector<int> active_clipping_planes;
  std::vector<int> active_edges;
  Scalar cal_cell_euler();
  Scalar cal_halfplane_facet_euler(int neigh_id);
  Scalar cal_edge_euler(aint2 neigh_min_max);
  void reload_active();

  // explicity compute a powercell's vertices coordinates and faces
  // use vertex as base element (instead of face)
  // [used for topology check and GUI]
  bool is_pc_explicit = false;
  std::vector<cfloat3> pc_points;
  std::vector<std::vector<int>> pc_local_active_faces;
  // map of lf -> active_lf; size = nb_p; -1 if not active
  std::vector<int> pc_lf2active_map;
  void reload_pc_explicit();
  // some vertices may be null
  bool is_vertex_null = false;

  // debug
  void print_cell_detail_euler() const;

  Status status;
  int thread_id;
  int voro_id;  // may not matching MedialSphere::id in all_medial_spheres,
                // if given partial spheres
  int tet_id;
  float weight = -1;  // sq_radius
  bool is_active = false;
  // [Option1]: let CPU compute explicit
  // (so that we can remove duplicated voro points)
  uchar nb_v, nb_p, nb_e;
  // format: [clip1, clip2, clip3, #adjacent_cells]
  cuchar4 ver_data_trans[_MAX_T_];
  // format: [a, b, c, d, #adjacent_cells]
  //         [x, y, z, w, h              ]
  cfloat5 clip_data_trans[_MAX_P_];
  // format: [fid1, fid2]
  // we store a unique id for clipping plane
  // fid matches sf_mesh::fid for faces on sf_mesh
  // 1. plane inherit from tets, we store [fid, -1]
  // 2. halfplane of two seeds, we store [seed_id_min, seed_id_max]
  cint2 clip_id2_data_trans[_MAX_P_];
  // format: [clip1, clip2, #adjacent_cells]
  cuchar3 edge_data[_MAX_E_];
  float euler = -1;
  float cell_vol = -1;

  // for check CC and Euler, debug and GUI
  // will be saved at PowerCell::cell_ids
  int id = UNK_INT;  // matching convex_cells_host_non_dup
};

bool is_convex_cell_valid(const ConvexCellHost& cell);
