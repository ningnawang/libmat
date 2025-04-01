#pragma once

// user-defined
struct Parameter {
  // original mesh
  float xmin_orig = -1., ymin_orig = -1., zmin_orig = -1., xmax_orig = -1.,
        ymax_orig = -1., zmax_orig = -1.;
  float bbox_diag_l_prev = 0.f;
  float scale_maxside_orig = 0.f;  // for normalization
  // normalized mesh
  float xmin = -1., ymin = -1., zmin = -1., xmax = -1., ymax = -1., zmax = -1.;
  float bbox_diag_l = 0.f;
  std::vector<float> bb_points;  // 8 bbox

  float scale_max = 1000.f;  // normalize models between [0,1000]^3

  // For fix_geo sampling step
  double geo_rel = 1. / 20.;
  double geo_error_bb_per = 0.4;

  // For finding sphere neighbors (using RT)
  // will be updated by function get_RT_vertex_neighbors()
  int site_k = 90;

  // For surface poisson disk sampling = sf_face_len * bbox_diag_l
  float sf_face_len = 1. / 20.;

  // For remove_close_medial_spheres()
  float clean_sphere_thres_rel = 1. / 300.f;

  // For ideal mface len = mface_rel_len * bbox_diag_l
  float mface_rel_len = 1. / 30.;

  // For adding spheres on feature edges if too long
  float fe_rel_len = 1. / 50.;
  // 10 for pipeline fix_extf
  // 30 for pipeline fix_geo

  // For measuring the Hausdorff distance
  // between medial slab and the surface, over bbox_diag_l
  float hd_rel_slab = 1. / 40.;
  // RPC volume over medial sphere volume, smaller the better
  float hd_rel_sphere = 2.5;

  // For adding medial spheres that tangent to concave lines using
  // ball-shrinking algo. We need to specify two parameters:
  // 1. cc_len_eps_rel:
  //    defines the length between to pin points on concave lines, relative to
  //    bbox_diag_l
  // 2. cc_normal_eps:
  //    define the angle between two normals given a pin point, scaled in
  //    [0,360]
  // Given a pin point and a normal, we can add a sphere using ball-shrinking
  // algo. Concave lines can be better preserved when cc_len_eps is smaller
  // (more medial spheres inserted), but it also means more processing time.
  double cc_len_eps_rel = 1. / 70.;
  double cc_normal_eps = 20;  // [no use]

  // Threshold for detecting sharp/concave edges and corners [0,180]
  // convex, bigger angle more convex, so smaller threshold more sensative
  double thres_convex = 55.;
  // concave, bigger angle more concave, so smaller threshold more sensative
  double thres_concave = 60.;
};

//----------------------------------------------------------------------------
// for sphere shriking, matching Parameter::scale_max
constexpr float INIT_RADIUS = 1000.f;
// for k-ring neighbors on sf_mesh
// each sphere's covered sf_mesh fids
constexpr int K_NEIGH = 3;

// // Feature Edges
// constexpr int SHARP_EDGE = 1;
// constexpr int CONCAVE_EDGE = 2;

// Topo fix
constexpr int TOPO_ITR_MAX = 31;
constexpr int ITR_MAX = 31;
// Geo fix
constexpr int GEO_ITR_MAX = 5;
constexpr unsigned int GEO_SAMPLE_MAX = 10000;
// Intf/Extf fix
// for non_cad models
constexpr int INTF_ITR_MAX = 5;
constexpr int EXTF_ITR_MAX = 5;

#define RAN_SEED 200                // random seed
#define SCALAR_FEATURE_RADIUS 1e-1  // for SE sphere
#define SCALAR_SE_MERGE_RADIUS 1    /* scaled to [0,1000]^3 */
#define SCALAR_CE_PIN_RADIUS 30
#define SQ_DIST_TO_CC 1e-2
