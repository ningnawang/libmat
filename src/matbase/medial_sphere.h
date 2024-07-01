#ifndef H_medial_sphere_H
#define H_medial_sphere_H

#include <geogram/mesh/mesh.h>

#include <vector>

#include "../src/rpd3d/voronoi_defs.h"
#include "common_geogram.h"
#include "input_types.h"

// for sphere shrinking
class ss_params {
 public:
  Vector3 p, p_normal, q, q_normal;
  // positive: matching GEO::Mesh facets (starts from 0)
  // negative: matching FeatureEdge:id (starts from -2)
  int q_fid = UNK_FACE;
  int p_fid = UNK_FACE;

 public:
  inline bool is_p_on_ce() { return p_fid <= -2; }
  inline bool is_q_on_ce() { return q_fid <= -2; }
  inline void set_p_fid(int _pfid, bool is_on_ce = false) {
    p_fid = is_on_ce ? -(_pfid + 2) : _pfid;
  }
  inline void set_q_fid(int _qfid, bool is_on_ce = false) {
    q_fid = is_on_ce ? -(_qfid + 2) : _qfid;
  }
  inline int get_p_fid() { return is_p_on_ce() ? -(p_fid + 2) : p_fid; }
  inline int get_q_fid() { return is_q_on_ce() ? -(q_fid + 2) : q_fid; }
};

class TangentPlane {
 public:
  TangentPlane(const SurfaceMesh& sf_mesh, const Vector3& _normal,
               const Vector3& _point, const int _fid);
  ~TangentPlane(){};

 public:
  void print_info() const;
  bool is_same_normal(const Vector3& bnormal,
                      const double eps_degree = EPS_DEGREE_30) const;
  static bool is_same_normal(const Vector3& anormal, const Vector3& bnormal,
                             const double eps_degree = EPS_DEGREE_30);
  bool is_same_tan_pl(const TangentPlane& tan_pl2, bool is_only_normal = false,
                      double eps_p_dist = SCALAR_1) const;
  static double get_energy_value(const Vector3& theta, const double& radius,
                                 const double alpha1, const double alpha2,
                                 const Vector3& tan_point,
                                 const Vector3& normal);
  double update_energy_value(const Vector3& theta, const double& radius,
                             const double alpha1, const double alpha2);
  void update_tan_point(const Vector3& _new_tan_point);
  void update_by_sf_mesh(const GEO::Mesh& sf_mesh,
                         const AABBWrapper& aabb_wrapper);
  bool update_covered_sf_fids(const SurfaceMesh& sf_mesh, int k = K_NEIGH);

 public:
  Vector3 normal;  // normal of tangent plane
  // std::vector<Vector3> points;  // points on plane
  Vector3 tan_point;  // tangent point
  int fid;            // corresponding fid from sf_mesh, -1 as default
  double energy;      // sphere to minimize this energy function,
                      // DBL_MAX default
  double energy_over_sq_radius;  // energy / sq_radius, used as break threshold

  bool is_deleted;

  std::set<int> sf_fids_covered;  // fid's k-ring sf_mesh neighbors
};

class TangentConcaveLine {
 public:
  TangentConcaveLine(const SurfaceMesh& sf_mesh, const int _id,
                     const FeatureEdge& fe);
  ~TangentConcaveLine(){};

 public:
  void print_info() const;
  // if concave line is a curve
  // then each line should cover more adjacent faces
  bool is_normal_covered_by_adj_fs(
      const Vector3& n, double esp_degree_given = EPS_DEGREE_20) const;
  bool purge_one_tan_plane(TangentPlane& one_tan_plane) const;
  void purge_tan_planes(std::vector<TangentPlane>& tan_planes,
                        bool is_debug = false) const;

  static double get_energy_value(const Vector3& theta, const double& radius,
                                 const double alpha3, const Vector3& tan_point,
                                 const Vector3& normal);
  double update_energy_value(const Vector3& theta, const double& radius,
                             const double alpha3);

  bool update_covered_sf_fids(const SurfaceMesh& sf_mesh, int k = K_NEIGH);

 public:
  // NOTE: keep it same as FeatureEdge::operator==()
  bool operator==(TangentConcaveLine const& b) const;

 public:
  int id;
  bool is_deleted;  // [no use] for now
  int id_fe;        // matching to TetMesh::feature_edges::id
  int id_fl;        // matching to FeatureLine::id
                    // (TetMesh::feature_edges::t2vs_group[3])

  aint2 t2vs_ids;        // format <tvid_min, tvid_max>
  avec2 t2vs_pos;        // 2 tet_vs indices aint2, order maps t2vs_ids
  aint2 adj_sf_fs_pair;  // mapping to GEO::Mesh sf_mesh
  avec2 adj_normals;     // order maps adj_sf_fs_pair

  Vector3 direction;  // (t2vs_pos[1] - t2vs_pos[0]).normalize()
  Vector3 tan_point;  // tangent point, init as center of t2vs_pos
  Vector3 normal;     // a random vector inside the range of 2 adjacent normals

  double energy;  // distance from sphere to concave line, updated by alpha_3
  double energy_over_sq_radius;  // used as break threshold

  std::vector<std::set<int>>
      sf_fids_covered_two;  // adj_sf_fs_pair's k-ring sf_mesh neighbors
};

enum Topo_Status {
  low_edge_euler = -4,
  low_facet_euler = -3,
  low_cell_euler = -2,  // ninwang: site_knn __K__ might be small (fixed)
  unkown = -1,          // newly added sphere to fix
  ok = 0,
  high_cell_cc = 1,  // ninwang: tet_knn tet_k might be small (fixed)
  high_facet_cc = 2,
  high_edge_cc = 3
};

// each PowerCell contains multiple ConvexCellHost
struct PowerCell {
  PowerCell(){};
  int voro_id = -1;

  // sum of all convex cell volumes
  float pcell_vol = -1;

  //  format: [num_itr_global, cell_id]
  std::set<int> cell_ids;  // matching ConvexCellHost::id
  std::set<int> tet_ids;   // matching ConvexCellHost::tet_id

  /** For PowerCell Cells **/
  // cell_id -> {neighboring cell ids}
  std::map<int, std::set<int>> cell_neighbors;
  // // original tet fid -> 1 or 2 cell ids
  // std::map<int, std::set<int>> tfid_to_cells;
  // cell id -> orignal tet fids
  // first #sf_mesh.facets.nb()-1 matches GEO::Mesh, later are unique fids
  // from orignal tet
  std::map<int, std::set<int>> cell_to_tfids;
  // pc_face (not surface triangle) centroids on surface mesh in
  // <pc_face_centroid, surf_fid> pair (<Vector3, int>)
  // (fid defines the same as cell_to_tfids, but we only care about the sf_mesh
  // fid that matches GEO::Mesh here).
  // cell id -> vector of <pc_face_centroid, surf_fid>
  std::map<int, std::vector<v2int>> cell_to_surfv2fid;

  std::vector<std::set<int>> cc_cells;  // grouped for each CC
  // std::vector<std::set<int>> cc_bfids;  // boundary fids, matching GEO::Mesh
  // or just unique id from original tet mesh

  std::vector<std::vector<v2int>>
      cc_surf_v2fids;  // surface vertex2fids, matching GEO::Mesh only

  /** For PowerCell Facets **/
  // halfplane seed_neigh_id -> list of cell ids
  // (all cells that clipped by the halfplane)
  std::map<int, std::set<int>> facet_neigh_to_cells;
  // store seed_neigh_id that needs to add new spheres around
  // neigh_id -> fixed (true) or not fixed (false)
  std::map<int, bool> facet_neigh_is_fixed;
  // neigh_id -> { set of cell_ids in one facet CC }
  // sorted based on the size of group (in set)
  std::map<int, std::vector<std::set<int>>> facet_cc_cells;
  // neigh_id -> { set of v2fids in one facet CC }
  // relates to cell_to_surfv2fid
  // here we only care about the sf_mesh fid that matches GEO::Mesh
  std::map<int, std::vector<std::vector<v2int>>> facet_cc_surf_v2fids;

  /** For PowerCell Edges **/
  // Each edge in powercell cur_id is uniquely defined by 3 halfplanes [cur_id,
  // neigh_id_min, neigh_id_max] and we use [neigh_id_min, neigh_id_max] to
  // represent this unique edge.
  //
  // We skip one poweredge if any neigh_id is -1 (boundary face), because there
  // might be multiple poweredges as [-1, -1] but not the same poweredge
  //
  // [neigh_id_min, neigh_id_max] -> list of cell ids
  // (all cells that clipped by 2 halfplanes)
  std::map<aint2, std::set<int>> e_to_cells;
  // store shared edge that needs to add new spheres around
  // [neigh_id_min, neigh_id_max] -> fixed (true) or not fixed (false)
  std::map<aint2, bool> e_is_fixed;
  // [neigh_id_min, neigh_id_max] -> { set of cell_ids in one edge CC }
  std::map<aint2, std::vector<std::set<int>>> edge_cc_cells;
  // Each powercell edge CC is dual to a medial face of 3 spheres
  // (each powercell may have multiple edge CC)
  // Here we store all convex edges on the powercell edges.
  // Note: convex edges are not grouped by powercell edge CC!!
  //
  // if any neighbor sphere is on external features, then neigh_id_min = -1
  // we do not store this powercell edge
  //
  // used for thinning, store dual medial face's importance
  // sf_fid = -1 if not on surface
  // std::map<aint2, std::vector<v2int>> edge_endpoints;
  //
  // [neigh_id_min, neigh_id_max] ->
  // {a set of convex edges [aint2, aint2] all powercell edges}
  // aint2 is [cell_id, lvid (from cc_trans.nb_v)]
  std::map<aint2, std::vector<std::array<aint2, 2>>> edge_2endvertices;

  /** For PowerCell Vertices **/
  // Each vertex in powercell cur_id is uniquely defined by
  // [cur_id, neigh_id1, neigh_id2, neigh_id3]
  // and we use sorted [neigh_id1, neigh_id2, neigh_id3] to represent this
  // unique vertex
  //
  // 1. If neigh_id1 is -1, then this powercell vertex is on surface
  // 2. If neigh_id1 and neigh_id2 are all -1, this powercell vertex is
  // either on sharp edge (external feature), or just a convex cell vertex,
  // we do not store this vertex
  // 3. If all neigh_id(i) > 0, then vid is inside the shape. If this vertex is
  // a powercell vertex (not just a convex cell vertex), then it is dual to a
  // medial tet
  //
  // ** not including sharp edge vertices
  //
  // [cell_id, lvid (from cc_trans.nb_v)] -> [neigh_id1, neigh_id2, neigh_id3]
  std::map<aint2, aint3> vertex_2id;
  // [cell_id, lvid (from cc_trans.nb_v)] -> <pos, sf_fid>
  // sf_fid can be -1 if vertex is not on surface
  // ** not including sharp edge vertices
  std::map<aint2, v2int> vertex_2pos;

  /** For External Edge Features **/
  // all touched sharp edges
  // stored in [cell_id, lvid1, lvid2, fe_line_id, fe_id]
  // fe_id matching FeatureEdge::id in TetMesh::feature_edges
  // se_line_id matching FeatureLine::id in TetMesh::se_lines
  std::set<aint5> se_covered_lvids;
  // each powercell's sharp line endpoint has a unique neigh_id
  // (maximum 2 endpoints for each powercell's sharp line)
  // [cell_id, lvid, neigh_id, se_line_id] -> pos
  // se_line_id matching FeatureLine::id in TetMesh::se_lines
  // used for checking and fixing geometry on sharp edges
  std::map<aint4, Vector3> se_line_endpos;

  /** For Concave Edges  **/
  // all touched concave edges
  // stored in [cell_id, lvid1, lvid2, fe_line_id, fe_id]
  // fe_id matching FeatureEdge::id in TetMesh::feature_edges
  // ce_line_id matching FeatureLine::id in TetMesh::ce_lines
  std::set<aint5> ce_covered_lvids;

  /** For Internal Features **/
  // grouped sf_mesh <pc_face_centroid, surf_fid> pairs (<Vector3, int>)
  // not crossing sharp edges / concave edges
  // only for surface fid this pcell is touched, no tangent info
  std::vector<std::vector<v2int>> surf_v2fid_in_groups;

  Topo_Status topo_status;
};

enum SphereType {
  T_1_N = -3,  // external feature, corners
  T_1_2 = -2,  // external feature, sharp edge
  T_UNK = -1,
  T_2 = 2,
  T_N = 3,       // including T_3, T_4 ...
  T_N_JUNC = 4,  // potential junctions T_4, T_5 ...
  T_c = 10,      // invisible pin sphere added on concave line to avoid RPD
                 // degeneration [no use?]
  T_2_c = 11,    // added for one concave line using sphere shrinking
  T_N_c = 12,    // added through internal feature preservation
  T_3_c = 13,    // added for concave corners (debug only) [no use]
  // todo: deprecate this
  T_X_c = 14  // added for one concave edge during topo fix
};

enum RtValidStates { Valid2Invalid = -1, NoChange = 0, Invalid2Valid = 1 };

class MedialSphere {
 public:
  MedialSphere(){};  // for mult-threads parallization
  MedialSphere(int _id, Vector3 _pin, Vector3 _pin_normal, int _pin_fid,
               SphereType _type = SphereType::T_UNK, int _itr = -1);
  MedialSphere(int _id, Vector3 _center, double _radius, SphereType _type,
               int _itr = -1);
  ~MedialSphere(){};

 public:
  void print_info() const;
  void print_ss_info() const;
  void print_tan_planes() const;
  void print_tan_cc_lines() const;
  int get_tan_element_size() const;
  int get_sf_covered_group_size() const;

  double update_all_energy_values(const double alpha1, const double alpha2,
                                  const double alpha3, bool is_debug);
  void save_old_center_radius(bool is_clear = true);
  void update_tan_planes_from_ss_params(const SurfaceMesh& sf_mesh,
                                        bool is_update_p, bool is_update_q);
  void new_tan_plane_no_dup(const SurfaceMesh& sf_mesh, const Vector3& _normal,
                            const Vector3& _point, const int _fid,
                            bool is_only_normal);
  void update_tan_cc_lines_from_ss_params(const SurfaceMesh& sf_mesh,
                                          const std::vector<FeatureEdge>& fe,
                                          bool is_update_p, bool is_update_q);
  bool new_cc_line_no_dup(const TangentConcaveLine& one_cc_line);
  void get_sphere_all_tangent_pairs_includes_cc_lines(
      std::vector<std::array<Vector3, 2>>& tan_pairs) const;
  void get_sphere_all_tangent_pairs_includes_cc_lines(
      std::vector<avec2int>& tan_pairs) const;
  void update_tan_planes_by_sf_mesh(const GEO::Mesh& sf_mesh,
                                    const AABBWrapper& aabb_wrapper);

  void purge_and_delete_tan_planes();
  void remove_deleted_tangents(bool is_run_cc);
  void dilate_sphere_radius();
  bool is_same_tangent_info(const MedialSphere& msphereB);
  bool is_tangent_info_covered_by_B(const MedialSphere& msphereB);
  // updating MedialSphere::covered_sf_fids_in_group
  void update_sphere_covered_sf_fids(const SurfaceMesh& sf_mesh, bool is_debug);
  // [no use] so far, use update_sphere_covered_sf_fids()
  void update_tangent_covered_fids_by_sf_mesh(const SurfaceMesh& sf_mesh,
                                              int k = K_NEIGH);

  bool is_on_se() const;
  bool is_on_corner() const;
  bool is_on_ce_pin() const;
  bool is_on_ce_pre() const;
  bool is_on_ce_pre_or_fix() const;
  bool is_on_ce() const;
  bool is_on_intf() const;
  bool is_on_extf() const;
  bool is_on_sheet() const;
  bool is_on_junction() const;

  void topo_clear();
  void pcell_insert(int cell_id);
  // We only consider Topo_Status::high_facet_cc as true
  // iff both spheres have Topo_Status::high_facet_cc
  bool fcc_is_to_fix(int neigh_id) const;
  void fcc_fixed(int neigh_id);

  bool operator<(const MedialSphere& m2) const;
  bool operator==(const MedialSphere& m2) const;
  bool operator!=(const MedialSphere& m2) const;
  bool is_sphere_too_close(const MedialSphere& m2, double threshold) const;
  void copy(const MedialSphere& b);

  // convert ss_params fid <-> eid
  inline static int convert_ss(const int _id) { return -(_id + 2); };

 public:
  int id;
  int site_id;  // for computing parital RPD, -1 as default
  int itr_cnt;  // count the number of iteration when this sphere is created
  bool is_deleted = false;
  int dup_cnt = 0;  // duplicated count, one matching a connected pcell

  Vector3 center;
  double radius;
  bool is_radius_dilated = false;
  SphereType type;
  ss_params ss;
  std::vector<TangentPlane> tan_planes;
  std::vector<TangentConcaveLine> tan_cc_lines;

  /* For RT related */
  /* mostly used for computing partial RPDs */
  bool is_rt_valid;  // valid in RT, reset for each RT, false as default
  bool is_rt_prev_valid;
  // -1 -> valid -> invalid
  // 0  -> no change
  // 1  -> invalid -> valid
  int rt_change_status;             // matching RtValidStates
  std::set<int> rt_neigh_ids_prev;  // matching MedialSphere::id, used for
                                    // computing partial RPDs
  std::set<int> rt_neigh_ids;       // matching MedialSphere::id,
                                    // used for computing partial RPDs

  /* For topology check */
  PowerCell pcell;  // one powercell per sphere
  // Euler characteristics for all cells
  double euler = 0;  // = euler_sum - num_cells
  double euler_sum = 0;
  uint num_cells = 0;
  int itr_topo_fix = 0;  // count the topo fix for this sphere

  /* For medial mesh */
  std::set<int> edges_;  // matches MedialMesh::edges
  std::set<int> faces_;  // matches MedialMesh::faces
  void clear_mm();

  // for relaxation (CVT-ish) and sphere iterating
  Vector3 old_center;
  double old_radius;

  /* For external feature */
  // matching FeatureEdge::id
  int se_edge_id = -1;
  // sharp line id, matching FeatureLine::id (TetMesh::se_lines)
  int se_line_id = -1;
  // corners feature line id, matching FeatureLine::id
  std::set<int> corner_fls;
  // matching FeatureEdge::id
  std::set<int> corner_fes;

  /* For concave lines */
  double min_sq_dist_2cc = -1;

  /* For internal feature */
  // collect surface regions the sphere is covered, in groups
  // use both 1. pcell, 2. tangent info [NO]
  // use only **tangent info** [YES]
  // format <pc_sf_centroid/sf_centroid, surf_fid> pairs (<Vector3, int>)
  //
  // updated by update_sphere_covered_sf_fids()
  std::vector<std::vector<v2int>> covered_sf_fids_in_group;
  void print_covered_sf_fids_in_group() const;

  // [no use]
  // updated by RPD3D_Wrapper::update_pcell_samples_sf_fids()
  std::vector<std::set<int>> pcell_samples_sf_fids_in_group;
};

void update_power_cells(const SurfaceMesh& sf_mesh,
                        std::vector<ConvexCellHost>& convex_cells_host,
                        std::vector<MedialSphere>& all_medial_spheres,
                        const std::map<aint3, aint3>& tet_es2fe_map,
                        bool is_debug);

bool validate_new_sphere(const std::vector<MedialSphere>& all_medial_spheres,
                         const MedialSphere& new_sphere,
                         bool is_small_threshold = false,
                         bool is_debug = false);
bool add_new_sphere_validate(std::vector<MedialSphere>& all_medial_spheres,
                             MedialSphere& new_sphere,
                             bool is_small_threshold = false,
                             bool is_debug = false);

int remove_duplicated_medial_spheres(
    std::vector<MedialSphere>& all_medial_spheres);

// no use
// replaced by function generate_RT_CGAL_and_purge_spheres()
void purge_deleted_medial_spheres(
    std::vector<MedialSphere>& all_medial_spheres);

bool is_two_mspheres_on_same_se(const MedialSphere& msphere1,
                                const MedialSphere& msphere2);

// must call after update_power_cells()
void update_se_tangent_planes(const SurfaceMesh& sf_mesh,
                              const std::vector<FeatureEdge>& feature_edges,
                              std::vector<MedialSphere>& all_medial_spheres,
                              bool is_clear_existing);

void copy_powercell_volume(const std::vector<float>& all_cell_vol,
                           std::vector<MedialSphere>& all_medial_spheres);

// load spheres that RT status has changed
// status of MedialSphere::rt_change_status (matches RtValidStates)
// deleted spheres are loaded in invalid_spheres, not in changed_spheres
void load_changed_spheres(const int num_itr_global,
                          const std::vector<MedialSphere> all_medial_spheres,
                          std::set<int>& changed_spheres,
                          std::set<int>& invalid_spheres, bool is_debug);

// convert tangent elements to covered surface fids
void A_B_spheres_common_diff_tangent_info_from_surface(
    const SurfaceMesh& sf_mesh, const MedialSphere& mat_A,
    const MedialSphere& mat_B, std::vector<TangentPlane>& common_tan_pls,
    std::vector<TangentPlane>& A_non_common_tan_pls,
    std::vector<TangentPlane>& B_non_common_tan_pls, bool is_debug = false);

// helper functions
v2int get_v2fid_max_to_point(const GEO::Mesh& sf_mesh,
                             const std::vector<v2int> surf_v2fids,
                             const Vector3& point);

v2int get_v2fid_min_to_point(const GEO::Mesh& sf_mesh,
                             const std::vector<v2int> surf_v2fids,
                             const Vector3& point);

// for degenerated medial spheres
int delete_degenerated_medial_spheres(
    std::vector<MedialSphere>& all_medial_spheres, bool is_debug);

#endif  // __H_medial_sphere_H__