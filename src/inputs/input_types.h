#pragma once
#include <assert.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_geometry.h>

#include <memory>

#include "common_geogram.h"
#include "mesh_AABB.h"
#include "params.h"

enum EdgeType {
  UE = -1,  // unkown edge
  SE = 1,   // convex sharp edge
  CE = 2    // concave edge
};

// either sharp edge or concave edge
class FeatureEdge {
 public:
  FeatureEdge(const int _id, const EdgeType &_t, const aint3 &_t2vs_group)
      : id(_id), type(_t), t2vs_group(_t2vs_group) {};
  ~FeatureEdge() {};
  void print_info() const;
  inline int get_fl_id() const { return t2vs_group[2]; }
  // NOTE: keep it same as TangentConcaveLine::operator==()
  //       used by shrink_sphere()
  inline bool operator==(FeatureEdge const &b) const {
    if (id == b.id) return true;
    // belongs to the same grouped concave line
    // this line can be a curve!!!
    if (get_fl_id() == b.get_fl_id()) {
      Vector3 direction = GEO::normalize(t2vs_pos[1] - t2vs_pos[0]);
      Vector3 b_direction = GEO::normalize(b.t2vs_pos[1] - b.t2vs_pos[0]);
      // we need to check the segment direction
      // if the direction deviate too much, then we consider them as different
      if (!is_vector_same_direction(direction, b_direction, EPS_DEGREE_30) &&
          !is_vector_oppo_direction(direction, b_direction, EPS_DEGREE_30))
        return false;
      return true;
    }
    return false;
  }

  int id;
  EdgeType type;
  // (store in tet mesh indices, matching TetMesh::tet_vertices)
  aint3 t2vs_group;  // edges in format <tvid_min, tvid_max,
                     // num_fe_group/FeatureLine::id>
  avec2 t2vs_pos;    // tet_vs indices of 2 endpoints, index matching t2vs_group
  aint2 adj_sf_fs_pair;  // mapping to GEO::Mesh sf_mesh
  avec2 adj_tan_points;  // centroid of adj_sf_fs_pair
  avec2 adj_normals;     // order maps adj_sf_fs_pair
  double length = 0.f;   // length of current feature edge
};

class ConcaveCorner {
 public:
  ConcaveCorner(const int _tvid) : tvid(_tvid) {};
  ~ConcaveCorner() {};

  void print_info() const;

  int tvid;                     // tet vid of this corner
  std::vector<int> adj_fe_ids;  // matching FeatureEdge::id
  std::vector<int> adj_fl_ids;  // matching FeatureLine::id

  // size matching adj_fe_ids.size()
  std::vector<Vector3> adj_fe_dir;  // from corner ->  another point on CE
  std::vector<double> adj_fe_len;   // maximum distance to perturb the corner
};

// sample point on FeatureLine
class FL_Sample {
 public:
  FL_Sample(const Vector3 &_p, const int &_fe) : point(_p), fe_id(_fe) {};
  ~FL_Sample() {};
  void print_info() const;

  Vector3 point;
  int fe_id;
};

// either sharp line or concave line
class FeatureLine {
 public:
  FeatureLine(const int _id, const EdgeType &_t) : id(_id), type(_t) {};
  ~FeatureLine() {};
  void print_info() const;

  int id;  // num_fe_group / t2vs_group[3]
  EdgeType type;
  bool is_loop = false;     // if feature line is a loop/cycle [no use]
  std::vector<int> fe_ids;  // sorted, matching FeatureEdge::id
  std::map<aint2, int> t2vs_to_feid;  // <tvid_min, tvid_max> -> FeatureEdge::id
  std::map<int, aint2> tvs_neighbors;  // tvid -> {tvid_left/-1, tvid_right/-1}
  double length = 0.f;                 // length of the whole feature line
  std::vector<FL_Sample> samples;      // sample points on feature line

  // feature edge fe: <tvid_min, tvid_max>
  // dir: 0 or 1, matching FeatureLine::tvs_neighbors::aint2
  // if dir = 0: return <pos_left, tvid_left>
  // if dir = 1: return <pos_right, tvid_right>
  v2int get_fe_endpoint_given_dir(const FeatureEdge &fe, int dir) const;

  // -1 if no next
  int get_next_fe_given_dir(int fe_cur, int dir);
};

class TetMesh {
 public:
  TetMesh(std::string path) : tet_path_with_ext(path) {};
  inline int get_fl_id(const int fe_id) const {
    assert(fe_id >= 0 && fe_id < feature_edges.size());
    return feature_edges[fe_id].get_fl_id();
  };

 public:
  bool is_non_cad = false;
  std::string tet_path_with_ext;
  std::vector<float> tet_vertices;  // tet vertices
  std::vector<int> tet_indices;     // tet 4 indices of vertices
  std::map<int, std::set<int>> v2tets;

  // mapping to SurfaceMesh
  std::map<int, std::set<int>> tet_vs2sf_fids;

  // adjacent info, used for RPD and topo fix
  std::vector<int> v_adjs;
  std::vector<int> e_adjs;
  std::vector<int> f_adjs;
  std::vector<int> f_ids;  // store tet's face ids, each id is unique

  // for both sharp edges and concave edges
  std::vector<FeatureEdge> feature_edges;

  // for feature edge (both sharp and concave)
  // mapping <tet_edge> -> <fe_info>
  // format: <tid, lfid_min, lfid_max> -> <fe_type, fe_id, fe_line_id>.
  // the mapping can be many-to-one (multiple tets share same sharp edge).
  //
  // fe_type:     matching EdgeType, 1 sharp edge, 2 concave edge
  // fe_id:       matching FeatureEdge::id in TetMesh::feature_edges.
  // fe_line_id:  a local id for sharp lines, matching FeatureLine::id in
  //              TetMesh::se_lines or TetMesh::ce_lines
  std::map<aint3, aint3> tet_es2fe_map;

  // for sharp edge
  // mapping vertex from <tid, lfid1, lfid2, lfid3> -> tvid
  std::map<aint4, int> tet_vs_lfs2tvs_map;
  // each feature line is grouped by feature edges
  // 1. index NOT matching FeatureLine::id
  // 2. index starts from 0
  // 3. index should have NO intersection with TetMesh::ce_lines
  std::vector<FeatureLine> se_lines;

  // for concave edges
  // each feature line is grouped by feature edges
  // 1. index NOT matching FeatureLine::id
  // 2. index starts from TetMesh::se_lines.size()
  // 3. index should have NO intersection with TetMesh::se_lines
  std::vector<FeatureLine> ce_lines;
  // allow fl to merge with other adjacent feature lines,
  // if the corner is only adjacent to concave lines, no convex lines,
  // include the id of itself
  //
  // updated by function update_allow_to_merge_ce_line_ids()
  //
  // used for clustering fids for msphere type
  // RPD3D_Wrapper::cluster_sphere_type_using_samples_fids()
  // [deprecated]
  std::map<int, std::set<int>> allow_to_merge_ce_line_ids;

  // not allow fl to merge with other adjacent feature lines,
  // controlled by TetMesh::corners_se_tet
  // updated by function update_not_allow_to_merge_ce_line_ids()
  std::map<int, std::set<int>> not_allow_to_merge_ce_line_ids;

  // for concave corners
  // adjacent to > 2 concave edges, no sharp edges
  // converted from TetMesh::corners_ce_tet
  std::vector<ConcaveCorner> cc_corners;

  // 1. regular se corner (connect to > 2 sharp edges)
  // 2. non-regular corner (> 0 sharp edges and > 0 concave edges)
  //    just to make sure we will add a zero-radius medial sphere
  std::set<int> corners_se_tet;
  // corner_tvs -> set of FeatureLine::id
  // defined by corners_se_tet
  // matching both TetMesh::se_lines and TetMesh::ce_lines
  std::map<int, std::set<int>> corner2fl;
  // corner_tvs -> set of FeatureEdge::id
  std::map<int, std::set<int>> corner2fe;
  // FeatureLine::id -> set of corner MedialSphere::isd
  // updated in init_corner_spheres()
  std::map<int, std::set<int>> fl2corner_sphere;
  // for concave corners (#adjacent_ce > 2 only)
  // convert to TetMesh::cc_corners using function store_concave_corners()
  std::set<int> corners_ce_tet;
  // not real corners
  // adjacent to > 2 (se or ce)
  // used to separate sharp edge (se) into different groups of regions
  // including all corners_se_tet and corners_ce_tet
  std::set<int> corner_fake_tet;
};

class AABBWrapper {
 private:
  GEO::Mesh se_mesh;  // for sharp edge lines
  bool is_se_mesh_exist = false;
  GEO::Mesh ce_mesh;  // for concave lines
  bool is_ce_mesh_exist = false;

  std::shared_ptr<GEO::MeshFacetsAABBLambda> sf_tree;
  std::shared_ptr<GEO::MeshFacetsAABBLambda> se_tree;  // sharp edge tree
  std::shared_ptr<GEO::MeshFacetsAABBLambda> ce_tree;  // convex edge tree

 public:
  AABBWrapper() {}

  bool init_mesh_from_edges(const std::vector<float> &input_vertices,
                            const std::vector<aint2> &edges,
                            const std::vector<int> &fe_ids, GEO::Mesh &mesh);

  // is_reorder == true => will reorder sf_mesh
  void init_sf_mesh_and_tree(GEO::Mesh &_sf_mesh, bool is_reorder = true);

  void init_feature_meshes_and_trees(
      const std::vector<float> &input_vertices,
      const std::vector<FeatureEdge> &feature_edges);

 public:
  inline int project_to_sf_get_nearest_face(Vector3 &p) const {
    Vector3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    int fidx = sf_tree->nearest_facet(p, nearest_p, sq_dist);
    p = nearest_p;
    return fidx;
  }

  inline int project_to_sf_get_nearest_face(Vector3 &p, double &sq_dist) const {
    Vector3 nearest_p;
    sq_dist = std::numeric_limits<double>::max();  //??
    int fidx = sf_tree->nearest_facet(p, nearest_p, sq_dist);
    p = nearest_p;
    return fidx;
  }

  inline int project_to_sf_get_nearest_face_with_lambdas(
      Vector3 &p, double &sq_dist, double &lambda0, double &lambda1,
      double &lambda2) const {
    Vector3 nearest_p;
    sq_dist = std::numeric_limits<double>::max();  //??
    int fidx = sf_tree->nearest_facet_lambdas(p, nearest_p, sq_dist, lambda0,
                                              lambda1, lambda2);
    // printf("lambda0 %f, lambda1 %f, lambda2 %f\n", lambda0, lambda1,
    // lambda2); Vector3 p0 = sf_tree->mesh()->vertices.point(
    //     sf_tree->mesh()->facets.vertex(fidx, 0));
    // Vector3 p1 = sf_tree->mesh()->vertices.point(
    //     sf_tree->mesh()->facets.vertex(fidx, 1));
    // Vector3 p2 = sf_tree->mesh()->vertices.point(
    //     sf_tree->mesh()->facets.vertex(fidx, 2));
    // Vector3 tmp = p0 * lambda0 + p1 * lambda1 + p2 * lambda2;
    // assert(GEO::Geom::distance(nearest_p, tmp) < 1e-6);

    p = nearest_p;
    return fidx;
  }

  inline bool get_ray_nearest_intersection(const Vector3 &orig,
                                           const Vector3 &dir, Vector3 &p,
                                           int &fid) const {
    GEO::Ray R(orig, dir);
    GEO::MeshFacetsAABBLambda::Intersection I;
    bool result = sf_tree->ray_nearest_intersection(R, I);
    p = I.p;
    fid = I.f;
    return result;
  }

  inline double project_to_sf(Vector3 &p) const {
    Vector3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    sf_tree->nearest_facet(p, nearest_p, sq_dist);
    p = nearest_p;
    return sq_dist;
  }

  inline int get_nearest_face_sf(const Vector3 &p) const {
    Vector3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    return sf_tree->nearest_facet(p, nearest_p, sq_dist);
  }

  inline int get_nearest_face_sf(const Vector3 &p, double &sq_dist) const {
    Vector3 nearest_p;
    sq_dist = std::numeric_limits<double>::max();  //??
    return sf_tree->nearest_facet(p, nearest_p, sq_dist);
  }

  inline int get_nearest_point_on_sf(const Vector3 &p, Vector3 &nearest_p,
                                     double &sq_dist) const {
    sq_dist = std::numeric_limits<double>::max();  //??
    int fidx = sf_tree->nearest_facet(p, nearest_p, sq_dist);
    return fidx;
  }

  inline double get_sq_dist_to_sf(const Vector3 &p) const {
    Vector3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    sf_tree->nearest_facet(p, nearest_p, sq_dist);
    return sq_dist;
  }

  inline double get_sq_dist_to_se(const Vector3 &p) const {
    if (!is_se_mesh_exist) return DBL_MAX;
    Vector3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    se_tree->nearest_facet(p, nearest_p, sq_dist);
    return sq_dist;
  }

  inline int project_to_se(Vector3 &p) const {
    if (!is_se_mesh_exist) return UNK_FACE;
    Vector3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    int eid = se_tree->nearest_facet(p, nearest_p, sq_dist);
    p = nearest_p;
    GEO::Attribute<int> attr_fe_ids(se_mesh.facets.attributes(), "fe_id");
    return attr_fe_ids[eid];
  }

  // inline double project_to_se(Vector3 &p) const {
  //   Vector3 nearest_p;
  //   double sq_dist = std::numeric_limits<double>::max();  //??
  //   se_tree->nearest_facet(p, nearest_p, sq_dist);
  //   p = nearest_p;
  //   return sq_dist;
  // }

  /////////////////////
  // for concave lines
  //
  // inline int project_to_ce_get_feid(Vector3 &p, double &sq_dist) const {
  //   if (!is_ce_mesh_exist) return -1;
  //   Vector3 nearest_p;
  //   sq_dist = std::numeric_limits<double>::max();  //??
  //   int eid = ce_tree->nearest_facet(p, nearest_p, sq_dist);
  //   p[0] = nearest_p[0];
  //   p[1] = nearest_p[1];
  //   p[2] = nearest_p[2];
  //   GEO::Attribute<int> attr_fe_ids(ce_mesh.facets.attributes(), "fe_id");
  //   return attr_fe_ids[eid];
  // }

  inline double get_nearest_point_on_ce(const Vector3 &p, Vector3 &nearest_p,
                                        double &sq_dist) const {
    if (!is_ce_mesh_exist) return UNK_FACE;
    sq_dist = std::numeric_limits<double>::max();  //??
    int eid = ce_tree->nearest_facet(p, nearest_p, sq_dist);
    GEO::Attribute<int> attr_fe_ids(ce_mesh.facets.attributes(), "fe_id");
    return attr_fe_ids[eid];
  }

  inline double get_sq_dist_to_ce(const Vector3 &p) const {
    if (!is_ce_mesh_exist) return DBL_MAX;
    Vector3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    ce_tree->nearest_facet(p, nearest_p, sq_dist);
    return sq_dist;
  }
};

class SurfaceMesh : public GEO::Mesh {
 public:
  void reload_sf2tet_vs_mapping();
  void collect_fid_centroids(const std::set<int> &given_fids,
                             std::vector<v2int> &one_group_fids) const;
  void cache_sf_fid_neighs();
  void cache_sf_fid_neighs_no_cross();
  void cache_sf_fid_krings_no_cross_se_only(const int k);
  void cache_sf_fid_to_fe_id(const std::vector<FeatureEdge> &feature_edges);
  void cache_fe_sf_fs_map();
  bool get_sf_fid_krings_no_cross_se_only(const int fid_given);
  bool get_sf_fid_krings_no_cross_se_only(const int fid_given,
                                          std::set<int> &kring_neighbors) const;
  bool collect_kring_neighbors_given_fid(const int k, int tan_fid,
                                         std::set<int> &kring_neighbors) const;
  bool collect_kring_neighbors_given_fid_se_only(
      const int k, int tan_fid, std::set<int> &kring_neighbors) const;
  void update_fe_sf_fs_pairs_to_ce_id(
      const std::vector<FeatureEdge> &feature_edges);

  // check if the projection point is on FE if any,
  // either CE or SE as requested in check_proj_type
  //
  // return: <SurfaceMesh::fid, FeatureEdge::id>
  // -1: no fid/FE found
  // >=0: id/FE id
  aint2 project_to_sf_and_get_FE_if_any(
      const std::vector<FeatureEdge> &feature_edges,
      const EdgeType &check_proj_type, Vector3 &p, double &sq_dist,
      bool is_debug) const;

 public:
  bool is_non_cad = false;
  AABBWrapper aabb_wrapper;
  std::vector<int> sf2tet_vs_mapping;  // matching TetMesh::tet_vertices
  // fid -> {neigh fids}
  // [no use]
  std::map<int, std::set<int>> sf_fid_neighs;
  // fid -> {neigh fids} but NOT cross fe_sf_fs_pairs,
  // updated in cache_sf_fid_neighs_no_cross(), after calling
  // detect_mark_sharp_features()
  // [no use]
  std::map<int, std::set<int>> sf_fid_neighs_no_cross;
  // fid -> {kring fids} but NOT cross fe_sf_fs_pairs_se_only,
  // updated in cache_sf_fid_krings_no_cross_se_only()
  std::map<int, std::set<int>> sf_fid_krings_no_cross_se_only;

  // for feature edges, both SE and CE
  // store <sf_fid_min, sf_fid_max> if on feature edge
  // updated in detect_mark_sharp_features()
  std::set<aint2> fe_sf_fs_pairs;
  // same as fe_sf_fs_pairs, but only for SE
  std::set<aint2> fe_sf_fs_pairs_se_only;
  // updated in cache_fe_sf_fs_map()
  // same info as fe_sf_fs_pairs_se_only()
  std::map<int, int> fe_sf_fs_map_se_only;
  // same as fe_sf_fs_pairs, but only for CE
  std::set<aint2> fe_sf_fs_pairs_ce_only;
  // updated in cache_fe_sf_fs_map()
  // same info as fe_sf_fs_pairs_ce_only
  std::map<int, int> fe_sf_fs_map_ce_only;
  // updated in update_fe_sf_fs_pairs_to_ce_id()
  std::map<aint2, int> fe_sf_fs_pairs_to_ce_id;  // mapping to FeatureEdge::id

  // surface samples
  std::vector<v2int> samples_fids;

  // fid -> FeatureEdge::id
  // updated by function SurfaceMesh::cache_sf_fid_to_fe_id()
  // used by function SurfaceMesh::project_to_sf_and_get_FE_if_any()
  std::map<int, int> sf_fid_to_fe_id;
};

void load_sf_tet_mapping(const GEO::Mesh &sf_mesh,
                         std::map<int, int> &vs_tet2sf_mapping,
                         std::map<int, std::set<int>> &vs2fids,
                         std::map<int, std::set<int>> &tet_vs2sf_fids);

void get_k_ring_neighbors_no_cross(const GEO::Mesh &sf_mesh,
                                   const std::set<aint2> &fe_sf_pairs_not_cross,
                                   const int fid_given, const int k,
                                   std::set<int> &k_ring_fids,
                                   bool is_clear_cur, bool is_debug);

void store_feature_line(const TetMesh &tet_mesh, const SurfaceMesh &sf_mesh,
                        const EdgeType &fe_type,
                        const std::vector<std::vector<aint2>> fe_tet_groups,
                        const int fl_starts_id,
                        std::vector<FeatureEdge> &feature_edges,
                        std::vector<FeatureLine> &feature_lines,
                        std::set<aint4> &fe_tet_info,
                        const std::set<int> &corners_se_tet,
                        std::map<int, std::set<int>> &corner2fl,
                        std::map<int, std::set<int>> &corner2fe, bool is_debug);

void sample_points_on_feature_line(
    const std::vector<FeatureEdge> &feature_edges, FeatureLine &fl_given,
    double esp_len, bool is_debug);
