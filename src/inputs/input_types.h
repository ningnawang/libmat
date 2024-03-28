#pragma once
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_geometry.h>

#include <memory>

#include "common_geogram.h"
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
      : id(_id), type(_t), t2vs_group(_t2vs_group){};
  ~FeatureEdge(){};
  void print_info() const;

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
  ConcaveCorner(const int _tvid) : tvid(_tvid){};
  ~ConcaveCorner(){};

  void print_info() const;

  int tvid;                     // tet vid of this corner
  std::vector<int> adj_fe_ids;  // matching FeatureEdge::id

  // size matching adj_fe_ids.size()
  std::vector<Vector3> adj_fe_dir;  // from corner ->  another point on CE
  std::vector<double> adj_fe_len;   // maximum distance to perturb the corner
};

// sample point on FeatureLine
class FL_Sample {
 public:
  FL_Sample(const Vector3 &_p, const int &_fe) : point(_p), fe_id(_fe){};
  ~FL_Sample(){};
  void print_info() const;

  Vector3 point;
  int fe_id;
};

// either sharp line or concave line
class FeatureLine {
 public:
  FeatureLine(const int _id, const EdgeType &_t) : id(_id), type(_t){};
  ~FeatureLine(){};
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
  TetMesh(std::string path) : tet_path_with_ext(path){};

 public:
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
  std::vector<FeatureLine> se_lines;

  // for concave edges
  // each feature line is grouped by feature edges
  std::vector<FeatureLine> ce_lines;

  // for concave corners
  std::vector<ConcaveCorner> cc_corners;

  // 1. regular se corner (connect to > 2 sharp edges)
  // 2. non-regular corner (> 0 sharp edges and > 0 concave edges)
  //    just to make sure we will add a zero-radius medial sphere
  std::set<int> corners_se_tet;
  // corner_tvs -> set of FeatureLine::id
  std::map<int, std::set<int>> corner2fl;
  // corner_tvs -> set of FeatureEdge::id
  std::map<int, std::set<int>> corner2fe;
  // FeatureLine::id -> set of corner MedialSphere::id
  // updated in init_corner_spheres()
  std::map<int, std::set<int>> fl2corner_sphere;
  // for concave corners (#adjacent_ce > 2 only)
  std::set<int> corners_ce_tet;
  // not real corners
  // adjacent to > 2 se and ce
  // used to separate sharp edge (se) into different groups of regions
  // including all corners_se_tet
  std::set<int> corner_fake_tet;
};

class AABBWrapper {
 private:
  GEO::Mesh se_mesh;  // for sharp edge lines
  bool is_se_mesh_exist = false;
  GEO::Mesh ce_mesh;  // for concave lines
  bool is_ce_mesh_exist = false;

  std::shared_ptr<GEO::MeshFacetsAABB> sf_tree;
  std::shared_ptr<GEO::MeshFacetsAABB> se_tree;  // sharp edge tree
  std::shared_ptr<GEO::MeshFacetsAABB> ce_tree;  // convex edge tree

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

  inline bool get_ray_nearest_intersection(const Vector3 &orig,
                                           const Vector3 &dir, Vector3 &p,
                                           int &fid) const {
    GEO::Ray R(orig, dir);
    GEO::MeshFacetsAABB::Intersection I;
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

  inline double project_to_ce(Vector3 &p) const {
    Vector3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    se_tree->nearest_facet(p, nearest_p, sq_dist);
    p = nearest_p;
    return sq_dist;
  }

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
  void reload_sf_fid_neighs();
  void reload_sf_fid_neighs_no_cross();  // must call after
  // detect_mark_sharp_features()
  void collect_fid_centroids(const std::set<int> &given_fids,
                             std::vector<v2int> &one_group_fids) const;
  bool collect_kring_neighbors_given_fid(const int k, int tan_fid,
                                         std::set<int> &kring_neighbors) const;

 public:
  AABBWrapper aabb_wrapper;
  std::vector<int> sf2tet_vs_mapping;          // matching TetMesh::tet_vertices
  std::map<int, std::set<int>> sf_fid_neighs;  // fid -> {neigh fids} [no use]

  // for feature edges, both SE and CE
  // store <sf_fid_min, sf_fid_max> if on feature edge
  // updated in detect_mark_sharp_features()
  std::set<aint2> fe_sf_fs_pairs;
  // fid -> {neigh fids} but NOT cross fe_sf_fs_pairs,
  // updated in reload_sf_fid_neighs_no_cross(), after calling
  // detect_mark_sharp_features()
  std::map<int, std::set<int>> sf_fid_neighs_no_cross;
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
                        std::vector<FeatureEdge> &feature_edges,
                        std::vector<FeatureLine> &feature_lines,
                        std::set<aint4> &fe_tet_info,
                        const std::set<int> &corners_se_tet,
                        std::map<int, std::set<int>> &corner2fl,
                        std::map<int, std::set<int>> &corner2fe, bool is_debug);

void sample_points_on_feature_line(
    const std::vector<FeatureEdge> &feature_edges, FeatureLine &fl_given,
    double esp_len, bool is_debug);
