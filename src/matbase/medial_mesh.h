#pragma once

#include "common_geogram.h"
#include "medial_primitives.h"
#include "medial_sphere.h"

enum MedialType {
  MUNKOWN = -1,
  SHEET = 0,
  SEAM = 1,     /*matching internal feature*/
  BOUNDARY = 2, /*matching external feature*/
  JUNCTION = 3
};

class MedialStruct {
 public:
  MedialStruct(int _id, MedialType _t) : id(_id), type(_t) {};
  ~MedialStruct() {};

 public:
  int id;
  MedialType type;
  std::set<int> m_face_ids;    // MedialFace::id
  std::set<int> m_edge_ids;    // MedialEdge::id
  std::set<int> m_sphere_ids;  // MedialSphere::id
};

class MedialEdge {
 public:
  int eid;
  int dup_cnt;  // duplicated count
  aint2 vertices_;
  std::set<int> faces_;  // triangle list
  bool is_deleted = false;
  bool is_intf = false;
  bool is_extf = false;
  bool is_on_same_sheet = false;  // including is_intf

  // for medial structure
  int mstruct_id = -1;  // MedialStruct::id
  std::vector<TangentPlane> common_tan_pls;
  // order may NOT matches MedialEdge::vertices_
  // we want to make sure v1 is on internal feature as possible
  bool is_vs_swapped = false;  // swap the order of vertices_
  std::vector<TangentPlane> v1_non_common_tan_pls;
  std::vector<TangentPlane> v2_non_common_tan_pls;

  bool HasVertex(int vid) {
    return ((vertices_[0] == vid) || (vertices_[1] == vid));
  }
  bool HasFace(int fid) { return (faces_.find(fid) != faces_.end()); }
};

class MedialFace {
 public:
  // setup in function get_triangles_from_three_spheres()
  enum Type {
    INVALID = 1,      // -- fail to generate 2 simple triangles
    TETRAHEDRON = 2,  // -- share an edge
    PRISM = 3         // -- share a point or nothing
  };
  int fid;
  int dup_cnt;  // duplicated count
  Type type;
  aint3 vertices_;
  aint3 edges_;
  std::set<int> tets_;  // neighboring tetrahedrons
  bool is_deleted = false;

  // for type
  bool is_on_same_sheet = false;

  // for medial structure
  int mstruct_id = -1;  // MedialStruct::id

  void print_medial_face() const;

  // v2int: endpoints on dual edges
  // and its surface fids (matching GEO::Mesh)
  //
  // generated from MedialSphere::pcell::edge_2endvertices in function
  // generate_medial_faces()
  //
  // NOTE: dual_edge_2endvertices size 'should be' 2 (1 dual_edge with 2
  // endpoints)
  std::vector<std::array<v2int, 2>>
      dual_edge_endpoints;  // dual edge on 3 powercells, may have multiple dual
                            // edges, one for each edge CC
  bool is_sorted_de = false;  // do not check geometry if false [no use]

  double importance = 2.;  // importance metric

  //////////////////////////////////////////////////////////////////////
  // There are 3 vertex pairs, each pairs share the same sphere
  // (must be on the opposite side of a same sphere)
  // the pairs are:
  // st0.v[0] - st1.v[0]
  // st0.v[1] - st1.v[2]
  // st0.v[2] - st1.v[1]
  //
  // Also 3 vertices in each SimpleTriangle are counter-clockwise
  //
  // Detail see Primitives::get_triangles_from_three_spheres()
  //////////////////////////////////////////////////////////////////////
  SimpleTriangle st[2];
  bool is_valid_st = false;

  // for relaxation
  double area;
  Vector3 circumcenter;
  Vector3 centroid;
  bool is_on_boundary = false;
};

// These tets need to be removed
class MedialTet {
 public:
  int tid;
  int dup_cnt;  // duplicated count

  // set has its own order
  std::set<int> vertices_;    // vertex list
  std::set<int> edges_;       // edge list
  std::set<int> faces_;       // facet list
  std::set<int> neigh_tets_;  // neighboring tetrahedrons
  bool is_deleted = false;

  Vector3 dual_vertex;
  bool is_dual_vertex_outside = false;
  double dual_sq_radius;
  std::array<int, 4> all_vs_tags;

  bool HasVertex(int vid) { return (vertices_.find(vid) != vertices_.end()); }
  bool HasEdge(int eid) { return (edges_.find(eid) != edges_.end()); }
  bool HasFace(int fid) { return (faces_.find(fid) != faces_.end()); }

  void print_medial_tet() const;
};

class MedialMesh {
 public:
  MedialMesh();
  MedialMesh(std::vector<MedialSphere>& all_medial_spheres);
  ~MedialMesh();

  // id matches MedialSphere::id
  std::vector<MedialSphere>* vertices = nullptr;
  std::vector<MedialEdge> edges;
  std::vector<MedialFace> faces;
  std::vector<MedialTet> tets;  // the one to remove!!!
  std::vector<MedialStruct> mstructure;

  // (mfid_min, mfid_max) -> int
  // -1: not visited
  // 0: not on the same sheet
  // 1: on the same sheet
  //
  // stored as a upper triangular matrix with size n,
  // given two sorted MedialFace:id indices (mfid_min, mfid_max),
  // we can get the index for the matrix by calling get_upper_tri_matrix_idx()
  // updated by RPD3D_Wrapper::cluster_mface_adj_type()
  std::vector<int> is_two_faces_on_the_same_sheet;

  int numSpheres_active;
  int numEdges_active;
  int numFaces_active;
  int numTets_active;

  // only for I/O
  int numEdges_no_dup;
  int numFaces_no_dup;

 public:
  void compute_face_simple_triangles_all();
  void compute_face_simple_triangles(int fid);  // helper
  void compute_faces_meta_data(int fid);
  void compute_faces_st_meta(
      const AABBWrapper& aabb_wrapper);  // for debug, fix_geo

 public:
  void clear();
  bool ValidVertex(int vid);
  bool Edge(int vid0, int vid1, int& eid);
  bool Face(const aint3& vset, int& fid);
  // bool Tet(const std::vector<int>& tvids, int& tid);

  void get_face_neighs(const int fid, std::set<int>& f_neighs) const;
  void get_edge_neighs(const int fid, std::set<int>& f_neighs) const;

  // for corner spheres
  // only create once
  bool is_corner_spheres_created = false;

  int create_edge(const int vid0, const int vid1,
                  int dup_cnt = 1 /*duplicated count*/,
                  int eid = -1 /*if not given then auto-increase*/);
  int create_face(const aint3& vvid, int dup_cnt = 1 /*duplicated count*/,
                  int fid = -1, /*if not given then auto-increase*/
                  bool is_auto_create_edge = true /*create edge if not exist*/);
  int create_tet(const std::array<int, 4>& tvids,
                 int dup_cnt = 1 /*duplicated count*/,
                 int tid = -1 /*if not given then auto-increase*/,
                 bool is_auto_create_face = false /*create face if not exist*/);
  void genearte_medial_spheres(std::vector<MedialSphere>& all_medial_spheres);
  void generate_medial_faces();
  void compute_common_diff_for_medge(const SurfaceMesh& sf_mesh,
                                     const std::pair<aint2, int>& me, int eid);
  void generate_medial_edges(const SurfaceMesh& sf_mesh,
                             bool is_update_tan_pls);
  void generate_medial_tets();

  // void sort_dual_edges_all(const GEO::Mesh& sf_mesh);
  // void sort_dual_edges(const GEO::Mesh& sf_mesh, MedialFace& mface);

  // for thinning
  int clear_all_tets();
  bool delete_vertex(const int vid);
  bool delete_edge(const int eid);
  bool delete_face(const int fid);
  bool delete_tet(const int tet_id);
  void check_and_store_unthin_tets_in_mat();

  // for medial structure
  void trace_medial_structure(bool is_debug);

  // for dup_cnt
  void validate_mmesh_dup_cnt();
  // For only edges with dup_cnt>1
  // edge must have dup_cnt > 1
  // face can be either dup_cnt > 1 or -1
  std::map<int, int> dup_e2f;
  // For only faces with dup_cnt>1
  // each dup_cntf(face) > 1 must has a dup_cnt(edge)>1
  // validated in validate_mmesh_dup_cnt()
  std::map<int, int> dup_f2e;
  void update_mmesh_dup_ef_map(bool is_debug = false);
  int get_edge_fid_min_importance_keep_connectivity(int eid);

  // for external/internal feature
  std::set<aint2> mat_intf_edges;
  std::set<aint2> mat_extf_edges;
};

int compute_Euler(const MedialMesh& mat);
