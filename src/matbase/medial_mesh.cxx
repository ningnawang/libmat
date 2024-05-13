#include "medial_mesh.h"

#include <assert.h>

void MedialFace::print_medial_face() const {
  printf("------ MAT face %d, is_deleted %d, is_on_boundary: %d\n", fid,
         is_deleted, is_on_boundary);
  // printf(
  //     "dual_segment len: %f\n",
  //     (dual_edge_endpoints[0].first -
  //     dual_edge_endpoints[1].first).length());
  printf("vertices_: (%d,%d,%d), importance: %f \n", vertices_[0], vertices_[1],
         vertices_[2], importance);
  if (!tets_.empty()) {
    printf("tets_: ");
    for (const auto tid : tets_) {
      printf("%d, ", tid);
    }
    printf("\n");
  }
}

void MedialTet::print_medial_tet() const {
  printf("------ MAT tet %d, is_deleted: %d\n", tid, is_deleted);
  if (!faces_.empty()) {
    printf("faces_: ");
    for (const auto fid : faces_) {
      printf("%d, ", fid);
    }
    printf("\n");
  }
}

/////////////////////////////////////////////
// MedialMesh Functions
/////////////////////////////////////////////
MedialMesh::MedialMesh() {
  numSpheres_active = numEdges_active = numFaces_active = numTets_active = 0;
  numEdges_no_dup = numFaces_no_dup = 0;
};

MedialMesh::MedialMesh(std::vector<MedialSphere>& all_medial_spheres) {
  vertices = &all_medial_spheres;
  numSpheres_active = numEdges_active = numFaces_active = numTets_active = 0;
  numEdges_no_dup = numFaces_no_dup = 0;
};

MedialMesh::~MedialMesh() { delete vertices; };

void MedialMesh::compute_face_simple_triangles_all() {
  for (int i = 0; i < faces.size(); i++) compute_face_simple_triangles(i);
}

void MedialMesh::clear() {
  if (vertices != nullptr) {
    for (uint vid = 0; vid < vertices->size(); vid++) {
      auto& mvertex = vertices->at(vid);
      mvertex.clear_mm();
    }
    vertices = nullptr;
  }
  edges.clear();
  faces.clear();
  tets.clear();
  mstructure.clear();
  numSpheres_active = numEdges_active = numFaces_active = numTets_active = 0;
  numEdges_no_dup = numFaces_no_dup = 0;

  // clear feature edges
  mat_intf_edges.clear();
  mat_extf_edges.clear();
}

int MedialMesh::clear_all_tets() {
  tets.clear();
  numTets_active = 0;

  // all adjacent info will also be cleaned
  for (int fid = 0; fid < faces.size(); fid++) {
    faces[fid].tets_.clear();
  }
}

bool MedialMesh::delete_vertex(const int vid) {
  assert(vertices != nullptr);
  if (vid < 0 || vid >= vertices->size()) return false;
  auto& mvertex = vertices->at(vid);
  if (mvertex.is_deleted) return true;
  mvertex.is_deleted = true;
  this->numSpheres_active -= mvertex.dup_cnt;

  // check all adjacencies are deleted
  assert(mvertex.edges_.empty());
  assert(mvertex.faces_.empty());
  return true;
}

bool MedialMesh::delete_edge(const int eid) {
  if (eid < 0 || eid >= edges.size()) return false;
  auto& edge = edges[eid];
  if (edge.is_deleted) return true;
  edge.is_deleted = true;
  numEdges_active -= edge.dup_cnt;
  numEdges_no_dup -= 1;
  // update vertices
  vertices->at(edge.vertices_[0]).edges_.erase(eid);
  if (vertices->at(edge.vertices_[0]).edges_.empty())
    delete_vertex(edge.vertices_[0]);
  vertices->at(edge.vertices_[1]).edges_.erase(eid);
  if (vertices->at(edge.vertices_[1]).edges_.empty())
    delete_vertex(edge.vertices_[1]);
  // check/update faces
  if (!edge.faces_.empty()) {
    printf(
        "[Delete_MAT_Edge] Edge %d should not be deleted,"
        "still exists in %ld faces\n ",
        eid, edge.faces_.size());
    assert(false);
  }
  return true;
}

bool MedialMesh::delete_face(const int fid) {
  if (fid < 0 || fid >= faces.size()) return false;
  auto& face = faces[fid];
  if (face.is_deleted) return true;
  face.is_deleted = true;
  numFaces_active -= face.dup_cnt;
  numFaces_no_dup -= face.dup_cnt;
  // update vertices
  for (const auto& vid : face.vertices_) vertices->at(vid).faces_.erase(fid);
  // update edges
  for (const auto& eid : face.edges_) {
    edges[eid].faces_.erase(fid);
    // if (edges[eid].faces_.empty()) delete_edge(eid);
  }
  // check/update tets
  if (!face.tets_.empty()) {
    printf(
        "[Delete_MAT_Face] face %d should not be deleted,"
        "still exists in %ld tet\n ",
        fid, face.tets_.size());
    assert(false);
  }
  return true;
}

bool MedialMesh::delete_tet(const int tet_id) {
  if (tet_id < 0 || tet_id > tets.size()) return false;

  auto& tet = tets[tet_id];
  if (tet.is_deleted) return true;
  tet.is_deleted = true;
  numTets_active -= tet.dup_cnt;
  // update faces
  for (const auto& fid : tet.faces_) faces[fid].tets_.erase(tet_id);
  return true;
}

void MedialMesh::compute_face_simple_triangles(int fid) {
  const int k = 3;  // every medial face is dim=3
  SimpleTriangle st0, st1;
  Vector3 pos[k];
  double radius[k];
  int count = 0;

  int i = 0;
  for (const auto vid : faces[fid].vertices_) {
    pos[i] = vertices->at(vid).center;
    radius[i] = vertices->at(vid).radius;
    i++;
  }

  int type = get_triangles_from_three_spheres(
      pos[0], radius[0], pos[1], radius[1], pos[2], radius[2], st0, st1);

  if (type == MedialFace::Type::INVALID) {
    faces[fid].is_valid_st = false;
    std::vector<int> fvids;
    for (const auto vid : faces[fid].vertices_) fvids.push_back(vid);
    // printf("[ERROR] MAT fid %d (%d,%d,%d) is invalid \n", fid, fvids[0],
    //        fvids[1], fvids[2]);
    // assert(false);
  } else {
    faces[fid].st[0] = st0;
    faces[fid].st[1] = st1;
    faces[fid].is_valid_st = true;
  }
}

// store circumcenter, area and centroid for medial faces
void MedialMesh::compute_faces_meta_data(int fid) {
  MedialFace& mface = faces.at(fid);
  std::array<Vector3, 3> mface_vs_pos;
  int extf_cnt = 0;
  FOR(i, 3) {
    const auto& ms = vertices->at(mface.vertices_[i]);
    mface_vs_pos[i] = ms.center;
    if (ms.is_on_extf()) extf_cnt++;
  }
  // update boundary
  if (extf_cnt > 1) mface.is_on_boundary = true;
  FOR(i, 3) {
    const auto& me = edges.at(mface.edges_[i]);
    if (me.faces_.size() == 1) mface.is_on_boundary = true;
  }
  mface.circumcenter = GEO::Geom::triangle_circumcenter(
      mface_vs_pos[0], mface_vs_pos[1], mface_vs_pos[2]);
  mface.area = GEO::Geom::triangle_area(mface_vs_pos[0], mface_vs_pos[1],
                                        mface_vs_pos[2]);
  mface.centroid =
      get_triangle_centroid(mface_vs_pos[0], mface_vs_pos[1], mface_vs_pos[2]);
}

void MedialMesh::compute_faces_st_meta(const AABBWrapper& aabb_wrapper) {
  for (auto& mface : faces) {
    if (!mface.is_valid_st) continue;
    for (int i = 0; i < 2; i++) {
      mface.st[i].centroid = get_triangle_centroid(
          mface.st[i].v[0], mface.st[i].v[1], mface.st[i].v[2]);

      // bool is_int = aabb_wrapper.get_ray_nearest_intersection(
      //     mface.st[i].centroid, mface.st[i].normal,
      //     mface.st[i].nearest_point, mface.st[i].nearest_sf_fid);
      // if (!is_int || mface.st[i].nearest_sf_fid == -1) {
      //   printf(
      //       "[SlabCheck] mface %d with mspheres (%d,%d,%d) has st %d has no "
      //       "intersection to surface??? check this \n",
      //       mface.fid, mface.vertices_[0], mface.vertices_[1],
      //       mface.vertices_[2], i);
      //   continue;
      // }
      // mface.st[i].dist_to_sf =
      //     (mface.st[i].centroid - mface.st[i].nearest_point).length();

      double sq_dist;
      mface.st[i].nearest_sf_fid = aabb_wrapper.get_nearest_point_on_sf(
          mface.st[i].centroid, mface.st[i].nearest_point, sq_dist);
      mface.st[i].dist_to_sf = std::sqrt(sq_dist);
    }  // for st
  }  // for mfaces
}

bool MedialMesh::ValidVertex(int vid) {
  if (vid > vertices->size()) return false;
  return true;
}

bool MedialMesh::Edge(int vid0, int vid1, int& eid) {
  if (!ValidVertex(vid0) || !ValidVertex(vid1)) return false;

  for (auto si = vertices->at(vid0).edges_.begin();
       si != vertices->at(vid0).edges_.end(); si++) {
    if (edges[*si].HasVertex(vid1)) {
      eid = *si;
      return true;
    }
  }
  return false;
}

bool MedialMesh::Face(const aint3& vset, int& fid) {
  for (auto si = vset.begin(); si != vset.end(); si++)
    if (!ValidVertex(*si)) return false;

  int vid0 = *(vset.begin());
  for (auto si = vertices->at(vid0).faces_.begin();
       si != vertices->at(vid0).faces_.end(); si++)
    if (faces[*si].vertices_ == vset) {
      fid = *si;
      return true;
    }
  return false;
}

void MedialMesh::get_face_neighs(const int fid, std::set<int>& f_neighs) const {
  f_neighs.clear();
  FOR(leid, 3) {
    int eid = faces.at(fid).edges_.at(leid);
    for (int nfid : edges.at(eid).faces_) {
      if (nfid == fid) continue;
      f_neighs.insert(nfid);
    }
  }
}

void MedialMesh::get_edge_neighs(const int eid, std::set<int>& e_neighs) const {
  e_neighs.clear();
  FOR(lvid, 2) {
    int vid = edges.at(eid).vertices_.at(lvid);
    for (int neid : vertices->at(vid).edges_) {
      if (neid == eid) continue;
      e_neighs.insert(neid);
    }
  }
}

int MedialMesh::create_edge(const int vid0, const int vid1, int dcnt,
                            int eid  // if not given then auto-increase
) {
  // edge already exists
  // if (Edge(vid0, vid1, eid)) return eid;
  assert(!Edge(vid0, vid1, eid));  // do not handle case when edge exists
  if (eid == -1) eid = edges.size();
  MedialEdge bep;
  bep.eid = eid;  // give each mat point an index
  bep.vertices_[0] = vid0;
  bep.vertices_[1] = vid1;
  bep.dup_cnt = dcnt;
  vertices->at(vid0).edges_.insert(eid);
  vertices->at(vid1).edges_.insert(eid);
  edges.push_back(bep);
  numEdges_active += dcnt;
  numEdges_no_dup += 1;
  return eid;
}

int MedialMesh::create_face(const aint3& vvid, int dcnt,
                            int fid,  // if not given then auto-increase
                            bool is_auto_create_edge) {
  aint3 vvid_sorted = get_sorted(vvid);
  // if (Face(vvid_sorted, fid)) return fid;
  assert(!Face(vvid_sorted, fid));  // do not handle case when face exists
  if (fid == -1) fid = faces.size();
  MedialFace bfp;
  bfp.fid = fid;
  bfp.vertices_ = vvid_sorted;
  std::set<int> eids_set;
  for (const auto& v1 : vvid_sorted) {
    for (const auto& v2 : vvid_sorted) {
      if (v1 == v2) continue;
      int one_eid = -1;
      if (Edge(v1, v2, one_eid)) {
        // nothing
      } else {
        if (!is_auto_create_edge) {
          // maybe degenerated case?
          printf(
              "[MAT] face %d not contain an edge with -1 index, vvid_sorted: "
              "(%d,%d,%d), edge (%d,%d) not exits\n",
              fid, vvid_sorted[0], vvid_sorted[1], vvid_sorted[2], v1, v2);
          assert(false);
          // return -1;
        }
        one_eid = create_edge(v1, v2, dcnt);
      }
      if (one_eid == -1) {
        printf(
            "[MAT] face not contain an edge with -1 index, vvid_sorted: "
            "(%d,%d,%d)\n",
            vvid_sorted[0], vvid_sorted[1], vvid_sorted[2]);
        assert(false);
      }
      eids_set.insert(one_eid);
    }
  }
  if (eids_set.size() != 3) {
    printf(
        "[MAT] face not contain 3 edges, vvid_sorted: (%d,%d,%d), edges size: "
        "%ld\n",
        vvid_sorted[0], vvid_sorted[1], vvid_sorted[2], eids_set.size());
    assert(false);
  }

  for (const auto& v : vvid_sorted) {
    vertices->at(v).faces_.insert(fid);
  }
  std::vector<int> eids;
  for (const auto& e : eids_set) {
    eids.push_back(e);
    edges[e].faces_.insert(fid);
  }
  bfp.edges_ = {{eids[0], eids[1], eids[2]}};
  bfp.dup_cnt = dcnt;
  faces.push_back(bfp);
  numFaces_active += dcnt;
  numFaces_no_dup += 1;
  return fid;
}

int MedialMesh::create_tet(const std::array<int, 4>& vids, int dcnt, int tid,
                           bool is_auto_create_face) {
  if (tid == -1) tid = tets.size();
  MedialTet tet;
  bool is_debug = false;
  tet.tid = tid;
  tet.dup_cnt = dcnt;
  std::vector<int> fids;
  // save faces
  for (int i = 0; i < 4; i++) {
    for (int j = i + 1; j < 4; j++) {
      for (int k = j + 1; k < 4; k++) {
        aint3 vvid;
        vvid[0] = vids[i];
        vvid[1] = vids[j];
        vvid[2] = vids[k];
        int fid = -1;
        if (Face(vvid, fid)) {  // check if mat face exist
          // nothing
        } else {
          if (!is_auto_create_face) {
            printf(
                "[MAT] tet not contain a face with -1 index, vids: "
                "(%d,%d,%d,%d), face (%d,%d,%d) not exits\n",
                vids[0], vids[1], vids[2], vids[3], vvid[0], vvid[1], vvid[2]);
            // assert(false);
            return -1;
          }
          // creating new face if not exist
          fid = create_face(vvid, -1, true);
        }
        const MedialFace& face = faces[fid];
        tet.vertices_.insert(face.vertices_.begin(), face.vertices_.end());
        tet.edges_.insert(face.edges_.begin(), face.edges_.end());
        tet.faces_.insert(fid);
        tet.neigh_tets_.insert(face.tets_.begin(), face.tets_.end());
        fids.push_back(fid);
      }  // k
    }  // j
  }  // i

  if (fids.size() != 4) {
    printf("[MAT TET] error creating tet %d, fids size not 4: %ld \n", tid,
           fids.size());
    assert(false);
  }
  if (is_debug)
    printf("tet %d has 4 vids: (%d,%d,%d,%d)\n", tid, vids[0], vids[1], vids[2],
           vids[3]);
  // only insert tid to face when tet is truly created.
  for (const auto& fid : fids) {
    MedialFace& face = faces[fid];
    face.tets_.insert(tid);
  }
  tets.push_back(tet);
  numTets_active += dcnt;
  return tid;
}

// dual to powercell edges
void MedialMesh::generate_medial_faces() {
  assert(vertices != nullptr);
  printf("[MedialMesh] generating medial faces ...\n");
  std::map<aint3, std::vector<std::array<v2int, 2>>>
      mfs_dual_edges;            // aint3 is sorted
  std::map<aint3, int> mfs_cnt;  // aint3 is sorted
  for (auto& msphere : *vertices) {
    if (msphere.is_deleted) continue;
    // [neigh_id_min, neigh_id_max] -> { set of cell_ids in one edge CC }
    for (const auto& pair : msphere.pcell.edge_cc_cells) {
      const aint2& neigh_min_max = pair.first;
      const std::vector<std::array<aint2, 2>>& edge_endvertices_vec =
          msphere.pcell.edge_2endvertices.at(neigh_min_max);
      aint3 mf = {{msphere.id, neigh_min_max[0], neigh_min_max[1]}};
      // // skip pin invisible spheres for concave edges
      // if (vertices->at(mf[1]).is_on_ce_pin() ||
      //     vertices->at(mf[2]).is_on_ce_pin())
      //   continue;
      std::sort(mf.begin(), mf.end());

      // three spheres must have same number of EdgeCC
      // after the topo fix!
      if (mfs_cnt.find(mf) != mfs_cnt.end()) {
        if (pair.second.size() != mfs_cnt[mf]) {
          printf("mf (%d,%d,%d) has msphere %d, edgeCC %d != %d \n", mf[0],
                 mf[1], mf[2], msphere.id, pair.second.size(), mfs_cnt[mf]);
          // for (const auto cids : pair.second) {
          //   printf("(");
          //   for (const auto cid : cids) printf("%d,", cid);
          //   printf("), ");
          // }
          // printf("]\n");
          assert(false);
        }
        continue;  // do nothing
      }

      // [s1, s2, s3] has been count, no need to recal
      if (mfs_dual_edges.find(mf) != mfs_dual_edges.end()) continue;
      for (const auto dual_edge : edge_endvertices_vec) {
        mfs_dual_edges[mf].push_back(
            {{msphere.pcell.vertex_2pos.at(dual_edge[0]),
              msphere.pcell.vertex_2pos.at(dual_edge[1])}});
      }
      // each edgeCC dual to one medial face
      mfs_cnt[mf] = pair.second.size();
    }
  }

  faces.clear();
  for (auto& pair : mfs_cnt) {
    const aint3& seed3 = pair.first;
    const int dup_cnt = pair.second;
    const std::vector<std::array<v2int, 2>>& edge_endpoints =
        mfs_dual_edges.at(seed3);
    // NOTE: must set is_auto_create_edge=true
    //       otherwise some edges may not have adjacent face
    //       which will make exporting PLY not happy
    int fid = create_face(seed3, dup_cnt, -1, true /*is_auto_create_edge*/);
    if (fid == -1) continue;
    faces[fid].dual_edge_endpoints = edge_endpoints;
    compute_face_simple_triangles(fid);
    compute_faces_meta_data(fid);
    // if (dup_cnt > 1) {
    //   printf("++++ medial face %d (%d,%d,%d) has dup_cnt %d\n", fid,
    //   seed3[0],
    //          seed3[1], seed3[2], dup_cnt);
    // }
  }
}

void MedialMesh::genearte_medial_spheres(
    std::vector<MedialSphere>& _all_medial_spheres) {
  this->vertices = &_all_medial_spheres;
  printf("[MedialMesh] generating medial spheres ...\n");
  for (auto& msphere : *vertices) {
    if (msphere.is_deleted) continue;
    msphere.dup_cnt = msphere.pcell.cc_cells.size();
    this->numSpheres_active += msphere.dup_cnt;
  }
}

// helper function for generate_medial_edges()
void MedialMesh::compute_common_diff_for_medge(const SurfaceMesh& sf_mesh,
                                               const std::pair<aint2, int>& me,
                                               int eid) {
  // update sf_mesh fids covered info
  auto& medge = this->edges.at(eid);
  const MedialSphere& A = this->vertices->at(me.first[0]);
  const MedialSphere& B = this->vertices->at(me.first[1]);
  medge.is_vs_swapped = false;  // reset
  A_B_spheres_common_diff_tangent_info_from_surface(
      sf_mesh, A, B, medge.common_tan_pls, medge.v1_non_common_tan_pls,
      medge.v2_non_common_tan_pls);

  // mark medial edge if on internal feature
  if ((A.is_on_intf() && B.is_on_intf()) &&
      (medge.v1_non_common_tan_pls.empty() ||
       medge.v2_non_common_tan_pls.empty())) {
    medge.is_intf = true;
    // printf("[MedialEdge] medge (%d,%d) is on intf\n", medge.vertices_[0],
    //        medge.vertices_[1]);
  } else {
    // printf("[MedialEdge] medge (%d,%d) is NOT on intf\n",
    // medge.vertices_[0],
    //        medge.vertices_[1]);
  }

  // mark medial edge if on external feature
  if (A.is_on_se() && B.is_on_se() && A.se_line_id == B.se_line_id &&
      (medge.v1_non_common_tan_pls.empty() ||
       medge.v2_non_common_tan_pls.empty())) {
    medge.is_extf = true;
  } else if (A.is_on_corner() && B.is_on_se() &&
             A.corner_fls.find(B.se_line_id) != A.corner_fls.end()) {
    medge.is_extf = true;
  } else if (B.is_on_corner() && A.is_on_se() &&
             B.corner_fls.find(A.se_line_id) != B.corner_fls.end()) {
    medge.is_extf = true;
  } else if (A.is_on_corner() && B.is_on_corner()) {
    medge.is_extf = true;
  }
}

void MedialMesh::generate_medial_edges(const SurfaceMesh& sf_mesh,
                                       bool is_compute_common_diff_sfids) {
  assert(vertices != nullptr);
  printf("[MedialMesh] generating medial edges ...\n");
  // aint2 is sorted -> num of created medial edge
  std::map<aint2, int> mes_cnt;
  for (auto& msphere : *vertices) {
    if (msphere.is_deleted) continue;
    // neigh_id -> { set of cell_ids in one facet CC }
    for (const auto& pair : msphere.pcell.facet_cc_cells) {
      const int& seed_neigh_id = pair.first;
      // if (vertices->at(seed_neigh_id).is_on_ce_pin()) continue;
      aint2 me = {{msphere.id, seed_neigh_id}};
      std::sort(me.begin(), me.end());

      // <seed_id, neigh_seed_id> has been checked, no need to count again
      if (mes_cnt.find(me) != mes_cnt.end()) {
        // two spheres must have same number of FacetCC
        // after the topo fix!
        if (pair.second.size() != mes_cnt[me]) {
          printf(
              "medial edge with [%d,%d] has sphere %d FaceCC size %d != %d\n",
              me[0], me[1], msphere.id, pair.second.size(), mes_cnt[me]);
          // for (const auto cids : pair.second) {
          //   printf("(");
          //   for (const auto cid : cids) printf("%d,", cid);
          //   printf(")");
          // }
          // printf("]\n");
          assert(false);
        }
        continue;
      }
      // each facetCC dual to one medial edge
      mes_cnt[me] = pair.second.size();
    }
  }

  // Note: cannot clear, generate_medial_faces() created some edges!
  // edges.clear();
  for (auto& me : mes_cnt) {
    // will abort and error if already exists
    int eid = create_edge(me.first[0], me.first[1], me.second);
    // if (me.second > 1) {
    //   printf("++++ medial edge %d (%d,%d) has dup_cnt %d\n", eid,
    //   me.first[0], me.first[1], me.second);
    // }

    if (is_compute_common_diff_sfids) {
      // for computing common/diff tangent planes
      compute_common_diff_for_medge(sf_mesh, me, eid);
    }
  }
}

void MedialMesh::generate_medial_tets() {
  assert(vertices != nullptr);
  printf("[MedialMesh] generating medial tets ...\n");
  std::map<aint4, int> mtets_cnt;  // aint4 is sorted
  for (auto& msphere : *vertices) {
    if (msphere.is_deleted) continue;
    // [cur_id, neigh_id1, neigh_id2, neigh_id3]
    // and we use sorted [neigh_id1, neigh_id2, neigh_id3] to represent this
    // unique vertex
    //
    // 1. If neigh_id1 is -1, then this powercell vertex is on surface
    // 2. If all neigh_id(i) > 0, then vid is inside the shape. If this vertex
    // is a powercell vertex (not just a convex cell vertex), then it is dual
    // to a medial tet
    for (const auto& pair : msphere.pcell.vertex_2id) {
      const aint3 neigh_ids = pair.second;
      if (neigh_ids[0] < 0 || neigh_ids[1] < 0 || neigh_ids[2] < 0) continue;
      // printf("[MedialMesh] processing msphere %d, neigh_ids: (%d,%d,%d)\n",
      //        msphere.id, neigh_ids[0], neigh_ids[1], neigh_ids[2]);
      // if (vertices->at(neigh_ids[0]).is_on_ce_pin()) continue;
      // if (vertices->at(neigh_ids[1]).is_on_ce_pin()) continue;
      // if (vertices->at(neigh_ids[2]).is_on_ce_pin()) continue;
      aint4 mtet = {{msphere.id, neigh_ids[0], neigh_ids[1], neigh_ids[2]}};
      std::sort(mtet.begin(), mtet.end());
      if (mtets_cnt.find(mtet) != mtets_cnt.end()) continue;
      mtets_cnt[mtet] += 1;
    }
  }

  // printf("mtets_cnt: %ld\n", mtets_cnt.size());
  tets.clear();
  for (auto& mtet : mtets_cnt) {
    // will skip if already exists
    create_tet(mtet.first, mtet.second, -1, false /*is_auto_create_face*/);
  }
}

// void MedialMesh::sort_dual_edges_all(const GEO::Mesh& sf_mesh) {
//   for (auto& mface : faces) {
//     if (!mface.is_valid_st) continue;
//     sort_dual_edges(sf_mesh, mface);
//   }
// }

// // sort the mface.dual_edge_endpoints based on the order of simple
// triangles
// // using normals
// void MedialMesh::sort_dual_edges(const GEO::Mesh& sf_mesh, MedialFace&
// mface)
// {
//   assert(mface.is_valid_st);
//   // must pass the topology check is_to_fix_edge_cc()
//   // only 1 dual_edge with 2 endpoints， at least 1 endpoint on surface
//   // assert(mface.dual_edge_endpoints.size() >= 1 &&
//   //        mface.dual_edge_endpoints.size() < 3);
//   if (mface.dual_edge_endpoints.size() < 1 &&
//       mface.dual_edge_endpoints.size() >= 3)
//     return;  // TODO: check this
//   std::vector<v2int> dual_edge_sorted(
//       2, std::make_pair<Vector3, int>(Vector3(0, 0, 0), -1));
//   for (auto& endpoint : mface.dual_edge_endpoints) {
//     int surf_fid = endpoint.second;
//     Vector3 point_normal = get_mesh_facet_normal(sf_mesh, surf_fid);
//     for (int i = 0; i < 2; i++) {
//       if (is_vector_same_direction(point_normal, mface.st[i].normal,
//                                    EPS_DEGREE_90)) {
//         dual_edge_sorted[i] = endpoint;
//         break;
//       }
//     }
//   }
//   mface.dual_edge_endpoints = dual_edge_sorted;
//   mface.is_sorted_de = true;
// }

// for debug
void MedialMesh::check_and_store_unthin_tets_in_mat() {
  clear_all_tets();

  // update mat_v_neighbors map
  std::map<int, std::unordered_set<int>> mat_v_neighbors;
  for (const auto face : faces) {
    aint3 f_vs = face.vertices_;
    for (int i = 0; i < 3; i++) {
      mat_v_neighbors[f_vs[i]].insert(f_vs[(i + 1) % 3]);
      mat_v_neighbors[f_vs[i]].insert(f_vs[(i + 2) % 3]);
    }
  }

  // now we have mat_v_neighbors map,
  // for every face with 3 vs, we check if another mat vertex
  // is adjacent to all 3 vs. If so, then tet exits in mat, not thin
  std::set<std::array<int, 4>> possible_tets;  // to avoid dup tets
  for (const auto face : faces) {
    std::vector<int> intersections;
    set_intersection<int>(mat_v_neighbors.at(face.vertices_[0]),
                          mat_v_neighbors.at(face.vertices_[1]),
                          mat_v_neighbors.at(face.vertices_[2]), intersections);
    std::vector<int> tet_vids;
    for (const auto& v_inter : intersections) {
      tet_vids.clear();
      tet_vids.push_back(v_inter);
      FOR(i, 3) tet_vids.push_back(face.vertices_[i]);
      std::sort(tet_vids.begin(), tet_vids.end());
      assert(tet_vids.size() == 4);
      possible_tets.insert(
          {{tet_vids[0], tet_vids[1], tet_vids[2], tet_vids[3]}});
    }  // for intersections
  }  // for faces

  // printf("[UNTHIN CHECK] found %ld possible_tets\n", possible_tets.size());

  for (const auto& one_tet_vids : possible_tets) {
    // will not create if any face not exists
    create_tet(one_tet_vids, -1, false);
  }
  printf("[UNTHIN CHECK] found %ld tets in MAT\n", tets.size());
}

bool is_add_neigh_eid(const MedialMesh& mmesh, int cur_eid, int n_eid) {
  // check if add n_eid
  const auto& cur_medge = mmesh.edges.at(cur_eid);
  const auto& n_medge = mmesh.edges.at(n_eid);

  // find current sphere and next sphere
  int common_vid = -1, nxt_vid = -1;
  for (int n_vid : n_medge.vertices_) {
    bool is_non_common = true;
    for (int cur_vid : cur_medge.vertices_) {
      if (cur_vid == n_vid) {
        common_vid = n_vid;
        is_non_common = false;
      }
    }
    if (!is_non_common) continue;  // found common
    nxt_vid = n_vid;               // non-common
  }
  assert(common_vid != -1 && nxt_vid != -1);
  const auto& cur_msphere = mmesh.vertices->at(common_vid);
  const auto& nxt_msphere = mmesh.vertices->at(nxt_vid);

  const std::vector<TangentPlane>& common_tan_pls = n_medge.common_tan_pls;
  const std::vector<TangentPlane>* A_non_common_tan_pls =
      &n_medge.v1_non_common_tan_pls;
  const std::vector<TangentPlane>* B_non_common_tan_pls =
      &n_medge.v2_non_common_tan_pls;
  if (common_vid == n_medge.vertices_[1]) {
    A_non_common_tan_pls = &(n_medge.v2_non_common_tan_pls);
    B_non_common_tan_pls = &(n_medge.v1_non_common_tan_pls);
  }

  printf("checking next edge %d (%d,%d)\n", n_eid, common_vid, nxt_vid);

  // check if nxt_msphere
  // 1. same type as cur_msphere
  if (cur_msphere.is_on_intf() && !nxt_msphere.is_on_intf()) return false;
  if (cur_msphere.is_on_extf() && !nxt_msphere.is_on_extf()) return false;
  if (cur_msphere.is_on_sheet() && !nxt_msphere.is_on_sheet()) return false;

  printf(
      "[MedialStruc] n_eid %d (%d,%d) has common_tan_pls %zu, "
      "A_non_common_tan_pls: %zu \n",
      n_eid, common_vid, nxt_vid, common_tan_pls.size(),
      A_non_common_tan_pls->size());
  printf("--- common\n");
  for (const auto& tan_pl : common_tan_pls) tan_pl.print_info();
  printf("--- A_non_common\n");
  for (const auto& tan_pl : *A_non_common_tan_pls) tan_pl.print_info();
  printf("--- B_non_common\n");
  for (const auto& tan_pl : *B_non_common_tan_pls) tan_pl.print_info();

  // 2. coveres cur_sphere (B covers A)
  if (!common_tan_pls.empty() && A_non_common_tan_pls->empty()) return true;
  return false;
}

bool is_add_neigh_fid(const MedialMesh& mmesh, int cur_fid, int nxt_fid) {
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

  // union
  // keep both f1 and f2
  auto check_and_union_two_v2fid_groups =
      [&](std::vector<std::vector<v2int>>& f1_v2fid_groups,
          const std::vector<std::vector<v2int>>& f2_v2fid_groups) {
        // fetch f1 fids and their groups
        std::map<int, int> f1_saved_fids;
        for (int group_id = 0; group_id < f1_v2fid_groups.size(); group_id++) {
          const auto& f1_v2fid_group = f1_v2fid_groups.at(group_id);
          for (const auto f1_v2fid : f1_v2fid_group)
            f1_saved_fids[f1_v2fid.second] = group_id;
        }

        // merge f2 to f1
        for (const auto& f2_v2fid_group : f2_v2fid_groups) {
          // if already saved, then get existing group
          int exist_f1_group_id = -1;
          for (const auto& f2_v2fid : f2_v2fid_group) {
            if (f1_saved_fids.find(f2_v2fid.second) == f1_saved_fids.end())
              continue;
            exist_f1_group_id = f1_saved_fids.at(f2_v2fid.second);
            break;
          }

          if (exist_f1_group_id != -1) {
            // merge with current f1_v2fid_group
            auto& f1_v2fid_group = f1_v2fid_groups.at(exist_f1_group_id);
            merge_two_v2fid_vec(f1_v2fid_group, f2_v2fid_group);
          } else {
            // push to f1_v2fid_groups
            f1_v2fid_groups.push_back(f2_v2fid_group);
            exist_f1_group_id = f1_v2fid_groups.size() - 1;
          }
          assert(exist_f1_group_id != -1);
        }
      };

  // join
  // only keep the (f1 join f2)'s group
  auto check_and_join_two_v2fid_groups =
      [&](const std::vector<std::vector<v2int>>& f1_v2fid_groups,
          const std::vector<std::vector<v2int>>& f2_v2fid_groups) {
        std::vector<std::vector<v2int>> join_v2fid_groups;
        // fetch f1 fids and their groups
        std::map<int, int> f1_saved_fids;
        for (int f1_group_id = 0; f1_group_id < f1_v2fid_groups.size();
             f1_group_id++) {
          const auto& f1_v2fid_group = f1_v2fid_groups.at(f1_group_id);
          for (const auto f1_v2fid : f1_v2fid_group)
            f1_saved_fids[f1_v2fid.second] = f1_group_id;
        }

        // join f2 with f1
        for (const auto& f2_v2fid_group : f2_v2fid_groups) {
          // if already saved, then get existing group
          int exist_f1_group_id = -1;
          for (const auto& f2_v2fid : f2_v2fid_group) {
            if (f1_saved_fids.find(f2_v2fid.second) == f1_saved_fids.end())
              continue;
            exist_f1_group_id = f1_saved_fids.at(f2_v2fid.second);
            break;
          }
          // skip if not in f1
          if (exist_f1_group_id == -1) continue;
          // merge with current f1, save to cur_v2fid_group
          auto cur_v2fid_group = f1_v2fid_groups.at(exist_f1_group_id);
          merge_two_v2fid_vec(cur_v2fid_group, f2_v2fid_group);
          join_v2fid_groups.push_back(cur_v2fid_group);
        }
        return join_v2fid_groups;
      };

  // join, not union
  auto get_mface_join_covered_fid_groups =
      [&](const int fid, std::vector<std::vector<v2int>>& cur_merged_fid_groups,
          bool is_skip_intf) {
        FOR(lv, 3) {
          const auto& msphere =
              mmesh.vertices->at(mmesh.faces.at(fid).vertices_[lv]);
          if (is_skip_intf && msphere.is_on_intf()) continue;
          if (cur_merged_fid_groups.empty()) {
            cur_merged_fid_groups = msphere.covered_sf_fids_in_group;
            continue;
          }
          // join
          cur_merged_fid_groups = check_and_join_two_v2fid_groups(
              cur_merged_fid_groups, msphere.covered_sf_fids_in_group);
        }
      };

  auto print_covered_sf_fids_in_group =
      [](const std::vector<std::vector<v2int>>& covered_sf_fids_in_group) {
        printf("covered_sf_fids_in_group: \n");
        for (const auto& one_fid_group : covered_sf_fids_in_group) {
          printf("(");
          for (const auto& fid : one_fid_group) {
            printf("%d ", fid.second);
          }
          printf(")\n");
        }
      };

  bool is_debug = true;
  if (is_debug) printf("checking cur_fid %d, nxt_fid %d\n", cur_fid, nxt_fid);
  // find common vs
  std::set<int> cur_vs, nxt_vs;
  std::vector<int> common_vs;
  cur_vs.insert(mmesh.faces.at(cur_fid).vertices_.begin(),
                mmesh.faces.at(cur_fid).vertices_.end());
  nxt_vs.insert(mmesh.faces.at(nxt_fid).vertices_.begin(),
                mmesh.faces.at(nxt_fid).vertices_.end());
  set_intersection(cur_vs, nxt_vs, common_vs);
  assert(common_vs.size() == 2);

  if (is_debug)
    printf("found common_vs: (%d,%d)\n", common_vs[0], common_vs[1]);

  ////////////////////////////////////
  // Step 1: check if common edge on intf
  // find common eid
  std::vector<int> common_edge;
  set_intersection(mmesh.vertices->at(common_vs[0]).edges_,
                   mmesh.vertices->at(common_vs[1]).edges_, common_edge);
  assert(common_edge.size() == 1);
  // if on intf, then skip
  const auto& medge = mmesh.edges.at(common_edge[0]);
  if (medge.is_intf) return false;

  ////////////////////////////////////
  // Step 2: check joint of three spheres
  // check common covered sf_mesh fids of cur_fid and nxt_fid
  std::vector<std::vector<v2int>> cur_merged_fid_groups, nxt_merged_fid_groups;
  get_mface_join_covered_fid_groups(cur_fid, cur_merged_fid_groups,
                                    false /*is_skip_intf*/);
  get_mface_join_covered_fid_groups(nxt_fid, nxt_merged_fid_groups,
                                    false /*is_skip_intf*/);

  if (is_debug)
    printf("cur_merged_fid_groups: %zu, nxt_merged_fid_groups: %zu \n",
           cur_merged_fid_groups.size(), nxt_merged_fid_groups.size());

  // dunno when
  if (cur_merged_fid_groups.empty() && nxt_merged_fid_groups.empty())
    return true;
  // dunno when
  if (cur_merged_fid_groups.empty() || nxt_merged_fid_groups.empty())
    return false;

  // if size not match, then false
  if (cur_merged_fid_groups.size() != nxt_merged_fid_groups.size())
    return false;

  // check if cur_merged_fid_groups and nxt_merged_fid_groups
  // share the same sf_mesh fids
  std::vector<std::vector<v2int>> tmp = cur_merged_fid_groups;
  check_and_union_two_v2fid_groups(cur_merged_fid_groups,
                                   nxt_merged_fid_groups);
  if (cur_merged_fid_groups.size() != tmp.size()) {
    if (is_debug) {
      print_covered_sf_fids_in_group(tmp);
      print_covered_sf_fids_in_group(cur_merged_fid_groups);
    }
    return false;
  }
  return true;
}

void MedialMesh::trace_medial_structure(bool is_debug) {
  auto get_medial_structure = [&](int fid) {
    const auto& msphere1 = this->vertices->at(this->faces.at(fid).vertices_[0]);
    const auto& msphere2 = this->vertices->at(this->faces.at(fid).vertices_[1]);
    const auto& msphere3 = this->vertices->at(this->faces.at(fid).vertices_[2]);

    if (msphere1.is_on_intf() && msphere2.is_on_intf() && msphere3.is_on_intf())
      return MedialType::SINGULAR;
    return MedialType::SHEET;
    // printf("ERROR: unkown msphere1 %d type %d, msphere2 %d type %d\n",
    //        msphere1.id, msphere1.type, msphere2.id, msphere2.type);
    // msphere1.print_info();
    // msphere2.print_info();
    // assert(false);
  };

  printf("[MedialStruc] start tracing medial structure ...\n");
  this->mstructure.clear();
  std::set<int> visited;
  std::queue<int> visit_queue;
  // for (const auto& mface : faces) {
  while (visited.size() != this->faces.size()) {
    int rand_fid = RANDOM_INT(0, this->faces.size());
    const auto& mface = this->faces.at(rand_fid);
    if (mface.is_deleted) continue;
    if (is_debug)
      printf("[MedialStruc] checking mface %d (%d,%d,%d), visited: %zu\n",
             mface.fid, mface.vertices_[0], mface.vertices_[1],
             mface.vertices_[2], visited.size());
    if (visited.find(mface.fid) != visited.end()) continue;

    // create a new MedialStruct
    MedialStruc mstruct(this->mstructure.size(), MedialType::MUNKOWN);
    mstruct.type = get_medial_structure(mface.fid);

    // expand a new medial structure
    visit_queue.push(mface.fid);
    while (!visit_queue.empty()) {
      int cur_fid = visit_queue.front();
      visit_queue.pop();
      assert(!this->faces.at(cur_fid).is_deleted);
      if (visited.find(cur_fid) != visited.end()) continue;
      visited.insert(cur_fid);
      mstruct.m_face_ids.insert(cur_fid);
      this->faces.at(cur_fid).mstruc_id = mstruct.id;

      // check if adding next edge to queue
      for (int eid : this->faces.at(cur_fid).edges_) {
        for (int nxt_fid : this->edges.at(eid).faces_) {
          if (nxt_fid == cur_fid) continue;
          if (this->faces.at(nxt_fid).is_deleted) continue;
          if (visited.find(nxt_fid) != visited.end()) continue;
          // check cur_fid and nxt_fid
          bool is_add = is_add_neigh_fid(*this, cur_fid, nxt_fid);
          if (is_debug)
            printf("[MedialStruc] cur_fid %d nxt_fid %d is_add %d\n", cur_fid,
                   nxt_fid, is_add);
          if (!is_add) continue;
          visit_queue.push(nxt_fid);
        }
      }
    }  // while visit_queu

    // print
    if (is_debug) {
      printf("[MedialStruc] saved mstruc %d type %d, fid: (", mstruct.id,
             mstruct.type);
      for (const auto fid : mstruct.m_face_ids) printf("%d ", fid);
      printf("), visited %zu/%zu \n", visited.size(), faces.size());
    }

    // save MedialStruc
    this->mstructure.push_back(mstruct);
  }  // for faces

  if (is_debug)
    printf("[MedialStruc] mstruc size %zu, visited: %zu, faces: %zu\n",
           this->mstructure.size(), visited.size(), faces.size());
}

// Validate:
// 1. each dup_cnt(face) contains > 1 dup_cnt(edge
void MedialMesh::validate_mmesh_dup_cnt() {
  std::set<int> all_dup_edges;

  // check mmesh edges dup_cnt
  int num_dup_edges = 0, num_dup_faces = 0;
  for (const auto& medge : this->edges) {
    if (medge.dup_cnt == 1) continue;
    num_dup_edges++;
    all_dup_edges.insert(medge.eid);

    // check if one neighboring mface dup_cnt
    int nfid_cnt = 0;
    for (const int nfid : medge.faces_) {
      if (this->faces.at(nfid).dup_cnt > 1) nfid_cnt++;
    }
    if (nfid_cnt < 1) {
      printf("[validate_mmesh] medge %d (%d,%d) has nfid_cnt: %d\n", medge.eid,
             medge.vertices_[0], medge.vertices_[1], nfid_cnt);
      // not true, some edges do not have neighboring face with dup_cnt > 1
      // assert(nfid_cnt > 0);
    }
  }

  // check mmesh faces dup_cnt
  for (const auto& mface : this->faces) {
    if (mface.dup_cnt == 1) continue;
    num_dup_faces++;

    int neid_cnt = 0;
    for (const int neid : mface.edges_) {
      if (this->edges.at(neid).dup_cnt > 1) {
        assert(all_dup_edges.find(neid) != all_dup_edges.end());
        neid_cnt++;
      }
    }
    if (neid_cnt < 1) {
      printf("[validate_mmesh] mface %d (%d,%d,%d) has neid_cnt: %d\n",
             mface.fid, mface.vertices_[0], mface.vertices_[1],
             mface.vertices_[2], neid_cnt);
      assert(neid_cnt > 0);  // must be true?
    }
  }

  printf("[validate_mesh] num_dup_edges: %d, num_dup_faces %d\n", num_dup_edges,
         num_dup_faces);
  assert(num_dup_edges >= num_dup_faces);
}

// validate_mmesh_dup_cnt() will makes sure:
// each dup_cnt(face) contains > 1 dup_cnt(edge)
void MedialMesh::update_mmesh_dup_ef_map(bool is_debug) {
  if (is_debug) printf("calling update_mmesh_dup_ef_map....\n");
  std::map<int, std::set<int>> map_edge2faces;
  for (const auto& medge : this->edges) {
    if (medge.dup_cnt == 1) continue;
    // keep the key, but make value empty
    map_edge2faces[medge.eid] = std::set<int>();
    for (const int nfid : medge.faces_) {
      if (this->faces.at(nfid).dup_cnt > 1)
        map_edge2faces[medge.eid].insert(nfid);
    }
  }

  // dup_e2f map for every dup_cnt(edge) > 1
  this->dup_e2f.clear();
  std::set<int> untouch_fs;
  // dup edge has only 1 dup face
  // update untouch_fs
  for (const auto& e2fs : map_edge2faces) {
    if (e2fs.second.size() != 1) continue;
    int eid = e2fs.first;
    int fid = *e2fs.second.begin();
    this->dup_e2f[eid] = fid;
    this->dup_f2e[fid] = eid;
    untouch_fs.insert(fid);
  }
  // dup edge has > 1 dup faces
  for (const auto& e2fs : map_edge2faces) {
    if (e2fs.second.size() <= 1) continue;
    assert(e2fs.second.size() > 1);
    int eid = e2fs.first;
    bool is_ef_exist = false;
    for (const auto& fid : e2fs.second) {
      if (untouch_fs.find(fid) != untouch_fs.end()) continue;
      this->dup_e2f[eid] = fid;
      this->dup_f2e[fid] = eid;
      is_ef_exist = true;
      break;
    }
    assert(is_ef_exist);
  }
  // dup edge has 0 dup faces
  // assign eid -> -1, let function handle_dupcnt_face_edge_pair() to handle
  // on-the-fly, to avoid disconnectivity
  for (const auto& e2fs : map_edge2faces) {
    if (e2fs.second.size() > 0) continue;
    assert(e2fs.second.size() == 0);
    int eid = e2fs.first;
    this->dup_e2f[eid] = -1;
  }
  // make sure each dup edge has a face assigned to it
  assert(map_edge2faces.size() == this->dup_e2f.size());
  // make sure each dup face has a edge assigned to it
  int num_dup_face = 0;
  for (const auto& mface : this->faces) {
    if (mface.dup_cnt > 1) num_dup_face++;
  }
  assert(this->dup_f2e.size() == num_dup_face);
}

int MedialMesh::get_edge_fid_min_importance_keep_connectivity(int eid) {
  // find neighboring fid with least importance
  double min_import = DBL_MAX;
  int min_fid = -1;
  for (const auto& fid : this->edges.at(eid).faces_) {
    assert(!this->faces.at(fid).is_deleted);
    assert(this->faces.at(fid).dup_cnt == 1);

    // check if he face has any adj_edge that is a simple pair
    bool is_good = true;
    for (const auto& neid : this->faces.at(fid).edges_) {
      if (this->edges.at(neid).faces_.size() == 1) {
        // (fid, neid) is a simple pair
        // we cannot assign fid to eid anymore!!!
        // because deleting fid will make neid diconnect!!!!
        is_good = false;
        break;
      }
    }
    if (!is_good) continue;
    if (this->faces.at(fid).importance < min_import) {
      min_fid = fid;
      min_import = this->faces.at(fid).importance;
    }
  }  // all neigh fids
  assert(min_fid != -1);
  return min_fid;
}

int compute_Euler(const MedialMesh& mat) {
  // V-E+F (=2)
  // V-E+F-T (=1)
  std::set<int> connected_vertices, connected_edges, connected_faces;
  int num_faces = 0, num_edges = 0, num_tets = 0;
  for (const auto& face : mat.faces) {
    if (face.is_deleted) continue;
    num_faces += face.dup_cnt;
    connected_faces.insert(face.fid);
    std::vector<int> f_vs;
    for (const auto v : face.vertices_) {
      if (mat.vertices->at(v).is_deleted) {
        printf("mat face %d has deleted vertex %d\n", face.fid, v);
        // assert(false);
      }
      // assert(!mat.vertices->at(v).is_deleted);
      connected_vertices.insert(v);
      f_vs.push_back(v);
    }
    if (f_vs.size() != 3) {
      printf("mat face has < 3 vertices: %ld \n", f_vs.size());
      assert(false);
    }
    // If edge not exist, then face should not exist
    for (const int e : face.edges_) {
      assert(!mat.edges.at(e).is_deleted);
      connected_edges.insert(e);
    }
  }

  for (const auto& edge : mat.edges) {
    if (edge.is_deleted) continue;
    num_edges += edge.dup_cnt;

    // // only for edges with adjacent faces
    // if (edge.faces_.empty()) continue;
    // aint2 tmp = edge.vertices_;
    // std::sort(tmp.begin(), tmp.end());
    // connected_edges.insert(tmp);
    // connected_vertices.insert(tmp[0]);
    // connected_vertices.insert(tmp[1]);
  }

  for (const auto& tet : mat.tets) {
    if (tet.is_deleted) continue;
    num_tets += tet.dup_cnt;
  }

  int num_vs = 0;
  for (const auto& msphere : *mat.vertices) {
    if (msphere.is_deleted) continue;
    // if (msphere.pcell.cell_ids.empty()) continue;
    if (msphere.edges_.empty() && msphere.faces_.empty()) continue;
    num_vs++;
  }

  printf(
      "connected_vertices: %zu, connected_edges: %zu, connected_faces: "
      "%zu\n",
      connected_vertices.size(), connected_edges.size(),
      connected_faces.size());

  printf("num_vs %d, num_edges: %d, num_faces: %d, num_tets %d \n", num_vs,
         num_edges, num_faces, num_tets);

  printf(
      "numSpheres_active: %d, numEdges_active: %d, numFaces_active: %d, "
      "numTets_active %d\n",
      mat.numSpheres_active, mat.numEdges_active, mat.numFaces_active,
      mat.numTets_active);

  int euler = num_vs - num_edges + num_faces - num_tets;
  int euler2 = connected_vertices.size() - connected_edges.size() +
               connected_faces.size() - num_tets;
  int euler3 = mat.numSpheres_active - mat.numEdges_active +
               mat.numFaces_active - mat.numTets_active;
  printf("Euler Characteristic: %d, euler2: %d. euler3 %d \n", euler, euler2,
         euler3);
  return euler;
}
