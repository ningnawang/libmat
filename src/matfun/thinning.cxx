#include "thinning.h"

// helper funtion to load mface's importance
void Thinning::compute_face_importance(const MedialMesh& mat, MedialFace& face,
                                       bool is_sort_randomly, bool is_debug) {
  if (is_sort_randomly) {
    face.importance = (double)std::rand() / (RAND_MAX);
  } else {
    // do not count feature spheres
    // they are making avg_radius smaller
    double avg_radius = 0;
    int nb_extf = 0;
    for (const auto& vid : face.vertices_) {
      auto& vertex = mat.vertices->at(vid);
      if (vertex.is_on_extf()) nb_extf++;
      avg_radius += vertex.radius;
    }
    int cnt = face.vertices_.size() - nb_extf;
    avg_radius /= std::max(cnt, 1);
    if (nb_extf > 2 || avg_radius <= SCALAR_ZERO_2) {
      // three vertices could all be feature vertices
      // this happends mostly for non-feature models
      face.importance = 0.01;
      if (is_debug)
        printf("[Prune] fid %d has nb_extf: %d, avg_radius: %f\n", face.fid,
               nb_extf, avg_radius);
    } else {
      double dual_edge_dist = 0;
      for (int i = 0; i < face.dual_edge_endpoints.size(); i++) {
        dual_edge_dist += (face.dual_edge_endpoints[i][0].first -
                           face.dual_edge_endpoints[i][1].first)
                              .length();
      }
      face.importance = dual_edge_dist / (avg_radius * 2);
    }
  }
}

// just store importance for each mat face
void Thinning::load_all_mat_face_importance_globally(MedialMesh& mat,
                                                     bool is_sort_randomly,
                                                     bool is_debug) {
  for (auto& face : mat.faces) {
    compute_face_importance(mat, face, is_sort_randomly, is_debug);
  }
  printf("Saved mat faces importances\n");
}

// priority queue is ordered by first element of the pair
void Thinning::sort_mat_face_importance_globally(
    MedialMesh& mat, std::set<Thinning::face_importance>& imp_queue,
    bool is_import_given, bool is_sort_randomly, bool is_debug) {
  imp_queue.clear();
  for (const auto& one_tet : mat.tets) {
    int tet_id = one_tet.tid;
    for (const auto& fid : one_tet.faces_) {
      auto& face = mat.faces[fid];
      if (face.tets_.empty()) {
        // face.importance = DBL_MAX;
        printf("face of tet cannot have empty tets_!!!\n");
        assert(false);
        continue;
      }

      if (!is_import_given)
        compute_face_importance(mat, face, is_sort_randomly, is_debug);
      // store to priority queue
      imp_queue.insert(std::pair<double, int>(face.importance, fid));
      if (is_debug)
        printf("[Prune] fid %d has importance: %f, pushed to priority queue \n",
               fid, face.importance);

    }  // for mat face
  }  // for mat tets
  if (is_debug)
    printf("[Prune] total %ld mat faces in priority queue \n",
           imp_queue.size());
}

void prune_one_simple_pair(const simple_pair& sp_del, MedialMesh& mat,
                           bool is_debug) {
  if (is_debug) {
    printf("[Prune] prune simple pair: \n");
    sp_del.print_info();
  }
  if (sp_del.type == 2) {  // tet-face pair
    int tet_id = sp_del.idx0;
    int fid = sp_del.idx1;
    mat.delete_tet(tet_id);
    mat.delete_face(fid);
    if (is_debug) {
      auto& face = mat.faces[fid];
      // face.delete_for_tet = tet_id;
      printf("[Prune] deleted tet %d and face %d\n", tet_id, fid);
      printf("[Prune] face %d have tets: %ld \n", fid, face.tets_.size());
    }
  } else if (sp_del.type == 1) {  // face-edge pair
    int fid = sp_del.idx0;
    int eid = sp_del.idx1;
    mat.delete_face(fid);
    mat.delete_edge(eid);
    if (is_debug) printf("[Prune] deleted face %d and edge %d \n", fid, eid);
  } else if (sp_del.type == 0) {  // edge-vertex pair
    int eid = sp_del.idx0;
    int vid = sp_del.idx1;
    mat.delete_edge(eid);
    mat.delete_vertex(vid);
    if (is_debug) printf("[Prune] deleted edge %d and vertex %d \n", eid, vid);
  } else {
    // error
    printf("[Prune] error, unkown type %d to prune \n", sp_del.type);
    assert(false);
  }
}

// face with dup_cnt>1 must have an adjacent edge with dup_cnt>1
void postprocess_dupcnt_tet_face(MedialMesh& mat, int tid, int fid,
                                 bool is_debug) {
  auto& mtet = mat.tets.at(tid);
  assert(mtet.dup_cnt == 1);
  auto& mface = mat.faces.at(fid);
  assert(mface.is_deleted);
  // case 1: tet-face-edge dup_cnt is 1-1-2: no nothing
  if (mface.dup_cnt < 2) return;  // do nothing
  // assert(mface.dup_cnt == 2);

  if (is_debug)
    printf("[DUP_CNT] tet %d has fid %d with dup_cnt %d\n", tid, fid,
           mface.dup_cnt);

  // case 2: tet-face-edge dup_cnt is 1-x-x, x>1
  // find adjacent edge with dup_cnt>1
  // given in mat.dup_f2e, updated through func update_mmesh_dup_ef_map()
  int neid = mat.dup_f2e.at(fid);
  assert(neid != -1);
  if (is_debug)
    printf(
        "[DUP_CNT] tet %d has fid %d dup(fid): %d, updating eid %d with "
        "dup(eid): %d\n",
        tid, fid, mface.dup_cnt, neid, mface.dup_cnt);

  // clear edge's dup_cnt
  auto& mat_e = mat.edges.at(neid);
  assert(mat_e.dup_cnt == mface.dup_cnt);  // 2-2 or 3-3 are all possible
  int diff = mat_e.dup_cnt - 1;
  mat_e.dup_cnt -= diff;
  mat.numEdges_active -= diff;
  assert(!mat_e.faces_.empty());
}

void Thinning::prune_tets_while_iteration(
    MedialMesh& mat, std::set<Thinning::face_importance>& imp_queue,
    bool is_dup_cnt, bool is_debug) {
  if (mat.tets.empty()) return;
  if (imp_queue.empty()) {
    printf("[Prune] imp_queue is not empty?? imp_queue %ld \n",
           imp_queue.size());
    return;
  }

  while (true) {
    if (is_debug) printf("still %d mat tets to prune \n", mat.numTets_active);
    for (const auto& imp_pair : imp_queue) {
      double f_imp = imp_pair.first;
      int fid = imp_pair.second;
      auto& face = mat.faces[fid];
      if (face.is_deleted) continue;
      if (face.tets_.empty() || face.tets_.size() > 1) continue;
      if (face.tets_.size() != 1) {
        printf("ERROR face has more than 1 tets \n");
        assert(false);
      }

      // everytime we delete a face, we break for loop of imp_queue
      // and checking again from face with smallest importance
      int tid = *(face.tets_.begin());
      // push <tet, face> to delete
      auto sp_del = simple_pair(2, tid, fid, tid);
      if (is_debug) {
        printf(
            "[Prune] add tet_face pair to prune: tid: %d, fid: %d, "
            "importance: %f\n",
            tid, fid, f_imp);
      }
      prune_one_simple_pair(sp_del, mat, is_debug);
      // taking care of dup_cnt
      if (is_dup_cnt)
        postprocess_dupcnt_tet_face(mat, tid, fid, is_debug /*is_debug*/);
      break;
    }  // for imp_queue

    // break while loop once we have all tets cleaned
    if (mat.numTets_active == 0) break;
  }  // while true
}

void Thinning::handle_dupcnt_face_edge_pair(MedialMesh& mat, bool is_debug) {
  for (const auto& ef_pair : mat.dup_e2f) {
    int eid = ef_pair.first;
    int fid = ef_pair.second;
    if (fid == -1) {  // given in MedialMesh::update_mmesh_dup_ef_map()
      // assign new fid on-the-fly, to avoid diconnectivity problem
      fid = mat.get_edge_fid_min_importance_keep_connectivity(eid);
    }

    auto& mface = mat.faces.at(fid);
    auto& medge = mat.edges.at(eid);
    if (mface.is_deleted && medge.is_deleted) continue;
    // handled by postprocess_dupcnt_tet_face()
    if (mface.is_deleted && mface.dup_cnt > 1 && !medge.is_deleted &&
        medge.dup_cnt == 1)
      continue;

    // check other unkown cases
    if ((mface.is_deleted && !medge.is_deleted) ||
        (!mface.is_deleted && medge.is_deleted)) {
      printf("[DUP_CNT] handling face %d dup_cnt %d, is_deleted %d\n", fid,
             mface.dup_cnt, mface.is_deleted);
      printf("[DUP_CNT] handling edge %d dup_cnt %d, is_deleted %d\n", eid,
             medge.dup_cnt, medge.is_deleted);
      assert(false);
    }

    if (is_debug) {
      printf("[DUP_CNT] handling face %d dup_cnt %d, is_deleted %d\n", fid,
             mface.dup_cnt, mface.is_deleted);
      printf("[DUP_CNT] handling edge %d dup_cnt %d, is_deleted %d\n", eid,
             medge.dup_cnt, medge.is_deleted);
    }
    assert(((medge.dup_cnt > 1 && mface.dup_cnt > 1) &&
            medge.dup_cnt == mface.dup_cnt) ||
           (medge.dup_cnt == 2 && mface.dup_cnt == 1));

    if ((medge.dup_cnt > 1 && mface.dup_cnt > 1) &&
        medge.dup_cnt == mface.dup_cnt) {
      // clear face's dup_cnt
      mat.faces.at(fid).dup_cnt -= (mface.dup_cnt - 1);
      mat.numFaces_active -= (mface.dup_cnt - 1);
      // clear edge's dup_cnt
      mat.edges.at(eid).dup_cnt -= (medge.dup_cnt - 1);
      mat.numEdges_active -= (medge.dup_cnt - 1);
    } else {  // case: medge.dup_cnt == 2 && mface.dup_cnt == 1
      // sanity check
      for (const auto& neid : mat.faces.at(fid).edges_) {
        if (mat.edges.at(neid).faces_.size() == 1) {
          printf(
              "-----ERROR: deleting face %d (%d,%d,%d,) edge %d (%d,%d) will "
              "make edge %d (%d,%d) disconnect!\n",
              fid, mat.faces.at(fid).vertices_[0],
              mat.faces.at(fid).vertices_[1], mat.faces.at(fid).vertices_[2],
              eid, mat.edges.at(eid).vertices_[0],
              mat.edges.at(eid).vertices_[1], neid,
              mat.edges.at(neid).vertices_[0], mat.edges.at(neid).vertices_[1]);
          assert(false);
        }
      }
      mat.delete_face(fid);
      // clear edge's dup_cnt
      mat.edges.at(eid).dup_cnt -= 1;
      mat.numEdges_active -= 1;
      assert(mat.edges.at(eid).faces_.size() > 0);
    }
  }
}

void Thinning::prune_faces_while_iteration(
    const std::vector<MedialSphere>& all_medial_spheres, MedialMesh& mat,
    std::set<Thinning::face_importance>& imp_queue, double imp_thres,
    bool is_debug) {
  if (imp_queue.empty()) {
    printf("[Prune] error, imp_queue is not empty?? imp_queue %ld\n",
           imp_queue.size());
    return;
  }
  auto is_edge_on_ext_feature = [&](const int eid) {
    const auto& medge = mat.edges.at(eid);
    if (medge.is_extf) return true;
    const auto& mv0 = mat.vertices->at(medge.vertices_.at(0));
    const auto& mv1 = mat.vertices->at(medge.vertices_.at(1));
    if (is_two_mspheres_on_same_sl_including_corners(mv0, mv1)) return true;
    return false;
  };

  auto is_skip_edge_if_on_boundary = [&](const int eid) {
    const auto& medge = mat.edges.at(eid);
    if (medge.is_deleted) return true;
    const auto& mv0 = mat.vertices->at(medge.vertices_.at(0));
    const auto& mv1 = mat.vertices->at(medge.vertices_.at(1));
    if (medge.is_on_same_sheet && (mv0.is_on_corner() || mv1.is_on_corner()))
      return true;
    return false;
  };

  int num_tet_faces_are_good = 0;
  while (true) {
    if (is_debug)
      printf("num_tet_faces_are_good: %d, imp_queue.size %ld \n",
             num_tet_faces_are_good, imp_queue.size());
    num_tet_faces_are_good = 0;
    for (const auto& imp_pair : imp_queue) {
      bool is_break_queue_loop = false;
      double f_imp = imp_pair.first;
      int fid = imp_pair.second;
      auto& face = mat.faces[fid];
      if (face.importance != f_imp) {
        printf("ERROR: face %d cannot have inconsistant importace: %f and %f\n",
               fid, face.importance, f_imp);
        assert(false);
      }
      // if (face.is_deleted) {
      if (face.is_deleted || face.importance >= imp_thres) {
        num_tet_faces_are_good++;
        continue;
      }
      if (!face.tets_.empty()) {
        printf("Tets of mat face cannot be empty\n");
        assert(false);
      }
      // everytime we delete a face, we break for loop of imp_queue
      // and checking again from face with smallest importance
      for (const auto& eid : face.edges_) {
        const auto edge = mat.edges[eid];
        if (edge.faces_.size() > 1) continue;  // check next edge

        // if edge is on sharp edge, then skip
        if (is_edge_on_ext_feature(eid)) continue;

        // if edge is on corner, and single face, then skip
        if (is_skip_edge_if_on_boundary(eid)) continue;

        // push <face, edge> to delete
        auto sp_del = simple_pair(1, fid, eid, -1);
        if (is_debug) {
          printf(
              "[Prune] add face_edge pair to prune: fid %d: (%d,%d,%d), eid "
              "%d: (%d,%d), importance: %f \n",
              fid, face.vertices_[0], face.vertices_[1], face.vertices_[2], eid,
              edge.vertices_[0], edge.vertices_[1], f_imp);
        }
        prune_one_simple_pair(sp_del, mat, is_debug);
        is_break_queue_loop = true;
        break;
      }  // for face.edges_
      if (is_break_queue_loop) break;
      // the face is good, every edge has 2 neighboring faces
      num_tet_faces_are_good++;
    }  // for imp_queue

    // break while loop once we have all tets faces cleaned
    if (num_tet_faces_are_good == imp_queue.size()) break;
  }  // while true

  return;
}

void Thinning::prune_edges_while_iteration(
    const std::vector<MedialSphere>& all_medial_spheres, MedialMesh& mat,
    bool is_debug) {
  int num_spheres_are_good;
  while (true) {
    num_spheres_are_good = 0;
    for (const auto& mvertex : all_medial_spheres) {
      if (mvertex.is_deleted) {
        num_spheres_are_good++;
        continue;
      }
      if (mvertex.edges_.size() > 1) {
        num_spheres_are_good++;
        continue;
      }
      int eid = *mvertex.edges_.begin();
      assert(mat.edges.at(eid).faces_.empty());
      // push <edge, vertex> to delete
      auto sp_del = simple_pair(0, eid, mvertex.id, -1);
      prune_one_simple_pair(sp_del, mat, is_debug);
    }
    if (num_spheres_are_good == all_medial_spheres.size()) break;
  }
}

void Thinning::assert_mmesh(MedialMesh& mat) {
  for (const auto& mv : *mat.vertices) {
    if (mv.is_deleted) continue;
    // assert(mv.dup_cnt == 1);
    assert(!mv.edges_.empty());
    if (mv.faces_.empty()) {
      printf("======= ERROR: msphere %d has edges [", mv.id);
      for (const auto& eid : mv.edges_) {
        printf("%d, ", eid);
      }
      printf("], but empty face\n");
    }
    // assert(!mv.faces_.empty());
  }
  for (const auto& medge : mat.edges) {
    if (medge.is_deleted) continue;
    assert(medge.dup_cnt == 1);
    // if (medge.faces_.empty())
    //   printf("===== ERROR: medge %d (%d,%d) has empty adjacent face\n",
    //          medge.eid, medge.vertices_[0], medge.vertices_[1]);
    // assert(!medge.faces_.empty());
  }
  for (const auto& mface : mat.faces) {
    if (mface.is_deleted) continue;
    assert(mface.dup_cnt == 1);
  }
}

// imp_thres: more faces will be prunes if bigger
void Thinning::prune(const std::vector<MedialSphere>& all_medial_spheres,
                     MedialMesh& mat, double imp_thres, bool is_dup_cnt,
                     bool is_import_given, bool is_sort_randomly,
                     bool is_debug) {
  printf("start thinning with importance threshold: %f is_debug: %d... \n",
         imp_thres, is_debug);
  // using set as priority queue
  // (importance, fid) starts with the smallest importance
  // tie-broken with smallest fid.
  std::set<Thinning::face_importance> imp_queue;  // in ascending order
  sort_mat_face_importance_globally(mat, imp_queue, is_import_given,
                                    is_sort_randomly, is_debug);

  if (is_dup_cnt) {
    mat.update_mmesh_dup_ef_map(is_debug);
    prune_tets_while_iteration(mat, imp_queue, is_dup_cnt, is_debug);
    // printf("done tet-face prune...\n");
    handle_dupcnt_face_edge_pair(mat, is_debug);
    // printf("done dupcnt_face_edge ...\n");
    prune_faces_while_iteration(all_medial_spheres, mat, imp_queue, imp_thres,
                                is_debug);
    prune_edges_while_iteration(all_medial_spheres, mat, is_debug);
    assert_mmesh(mat);
  } else {
    prune_tets_while_iteration(mat, imp_queue, is_dup_cnt, is_debug);
    prune_faces_while_iteration(all_medial_spheres, mat, imp_queue, imp_thres,
                                is_debug);
  }
}