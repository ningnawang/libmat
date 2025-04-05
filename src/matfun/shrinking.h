#ifndef H_SHRINKING_H
#define H_SHRINKING_H

#include <assert.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "common_geogram.h"
#include "input_types.h"
#include "medial_sphere.h"
#include "params.h"

void pre_and_init_aabb(GEO::Mesh &sf_mesh, AABBWrapper &aabb_wrapper,
                       bool is_reorder);
void pre_and_init_feature_aabb(const TetMesh &tet_mesh,
                               AABBWrapper &aabb_wrapper);

// will update sf_mesh
void init_and_shrink(const SurfaceMesh &sf_mesh, const TetMesh &tet_mesh,
                     std::vector<MedialSphere> &all_medial_spheres,
                     int num_init_spheres, int itr_limit = -1,
                     bool is_debug = false);
void init_and_shrink_given_pins(const SurfaceMesh &sf_mesh,
                                const TetMesh &tet_mesh,
                                std::vector<MedialSphere> &all_medial_spheres,
                                const std::vector<v2int> &surface_pins,
                                int itr_limit = -1, bool is_debug = false);

// shrink one sphere
bool shrink_sphere(const SurfaceMesh &sf_mesh, const AABBWrapper &aabb_wrapper,
                   const std::set<aint2> &fe_sf_fs_pairs,
                   const std::vector<FeatureEdge> &feature_edges,
                   MedialSphere &msphere, int itr_limit, bool is_del_near_ce,
                   bool is_del_near_se, bool is_debug);

// will update TetMesh::fl2corner_sphere
int init_corner_spheres(const int num_itr_global, TetMesh &tet_mesh,
                        std::vector<MedialSphere> &all_medial_spheres);

void insert_spheres_for_concave_lines_new(
    const SurfaceMesh &sf_mesh, const std::vector<ConcaveCorner> &cc_corners,
    const std::vector<FeatureEdge> &feature_edges,
    std::vector<FeatureLine> &ce_lines,
    std::vector<MedialSphere> &all_medial_spheres,
    const double cc_len_eps /*=length, scaled in [0, Parameter::scale_max]*/,
    bool is_debug);

bool update_new_concave_sphere(const SurfaceMesh &sf_mesh,
                               const std::vector<FeatureEdge> &feature_edges,
                               const Vector3 &pin_point, const int fe_id,
                               const int sphere_type, MedialSphere &new_sphere,
                               bool is_debug);

bool update_msphere_given_v2fid(const SurfaceMesh &sf_mesh,
                                const TetMesh &tet_mesh,
                                const v2int v2fid_chosen,
                                const int new_sphere_id,
                                MedialSphere &new_msphere,
                                const bool is_merge_to_ce, bool is_debug);

bool insert_new_sphere_given_v2fid(
    const int num_itr_global, const SurfaceMesh &sf_mesh,
    const TetMesh &tet_mesh, const v2int v2fid_chosen,
    std::vector<MedialSphere> &all_medial_spheres, const bool is_merge_to_ce,
    bool is_debug);

#endif  // __SHRINKING_H__
