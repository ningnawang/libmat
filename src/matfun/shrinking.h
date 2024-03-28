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

void pre_and_init_aabb(GEO::Mesh &sf_mesh, AABBWrapper &aabb_wrapper);
void pre_and_init_feature_aabb(const TetMesh &tet_mesh,
                               AABBWrapper &aabb_wrapper);

// will update sf_mesh
void init_and_shrink(const SurfaceMesh &sf_mesh, const TetMesh &tet_mesh,
                     std::vector<MedialSphere> &all_medial_spheres,
                     int num_init_spheres, int itr_limit = -1,
                     bool is_debug = false);

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

int create_new_concave_sphere_given_pin(
    const SurfaceMesh &sf_mesh, const std::vector<FeatureEdge> &feature_edges,
    const Vector3 &pin_point, const int fe_id,
    std::vector<MedialSphere> &all_medial_spheres, int sphere_type,
    bool is_debug = false);

// create CC spheres in a while loop, try 10 times
// 1. init concave spheres
// 2. topo fix
int create_new_concave_sphere_given_pin_wrapper(
    const SurfaceMesh &sf_mesh, const std::vector<FeatureEdge> &feature_edges,
    const Vector3 &pin_point, const int fe_id,
    std::vector<MedialSphere> &all_medial_spheres, int sphere_type,
    bool is_debug);

#endif  // __SHRINKING_H__
