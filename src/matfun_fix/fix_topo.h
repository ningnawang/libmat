#ifndef H_TOPO_H
#define H_TOPO_H

#include <assert.h>
// #include <cuda_runtime.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "common_geogram.h"
#include "fix_common.h"
#include "medial_sphere.h"
#include "params.h"
#include "shrinking.h"
#include "voronoi_defs.h"

void check_cc_and_euler(std::vector<ConvexCellHost>& convex_cells_host,
                        std::vector<MedialSphere>& all_medial_spheres,
                        std::set<int>& spheres_to_fix, bool is_debug);
int fix_topo_by_adding_new_sphere(
    const int num_itr_global, const int num_topo_itr,
    const std::vector<ConvexCellHost>& convex_cells_host,
    const SurfaceMesh& sf_mesh, const TetMesh& tet_mesh,
    const std::set<int>& spheres_to_fix,
    std::vector<MedialSphere>& all_medial_spheres, bool is_debug);

v2int get_v2fid_max_to_point(const GEO::Mesh& sf_mesh,
                             const std::vector<v2int> surf_v2fids,
                             const Vector3& point);

v2int get_v2fid_min_to_point(const GEO::Mesh& sf_mesh,
                             const std::vector<v2int> surf_v2fids,
                             const Vector3& point);

void get_facet_CC_surf_v2fids(
    const std::map<int, std::vector<v2int>>& cell_to_surfv2fid,
    const std::map<int, std::vector<std::set<int>>>& facet_cc_cells,
    std::map<int, std::vector<std::vector<v2int>>>& facet_cc_surf_v2fids);

#endif  // __TOPO_H__
