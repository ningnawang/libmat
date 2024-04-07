#ifndef H_IO_H
#define H_IO_H

#include <assert.h>
#include <geogram/mesh/mesh.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "../src/external/Predicates.hpp"
#include "common_geogram.h"
#include "medial_mesh.h"
#include "medial_sphere.h"
#include "params.h"
#include "voronoi_defs.h"

bool is_inverted(const Vector3& v0, const Vector3& v1, const Vector3& v2,
                 const Vector3& v3);

void get_bbox(const std::vector<float>& vertices, float& xmin, float& ymin,
              float& zmin, float& xmax, float& ymax, float& zmax,
              float& bbox_diag_l);

bool load_tet(const std::string& filename, std::vector<float>& vertices,
              std::vector<int>& indices, bool normalize, Parameter& params);

void load_tet_adj_info(const std::map<int, std::set<int>>& v2tets,
                       const std::vector<int>& tet_indices,
                       const std::map<int, std::set<int>>& tet_vs2sf_fids,
                       const int n_sf_facets, std::vector<int>& v_adjs,
                       std::vector<int>& e_adjs, std::vector<int>& f_adjs,
                       std::vector<int>& f_ids);
void load_spheres_to_sites_given(std::vector<MedialSphere>& all_medial_spheres,
                                 bool& site_is_transposed,
                                 std::vector<float>& site,
                                 std::vector<float>& site_weights, int& n_site);

void load_spheres_from_file(const char* filename,
                            std::vector<MedialSphere>& all_medial_spheres,
                            bool is_load_type, bool is_load_deleted);

void save_spheres_file(const std::vector<MedialSphere>& all_medial_spheres,
                       const std::string filename, bool is_save_type,
                       bool is_load_deleted, int num_rpd_itr = -1);

void load_v2tets(const std::vector<float>& vertices,
                 const std::vector<int>& indices,
                 std::map<int, std::set<int>>& v2tets);

void load_surface_vertices(const std::vector<float>& vertices,
                           const std::vector<int>& indices,
                           std::vector<float>& surf_vertices);

bool load_surface_mesh(const std::string& path, GEO::Mesh& input);
bool load_surface_mesh_geogram(const std::string& path, GEO::Mesh& input);

void load_sf_mesh_from_internal(const std::vector<std::array<float, 3>>& points,
                                const std::vector<std::array<int, 3>>& faces,
                                const std::vector<int>& sf2tet_vs_mapping,
                                GEO::Mesh& input);

// void write_convex_cells(std::vector<float3>& voro_points,
//                         std::vector<std::array<int, 4>>& voro_tets,
//                         std::vector<int>& voro_tets_sites);

void get_surface_from_tet(const std::vector<float>& tet_vertices,
                          const std::vector<int>& tet_indices,
                          std::vector<std::array<float, 3>>& surf_vertices,
                          std::vector<std::array<int, 3>>& surf_faces,
                          std::vector<int>& sf_vs_2_tet_vs);

bool save_sf_mesh(const std::string sf_path, const GEO::Mesh& sf_mesh);
bool save_sf_mesh_with_extf(const std::string sf_path,
                            const GEO::Mesh& sf_mesh);
bool save_sf_mesh_scaled(const std::string sf_path_scaled,
                         const GEO::Mesh& sf_mesh, const Parameter& param);
bool save_sf_mesh_geogram(const std::string sf_path, GEO::Mesh& sf_mesh);

bool is_slice_by_plane(const Vector3& bary, const Parameter& params);
Vector3 compute_cell_barycenter(const ConvexCellHost& cc_trans);

bool save_convex_cells_houdini(
    const Parameter params, const std::vector<MedialSphere>& all_medial_spheres,
    const std::vector<ConvexCellHost>& convex_cells_returned,
    std::string rpd_name, const int max_sf_fid, const bool is_boundary_only,
    bool is_slice_plane = false);

// load ma
void unnormalize_matfp(const Parameter& params, MedialMesh& mmesh_matfp);
void renormalize_matfp(const Parameter& params, MedialMesh& mmesh_matfp);
void load_matfp(const std::string& ma_path,
                std::vector<MedialSphere>& all_medial_spheres, MedialMesh& mat);

// save ma
void export_ma(const std::string& maname, const MedialMesh& mat);
void write_ma_ply(const std::string& maname, const MedialMesh& mat);

#endif  // __IO_H__
