#include "io_cuda.h"

#include "io.h"
#include "io_utils.hpp"

Vector3 compute_cell_barycenter(const ConvexCellHost& cc_trans) {
  cfloat3 bary = {0.0f, 0.0f, 0.0f};  // barycenter of cell
  FOR(i, cc_trans.nb_v) {
    cfloat4 voro_vertex = cc_trans.compute_vertex_coordinates(
        cmake_uchar3(cc_trans.ver_trans_const(i)));
    bary =
        cplus3(cmake_float3(voro_vertex.x, voro_vertex.y, voro_vertex.z), bary);
  }

  bary = cdivide3(bary, cc_trans.nb_v);
  return Vector3(bary.x, bary.y, bary.z);
}

// return 3d triangle mesh if is_triangle is true
void get_one_convex_cell_faces(
    const ConvexCellHost& cc_trans, std::vector<float3>& voro_points,
    std::vector<std::vector<unsigned>>& one_voro_cell_faces,
    std::vector<int>& voro_faces_sites, bool is_triangle, int max_sf_fid,
    bool is_boundary_only) {
  int row = voro_points.size();  // index from 0
  assert(one_voro_cell_faces.size() == voro_faces_sites.size());  // not clean!!

  // save all vertices
  FOR(i, cc_trans.nb_v) {
    cfloat4 voro_vertex = cc_trans.compute_vertex_coordinates(
        cmake_uchar3(cc_trans.ver_trans_const(i)));
    // ninwang: TODO check this
    if (std::isnan(voro_vertex.x) || std::isnan(voro_vertex.y) ||
        std::isnan(voro_vertex.z)) {
      printf(
          "ERROR: cell of site %d has vertex is NaN (%f, %f, %f), dual "
          "triangle (%d, %d, %d)\n",
          cc_trans.voro_id, voro_vertex.x, voro_vertex.y, voro_vertex.z,
          cc_trans.ver_trans_const(i).x, cc_trans.ver_trans_const(i).y,
          cc_trans.ver_trans_const(i).z);
      cc_trans.print_info();
      return;
    }
    voro_points.push_back(
        make_float3(voro_vertex.x, voro_vertex.y, voro_vertex.z));
  }

  // some clipping planes may not exist in tri but
  // we still sotre it, here is to filter those planes
  // if (!cc_trans.is_active_updated) cc_trans.reload_active();
  assert(cc_trans.is_active_updated);
  const std::vector<int>& active_clipping_planes =
      cc_trans.active_clipping_planes;

  FOR(plane, cc_trans.nb_p) {
    if (active_clipping_planes[plane] <= 0) continue;
    cint2 hp = cc_trans.clip_id2_const(plane);

    // check if only shown boundary or not
    // 1. halfplane; 2. on sf_mesh fid
    bool is_shown = true;
    if (is_boundary_only) {
      is_shown = false;
      if ((hp.y != -1) || (hp.y == -1 && hp.x < max_sf_fid)) {
        is_shown = true;
      }
    }
    if (!is_shown) continue;

    std::vector<int> tab_v;   // index of dual vertex
    std::vector<int> tab_lp;  // local index of dual vertex in dual triangle
    // for each dual triangle
    FOR(t, cc_trans.nb_v) {
      // store info of dual vertex
      if ((int)cc_trans.ver_trans_const(t).x == plane) {
        tab_v.push_back(t);
        tab_lp.push_back(0);
      } else if ((int)cc_trans.ver_trans_const(t).y == plane) {
        tab_v.push_back(t);
        tab_lp.push_back(1);
      } else if ((int)cc_trans.ver_trans_const(t).z == plane) {
        tab_v.push_back(t);
        tab_lp.push_back(2);
      }
    }

    if (tab_lp.size() <= 2) {
      std::cout << (int)plane << std::endl;
    }

    int i = 0;
    int j = 0;
    int lf = 0;  // local fid
    std::vector<std::vector<int>> voro_local_faces;
    voro_local_faces.push_back(std::vector<int>(0));

    while (voro_local_faces[lf].size() < tab_lp.size()) {
      int ind_i = (tab_lp[i] + 1) % 3;
      bool temp = false;
      j = 0;
      while (temp == false) {
        int ind_j = (tab_lp[j] + 2) % 3;
        if ((int)cc_trans.ith_plane_const(tab_v[i], ind_i) ==
            (int)cc_trans.ith_plane_const(tab_v[j], ind_j)) {
          voro_local_faces[lf].push_back(tab_v[i]);
          temp = true;
          i = j;
        }
        j++;
      }
    }

    int nb_pts = voro_local_faces[lf].size();
    std::vector<unsigned> one_face;
    if (is_triangle) {
      // triangulate face, then make a tet
      for (uint p = 1; p < nb_pts - 1; p++) {
        // (0, p, p+1) is a triangulation of face
        one_face.clear();
        one_face.push_back(row + voro_local_faces[lf][0]);
        one_face.push_back(row + voro_local_faces[lf][p]);
        one_face.push_back(row + voro_local_faces[lf][(p + 1) % nb_pts]);
        one_voro_cell_faces.push_back(one_face);
      }
    } else {
      one_face.clear();
      for (uint p = 0; p < nb_pts; p++) {
        one_face.push_back(row + voro_local_faces[lf][p]);
      }
      one_voro_cell_faces.push_back(one_face);
    }

    voro_faces_sites.push_back(cc_trans.voro_id);
    lf++;
  }
}

// shibo: houdini output test (file .geo)
// TODO: customized file name
bool save_convex_cells_houdini(
    const Parameter params, const std::vector<MedialSphere>& all_medial_spheres,
    const std::vector<ConvexCellHost>& convex_cells_returned,
    std::string rpd_name, const int max_sf_fid, const bool is_boundary_only,
    bool is_slice_plane) {
  IO::Geometry geometry;
  IO::GeometryWriter geometry_writer(".");
  std::string rpd_path = "../out/rpd/rpd_" + rpd_name + "_" + get_timestamp();

  std::vector<float3> voro_points;
  std::vector<std::vector<unsigned>> all_voro_faces;
  std::vector<int> voro_faces_sites;

  for (const auto& cc_trans : convex_cells_returned) {
    // filter spheres by slicing plane
    Vector3 last_bary = compute_cell_barycenter(cc_trans);
    if (is_slice_plane && is_slice_by_plane(last_bary, params)) continue;
    get_one_convex_cell_faces(cc_trans, voro_points, all_voro_faces,
                              voro_faces_sites, false /*is_triangle*/,
                              max_sf_fid, is_boundary_only);
  }
  geometry.AddParticleAttribute("P", voro_points);
  geometry.AddPolygon(all_voro_faces);
  geometry.AddPrimitiveAttribute("PrimAttr", voro_faces_sites);
  geometry_writer.OutputGeometry(rpd_path, geometry);

  printf("saved houdini file: %s \n", rpd_path.c_str());
  return true;
}