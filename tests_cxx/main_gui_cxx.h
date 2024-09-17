#ifndef H_MAIN_GUI_H
#define H_MAIN_GUI_H

#include <geogram/mesh/mesh.h>

#include <cmath>

#include "../src/matbase/medial_mesh.h"
#include "../src/matbase/medial_sphere.h"

class MainGuiWindow {
 private:
  static MainGuiWindow* instance_;
  static constexpr double point_radius_rel = 0.01;
  static constexpr double edge_radius_rel = 0.0008;
  using timer = std::chrono::system_clock;

 private:
  // a pointer
  std::string file_name = "";
  TetMesh* tet_mesh = nullptr;
  SurfaceMesh* sf_mesh = nullptr;
  std::vector<MedialSphere>* all_medial_spheres = nullptr;
  Parameter* params = nullptr;  // not const, update site_k
  MedialMesh* mmesh = nullptr;  // not const

 public:
  void set_file_name(const std::string& _file_name);
  void set_params(Parameter& _params);
  void set_tet_mesh(TetMesh& _tet_mesh);
  void set_sf_mesh(SurfaceMesh& _sf_mesh);
  void set_all_medial_spheres(std::vector<MedialSphere>& _all_medial_spheres);
  void set_medial_mesh(MedialMesh& _mmesh);

 public:
  MainGuiWindow();
  ~MainGuiWindow();  // implemented in main_gui.cpp

  // show polyscope window
  void show();
  static void callbacks();

  // Tetmesh
  void show_tet_mesh(const TetMesh& tet_mesh, bool is_shown_meta = false);

  // Surface
  int given_sf_face_id = -1;
  void show_surface_mesh(const GEO::Mesh& sf_mesh, const int& given_sf_face_id);

  int given_sphere_id = -1;
  void show_pin_points(const std::vector<MedialSphere>& all_medial_spheres,
                       const int given_sphere_id);
  void show_all_medial_spheres(
      const std::vector<MedialSphere>& all_medial_spheres);
  void show_one_sphere(const std::vector<MedialSphere>& all_medial_spheres,
                       const int given_sphere_id,
                       const bool is_show_radius = false,
                       const bool is_clear_all = false);

  // Medial Mesh
  int given_mat_face_id = -1;
  void show_medial_mesh(const MedialMesh& mmesh, int given_mat_face_id = -1,
                        std::string mmname = "MedialMesh");
  double given_thinning_thres = 0.3;
  void show_medial_edges(const MedialMesh& mmesh);
};

#endif  // __H_MAIN_GUI_H__
