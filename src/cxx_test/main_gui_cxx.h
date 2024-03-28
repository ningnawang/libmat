#ifndef H_MAIN_GUI_H
#define H_MAIN_GUI_H

#include <geogram/mesh/mesh.h>

#include <cmath>

#include "../medial_sphere.h"

class MainGuiWindow {
 private:
  static MainGuiWindow* instance_;
  static constexpr double point_radius_rel = 0.01;

 private:
  // a pointer
  GEO::Mesh* sf_mesh = nullptr;
  std::vector<MedialSphere>* all_medial_spheres = nullptr;
  //   std::vector<float>* tet_vertices = nullptr;
  //   std::vector<int>* tet_indices = nullptr;
  //   std::vector<std::array<float, 3>>* sf_vertices = nullptr;
  //   std::vector<std::array<int, 3>>* sf_faces = nullptr;

 public:
  void set_sf_mesh(GEO::Mesh& _sf_mesh);
  void set_all_medial_spheres(std::vector<MedialSphere>& _all_medial_spheres);
  //   void set_tet_mesh(std::vector<float>& _tet_vertices,
  //                     std::vector<int>& _tet_indices);
  //   void set_sf_mesh_internal(std::vector<std::array<float, 3>>&
  //   _sf_vertices,
  //                             std::vector<std::array<int, 3>>& _sf_fases);

 public:
  MainGuiWindow();
  ~MainGuiWindow();  // implemented in main_gui.cpp

  // show polyscope window
  void show();
  static void callbacks();

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
};

#endif  // __H_MAIN_GUI_H__
