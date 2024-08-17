#pragma once

#include <geogram/mesh/mesh.h>

#include <chrono>
#include <cmath>

#include "input_types.h"
#include "medial_mesh.h"
#include "medial_sphere.h"
#include "rpd_api.h"
#include "triangulation.h"
#include "voronoi_defs.h"
class MainGuiWindow {
 private:
  static MainGuiWindow* instance_;
  static constexpr double point_radius_rel = 0.002;
  static constexpr double edge_radius_rel = 0.0008;
  using timer = std::chrono::system_clock;

 private:
  // pointers
  Parameter* params = nullptr;  // not const, update site_k
  std::vector<MedialSphere>* all_medial_spheres = nullptr;  // not const
  MedialMesh* mmesh = nullptr;                              // not const

 public:
  void set_params(Parameter& _params);
  // void set_tet_mesh(TetMesh& _tet_mesh);
  // void set_sf_mesh(const SurfaceMesh& _sf_mesh);
  void set_all_medial_spheres(std::vector<MedialSphere>& _all_medial_spheres);
  void set_medial_mesh(MedialMesh& _mmesh);

 public:
  MainGuiWindow();
  ~MainGuiWindow();  // implemented in main_gui.cpp
  static void callbacks();
  void show();
};