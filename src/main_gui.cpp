#include "main_gui.h"

#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/volume_mesh.h>

#include "fix_extf.h"
#include "fix_geo.h"
#include "fix_geo_error.h"
#include "fix_intf.h"

MainGuiWindow* MainGuiWindow::instance_ = nullptr;
MainGuiWindow::MainGuiWindow() {
  if (instance_ != nullptr) {
    printf("ERROR: GuiWindow instance is not nullptr!!");
    exit(1);
  }
  instance_ = this;
}

MainGuiWindow::~MainGuiWindow() {
  // we delete here since we allocate memory
  // (called new) for each of these variables
  delete params;
  //   delete tet_mesh;
  //   delete sf_mesh;
  delete all_medial_spheres;
  delete mmesh;
  //   delete mmesh_gt;
  //   delete spheres_to_fix;
  //   delete rt;
  //   delete rpd3d;
  //   delete opt_rpd;

  this->instance_ = nullptr;
}

void MainGuiWindow::set_params(Parameter& _params) { params = &_params; }
void MainGuiWindow::set_all_medial_spheres(
    std::vector<MedialSphere>& _all_medial_spheres) {
  all_medial_spheres = &_all_medial_spheres;
}
void MainGuiWindow::set_medial_mesh(MedialMesh& _mmesh) { mmesh = &_mmesh; }

void MainGuiWindow::show() {
  // a few camera options
  polyscope::view::upDir = polyscope::UpDir::ZUp;

  // Initialize Polyscope
  polyscope::init();

  // Add the callback
  polyscope::state::userCallback = MainGuiWindow::callbacks;

  // Show the GUI
  polyscope::show();
}

void MainGuiWindow::callbacks() {
  ImGui::PushItemWidth(100);

  // // Surface
  // if (ImGui::Button("show surface mesh")) {
  //   instance_->show_surface_mesh(*(instance_->sf_mesh),
  //                                -1 /*given_sf_face_id*/);
  // }
  // ImGui::SameLine();
  if (ImGui::Button("show pins")) {
    instance_->show_pin_points(*(instance_->all_medial_spheres),
                               -1 /*given_sphere_id*/);
  }

  ImGui::InputInt("#face_id", &(instance_->given_sf_face_id));
  ImGui::SameLine();
  if (ImGui::Button("show")) {
    instance_->show_surface_mesh(*(instance_->sf_mesh),
                                 instance_->given_sf_face_id);
  }

  // Medial spheres
  if (ImGui::SmallButton("show all mat centers")) {
    instance_->show_all_medial_spheres(*(instance_->all_medial_spheres));
  }

  ImGui::InputInt("#sphere_id", &(instance_->given_sphere_id));
  ImGui::SameLine();
  if (ImGui::SmallButton("show/clear")) {
    instance_->show_one_sphere(*(instance_->all_medial_spheres),
                               instance_->given_sphere_id,
                               true, /*is_show_radius*/
                               true  /*is_clear_all*/
    );
    instance_->show_pin_points(*(instance_->all_medial_spheres),
                               instance_->given_sphere_id);
  }
  ImGui::SameLine();
  if (ImGui::SmallButton("add sphere no R")) {
    instance_->show_one_sphere(*(instance_->all_medial_spheres),
                               instance_->given_sphere_id,
                               false, /*is_show_radius*/
                               false  /*is_clear_all*/
    );
  }
  ImGui::SameLine();
  if (ImGui::SmallButton("with R")) {
    instance_->show_one_sphere(*(instance_->all_medial_spheres),
                               instance_->given_sphere_id,
                               true, /*is_show_radius*/
                               false /*is_clear_all*/
    );
  }

  //   ImGui::SameLine();
  ImGui::PopItemWidth();
}