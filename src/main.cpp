#include ""

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include "IO/IO_CXX/io.h"
#include "input_types.h"
#include "main_gui.h"
#include "medial_mesh.h"
#include "params.h"
#include "sharp_feature_detection.h"
#include "shrinking.h"
#include "triangulation.h"

int main(int argc, char** argv) {
  if (1 > argc) {
    std::cerr << "Usage: " << argv[0]
              << "./bin/MAT_MODULES <tet_mesh.msh/medial_mesh.ma>" << std::endl;
    return 1;
  }

  std::srand(RAN_SEED);  // set random seed
  std::string mesh_file_path = argv[1];

  // loading extra file
  // 0: no file
  // 1: sph file for spheres
  // 2: ma file for medial mesh
  // 3: non_cad model
  std::string file_ext = get_file_ext(mesh_file_path);
  int load_file_flag = 0;
  if (file_ext == "msh") {
    return 1;
  } else if (file_ext == "sph") {
    printf("loading sph file: %s\n", mesh_file_path.c_str());
    load_file_flag = 1;
  } else if (file_ext == "ma") {
    printf("loading ma file: %s\n", mesh_file_path.c_str());
    load_file_flag = 2;
  }

  Parameter params;
  std::vector<MedialSphere> all_medial_spheres;
  MedialMesh mmesh;

  // for Gui
  MainGuiWindow main_gui;
  main_gui.set_params(params);
  main_gui.set_all_medial_spheres(all_medial_spheres);
  main_gui.set_medial_mesh(mmesh);

  if (load_file_flag = 2) {
    load_matfp(mesh_file_path, all_medial_spheres, mmesh);
    unnormalize_matfp(params, mmesh);
    renormalize_matfp(params, mmesh);
  }
  main_gui.show();
}