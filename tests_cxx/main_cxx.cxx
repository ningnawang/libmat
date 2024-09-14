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

#include "../src/IO/IO_CXX/io.h"
#include "../src/inputs/params.h"
#include "../src/matfun/shrinking.h"
#include "../src/matfun/updating.h"
#include "main_gui_cxx.h"

// Scale matfp to mattopo, without scaling radius
//  $ ./bin/LIBMAT_TEST_CXX ../data/mattopo/input/bug.off_.msh
//  ../data/matfp/mat_mesh_bug.off__sf_2024-08-17_15:06:47.ma
int main(int argc, char** argv) {
  if (3 != argc) {
    std::cerr << "Usage: " << argv[0]
              << "./bin/LIBMAT_TEST <tet_mesh.msh> "
                 "<matfp_medial_mesh.ma/medial_spheres.sph>"
              << std::endl;
    return 1;
  }

  std::srand(RAN_SEED);  // set random seed
  std::string tet_mesh_file = argv[1];
  std::string other_mesh_file = argv[2];

  // loading extra file
  // 0: no file
  // 1: sph file for spheres
  // 2: ma file for medial mesh
  // 3: non_cad model
  std::string file_ext = get_file_ext(other_mesh_file);
  std::string file_name_no_ext =
      get_only_file_name(other_mesh_file, false /*withExtension*/);
  int load_file_flag = 0;
  if (file_ext == "msh") {
    return 1;
  } else if (file_ext == "sph") {
    printf("loading sph file: %s\n", other_mesh_file.c_str());
    load_file_flag = 1;
  } else if (file_ext == "ma") {
    printf("loading ma file: %s\n", other_mesh_file.c_str());
    load_file_flag = 2;
  }

  // Load tet from file
  TetMesh tet_mesh(tet_mesh_file);
  Parameter params;
  std::vector<MedialSphere> all_medial_spheres;
  MedialMesh mmesh;

  // please set normalize to true for model joint.tet/.xyz
  if (!load_tet(tet_mesh_file, tet_mesh.tet_vertices, tet_mesh.tet_indices,
                true /*normalize*/, params)) {
    std::cerr << tet_mesh_file << ": could not load file" << std::endl;
    return 1;
  }

  if (load_file_flag = 2) {
    // load_matfp(mesh_file_path, all_medial_spheres, mmesh);
    load_mat_clean(other_mesh_file, all_medial_spheres, mmesh);
  }

  // for Gui
  MainGuiWindow main_gui;
  main_gui.set_file_name(file_name_no_ext);
  main_gui.set_params(params);
  main_gui.set_tet_mesh(tet_mesh);
  main_gui.set_all_medial_spheres(all_medial_spheres);
  main_gui.set_medial_mesh(mmesh);

  main_gui.show();
}