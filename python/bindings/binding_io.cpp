#include "common.h"
#include "libmat/read_ma.h"
#include "libmat/read_surface_mesh.h"
#include "libmat/read_tet_mesh.h"
#include "libmat/write_ma.h"

void bind_io(py::module_& m) {
  m.def("read_tet_mesh", [](const std::string& path) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi T;
    bool ok;
    {
      libmat_py::CallGuard guard;
      ok = libmat::read_tet_mesh(path, V, T);
    }
    if (!ok) throw py::value_error("could not read tet mesh: " + path);
    return py::make_tuple(V, T);
  });

  m.def("read_surface_mesh", [](const std::string& path) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    bool ok;
    {
      libmat_py::CallGuard guard;
      ok = libmat::read_surface_mesh(path, V, F);
    }
    if (!ok) throw py::value_error("could not read surface mesh: " + path);
    return py::make_tuple(V, F);
  });

  m.def("read_ma", [](const std::string& path) {
    Eigen::MatrixXd V;
    Eigen::VectorXd R;
    Eigen::MatrixXi E, F;
    bool ok;
    {
      libmat_py::CallGuard guard;
      ok = libmat::read_ma(path, V, R, E, F);
    }
    if (!ok) throw py::value_error("could not read .ma file: " + path);
    return py::make_tuple(V, R, E, F);
  });

  m.def("write_ma", [](const std::string& path, const Eigen::MatrixXd& V,
                       const Eigen::VectorXd& R, const Eigen::MatrixXi& E,
                       const Eigen::MatrixXi& F) {
    bool ok;
    {
      libmat_py::CallGuard guard;
      ok = libmat::write_ma(path, V, R, E, F);
    }
    if (!ok) throw py::value_error("could not write .ma file: " + path);
  });
}
