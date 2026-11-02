#include "voronoi_fallback.h"

#include "../inputs/params.h"
#include "Delaunay_psm.h"
#include "convex_cell.h"

void initialize_geogram(int argc, char** argv) {
  GEO::initialize();
  GEO::CmdLine::import_arg_group("standard");
  GEO::CmdLine::import_arg_group("algo");
  GEO::CmdLine::set_arg("log:quiet", true);
  GEO::CmdLine::set_arg("algo:delaunay", "BDEL");
  GEO::CmdLine::parse(argc, argv);
  // Use lexicographic order for symbolic perturbations.
  // We do that instead of point address, because our
  // points are local variables, generated on the fly,
  // with randomly varying addresses.
  GEO::PCK::set_SOS_mode(GEO::PCK::SOS_LEXICO);
}

void fallback_voro_diagram_CPU(

) {}