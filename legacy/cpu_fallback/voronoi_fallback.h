#ifndef CPU_FALLBACK
#define CPU_FALLBACK

#include <vector>

#include "../rpd3d/voronoi_defs.h"

void initialize_geogram(int argc, char** argv);
void terminate_geogram();

/*
 * Fallback CPU implementation for failed vertices.
 * Recomputes Voronoi cell v for each v such that stat[v] != success.
 */
void fallback_voro_diagram_CPU(std::vector<ConvexCellHost>&);

#endif
