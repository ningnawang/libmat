This directory contains the CPU fallback code that recomputes the
individual Voronoi cells that could not be computed on the GPU (if any).

It provides the function: fallback_voro_diagram_CPU(pts, stat, bary, KNN)

See main() in test_voronoi.cu

