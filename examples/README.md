# LibMAT examples

Python scripts showing how to use `libmat`. Run them from the repository root
after installing the package (`pip install .`):

```bash
python examples/03_shrink_spheres.py
```

Every script opens a [polyscope](https://polyscope.run/py/) viewer if it is
installed (`pip install polyscope`); pass `--no-gui` to only print results.
Outputs are written to `examples/output/`.

| script | what it shows |
| --- | --- |
| `01_read_and_extract_surface.py` | read a tet mesh, extract its boundary |
| `02_detect_sharp_features.py` | sharp / concave edges, feature lines, corners |
| `03_shrink_spheres.py` | medial spheres with the ball-shrinking algorithm |
| `04_regular_triangulation.py` | neighbors of medial spheres (power diagram dual) |
| `05_medial_mesh_visualization.py` | read a `.ma` file, medial slabs envelope |
| `06_distance_to_medial_mesh.py` | approximation error of a medial mesh (CPU / CUDA) |
| `07_thin_medial_mesh.py` | topology-preserving thinning |
| `08_medial_mesh_from_spheres_cuda.py` | full pipeline with the CUDA restricted power diagram |
