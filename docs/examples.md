---
title: Examples
---

# Examples

The [`examples/`](https://github.com/ningnawang/libmat/tree/main/examples) folder contains one Python script per topic. Run them from the repository root; each one opens a [polyscope](https://polyscope.run/py/) viewer if it is installed (`pip install polyscope`), or pass `--no-gui`:

```bash
python examples/03_shrink_spheres.py
```

## Read a tet mesh and extract its boundary surface

`examples/01_read_and_extract_surface.py`

```python
--8<-- "examples/01_read_and_extract_surface.py"
```

## Detect sharp edges, concave edges and corners of a CAD-like model

`examples/02_detect_sharp_features.py`

```python
--8<-- "examples/02_detect_sharp_features.py"
```

## Sample medial spheres with the ball-shrinking algorithm

`examples/03_shrink_spheres.py`

```python
--8<-- "examples/03_shrink_spheres.py"
```

## Neighborhood of medial spheres: regular triangulation (power diagram dual)

`examples/04_regular_triangulation.py`

```python
--8<-- "examples/04_regular_triangulation.py"
```

## Visualize a medial mesh (.ma) and the shape it represents

`examples/05_medial_mesh_visualization.py`

```python
--8<-- "examples/05_medial_mesh_visualization.py"
```

## Measure how well a medial mesh approximates a surface

`examples/06_distance_to_medial_mesh.py`

```python
--8<-- "examples/06_distance_to_medial_mesh.py"
```

## Topology-preserving thinning of a medial mesh (as in MATTopo)

`examples/07_thin_medial_mesh.py`

```python
--8<-- "examples/07_thin_medial_mesh.py"
```

## Full medial mesh pipeline on the GPU (requires libmat built with CUDA)

`examples/08_medial_mesh_from_spheres_cuda.py`

```python
--8<-- "examples/08_medial_mesh_from_spheres_cuda.py"
```
