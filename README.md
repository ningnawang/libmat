# LibMAT
Some commonly used libraries for computing the 3D medial axis transform (MAT) given a 3D triangle mesh (or tetrahedral mesh).

Starter code is from the paper "[MATFP: Computing Medial Axis Transform with Feature Preservation via Restricted Power Diagram](https://github.com/ningnawang/MATFP)".

The extended work "[MATTopo: Topology-preserving Medial Axis Transform with Restricted Power Diagram](https://github.com/ningnawang/mattopo)" heavily replies on this repo using tag **v0.0.1**.


## Attribution
If you use our library in your research paper, please cite us! You can use the bibtex block below:

```
@misc{gpytoolbox,
  title = {{LibMAT}: A C++ Library for Medial Axis Transform},
  author = {Ningna Wang},
  note = {https://github.com/ningnawang/libmat},
  year = {2024}
}
```

## 1. Lib using GPU (set option **LIBMAT_IS_USE_GPU** as ON)

### 1.1. dist2mat 
Given a 3D sample, compute its closest medial element (sphere/cone/slab) on the given medial mesh.

### 1.2. rpd3d & rpd3d_api
Libs and APIs for computing 3D RPD using CUDA.

### 1.3. matfun_fix
- fix_common
- fix_topo
- fix_extf
- fix_intf
- fix_geo

### 1.4. IO_CUDA
Some IO related libs related to CUDA 3D RPD output.

## 2. Lib NOT using GPU (set option **LIBMAT_IS_USE_GPU** as OFF)
### 2.1. inputs
- SurfaceMesh (extends GEO::Mesh)
- TetMesh
- AABBWrapper
- Sharp/Concave Feature Detection

### 2.2. matbase
- medial_sphere
- medial_primitives
- medial_mesh


### 2.4. matfun
- shrinking
- updating (sphere-optimization)
- thinning

### 2.3. IO_CXX
Some IO related libs.

### 2.5. rpd3d_base
- 3D triangulation using CGAL
- ConvexCellHost

## TODO list:
1. ~~Update CMakeLists.txt to remove some unused external dependencies, such as `polyscope`.~~
2. write more docs 
