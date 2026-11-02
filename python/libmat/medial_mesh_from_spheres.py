from ._libmat import has_cuda
from ._libmat import medial_mesh_from_spheres as _medial_mesh_from_spheres
from ._utils import float_matrix, float_vector, int_matrix


def medial_mesh_from_spheres(V, T, C, R):
    """Connect medial spheres into a medial mesh with the restricted power diagram (CUDA).

    The power diagram of the spheres is clipped by the input tet mesh on the
    GPU (restricted power diagram, RPD). Two spheres are connected by a medial
    edge when their restricted power cells share a face, and three spheres
    span a medial face when their cells share an edge.

    **Requires a build with CUDA** (`LIBMAT_WITH_CUDA=ON`).

    Parameters
    ----------
    V : (n,3) numpy double array
        Tet mesh vertex positions.
    T : (m,4) numpy int array
        Tetrahedra.
    C : (s,3) numpy double array
        Medial sphere centers (same coordinates as `V`).
    R : (s,) numpy double array
        Medial sphere radii.

    Returns
    -------
    E : (e,2) numpy int array
        Medial edges (indices into `C`).
    F : (f,3) numpy int array
        Medial faces (indices into `C`).
    importance : (f,) numpy double array
        Face importance (dual RPD edge length over the average sphere
        diameter), to be used with `thin_medial_mesh`.

    Raises
    ------
    RuntimeError
        If libmat was built without CUDA.

    Notes
    -----
    C++: `void libmat::medial_mesh_from_spheres(V, T, C, R, E, F,
    importance)` in `<libmat/medial_mesh_from_spheres.h>`, implemented in
    `src/cuda/api/`.

    Examples
    --------
    ```python
    import libmat
    V, T = libmat.read_tet_mesh("data/lblock.tet")
    SV, SF, _ = libmat.extract_surface(V, T)
    C, R, _ = libmat.shrink_spheres(SV, SF, num_spheres=1000)
    if libmat.has_cuda():
        E, F, importance = libmat.medial_mesh_from_spheres(V, T, C, R)
        E, F = libmat.thin_medial_mesh(C, R, E, F, importance, 0.3)
    ```
    """
    if not has_cuda():
        raise RuntimeError(
            "medial_mesh_from_spheres requires CUDA: rebuild libmat with "
            "LIBMAT_WITH_CUDA=ON")
    V = float_matrix(V, 3, "V")
    T = int_matrix(T, 4, "T")
    C = float_matrix(C, 3, "C")
    R = float_vector(R, C.shape[0], "R")
    return _medial_mesh_from_spheres(V, T, C, R)
