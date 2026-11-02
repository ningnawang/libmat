from ._libmat import distance_to_medial_mesh as _distance_to_medial_mesh
from ._libmat import has_cuda
from ._utils import float_matrix, float_vector, int_matrix


def distance_to_medial_mesh(P, C, R, E=None, F=None, device="auto"):
    """Signed distance from points to the volume of a medial mesh.

    A medial mesh represents the union of its medial spheres, medial cones
    (convex hulls of the two spheres of an edge) and medial slabs (convex
    hulls of the three spheres of a face). This function returns, for every
    query point, the distance to that union and the closest primitive. It is
    the building block for measuring how well a medial mesh approximates a
    shape (e.g. the one-sided Hausdorff distance from surface samples).

    Parameters
    ----------
    P : (p,3) numpy double array
        Query points.
    C : (n,3) numpy double array
        Sphere centers.
    R : (n,) numpy double array
        Sphere radii.
    E : (e,2) numpy int array, optional (default None)
        Medial edges.
    F : (f,3) numpy int array, optional (default None)
        Medial faces.
    device : {"auto", "cpu", "cuda"}, optional (default "auto")
        `"cuda"` runs the CUDA kernel (single precision, requires a build
        with CUDA), `"cpu"` runs a multi-threaded double precision
        implementation, `"auto"` uses CUDA when available.

    Returns
    -------
    D : (p,) numpy double array
        Signed distances: positive outside the union, negative inside (inside
        it is the largest penetration depth over the primitives, not an exact
        signed distance).
    prim_type : (p,) numpy int array
        Type of the closest primitive: `0` sphere, `1` cone, `2` slab.
    prim_id : (p,) numpy int array
        Index of the closest primitive in `C`, `E` or `F` (by `prim_type`).

    Notes
    -----
    The cost is `O(p * (n + e + f))`: every point is tested against every
    primitive.

    C++: `void libmat::distance_to_medial_mesh(P, C, R, E, F, use_cuda, D,
    prim_type, prim_id)` in `<libmat/distance_to_medial_mesh.h>`.

    Examples
    --------
    ```python
    import numpy as np
    import libmat
    SV, SF = libmat.read_surface_mesh("data/bug.obj")
    V, R, E, F = libmat.read_ma("data/bug.ma")
    D, prim_type, prim_id = libmat.distance_to_medial_mesh(SV, V, R, E, F)
    print("max |distance| surface -> MAT:", np.abs(D).max())
    ```
    """
    if device not in ("auto", "cpu", "cuda"):
        raise ValueError('device must be "auto", "cpu" or "cuda"')
    if device == "cuda" and not has_cuda():
        raise RuntimeError("libmat was built without CUDA, use device='cpu'")
    use_cuda = device == "cuda" or (device == "auto" and has_cuda())
    P = float_matrix(P, 3, "P")
    C = float_matrix(C, 3, "C")
    R = float_vector(R, C.shape[0], "R")
    E = int_matrix(E, 2, "E", allow_none=True)
    F = int_matrix(F, 3, "F", allow_none=True)
    return _distance_to_medial_mesh(P, C, R, E, F, use_cuda)
