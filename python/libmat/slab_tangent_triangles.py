from ._libmat import slab_tangent_triangles as _slab_tangent_triangles
from ._utils import float_matrix, float_vector, int_matrix


def slab_tangent_triangles(C, R, F):
    """Triangles tangent to the three spheres of each medial slab.

    A medial face spans a *medial slab*, the convex hull of its three
    spheres. Its boundary consists of two triangles tangent to all three
    spheres (one on each side of the face) and of pieces of the medial cones
    of its edges. Rendering these triangles is a quick way to visualize the
    shape represented by a medial mesh.

    Parameters
    ----------
    C : (n,3) numpy double array
        Sphere centers.
    R : (n,) numpy double array
        Sphere radii.
    F : (f,3) numpy int array
        Medial faces.

    Returns
    -------
    T0 : (f,3,3) numpy double array
        Tangent triangle on the side of the face normal (right-hand rule over
        `F[i]`). Corner `T0[i, k]` touches sphere `F[i, k]`, so `T0[i]` is
        oriented like `F[i]`, with its normal pointing away from the slab.
    T1 : (f,3,3) numpy double array
        Tangent triangle on the other side. Corner `T1[i, k]` touches sphere
        `F[i, k]`; use `T1[:, [0, 2, 1]]` for normals pointing away from the
        slab.
    valid : (f,) numpy bool array
        `False` where the tangent triangles do not exist (e.g. coincident
        spheres); those rows are zero.

    Notes
    -----
    C++: `void libmat::slab_tangent_triangles(C, R, F, T0, T1, valid)` in
    `<libmat/slab_tangent_triangles.h>`, with `T0`, `T1` stored as `#F by 9`.

    Examples
    --------
    ```python
    import numpy as np
    import libmat
    V, R, E, F = libmat.read_ma("data/bug.ma")
    T0, T1, valid = libmat.slab_tangent_triangles(V, R, F)
    # as an outward oriented triangle soup
    tris = np.concatenate([T0[valid], T1[valid][:, [0, 2, 1]]])
    SV = tris.reshape(-1, 3)
    SF = np.arange(len(SV)).reshape(-1, 3)
    ```
    """
    C = float_matrix(C, 3, "C")
    R = float_vector(R, C.shape[0], "R")
    F = int_matrix(F, 3, "F")
    T0, T1, valid = _slab_tangent_triangles(C, R, F)
    return T0.reshape(-1, 3, 3), T1.reshape(-1, 3, 3), valid.astype(bool)
