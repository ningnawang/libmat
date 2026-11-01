from ._libmat import regular_triangulation as _regular_triangulation
from ._utils import float_matrix, float_vector


def regular_triangulation(C, R):
    """Regular triangulation of a set of spheres (dual of the power diagram).

    The sphere with center `c` and radius `r` is a weighted point with weight
    `r^2`. Two spheres are connected in the regular triangulation iff their
    power cells share a face. A sphere is *hidden* (its power cell is empty)
    if it does not appear as a vertex. Computed with CGAL.

    Parameters
    ----------
    C : (n,3) numpy double array
        Sphere centers.
    R : (n,) numpy double array
        Sphere radii.

    Returns
    -------
    T : (t,4) numpy int array
        Finite tetrahedra, as indices into `C`.
    E : (e,2) numpy int array
        Finite edges (sorted pairs of neighboring spheres).
    valid : (n,) numpy bool array
        `False` for hidden spheres.

    Notes
    -----
    MATFP/MATTopo use this triangulation to find the neighbors of each
    medial sphere before clipping its power cell against the input shape.

    C++: `void libmat::regular_triangulation(C, R, T, E, valid)` in
    `<libmat/regular_triangulation.h>`.

    Examples
    --------
    ```python
    import libmat
    V, F = libmat.read_surface_mesh("data/lblock.obj")
    C, R, _ = libmat.shrink_spheres(V, F, num_spheres=1000)
    T, E, valid = libmat.regular_triangulation(C, R)
    print(len(E), "neighboring pairs,", (~valid).sum(), "hidden spheres")
    ```
    """
    C = float_matrix(C, 3, "C")
    R = float_vector(R, C.shape[0], "R")
    T, E, valid = _regular_triangulation(C, R)
    return T, E, valid.astype(bool)
