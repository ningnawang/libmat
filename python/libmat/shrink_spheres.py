from ._libmat import shrink_spheres as _shrink_spheres
from ._utils import float_matrix, int_matrix

#: Names of the sphere types returned by `shrink_spheres` (`SphereType` in
#: libmat/matbase/medial_sphere.h).
SPHERE_TYPES = {
    -3: "T_1_N (corner)",
    -2: "T_1_2 (sharp edge)",
    -1: "T_UNK",
    2: "T_2 (medial sheet)",
    3: "T_N (seam)",
    4: "T_N_JUNC (junction)",
    10: "T_c (concave line pin)",
    11: "T_2_c (concave line)",
    12: "T_N_c (internal feature)",
    13: "T_3_c (concave corner)",
    14: "T_X_c (concave edge)",
}


def shrink_spheres(V, F, num_spheres=1000, preserve_features=True,
                   thres_concave=0.2, thres_convex=40.0, max_iter=-1, seed=200):
    """Sample medial spheres with the ball-shrinking algorithm.

    For each random surface sample `p` with inward normal `n`, a large
    sphere tangent to the surface at `p` is shrunk until it touches the
    surface at a second point and becomes maximal inside the shape
    [Ma et al. 2012]. Its center lies (approximately) on the medial axis.
    With `preserve_features=True`, sharp and concave features are detected
    first and spheres stop at concave edges.

    Parameters
    ----------
    V : (n,3) numpy double array
        Vertex positions of a closed, manifold triangle mesh.
    F : (m,3) numpy int array
        Triangles.
    num_spheres : int, optional (default 1000)
        Number of surface samples. At most one sample is drawn per triangle,
        so at most `m` spheres are created; spheres that fail to converge or
        duplicate an existing sphere are discarded.
    preserve_features : bool, optional (default True)
        Detect sharp/concave features and take them into account.
    thres_concave : float, optional (default 0.2)
        Concave edge threshold, see `detect_sharp_features`.
    thres_convex : float, optional (default 40.0)
        Sharp edge threshold in degrees, see `detect_sharp_features`.
    max_iter : int, optional (default -1)
        Maximum number of shrinking iterations per sphere (-1: no limit).
    seed : int, optional (default 200)
        Seed of the random surface sampling.

    Returns
    -------
    C : (s,3) numpy double array
        Sphere centers, in the coordinates of `V`.
    R : (s,) numpy double array
        Sphere radii.
    types : (s,) numpy int array
        Sphere types, see `libmat.SPHERE_TYPES`.

    See Also
    --------
    regular_triangulation : neighborhood of the spheres (power diagram dual).
    medial_mesh_from_spheres : connect the spheres into a medial mesh (CUDA).

    Notes
    -----
    Internally the mesh is normalized to `[0, 1000]^3`, the scale the default
    parameters are tuned for; outputs are mapped back to the input
    coordinates.

    C++: `void libmat::shrink_spheres(V, F, options, C, R, types)` in
    `<libmat/shrink_spheres.h>`.

    Examples
    --------
    ```python
    import libmat
    V, F = libmat.read_surface_mesh("data/mesh_69058.obj")
    C, R, types = libmat.shrink_spheres(V, F, num_spheres=1000)
    print(C.shape, R.min(), R.max())
    ```
    """
    V = float_matrix(V, 3, "V")
    F = int_matrix(F, 3, "F")
    return _shrink_spheres(V, F, int(num_spheres), bool(preserve_features),
                           float(thres_concave), float(thres_convex),
                           int(max_iter), int(seed))
