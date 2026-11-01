from ._libmat import thin_medial_mesh as _thin_medial_mesh
from ._utils import float_matrix, float_vector, int_matrix


def thin_medial_mesh(C, R, E, F, importance, threshold):
    """Remove medial tets and unimportant faces by topology-preserving thinning.

    Medial meshes computed from a power diagram may contain *medial tets*:
    four spheres whose four triangles are all medial faces. Following
    MATTopo, each tet is removed by collapsing a tet-face simple pair (the
    face with the smallest importance first); then, among the faces of those
    tets, faces with `importance < threshold` are pruned through their free
    edges (face-edge simple pairs). Collapsing simple pairs does not change
    the topology of the medial mesh.

    Parameters
    ----------
    C : (n,3) numpy double array
        Sphere centers.
    R : (n,) numpy double array
        Sphere radii.
    E : (e,2) numpy int array
        Medial edges.
    F : (f,3) numpy int array
        Medial faces.
    importance : (f,) numpy double array
        Importance of each face. MATTopo uses the length of the dual edge of
        the restricted power diagram over the average sphere diameter, as
        returned by `medial_mesh_from_spheres`.
    threshold : float
        Faces of medial tets with `importance >= threshold` are kept.

    Returns
    -------
    E_out : (e',2) numpy int array
        Remaining medial edges (indices into `C`).
    F_out : (f',3) numpy int array
        Remaining medial faces (indices into `C`).

    Notes
    -----
    Only faces that belong to a medial tet are candidates for pruning; a
    medial mesh without tets is returned unchanged.

    C++: `void libmat::thin_medial_mesh(C, R, E, F, importance, threshold,
    E_out, F_out)` in `<libmat/thin_medial_mesh.h>`.

    Examples
    --------
    ```python
    import numpy as np
    import libmat
    # four spheres forming a medial tet, plus one extra face
    C = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]])
    R = np.full(5, 0.3)
    F = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3], [1, 2, 4]])
    importance = np.array([0.9, 0.1, 0.2, 0.8, 1.0])
    E_out, F_out = libmat.thin_medial_mesh(C, R, None, F, importance, 0.5)
    print(F_out)   # the low-importance faces of the tet are gone
    ```
    """
    C = float_matrix(C, 3, "C")
    R = float_vector(R, C.shape[0], "R")
    E = int_matrix(E, 2, "E", allow_none=True)
    F = int_matrix(F, 3, "F")
    importance = float_vector(importance, F.shape[0], "importance")
    return _thin_medial_mesh(C, R, E, F, importance, float(threshold))
