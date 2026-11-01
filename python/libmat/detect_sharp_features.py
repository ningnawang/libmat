from ._libmat import detect_sharp_features as _detect_sharp_features
from ._utils import float_matrix, int_matrix


def detect_sharp_features(V, F, thres_concave=0.2, thres_convex=40.0):
    """Detect sharp edges, concave edges and corners of a triangle mesh.

    An edge is *sharp* (convex) if the dihedral angle between its two faces
    deviates from flat by more than `thres_convex` degrees, and *concave* if
    the surface folds inwards along it by more than `thres_concave`. Feature
    edges are grouped into feature lines, split at corners. These are the
    features that a feature-preserving medial axis has to keep.

    Parameters
    ----------
    V : (n,3) numpy double array
        Vertex positions of a closed, manifold triangle mesh.
    F : (m,3) numpy int array
        Triangles.
    thres_concave : float, optional (default 0.2)
        Threshold for concave edges (smaller is more sensitive).
    thres_convex : float, optional (default 40.0)
        Dihedral angle threshold in degrees for sharp edges (smaller is more
        sensitive).

    Returns
    -------
    E : (k,2) numpy int array
        Feature edges, as indices into `V`.
    E_type : (k,) numpy int array
        `1` for a sharp (convex) edge, `2` for a concave edge.
    E_line : (k,) numpy int array
        Id of the feature line containing each edge.
    corners : (c,) numpy int array
        Corner vertices, as indices into `V`.
    corner_type : (c,) numpy int array
        `1` for a sharp corner (adjacent to sharp edges), `2` for a concave
        corner (adjacent to more than two concave edges only).

    Notes
    -----
    The mesh must be closed (no boundary edges).

    C++: `void libmat::detect_sharp_features(V, F, thres_concave,
    thres_convex, E, E_type, E_line, corners, corner_type)` in
    `<libmat/detect_sharp_features.h>`.

    Examples
    --------
    ```python
    import numpy as np
    import libmat
    V, F = libmat.read_surface_mesh("data/mesh_69058.obj")
    E, E_type, E_line, corners, corner_type = libmat.detect_sharp_features(V, F)
    print("sharp edges:  ", np.sum(E_type == 1))
    print("concave edges:", np.sum(E_type == 2))
    print("feature lines:", len(np.unique(E_line)))
    ```
    """
    V = float_matrix(V, 3, "V")
    F = int_matrix(F, 3, "F")
    return _detect_sharp_features(V, F, float(thres_concave), float(thres_convex))
