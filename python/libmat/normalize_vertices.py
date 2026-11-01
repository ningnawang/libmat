from ._libmat import normalize_vertices as _normalize_vertices
from ._utils import float_matrix


def normalize_vertices(V, scale=1000.0):
    """Normalize vertices into the box `[0, scale]^3`.

    The bounding box is translated to the origin and uniformly scaled so that
    its longest side has length `scale`:

        VN = (V - offset) * factor

    This is the normalization of MATFP/MATTopo (`scale = 1000`). The
    functions of libmat normalize their inputs internally and return results
    in the input coordinates, so you only need this function to compare with
    results produced by those pipelines (e.g. `.ma` files).

    Parameters
    ----------
    V : (n,3) numpy double array
        Vertex positions.
    scale : float, optional (default 1000.0)
        Length of the longest side of the bounding box after normalization.

    Returns
    -------
    VN : (n,3) numpy double array
        Normalized vertex positions.
    offset : (3,) numpy double array
        Minimum corner of the bounding box of `V`.
    factor : float
        Uniform scaling factor.

    Notes
    -----
    C++: `void libmat::normalize_vertices(V, scale, VN, offset, factor)` in
    `<libmat/normalize_vertices.h>`.

    Examples
    --------
    ```python
    import libmat
    V, T = libmat.read_tet_mesh("data/lblock.tet")
    VN, offset, factor = libmat.normalize_vertices(V)
    V_back = VN / factor + offset
    ```
    """
    VN, offset, factor = _normalize_vertices(float_matrix(V, 3, "V"), float(scale))
    return VN, offset.reshape(-1), factor
