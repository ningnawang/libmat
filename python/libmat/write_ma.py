from ._libmat import write_ma as _write_ma
from ._utils import float_matrix, float_vector, int_matrix


def write_ma(path, V, R, E=None, F=None):
    """Write a medial mesh to a `.ma` file.

    See `read_ma` for a description of the format.

    Parameters
    ----------
    path : str
        Output path (overwritten if it exists).
    V : (n,3) numpy double array
        Medial sphere centers.
    R : (n,) numpy double array
        Medial sphere radii.
    E : (e,2) numpy int array, optional (default None)
        Medial edges.
    F : (f,3) numpy int array, optional (default None)
        Medial faces.

    Raises
    ------
    ValueError
        If the inputs are inconsistent or the file cannot be written.

    Notes
    -----
    C++: `bool libmat::write_ma(path, V, R, E, F)` in `<libmat/write_ma.h>`.

    Examples
    --------
    ```python
    import libmat
    V, F = libmat.read_surface_mesh("data/lblock.obj")
    C, R, _ = libmat.shrink_spheres(V, F, num_spheres=500)
    libmat.write_ma("spheres.ma", C, R)   # spheres only, no edges/faces
    ```
    """
    V = float_matrix(V, 3, "V")
    R = float_vector(R, V.shape[0], "R")
    E = int_matrix(E, 2, "E", allow_none=True)
    F = int_matrix(F, 3, "F", allow_none=True)
    _write_ma(str(path), V, R, E, F)
