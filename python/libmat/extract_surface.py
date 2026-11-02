from ._libmat import extract_surface as _extract_surface
from ._utils import float_matrix, int_matrix


def extract_surface(V, T):
    """Extract the boundary triangle mesh of a tetrahedral mesh.

    Boundary faces are the tet faces that belong to a single tetrahedron.
    They are oriented so that their normals point outwards.

    Parameters
    ----------
    V : (n,3) numpy double array
        Tet mesh vertex positions.
    T : (m,4) numpy int array
        Tetrahedra.

    Returns
    -------
    SV : (k,3) numpy double array
        Surface vertex positions.
    SF : (f,3) numpy int array
        Surface triangles, as indices into `SV`.
    sf2tet : (k,) numpy int array
        Index of each surface vertex in `V`, i.e. `SV == V[sf2tet]`.

    Notes
    -----
    C++: `void libmat::extract_surface(V, T, SV, SF, sf2tet)` in
    `<libmat/extract_surface.h>`.

    Examples
    --------
    ```python
    import libmat
    V, T = libmat.read_tet_mesh("data/lblock.tet")
    SV, SF, sf2tet = libmat.extract_surface(V, T)
    print(SF.shape)   # (2800, 3)
    ```
    """
    return _extract_surface(float_matrix(V, 3, "V"), int_matrix(T, 4, "T"))
