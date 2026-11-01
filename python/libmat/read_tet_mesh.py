from ._libmat import read_tet_mesh as _read_tet_mesh


def read_tet_mesh(path):
    """Read a tetrahedral mesh from a .msh, .tet or .vtk file.

    Supports gmsh `.msh` files (e.g. produced by
    [fTetWild](https://github.com/wildmeshing/fTetWild)), the plain-text
    `.tet` format and legacy `.vtk` unstructured grids.

    Parameters
    ----------
    path : str
        Path to the tet mesh file.

    Returns
    -------
    V : (n,3) numpy double array
        Vertex positions.
    T : (m,4) numpy int array
        Tetrahedra, as 0-based indices into `V`.

    Raises
    ------
    ValueError
        If the file cannot be read.

    See Also
    --------
    extract_surface : boundary triangle mesh of a tet mesh.

    Notes
    -----
    C++: `bool libmat::read_tet_mesh(path, V, T)` in `<libmat/read_tet_mesh.h>`.

    Examples
    --------
    ```python
    import libmat
    V, T = libmat.read_tet_mesh("data/lblock.tet")
    print(V.shape, T.shape)   # (3751, 3) (18000, 4)
    ```
    """
    return _read_tet_mesh(str(path))
