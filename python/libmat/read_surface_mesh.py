from ._libmat import read_surface_mesh as _read_surface_mesh


def read_surface_mesh(path):
    """Read a triangle mesh from a file.

    Any format supported by geogram can be read (`.obj`, `.off`, `.ply`,
    `.stl`, `.mesh`, `.geogram`, ...). Polygonal faces are triangulated.

    Parameters
    ----------
    path : str
        Path to the surface mesh file.

    Returns
    -------
    V : (n,3) numpy double array
        Vertex positions.
    F : (m,3) numpy int array
        Triangles, as 0-based indices into `V`.

    Raises
    ------
    ValueError
        If the file cannot be read.

    Notes
    -----
    C++: `bool libmat::read_surface_mesh(path, V, F)` in
    `<libmat/read_surface_mesh.h>`.

    Examples
    --------
    ```python
    import libmat
    V, F = libmat.read_surface_mesh("data/lblock.obj")
    print(V.shape, F.shape)   # (1402, 3) (2800, 3)
    ```
    """
    return _read_surface_mesh(str(path))
