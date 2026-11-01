from ._libmat import read_ma as _read_ma


def read_ma(path):
    """Read a medial mesh from a `.ma` file.

    The `.ma` format is used by MATFP, MATTopo and Q-MAT:

        #V #E #F
        v x y z r      (one line per medial sphere)
        e i j          (one line per medial edge / cone)
        f i j k        (one line per medial face / slab)

    Parameters
    ----------
    path : str
        Path to the `.ma` file.

    Returns
    -------
    V : (n,3) numpy double array
        Medial sphere centers.
    R : (n,) numpy double array
        Medial sphere radii.
    E : (e,2) numpy int array
        Medial edges (each one spans a medial cone).
    F : (f,3) numpy int array
        Medial faces (each one spans a medial slab).

    Raises
    ------
    ValueError
        If the file cannot be read or is malformed.

    See Also
    --------
    write_ma : write a medial mesh.

    Notes
    -----
    C++: `bool libmat::read_ma(path, V, R, E, F)` in `<libmat/read_ma.h>`.

    Examples
    --------
    ```python
    import libmat
    V, R, E, F = libmat.read_ma("data/bug.ma")
    print(V.shape, E.shape, F.shape)   # (5826, 3) (19797, 2) (13972, 3)
    ```
    """
    return _read_ma(str(path))
