import numpy as np
import pytest

import libmat


def test_read_tet_mesh(data_path):
    V, T = libmat.read_tet_mesh(data_path("lblock.tet"))
    assert V.shape == (3751, 3) and T.shape == (18000, 4)
    assert T.min() == 0 and T.max() == len(V) - 1
    # L-block volume = 3 (positively oriented tets)
    a, b, c, d = (V[T[:, i]] for i in range(4))
    vol = np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6
    assert np.all(vol > 0)
    assert vol.sum() == pytest.approx(3.0)


def test_read_surface_mesh(data_path):
    V, F = libmat.read_surface_mesh(data_path("lblock.obj"))
    assert V.shape == (1402, 3) and F.shape == (2800, 3)


def test_read_missing_file():
    with pytest.raises(ValueError):
        libmat.read_tet_mesh("does_not_exist.msh")
    with pytest.raises(ValueError):
        libmat.read_surface_mesh("does_not_exist.obj")
    with pytest.raises(ValueError):
        libmat.read_ma("does_not_exist.ma")


def test_ma_roundtrip(tmp_path, data_path):
    V, R, E, F = libmat.read_ma(data_path("bug.ma"))
    assert V.shape == (5826, 3) and R.shape == (5826,)
    assert E.shape == (19797, 2) and F.shape == (13972, 3)
    out = tmp_path / "out.ma"
    libmat.write_ma(out, V, R, E, F)
    V2, R2, E2, F2 = libmat.read_ma(out)
    assert np.allclose(V, V2) and np.allclose(R, R2)
    assert np.array_equal(E, E2) and np.array_equal(F, F2)


def test_write_ma_spheres_only(tmp_path):
    V = np.random.default_rng(0).random((10, 3))
    R = np.full(10, 0.1)
    libmat.write_ma(tmp_path / "s.ma", V, R)
    V2, R2, E2, F2 = libmat.read_ma(tmp_path / "s.ma")
    assert np.allclose(V, V2) and E2.shape == (0, 2) and F2.shape == (0, 3)


def test_write_ma_bad_input(tmp_path):
    with pytest.raises(ValueError):
        libmat.write_ma(tmp_path / "x.ma", np.zeros((3, 3)), np.zeros(2))
    with pytest.raises(ValueError):
        libmat.write_ma(tmp_path / "x.ma", np.zeros((3, 3)), np.zeros(3),
                        F=np.array([[0, 1, 5]]))
