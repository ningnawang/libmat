import numpy as np

import libmat
from conftest import box


def test_normalize_vertices():
    V = np.random.default_rng(1).random((50, 3)) * [3, 1, 2] + 5
    VN, offset, factor = libmat.normalize_vertices(V, 1000)
    assert np.allclose(VN.min(0), 0)
    assert np.isclose(VN.max(), 1000)
    assert np.allclose(VN / factor + offset, V)


def test_extract_surface(data_path):
    V, T = libmat.read_tet_mesh(data_path("lblock.tet"))
    SV, SF, sf2tet = libmat.extract_surface(V, T)
    assert SF.shape == (2800, 3)
    assert np.allclose(SV, V[sf2tet])
    # closed and outward oriented: divergence theorem gives the volume 3
    a, b, c = SV[SF[:, 0]], SV[SF[:, 1]], SV[SF[:, 2]]
    vol = np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6
    assert np.isclose(vol, 3.0)


def test_detect_sharp_features_box():
    V, F = box()
    E, E_type, E_line, corners, corner_type = libmat.detect_sharp_features(V, F)
    assert E.shape == (12, 2)
    assert np.all(E_type == 1)
    assert len(np.unique(E_line)) == 12
    assert sorted(corners) == list(range(8)) and np.all(corner_type == 1)


def test_detect_sharp_features_lblock(data_path):
    V, F = libmat.read_surface_mesh(data_path("lblock.obj"))
    E, E_type, E_line, corners, corner_type = libmat.detect_sharp_features(V, F)
    concave = E[E_type == 2]
    # the concave edge is the segment x = y = 1, split into 10 grid edges
    assert len(concave) == 10
    assert np.allclose(V[concave.ravel()][:, :2], 1.0)
    # sharp and concave lines have different ids
    assert not set(E_line[E_type == 1]) & set(E_line[E_type == 2])
