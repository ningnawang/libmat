import numpy as np
import pytest

import libmat
from conftest import box


def inside_lblock(P, eps=1e-9):
    x, y, z = P.T
    in_a = (x > -eps) & (x < 2 + eps) & (y > -eps) & (y < 1 + eps)
    in_b = (x > -eps) & (x < 1 + eps) & (y > -eps) & (y < 2 + eps)
    return (in_a | in_b) & (z > -eps) & (z < 1 + eps)


def test_shrink_spheres_lblock(data_path):
    V, F = libmat.read_surface_mesh(data_path("lblock.obj"))
    C, R, types = libmat.shrink_spheres(V, F, num_spheres=500)
    assert 0 < len(C) <= 500
    assert C.shape == (len(R), 3) and types.shape == R.shape
    assert np.all(R > 0)
    assert np.all(inside_lblock(C))
    # the block is 1 thick: no inscribed sphere has radius > 0.5
    assert R.max() <= 0.5 + 1e-6
    assert set(np.unique(types)) <= set(libmat.SPHERE_TYPES)


def lblock_boundary_distance(P):
    """Exact distance to the boundary of the L-block for points inside."""
    x, y, z = P.T
    dz = np.minimum(z, 1 - z)
    dx = np.minimum(x, np.where(y > 1, 1 - x, 2 - x))
    dy = np.minimum(y, np.where(x > 1, 1 - y, 2 - y))
    d = np.minimum.reduce([dx, dy, dz])
    # in the square shared by both boxes, the concave edge x = y = 1 is close
    in_corner = (x < 1) & (y < 1)
    return np.where(in_corner, np.minimum(d, np.hypot(1 - x, 1 - y)), d)


def test_shrink_spheres_is_maximal(data_path):
    # spheres touch the surface: distance from center to the surface = radius
    V, F = libmat.read_surface_mesh(data_path("lblock.obj"))
    C, R, _ = libmat.shrink_spheres(V, F, num_spheres=300)
    assert np.allclose(lblock_boundary_distance(C), R, atol=1e-3)


def test_shrink_spheres_deterministic(data_path):
    V, F = libmat.read_surface_mesh(data_path("lblock.obj"))
    C1, R1, _ = libmat.shrink_spheres(V, F, num_spheres=200, seed=7)
    C2, R2, _ = libmat.shrink_spheres(V, F, num_spheres=200, seed=7)
    assert np.array_equal(C1, C2) and np.array_equal(R1, R2)


def test_shrink_spheres_scale_invariant(data_path):
    V, F = libmat.read_surface_mesh(data_path("lblock.obj"))
    C1, R1, _ = libmat.shrink_spheres(V, F, num_spheres=200)
    C2, R2, _ = libmat.shrink_spheres(V * 10 + 3, F, num_spheres=200)
    assert np.allclose(C1 * 10 + 3, C2, atol=1e-3)
    assert np.allclose(R1 * 10, R2, atol=1e-3)


def test_shrink_spheres_bad_input():
    V, F = box()
    with pytest.raises(ValueError):
        libmat.shrink_spheres(V[:, :2], F)
    with pytest.raises(ValueError):
        libmat.shrink_spheres(V, F + 100)


def test_regular_triangulation_equal_radii():
    # equal radii: regular triangulation == Delaunay triangulation
    rng = np.random.default_rng(3)
    C = rng.random((200, 3))
    R = np.full(len(C), 0.05)
    T, E, valid = libmat.regular_triangulation(C, R)
    assert valid.all()
    assert T.min() >= 0 and T.max() < len(C)
    # Euler characteristic of a triangulated ball: V - E + F - T = 1
    faces = np.unique(np.sort(np.concatenate(
        [T[:, [0, 1, 2]], T[:, [0, 1, 3]], T[:, [0, 2, 3]], T[:, [1, 2, 3]]]),
        axis=1), axis=0)
    assert len(C) - len(E) + len(faces) - len(T) == 1


def test_regular_triangulation_hidden_sphere():
    C = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0.2, 0.2, 0.2]])
    R = np.array([1.0, 1.0, 1.0, 1.0, 0.01])  # small sphere is hidden
    T, E, valid = libmat.regular_triangulation(C, R)
    assert list(valid) == [True, True, True, True, False]
    assert len(T) == 1 and 4 not in E
