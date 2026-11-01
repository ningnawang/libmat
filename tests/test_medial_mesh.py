import numpy as np
import pytest

import libmat


def brute_force_distance(p, C, R, E, F, n=200):
    best = np.min(np.linalg.norm(C - p, axis=1) - R)
    t = np.linspace(0, 1, 2001)[:, None]
    for i, j in E:
        c = C[i] * (1 - t) + C[j] * t
        r = R[i] * (1 - t[:, 0]) + R[j] * t[:, 0]
        best = min(best, np.min(np.linalg.norm(c - p, axis=1) - r))
    a, b = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n))
    a, b = a.ravel(), b.ravel()
    keep = a + b <= 1
    a, b = a[keep], b[keep]
    w = np.stack([a, b, 1 - a - b], axis=1)
    for f in F:
        c = w @ C[f]
        r = w @ R[f]
        best = min(best, np.min(np.linalg.norm(c - p, axis=1) - r))
    return best


def test_distance_to_spheres_only():
    rng = np.random.default_rng(0)
    C = rng.random((20, 3))
    R = rng.random(20) * 0.1
    P = rng.random((100, 3)) * 2 - 0.5
    D, prim_type, prim_id = libmat.distance_to_medial_mesh(P, C, R, device="cpu")
    expected = np.linalg.norm(P[:, None] - C[None], axis=2) - R
    assert np.allclose(D, expected.min(1))
    assert np.all(prim_type == 0)
    assert np.array_equal(prim_id, expected.argmin(1))


@pytest.mark.parametrize("seed", range(5))
def test_distance_matches_brute_force(seed):
    rng = np.random.default_rng(seed)
    C = rng.normal(size=(5, 3))
    R = rng.uniform(0.2, 0.8, size=5)
    E = np.array([[0, 1], [1, 2], [3, 4]])
    F = np.array([[0, 1, 2], [1, 3, 4]])
    P = rng.normal(size=(30, 3)) * 2
    D, _, _ = libmat.distance_to_medial_mesh(P, C, R, E, F, device="cpu")
    expected = np.array([brute_force_distance(p, C, R, E, F) for p in P])
    # brute force is a sampled upper bound
    assert np.all(D <= expected + 1e-9)
    assert np.allclose(D, expected, atol=2e-3)


def test_distance_surface_to_mat(data_path):
    SV, SF = libmat.read_surface_mesh(data_path("bug.obj"))
    V, R, E, F = libmat.read_ma(data_path("bug.ma"))
    D, prim_type, prim_id = libmat.distance_to_medial_mesh(SV, V, R, E, F)
    # the MATTopo result approximates the surface within 1% of the bbox
    assert np.abs(D).max() < 10.0
    assert set(np.unique(prim_type)) <= {0, 1, 2}


def test_distance_cuda_unavailable():
    if libmat.has_cuda():
        pytest.skip("built with CUDA")
    with pytest.raises(RuntimeError):
        libmat.distance_to_medial_mesh(np.zeros((1, 3)), np.zeros((1, 3)),
                                       np.ones(1), device="cuda")


def test_slab_tangent_triangles():
    C = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0]])
    R = np.array([0.1, 0.1, 0.1])
    T0, T1, valid = libmat.slab_tangent_triangles(C, R, np.array([[0, 1, 2]]))
    assert valid.all()
    # equal radii: the tangent triangles are the face offset by +-r along the
    # face normal (+z for this winding)
    assert np.allclose(T0[0], C + [0, 0, 0.1])
    assert np.allclose(T1[0], C - [0, 0, 0.1])


def test_slab_tangent_triangles_tangent(data_path):
    V, R, E, F = libmat.read_ma(data_path("bug.ma"))
    T0, T1, valid = libmat.slab_tangent_triangles(V, R, F[:500])
    assert valid.mean() > 0.9
    # every corner of a tangent triangle lies on the sphere of its vertex
    Fv = F[:500][valid]
    for T in (T0[valid], T1[valid]):
        d = np.linalg.norm(T - V[Fv], axis=2)
        assert np.allclose(d, R[Fv], rtol=1e-6, atol=1e-6)
        # and its plane is tangent to all three spheres
        nT = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
        nT /= np.linalg.norm(nT, axis=1, keepdims=True)
        h = np.abs(np.einsum("fkd,fd->fk", V[Fv] - T[:, :1], nT))
        assert np.allclose(h, R[Fv], rtol=1e-5, atol=1e-5)
    # T0 is on the side of the face normal, T1 on the other side
    n = np.cross(V[Fv[:, 1]] - V[Fv[:, 0]], V[Fv[:, 2]] - V[Fv[:, 0]])
    center = V[Fv].mean(1)
    assert np.all(np.einsum("ij,ij->i", T0[valid].mean(1) - center, n) >= -1e-9)
    assert np.all(np.einsum("ij,ij->i", T1[valid].mean(1) - center, n) <= 1e-9)


def test_thin_medial_mesh():
    C = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]])
    R = np.full(5, 0.3)
    F = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3], [1, 2, 4]])
    importance = np.array([0.9, 0.1, 0.2, 0.8, 1.0])
    E_out, F_out = libmat.thin_medial_mesh(C, R, None, F, importance, 0.5)
    assert sorted(map(tuple, F_out)) == [(0, 1, 2), (1, 2, 3), (1, 2, 4)]
    # remaining edges are exactly the edges of the remaining faces
    edges = {tuple(sorted(e)) for f in F_out for e in [f[[0, 1]], f[[1, 2]], f[[0, 2]]]}
    assert {tuple(sorted(e)) for e in E_out} == edges


def test_thin_medial_mesh_without_tets(data_path):
    V, R, E, F = libmat.read_ma(data_path("bug.ma"))
    E_out, F_out = libmat.thin_medial_mesh(V, R, E, F, np.zeros(len(F)), 0.5)
    # the MATTopo result is already thin: nothing to remove
    assert len(F_out) == len(F) and len(E_out) == len(E)


def test_thin_closed_tet_complex():
    # boundary of a 4-simplex: 5 tets, no free face; must not hang
    C = np.random.default_rng(0).random((5, 3))
    F = np.array([f for f in __import__("itertools").combinations(range(5), 3)])
    E_out, F_out = libmat.thin_medial_mesh(C, np.full(5, 0.1), None, F,
                                           np.ones(len(F)), 0.5)
    assert len(F_out) <= len(F)


def test_medial_mesh_from_spheres_requires_cuda(data_path):
    if libmat.has_cuda():
        pytest.skip("built with CUDA")
    V, T = libmat.read_tet_mesh(data_path("lblock.tet"))
    with pytest.raises(RuntimeError):
        libmat.medial_mesh_from_spheres(V, T, np.zeros((1, 3)), np.ones(1))
