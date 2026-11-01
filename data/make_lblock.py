"""Generate the L-shaped block used by the examples and tests.

    python data/make_lblock.py

writes data/lblock.tet (tet mesh) and data/lblock.obj (its boundary). The
block is the union of [0,2]x[0,1]x[0,1] and [0,1]x[1,2]x[0,1]: it has sharp
convex edges, one concave edge (x=y=1) and sharp/concave corners.
Each grid cube is split into 6 tets along its main diagonal (Freudenthal),
which gives a conforming tet mesh.
"""
import itertools
import os

import numpy as np

H = 0.1  # grid spacing


def lblock_tets(h=H):
    n = int(round(2 / h))
    m = int(round(1 / h))
    inside = lambda i, j: (j < m) or (i < m)  # cube (i, j) in the L shape
    vid = {}
    V, T = [], []

    def vertex(i, j, k):
        key = (i, j, k)
        if key not in vid:
            vid[key] = len(V)
            V.append((i * h, j * h, k * h))
        return vid[key]

    for i, j, k in itertools.product(range(n), range(n), range(m)):
        if not inside(i, j):
            continue
        for perm in itertools.permutations(range(3)):
            p = [i, j, k]
            tet = [vertex(*p)]
            for axis in perm:
                p[axis] += 1
                tet.append(vertex(*p))
            T.append(tet)
    V, T = np.array(V), np.array(T)
    # positive orientation
    a, b, c, d = (V[T[:, i]] for i in range(4))
    vol = np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a)
    T[vol < 0] = T[vol < 0][:, [0, 2, 1, 3]]
    return V, T


def boundary(V, T):
    # faces opposite to vertex 3, 2, 1, 0, oriented outwards for positive tets
    faces = np.concatenate([T[:, [0, 2, 1]], T[:, [0, 1, 3]],
                            T[:, [0, 3, 2]], T[:, [1, 2, 3]]])
    key = np.sort(faces, axis=1)
    _, inv, cnt = np.unique(key, axis=0, return_inverse=True,
                            return_counts=True)
    F = faces[cnt[inv.ravel()] == 1]
    used = np.unique(F)
    remap = -np.ones(len(V), dtype=int)
    remap[used] = np.arange(len(used))
    return V[used], remap[F]


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    V, T = lblock_tets()
    with open(os.path.join(here, "lblock.tet"), "w") as f:
        f.write(f"{len(V)} {len(T)}\n")
        for v in V:
            f.write(f"{v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        for t in T:
            f.write(f"4 {t[0]} {t[1]} {t[2]} {t[3]}\n")
    SV, SF = boundary(V, T)
    with open(os.path.join(here, "lblock.obj"), "w") as f:
        for v in SV:
            f.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        for t in SF + 1:
            f.write(f"f {t[0]} {t[1]} {t[2]}\n")
    print(f"lblock.tet: {len(V)} vertices, {len(T)} tets")
    print(f"lblock.obj: {len(SV)} vertices, {len(SF)} triangles")
