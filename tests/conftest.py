import os

import pytest

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


@pytest.fixture
def data_path():
    return lambda name: os.path.join(DATA, name)


def box(lo=(0.0, 0.0, 0.0), hi=(2.0, 1.0, 1.0)):
    """12-triangle box, outward oriented."""
    import numpy as np
    (x0, y0, z0), (x1, y1, z1) = lo, hi
    V = np.array([[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                  [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]])
    F = np.array([[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5],
                  [0, 5, 4], [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6],
                  [3, 0, 4], [3, 4, 7]])
    return V, F
