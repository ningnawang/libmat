"""Input validation helpers shared by the libmat python wrappers."""
import numpy as np


def float_matrix(X, cols, name):
    X = np.asarray(X, dtype=np.float64)
    if X.ndim != 2 or X.shape[1] != cols:
        raise ValueError(f"{name} must be a (n, {cols}) array, got shape {X.shape}")
    return np.ascontiguousarray(X)


def int_matrix(X, cols, name, allow_none=False):
    if X is None:
        if allow_none:
            return np.zeros((0, cols), dtype=np.int32)
        raise ValueError(f"{name} is required")
    X = np.asarray(X)
    if X.size == 0:
        return np.zeros((0, cols), dtype=np.int32)
    if X.ndim != 2 or X.shape[1] != cols:
        raise ValueError(f"{name} must be a (n, {cols}) array, got shape {X.shape}")
    if not np.issubdtype(X.dtype, np.integer):
        raise ValueError(f"{name} must be an integer array")
    return np.ascontiguousarray(X, dtype=np.int32)


def float_vector(x, n, name):
    x = np.asarray(x, dtype=np.float64).reshape(-1)
    if n is not None and x.shape[0] != n:
        raise ValueError(f"{name} must have {n} entries, got {x.shape[0]}")
    return x
