"""Stage 2: GF(2) linear algebra matches a NumPy oracle; Pallas == jnp."""

import numpy as np
import pytest

from ic.linalg import (
    gf2_matmul,
    gf2_nullspace_vector,
    gf2_rank,
    gf2_solve,
)
from ic.pallas_kernels import bitmatmul_mod2


def _np_rank_gf2(A):
    A = A.copy().astype(np.int32) % 2
    R, C = A.shape
    r = 0
    for c in range(C):
        piv = None
        for i in range(r, R):
            if A[i, c]:
                piv = i
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        for i in range(R):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        r += 1
    return r


@pytest.mark.parametrize("seed", range(6))
def test_gf2_matmul_matches_numpy(seed):
    rng = np.random.default_rng(seed)
    A = rng.integers(0, 2, size=(17, 23), dtype=np.int32)
    B = rng.integers(0, 2, size=(23, 13), dtype=np.int32)
    got = np.asarray(gf2_matmul(A, B))
    assert np.array_equal(got, (A @ B) % 2)


@pytest.mark.parametrize("seed", range(6))
def test_gf2_rank_matches_numpy(seed):
    rng = np.random.default_rng(100 + seed)
    A = rng.integers(0, 2, size=(20, 16), dtype=np.int32)
    assert gf2_rank(A) == _np_rank_gf2(A)


@pytest.mark.parametrize("seed", range(6))
def test_gf2_solve_roundtrip(seed):
    rng = np.random.default_rng(200 + seed)
    # build a consistent system from a known x
    A = rng.integers(0, 2, size=(24, 16), dtype=np.int32)
    x = rng.integers(0, 2, size=(16,), dtype=np.int32)
    b = (A @ x) % 2
    sol = gf2_solve(A, b)
    assert sol is not None
    assert np.array_equal((A @ sol) % 2, b)


def test_gf2_solve_detects_inconsistent():
    A = np.array([[1, 0], [1, 0]], dtype=np.int32)
    b = np.array([0, 1], dtype=np.int32)  # x0 = 0 and x0 = 1
    assert gf2_solve(A, b) is None


@pytest.mark.parametrize("seed", range(6))
def test_nullspace_vector(seed):
    rng = np.random.default_rng(300 + seed)
    # wide matrix -> guaranteed nontrivial kernel
    A = rng.integers(0, 2, size=(8, 20), dtype=np.int32)
    v = gf2_nullspace_vector(A)
    assert v is not None
    assert np.any(v)
    assert np.array_equal((A @ v) % 2, np.zeros(A.shape[0], dtype=np.int32))


def test_pallas_bitmatmul_matches_jnp():
    rng = np.random.default_rng(7)
    A = rng.integers(0, 2, size=(128, 96), dtype=np.int32)
    B = rng.integers(0, 2, size=(96, 64), dtype=np.int32)
    got = np.asarray(bitmatmul_mod2(A, B, block_m=64, block_n=32))
    assert np.array_equal(got, (A @ B) % 2)
