"""Bit-exact agreement with the Rust library's packing and filter hash.

These two functions are the seam at which a table built by this backend
has to be readable by the Rust/CUDA backends, so they are pinned to
explicit values and to the formulas transcribed in :mod:`tpu.ic.reference`.
"""

from ic.reference import pack, pair_filter_hash

U64 = (1 << 64) - 1


def test_pack_formula():
    # O -> 0
    assert pack((0, 0, True)) == 0
    # (x, y): 2(x+1) + sign, sign = (y > x^y)
    # x=5 (0b101), y=2 (0b010): x^y = 7, y=2 < 7 -> sign 0
    assert pack((5, 2, False)) == ((5 + 1) << 1) | 0
    # x=4, y=7: x^y = 3, y=7 > 3 -> sign 1
    assert pack((4, 7, False)) == ((4 + 1) << 1) | 1
    # real points never collide with the O sentinel (0): smallest is 2
    assert pack((0, 1, False)) >= 2


def test_pair_filter_hash_reference_vectors():
    # recomputed directly from the two constants in koblitz_index_calculus.rs
    def ref(key):
        h = (key * 0xFF51AFD7ED558CCD) & U64
        h ^= h >> 33
        return (h * 0xC4CEB9FE1A85EC53) & U64

    for key in [0, 1, 2, 3, 42, 1 << 40, (1 << 62) + 7, U64]:
        assert pair_filter_hash(key) == ref(key), key


def test_pack_sign_separates_P_from_negP():
    # -P = (x, x ^ y) on a binary curve; pack must separate them when they
    # differ, i.e. whenever x != 0.
    x, y = 6, 3
    negy = x ^ y
    assert pack((x, y, False)) != pack((x, negy, False)) or y == negy
