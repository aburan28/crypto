"""Stage 1: the pair table is built correctly and collects real relations."""

import numpy as np
import pytest

from ic.curve import CurveJax
from ic.instances import factor_base, find_generator, make_curve
from ic.pairtable import PairTable, collect_m3, verify_relation_m3
from ic.reference import INF, pack


def _oracle_table_keys(curve, base):
    """Every P_i + P_j (i<=j) packed, via the scalar oracle."""
    keys = []
    F = len(base)
    for i in range(F):
        for j in range(i, F):
            keys.append(pack(curve.add(base[i], base[j])))
    return sorted(keys)


@pytest.mark.parametrize("n", [9, 11])
def test_table_keys_match_oracle(n):
    curve = make_curve(n, a=0)
    cj = CurveJax(curve)
    base = factor_base(curve, size=24, seed=3)
    table = PairTable.build(cj, base)
    assert sorted(int(k) for k in table.keys) == _oracle_table_keys(curve, base)


@pytest.mark.parametrize("n", [9, 11])
def test_recovered_pairs_sum_correctly(n):
    curve = make_curve(n, a=0)
    cj = CurveJax(curve)
    base = factor_base(curve, size=24, seed=5)
    table = PairTable.build(cj, base)
    # every stored key really is the sum of its recorded pair
    for pos in range(len(table.keys)):
        i, j = table.pairs[pos]
        s = curve.add(base[int(i)], base[int(j)])
        assert pack(s) == int(table.keys[pos])


def test_m3_collection_end_to_end():
    # small curve we can enumerate and solve
    n = 11
    curve = make_curve(n, a=0)
    cj = CurveJax(curve)
    g, order = find_generator(curve)
    base = factor_base(curve, size=40, seed=7)
    table = PairTable.build(cj, base)

    # probe targets R = [a]G for a spread of scalars
    scalars = [3, 17, 29, 44, 58, 71, 90, 101, 120, 133]
    targets = [curve.mul_scalar(g, a) for a in scalars]

    rels = collect_m3(cj, table, base, targets, scalars)
    assert rels, "expected at least one 3-summand relation on this base"
    for rec, a in zip(rels, [r["a"] for r in rels]):
        # map a -> its target
        target = curve.mul_scalar(g, a)
        assert verify_relation_m3(curve, base, rec["points"], target), rec
