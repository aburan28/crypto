#!/usr/bin/env python3
"""Independent bit-serial/Fermat replay of the direct public-point fixture."""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / "rotated_subspace_support_20260925/verify.py"
P = (0x051C99BFA6F18DE467C80C23B98C7994AA,
     0x042EA2D112ECEC71FCF7E000D7EFC978BD)
Q = (0x06C997F3E7F2C66A4A5D2FDA13756A37B1,
     0x04A38D11829D32D347BD0C0F584D546E9A)
POLY = (1 << 131) | (1 << 13) | (1 << 2) | (1 << 1) | 1


def prior():
    spec = importlib.util.spec_from_file_location("challenge_import_independent", PRIOR)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def point(raw):
    return None if raw is None else tuple(int(c, 16) for c in raw)


def rank(vectors):
    pivots = {}
    for value in vectors:
        while value:
            bit = value.bit_length() - 1
            if bit in pivots:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                break
    return len(pivots)


def verify(path: Path):
    fixture = json.loads(path.read_text())
    mod = prior()
    assert fixture["schema"] == "ecc2k130_direct_challenge_point_import_v1"
    assert fixture["source"] == {
        "relations_path": "research/ecc2k130_relations/relations.py",
        "relations_git_blob": "a1df88fbf76fd14e565a5ec6b7d0ea4402f954ed",
        "fastfield_path": "research/ecc2k130_relations/fastfield.py",
        "fastfield_git_blob": "bb57818dd85db7bff049e6f590d818616c5164e6"}
    assert fixture["field"] == {"degree": 131, "polynomial": hex(POLY),
                                "encoding": "polynomial coefficient bit i is z^i",
                                "normal_beta": 3}
    q, lam = mod.Q131, mod.LAMBDA131
    assert q == 680564733841876926932320129493409985129
    assert lam == 196511074115861092422032515080945363956
    assert fixture["curve"] == {"equation": "y^2+x*y=x^3+1", "cofactor": 4,
                                "order": str(4 * q), "subgroup_order": str(q),
                                "frobenius_scalar": str(lam)}
    assert mod.P131 == POLY and mod.source_order_by_recurrence(131) == 4 * q
    field = mod.GF(131, POLY)
    curve = mod.E(field)
    conjugates = []
    beta = 3
    for _ in range(131):
        conjugates.append(beta)
        beta = field.square(beta)
    assert beta == 3 and rank(conjugates) == 131
    assert field.trace(3) == 1

    def coordinates(entry, expected):
        actual = point(entry["polynomial"])
        assert actual == expected and actual is not None
        assert all(0 <= c < 1 << 131 for c in actual)
        assert curve.on(actual)
        normal_masks = tuple(int(c, 16) for c in entry["normal_beta3"])
        assert all(0 <= c < 1 << 131 for c in normal_masks)
        for poly, mask in zip(actual, normal_masks):
            # An independently checked full-rank basis makes this encoding unique.
            decoded = 0
            for i, basis in enumerate(conjugates):
                if (mask >> i) & 1:
                    decoded ^= basis
            assert decoded == poly
        assert entry["trace_x"] == field.trace(actual[0])
        assert entry["trace_x"] == (normal_masks[0].bit_count() & 1)
        return actual

    p = coordinates(fixture["P"], P)
    target = coordinates(fixture["Q"], Q)
    assert p != target and p != (target[0], target[0] ^ target[1])
    for pt in (p, target):
        assert curve.scalar(pt, q) is None
        assert curve.tau(pt) == curve.scalar(pt, lam)
        assert curve.tau(pt) != pt  # n=131 prime, so orbit length is 131.
    assert pow(lam, 131, q) == 1 and lam != 1
    four_q = coordinates(fixture["four_Q"], curve.scalar(target, 4))

    torsion = [None, (0, 1), (1, 0), (1, 1)]
    assert len(fixture["Q_torsion_translates"]) == 4
    for i, record in enumerate(fixture["Q_torsion_translates"]):
        t = torsion[i]
        assert record["torsion_index"] == i and point(record["torsion"]) == t
        assert curve.on(t) and curve.scalar(t, 4) is None
        translated = coordinates(record["Q_plus_T"], curve.add(target, t))
        assert curve.scalar(translated, 4) == four_q
        assert curve.tau(translated) == curve.add(curve.tau(target), t)
    return {"status": "PASS", "field_model_identity": True,
            "normal_basis_rank": 131, "public_points_verified": 2,
            "torsion_translates_verified": 4,
            "all_group_law_and_normal_coordinate_checks": True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixture", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    assert not args.report.exists()
    args.report.write_text(json.dumps(verify(args.fixture),
                                      sort_keys=True, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
