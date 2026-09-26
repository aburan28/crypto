#!/usr/bin/env python3
"""Make a literal, directly encoded public ECC2K-130 P/Q input fixture."""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / "rotated_subspace_support_20260925/gate.py"
P = (0x051C99BFA6F18DE467C80C23B98C7994AA,
     0x042EA2D112ECEC71FCF7E000D7EFC978BD)
Q = (0x06C997F3E7F2C66A4A5D2FDA13756A37B1,
     0x04A38D11829D32D347BD0C0F584D546E9A)
POLY = (1 << 131) | (1 << 13) | (1 << 2) | (1 << 1) | 1
SOURCE_RELATIONS_BLOB = "a1df88fbf76fd14e565a5ec6b7d0ea4402f954ed"
SOURCE_FASTFIELD_BLOB = "bb57818dd85db7bff049e6f590d818616c5164e6"


def prior():
    spec = importlib.util.spec_from_file_location("challenge_import_prior", PRIOR)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def masks_for_basis(vectors):
    """Leftmost-high-bit F2 elimination, carrying normal coordinate masks."""
    pivots = {}
    for i, item in enumerate(vectors):
        value, mask = item, 1 << i
        while value:
            bit = value.bit_length() - 1
            if bit in pivots:
                old_value, old_mask = pivots[bit]
                value ^= old_value
                mask ^= old_mask
            else:
                pivots[bit] = (value, mask)
                break
        assert value, "normal conjugates are dependent"
    assert len(pivots) == 131

    def encode(value):
        mask = 0
        while value:
            bit = value.bit_length() - 1
            vector, coefficient = pivots[bit]
            value ^= vector
            mask ^= coefficient
        return mask

    return encode


def hex_point(p):
    return None if p is None else [hex(p[0]), hex(p[1])]


def make():
    mod = prior()
    field = mod.Field(131, mod.MODELS[131]["low"])
    curve = mod.Curve(field)
    assert field.poly == POLY and mod.Q131 == 680564733841876926932320129493409985129
    assert mod.source_group_order(131) == 4 * mod.Q131
    assert mod.LAMBDA131 == 196511074115861092422032515080945363956
    assert pow(mod.LAMBDA131, 131, mod.Q131) == 1
    assert mod.LAMBDA131 != 1
    for point in (P, Q):
        assert all(0 <= c <= field.mask for c in point)
        assert curve.on_curve(point) and curve.scalar(point, mod.Q131) is None
        assert curve.tau(point) == curve.scalar(point, mod.LAMBDA131)
    assert P != Q and P != curve.neg(Q)
    conjugates = mod.normal_conjugates(field, 3)
    encode = masks_for_basis(conjugates)

    def coordinates(p):
        if p is None:
            return None
        mask_x, mask_y = encode(p[0]), encode(p[1])
        assert mod.x_from_mask(conjugates, mask_x) == p[0]
        assert mod.x_from_mask(conjugates, mask_y) == p[1]
        return {"polynomial": hex_point(p), "normal_beta3": [hex(mask_x), hex(mask_y)],
                "trace_x": field.trace(p[0])}

    torsion = [None, (0, 1), (1, 0), (1, 1)]
    assert all(curve.on_curve(t) and curve.scalar(t, 4) is None for t in torsion)
    assert len(set(torsion)) == 4
    four_q = curve.scalar(Q, 4)
    translates = []
    for i, t in enumerate(torsion):
        target = curve.add(Q, t)
        assert target is not None and curve.scalar(target, 4) == four_q
        assert curve.tau(target) == curve.add(curve.tau(Q), t)
        translates.append({"torsion_index": i, "torsion": hex_point(t),
                           "Q_plus_T": coordinates(target)})
    return {
        "schema": "ecc2k130_direct_challenge_point_import_v1",
        "source": {"relations_path": "research/ecc2k130_relations/relations.py",
                   "relations_git_blob": SOURCE_RELATIONS_BLOB,
                   "fastfield_path": "research/ecc2k130_relations/fastfield.py",
                   "fastfield_git_blob": SOURCE_FASTFIELD_BLOB},
        "field": {"degree": 131, "polynomial": hex(POLY),
                  "encoding": "polynomial coefficient bit i is z^i",
                  "normal_beta": 3},
        "curve": {"equation": "y^2+x*y=x^3+1", "cofactor": 4,
                  "order": str(4 * mod.Q131), "subgroup_order": str(mod.Q131),
                  "frobenius_scalar": str(mod.LAMBDA131)},
        "P": coordinates(P), "Q": coordinates(Q),
        "four_Q": coordinates(four_q),
        "Q_torsion_translates": translates,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    assert not args.out.exists()
    args.out.write_text(json.dumps(make(), sort_keys=True, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
