"""Exact two-fiber character cells for 141 line-derived quotient classes."""

import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, Integer, prod

import type_i_divisor_transfer as dt

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
source_path = HERE / "type_i_line_guided_collector.json"
source = json.loads(source_path.read_text())
assert source["status"] == "verified_line_base_known_scalar_query_law_comparison"
assert source["base_classes"] == 141
p, r, L = dt.p, dt.r, dt.L
E = EllipticCurve(L, [L(dt.gm.a), L(dt.gm.b)])
O = E(0)
inv2 = Integer(pow(2, -1, int(r)))


def fiber_character(P):
    if P == O or P[0] == dt.c0:
        return "exceptional"
    X = dt.inf + dt.d0 / (P[0] - dt.c0)
    f = prod(X - root ** p for root in dt.quartic_roots)
    assert f ** (p ** 3) == f
    if f == 0:
        return "exceptional"
    return "+" if f ** ((p ** 3 - 1) // 2) == 1 else "-"


cells = {"++": 0, "+-": 0, "-+": 0, "--": 0, "exceptional": 0}
rows = []
for i, (x, y) in enumerate(source["base_points_fp"]):
    Q = E(L(Integer(x)), L(Integer(y)))
    assert Q != O and r * Q == O
    H = inv2 * Q
    assert 2 * H == Q
    first, second = fiber_character(Q), fiber_character(H)
    cell = first + second if first != "exceptional" and second != "exceptional" else "exceptional"
    cells[cell] += 1
    rows.append({"index": i, "cell": cell})
assert sum(cells.values()) == 141
out = {
    "status": "verified_line_quotient_two_fiber_character_census",
    "field_prime": str(p), "subgroup_order": str(r),
    "input_classes": 141, "cells": cells, "rows": rows,
    "fiber_rule": "Q and unique r-torsion half H=2^-1 Q; quartic RHS at inverse x is nonzero square in F_(p^3) iff character +",
    "source_sha256": sha(__file__),
    "input_receipt_sha256": sha(source_path),
    "runtime_sha256": sha(HERE / "type_i_line_fiber_runtime_info.json"),
}
(HERE / "type_i_line_fiber_eligibility.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("LINE_QUOTIENT_FIBER_PASS", json.dumps(cells, sort_keys=True))
