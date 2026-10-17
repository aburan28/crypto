"""Check split-line quotient identities on frozen P-192 Type I points.

This is a correctness probe over the first 32 split lines. Run through the
repository's checked Sage launcher after saving its runtime-info receipt.
"""

import hashlib
import json
import random
from pathlib import Path

from sage.all import EllipticCurve, Integer, prod, vector

import type_i_divisor_transfer as dt
from plane_point_inverse import PlanePointInverse, elliptic_quotient_image

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
prior_path = HERE / "type_i_line_guided_collector.json"
projection_path = HERE / "type_i_line_projection.json"
base_receipt_path = HERE / "type_i_quotient_base_probe.json"
prior = json.loads(prior_path.read_text())
projection = json.loads(projection_path.read_text())
base_receipt = json.loads(base_receipt_path.read_text())
assert len(prior["base_points_fp"]) == 141
assert projection["new_projected_distinct"] == 9784
assert projection["base141_signed_matches"] == 0

p, r, L = dt.p, dt.r, dt.L
F = dt.gm.F
E = EllipticCurve(L, [L(dt.gm.a), L(dt.gm.b)])
O = E(0)
h = Integer(base_receipt["candidate_cofactor"])
assert h * r == Integer(base_receipt["extension_elliptic_order"])


def key(P):
    if P == O:
        return "identity"
    assert P[0] ** p == P[0] and P[1] ** p == P[1]
    return f"{int(P[0].polynomial()[0])}:{int(P[1].polynomial()[0])}"


def cover_from_x(x):
    X = L(x)
    U = (X - dt.zeta) * (X - dt.beta)
    y = U ** ((dt.q + 1) // 2)
    norm = prod(X - root ** (p ** j)
                for root in (dt.zeta, dt.beta) for j in range(3))
    if y * y ** p * y ** (p * p) == -norm:
        y = -y
    assert y * y ** p * y ** (p * p) == norm
    return (X, *dt.Vandermonde.solve_right(
        vector(L, [y, y ** p, y ** (p * p)])))


rng = random.Random(base_receipt["pilot_seed"])
anchor_x = [Integer(rng.randrange(int(p)))
            for _ in range(base_receipt["anchor_pilot_index"] + 1)][-1]
anchor = cover_from_x(anchor_x)
inverse = PlanePointInverse(dt)
anchor_image = elliptic_quotient_image(dt, E, anchor, 0)


def project(u, v):
    point = inverse.invert_affine(F(u), F(v))
    Q = h * (elliptic_quotient_image(dt, E, point, 0) - anchor_image)
    assert Q != O and r * Q == O
    return Q


base = [E(L(Integer(x)), L(Integer(y))) for x, y in prior["base_points_fp"]]
assert len({key(Q) for Q in base}) == 141
new_rows = [tuple(map(int, line.split())) for line in
            (HERE / "type_i_line_new_points.tsv").read_text().splitlines()]
pair_rows = [tuple(map(int, line.split())) for line in
             (HERE / "type_i_line_new_point_pairs.tsv").read_text().splitlines()]
assert len(new_rows) == len(pair_rows) == 9784
assert all(pair_rows[2 * k] == pair_rows[2 * k + 1] for k in range(4892))

line_rows = []
constant = None
fresh_keys = set()
for k in range(32):
    i, j = pair_rows[2 * k]
    A = project(*new_rows[2 * k])
    B = project(*new_rows[2 * k + 1])
    assert key(A) != key(B)
    fresh_keys.update((key(A), key(B)))
    total = base[i - 1] + base[j - 1] + A + B
    if constant is None:
        constant = total
    assert total == constant
    line_rows.append({"base_indices_1based": [i, j],
                      "fresh_points_fp": [[int(P[0].polynomial()[0]),
                                           int(P[1].polynomial()[0])]
                                          for P in (A, B)],
                      "constant_key": key(total)})
assert len(fresh_keys) == 64
out = {
    "status": "verified_first32_split_line_quotient_constant",
    "line_count": 32,
    "fresh_point_count": len(fresh_keys),
    "full_split_line_count": 4892,
    "full_fresh_point_count": 9784,
    "base_point_count": 141,
    "constant_key": key(constant),
    "constant_is_identity": constant == O,
    "line_rows": line_rows,
    "input_sha256": {name: sha(HERE / name) for name in
                     ("type_i_line_new_points.tsv",
                      "type_i_line_new_point_pairs.tsv")},
    "prior_sha256": sha(prior_path),
    "projection_sha256": sha(projection_path),
    "source_sha256": sha(__file__),
    "runtime_sha256": sha(HERE / "type_i_line_leaf_runtime_info.json"),
    "boundary": "First 32 of 4892 split lines checked under Sage group law; full line incidence and distinctness come from frozen GP and Sage projection receipts.",
}
(HERE / "type_i_line_leaf_relation_probe.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("LINE_LEAF_RELATION_PASS", json.dumps({
    "line_count": out["line_count"], "fresh_point_count": len(fresh_keys),
    "constant_is_identity": constant == O, "constant_key": key(constant)
}, sort_keys=True), flush=True)
