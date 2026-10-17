"""Replay every split-line quotient equation and center its common constant.

Run with the checked repository Sage launcher after saving --runtime-info.
"""

import hashlib
import json
import random
import time
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
assert prior["base_classes"] == 141
assert projection["new_projected_distinct"] == 9784

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


def coords(P):
    assert P != O
    return [int(P[0].polynomial()[0]), int(P[1].polynomial()[0])]


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
anchor_images = [elliptic_quotient_image(dt, E, anchor, j)
                 for j in range(3)]
old = [E(L(Integer(x)), L(Integer(y))) for x, y in prior["base_points_fp"]]
assert len(old) == len({key(Q) for Q in old}) == 141
assert hashlib.sha256("\n".join(key(Q) for Q in old).encode()).hexdigest() \
       == projection["base_key_sha256"]

new_rows = [tuple(map(int, row.split())) for row in
            (HERE / "type_i_line_new_points.tsv").read_text().splitlines()]
pair_rows = [tuple(map(int, row.split())) for row in
             (HERE / "type_i_line_new_point_pairs.tsv").read_text().splitlines()]
assert len(new_rows) == len(pair_rows) == 9784
assert all(pair_rows[2 * k] == pair_rows[2 * k + 1]
           for k in range(4892))
start = time.monotonic_ns()
new = []
new_keys = []
seen = set()
seen_signed = set()
old_signed = {min(key(Q), key(-Q)) for Q in old}
constant = None
center = None
for k in range(4892):
    projected = []
    for index in (2 * k, 2 * k + 1):
        u, v = new_rows[index]
        point = inverse.invert_affine(F(u), F(v))
        Q = h * (elliptic_quotient_image(dt, E, point, 0)
                 - anchor_images[0])
        assert Q != O and r * Q == O
        if index < 64:
            assert all(h * (elliptic_quotient_image(dt, E, point, j)
                            - anchor_images[j]) == Q for j in range(1, 3))
        current = key(Q)
        signed = min(current, key(-Q))
        assert current not in seen and signed not in seen_signed
        assert signed not in old_signed
        seen.add(current)
        seen_signed.add(signed)
        new_keys.append(current)
        new.append(coords(Q))
        projected.append(Q)
    i, j = pair_rows[2 * k]
    assert 1 <= i < j <= 141
    total = old[i - 1] + old[j - 1] + projected[0] + projected[1]
    if constant is None:
        constant = total
        assert constant != O and r * constant == O
        center = Integer(pow(4, -1, int(r))) * constant
        assert 4 * center == constant
    assert total == constant
    assert (old[i - 1] - center) + (old[j - 1] - center) \
           + (projected[0] - center) + (projected[1] - center) == O
    if (k + 1) % 500 == 0:
        print(f"LINE_FULL_CONSTANT_PROGRESS {k + 1}/4892", flush=True)

assert len(new) == len(seen) == len(seen_signed) == 9784
assert hashlib.sha256("\n".join(new_keys).encode()).hexdigest() \
       == projection["new_key_sha256"]
new_path = HERE / "type_i_line_full_new_quotient_points.json"
new_path.write_text(json.dumps(new, separators=(",", ":")) + "\n")
out = {
    "status": "verified_all_4892_split_line_quotient_constants",
    "field_prime": str(p), "subgroup_order": str(r),
    "old_points": len(old), "new_points": len(new), "line_cuts": 4892,
    "constant_fp": coords(constant), "center_fp": coords(center),
    "constant_is_identity": False, "center_fourfold_is_constant": True,
    "all_centered_line_sums_zero": True,
    "new_key_sha256": projection["new_key_sha256"],
    "new_coordinate_sha256": sha(new_path),
    "source_sha256": sha(__file__),
    "input_sha256": {name: sha(HERE / name) for name in
                     ("type_i_line_new_points.tsv",
                      "type_i_line_new_point_pairs.tsv")},
    "prior_sha256": sha(prior_path),
    "projection_sha256": sha(projection_path),
    "runtime_sha256": sha(HERE / "type_i_line_full_constant_runtime_info.json"),
    "elapsed_ns_unisolated": time.monotonic_ns() - start,
    "boundary": "Complete first-quotient replay of the 4892 frozen split lines; first 64 new points checked across all three quotient components. The centered coordinates do not themselves provide logs or a closed relation core.",
}
(HERE / "type_i_line_full_constant.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("LINE_FULL_CONSTANT_PASS", json.dumps({
    "line_cuts": 4892, "new_points": len(new),
    "constant_fp": coords(constant), "center_fp": coords(center)
}, sort_keys=True), flush=True)
