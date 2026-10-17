"""Project the frozen rational line-support points into P-192 r-torsion.

Run with the repository's checked Sage launcher. Elapsed values are local
correctness-run diagnostics, not isolated wall-time comparisons.
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
base_receipt = json.loads((HERE / "type_i_quotient_base_probe.json").read_bytes())
plane_receipt = json.loads((HERE / "type_i_plane_quartic.json").read_bytes())
fixture = json.loads((HERE / "type_i_base_points_legacy_fixture.json").read_bytes())
gp_log = (HERE / "type_i_line_support_screen.log").read_text().splitlines()
assert gp_log[0].startswith("LINE_RESIDUAL_GP_PASS points=141 pairs=9870 ")
assert "nonsplit_split=[4978, 4892]" in gp_log[0]
assert "new_rational_points=9784 all_points_distinct=9925" in gp_log[0]
assert len(gp_log) == 2 and not any("***" in line for line in gp_log)
assert base_receipt["p"] == str(dt.p) and base_receipt["r"] == str(dt.r)
assert fixture["field_prime"] == str(dt.p) and len(fixture["points"]) == 16

F, L, p = dt.gm.F, dt.L, dt.p
inverse = PlanePointInverse(dt)
E = EllipticCurve(L, [L(dt.gm.a), L(dt.gm.b)])
O = E(0)
h = Integer(base_receipt["candidate_cofactor"])
assert h * dt.r == Integer(base_receipt["extension_elliptic_order"])

def cover_from_x(x):
    X = L(x)
    U = (X - dt.zeta) * (X - dt.beta)
    y = U ** ((dt.q + 1) // 2)
    norm = prod(X - root ** (p ** j)
                for root in (dt.zeta, dt.beta) for j in range(3))
    if y * y ** p * y ** (p * p) == -norm:
        y = -y
    assert y * y ** p * y ** (p * p) == norm
    Y = dt.Vandermonde.solve_right(vector(L, [y, y ** p, y ** (p * p)]))
    point = (X, *Y)
    dt.check_model(point)
    return point


def images(point):
    return [elliptic_quotient_image(dt, E, point, j) for j in range(3)]


def point_key(Q):
    if Q == O:
        return "identity"
    assert Q[0] ** p == Q[0] and Q[1] ** p == Q[1]
    return f"{int(Q[0].polynomial()[0])}:{int(Q[1].polynomial()[0])}"


rng = random.Random(base_receipt["pilot_seed"])
anchor_index = base_receipt["anchor_pilot_index"]
anchor_x = [Integer(rng.randrange(int(p))) for _ in range(anchor_index + 1)][-1]
anchor_images = images(cover_from_x(anchor_x))
fixture_points = [E(L(int(row["x"])), L(int(row["y"])))
                  for row in fixture["points"]]
fixture_keys = {point_key(Q) for Q in fixture_points}
negative_fixture_keys = {point_key(-Q) for Q in fixture_points}

# Positive control: invert all 16 old plane points and recover their exact
# projected quotient coordinates, not merely their support abscissae.
controls = [row for row in plane_receipt["known_points_verified"]
            if row["source"] == "base"]
assert len(controls) == 16
for row in controls:
    index = row["index"]
    c = [F(int(value)) for value in row["coordinates_Fp_decimal"]]
    point = inverse.invert_affine(c[0] / c[2], c[1] / c[2])
    assert point[0] == L(F(int(base_receipt["base_selection"][index]["x"])))
    quotient = images(point)
    projected = [h * (quotient[j] - anchor_images[j]) for j in range(3)]
    assert projected[0] == projected[1] == projected[2] == fixture_points[index]

base_path = HERE / "type_i_line_base_points.tsv"
new_path = HERE / "type_i_line_new_points.tsv"
pair_path = HERE / "type_i_line_new_point_pairs.tsv"
base_rows = [tuple(map(int, line.split())) for line in base_path.read_text().splitlines()]
new_rows = [tuple(map(int, line.split())) for line in new_path.read_text().splitlines()]
pair_rows = [tuple(map(int, line.split())) for line in pair_path.read_text().splitlines()]
assert len(base_rows) == 141 and len(new_rows) == len(pair_rows) == 9784
assert all(1 <= index <= 160 and 0 <= u < p and 0 <= v < p
           for index, u, v in base_rows)
assert all(0 <= u < p and 0 <= v < p for u, v in new_rows)

base_start = time.monotonic_ns()
base_points = []
for _, u, v in base_rows:
    point = inverse.invert_affine(F(u), F(v))
    quotient = images(point)
    projected = [h * (quotient[j] - anchor_images[j]) for j in range(3)]
    assert projected[0] == projected[1] == projected[2]
    assert projected[0] != O and dt.r * projected[0] == O
    base_points.append(projected[0])
base_elapsed = time.monotonic_ns() - base_start
base_keys = [point_key(Q) for Q in base_points]
base_sign_keys = [min(point_key(Q), point_key(-Q)) for Q in base_points]
assert len(set(base_keys)) == len(set(base_sign_keys)) == 141
assert all(key not in fixture_keys | negative_fixture_keys for key in base_keys)

pair_start = time.monotonic_ns()
pair_keys = [point_key(base_points[i] + base_points[j])
             for i in range(141) for j in range(i + 1, 141)]
pair_elapsed = time.monotonic_ns() - pair_start
assert len(pair_keys) == len(set(pair_keys)) == 9870

new_start = time.monotonic_ns()
seen = set()
seen_sign = set()
base_sign_set = set(base_sign_keys)
new_keys = []
new_zero = old_matches = base_matches = new_repeats = sign_repeats = 0
for index, (u, v) in enumerate(new_rows):
    point = inverse.invert_affine(F(u), F(v))
    quotient0 = elliptic_quotient_image(dt, E, point, 0)
    Q = h * (quotient0 - anchor_images[0])
    assert dt.r * Q == O
    key = point_key(Q)
    signed = min(key, point_key(-Q))
    if index < 64:
        quotient = images(point)
        projected = [h * (quotient[j] - anchor_images[j]) for j in range(3)]
        assert projected[0] == projected[1] == projected[2] == Q
    new_zero += Q == O
    old_matches += key in fixture_keys | negative_fixture_keys
    base_matches += signed in base_sign_set
    new_repeats += key in seen
    sign_repeats += signed in seen_sign
    seen.add(key)
    seen_sign.add(signed)
    new_keys.append(key)
    if (index + 1) % 1000 == 0:
        print(f"LINE_PROJECTION_PROGRESS {index + 1}/{len(new_rows)}", flush=True)
new_elapsed = time.monotonic_ns() - new_start

out = {
    "status": "verified_line_support_first_quotient_projection_screen",
    "candidate_id": None,
    "field_prime": str(p), "subgroup_order": str(dt.r),
    "base_rational_points": 141, "new_line_rational_points": len(new_rows),
    "pair_lines": 9870, "distinct_base_pair_sums": len(set(pair_keys)),
    "base_projected_distinct": len(set(base_keys)),
    "new_projected_distinct": len(seen),
    "new_projected_sign_distinct": len(seen_sign),
    "new_zero": int(new_zero), "old16_signed_matches": int(old_matches),
    "base141_signed_matches": int(base_matches),
    "new_exact_repeats": int(new_repeats), "new_signed_repeats": int(sign_repeats),
    "old16_controls": 16, "new_three_component_controls": 64,
    "base_projection_ns_unisolated": base_elapsed,
    "pair_sum_ns_unisolated": pair_elapsed,
    "new_projection_ns_unisolated": new_elapsed,
    "input_sha256": {name: sha(path) for name, path in
                     (("base_points", base_path), ("new_points", new_path),
                      ("new_point_pairs", pair_path))},
    "source_sha256": sha(__file__),
    "inverse_source_sha256": sha(HERE / "plane_point_inverse.py"),
    "gp_source_sha256": sha(HERE / "line_support_screen.gp"),
    "gp_log_sha256": sha(HERE / "type_i_line_support_screen.log"),
    "runtime_sha256": sha(HERE / "type_i_line_sage_runtime_info.json"),
    "base_key_sha256": hashlib.sha256("\n".join(base_keys).encode()).hexdigest(),
    "new_key_sha256": hashlib.sha256("\n".join(new_keys).encode()).hexdigest(),
    "boundary": "First quotient for every new point; all 141 old and 64 new points checked in all three quotient coordinates. Local elapsed diagnostics are not isolated performance results. A full relation collector and target cost remain separate gates.",
}
(HERE / "type_i_line_projection.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print(json.dumps(out, sort_keys=True), flush=True)
