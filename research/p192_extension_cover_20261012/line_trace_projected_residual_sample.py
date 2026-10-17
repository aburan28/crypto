"""Exact old-anchored residual columns using trace and the line-section sum.

Enumerate the first N split non-source cuts in canonical pivot/neighbor
order. Project one residual per cut by Frobenius trace; obtain its partner
from the independently replayed constant line-section quotient sum. Direct
projection of the first 128 partners checks the second-column shortcut.
"""

import argparse
import hashlib
import json
import random
import time
from collections import Counter
from math import comb
from pathlib import Path

from sage.all import Integer, prod, vector

import type_i_divisor_transfer as dt
from frobenius_trace_projection import FrobeniusTraceProjector
from plane_point_inverse_cached import CachedPlanePointInverse

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--split-lines", type=int, default=4096)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
assert 1 <= args.split_lines <= 680474
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
inputs = {
    "plane": HERE / "type_i_plane_quartic.json",
    "old_quartic": HERE / "type_i_line_base_points.tsv",
    "new_quartic": HERE / "type_i_line_new_points.tsv",
    "collinearity": HERE / "type_i_line_full_collinearity.json",
    "old_quotient": HERE / "type_i_line_guided_collector.json",
    "new_quotient": HERE / "type_i_line_full_new_quotient_points.json",
    "constant": HERE / "type_i_line_full_constant.json",
    "base": HERE / "type_i_quotient_base_probe.json",
    "full_screen": HERE / "type_i_line_nonsource_reuse_old141.json",
}
plane = json.loads(inputs["plane"].read_text())
collinearity = json.loads(inputs["collinearity"].read_text())
old = [tuple(map(int, line.split()[1:])) for line in
       inputs["old_quartic"].read_text().splitlines()]
new = [tuple(map(int, line.split())) for line in
       inputs["new_quartic"].read_text().splitlines()]
points = old + new
assert len(old) == 141 and len(points) == len(set(points)) == 9925
old_q = json.loads(inputs["old_quotient"].read_text())["base_points_fp"]
new_q = json.loads(inputs["new_quotient"].read_text())
assert len(old_q) == 141 and len(new_q) == 9784
full_screen = json.loads(inputs["full_screen"].read_text())
assert full_screen["old_pivots"] == 141
source_neighbors = [set() for _ in range(141)]
for row in collinearity["line_records"]:
    indices = list(map(int, row["point_indices_zero_based"]))
    for i in indices:
        if i < 141:
            source_neighbors[i].update(j for j in indices if j > i)

p = int(dt.p)
assert p == int(plane["p"]) and p % 4 == 3
terms = [(tuple(map(int, label.split(","))), int(value) % p)
         for label, value in plane["quartic_coefficients"].items()]


def taylor_at(x, y):
    xp, yp = [pow(x, i, p) for i in range(5)], [pow(y, i, p) for i in range(5)]
    out = [[0] * (k + 1) for k in range(5)]
    for (a, b, _), value in terms:
        for ix in range(a + 1):
            for iy in range(b + 1):
                k = ix + iy
                out[k][iy] = (out[k][iy] + value * comb(a, ix) * comb(b, iy)
                               * xp[a - ix] * yp[b - iy]) % p
    assert out[0][0] == 0
    return out


def evaluate(coefficients, z):
    value = 0
    for c in reversed(coefficients):
        value = (value * z + c) % p
    return value


def cover_from_x(value):
    X = dt.L(value)
    u = (X - dt.zeta) * (X - dt.beta)
    y = u ** ((dt.q + 1) // 2)
    norm = prod(X - root ** (dt.p ** j)
                for root in (dt.zeta, dt.beta) for j in range(3))
    if y * y ** dt.p * y ** (dt.p * dt.p) == -norm:
        y = -y
    assert y * y ** dt.p * y ** (dt.p * dt.p) == norm
    Y = dt.Vandermonde.solve_right(vector(dt.L, [y, y ** dt.p, y ** (dt.p * dt.p)]))
    point = (X, *Y)
    dt.check_model(point)
    return point


base = json.loads(inputs["base"].read_text())
h = Integer(base["candidate_cofactor"])
projector = FrobeniusTraceProjector(dt.gm.F, dt.L, dt.gm.a, dt.gm.b, dt.r, h)
assert h * dt.r == Integer(base["extension_elliptic_order"])
inverse = CachedPlanePointInverse(dt)
rng = random.Random(base["pilot_seed"])
anchor_x = [Integer(rng.randrange(p))
            for _ in range(base["anchor_pilot_index"] + 1)][-1]
anchor_point = cover_from_x(anchor_x)
anchor_y = dt.check_model(anchor_point)[0]
anchor = projector.E6(dt.c0 + dt.d0 / (anchor_point[0] - dt.inf),
                      dt.kappa * anchor_y / (anchor_point[0] - dt.inf) ** 2)
anchor_trace = projector.trace(anchor)
E = projector.E0
selected_q = [E(dt.gm.F(int(x)), dt.gm.F(int(y))) for x, y in old_q + new_q]
selected_keys = {projector.key(Q) for Q in selected_q}
constant_data = json.loads(inputs["constant"].read_text())
S = E(*[dt.gm.F(int(x)) for x in constant_data["constant_fp"]])
assert len(selected_keys) == 9925 and dt.r * S == E(0)
assert constant_data["new_coordinate_sha256"] == sha(inputs["new_quotient"])

rows = []
seen = {}
signed_seen = {}
selected_hits = []
repeats = []
signed_repeats = []
first128 = full_screen["first_128_split_cases"]
counts = Counter()
projection_ns = 0
start_all = time.monotonic_ns()
for i in range(141):
    if len(rows) == args.split_lines:
        break
    x, y = points[i]
    taylor = taylor_at(x, y)
    for j in range(i + 1, len(points)):
        if j in source_neighbors[i]:
            continue
        u, v = points[j]
        d = (u - x) % p
        m = ((v - y) * pow(d, -1, p)) % p
        c1, c2, c3, A = (evaluate(taylor[k], m) for k in range(1, 5))
        B = (c3 + d * A) % p
        C = (c2 + d * B) % p
        assert A and C and (c1 + d * C) % p == 0
        discr = (B * B - 4 * A * C) % p
        assert discr
        if pow(discr, (p - 1) // 2, p) != 1:
            continue
        sqrt = pow(discr, (p + 1) // 4, p)
        den = pow(2 * A, -1, p)
        residuals = [((x + root) % p, (y + m * root) % p)
                     for root in (((-B + sqrt) * den) % p,
                                  ((-B - sqrt) * den) % p)]
        assert residuals[0] != residuals[1]
        row_index = len(rows)
        if row_index < 128:
            frozen = first128[row_index]
            assert frozen["selected_pair"] == [i, j]
            assert residuals == [tuple(map(int, q)) for q in frozen["residual_points"]]
        start = time.monotonic_ns()
        cover, ys = inverse.invert_affine(*[dt.gm.F(z) for z in residuals[0]])
        image = inverse.quotient_image(projector.E6, cover, ys, 0)
        Q1 = projector.project_delta(image, anchor_trace)
        Q2 = S - selected_q[i] - selected_q[j] - Q1
        projection_ns += time.monotonic_ns() - start
        if row_index < 128 or row_index % 257 == 0:
            cover2, ys2 = inverse.invert_affine(*[dt.gm.F(z) for z in residuals[1]])
            image2 = inverse.quotient_image(projector.E6, cover2, ys2, 0)
            assert projector.key(projector.project_delta(image2, anchor_trace)) == projector.key(Q2)
            counts["direct_partner_controls"] += 1
        keys = [projector.key(Q1), projector.key(Q2)]
        for side, (Q, key) in enumerate(zip((Q1, Q2), keys)):
            counts["projected_columns"] += 1
            if key is None:
                counts["identity_columns"] += 1
            if key in selected_keys and len(selected_hits) < 100:
                selected_hits.append([row_index, side, key])
            if key in seen and len(repeats) < 100:
                repeats.append([row_index, side, seen[key]])
            seen.setdefault(key, [row_index, side])
            negative_key = projector.key(-Q)
            signed_key = min(key, negative_key) if key is not None else None
            if signed_key in signed_seen and len(signed_repeats) < 100:
                signed_repeats.append([row_index, side, signed_seen[signed_key]])
            signed_seen.setdefault(signed_key, [row_index, side])
        rows.append({
            "selected_pair": [i, j],
            "residual_points": [[str(z) for z in point] for point in residuals],
            "projected_residuals": [[str(z) for z in key] if key is not None else None
                                    for key in keys],
        })
        if len(rows) % 512 == 0:
            print("TRACE_RESIDUAL_PROGRESS", len(rows), args.split_lines,
                  "repeats", counts["projected_columns"] - len(seen), flush=True)
        if len(rows) == args.split_lines:
            break

assert len(rows) == args.split_lines
assert counts["projected_columns"] == 2 * args.split_lines
out = {
    "status": "verified_p192_old_anchor_trace_projected_residual_sample",
    "split_lines": len(rows),
    "projected_columns": counts["projected_columns"],
    "distinct_projected_columns": len(seen),
    "distinct_signed_projected_columns": len(signed_seen),
    "selected_column_hits_first_100": selected_hits,
    "repeated_columns_first_100": repeats,
    "signed_repeated_columns_first_100": signed_repeats,
    "identity_columns": counts["identity_columns"],
    "direct_partner_controls": counts["direct_partner_controls"],
    "projection_ns_unisolated": projection_ns,
    "elapsed_ns_unisolated": time.monotonic_ns() - start_all,
    "rows": rows,
    "source_sha256": sha(__file__),
    "projector_source_sha256": sha(HERE / "frobenius_trace_projection.py"),
    "inverse_source_sha256": sha(HERE / "plane_point_inverse_cached.py"),
    "input_sha256": {name: sha(path) for name, path in inputs.items()},
    "runtime_info_sha256": sha(HERE / "type_i_line_trace_projection_runtime_info.json"),
    "boundary": "First N rational old-anchored non-source cuts in canonical pivot order; first residual directly projected, second inferred from verified common line-section constant with direct controls. Local costs are unisolated diagnostics.",
}
args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print("TRACE_RESIDUAL_PASS", json.dumps({
    "split_lines": len(rows),
    "projected_columns": counts["projected_columns"],
    "distinct_columns": len(seen),
    "selected_hits": len(selected_hits),
    "repeated_columns": counts["projected_columns"] - len(seen),
    "direct_partner_controls": counts["direct_partner_controls"],
}, sort_keys=True))
