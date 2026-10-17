"""Exact projected line fan from the first non-source residual quartic point.

This is a second-generation pilot.  It retains the full selected point set,
derives every rational split cut through one frozen residual pivot, projects
its quotient images by three-step Frobenius trace, and streams the complete
old-anchored census to test subgroup-column reuse.  Sparse six-step controls
check the optimized projection independently.
"""

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
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
input_names = (
    "type_i_plane_quartic.json",
    "type_i_line_base_points.tsv",
    "type_i_line_new_points.tsv",
    "type_i_line_guided_collector.json",
    "type_i_line_full_new_quotient_points.json",
    "type_i_line_full_constant.json",
    "type_i_quotient_base_probe.json",
    "type_i_line_trace_full_chunks/pivot-000.json",
    "type_i_line_trace_full_verified_20261016T1849Z.json",
)
inputs = {name: sha(HERE / name) for name in input_names}
runtime = HERE / "type_i_residual_pivot000_three_runtime_info.json"
assert runtime.exists()
plane = json.loads((HERE / input_names[0]).read_text())
p = int(plane["p"])
assert p == int(dt.p) and p % 4 == 3
old_points = [tuple(map(int, row.split()[1:])) for row in
              (HERE / input_names[1]).read_text().splitlines()]
new_points = [tuple(map(int, row.split())) for row in
              (HERE / input_names[2]).read_text().splitlines()]
points = old_points + new_points
assert len(old_points) == 141 and len(points) == len(set(points)) == 9925
first = json.loads((HERE / input_names[7]).read_text())["rows"][0]
assert first["selected_pair"] == [0, 285]
pivot = tuple(map(int, first["residual_points"][0]))
partner = tuple(map(int, first["residual_points"][1]))
assert pivot not in set(points) and partner not in set(points)

terms = [(tuple(map(int, label.split(","))), int(value) % p)
         for label, value in plane["quartic_coefficients"].items()]


def taylor_at(x, y):
    xp, yp = [pow(x, i, p) for i in range(5)], [pow(y, i, p) for i in range(5)]
    out = [[0] * (k + 1) for k in range(5)]
    for (a, b, _), value in terms:
        for ix in range(a + 1):
            for iy in range(b + 1):
                out[ix + iy][iy] = (out[ix + iy][iy]
                                    + value * comb(a, ix) * comb(b, iy)
                                    * xp[a - ix] * yp[b - iy]) % p
    assert out[0][0] == 0
    return out


def evaluate(coefficients, z):
    value = 0
    for c in reversed(coefficients):
        value = (value * z + c) % p
    return value


base = json.loads((HERE / input_names[6]).read_text())
projector = FrobeniusTraceProjector(
    dt.gm.F, dt.L, dt.gm.a, dt.gm.b, dt.r,
    Integer(base["candidate_cofactor"]))
E6, E0 = projector.E6, projector.E0
assert not E0.a1() and not E0.a3()
inverse = CachedPlanePointInverse(dt)


def trace3(point):
    P = E6(point)
    if not P:
        return E0(0)
    assert P[0] ** (p ** 3) == P[0] and P[1] ** (p ** 3) == P[1]
    total = E6(0)
    x, y = P[0], P[1]
    for _ in range(3):
        total += E6(x, y)
        x, y = x ** p, y ** p
    if not total:
        return E0(0)
    assert total[0] ** p == total[0] and total[1] ** p == total[1]
    return E0(dt.gm.F(total[0]), dt.gm.F(total[1]))


def cover_from_x(value):
    X = dt.L(value)
    y = ((X - dt.zeta) * (X - dt.beta)) ** ((dt.q + 1) // 2)
    norm = prod(X - root ** (dt.p ** j)
                for root in (dt.zeta, dt.beta) for j in range(3))
    if y * y ** dt.p * y ** (dt.p * dt.p) == -norm:
        y = -y
    assert y * y ** dt.p * y ** (dt.p * dt.p) == norm
    Y = dt.Vandermonde.solve_right(
        vector(dt.L, [y, y ** dt.p, y ** (dt.p * dt.p)]))
    cover = (X, *Y)
    dt.check_model(cover)
    return cover


rng = random.Random(base["pilot_seed"])
anchor_x = [Integer(rng.randrange(p))
            for _ in range(base["anchor_pilot_index"] + 1)][-1]
anchor_cover = cover_from_x(anchor_x)
anchor_y = dt.check_model(anchor_cover)[0]
anchor = E6(dt.c0 + dt.d0 / (anchor_cover[0] - dt.inf),
            dt.kappa * anchor_y / (anchor_cover[0] - dt.inf) ** 2)
anchor_trace3 = trace3(anchor)
assert projector.trace(anchor) == 2 * anchor_trace3
scalar3 = Integer(2 * projector.scalar % dt.r)
selected_old = json.loads((HERE / input_names[3]).read_text())["base_points_fp"]
selected_new = json.loads((HERE / input_names[4]).read_text())
selected_q = [E0(dt.gm.F(int(x)), dt.gm.F(int(y)))
              for x, y in selected_old + selected_new]
assert len(selected_q) == 9925
pivot_q = E0(*[dt.gm.F(int(z)) for z in first["projected_residuals"][0]])
partner_q = E0(*[dt.gm.F(int(z)) for z in first["projected_residuals"][1]])
constant = json.loads((HERE / input_names[5]).read_text())
S = E0(*[dt.gm.F(int(z)) for z in constant["constant_fp"]])
assert selected_q[0] + selected_q[285] + pivot_q + partner_q == S


def signed_key(key):
    if key is None:
        return None
    x, y = key
    return min(key, (x, -y % p))


selected_keys = {projector.key(Q) for Q in selected_q}
selected_signed = {signed_key(k) for k in selected_keys}
selected_geometry = set(points)
seen_geometry, seen_exact, seen_signed = {}, {}, {}
rows, counts = [], Counter()
sample_controls = []
taylor = taylor_at(*pivot)
start_ns = time.monotonic_ns()
projection_ns = 0
for j, (u, v) in enumerate(points):
    counts["selected_pairs"] += 1
    d = (u - pivot[0]) % p
    if not d:
        counts["vertical"] += 1
        continue
    m = (v - pivot[1]) * pow(d, -1, p) % p
    c1, c2, c3, A = (evaluate(taylor[k], m) for k in range(1, 5))
    if not A:
        counts["degree_below_four"] += 1
        continue
    B, C = (c3 + d * A) % p, (c2 + d * (c3 + d * A)) % p
    if not C or (c1 + d * C) % p:
        counts["tangent_or_invalid"] += 1
        continue
    discr = (B * B - 4 * A * C) % p
    if not discr:
        counts["double_residual"] += 1
        continue
    if pow(discr, (p - 1) // 2, p) == p - 1:
        counts["nonsplit"] += 1
        continue
    counts["split"] += 1
    sq, den = pow(discr, (p + 1) // 4, p), pow(2 * A, -1, p)
    roots = ((-B + sq) * den % p, (-B - sq) * den % p)
    residuals = [((pivot[0] + z) % p, (pivot[1] + m * z) % p)
                 for z in roots]
    assert residuals[0] != residuals[1]
    for side, point in enumerate(residuals):
        if point in selected_geometry:
            counts["geometric_selected_hits"] += 1
        if point in seen_geometry:
            counts["geometric_repeats"] += 1
        seen_geometry.setdefault(point, [j, side])
    start = time.monotonic_ns()
    cover, ys = inverse.invert_affine(*[dt.gm.F(z) for z in residuals[0]])
    image = inverse.quotient_image(E6, cover, ys, 0)
    Q1 = scalar3 * (trace3(image) - anchor_trace3)
    Q2 = S - pivot_q - selected_q[j] - Q1
    projection_ns += time.monotonic_ns() - start
    assert pivot_q + selected_q[j] + Q1 + Q2 == S
    if len(rows) % 256 == 0:
        assert projector.project_delta(image, projector.trace(anchor)) == Q1
        cover2, ys2 = inverse.invert_affine(*[dt.gm.F(z) for z in residuals[1]])
        image2 = inverse.quotient_image(E6, cover2, ys2, 0)
        assert projector.project_delta(image2, projector.trace(anchor)) == Q2
        sample_controls.append(len(rows))
    keys = [projector.key(Q1), projector.key(Q2)]
    for side, key in enumerate(keys):
        skey = signed_key(key)
        if key is None:
            counts["identity_columns"] += 1
        if key in selected_keys:
            counts["selected_hits"] += 1
        if skey in selected_signed:
            counts["signed_selected_hits"] += 1
        if key in seen_exact:
            counts["exact_repeats"] += 1
        if skey in seen_signed:
            counts["signed_repeats"] += 1
        seen_exact.setdefault(key, [len(rows), side])
        seen_signed.setdefault(skey, [len(rows), side])
    rows.append({"selected_index": j,
                 "residual_points": [list(map(str, point)) for point in residuals],
                 "projected_residuals": [list(map(str, key)) if key is not None
                                         else None for key in keys]})
    if len(rows) % 1024 == 0:
        print("RESIDUAL_PIVOT_THREE_PROGRESS", len(rows), j, flush=True)

assert len(rows) == 5037 and len(seen_geometry) == 10073
assert len(sample_controls) == 20
old_exact_hits, old_signed_hits = [], []
old_columns = 0
for i in range(141):
    chunk = json.loads((HERE / f"type_i_line_trace_full_chunks/pivot-{i:03d}.json").read_text())
    for old_row_index, old_row in enumerate(chunk["rows"]):
        for side, pair in enumerate(old_row["projected_residuals"]):
            key = tuple(map(int, pair))
            old_columns += 1
            if key in seen_exact:
                old_exact_hits.append([i, old_row_index, side, seen_exact[key]])
            if signed_key(key) in seen_signed:
                old_signed_hits.append([i, old_row_index, side,
                                        seen_signed[signed_key(key)]])
assert old_columns == 1360948

receipt = {
    "status": "p192_residual_pivot000_three_step_projected_fan",
    "boundary": "One frozen residual pivot through every selected quartic point; full old-anchored subgroup-column membership scan. Six-step direct controls at every 256th split row; other second residuals use the checked line-section equation.",
    "field_prime": str(p), "subgroup_order": str(dt.r),
    "source_line_selected_pair": first["selected_pair"],
    "pivot_quartic_point": list(map(str, pivot)),
    "pivot_quotient_point": list(map(str, projector.key(pivot_q))),
    "counts": dict(counts),
    "distinct_geometry": len(seen_geometry),
    "distinct_exact_columns": len(seen_exact),
    "distinct_signed_columns": len(seen_signed),
    "old_family_columns_scanned": old_columns,
    "old_exact_hits": old_exact_hits,
    "old_signed_hits": old_signed_hits,
    "six_step_control_row_indices": sample_controls,
    "projection_ns_unisolated": projection_ns,
    "elapsed_ns_unisolated": time.monotonic_ns() - start_ns,
    "source_sha256": sha(__file__),
    "input_sha256": inputs,
    "old_family_aggregate_sha256": inputs[input_names[8]],
    "runtime_info_sha256": sha(runtime),
    "projector_source_sha256": sha(HERE / "frobenius_trace_projection.py"),
    "inverse_source_sha256": sha(HERE / "plane_point_inverse_cached.py"),
    "rows": rows,
}
out = HERE / "type_i_residual_pivot000_three.json"
out.write_text(json.dumps(receipt, sort_keys=True, separators=(",", ":")) + "\n")
print("RESIDUAL_PIVOT_THREE_DONE", dict(counts), len(old_exact_hits),
      len(old_signed_hits), flush=True)
