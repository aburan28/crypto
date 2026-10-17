"""Replay a three-step Frobenius trace on frozen P-192 residual columns.

The inverse formulas use only constants fixed by p^3-Frobenius and rational
plane inputs. A resulting quotient image is therefore defined over F_(p^3).
This script checks the constants, image field membership, both trace laws,
and the frozen subgroup column on both residuals of 128 exact quartic cuts.
"""

import hashlib
import json
import random
import time
from pathlib import Path

from sage.all import Integer, prod, vector

import type_i_divisor_transfer as dt
from frobenius_trace_projection import FrobeniusTraceProjector
from plane_point_inverse_cached import CachedPlanePointInverse


HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
runtime_path = HERE / "type_i_three_trace_runtime_info.json"
sample_path = HERE / "type_i_line_trace_residual_4096.json"
base_path = HERE / "type_i_quotient_base_probe.json"
base = json.loads(base_path.read_text())
sample = json.loads(sample_path.read_text())["rows"][:128]
assert len(sample) == 128
p = int(dt.p)
projector = FrobeniusTraceProjector(
    dt.gm.F, dt.L, dt.gm.a, dt.gm.b, dt.r,
    Integer(base["candidate_cofactor"]))
inverse = CachedPlanePointInverse(dt)
E6, E0 = projector.E6, projector.E0

constants = {
    "zeta": dt.zeta, "beta": dt.beta, "inf": dt.inf,
    "c0": dt.c0, "d0": dt.d0, "kappa": dt.kappa,
}
constants.update({f"quartic_root_{i}": x
                  for i, x in enumerate(dt.quartic_roots)})
constant_fixed = {name: bool(x ** (p ** 3) == x)
                  for name, x in constants.items()}
assert all(constant_fixed.values())


def trace3(point):
    P = E6(point)
    # Field-membership tests are made by each caller outside the trace timer.
    total = E6(0)
    x, y = P[0], P[1]
    for _ in range(3):
        total += E6(x, y)
        x, y = x ** p, y ** p
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
    point = (X, *Y)
    dt.check_model(point)
    return point


rng = random.Random(base["pilot_seed"])
anchor_x = [Integer(rng.randrange(p))
            for _ in range(base["anchor_pilot_index"] + 1)][-1]
anchor_cover = cover_from_x(anchor_x)
anchor_y = dt.check_model(anchor_cover)[0]
anchor = E6(dt.c0 + dt.d0 / (anchor_cover[0] - dt.inf),
            dt.kappa * anchor_y / (anchor_cover[0] - dt.inf) ** 2)
anchor_fixed = bool(anchor[0] ** (p ** 3) == anchor[0]
                    and anchor[1] ** (p ** 3) == anchor[1])
anchor_trace6 = projector.trace(anchor)
anchor_trace3 = trace3(anchor) if anchor_fixed else None
if anchor_fixed:
    assert anchor_trace6 == 2 * anchor_trace3
scalar3 = Integer(2 * projector.scalar % dt.r)

elapsed = {"inverse_ns": 0, "quotient_ns": 0,
           "trace3_ns": 0, "trace6_ns": 0}
column_keys = []
for row in sample:
    for residual, frozen in zip(row["residual_points"],
                                row["projected_residuals"]):
        u, v = map(int, residual)
        start = time.monotonic_ns()
        cover, ys = inverse.invert_affine(dt.gm.F(u), dt.gm.F(v))
        elapsed["inverse_ns"] += time.monotonic_ns() - start
        start = time.monotonic_ns()
        image = inverse.quotient_image(E6, cover, ys, 0)
        elapsed["quotient_ns"] += time.monotonic_ns() - start
        assert image[0] ** (p ** 3) == image[0]
        assert image[1] ** (p ** 3) == image[1]
        order = ("trace3", "trace6") if len(column_keys) % 2 else (
            "trace6", "trace3")
        traces = {}
        for variant in order:
            start = time.monotonic_ns()
            traces[variant] = (trace3(image) if variant == "trace3"
                               else projector.trace(image))
            elapsed[variant + "_ns"] += time.monotonic_ns() - start
        assert traces["trace6"] == 2 * traces["trace3"]
        original = projector.scalar * (traces["trace6"] - anchor_trace6)
        if anchor_fixed:
            optimized = scalar3 * (traces["trace3"] - anchor_trace3)
        else:
            optimized = projector.scalar * (
                2 * traces["trace3"] - anchor_trace6)
        expected = E0(*[dt.gm.F(int(x)) for x in frozen])
        assert original == optimized == expected
        column_keys.append(projector.key(optimized))

out = {
    "status": "verified_p192_three_trace_subfield_256_columns",
    "rational_cuts": len(sample),
    "residual_images": len(column_keys),
    "all_model_constants_fixed_by_frobenius_cube": constant_fixed,
    "anchor_fixed_by_frobenius_cube": anchor_fixed,
    "trace6_equals_twice_trace3": True,
    "three_trace_projection_matches_frozen_columns": True,
    "subfield_membership_checks_outside_trace_timers": True,
    "column_digest_sha256": hashlib.sha256(json.dumps(
        column_keys, separators=(",", ":")).encode()).hexdigest(),
    "alternating_unisolated_elapsed_ns": elapsed,
    "source_sha256": sha(__file__),
    "runtime_info_sha256": sha(runtime_path),
    "input_sha256": {path.name: sha(path) for path in
                     (sample_path, base_path)},
    "projection_source_sha256": sha(HERE / "frobenius_trace_projection.py"),
    "inverse_source_sha256": sha(HERE / "plane_point_inverse_cached.py"),
    "boundary": "Exact equivalence on both residuals of 128 frozen cuts. The field closure and trace identity hold on the inverse open locus; trace timers exclude separately performed subfield-membership checks and are exploratory with concurrent producers.",
}
(HERE / "type_i_three_trace_subfield_256.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("THREE_TRACE_SUBFIELD_PASS", len(column_keys), anchor_fixed, flush=True)
