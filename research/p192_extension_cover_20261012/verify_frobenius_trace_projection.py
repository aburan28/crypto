"""Paired direct/trace subgroup projection on 256 frozen P-192 residuals."""

import argparse
import hashlib
import json
import random
import time
from pathlib import Path

from sage.all import EllipticCurve, Integer, prod, vector

import type_i_divisor_transfer as dt
from frobenius_trace_projection import FrobeniusTraceProjector
from plane_point_inverse_cached import CachedPlanePointInverse

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
screen_path = HERE / "type_i_line_nonsource_reuse_old141.json"
base_path = HERE / "type_i_quotient_base_probe.json"
screen = json.loads(screen_path.read_text())
base = json.loads(base_path.read_text())
assert screen["old_pivots"] == 141 and len(screen["first_128_split_cases"]) == 128
h = Integer(base["candidate_cofactor"])
projector = FrobeniusTraceProjector(dt.gm.F, dt.L, dt.gm.a, dt.gm.b, dt.r, h)
assert projector.E6.cardinality() == h * dt.r
inverse = CachedPlanePointInverse(dt)

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

rng = random.Random(base["pilot_seed"])
anchor_x = [Integer(rng.randrange(int(dt.p)))
            for _ in range(base["anchor_pilot_index"] + 1)][-1]
anchor_point = cover_from_x(anchor_x)
anchor_y = dt.check_model(anchor_point)[0]
E6 = projector.E6
anchor = E6(dt.c0 + dt.d0 / (anchor_point[0] - dt.inf),
            dt.kappa * anchor_y / (anchor_point[0] - dt.inf) ** 2)
anchor_trace = projector.trace(anchor)
direct_ns = trace_ns = inverse_ns = 0
seen = set()
for case_index, case in enumerate(screen["first_128_split_cases"]):
    for u, v in case["residual_points"]:
        start = time.monotonic_ns()
        cover, ys = inverse.invert_affine(dt.gm.F(int(u)), dt.gm.F(int(v)))
        image = inverse.quotient_image(E6, cover, ys, 0)
        inverse_ns += time.monotonic_ns() - start
        start = time.monotonic_ns()
        direct = projector.direct_delta(image, anchor)
        direct_ns += time.monotonic_ns() - start
        start = time.monotonic_ns()
        trace = projector.project_delta(image, anchor_trace)
        trace_ns += time.monotonic_ns() - start
        assert projector.key(direct) == projector.key(trace)
        assert projector.key(trace) is not None
        seen.add(projector.key(trace))
    if (case_index + 1) % 16 == 0:
        print("TRACE_PROJECTION_PROGRESS", case_index + 1, 128, flush=True)
assert len(seen) == 256
out = {
    "status": "verified_p192_six_frobenius_trace_projection",
    "residual_points": 256,
    "paired_direct_equals_trace": 256,
    "distinct_projected_classes": len(seen),
    "field_prime": str(dt.p),
    "subgroup_order": str(dt.r),
    "cofactor": str(h),
    "trace_scalar": str(projector.scalar),
    "inverse_ns_unisolated": inverse_ns,
    "direct_projection_ns_unisolated": direct_ns,
    "trace_projection_ns_unisolated": trace_ns,
    "source_sha256": sha(__file__),
    "projector_source_sha256": sha(HERE / "frobenius_trace_projection.py"),
    "screen_sha256": sha(screen_path),
    "base_sha256": sha(base_path),
    "inverse_source_sha256": sha(HERE / "plane_point_inverse_cached.py"),
    "boundary": "Paired correctness and local component cost on exactly the first 256 residual points; no full outside-point projection or isolated speed claim.",
}
args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print("TRACE_PROJECTION_PASS", json.dumps({
    "paired_points": 256,
    "distinct_projected_classes": len(seen),
    "direct_seconds": direct_ns / 1e9,
    "trace_seconds": trace_ns / 1e9,
}, sort_keys=True))
