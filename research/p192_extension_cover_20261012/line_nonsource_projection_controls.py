"""Project the first 128 new rational line cuts into the P-192 subgroup."""

import argparse
import hashlib
import json
import random
import time
from pathlib import Path

from sage.all import EllipticCurve, Integer, prod, vector

import type_i_divisor_transfer as dt
from plane_point_inverse_cached import CachedPlanePointInverse
from plane_point_inverse import PlanePointInverse, elliptic_quotient_image

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--screen", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
screen = json.loads(args.screen.read_text())
assert screen["status"] == "verified_nonsource_old_pivot_residual_reuse_screen"
assert len(screen["first_128_split_cases"]) == 128
base_path = HERE / "type_i_quotient_base_probe.json"
old_path = HERE / "type_i_line_guided_collector.json"
new_path = HERE / "type_i_line_full_new_quotient_points.json"
base = json.loads(base_path.read_text())
old = json.loads(old_path.read_text())
new = json.loads(new_path.read_text())
p = dt.p
E = EllipticCurve(dt.L, [dt.L(dt.gm.a), dt.L(dt.gm.b)])
O = E(0)
h = Integer(base["candidate_cofactor"])
assert h * dt.r == Integer(base["extension_elliptic_order"])
inverse = CachedPlanePointInverse(dt)
reference_inverse = PlanePointInverse(dt)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cover_from_x(value):
    X = dt.L(value)
    u = (X - dt.zeta) * (X - dt.beta)
    y = u ** ((dt.q + 1) // 2)
    norm = prod(X - root ** (p ** j)
                for root in (dt.zeta, dt.beta) for j in range(3))
    if y * y ** p * y ** (p * p) == -norm:
        y = -y
    assert y * y ** p * y ** (p * p) == norm
    Y = dt.Vandermonde.solve_right(vector(dt.L, [y, y ** p, y ** (p * p)]))
    point = (X, *Y)
    dt.check_model(point)
    return point


rng = random.Random(base["pilot_seed"])
anchor_index = base["anchor_pilot_index"]
anchor_x = [Integer(rng.randrange(int(p))) for _ in range(anchor_index + 1)][-1]
anchor_point = cover_from_x(anchor_x)
anchor_y = dt.check_model(anchor_point)[0]
anchor = E(dt.c0 + dt.d0 / (anchor_point[0] - dt.inf),
           dt.kappa * anchor_y / (anchor_point[0] - dt.inf) ** 2)


def key(Q):
    assert Q != O and Q[0] ** p == Q[0] and Q[1] ** p == Q[1]
    return (int(Q[0].polynomial()[0]), int(Q[1].polynomial()[0]))


selected = {tuple(map(int, point)) for point in old["base_points_fp"]}
selected.update(tuple(map(int, point)) for point in new)
assert len(selected) == 9925
seen = {}
matches = []
start_ns = time.monotonic_ns()
for case_index, case in enumerate(screen["first_128_split_cases"]):
    for point_index, (u, v) in enumerate(case["residual_points"]):
        cover, ys = inverse.invert_affine(dt.gm.F(int(u)), dt.gm.F(int(v)))
        image = inverse.quotient_image(E, cover, ys, 0)
        if case_index < 16:
            reference_cover = reference_inverse.invert_affine(
                dt.gm.F(int(u)), dt.gm.F(int(v)))
            assert reference_cover == cover
            assert elliptic_quotient_image(dt, E, reference_cover, 0) == image
        Q = h * (image - anchor)
        assert dt.r * Q == O
        point_key = key(Q)
        if point_key in selected or point_key in seen:
            matches.append({
                "case_index": case_index,
                "point_index": point_index,
                "projected": [str(value) for value in point_key],
                "selected_match": point_key in selected,
                "prior_case": seen.get(point_key),
            })
        seen.setdefault(point_key, [case_index, point_index])
    if (case_index + 1) % 16 == 0:
        print("NONSOURCE_PROJECTION_PROGRESS", case_index + 1, 128,
              len(matches), flush=True)

out = {
    "status": "verified_nonsource_residual_projection_controls",
    "line_cuts": 128,
    "residual_points": 256,
    "reference_inverse_and_quotient_controls": 32,
    "distinct_projected_classes": len(seen),
    "selected_or_repeated_matches": matches,
    "elapsed_ns_unisolated": time.monotonic_ns() - start_ns,
    "source_sha256": sha(Path(__file__)),
    "screen_sha256": sha(args.screen),
    "base_sha256": sha(base_path),
    "old_sha256": sha(old_path),
    "new_sha256": sha(new_path),
    "inverse_source_sha256": sha(HERE / "plane_point_inverse_cached.py"),
    "reference_inverse_source_sha256": sha(HERE / "plane_point_inverse.py"),
    "runtime_info_sha256": sha(HERE / "type_i_line_new64_runtime_info.json"),
    "boundary": "First 128 rational non-source cuts only; exact point projection and subgroup checks for their 256 residuals.",
}
args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print("NONSOURCE_PROJECTION_PASS", json.dumps({
    "residual_points": 256,
    "distinct_projected_classes": len(seen),
    "matches": len(matches),
}, sort_keys=True))
