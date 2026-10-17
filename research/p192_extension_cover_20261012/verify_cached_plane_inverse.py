"""Exact paired replay of cached and established Type I point inverses.

The alternating local durations are diagnostics only; this script does not
claim a host-isolated wall-time improvement.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import EllipticCurve

import type_i_divisor_transfer as dt
from plane_point_inverse import PlanePointInverse, elliptic_quotient_image
from plane_point_inverse_cached import CachedPlanePointInverse

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
F, L, p, r = dt.gm.F, dt.L, dt.p, dt.r
E = EllipticCurve(L, [L(dt.gm.a), L(dt.gm.b)])
legacy = PlanePointInverse(dt)
cached = CachedPlanePointInverse(dt)
base_rows = [tuple(map(int, line.split()[1:])) for line in
             (HERE / "type_i_line_base_points.tsv").read_text().splitlines()]
new_rows = [tuple(map(int, line.split())) for line in
            (HERE / "type_i_line_new_points.tsv").read_text().splitlines()]
rows = [("old", i, *row) for i, row in enumerate(base_rows[:16])]
rows += [("new", i, *row) for i, row in enumerate(new_rows[:64])]
assert len(rows) == 80
elapsed = {"legacy": 0, "cached": 0}
controls = []
for n, (family, index, u, v) in enumerate(rows):
    results = {}
    order = ("cached", "legacy") if n % 2 else ("legacy", "cached")
    for variant in order:
        start = time.monotonic_ns()
        if variant == "legacy":
            point = legacy.invert_affine(F(u), F(v))
            images = [elliptic_quotient_image(dt, E, point, j)
                      for j in range(3)]
        else:
            point, y = cached.invert_affine(F(u), F(v))
            images = [cached.quotient_image(E, point, y, j)
                      for j in range(3)]
        elapsed[variant] += time.monotonic_ns() - start
        results[variant] = (point, images)
    assert results["legacy"] == results["cached"]
    point, images = results["legacy"]
    controls.append({"family": family, "index": index,
                     "point_x": str(point[0]),
                     "quotient0_x": str(images[0][0]),
                     "quotient0_y": str(images[0][1])})
out = {
    "status": "CACHED_PLANE_INVERSE_EQUIVALENT_80_POINTS",
    "base_controls": 16, "new_controls": 64,
    "all_three_quotient_images_equal": True,
    "alternating_local_elapsed_ns_unisolated": elapsed,
    "controls_sha256": hashlib.sha256(
        json.dumps(controls, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest(),
    "source_sha256": sha(__file__),
    "cached_source_sha256": sha(HERE / "plane_point_inverse_cached.py"),
    "legacy_source_sha256": sha(HERE / "plane_point_inverse.py"),
    "runtime_sha256": sha(HERE / "type_i_cached_inverse_runtime_info.json"),
    "input_sha256": {name: sha(HERE / name) for name in
                     ("type_i_line_base_points.tsv",
                      "type_i_line_new_points.tsv")},
    "boundary": "Exact point and all three quotient-image replay on 80 frozen inputs; local alternating elapsed times are exploratory diagnostics without isolation or A/A noise gates.",
}
(HERE / "type_i_cached_inverse_verification.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("CACHED_PLANE_INVERSE_PASS", len(rows), flush=True)
