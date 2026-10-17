"""Independent reverse-pivot replay of all old-anchored residual points."""

import hashlib
import json
import time
from collections import Counter
from itertools import combinations
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
producer_path = HERE / "type_i_line_nonsource_reuse_old141.json"
producer = json.loads(producer_path.read_text())
plane_path = HERE / "type_i_plane_quartic.json"
old_path = HERE / "type_i_line_base_points.tsv"
new_path = HERE / "type_i_line_new_points.tsv"
lines_path = HERE / "type_i_line_full_collinearity.json"
plane = json.loads(plane_path.read_text())
lines = json.loads(lines_path.read_text())
p = int(plane["p"])
old = [tuple(map(int, line.split()[1:])) for line in old_path.read_text().splitlines()]
new = [tuple(map(int, line.split())) for line in new_path.read_text().splitlines()]
points = old + new
selected = set(points)
assert len(old) == 141 and len(new) == 9784 and len(selected) == len(points) == 9925
assert producer["old_pivots"] == 141


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


assert producer["input_sha256"]["quartic"] == sha(plane_path)
assert producer["input_sha256"]["old"] == sha(old_path)
assert producer["input_sha256"]["new"] == sha(new_path)
assert producer["input_sha256"]["collinearity"] == sha(lines_path)

source_pairs = set()
for record in lines["line_records"]:
    vertices = list(map(int, record["point_indices_zero_based"]))
    assert len(vertices) == 4
    source_pairs.update((a, b) for a, b in combinations(sorted(vertices), 2)
                        if a < 141)

coefficients = [(tuple(map(int, label.split(","))), int(value) % p)
                for label, value in plane["quartic_coefficients"].items()]


def directional_coefficients(x, y):
    """Coefficient of t^k m^s in F(x+t,y+mt,1)."""
    by_degree = [[0] * (k + 1) for k in range(5)]
    for (a, b, _), value in coefficients:
        for k in range(a + b + 1):
            for s in range(max(0, k - a), min(b, k) + 1):
                by_degree[k][s] += (value * comb(a, k - s) * comb(b, s)
                                     * pow(x, a - k + s, p) * pow(y, b - s, p))
    assert by_degree[0][0] % p == 0
    return [[c % p for c in degree] for degree in by_degree]


def horner(values, point):
    result = 0
    for value in reversed(values):
        result = (result * point + value) % p
    return result


def batch_inverse(values):
    if not values:
        return []
    prefix = [1]
    for value in values:
        assert value % p
        prefix.append(prefix[-1] * value % p)
    inv = pow(prefix[-1], -1, p)
    answer = [0] * len(values)
    for k in range(len(values) - 1, -1, -1):
        answer[k] = inv * prefix[k] % p
        inv = inv * values[k] % p
    assert inv == 1
    return answer


def check_quartic(point):
    x, y = point
    return sum(value * pow(x, a, p) * pow(y, b, p)
               for (a, b, _), value in coefficients) % p == 0


counts = Counter()
seen = set()
first_repeats = []
start_ns = time.monotonic_ns()
for i in range(140, -1, -1):
    x, y = points[i]
    coeff = directional_coefficients(x, y)
    for start in range(i + 1, len(points), 1024):
        stop = min(len(points), start + 1024)
        neighbors = list(range(start, stop))
        differences = [(points[j][0] - x) % p for j in neighbors]
        slopes_denominators = batch_inverse(differences)
        split_cases = []
        for j, d, inv_d in zip(neighbors, differences, slopes_denominators):
            counts["selected_pairs"] += 1
            if (i, j) in source_pairs:
                counts["source_pair_presentations"] += 1
                continue
            slope = (points[j][1] - y) * inv_d % p
            c1, c2, c3, c4 = (horner(coeff[k], slope) for k in range(1, 5))
            A = c4
            B = (c3 + d * A) % p
            C = (c2 + d * B) % p
            assert A and C and (c1 + d * C) % p == 0
            assert (A * d * d + B * d + C) % p
            discriminant = (B * B - 4 * A * C) % p
            assert discriminant
            symbol = pow(discriminant, (p - 1) // 2, p)
            if symbol == p - 1:
                counts["nonsplit_lines"] += 1
                continue
            assert symbol == 1
            root = pow(discriminant, (p + 1) // 4, p)
            assert root * root % p == discriminant
            split_cases.append((j, d, slope, A, B, C, root))
            counts["split_lines"] += 1
        denominator_inverses = batch_inverse([2 * row[3] % p for row in split_cases])
        for (j, d, slope, A, B, C, root), inv in zip(split_cases, denominator_inverses):
            for sign in (-1, 1):
                t = (-B + sign * root) * inv % p
                assert t not in (0, d) and (A * t * t + B * t + C) % p == 0
                point = ((x + t) % p, (y + slope * t) % p)
                assert point not in selected
                key = point[0].to_bytes(24, "big") + point[1].to_bytes(24, "big")
                if key in seen:
                    counts["reused_residual_presentations"] += 1
                    if len(first_repeats) < 100:
                        first_repeats.append({
                            "selected_pair": [i, j],
                            "point": [str(value) for value in point],
                        })
                else:
                    seen.add(key)
                counts["residual_presentations"] += 1
                if counts["residual_presentations"] % 3001 == 0:
                    assert check_quartic(point)
                    counts["direct_quartic_controls"] += 1
    if i % 8 == 0:
        print("NONSOURCE_REVERSE_PROGRESS", 141 - i, 141,
              counts["split_lines"],
              counts.get("reused_residual_presentations", 0), flush=True)

assert counts["selected_pairs"] == 1389414
assert counts["source_pair_presentations"] == 24460
assert counts["split_lines"] + counts["nonsplit_lines"] == 1389414 - 24460
assert counts["residual_presentations"] == 2 * counts["split_lines"]
assert len(seen) + counts.get("reused_residual_presentations", 0) == counts["residual_presentations"]
assert counts.get("reused_residual_presentations", 0) == 0
for key in ("selected_pairs", "source_pair_presentations", "split_lines",
            "nonsplit_lines", "residual_presentations"):
    assert counts[key] == producer["counts"][key], key
assert len(seen) == producer["distinct_outside_residual_points"]

out = {
    "status": "verified_reverse_old_pivot_residual_reuse_replay",
    "counts": dict(counts),
    "distinct_outside_residual_points": len(seen),
    "reuse_witnesses_first_100": first_repeats,
    "producer_sha256": sha(producer_path),
    "source_sha256": sha(Path(__file__)),
    "input_sha256": {"quartic": sha(plane_path), "old": sha(old_path),
                     "new": sha(new_path), "collinearity": sha(lines_path)},
    "elapsed_ns_unisolated": time.monotonic_ns() - start_ns,
    "boundary": "Complete independent old-pivot replay in descending pivot and ascending block order, with separate 1,024-line inversions and residual-point set.",
}
path = HERE / "type_i_line_nonsource_reuse_old141_reverse.json"
path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print("NONSOURCE_REVERSE_PASS", json.dumps({
    "split_lines": counts["split_lines"],
    "outside_residual_points": len(seen),
    "reused_residual_presentations": counts.get("reused_residual_presentations", 0),
}, sort_keys=True))
