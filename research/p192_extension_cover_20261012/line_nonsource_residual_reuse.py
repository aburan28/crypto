"""Exact rational-residual reuse on frozen old-point pivot line cuts.

The complete collinearity census identifies every source line. A non-source
line has two selected endpoints and, by the transverse-cut certificate, two
distinct residual points outside the frozen 9,925-point set. We compute
those residuals and test whether independent cuts reuse an outside point.
"""

import argparse
import hashlib
import json
import time
from collections import Counter
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--old-pivots", type=int, default=64)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
assert 1 <= args.old_pivots <= 141
inputs = {
    "quartic": HERE / "type_i_plane_quartic.json",
    "old": HERE / "type_i_line_base_points.tsv",
    "new": HERE / "type_i_line_new_points.tsv",
    "collinearity": HERE / "type_i_line_full_collinearity.json",
    "transverse": HERE / "type_i_line_full_residual_discriminant.json",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


plane = json.loads(inputs["quartic"].read_text())
collinearity = json.loads(inputs["collinearity"].read_text())
transverse = json.loads(inputs["transverse"].read_text())
p = int(plane["p"])
assert p % 4 == 3
old = [tuple(map(int, row.split()[1:])) for row in
       inputs["old"].read_text().splitlines()]
new = [tuple(map(int, row.split())) for row in
       inputs["new"].read_text().splitlines()]
points = old + new
n = len(points)
selected = set(points)
assert len(old) == 141 and len(new) == 9784 and len(selected) == n == 9925
assert collinearity["source_quartets"] == 4892
assert collinearity["novel_three_point_lines"] == 0
assert collinearity["novel_four_point_lines"] == 0
assert transverse["counts"] == {
    "pairs": 49247850, "vertical_pairs": 0, "degree_below_four": 0,
    "residual_double_root": 0, "tangent_at_first": 0,
    "tangent_at_second": 0,
}
for name in ("quartic", "old", "new"):
    assert collinearity["input_sha256"][inputs[name].name] == sha(inputs[name])
    assert transverse["input_sha256"][inputs[name].name] == sha(inputs[name])

source_neighbors = [set() for _ in range(args.old_pivots)]
for record in collinearity["line_records"]:
    indices = list(map(int, record["point_indices_zero_based"]))
    assert len(indices) == 4
    for i in indices:
        if i < args.old_pivots:
            source_neighbors[i].update(j for j in indices if j > i)

terms = [(tuple(map(int, label.split(","))), int(value) % p)
         for label, value in plane["quartic_coefficients"].items()]
assert all(sum(exponents) == 4 for exponents, _ in terms)


def taylor_at(x, y):
    xp = [pow(x, i, p) for i in range(5)]
    yp = [pow(y, i, p) for i in range(5)]
    out = [[0] * (k + 1) for k in range(5)]
    for (a, b, _), value in terms:
        for ix in range(a + 1):
            for iy in range(b + 1):
                k = ix + iy
                out[k][iy] = (out[k][iy] + value * comb(a, ix)
                               * comb(b, iy) * xp[a - ix] * yp[b - iy]) % p
    assert out[0][0] == 0
    return out


def evaluate(coefficients, x):
    result = 0
    for coefficient in reversed(coefficients):
        result = (result * x + coefficient) % p
    return result


def on_quartic(x, y):
    return sum(value * pow(x, a, p) * pow(y, b, p)
               for (a, b, _), value in terms) % p == 0


counts = Counter()
first_split_cases = []
reuse_witnesses = []
seen = {}
residual_stream = hashlib.sha256()
start_ns = time.monotonic_ns()
for i in range(args.old_pivots):
    x, y = points[i]
    t = taylor_at(x, y)
    js = list(range(i + 1, n))
    prefix = []
    product = 1
    for j in js:
        d = (points[j][0] - x) % p
        assert d != 0
        product = product * d % p
        prefix.append(product)
    inverse_product = pow(product, -1, p)
    split = []
    for k in range(len(js) - 1, -1, -1):
        j = js[k]
        u, v = points[j]
        d = (u - x) % p
        inverse_d = inverse_product * (prefix[k - 1] if k else 1) % p
        inverse_product = inverse_product * d % p
        counts["selected_pairs"] += 1
        if j in source_neighbors[i]:
            counts["source_pair_presentations"] += 1
            continue
        m = (v - y) * inverse_d % p
        c1, c2, c3, c4 = (evaluate(t[z], m) for z in range(1, 5))
        A = c4
        B = (c3 + d * A) % p
        C = (c2 + d * B) % p
        assert A != 0 and C != 0
        assert (c1 + d * C) % p == 0
        assert (A * d * d + B * d + C) % p != 0
        discriminant = (B * B - 4 * A * C) % p
        assert discriminant != 0
        symbol = pow(discriminant, (p - 1) // 2, p)
        if symbol == p - 1:
            counts["nonsplit_lines"] += 1
            continue
        assert symbol == 1
        square_root = pow(discriminant, (p + 1) // 4, p)
        assert square_root * square_root % p == discriminant
        split.append((j, m, d, A, B, C, square_root))
        counts["split_lines"] += 1
    assert inverse_product == 1

    # One inversion for all accepted denominators at this pivot.
    denominators = []
    denominator_prefix = []
    denominator_product = 1
    for _, _, _, A, _, _, _ in split:
        denominator = 2 * A % p
        denominator_product = denominator_product * denominator % p
        denominators.append(denominator)
        denominator_prefix.append(denominator_product)
    inverse_product = pow(denominator_product, -1, p)
    for k in range(len(split) - 1, -1, -1):
        j, m, d, A, B, C, square_root = split[k]
        inverse = inverse_product * (denominator_prefix[k - 1] if k else 1) % p
        inverse_product = inverse_product * denominators[k] % p
        roots = [((-B + sign * square_root) * inverse) % p
                 for sign in (1, -1)]
        assert roots[0] != roots[1]
        assert all(root not in (0, d) for root in roots)
        residuals = []
        for root in roots:
            assert (A * root * root + B * root + C) % p == 0
            point = ((x + root) % p, (y + m * root) % p)
            assert point not in selected
            key = point[0].to_bytes(24, "big") + point[1].to_bytes(24, "big")
            previous = seen.get(key)
            if previous is not None:
                counts["reused_residual_presentations"] += 1
                if len(reuse_witnesses) < 100:
                    reuse_witnesses.append({
                        "point": [str(value) for value in point],
                        "first_pair": previous,
                        "second_pair": [i, j],
                    })
            else:
                seen[key] = [i, j]
            residual_stream.update(key)
            residuals.append(point)
            counts["residual_presentations"] += 1
            if counts["residual_presentations"] % 2003 == 0:
                assert on_quartic(*point)
                counts["direct_quartic_controls"] += 1
        if len(first_split_cases) < 128:
            first_split_cases.append({
                "selected_pair": [i, j],
                "residual_points": [[str(v) for v in point] for point in residuals],
            })
    assert inverse_product == 1
    print("NONSOURCE_REUSE_PROGRESS", i + 1, args.old_pivots,
          counts["split_lines"], counts["reused_residual_presentations"],
          flush=True)

assert counts["selected_pairs"] == sum(n - i - 1 for i in range(args.old_pivots))
assert counts["selected_pairs"] == (counts["source_pair_presentations"] +
                                    counts["split_lines"] +
                                    counts["nonsplit_lines"])
assert counts["residual_presentations"] == 2 * counts["split_lines"]
assert len(seen) + counts["reused_residual_presentations"] == counts["residual_presentations"]
out = {
    "status": "verified_nonsource_old_pivot_residual_reuse_screen",
    "field_prime": str(p),
    "selected_point_count": n,
    "old_pivots": args.old_pivots,
    "counts": dict(counts),
    "distinct_outside_residual_points": len(seen),
    "first_128_split_cases": first_split_cases,
    "reuse_witnesses_first_100": reuse_witnesses,
    "residual_stream_sha256": residual_stream.hexdigest(),
    "input_sha256": {name: sha(path) for name, path in inputs.items()},
    "source_sha256": sha(Path(__file__)),
    "elapsed_ns_unisolated": time.monotonic_ns() - start_ns,
    "boundary": "Old-pivot selected-pair lines only. The full selected-pair transverse and collinearity certificates establish non-source status and distinct outside residuals; this screen measures rational splitting and point reuse on the declared pivot subset.",
}
args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print("NONSOURCE_REUSE_PASS", json.dumps({
    "old_pivots": args.old_pivots,
    "selected_pairs": counts["selected_pairs"],
    "split_lines": counts["split_lines"],
    "distinct_outside_residual_points": len(seen),
    "reused_residual_presentations": counts["reused_residual_presentations"],
}, sort_keys=True), flush=True)
