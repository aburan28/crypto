"""Census lines through at least three of all 9,925 frozen quartic points.

For each pivot, batch-invert nonzero abscissa differences and group equal
slopes. This visits every unordered pair, using O(number of points) working
memory per pivot rather than retaining all 49 million pair records.
"""

import argparse
import hashlib
import json
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUTS = ("type_i_plane_quartic.json", "type_i_line_base_points.tsv",
          "type_i_line_new_points.tsv", "type_i_line_new_point_pairs.tsv")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_points():
    old = [tuple(map(int, row.split()[1:])) for row in
           (HERE / INPUTS[1]).read_text().splitlines()]
    new = [tuple(map(int, row.split())) for row in
           (HERE / INPUTS[2]).read_text().splitlines()]
    pairs = [tuple(map(int, row.split())) for row in
             (HERE / INPUTS[3]).read_text().splitlines()]
    assert len(old) == 141 and len(new) == len(pairs) == 9784
    assert len(set(old + new)) == 9925
    return old, new, pairs


def census(points, p, pivot_limit):
    n = len(points)
    lines = {}
    pair_count = 0
    start = time.monotonic_ns()
    for i in range(pivot_limit):
        x0, y0 = points[i]
        js, prefixes, vertical = [], [], []
        product = 1
        for j in range(i + 1, n):
            xj = points[j][0]
            dx = (xj - x0) % p
            if dx == 0:
                vertical.append(j)
            else:
                product = product * dx % p
                js.append(j)
                prefixes.append(product)
            pair_count += 1
        if len(vertical) >= 2:
            lines.setdefault((-1, x0), set()).update((i, *vertical))
        if js:
            inverse_product = pow(product, -1, p)
            first, groups = {}, {}
            for k in range(len(js) - 1, -1, -1):
                j = js[k]
                dx = (points[j][0] - x0) % p
                inv_dx = inverse_product * (prefixes[k - 1] if k else 1) % p
                inverse_product = inverse_product * dx % p
                slope = (points[j][1] - y0) * inv_dx % p
                if slope in first:
                    groups.setdefault(slope, {first[slope]}).add(j)
                else:
                    first[slope] = j
            assert inverse_product == 1
            for slope, rest in groups.items():
                intercept = (y0 - slope * x0) % p
                lines.setdefault((slope, intercept), set()).update(rest | {i})
        if (i + 1) % 250 == 0:
            print("FULL_COLLINEARITY_PROGRESS", i + 1, pivot_limit,
                  pair_count, len(lines), flush=True)
    return lines, pair_count, time.monotonic_ns() - start


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pivot-limit", type=int)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    plane = json.loads((HERE / INPUTS[0]).read_text())
    p = int(plane["p"])
    old, new, pairs = load_points()
    points = old + new
    n = len(points)
    pivot_limit = n if args.pivot_limit is None else args.pivot_limit
    assert 1 <= pivot_limit <= n
    lines, pair_count, elapsed_ns = census(points, p, pivot_limit)
    if pivot_limit == n:
        assert pair_count == n * (n - 1) // 2
    records = []
    for (slope, intercept), members in sorted(lines.items()):
        indices = sorted(members)
        assert 3 <= len(indices) <= 4
        for index in indices:
            x, y = points[index]
            assert ((x - intercept) % p == 0 if slope == -1 else
                    (y - slope * x - intercept) % p == 0)
        records.append({"line": [str(slope), str(intercept)],
                        "point_indices_zero_based": indices})
    source = set()
    for k in range(4892):
        i, j = pairs[2 * k]
        assert pairs[2 * k + 1] == (i, j)
        source.add(tuple(sorted((i - 1, j - 1,
                                 len(old) + 2 * k, len(old) + 2 * k + 1))))
    assert len(source) == 4892
    discovered = {tuple(row["point_indices_zero_based"])
                  for row in records}
    if pivot_limit == n:
        assert source <= discovered
    degree_counts = Counter(len(row["point_indices_zero_based"])
                            for row in records)
    families = Counter((sum(i < len(old) for i in row["point_indices_zero_based"]),
                        sum(i >= len(old) for i in row["point_indices_zero_based"]))
                       for row in records)
    out = {
        "status": ("FULL_9925_POINT_COLLINEARITY_CENSUS" if pivot_limit == n
                   else "PARTIAL_PIVOT_COLLINEARITY_DIAGNOSTIC"),
        "field_prime": str(p), "point_count": n, "old_points": len(old),
        "new_points": len(new), "pivot_limit": pivot_limit,
        "unordered_pairs_visited": pair_count,
        "source_quartets": len(source),
        "source_quartets_found": len(source & discovered),
        "novel_three_point_lines": degree_counts[3],
        "novel_four_point_lines": len({row for row in discovered
                                       if len(row) == 4} - source),
        "selected_point_count_distribution": {str(k): v
                                               for k, v in sorted(degree_counts.items())},
        "old_new_count_distribution": {f"{a}+{b}": v
                                       for (a, b), v in sorted(families.items())},
        "line_records": records,
        "source_sha256": sha(__file__),
        "input_sha256": {name: sha(HERE / name) for name in INPUTS},
        "elapsed_ns_unisolated": elapsed_ns,
        "boundary": "Exact collinearity census on the frozen point set; line incidence alone is not a known-scalar relation or a measured index-calculus cost.",
    }
    path = args.output or (HERE / "type_i_line_full_collinearity.json")
    path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print("FULL_COLLINEARITY_PASS", out["status"], pair_count,
          len(records), dict(degree_counts), flush=True)


if __name__ == "__main__":
    main()
