"""Certify the 32-line common-point pencil inside the 205-class base."""

import collections
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from line_full_character_support import rank_mod  # noqa: E402
from verify_adaptive_private_columns import add, on_curve  # noqa: E402


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


paths = {
    "old_quotient": HERE / "type_i_line_guided_collector.json",
    "new_quotient": HERE / "type_i_line_new64_fixed_origin.json",
    "old_quartic": HERE / "type_i_line_base_points.tsv",
    "new_quartic": HERE / "type_i_line_new_points.tsv",
    "source_pairs": HERE / "type_i_line_new_point_pairs.tsv",
}
old = json.loads(paths["old_quotient"].read_text())
new = json.loads(paths["new_quotient"].read_text())
p, r = int(old["field_prime"]), int(old["subgroup_order"])
old_q = [tuple(map(int, point)) for point in old["base_points_fp"]]
new_q = [tuple(map(int, row["quotient"])) for row in new["rows"]]
assert len(old_q) == 141 and len(new_q) == 64
assert len(set(old_q + new_q)) == 205
assert all(on_curve(point, p) for point in old_q + new_q)
old_c = [tuple(map(int, line.split()[1:])) for line in
         paths["old_quartic"].read_text().splitlines()]
new_c = [tuple(map(int, line.split())) for line in
         paths["new_quartic"].read_text().splitlines()[:64]]
pairs = [tuple(map(int, line.split())) for line in
         paths["source_pairs"].read_text().splitlines()[:64]]
assert len(old_c) == 141 and len(new_c) == len(pairs) == 64
assert all(pairs[2 * k] == pairs[2 * k + 1] for k in range(32))
assert {pairs[2 * k][0] for k in range(32)} == {1}


def collinear(A, B, C):
    return ((B[0] - A[0]) * (C[1] - A[1]) -
            (B[1] - A[1]) * (C[0] - A[0])) % p == 0


triples = []
targets = []
for k in range(32):
    i, j = pairs[2 * k]
    assert i == 1 and 1 <= j <= 141 and j != i
    A, B = old_c[i - 1], old_c[j - 1]
    X, Y = new_c[2 * k], new_c[2 * k + 1]
    assert X != Y and all(collinear(A, B, point) for point in (X, Y))
    triple = (j - 1, 141 + 2 * k, 141 + 2 * k + 1)
    triples.append(triple)
    targets.append(add(add(old_q[j - 1], new_q[2 * k], p), new_q[2 * k + 1], p))
assert len(set(triples)) == 32 and len(set(targets)) == 1
reference = collections.Counter(triples[0])
rows = []
for triple in triples[1:]:
    row = collections.Counter(triple)
    row.subtract(reference)
    rows.append(dict(row))
assert rank_mod(rows, r) == 31
assert all(row[141 + 2 * k] == row[141 + 2 * k + 1] == 1
           for k, row in enumerate(rows, start=1))
out = {
    "status": "verified_32_line_pencil_three_sum_collision",
    "field_prime": str(p), "subgroup_order": str(r),
    "selected_old_points": 141, "selected_new_points": 64,
    "source_lines": 32, "common_old_index_one_based": 1,
    "distinct_three_sum_decompositions": 32,
    "independent_collision_rows": 31,
    "common_target": [str(value) for value in targets[0]],
    "triple_indices_zero_based": triples,
    "input_sha256": {name: sha(path) for name, path in paths.items()},
    "source_sha256": sha(Path(__file__)),
    "boundary": "The 31 collision rows contain private new-point columns; the old-column projection remains unconstrained by this pencil alone.",
}
path = HERE / "type_i_line_pencil_three_sum.json"
path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print("PENCIL_THREE_SUM_PASS", json.dumps({
    "source_lines": 32, "distinct_decompositions": 32,
    "collision_rank": 31,
}, sort_keys=True))
