"""Independent prime-field replay and capacity audit for split-line rows."""

import hashlib
import json
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
receipt_path = HERE / "type_i_line_leaf_relation_probe.json"
receipt = json.loads(receipt_path.read_text())
prior = json.loads((HERE / "type_i_line_guided_collector.json").read_text())
projection = json.loads((HERE / "type_i_line_projection.json").read_text())
p = int(prior["field_prime"])
r = int(prior["subgroup_order"])
a = p - 3
b = int("64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1", 16)
assert receipt["source_sha256"] == sha(HERE / "line_leaf_relation_probe.py")
assert receipt["runtime_sha256"] == sha(HERE / "type_i_line_leaf_runtime_info.json")
assert receipt["prior_sha256"] == sha(HERE / "type_i_line_guided_collector.json")
assert receipt["projection_sha256"] == sha(HERE / "type_i_line_projection.json")
for name, digest in receipt["input_sha256"].items():
    assert sha(HERE / name) == digest
assert projection["new_projected_distinct"] == 9784
assert projection["new_projected_sign_distinct"] == 9784
assert projection["base141_signed_matches"] == 0
assert len(receipt["line_rows"]) == receipt["line_count"] == 32


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    u, v = Q
    if x == u and (y + v) % p == 0:
        return None
    slope = ((3 * x * x + a) * pow(2 * y, -1, p) if P == Q
             else (v - y) * pow((u - x) % p, -1, p)) % p
    z = (slope * slope - x - u) % p
    return z, (slope * (x - z) - y) % p


def point(xy):
    x, y = map(int, xy)
    assert 0 <= x < p and 0 <= y < p
    assert (y * y - x * x * x - a * x - b) % p == 0
    return x, y


base = [point(xy) for xy in prior["base_points_fp"]]
assert len(base) == len(set(base)) == 141
pair_rows = [tuple(map(int, row.split())) for row in
             (HERE / "type_i_line_new_point_pairs.tsv").read_text().splitlines()]
assert len(pair_rows) == 9784
assert all(pair_rows[2 * k] == pair_rows[2 * k + 1]
           for k in range(4892))
assert len(set(pair_rows[::2])) == 4892
fresh = set()
constant = None
for k, row in enumerate(receipt["line_rows"]):
    i, j = row["base_indices_1based"]
    assert (i, j) == pair_rows[2 * k]
    A, B = [point(xy) for xy in row["fresh_points_fp"]]
    assert A != B and A not in fresh and B not in fresh
    fresh.update((A, B))
    total = add(add(base[i - 1], base[j - 1]), add(A, B))
    if constant is None:
        constant = total
    assert total == constant
    assert row["constant_key"] == f"{total[0]}:{total[1]}"
assert len(fresh) == 64
assert receipt["constant_key"] == f"{constant[0]}:{constant[1]}"
assert not receipt["constant_is_identity"]

B = 141 + 9784
capacity = comb(B + 2, 3)
assert B == 9925 and capacity == 162993977475
assert 64 * capacity < r - 1
out = {
    "status": "LINE_LEAF_INDEPENDENT_PASS",
    "sample_lines_replayed": 32,
    "sample_fresh_points_replayed": len(fresh),
    "full_line_count_from_incidence": 4892,
    "full_distinct_fresh_from_projection": 9784,
    "full_base_size": B,
    "three_sum_capacity_ceiling": capacity,
    "uniform_64_query_union_bound_num": 64 * capacity,
    "uniform_64_query_union_bound_den": r - 1,
    "receipt_sha256": sha(receipt_path),
    "source_sha256": sha(__file__),
}
(HERE / "type_i_line_leaf_relation_verification.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("LINE_LEAF_INDEPENDENT_PASS", json.dumps({
    "sample_lines": 32, "full_base_size": B,
    "three_sum_capacity_ceiling": capacity}, sort_keys=True), flush=True)
