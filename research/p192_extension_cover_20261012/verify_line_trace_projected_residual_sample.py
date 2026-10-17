"""Independent prime-field replay of the trace residual-column sample.

Use quartic evaluations at four fixed offsets to recover each line-cut
quadratic, rather than the producer's directional Taylor expansion.  Replay
the complete prefix of candidate pairs, all reported residual coordinates,
elliptic quotient line identities, and exact/sign-folded column identities.
"""

import argparse
import hashlib
import json
from pathlib import Path

from verify_adaptive_private_columns import add, neg, on_curve

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--sample", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
sample = json.loads(args.sample.read_text())
assert sample["status"] in {
    "verified_p192_old_anchor_trace_projected_residual_sample",
    "verified_p192_all_old_pivots_trace_projected_residual_stratified",
}
stratified = sample["status"].endswith("_stratified")
if stratified:
    assert sample["selection_policy"] == "first_split_nonsource_lines_per_old_pivot"
    assert sample["old_pivots"] == 141
    assert sample["split_lines"] == 141 * sample["per_pivot"]
paths = {key: HERE / name for key, name in {
    "plane": "type_i_plane_quartic.json",
    "old_quartic": "type_i_line_base_points.tsv",
    "new_quartic": "type_i_line_new_points.tsv",
    "collinearity": "type_i_line_full_collinearity.json",
    "old_quotient": "type_i_line_guided_collector.json",
    "new_quotient": "type_i_line_full_new_quotient_points.json",
    "constant": "type_i_line_full_constant.json",
    "base": "type_i_quotient_base_probe.json",
    "full_screen": "type_i_line_nonsource_reuse_old141.json",
}.items()}
assert sample["input_sha256"] == {key: sha(path) for key, path in paths.items()}
plane = json.loads(paths["plane"].read_text())
p = int(plane["p"])
inv2, inv3, inv4, inv8 = (pow(k, -1, p) for k in (2, 3, 4, 8))
terms = [(tuple(map(int, label.split(","))), int(value) % p)
         for label, value in plane["quartic_coefficients"].items()]
old = [tuple(map(int, line.split()[1:])) for line in paths["old_quartic"].read_text().splitlines()]
new = [tuple(map(int, line.split())) for line in paths["new_quartic"].read_text().splitlines()]
points = old + new
assert len(old) == 141 and len(points) == len(set(points)) == 9925
selected_quartic = set(points)
source = set()
for line in json.loads(paths["collinearity"].read_text())["line_records"]:
    ids = list(map(int, line["point_indices_zero_based"]))
    source.update((i, j) for i in ids for j in ids if i < 141 and i < j)
old_q = [tuple(map(int, x)) for x in
         json.loads(paths["old_quotient"].read_text())["base_points_fp"]]
new_q = [tuple(map(int, x)) for x in json.loads(paths["new_quotient"].read_text())]
selected = old_q + new_q
selected_keys = set(selected)
assert len(selected) == len(selected_keys) == 9925
S = tuple(map(int, json.loads(paths["constant"].read_text())["constant_fp"]))
assert on_curve(S, p)


def quartic(x, y):
    return sum(c * pow(x, a, p) * pow(y, b, p)
               for (a, b, _), c in terms) % p


def canonical(point):
    return min(point, neg(point, p)) if point is not None else None


signed_selected_keys = {canonical(q) for q in selected_keys}
rows = sample["rows"]
seen, signed_seen = {}, {}
selected_hits = []
signed_selected_hits = []
repeats = []
signed_repeats = []
next_index = 0
nonsplit = 0
last_pair = rows[-1]["selected_pair"]
for i in range(141 if stratified else last_pair[0] + 1):
    x, y = points[i]
    start_index = next_index
    last_j = (len(points) - 1 if stratified or i != last_pair[0]
              else last_pair[1])
    for j in range(i + 1, last_j + 1):
        if (i, j) in source:
            continue
        u, v = points[j]
        d = (u - x) % p
        assert d
        m = (v - y) * pow(d, -1, p) % p
        def at(t):
            return quartic((x + t) % p, (y + m * t) % p)
        fp1, fm1, fp2, fm2 = at(1), at(-1), at(2), at(-2)
        odd1 = (fp1 - fm1) * inv2 % p
        even1 = (fp1 + fm1) * inv2 % p
        odd2 = (fp2 - fm2) * inv4 % p
        even2 = (fp2 + fm2) * inv8 % p
        c3 = (odd2 - odd1) * inv3 % p
        c1 = (odd1 - c3) % p
        A = (even2 - even1) * inv3 % p
        c2 = (even1 - A) % p
        B = (c3 + d * A) % p
        C = (c2 + d * B) % p
        assert A and C and (c1 + d * C) % p == 0
        discr = (B * B - 4 * A * C) % p
        assert discr
        symbol = pow(discr, (p - 1) // 2, p)
        if symbol == p - 1:
            nonsplit += 1
            continue
        assert symbol == 1
        assert next_index < len(rows)
        row = rows[next_index]
        assert row["selected_pair"] == [i, j]
        sq = pow(discr, (p + 1) // 4, p)
        den = pow(2 * A, -1, p)
        roots = [(-B + sq) * den % p, (-B - sq) * den % p]
        expected = [((x + root) % p, (y + m * root) % p) for root in roots]
        residuals = [tuple(map(int, point)) for point in row["residual_points"]]
        assert expected == residuals and len(set(residuals)) == 2
        assert all(point not in selected_quartic for point in residuals)
        assert all(quartic(*point) == 0 for point in residuals)
        projected = [tuple(map(int, point)) if point is not None else None
                     for point in row["projected_residuals"]]
        assert all(on_curve(point, p) for point in projected)
        total = add(add(selected[i], selected[j], p), projected[0], p)
        assert add(total, projected[1], p) == S
        for side, point in enumerate(projected):
            if point in selected_keys:
                selected_hits.append([next_index, side, point])
            if canonical(point) in signed_selected_keys:
                signed_selected_hits.append([next_index, side, point])
            if point in seen:
                repeats.append([next_index, side, seen[point]])
            seen.setdefault(point, [next_index, side])
            signed = canonical(point)
            if signed in signed_seen:
                signed_repeats.append([next_index, side, signed_seen[signed]])
            signed_seen.setdefault(signed, [next_index, side])
        next_index += 1
        if stratified and next_index - start_index == sample["per_pivot"]:
            break
    if stratified:
        assert next_index - start_index == sample["per_pivot"]

assert next_index == len(rows) == sample["split_lines"]
assert len(seen) == sample["distinct_projected_columns"]
assert len(signed_seen) == sample["distinct_signed_projected_columns"]
assert selected_hits == sample["selected_column_hits_first_100"]
assert repeats[:100] == sample["repeated_columns_first_100"]
assert signed_repeats[:100] == sample["signed_repeated_columns_first_100"]
assert sample["identity_columns"] == sum(point is None for point in seen)
assert sample["direct_partner_controls"] >= 128
out = {
    "status": "verified_independent_trace_residual_sample_replay",
    "split_lines": next_index,
    "nonsplit_prefix_lines": nonsplit,
    "projected_columns": 2 * next_index,
    "distinct_projected_columns": len(seen),
    "distinct_signed_projected_columns": len(signed_seen),
    "selected_hits": len(selected_hits),
    "signed_selected_hits": len(signed_selected_hits),
    "repeated_columns": len(repeats),
    "signed_repeated_columns": len(signed_repeats),
    "sample_sha256": sha(args.sample),
    "source_sha256": sha(__file__),
    "boundary": "Independent interpolation and prime-field group replay of the complete frozen sample prefix; subgroup projection equality relies on the paired native trace/direct controls.",
}
args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print("TRACE_RESIDUAL_REPLAY_PASS", json.dumps({
    "split_lines": next_index,
    "projected_columns": 2 * next_index,
    "selected_hits": len(selected_hits),
    "signed_selected_hits": len(signed_selected_hits),
    "repeated_columns": len(repeats),
    "signed_repeated_columns": len(signed_repeats),
}, sort_keys=True))
