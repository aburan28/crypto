"""Independent prime-field replay of one complete projected pivot chunk."""

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from verify_adaptive_private_columns import add, neg, on_curve

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--chunk", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
chunk = json.loads(args.chunk.read_text())
assert chunk["status"] == "verified_p192_full_old_pivot_trace_projection_chunk"
i = chunk["pivot"]
assert 0 <= i < 141
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
assert chunk["input_sha256"] == {key: sha(path) for key, path in paths.items()}
plane = json.loads(paths["plane"].read_text())
p = int(plane["p"])
assert p == int(chunk["field_prime"]) and p % 4 == 3
inv2, inv3, inv4, inv8 = (pow(k, -1, p) for k in (2, 3, 4, 8))
terms = [(tuple(map(int, label.split(","))), int(value) % p)
         for label, value in plane["quartic_coefficients"].items()]
old = [tuple(map(int, line.split()[1:])) for line in paths["old_quartic"].read_text().splitlines()]
new = [tuple(map(int, line.split())) for line in paths["new_quartic"].read_text().splitlines()]
points = old + new
selected_geom = set(points)
assert len(old) == 141 and len(points) == len(selected_geom) == 9925
source = set()
for record in json.loads(paths["collinearity"].read_text())["line_records"]:
    ids = list(map(int, record["point_indices_zero_based"]))
    source.update(j for x in ids for j in ids if x == i and j > i)
old_q = [tuple(map(int, x)) for x in
         json.loads(paths["old_quotient"].read_text())["base_points_fp"]]
new_q = [tuple(map(int, x)) for x in json.loads(paths["new_quotient"].read_text())]
selected = old_q + new_q
selected_keys = set(selected)
signed_selected = {min(Q, neg(Q, p)) for Q in selected}
assert len(selected_keys) == len(selected) == 9925
S = tuple(map(int, json.loads(paths["constant"].read_text())["constant_fp"]))
assert on_curve(S, p)


def quartic(x, y):
    return sum(c * pow(x, a, p) * pow(y, b, p)
               for (a, b, _), c in terms) % p


rows = chunk["rows"]
counts = Counter()
seen = {}
signed_seen = {}
geom_seen = set()
stream = hashlib.sha256()
next_index = 0
x, y = points[i]
assert quartic(x, y) == 0
for j in range(i + 1, len(points)):
    counts["selected_pairs"] += 1
    if j in source:
        counts["source_presentations"] += 1
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
    if pow(discr, (p - 1) // 2, p) == p - 1:
        counts["nonsplit_lines"] += 1
        continue
    assert pow(discr, (p - 1) // 2, p) == 1
    assert next_index < len(rows)
    row = rows[next_index]
    assert row["selected_pair"] == [i, j]
    sq = pow(discr, (p + 1) // 4, p)
    den = pow(2 * A, -1, p)
    roots = [(-B + sq) * den % p, (-B - sq) * den % p]
    residuals = [((x + root) % p, (y + m * root) % p) for root in roots]
    reported = [tuple(map(int, q)) for q in row["residual_points"]]
    assert reported == residuals and len(set(reported)) == 2
    assert all(Q not in selected_geom and quartic(*Q) == 0 for Q in reported)
    projected = [tuple(map(int, Q)) if Q is not None else None
                 for Q in row["projected_residuals"]]
    assert all(on_curve(Q, p) for Q in projected)
    assert add(add(add(selected[i], selected[j], p), projected[0], p),
               projected[1], p) == S
    for point, key in zip(reported, projected):
        encoded = point[0].to_bytes(24, "big") + point[1].to_bytes(24, "big")
        stream.update(encoded)
        assert encoded not in geom_seen
        geom_seen.add(encoded)
        counts["projected_columns"] += 1
        counts["identity_columns"] += key is None
        counts["selected_hits"] += key in selected_keys
        signed = min(key, neg(key, p)) if key is not None else None
        counts["signed_selected_hits"] += signed in signed_selected
        counts["repeats"] += key in seen
        seen.setdefault(key, next_index)
        counts["signed_repeats"] += signed in signed_seen
        signed_seen.setdefault(signed, next_index)
    next_index += 1
    counts["split_lines"] += 1

assert next_index == len(rows) == chunk["counts"]["split_lines"]
assert all(chunk["counts"].get(key, 0) == value for key, value in counts.items())
assert chunk["counts"]["direct_partner_controls"] == (next_index - 1) // 1024 + 1
assert chunk["distinct_projected_columns"] == len(seen)
assert chunk["distinct_signed_projected_columns"] == len(signed_seen)
assert chunk["residual_stream_sha256"] == stream.hexdigest()
out = {
    "status": "verified_independent_p192_trace_projected_pivot_replay",
    "pivot": i,
    "counts": dict(counts),
    "distinct_projected_columns": len(seen),
    "distinct_signed_projected_columns": len(signed_seen),
    "residual_stream_sha256": stream.hexdigest(),
    "chunk_sha256": sha(args.chunk),
    "source_sha256": sha(Path(__file__)),
    "boundary": "Complete independent prime-field quartic interpolation and line-group replay for one pivot; subgroup projection correctness rests on the paired native trace and partner controls.",
}
args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print("TRACE_PIVOT_REPLAY_PASS", json.dumps({
    "pivot": i,
    "split_lines": next_index,
    "projected_columns": counts["projected_columns"],
    "repeats": counts["repeats"],
    "signed_repeats": counts["signed_repeats"],
}, sort_keys=True))
