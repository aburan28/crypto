"""Independent prime-field replay of the complete residual-pivot pilot."""

import hashlib
import json
from collections import Counter
from pathlib import Path

from verify_adaptive_private_columns import add, neg, on_curve


HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
receipt_path = HERE / "type_i_residual_pivot000_three.json"
receipt = json.loads(receipt_path.read_text())
assert receipt["status"] == "p192_residual_pivot000_three_step_projected_fan"
assert receipt["source_sha256"] == sha(HERE / "line_residual_pivot_three.py")
assert all(sha(HERE / name) == digest
           for name, digest in receipt["input_sha256"].items())
assert receipt["runtime_info_sha256"] == sha(
    HERE / "type_i_residual_pivot000_three_runtime_info.json")
plane = json.loads((HERE / "type_i_plane_quartic.json").read_text())
p = int(plane["p"])
assert p == int(receipt["field_prime"])
terms = [(tuple(map(int, label.split(","))), int(value) % p)
         for label, value in plane["quartic_coefficients"].items()]
old = [tuple(map(int, line.split()[1:])) for line in
       (HERE / "type_i_line_base_points.tsv").read_text().splitlines()]
new = [tuple(map(int, line.split())) for line in
       (HERE / "type_i_line_new_points.tsv").read_text().splitlines()]
points = old + new
assert len(points) == len(set(points)) == 9925
selected_geometry = set(points)
selected = ([tuple(map(int, x)) for x in json.loads(
    (HERE / "type_i_line_guided_collector.json").read_text())["base_points_fp"]]
    + [tuple(map(int, x)) for x in json.loads(
        (HERE / "type_i_line_full_new_quotient_points.json").read_text())])
assert len(selected) == 9925 and all(on_curve(q, p) for q in selected)
selected_keys = set(selected)
selected_signed = {min(q, neg(q, p)) for q in selected}
S = tuple(map(int, json.loads(
    (HERE / "type_i_line_full_constant.json").read_text())["constant_fp"]))
pivot = tuple(map(int, receipt["pivot_quartic_point"]))
pivot_q = tuple(map(int, receipt["pivot_quotient_point"]))
assert on_curve(S, p) and on_curve(pivot_q, p)
assert pivot not in selected_geometry
inv2, inv3, inv4, inv8 = (pow(k, -1, p) for k in (2, 3, 4, 8))


def quartic(x, y):
    return sum(c * pow(x, a, p) * pow(y, b, p)
               for (a, b, _), c in terms) % p


assert quartic(*pivot) == 0
rows = receipt["rows"]
next_index = 0
counts = Counter()
seen_geometry, seen_exact, seen_signed = {}, {}, {}
for j, (u, v) in enumerate(points):
    counts["selected_pairs"] += 1
    d = (u - pivot[0]) % p
    assert d
    m = (v - pivot[1]) * pow(d, -1, p) % p

    def at(t):
        return quartic((pivot[0] + t) % p, (pivot[1] + m * t) % p)

    fp1, fm1, fp2, fm2 = at(1), at(-1), at(2), at(-2)
    odd1, even1 = (fp1 - fm1) * inv2 % p, (fp1 + fm1) * inv2 % p
    odd2, even2 = (fp2 - fm2) * inv4 % p, (fp2 + fm2) * inv8 % p
    c3 = (odd2 - odd1) * inv3 % p
    c1 = (odd1 - c3) % p
    A = (even2 - even1) * inv3 % p
    c2 = (even1 - A) % p
    B, C = (c3 + d * A) % p, (c2 + d * (c3 + d * A)) % p
    assert A and C and (c1 + d * C) % p == 0
    discr = (B * B - 4 * A * C) % p
    assert discr
    if pow(discr, (p - 1) // 2, p) == p - 1:
        counts["nonsplit"] += 1
        continue
    assert pow(discr, (p - 1) // 2, p) == 1
    assert next_index < len(rows)
    row = rows[next_index]
    assert row["selected_index"] == j
    sq, den = pow(discr, (p + 1) // 4, p), pow(2 * A, -1, p)
    roots = ((-B + sq) * den % p, (-B - sq) * den % p)
    expected = [((pivot[0] + z) % p, (pivot[1] + m * z) % p)
                for z in roots]
    reported = [tuple(map(int, q)) for q in row["residual_points"]]
    assert reported == expected and reported[0] != reported[1]
    assert all(quartic(*q) == 0 for q in reported)
    projected = [tuple(map(int, q)) if q is not None else None
                 for q in row["projected_residuals"]]
    assert all(on_curve(q, p) for q in projected)
    total = add(add(add(pivot_q, selected[j], p), projected[0], p),
                projected[1], p)
    assert total == S
    for side, (point, key) in enumerate(zip(reported, projected)):
        if point in selected_geometry:
            counts["geometric_selected_hits"] += 1
        if point in seen_geometry:
            counts["geometric_repeats"] += 1
        seen_geometry.setdefault(point, [next_index, side])
        if key is None:
            counts["identity_columns"] += 1
        if key in selected_keys:
            counts["selected_hits"] += 1
        signed = min(key, neg(key, p)) if key is not None else None
        if signed in selected_signed:
            counts["signed_selected_hits"] += 1
        if key in seen_exact:
            counts["exact_repeats"] += 1
        if signed in seen_signed:
            counts["signed_repeats"] += 1
        seen_exact.setdefault(key, [next_index, side])
        seen_signed.setdefault(signed, [next_index, side])
    counts["split"] += 1
    next_index += 1

assert next_index == len(rows) == 5037
assert dict(counts) == receipt["counts"]
assert len(seen_geometry) == receipt["distinct_geometry"]
assert len(seen_exact) == receipt["distinct_exact_columns"]
assert len(seen_signed) == receipt["distinct_signed_columns"]
assert receipt["six_step_control_row_indices"] == list(range(0, 5037, 256))
old_exact_hits, old_signed_hits = [], []
old_columns = 0
for i in range(141):
    old_chunk = json.loads((HERE / f"type_i_line_trace_full_chunks/pivot-{i:03d}.json").read_text())
    for row_index, row in enumerate(old_chunk["rows"]):
        for side, pair in enumerate(row["projected_residuals"]):
            key = tuple(map(int, pair))
            old_columns += 1
            if key in seen_exact:
                old_exact_hits.append([i, row_index, side, seen_exact[key]])
            signed = min(key, neg(key, p))
            if signed in seen_signed:
                old_signed_hits.append([i, row_index, side, seen_signed[signed]])
assert old_columns == receipt["old_family_columns_scanned"] == 1360948
assert old_exact_hits == receipt["old_exact_hits"]
assert old_signed_hits == receipt["old_signed_hits"]
assert old_exact_hits == old_signed_hits == [[0, 0, 1, [0, 0]]]

out = {
    "status": "verified_independent_residual_pivot000_full_fan",
    "split_lines": next_index,
    "distinct_projected_columns": len(seen_exact),
    "old_family_exact_hits": old_exact_hits,
    "old_family_signed_hits": old_signed_hits,
    "counts": dict(counts),
    "receipt_sha256": sha(receipt_path),
    "source_sha256": sha(Path(__file__)),
    "boundary": "Independent quartic intersections, subgroup curve equations, line sums, and old-family column membership for the complete one-pivot fan; native trace controls establish sampled projection equivalence.",
}
target = HERE / "type_i_residual_pivot000_three_replay.json"
target.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print("RESIDUAL_PIVOT_THREE_REPLAY_PASS", next_index, old_exact_hits, flush=True)
