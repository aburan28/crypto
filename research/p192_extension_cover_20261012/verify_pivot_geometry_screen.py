"""Independent polynomial-expansion controls for the pivot geometry screen."""

import argparse
import hashlib
import json
from collections import Counter
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--pivots", type=int, default=32)
args = parser.parse_args()
screen_path = HERE / f"type_i_pivot_geometry_screen_{args.pivots:03d}.json"
screen = json.loads(screen_path.read_text())
assert screen["status"] == "complete_declared_pivot_residual_geometry_screen"
assert screen["pivot_fans"] == len(screen["rows"]) == args.pivots
source = HERE / "line_all_pivot_geometry_screen.py"
assert screen["source_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
for name, digest in screen["input_sha256"].items():
    assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == digest
for name, digest in screen["chunk_input_sha256"].items():
    assert hashlib.sha256((HERE / "type_i_line_trace_full_chunks" / name).read_bytes()).hexdigest() == digest

plane = json.loads((HERE / "type_i_plane_quartic.json").read_text())
p = int(plane["p"])
assert p == int(screen["field_prime"]) and p % 4 == 3
terms = [(tuple(map(int, key.split(","))), int(value) % p)
         for key, value in plane["quartic_coefficients"].items()]
selected = [tuple(map(int, row.split()[1:])) for row in
            (HERE / "type_i_line_base_points.tsv").read_text().splitlines()]
selected += [tuple(map(int, row.split())) for row in
             (HERE / "type_i_line_new_points.tsv").read_text().splitlines()]
assert len(selected) == screen["selected_classes"] == 9925
selected_index_by_point = {point: index for index, point in enumerate(selected)}
assert len(selected_index_by_point) == len(selected)


def quartic(x, y):
    return sum(coefficient * pow(x, a, p) * pow(y, b, p)
               for (a, b, _), coefficient in terms) % p


def line_coefficients(pivot, slope):
    """Expand F(P+z(1,m)) directly by multiplying two binomial arrays."""
    x, y = pivot
    out = [0] * 5
    for (a, b, _), coefficient in terms:
        left = [comb(a, i) * pow(x, a-i, p) % p for i in range(a+1)]
        right = [comb(b, j) * pow(y, b-j, p) * pow(slope, j, p) % p
                 for j in range(b+1)]
        for i, li in enumerate(left):
            for j, rj in enumerate(right):
                out[i+j] = (out[i+j] + coefficient * li * rj) % p
    return out


def case(pivot_index, selected_index, side=None):
    pivot = tuple(map(int, screen["rows"][pivot_index]["pivot"]))
    u, v = selected[selected_index]
    d = (u-pivot[0]) % p
    if not d:
        return "vertical", None
    slope = (v-pivot[1]) * pow(d, -1, p) % p
    c = line_coefficients(pivot, slope)
    assert c[0] == 0 and sum(c[i]*pow(d, i, p) for i in range(5)) % p == 0
    A = c[4]
    if not A:
        return "degree_below_four", None
    B = (c[3]+d*A) % p
    C = (c[2]+d*B) % p
    if not C or (c[1]+d*C) % p:
        return "tangent_or_invalid", None
    discriminant = (B*B-4*A*C) % p
    if not discriminant:
        return "double_residual", None
    if pow(discriminant, (p-1)//2, p) == p-1:
        return "nonsplit", None
    root = pow(discriminant, (p+1)//4, p)
    assert root*root % p == discriminant
    assert side in (None, 0, 1)
    roots = [(-B+root) * pow(2*A, -1, p) % p,
             (-B-root) * pow(2*A, -1, p) % p]
    points = [((pivot[0]+z) % p, (pivot[1]+slope*z) % p) for z in roots]
    assert all(quartic(*point) == 0 for point in points)
    return "split", points if side is None else points[side]


totals = Counter()
source_selected_indices = []
sampled_line_cases = 0
for row in screen["rows"]:
    pivot_index = row["pivot_index"]
    assert 0 <= pivot_index < args.pivots
    pivot = tuple(map(int, row["pivot"]))
    partner = tuple(map(int, row["partner"]))
    assert quartic(*pivot) == quartic(*partner) == 0
    source_chunk = json.loads((HERE / "type_i_line_trace_full_chunks" /
                               f"pivot-{pivot_index:03d}.json").read_text())
    source_i, source_j = source_chunk["rows"][0]["selected_pair"]
    assert source_i == pivot_index
    source_selected_indices.append(source_j)
    assert row["counts"].get("source_fiber_presentations") == 2
    assert row["counts"].get("within_fan_repeats") == 1
    assert set(case(pivot_index, source_i)[1]) == {partner, selected[source_j]}
    assert set(case(pivot_index, source_j)[1]) == {partner, selected[source_i]}
    assert sum(row["counts"].get(key, 0) for key in
               ("vertical", "degree_below_four", "tangent_or_invalid",
                "double_residual", "nonsplit", "split")) == 9925
    totals.update(row["counts"])
    sample_indices = {0, 1, 2, 17, 127, 1024, 4096, 8192, 9924,
                      (pivot_index*71) % 9925}
    sampled_line_cases += len(sample_indices)
    for j in sample_indices:
        status, _ = case(pivot_index, j)
        assert status in row["counts"]
assert dict(totals) == screen["counts"]
assert (screen["distinct_residual_points"] ==
        2*totals["split"]-totals["within_fan_repeats"]-totals["cross_fan_repeats"])
assert len(screen["cross_fan_witnesses_first_1024"]) == min(
    totals["cross_fan_repeats"], 1024)
for witness in screen["cross_fan_witnesses_first_1024"]:
    point = tuple(map(int, witness["point"]))
    assert quartic(*point) == 0
    for pivot_index, selected_index, side in (witness["first"], witness["second"]):
        assert case(pivot_index, selected_index, side) == ("split", point)
    assert witness["first"][0] != witness["second"][0]
    assert selected_index_by_point[point] in source_selected_indices
    for pivot_index, selected_i, _ in (witness["first"], witness["second"]):
        assert selected_i == pivot_index
        assert selected_index_by_point[point] == source_selected_indices[pivot_index]

assert totals["within_fan_repeats"] == args.pivots
assert totals["source_fiber_presentations"] == 2*args.pivots
assert totals["cross_fan_repeats"] == len(screen["cross_fan_witnesses_first_1024"])

out = {
    "status": "independent_pivot_geometry_polynomial_controls_pass",
    "pivot_fans": args.pivots,
    "sampled_line_cases": sampled_line_cases,
    "cross_fan_witnesses_checked": len(screen["cross_fan_witnesses_first_1024"]),
    "cross_fan_repeats_on_original_source_lines": totals["cross_fan_repeats"],
    "cross_fan_repeated_selected_indices": sorted({
        selected_index_by_point[tuple(map(int, w["point"]))]
        for w in screen["cross_fan_witnesses_first_1024"]}),
    "non_source_geometric_repeats": 0,
    "counts": dict(totals),
    "screen_sha256": hashlib.sha256(screen_path.read_bytes()).hexdigest(),
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "boundary": "Independent quartic polynomial expansion for stratified line cases and every stored cross-fan witness; all detected repeats are source-partner or original selected-source incidences. This is a geometric statement; projected subgroup collisions require separate tests. The full ordered stream digest is producer-only.",
}
path = HERE / f"type_i_pivot_geometry_screen_{args.pivots:03d}_verification.json"
path.write_text(json.dumps(out, indent=2, sort_keys=True)+"\n")
print("PIVOT_GEOMETRY_REPLAY_PASS", args.pivots,
      out["cross_fan_witnesses_checked"], flush=True)
