"""Screen all 141 residual-pivot fans for repeated quartic geometry."""

import hashlib
import json
import time
import argparse
from collections import Counter
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--pivots", type=int, default=141)
args = parser.parse_args()
assert 1 <= args.pivots <= 141
plane_path = HERE / "type_i_plane_quartic.json"
old_path = HERE / "type_i_line_base_points.tsv"
new_path = HERE / "type_i_line_new_points.tsv"
chunks = HERE / "type_i_line_trace_full_chunks"
plane = json.loads(plane_path.read_text())
p = int(plane["p"])
assert p % 4 == 3
selected = [tuple(map(int, row.split()[1:])) for row in old_path.read_text().splitlines()]
selected += [tuple(map(int, row.split())) for row in new_path.read_text().splitlines()]
assert len(selected) == len(set(selected)) == 9925
terms = [(tuple(map(int, label.split(","))), int(value) % p)
         for label, value in plane["quartic_coefficients"].items()]


def taylor_at(x, y):
    xp = [pow(x, i, p) for i in range(5)]
    yp = [pow(y, i, p) for i in range(5)]
    out = [[0]*(k+1) for k in range(5)]
    for (a, b, _), value in terms:
        for ix in range(a+1):
            for iy in range(b+1):
                out[ix+iy][iy] = (out[ix+iy][iy] + value*comb(a, ix)*comb(b, iy)
                                   * xp[a-ix]*yp[b-iy]) % p
    assert out[0][0] == 0
    return out


def eval_poly(coefficients, x):
    result = 0
    for value in reversed(coefficients):
        result = (result*x + value) % p
    return result


start = time.monotonic_ns()
seen = {}
rows = []
cross_witnesses = []
totals = Counter()
geometry_hash = hashlib.sha256()
source_hashes = {}
for pivot_index in range(args.pivots):
    path = chunks / f"pivot-{pivot_index:03d}.json"
    source_hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    first = json.loads(path.read_text())["rows"][0]
    assert first["selected_pair"][0] == pivot_index
    pivot = tuple(map(int, first["residual_points"][0]))
    partner = tuple(map(int, first["residual_points"][1]))
    taylor = taylor_at(*pivot)
    count = Counter()
    denominators = [(u-pivot[0]) % p for u, _ in selected]
    nonzero = [j for j, d in enumerate(denominators) if d]
    prefixes = []
    product = 1
    for j in nonzero:
        product = product * denominators[j] % p
        prefixes.append(product)
    inverse_product = pow(product, -1, p)
    inverse_denominators = {}
    for k in range(len(nonzero)-1, -1, -1):
        j = nonzero[k]
        inverse_denominators[j] = inverse_product * (prefixes[k-1] if k else 1) % p
        inverse_product = inverse_product * denominators[j] % p
    assert inverse_product == 1
    for selected_index, (u, v) in enumerate(selected):
        d = denominators[selected_index]
        if not d:
            count["vertical"] += 1
            continue
        m = (v-pivot[1])*inverse_denominators[selected_index] % p
        c1, c2, c3, A = (eval_poly(taylor[k], m) for k in range(1, 5))
        if not A:
            count["degree_below_four"] += 1
            continue
        B, C = (c3+d*A) % p, (c2+d*(c3+d*A)) % p
        if not C or (c1+d*C) % p:
            count["tangent_or_invalid"] += 1
            continue
        disc = (B*B-4*A*C) % p
        if not disc:
            count["double_residual"] += 1
            continue
        if pow(disc, (p-1)//2, p) == p-1:
            count["nonsplit"] += 1
            continue
        count["split"] += 1
        sq, den = pow(disc, (p+1)//4, p), pow(2*A, -1, p)
        for side, z in enumerate(((-B+sq)*den % p, (-B-sq)*den % p)):
            point = (pivot[0]+z) % p, (pivot[1]+m*z) % p
            label = (pivot_index, selected_index, side)
            geometry_hash.update(pivot_index.to_bytes(2, "big"))
            geometry_hash.update(selected_index.to_bytes(2, "big"))
            geometry_hash.update(bytes([side]))
            geometry_hash.update(point[0].to_bytes(24, "big"))
            geometry_hash.update(point[1].to_bytes(24, "big"))
            prior = seen.get(point)
            if prior is None:
                seen[point] = label
            elif prior[0] == pivot_index:
                count["within_fan_repeats"] += 1
            else:
                count["cross_fan_repeats"] += 1
                if len(cross_witnesses) < 1024:
                    cross_witnesses.append({
                        "point": [str(point[0]), str(point[1])],
                        "first": list(prior), "second": list(label),
                    })
            if point in (pivot, partner):
                count["source_fiber_presentations"] += 1
    assert sum(count[k] for k in ("vertical", "degree_below_four",
                                  "tangent_or_invalid", "double_residual",
                                  "nonsplit", "split")) == len(selected)
    rows.append({"pivot_index": pivot_index, "pivot": list(map(str, pivot)),
                 "partner": list(map(str, partner)), "counts": dict(count)})
    totals.update(count)
    if (pivot_index+1) % 10 == 0:
        print("ALL_PIVOT_GEOMETRY_PROGRESS", pivot_index+1, args.pivots,
              "cross", totals["cross_fan_repeats"], flush=True)

output = {
    "status": "complete_declared_pivot_residual_geometry_screen",
    "field_prime": str(p), "selected_classes": len(selected),
    "pivot_fans": len(rows), "counts": dict(totals),
    "distinct_residual_points": len(seen),
    "cross_fan_witnesses_first_1024": cross_witnesses,
    "rows": rows,
    "geometry_presentation_sha256": geometry_hash.hexdigest(),
    "elapsed_ns_unisolated": time.monotonic_ns()-start,
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "input_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                     for path in (plane_path, old_path, new_path)},
    "chunk_input_sha256": source_hashes,
    "boundary": "Complete rational split-line geometry for the declared pivot prefix. Projected subgroup-column collisions and useful rank require separate checks; geometry alone cannot certify their absence.",
}
path = HERE / f"type_i_pivot_geometry_screen_{args.pivots:03d}.json"
path.write_text(json.dumps(output, indent=2, sort_keys=True)+"\n")
print("ALL_PIVOT_GEOMETRY_PASS", len(rows), len(seen),
      totals["cross_fan_repeats"], flush=True)
