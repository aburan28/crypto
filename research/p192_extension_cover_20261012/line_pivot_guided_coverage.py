"""Exact three-sum coverage for reduced pivot classes and known-scalar laws."""

import hashlib
import json
import time
from collections import Counter
from pathlib import Path

from sage.all import EllipticCurve, GF, Integer, matrix

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


reduction_path = HERE / "type_i_two_pivot_reduced64.json"
fixture_path = HERE / "type_i_base_points_legacy_fixture.json"
guide_path = HERE / "type_i_line_guided_collector.json"
density_path = HERE / "type_i_character_density.json"
reduction = json.loads(reduction_path.read_text())
fixture = json.loads(fixture_path.read_text())
guide = json.loads(guide_path.read_text())
density = json.loads(density_path.read_text())
p, r = Integer(reduction["field_prime"]), Integer(reduction["subgroup_order"])
assert p == Integer(fixture["field_prime"]) == Integer(guide["field_prime"])
assert r == Integer(fixture["subgroup_order"]) == Integer(guide["subgroup_order"])
assert len(reduction["representatives"]) == 64
assert all(row["status"] == "verified_fixed_origin_degree_three"
           for row in reduction["representatives"])
F = GF(p)
E = EllipticCurve(F, [F(-3), F(Integer(
    "64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1", 16))])
G = E(F(Integer("188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012", 16)),
      F(Integer("07192B95FFC8DA78631011ED6B24CDD573F977A11E794811", 16)))
O = E(0)
assert r * G == O


def checked_point(x, y):
    x, y = int(x), int(y)
    assert (y*y - x*x*x + 3*x - int(E.a6())) % int(p) == 0
    return E.point([F(x), F(y), F(1)], check=False)


old = [checked_point(item["x"], item["y"]) for item in fixture["points"]]
fan = [checked_point(*row["quotient"]) for row in reduction["rows"]]
assert len(old) == 16 and len(fan) == 64
assert len(set(old + fan)) == 80
assert all(r * point == O for point in fan)
coefficients = [F(int(value)) for value in density["sextic_coefficients_decimal"]]
inv2 = Integer(pow(2, -1, int(r)))


def char(point):
    assert point != O
    x = point[0]
    value = F(0)
    for coefficient in reversed(coefficients):
        value = value*x + coefficient
    assert value != 0
    return "+" if value.is_square() else "-"


def cell(query):
    pre = -inv2*query
    assert -2*pre == query and pre != O
    return char(pre) + char(2*pre)


def key(point):
    return "O" if point == O else f"{int(point[0])}:{int(point[1])}"


def support(base):
    first = {}
    collisions = []
    for i in range(len(base)):
        for j in range(i, len(base)):
            pair = base[i] + base[j]
            for k in range(j, len(base)):
                target = pair + base[k]
                name = key(target)
                if name in first:
                    collisions.append((first[name][0], (i, j, k)))
                else:
                    first[name] = ((i, j, k), target)
    return first, collisions


start = time.monotonic_ns()
cases = {}
for name, base in (
    ("old16", old),
    ("old16_fan000_16", old + fan[:16]),
    ("old16_two_fans_16_each", old + fan[:16] + fan[32:48]),
):
    targets, collisions = support(base)
    cells = Counter()
    for witness, target in targets.values():
        if target != O:
            cells[cell(target)] += 1
    equations = []
    for first, second in collisions:
        row = [0] * len(base)
        for index in first:
            row[index] += 1
        for index in second:
            row[index] -= 1
        equations.append(row)
    rank = int(matrix(GF(r), equations).rank()) if equations else 0
    cases[name] = {
        "base_size": len(base),
        "unordered_triples": len(targets) + len(collisions),
        "distinct_targets": len(targets),
        "duplicates": len(collisions),
        "collision_rank": rank,
        "character_cells": dict(sorted(cells.items())),
        "identity_target": key(O) in targets,
        "base_quotient_points": [[str(int(Q[0])), str(int(Q[1]))] for Q in base],
        "target_witness_by_key": {name: list(witness)
                                  for name, (witness, _) in targets.items()},
    }
    print("PIVOT_SUPPORT_PROGRESS", name, len(targets), dict(cells), flush=True)

query_laws = {}
for name, scalars in (
    ("uniform", guide["uniform_scalars"]),
    ("character_guided", guide["guided_scalars"]),
):
    hits = {}
    observed_cells = Counter()
    for index, encoded in enumerate(scalars):
        a = Integer(encoded)
        assert 0 < a < r
        target = -2*a*G
        assert target != O and target == a*(-2*G)
        observed_cells[cell(target)] += 1
        for case_name, case in cases.items():
            witness = case["target_witness_by_key"].get(key(target))
            if witness is not None:
                hits.setdefault(case_name, []).append({
                    "query_index": index, "scalar": str(a),
                    "witness": witness, "target_key": key(target),
                })
    if name == "character_guided":
        assert observed_cells == {"++": 64}
    query_laws[name] = {"queries": len(scalars),
                        "cells": dict(sorted(observed_cells.items())),
                        "hits": hits, "scalar_sha256": hashlib.sha256(
                            "\n".join(map(str, scalars)).encode()).hexdigest()}

for case in cases.values():
    del case["target_witness_by_key"]
output = {
    "status": "verified_two_pivot_geometry_guided_three_sum_coverage",
    "candidate_id": None,
    "field_prime": str(p), "subgroup_order": str(r),
    "cases": cases, "known_scalar_query_laws": query_laws,
    "elapsed_ns_unisolated": time.monotonic_ns() - start,
    "source_sha256": sha(__file__),
    "input_sha256": {path.name: sha(path) for path in
                     (reduction_path, fixture_path, guide_path, density_path)},
    "boundary": "Exact three-sum support and character cells for the declared 16/32/48-class bases. Frozen known-scalar query samples retain q=a(-2G); full pivot fans, independent divisor replay, and relation-stage costs are separate gates.",
}
path = HERE / "type_i_two_pivot_guided_coverage.json"
path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print("PIVOT_GUIDED_COVERAGE_PASS", {name: {
    "support": case["distinct_targets"], "admitted": case["character_cells"].get("++", 0),
    "rank": case["collision_rank"]} for name, case in cases.items()}, flush=True)
