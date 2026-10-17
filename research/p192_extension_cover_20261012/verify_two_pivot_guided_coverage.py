"""Independent integer-curve replay of two-pivot support and scalar queries."""

import hashlib
import json
from collections import Counter
from pathlib import Path

from verify_adaptive_private_columns import (
    GX, GY, add, cell, key, mul as affine_mul, neg, on_curve,
)

HERE = Path(__file__).resolve().parent
path = HERE / "type_i_two_pivot_guided_coverage.json"
guide_path = HERE / "type_i_line_guided_collector.json"
density_path = HERE / "type_i_character_density.json"
result = json.loads(path.read_text())
guide = json.loads(guide_path.read_text())
density = json.loads(density_path.read_text())
p, r = int(result["field_prime"]), int(result["subgroup_order"])
assert p == int(guide["field_prime"]) == int(density["p"])
assert r == int(guide["subgroup_order"]) == int(density["r"])
assert result["source_sha256"] == hashlib.sha256(
    (HERE / "line_pivot_guided_coverage.py").read_bytes()).hexdigest()
coefficients = [int(value) for value in density["sextic_coefficients_decimal"]]
generator = GX, GY
assert on_curve(generator, p) and affine_mul(r, generator, p) is None
inv2 = pow(2, -1, r)


def mul(scalar, point, modulus):
    """Independent mixed-Jacobian scalar multiply with one final inversion."""
    if point is None or scalar == 0:
        return None
    x, y = point
    acc = None
    for digit in bin(scalar)[2:]:
        if acc is not None:
            X, Y, Z = acc
            if Y == 0:
                acc = None
            else:
                y2 = Y*Y % modulus
                z2 = Z*Z % modulus
                s = 4*X*y2 % modulus
                m = 3*(X-z2)*(X+z2) % modulus
                acc = ((m*m-2*s) % modulus,
                       (m*(s-(m*m-2*s))-8*y2*y2) % modulus,
                       2*Y*Z % modulus)
        if digit == "1":
            if acc is None:
                acc = (x, y, 1)
            else:
                X, Y, Z = acc
                z2 = Z*Z % modulus
                u = x*z2 % modulus
                v = y*Z*z2 % modulus
                h, rr = (u-X) % modulus, (v-Y) % modulus
                if h == 0:
                    if rr == 0:
                        y2 = Y*Y % modulus
                        s = 4*X*y2 % modulus
                        m = 3*(X-z2)*(X+z2) % modulus
                        acc = ((m*m-2*s) % modulus,
                               (m*(s-(m*m-2*s))-8*y2*y2) % modulus,
                               2*Y*Z % modulus)
                    else:
                        acc = None
                else:
                    hh = h*h % modulus
                    hhh = h*hh % modulus
                    vv = X*hh % modulus
                    xx = (rr*rr-hhh-2*vv) % modulus
                    acc = (xx, (rr*(vv-xx)-Y*hhh) % modulus, Z*h % modulus)
    if acc is None:
        return None
    X, Y, Z = acc
    inverse = pow(Z, -1, modulus)
    inverse2 = inverse*inverse % modulus
    return X*inverse2 % modulus, Y*inverse2*inverse % modulus


for test_scalar in (1, 2, 3, 17, 123456789, int(r)):
    assert mul(test_scalar, generator, p) == affine_mul(test_scalar, generator, p)


def rank(rows, n):
    pivots = {}
    for row in rows:
        v = row[:]
        for j in sorted(pivots):
            if v[j]:
                factor = v[j]
                v = [(x - factor*y) % r for x, y in zip(v, pivots[j])]
        lead = next((j for j, x in enumerate(v) if x), None)
        if lead is not None:
            inverse = pow(v[lead], -1, r)
            pivots[lead] = [x*inverse % r for x in v]
    return len(pivots)


support_by_case = {}
for name, case in result["cases"].items():
    base = [tuple(map(int, point)) for point in case["base_quotient_points"]]
    assert len(base) == case["base_size"] and len(set(base)) == len(base)
    assert all(on_curve(point, p) for point in base)
    first = {}
    collisions = []
    for i in range(len(base)):
        for j in range(i, len(base)):
            pair = add(base[i], base[j], p)
            for k in range(j, len(base)):
                target = add(pair, base[k], p)
                name_key = key(target)
                if name_key in first:
                    collisions.append((first[name_key][0], (i, j, k)))
                else:
                    first[name_key] = ((i, j, k), target)
    assert len(first) == case["distinct_targets"]
    assert len(collisions) == case["duplicates"]
    assert len(first) + len(collisions) == case["unordered_triples"]
    assert ("O" in first) == case["identity_target"]
    equations = []
    for a, b in collisions:
        row = [0]*len(base)
        for j in a:
            row[j] += 1
        for j in b:
            row[j] -= 1
        equations.append(row)
    assert rank(equations, len(base)) == case["collision_rank"]
    cells = Counter()
    for _, target in first.values():
        if target is not None:
            pre = neg(mul(inv2, target, p), p)
            assert neg(add(pre, pre, p), p) == target
            cells[cell(pre, coefficients, p)] += 1
    assert dict(sorted(cells.items())) == case["character_cells"]
    support_by_case[name] = first
    print("PIVOT_COVERAGE_REPLAY", name, len(first), dict(cells), flush=True)

for name, law in result["known_scalar_query_laws"].items():
    scalars = guide["uniform_scalars" if name == "uniform" else "guided_scalars"]
    assert law["queries"] == len(scalars) == 64
    assert law["scalar_sha256"] == hashlib.sha256(
        "\n".join(map(str, scalars)).encode()).hexdigest()
    observed = Counter()
    hits = {case_name: [] for case_name in result["cases"]}
    for index, encoded in enumerate(scalars):
        scalar = int(encoded)
        assert 0 < scalar < r
        pre = mul(scalar, generator, p)
        target = neg(add(pre, pre, p), p)
        assert target == mul(scalar, neg(add(generator, generator, p), p), p)
        observed[cell(pre, coefficients, p)] += 1
        for case_name, support in support_by_case.items():
            if key(target) in support:
                hits[case_name].append((index, list(support[key(target)][0])))
    assert dict(sorted(observed.items())) == law["cells"]
    for case_name, matches in hits.items():
        claimed = law["hits"].get(case_name, [])
        assert [(item["query_index"], item["witness"]) for item in claimed] == matches
    print("PIVOT_QUERY_REPLAY", name, dict(observed),
          {case_name: len(matches) for case_name, matches in hits.items()}, flush=True)

output = {
    "status": "independent_two_pivot_coverage_replay_pass",
    "producer_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "cases_replayed": list(result["cases"]),
    "query_laws_replayed": list(result["known_scalar_query_laws"]),
    "boundary": "Integer-curve replay of exact support, character cells, collision rank, and known-scalar query hits for the frozen bases.",
}
(HERE / "type_i_two_pivot_guided_coverage_replay.json").write_text(
    json.dumps(output, indent=2, sort_keys=True) + "\n")
