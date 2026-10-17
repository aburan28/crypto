"""Exact 205-class three-sum support and known-scalar geometry-guided laws."""

import collections
import hashlib
import json
import math
import sys
import time
from pathlib import Path

from cypari2 import Pari

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from line_full_character_support import rank_mod  # noqa: E402
from verify_adaptive_private_columns import add, cell, mul, neg, on_curve  # noqa: E402


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


paths = {
    "old": HERE / "type_i_line_guided_collector.json",
    "new": HERE / "type_i_line_new64_fixed_origin.json",
    "density": HERE / "type_i_character_density.json",
    "collector": HERE / "type_i_line_all_cell_collector.json",
}
data = {name: json.loads(path.read_text()) for name, path in paths.items()}
p = int(data["old"]["field_prime"])
r = int(data["old"]["subgroup_order"])
assert data["new"]["reduced_classes"] == 64
old = [tuple(map(int, item)) for item in data["old"]["base_points_fp"]]
new = [tuple(map(int, row["quotient"])) for row in data["new"]["rows"]]
base = old + new
n = len(base)
assert len(old) == 141 and len(new) == 64 and n == 205
assert len(set(base)) == n and all(on_curve(point, p) for point in base)
coefficients = [int(v) for v in data["density"]["sextic_coefficients_decimal"]]
pari = Pari()
assert int(pari.isprime(p)) == 1
pari_p = pari(p)
inv2 = pow(2, -1, r)
halves = [neg(mul(inv2, point, p), p) for point in base]


def chi(point):
    if point is None:
        return "I"
    value = 0
    for c in reversed(coefficients):
        value = (value * point[0] + c) % p
    if value == 0:
        return "0"
    return "+" if int(pari.kronecker(value, pari_p)) == 1 else "-"


def character(half):
    return chi(add(half, half, p)) + chi(half)


for i, row in enumerate(data["new"]["rows"]):
    half = neg(halves[141 + i], p)
    assert row["cell"] == character(half) == cell(half, coefficients, p)[::-1]

start_ns = time.monotonic_ns()
support = {}
repeated = collections.defaultdict(list)
counts = collections.Counter()
first16 = collections.Counter()
euler_controls = 0
collision_rows = []
triple_count = 0
for i in range(n):
    for j in range(i, n):
        pair = add(base[i], base[j], p)
        pair_half = add(halves[i], halves[j], p)
        for k in range(j, n):
            target = add(pair, base[k], p)
            half = add(pair_half, halves[k], p)
            category = character(half)
            triple_count += 1
            if triple_count % 2003 == 0:
                assert category == cell(half, coefficients, p)[::-1]
                euler_controls += 1
            triple = (i, j, k)
            if target in support:
                repeated[target].append(triple)
                row = collections.Counter(triple)
                row.subtract(support[target])
                collision_rows.append(dict(row))
            else:
                support[target] = triple
                counts[category] += 1
            if k < 16:
                first16[category] += 1
    if (i + 1) % 16 == 0:
        print("NEW64_CHARACTER_PROGRESS", i + 1, n, len(support), flush=True)

assert triple_count == math.comb(n + 2, 3) == 1456935
assert len(support) + len(collision_rows) == triple_count
assert first16 == collections.Counter({"++": 206, "+-": 218,
                                       "-+": 199, "--": 193})
assert set(counts) == {"++", "+-", "-+", "--"}

generator = (int("188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012", 16),
             int("07192B95FFC8DA78631011ED6B24CDD573F977A11E794811", 16))
arms = {}
for name in ("control", "balanced"):
    queries = data["collector"][name + "_queries"]
    assert len(queries) == 64
    hits = []
    for query in queries:
        a = int(query["scalar"])
        target = neg(mul(2 * a, generator, p), p)
        half = mul(a, generator, p)
        assert add(half, half, p) == neg(target, p)
        assert character(half) == data["collector"]["trial_rows"][query["trial_index"]]["cell"]
        if target in support:
            for triple in (support[target], *repeated.get(target, [])):
                hits.append({"scalar": str(a), "indices": triple})
    arms[name] = {"queries": len(queries), "hits": hits,
                  "rank": rank_mod([hit["indices"] for hit in hits], r)}

out = {
    "status": "verified_205_class_character_three_sum_collector",
    "field_prime": str(p), "subgroup_order": str(r),
    "base_classes": n, "old_classes": len(old), "new_classes": len(new),
    "unordered_triples": triple_count,
    "distinct_three_sum_targets": len(support),
    "duplicate_decompositions": len(collision_rows),
    "collision_relation_rank": rank_mod(collision_rows, r),
    "support_by_Q_H_cell": dict(sorted(counts.items())),
    "first16_by_Q_H_cell": dict(sorted(first16.items())),
    "plus_plus_only_coverage": counts["++"],
    "four_cell_coverage": len(support),
    "euler_character_controls": euler_controls,
    "known_scalar_arms": arms,
    "elapsed_ns_unisolated": time.monotonic_ns() - start_ns,
    "source_sha256": sha(Path(__file__)),
    "input_sha256": {name: sha(path) for name, path in paths.items()},
    "boundary": "Exact support and two frozen known-scalar query laws on 205 reduced quotient classes; collision rank concerns only equal three-sum targets, and one-target pipeline cost is separate.",
}
path = HERE / "type_i_line_new64_character_collector.json"
path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print("NEW64_CHARACTER_COLLECTOR_PASS", json.dumps({
    "base_classes": n, "support": len(support),
    "duplicate_decompositions": len(collision_rows),
    "collision_rank": out["collision_relation_rank"],
    "support_by_cell": out["support_by_Q_H_cell"],
    "control_hits": len(arms["control"]["hits"]),
    "balanced_hits": len(arms["balanced"]["hits"]),
}, sort_keys=True), flush=True)
