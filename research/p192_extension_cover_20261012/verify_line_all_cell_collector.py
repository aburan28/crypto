"""Independent prime-field group replay of the 141-class three-sum lookup."""

import hashlib
import json
import random
from collections import Counter
from pathlib import Path

from sage.all import EllipticCurve, GF, Integer, matrix

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
path = HERE / "type_i_line_all_cell_collector.json"
v = json.loads(path.read_text())
assert v["source_sha256"] == sha(HERE / "line_all_cell_collector.py")
assert v["prior_receipt_sha256"] == sha(HERE / "type_i_line_guided_collector.json")
assert v["eligibility_receipt_sha256"] == sha(HERE / "type_i_line_fiber_eligibility.json")
assert v["reduction_receipt_sha256"] == sha(HERE / "type_i_line_fixed_origin_reduction.json")
assert v["runtime_sha256"] == sha(HERE / "type_i_line_all_cell_runtime_info.json")
prior = json.loads((HERE / "type_i_line_guided_collector.json").read_text())
p, r = Integer(v["field_prime"]), Integer(v["subgroup_order"])
assert p == Integer("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF", 16)
F = GF(p)
curve_a = p - 3
curve_b = Integer("64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1", 16)
gx = Integer("188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012", 16)
gy = Integer("07192B95FFC8DA78631011ED6B24CDD573F977A11E794811", 16)
E = EllipticCurve(F, [F(curve_a), F(curve_b)])
G = E(F(gx), F(gy))
O = E(0)
assert r * G == O and G != O
base = [E(F(x), F(y)) for x, y in prior["base_points_fp"]]
assert len(base) == 141 and len(set(base)) == 141
assert all(r * Q == O for Q in base)
pair_sums = {}
for i in range(141):
    for j in range(i, 141):
        pair_sums.setdefault(base[i] + base[j], []).append((i, j))
assert sum(map(len, pair_sums.values())) == 10011
rng = random.Random(v["seed"])
for trial in v["trial_rows"]:
    assert int(trial["scalar"]) == rng.randrange(1, int(r))
assert v["first16_distinct_three_sum_targets"] == 816
assert sum(v["first16_support_by_cell"].values()) + v["first16_exceptional_support"] == 816
assert v["first16_support_by_cell"]["++"] == 206
for arm, expected_cells in (("control", {"++": 64}),
                            ("balanced", {"++": 16, "+-": 16, "-+": 16, "--": 16})):
    queries = v[f"{arm}_queries"]
    assert len(queries) == 64
    assert len({q["trial_index"] for q in queries}) == 64
    for q in queries:
        assert q["scalar"] == v["trial_rows"][q["trial_index"]]["scalar"]
        assert q.get("cell", v["trial_rows"][q["trial_index"]]["cell"]) == v["trial_rows"][q["trial_index"]]["cell"]
    assert Counter(v["trial_rows"][q["trial_index"]]["cell"] for q in queries) == expected_cells
    found = []
    for qi, q in enumerate(queries):
        a = Integer(q["scalar"])
        T = -2 * a * G
        for k in range(141):
            for i, j in pair_sums.get(T - base[k], []):
                if j <= k:
                    assert base[i] + base[j] + base[k] == T
                    found.append((qi, q["trial_index"], str(a), i, j, k))
    reported = [(h["query_index"], h["trial_index"], h["scalar"], *h["indices"])
                for h in v[f"{arm}_hits"]]
    assert sorted(found) == sorted(reported)
    rows = []
    for _, _, _, i, j, k in found:
        row = [0] * 141
        for t in (i, j, k):
            row[t] += 1
        rows.append(row)
    rank = int(matrix(GF(r), rows).rank()) if rows else 0
    assert v[f"{arm}_rank"] == rank
print("LINE_ALL_CELL_INDEPENDENT_PASS", json.dumps({
    "classes": 141, "queries_per_arm": 64,
    "control_hits": len(v["control_hits"]),
    "balanced_hits": len(v["balanced_hits"])}, sort_keys=True))
