"""Known-scalar four-cell query law on the frozen P-192 line quotient base."""

import hashlib
import json
import random
import time
from collections import Counter
from pathlib import Path

from sage.all import EllipticCurve, GF, Integer, matrix, prod
import type_i_divisor_transfer as dt

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
prior_path = HERE / "type_i_line_guided_collector.json"
eligibility_path = HERE / "type_i_line_fiber_eligibility.json"
reduction_path = HERE / "type_i_line_fixed_origin_reduction.json"
prior = json.loads(prior_path.read_text())
eligibility = json.loads(eligibility_path.read_text())
reduction = json.loads(reduction_path.read_text())
assert prior["base_classes"] == reduction["reduced_classes"] == 141
assert eligibility["cells"] == {"++": 30, "+-": 26, "-+": 40, "--": 45, "exceptional": 0}

p, r, L = dt.p, dt.r, dt.L
E = EllipticCurve(L, [L(dt.gm.a), L(dt.gm.b)])
G = E(L(dt.gm.gx), L(dt.gm.gy))
O = E(0)
inv2 = Integer(pow(2, -1, int(r)))
assert r * G == O and r.is_prime()

def key(P):
    if P == O:
        return "identity"
    assert P[0] ** p == P[0] and P[1] ** p == P[1]
    return f"{int(P[0].polynomial()[0])}:{int(P[1].polynomial()[0])}"

def fiber_character(P):
    if P == O or P[0] == dt.c0:
        return "exceptional"
    X = dt.inf + dt.d0 / (P[0] - dt.c0)
    rhs = prod(X - root ** p for root in dt.quartic_roots)
    assert rhs ** (p ** 3) == rhs
    if rhs == 0:
        return "exceptional"
    return "+" if rhs ** ((p ** 3 - 1) // 2) == 1 else "-"

def cell_for_target(T):
    # For a query target T=-2H, the geometric fiber pair is Q=2H,H.
    H = -inv2 * T
    Q = 2 * H
    first, second = fiber_character(Q), fiber_character(H)
    return first + second if first != "exceptional" and second != "exceptional" else "exceptional"

start = time.monotonic_ns()
base = [E(L(Integer(x)), L(Integer(y))) for x, y in prior["base_points_fp"]]
assert len(base) == 141 and len({key(P) for P in base}) == 141
assert all(r * P == O for P in base)
for index, P in enumerate(base):
    assert cell_for_target(-P) == eligibility["rows"][index]["cell"]

pairs = {}
for i in range(141):
    for j in range(i, 141):
        pairs.setdefault(key(base[i] + base[j]), []).append((i, j))
assert sum(map(len, pairs.values())) == 10011

support = {}
for i in range(16):
    for j in range(i, 16):
        for k in range(j, 16):
            T = base[i] + base[j] + base[k]
            support.setdefault(key(T), T)
assert len(support) == prior["first16_distinct_three_sum_targets"]
support_cells = Counter(cell_for_target(T) for T in support.values())
assert support_cells["++"] == prior["first16_character_admitted_three_sum_targets"]
assert sum(support_cells.values()) == len(support)

# One frozen scalar stream feeds both arms. Control: first 64 ++.
# Broadened policy: first 16 from each of the four cells.
seed = 0x19220261015
rng = random.Random(seed)
trial_rows = []
by_cell = {name: [] for name in ("++", "+-", "-+", "--")}
until_control = None
until_balanced = None
while until_control is None or until_balanced is None:
    a = Integer(rng.randrange(1, int(r)))
    T = -2 * a * G
    assert r * T == O
    cell = cell_for_target(T)
    trial_rows.append({"scalar": str(a), "cell": cell})
    if cell in by_cell:
        by_cell[cell].append((a, T, len(trial_rows) - 1))
    if until_control is None and len(by_cell["++"]) >= 64:
        until_control = len(trial_rows)
    if until_balanced is None and all(len(by_cell[c]) >= 16 for c in by_cell):
        until_balanced = len(trial_rows)
    assert len(trial_rows) < 4096

def collect(items):
    hits = []
    for query_index, (a, T, trial_index) in enumerate(items):
        assert T == -2 * a * G
        for k in range(141):
            for i, j in pairs.get(key(T - base[k]), []):
                if j > k:
                    continue
                assert base[i] + base[j] + base[k] == T
                hits.append({"query_index": query_index, "trial_index": trial_index,
                             "scalar": str(a), "indices": [i, j, k]})
    rows = []
    for hit in hits:
        row = [0] * 141
        for i in hit["indices"]:
            row[i] += 1
        rows.append(row)
    return hits, int(matrix(GF(r), rows).rank()) if rows else 0

control = by_cell["++"][:64]
balanced = sorted((item for cell in by_cell for item in by_cell[cell][:16]),
                  key=lambda item: item[2])
assert len(control) == len(balanced) == 64
assert Counter(cell_for_target(item[1]) for item in balanced) == dict.fromkeys(by_cell, 16)
control_hits, control_rank = collect(control)
balanced_hits, balanced_rank = collect(balanced)
out = {
    "status": "verified_line_all_cell_known_scalar_query_law_comparison",
    "candidate_id": None,
    "field_prime": str(p), "subgroup_order": str(r),
    "base_classes": 141, "pair_index_entries": 10011,
    "first16_distinct_three_sum_targets": len(support),
    "first16_support_by_cell": {c: support_cells[c] for c in by_cell},
    "first16_exceptional_support": support_cells["exceptional"],
    "seed": seed, "trial_count": len(trial_rows),
    "control_ready_after_trials": until_control,
    "balanced_ready_after_trials": until_balanced,
    "trial_rows": trial_rows,
    "control_queries": [{"scalar": str(a), "trial_index": i} for a, _, i in control],
    "balanced_queries": [{"scalar": str(a), "trial_index": i,
                           "cell": cell_for_target(T)} for a, T, i in balanced],
    "control_hits": control_hits, "balanced_hits": balanced_hits,
    "control_rank": control_rank, "balanced_rank": balanced_rank,
    "elapsed_ns_unisolated": time.monotonic_ns() - start,
    "source_sha256": sha(__file__),
    "prior_receipt_sha256": sha(prior_path),
    "eligibility_receipt_sha256": sha(eligibility_path),
    "reduction_receipt_sha256": sha(reduction_path),
    "runtime_sha256": sha(HERE / "type_i_line_all_cell_runtime_info.json"),
    "boundary": "Exact first16 target-support coverage and 64-query exact three-sum lookup on 141 P-192 quotient classes. Both laws select known scalars from one deterministic stream; local time is diagnostic.",
}
(HERE / "type_i_line_all_cell_collector.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("LINE_ALL_CELL_COLLECTOR_PASS", json.dumps({
    "support_cells": out["first16_support_by_cell"],
    "control_hits": len(control_hits), "balanced_hits": len(balanced_hits),
    "control_rank": control_rank, "balanced_rank": balanced_rank,
    "trials": len(trial_rows)}, sort_keys=True), flush=True)
