#!/usr/bin/env python3
"""Exact target-blind third-factor probe-order screen on the archived m3 bases."""
from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
from itertools import permutations
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
M3 = ROOT / "research/notes/ecc2k130/m3_four_policy_20260930"
sys.path.insert(0, str(M3))
import verify as m3_verify  # noqa: E402

EVIDENCE = M3 / "evidence_run_36722040881"
POLICIES = ("original", "transported", "descendant_native", "pullback")
SCORE_GROUP_ADDITIONS = 8 * 36  # All 36 unordered pairs plus each of 8 thirds.


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_order_bound(curve, base: list[tuple[int, int]], targets: set) -> dict:
    assert len(base) == 8 and len(set(base)) == 8
    pairs = [curve.add(base[i], base[j]) for i in range(8) for j in range(i, 8)]
    assert len(pairs) == 36
    supports = [{curve.add(pair, base[k]) for pair in pairs} for k in range(8)]
    unions = [set() for _ in range(1 << 8)]
    for mask in range(1, 1 << 8):
        bit = mask & -mask
        k = bit.bit_length() - 1
        unions[mask] = unions[mask ^ bit] | supports[k]
    n_targets = len(targets)
    assert n_targets in (210, 420)
    unseen = [n_targets - len(union & targets) for union in unions]

    def order_cost(order) -> int:
        mask = 0
        cost = 0
        for k in order:
            cost += unseen[mask]
            mask |= 1 << k
        assert mask == 255
        return cost

    @lru_cache(None)
    def optimum(mask: int) -> tuple[int, tuple[int, ...]]:
        if mask == 255:
            return 0, ()
        choices = []
        for k in range(8):
            if not mask & (1 << k):
                tail_cost, tail_order = optimum(mask | (1 << k))
                choices.append((unseen[mask] + tail_cost, (k,) + tail_order))
        return min(choices)

    best_total, best_order = optimum(0)
    baseline_total = order_cost(range(8))
    # This second search checks the dynamic program against all 8! orders.
    brute_total, brute_order = min((order_cost(order), order)
                                   for order in permutations(range(8)))
    assert (best_total, best_order) == (brute_total, brute_order)
    assert best_total <= baseline_total
    saved_total = baseline_total - best_total
    break_even = (math.ceil(SCORE_GROUP_ADDITIONS * n_targets / saved_total)
                  if saved_total else None)
    return {
        "target_universe_size": n_targets,
        "distinct_three_sum_support": len(unions[255] & targets),
        "baseline_total_probes": baseline_total,
        "optimal_total_probes": best_total,
        "baseline_mean_probes": baseline_total / n_targets,
        "optimal_mean_probes": best_total / n_targets,
        "saved_probes_per_target": saved_total / n_targets,
        "optimal_order": list(best_order),
        "score_group_additions": SCORE_GROUP_ADDITIONS,
        "optimistic_break_even_targets_for_score": break_even,
        "saved_group_additions_at_L1024": 1024 * saved_total / n_targets,
    }


def derive() -> dict:
    m3_verify.checked_config()
    lock = json.loads((M3 / "FROZEN_REPLAY.json").read_text())
    assert sha(EVIDENCE / "result.json") == lock["result_sha256"]
    replay = json.loads((EVIDENCE / "replay_repaired.json").read_text())
    assert replay["status"] == "PASS" and replay["cells_replayed"] == 16
    assert replay["cases_replayed"] == 8192
    result = json.loads((EVIDENCE / "result.json").read_text())
    assert result["status"] == "PASS_PANEL" and result["subgroup_order"] == 421

    m3_verify.restore_bare_curve()
    field = m3_verify.pilot.FastGF2m(21, m3_verify.pilot.IRR)
    source = m3_verify.pilot.Koblitz(field, 0, 1)
    twist = m3_verify.pilot.Koblitz(field, 1, 1)
    lines, selected, _ = m3_verify.pilot.order_seven_lines(field, twist)
    forward = m3_verify.pilot.BinaryVeluMap.from_generator(
        source, twist, lines[selected]["generator"], 7)
    leaf = forward.codomain
    generator = tuple(result["challenge"]["G"])
    reps, by_point, _ = m3_verify.independent_orbits(source, generator)
    assert len(by_point) == 420 and len(reps) == 10
    source_targets = {
        "all": set(by_point),
        "A": {point for point, orbit in by_point.items() if orbit in set(reps[:5])},
        "B": {point for point, orbit in by_point.items() if orbit in set(reps[5:])},
    }
    assert len(source_targets["A"]) == len(source_targets["B"]) == 210
    leaf_targets = {label: {forward(point) for point in points}
                    for label, points in source_targets.items()}
    assert all(len(leaf_targets[label]) == len(source_targets[label])
               for label in source_targets)

    rows = []
    for seed in sorted(result["variants"]):
        for policy in POLICIES:
            is_leaf = policy in ("transported", "descendant_native")
            curve = leaf if is_leaf else source
            targets_by_label = leaf_targets if is_leaf else source_targets
            base = [tuple(point) for point in result["bases"][seed][policy]]
            for label in ("all", "A", "B"):
                row = {"seed": int(seed), "policy": policy, "universe": label,
                       **exact_order_bound(curve, base, targets_by_label[label])}
                rows.append(row)
    keyed = {(row["seed"], row["policy"], row["universe"]): row for row in rows}
    covariance_fields = ("target_universe_size", "distinct_three_sum_support",
                         "baseline_total_probes", "optimal_total_probes")
    for seed in sorted(result["variants"]):
        for left, right in (("original", "transported"),
                            ("pullback", "descendant_native")):
            for label in ("all", "A", "B"):
                first = keyed[int(seed), left, label]
                second = keyed[int(seed), right, label]
                assert all(first[field] == second[field] for field in covariance_fields)

    original_cold = []
    for seed in sorted(result["variants"]):
        for label in ("A", "B"):
            original_cold.append({
                "seed": int(seed), "holdout": label,
                "observed_source_base_field_mul": result["phase_costs"][
                    f"source_base_{seed}"]["mul"],
                "observed_original_scan_to_rank_field_mul": result["variants"][
                    seed][label]["original"]["costs"]["scan_to_rank"]["mul"],
            })
    return {"schema": "ecc2k130-m3-probe-order-bound-v1",
            "classification": "retrospective_probe_only_screen",
            "m3_result_sha256": lock["result_sha256"],
            "m3_repaired_replay_sha256": sha(EVIDENCE / "replay_repaired.json"),
            "pair_table_entries": 36, "third_factors": 8,
            "score_group_additions": SCORE_GROUP_ADDITIONS,
            "rows": rows, "original_cold_observations": original_cold,
            "method_speedup": None, "n131_transfer": None}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite a derived receipt"
    result = derive()
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"schema": result["schema"], "rows": len(result["rows"]),
                      "source": result["m3_result_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    main()
