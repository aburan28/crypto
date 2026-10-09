"""Validate one-target S3/rho blocks and summarize exploratory paired ratios."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import random
import statistics


HERE = Path(__file__).resolve().parent
CROSSCURVE = HERE.parent / "koblitz-s3-crosscurve-20261009"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")
Z95_ONE_SIDED = 1.6448536269514722
Z80 = 0.8416212335729143


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def bootstrap_geomean_interval(log_ratios: list[float]) -> list[float] | None:
    if len(log_ratios) < 4:
        return None
    rng = random.Random(0x533352484F)
    values = sorted(math.exp(statistics.mean(rng.choices(log_ratios,
                                                        k=len(log_ratios))))
                    for _ in range(10000))
    return [values[249], values[9749]]


def wilson_interval(wins: int, total: int) -> list[float] | None:
    if not total:
        return None
    z = 1.959963984540054
    p = wins / total
    den = 1 + z * z / total
    center = (p + z * z / (2 * total)) / den
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / den
    return [center - half, center + half]


def analyze(ledgers: list[Path], output: Path) -> dict:
    freeze_bytes = (HERE / "FREEZE.json").read_bytes()
    freeze = json.loads(freeze_bytes)
    plan_bytes = (CROSSCURVE / "FREEZE.json").read_bytes()
    plan = json.loads(plan_bytes)
    assert sha(plan_bytes) == freeze["crosscurve_plan_sha256"]
    assert ledgers
    ledger_bytes = [path.read_bytes() for path in ledgers]
    rows = [json.loads(line) for data in ledger_bytes for line in data.splitlines()]
    assert rows and len(rows) % 8 == 0
    targets = {(t["n"], t["label"]): t for t in plan["targets"]}
    grouped = defaultdict(list)
    binary_hashes = defaultdict(set)
    run_ids = set()
    for row in rows:
        n, label, arm = row["n"], row["label"], row["arm"]
        assert n in (41, 53) and arm in ("s3", "rho")
        target = targets[(n, label)]
        assert row["kind"] == "s3_rho_one_target_run_v1"
        assert row["target_sha256"] == target["point_sha256"]
        assert row["workload_id"] == target["workload_id"]
        expected_run = (f"{plan['curves'][str(n)]['rank243_candidate_id']}"
                        if arm == "s3" else f"RHO1N{n}Csignedfrobenius")
        expected_run += f"W{target['workload_id']}R{row['run']}"
        assert row["run_id"] == expected_run and row["run_id"] not in run_ids
        run_ids.add(row["run_id"])
        assert row["run"] in (1, 2)
        assert row["candidate_id"] == (plan["curves"][str(n)]["rank243_candidate_id"]
                                         if arm == "s3" else None)
        expected_hash = (plan["partial_binary_sha256_local"] if arm == "s3" else
                         freeze["rho_binary_sha256_local"])
        assert row["binary_sha256"] == expected_hash
        binary_hashes[arm].add(expected_hash)
        if arm == "s3":
            report = row["report"]
            if report:
                timing = report.get("timing_ms", {})
                actual_phase_ok = (all(isinstance(timing.get(k), (int, float))
                                       and timing[k] >= 0 for k in PHASES)
                                   and abs(sum(timing[k] for k in PHASES) -
                                           timing["target_online_after_reusable_setup"]) < .02)
                assert row["phase_ok"] == actual_phase_ok
                if row["verified_native_and_fixture"]:
                    assert report["recovered_scalar"] == target["fixture_scalar"]
                    assert report["group_verified"] and report["target_count"] == 1
                    assert report["target"] == target["public_point"]
                    assert report["prepared_rank_reached"] is True
            else:
                assert not row["verified_native_and_fixture"]
        else:
            report = row["report"]
            fixture, summary = report.get("fixture"), report.get("summary")
            if fixture and summary:
                actual_phase_ok = (fixture["online_ns"] == fixture["walk_ns"] +
                                   fixture["collision_ns"] + fixture["recovery_check_ns"])
                assert row["phase_ok"] == actual_phase_ok
                if row["verified_native_and_fixture"]:
                    assert fixture["recovered_fixture_scalar"] == target["fixture_scalar"]
                    assert fixture["published_q"] == target["public_point"]
                    assert fixture["table_entries_before"] == 0
                    assert fixture["verified"] and summary["all_verified"]
                    assert fixture["online_stop_event"] == "recovery_check_true"
                    assert summary["fixtures"] == 1 and summary["target_source"] == "public_point"
            else:
                assert not row["verified_native_and_fixture"]
        assert not row["verified_native_and_fixture"] or row["exit_code"] == 0
        grouped[(n, label)].append(row)
    indices = {int(label[1:]) for _, label in grouped}
    assert set(grouped) == {(n, f"T{i:03}") for i in indices for n in (41, 53)}
    blocks = []
    for (n, label), block in sorted(grouped.items()):
        assert len(block) == 4 and Counter(row["arm"] for row in block) == {
            "s3": 2, "rho": 2}
        expected_order = (("rho", "s3", "s3", "rho")
                          if (int(label[1:]) + (n == 53)) % 2 else
                          ("s3", "rho", "rho", "s3"))
        assert tuple(row["arm"] for row in block) == expected_order
        good = all(row["verified_native_and_fixture"] and row["phase_ok"] for row in block)
        item = {"n": n, "label": label, "workload_id": targets[(n, label)]["workload_id"],
                "all_verified": good, "run_ids": [r["run_id"] for r in block],
                "failed_run_ids": [r["run_id"] for r in block
                                   if not r["verified_native_and_fixture"]],
                "exploratory_rho_over_s3_online_ratio": None,
                "online_ms_by_arm": None, "aa_abs_log_ratio_by_arm": None,
                "rho_walk_steps": None, "rho_dp_entries": None,
                "s3_online_root_calls": None, "s3_fallback_attempts": None}
        if good:
            by_arm = {arm: [r for r in block if r["arm"] == arm]
                      for arm in ("s3", "rho")}
            times = {arm: [r["online_ms"] for r in by_arm[arm]]
                     for arm in ("s3", "rho")}
            assert all(isinstance(v, (int, float)) and v > 0
                       for vv in times.values() for v in vv)
            item.update({
                "exploratory_rho_over_s3_online_ratio":
                    math.sqrt(math.prod(times["rho"]) / math.prod(times["s3"])),
                "online_ms_by_arm": times,
                "aa_abs_log_ratio_by_arm": {
                    arm: abs(math.log(times[arm][0] / times[arm][1]))
                    for arm in ("s3", "rho")},
                "rho_walk_steps": [r["report"]["fixture"]["walk_steps"]
                                   for r in by_arm["rho"]],
                "rho_dp_entries": [r["report"]["fixture"]["table_entries_after"]
                                   for r in by_arm["rho"]],
                "s3_online_root_calls": [r["report"]["fallback_s3_calls"] +
                                         r["report"]["target_s3_calls"]
                                         for r in by_arm["s3"]],
                "s3_fallback_attempts": [r["report"]["fallback_attempts"]
                                          for r in by_arm["s3"]],
            })
        blocks.append(item)
    curves = {}
    for n in (41, 53):
        selected = [b for b in blocks if b["n"] == n]
        valid = [b for b in selected if b["all_verified"]]
        ratios = [b["exploratory_rho_over_s3_online_ratio"] for b in valid]
        logs = [math.log(r) for r in ratios]
        wins = sum(r >= 2 for r in ratios)
        aa = {arm: sorted(b["aa_abs_log_ratio_by_arm"][arm] for b in valid)
              for arm in ("s3", "rho")}
        margins = [value - math.log(2) for value in logs]
        mean_margin = statistics.mean(margins) if margins else None
        needed = (max(24, math.ceil(((Z95_ONE_SIDED + Z80) *
                                     statistics.stdev(margins) / mean_margin) ** 2))
                  if len(margins) >= 4 and mean_margin > 0 else None)
        curves[str(n)] = {
            "target_count": len(selected), "verified_target_count": len(valid),
            "failed_run_count": sum(len(b["failed_run_ids"]) for b in selected),
            "at_least_2_online_ratio_count": wins,
            "at_least_2_fraction_wilson_95": wilson_interval(wins, len(valid)),
            "exploratory_rho_over_s3_geomean": math.exp(statistics.mean(logs)) if logs else None,
            "exploratory_rho_over_s3_geomean_bootstrap_95": bootstrap_geomean_interval(logs),
            "exploratory_rho_over_s3_min": min(ratios) if ratios else None,
            "aa_p95_abs_log_ratio": {
                arm: (values[max(0, math.ceil(.95 * len(values)) - 1)] if values else None)
                for arm, values in aa.items()},
            "s3_fallback_attempts_total": sum(sum(b["s3_fallback_attempts"]) for b in valid),
            "provisional_pairs_for_geomean_above_2_one_sided_95_80_percent_power": needed,
        }
    result = {
        "kind": "s3_rho_one_target_local_pilot_analysis_v1",
        "freeze_sha256": sha(freeze_bytes),
        "crosscurve_plan_sha256": sha(plan_bytes),
        "ledger_sha256": [sha(data) for data in ledger_bytes],
        "run_count": len(rows), "matched_index_count": len(indices),
        "binary_sha256": {arm: next(iter(values)) for arm, values in binary_hashes.items()},
        "rho_threads": freeze["rho_threads"], "s3_workers": freeze["s3_workers"],
        "host_isolation_receipt": None, "controlled_speedup": None,
        "sample_size_estimates_are_exploratory": True,
        "curves": curves, "blocks": blocks,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(output), "runs": len(rows), "curves": curves},
                     sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.ledger, args.output)
