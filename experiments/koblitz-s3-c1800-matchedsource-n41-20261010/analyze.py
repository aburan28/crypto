"""Audit 24 N41 same-source C1200/C1800 pairs and their online ratio."""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import random
import statistics


HERE = Path(__file__).resolve().parent


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def geometric(values: list[float]) -> float:
    assert values and all(value > 0 for value in values)
    return math.exp(statistics.mean(math.log(value) for value in values))


def analyze(ledger: Path, output: Path) -> dict:
    freeze_path = HERE / "FREEZE.json"
    freeze = json.loads(freeze_path.read_bytes())
    raw = ledger.read_bytes() if ledger.exists() else b""
    rows = [json.loads(line) for line in raw.splitlines()]
    expected = [(entry["label"], entry["arm"], entry["repetition"])
                for entry in freeze["run_order"]]
    actual = [(row["label"], row["arm"], row["repetition"]) for row in rows]
    order_prefix_ok = actual == expected[:len(actual)]
    targets = {target["label"]: target for target in freeze["selected_targets"]}
    audited = []
    for row in rows:
        target = targets[row["label"]]
        frozen = freeze["arms"][row["arm"]]
        stem = f"{row['label']}_{row['arm']}_R{row['repetition']}"
        report_file = ledger.parent / f"{stem}.json"
        stdout_file = ledger.parent / f"{stem}.stdout"
        stderr_file = ledger.parent / f"{stem}.stderr"
        report_raw = report_file.read_bytes() if report_file.exists() else b""
        stdout_raw = stdout_file.read_bytes() if stdout_file.exists() else b""
        stderr_raw = stderr_file.read_bytes() if stderr_file.exists() else b""
        try:
            report = json.loads(report_raw) if report_raw else None
        except ValueError:
            report = None
        valid = bool(row["verified_native_and_fixture"] and report
                     and row["phase_ok"] and row["exit_code"] == 0
                     and row["error"] is None
                     and row["candidate_id"] == frozen["candidate_id"]
                     and row["run_id"] ==
                     f"{frozen['candidate_id']}W{target['workload_id']}R{row['repetition']}"
                     and row["binary_sha256"] == frozen["binary_sha256"]
                     and row["source_sha256"] == frozen["source_sha256"]
                     and row["workload_id"] == target["workload_id"]
                     and row["point_sha256"] == target["point_sha256"]
                     and row["report_sha256"] == sha(report_raw)
                     and row["stdout_sha256"] == sha(stdout_raw)
                     and row["stderr_sha256"] == sha(stderr_raw)
                     and stdout_raw == report_raw
                     and report["target"] == target["public_point"]
                     and report["recovered_scalar"] == target["fixture_scalar"]
                     and report["group_verified"] is True
                     and report["factor_base_points"] == frozen["factor_base_points"]
                     and report["factor_base_digest"] == frozen["factor_base_digest"]
                     and report["prepared_rank"] == frozen["prepare_rank"]
                     and report["root_table_slots"] == frozen["root_table_slots"]
                     and row["online_ms"] ==
                     report["timing_ms"]["target_online_after_reusable_setup"]
                     and row["target_pdp_charged_ms"] ==
                     report["timing_ms"]["target_pdp_charged"]
                     and row["target_s3_calls"] == report["target_s3_calls"]
                     and row["reusable_setup_ms"] ==
                     report["timing_ms"]["reusable_setup_total"]
                     and row["peak_rss_bytes"] <=
                     freeze["observed_peak_rss_envelope_bytes"])
        audited.append({"label": row["label"], "arm": row["arm"],
                        "valid": valid, "report_sha256": row["report_sha256"]})
    complete = bool(order_prefix_ok and len(rows) == len(expected)
                    and all(entry["valid"] for entry in audited))
    result = {
        "kind": "s3_n41_c1200_c1800_matchedsource24_analysis_v1",
        "freeze_sha256": sha(freeze_path.read_bytes()),
        "ledger_sha256": sha(raw), "planned_runs": len(expected),
        "observed_runs": len(rows),
        "verified_runs": sum(entry["valid"] for entry in audited),
        "order_prefix_ok": order_prefix_ok, "complete": complete,
        "audited_rows": audited, "target_results": [],
        "online_geometric_baseline_over_expanded": None,
        "online_target_bootstrap_95pct": None,
        "s3_calls_geometric_baseline_over_expanded": None,
        "development_gate_at_least_twofold": None,
        "matched_source_gate_bootstrap_lower_at_least_twofold": None,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    if complete:
        by_label = defaultdict(dict)
        for row in rows:
            by_label[row["label"]][row["arm"]] = row
        online_ratios = []
        call_ratios = []
        for label in freeze["selected_labels"]:
            pair = by_label[label]
            assert set(pair) == {"baseline", "expanded"}
            base, expanded = pair["baseline"], pair["expanded"]
            online = base["online_ms"] / expanded["online_ms"]
            calls = base["target_s3_calls"] / expanded["target_s3_calls"]
            online_ratios.append(online)
            call_ratios.append(calls)
            result["target_results"].append({
                "label": label, "workload_id": targets[label]["workload_id"],
                "baseline_online_ms": base["online_ms"],
                "expanded_online_ms": expanded["online_ms"],
                "baseline_target_pdp_charged_ms": base["target_pdp_charged_ms"],
                "expanded_target_pdp_charged_ms": expanded["target_pdp_charged_ms"],
                "baseline_target_s3_calls": base["target_s3_calls"],
                "expanded_target_s3_calls": expanded["target_s3_calls"],
                "online_ratio": online, "s3_calls_ratio": calls,
            })
        rng = random.Random(20261010)
        boot = sorted(geometric([online_ratios[rng.randrange(len(online_ratios))]
                                 for _ in online_ratios]) for _ in range(20000))
        result["online_geometric_baseline_over_expanded"] = geometric(online_ratios)
        result["online_target_bootstrap_95pct"] = [boot[499], boot[19499]]
        result["s3_calls_geometric_baseline_over_expanded"] = geometric(call_ratios)
        result["development_gate_at_least_twofold"] = (
            result["online_geometric_baseline_over_expanded"] >= 2.0)
        result["matched_source_gate_bootstrap_lower_at_least_twofold"] = (
            result["online_target_bootstrap_95pct"][0] >= 2.0)
        result["median_online_ms"] = {
            arm: statistics.median(row["online_ms"] for row in rows if row["arm"] == arm)
            for arm in ("baseline", "expanded")}
        result["median_reusable_setup_ms"] = {
            arm: statistics.median(row["reusable_setup_ms"] for row in rows
                                   if row["arm"] == arm)
            for arm in ("baseline", "expanded")}
        result["max_peak_rss_bytes"] = max(row["peak_rss_bytes"] for row in rows)
        result["total_major_page_faults"] = sum(
            row["resource_counters"]["major_page_faults"] for row in rows)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"complete": complete, "verified_runs": result["verified_runs"],
                      "online_ratio": result["online_geometric_baseline_over_expanded"],
                      "call_ratio": result["s3_calls_geometric_baseline_over_expanded"],
                      "development_gate": result["development_gate_at_least_twofold"]},
                     sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyze(args.ledger.resolve(), args.output.resolve())
