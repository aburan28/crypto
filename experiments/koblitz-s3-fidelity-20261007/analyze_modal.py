"""Audit the exact-output Modal S3 panel and summarize exploratory paired costs.

Usage: python3 analyze_modal.py modal/pilot_raw_v2.json.gz modal/analysis_raw_v2.json
The twelve points per curve are the units of analysis. Repeating a panel on
those same points is a custody/noise replication, not 24 new targets.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys

from scipy import stats


PHASES = (
    "target_query", "target_pdp_charged", "target_relation_check",
    "target_descent", "target_recovery_check",
)


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def t_interval(values: list[float], confidence: float = 0.95) -> dict:
    mean = statistics.mean(values)
    sd = statistics.stdev(values)
    half = stats.t.ppf((1 + confidence) / 2, len(values) - 1) * sd / math.sqrt(len(values))
    return {"mean_log": mean, "sample_sd_log": sd,
            "geometric_mean_ratio": math.exp(mean),
            "descriptive_ratio_interval": [math.exp(mean - half), math.exp(mean + half)],
            "confidence": confidence}


def interaction(a: list[float], b: list[float]) -> dict:
    """Return N53 minus N41 on paired target-level log speedups."""
    va = statistics.variance(a) / len(a)
    vb = statistics.variance(b) / len(b)
    se = math.sqrt(va + vb)
    df = (va + vb) ** 2 / (va * va / (len(a) - 1) + vb * vb / (len(b) - 1))
    diff = statistics.mean(b) - statistics.mean(a)
    half = stats.t.ppf(0.975, df) * se
    return {"contrast": "N53 minus N41", "mean_log_difference": diff,
            "descriptive_ratio_of_ratios": math.exp(diff),
            "descriptive_95pct_interval": [math.exp(diff - half), math.exp(diff + half)],
            "welch_df": df, "descriptive_two_sided_p": 2 * stats.t.sf(abs(diff / se), df)}


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: python3 analyze_modal.py INPUT.json.gz OUTPUT.json")
    source, output = map(Path, sys.argv[1:])
    if output.exists():
        raise FileExistsError(output)
    compressed = source.read_bytes()
    panel = json.loads(gzip.decompress(compressed))
    runs = panel["runs"]
    blocks = panel["blocks"]
    by_tag = {run["tag"]: run for run in runs}
    if len(by_tag) != len(runs):
        raise ValueError("duplicate run tag")
    audited = []
    for run in runs:
        raw = run.get("raw_text")
        report = run.get("report")
        raw_matches = bool(raw is not None and sha(raw.encode()) == run.get("raw_sha256")
                           and json.loads(raw) == report)
        binary_stable = run.get("binary_sha256_before") == run.get("binary_sha256_after")
        timing = (report or {}).get("timing_ms", {})
        phases_complete = all(isinstance(timing.get(k), (int, float)) for k in PHASES)
        phases_sum = phases_complete and abs(
            sum(timing[k] for k in PHASES) - timing["target_online_after_reusable_setup"]
        ) < 0.02
        audited.append({"tag": run["tag"], "n": run["n"], "arm": run["arm"],
                        "exit_code": run["exit_code"],
                        "raw_hash_matches": raw_matches,
                        "stdout_matches_raw": run.get("stdout_matches_raw") is True,
                        "binary_stable": binary_stable,
                        "native_and_fixture_verified": run["verified_native_and_fixture"],
                        "exclusive_phases_sum": bool(phases_sum)})
    curve_results = {}
    for n in (41, 53):
        targets = []
        log_ratios = []
        aa_log_noise = []
        for block in (block for block in blocks if block["n"] == n):
            pair = [by_tag[tag] for tag in block["run_tags"][:4]]
            baseline = [r["report"]["timing_ms"]["target_online_after_reusable_setup"]
                        for r in pair if r["arm"] == "baseline"]
            candidate = [r["report"]["timing_ms"]["target_online_after_reusable_setup"]
                         for r in pair if r["arm"] == "candidate"]
            if len(baseline) != 2 or len(candidate) != 2:
                raise ValueError(f"incomplete pair: {block['label']}")
            d = math.log(math.sqrt(baseline[0] * baseline[1] /
                                   (candidate[0] * candidate[1])))
            aa = [by_tag[tag]["report"]["timing_ms"]["target_online_after_reusable_setup"]
                  for tag in block["run_tags"][4:]]
            if len(aa) != 2:
                raise ValueError(f"missing A/A control: {block['label']}")
            aa_d = abs(math.log(aa[0] / aa[1]))
            log_ratios.append(d)
            aa_log_noise.append(aa_d)
            targets.append({"label": block["label"], "public_point": block["point"],
                            "run_tags": block["run_tags"], "order": block["order"],
                            "baseline_online_ms": baseline,
                            "candidate_online_ms": candidate,
                            "aa_online_ms": aa, "exploratory_paired_ratio": math.exp(d),
                            "aa_abs_log_ratio": aa_d,
                            "all_native_verified": block["all_native_verified"],
                            "semantic_equal": block["semantic_equal"]})
        p95 = sorted(aa_log_noise)[math.ceil(0.95 * len(aa_log_noise)) - 1]
        curve_results[str(n)] = {"independent_target_count": len(targets),
                                 "faster_target_count": sum(x > 0 for x in log_ratios),
                                 "paired_log_ratio": log_ratios,
                                 "descriptive_t_interval": t_interval(log_ratios),
                                 "aa_p95_abs_log_ratio": p95,
                                 "aa_p95_ratio": math.exp(p95),
                                 "aa_noise_gate_pass": p95 < math.log(1.05),
                                 "factor_base_points": sorted({r["report"]["factor_base_points"]
                                                               for r in runs if r["n"] == n}),
                                 "folded_columns": sorted({r["report"]["orbit_columns"]
                                                            for r in runs if r["n"] == n}),
                                 "factor_base_digests": sorted({r["report"]["factor_base_digest"]
                                                                for r in runs if r["n"] == n}),
                                 "targets": targets}
    selection = panel["host"]["selection"]
    cgroup = panel["host"]["cgroup"]
    fidelity = {
        "strict_host_isolation_receipt_present": panel["strict_host_isolation_receipt"] is not None,
        "physical_smt_topology_known": selection["topology_known"],
        "numa_node_known": selection["node"] is not None,
        "exclusive_cpuset_visible": bool(cgroup.get("cpuset.cpus.exclusive.effective")),
        "cpuset_partition_visible": bool(cgroup.get("cpuset.cpus.partition")),
        "cpu_psi_visible": bool(panel["host"]["psi"].get("cpu")),
        "memory_psi_visible": bool(panel["host"]["psi"].get("memory")),
        "worker_thread_affinities_visible": all(bool(r["thread_masks_seen"]) for r in runs),
        "perf_whole_process_available": all(p["result"].get("exit_code") == 0
                                            for p in panel["perf_whole_process"]),
        "aa_noise_gate_pass_both_curves": all(curve_results[str(n)]["aa_noise_gate_pass"]
                                              for n in (41, 53)),
    }
    all_custody = all(all(value for key, value in row.items() if key in (
        "raw_hash_matches", "stdout_matches_raw", "binary_stable",
        "native_and_fixture_verified", "exclusive_phases_sum")) for row in audited)
    result = {
        "kind": "modal_s3_two_curve_pilot_audit_v2",
        "source_artifact": str(source), "source_artifact_sha256": sha(compressed),
        "source_commit": panel["source_commit"], "runner_sha256": panel["runner_sha256"],
        "modal_image_id": panel["host"].get("modal_image_id"),
        "binary_sha256": panel["host"]["binary_sha256"],
        "raw_run_count": len(runs), "target_block_count": len(blocks),
        "run_audit": audited, "all_run_custody_and_native_checks_pass": all_custody,
        "all_semantic_equal": all(b["semantic_equal"] for b in blocks),
        "fidelity": fidelity, "curves": curve_results,
        "interaction": interaction(curve_results["41"]["paired_log_ratio"],
                                   curve_results["53"]["paired_log_ratio"]),
        "confirmatory_targets_per_curve": None,
        "controlled_speedup": None,
        "claim_status": "exploratory: no physical-host isolation receipt; A/A gate must pass on both curves",
        "statistical_intervals_status": "descriptive only; not claim-grade because predeclared fidelity gate failed",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(output), "raw_runs": len(runs), "targets": len(blocks),
                      "all_custody": all_custody, "fidelity": fidelity,
                      "curve_ratios": {n: curve_results[n]["descriptive_t_interval"]["geometric_mean_ratio"]
                                       for n in ("41", "53")}}, sort_keys=True))


if __name__ == "__main__":
    main()
