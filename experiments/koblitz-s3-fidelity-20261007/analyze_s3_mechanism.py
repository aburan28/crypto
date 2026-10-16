"""Attribute the archived Modal VM S3 timing diagnostic to exclusive phases.

This reads the frozen paired ledger and its custody audit. It produces a
descriptive stage diagnostic, not a controlled speedup or a new target sample.

Usage: python3 analyze_s3_mechanism.py RAW.json.gz AUDIT.json OUTPUT.json
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
ONLINE = "target_online_after_reusable_setup"
NESTED = ("rank_pdp_wall", "rank_query_generation", "rank_relation_check",
          "rank_linear_algebra_incremental", "index_build")
COUNTERS = ("rank_state_probes", "rank_target_s3_calls", "target_state_probes",
            "target_s3_calls")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def paired_mean_interval(values: list[float]) -> dict:
    if len(values) < 2:
        raise ValueError("at least two independent target blocks required")
    mean = statistics.mean(values)
    sd = statistics.stdev(values)
    half = stats.t.ppf(0.975, len(values) - 1) * sd / math.sqrt(len(values))
    return {"mean_baseline_minus_candidate_ms": mean,
            "sample_sd_paired_delta_ms": sd,
            "descriptive_95pct_interval_ms": [mean - half, mean + half],
            "target_count": len(values)}


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit("usage: python3 analyze_s3_mechanism.py RAW.json.gz AUDIT.json OUTPUT.json")
    raw_path, audit_path, output = map(Path, sys.argv[1:])
    if output.exists():
        raise FileExistsError(output)
    compressed = raw_path.read_bytes()
    panel = json.loads(gzip.decompress(compressed))
    audit_bytes = audit_path.read_bytes()
    audit = json.loads(audit_bytes)
    if audit["source_artifact_sha256"] != sha(compressed):
        raise ValueError("custody audit is for a different raw panel")
    if not audit["all_run_custody_and_native_checks_pass"] or not audit["all_semantic_equal"]:
        raise ValueError("raw panel did not pass custody and semantic checks")
    if len(panel["runs"]) != 146 or len(panel["blocks"]) != 24:
        raise ValueError("incomplete frozen VM panel")

    source_dir = Path(__file__).resolve().parents[1] / "koblitz-s3-pair-query-20261007-v4"
    for arm in ("baseline", "candidate"):
        if sha((source_dir / f"{arm}.rs").read_bytes()) != panel["host"]["source_sha256"][arm]:
            raise ValueError(f"{arm} source no longer matches VM panel")
    by_tag = {run["tag"]: run for run in panel["runs"]}
    if len(by_tag) != len(panel["runs"]):
        raise ValueError("duplicate run tag")

    curves = {}
    for n in (41, 53):
        target_rows = []
        for block in (b for b in panel["blocks"] if b["n"] == n):
            if not block["all_native_verified"] or not block["semantic_equal"]:
                raise ValueError(f"unverified target block {n}/{block['label']}")
            pair = [by_tag[tag] for tag in block["run_tags"][:4]]
            arms = {arm: [run["report"] for run in pair if run["arm"] == arm]
                    for arm in ("baseline", "candidate")}
            if any(len(reports) != 2 for reports in arms.values()):
                raise ValueError(f"incomplete ABBA/BAAB pair {n}/{block['label']}")
            row = {"label": block["label"], "run_tags": block["run_tags"][:4],
                   "timing_ms": {}, "counters": {}}
            for key in (ONLINE,) + PHASES + NESTED:
                means = {arm: statistics.mean(report["timing_ms"][key]
                                               for report in reports)
                         for arm, reports in arms.items()}
                row["timing_ms"][key] = {
                    "baseline": means["baseline"], "candidate": means["candidate"],
                    "baseline_minus_candidate": means["baseline"] - means["candidate"]}
            for key in COUNTERS:
                means = {arm: statistics.mean(report[key] for report in reports)
                         for arm, reports in arms.items()}
                row["counters"][key] = means
            online_delta = row["timing_ms"][ONLINE]["baseline_minus_candidate"]
            phase_delta = sum(row["timing_ms"][key]["baseline_minus_candidate"]
                              for key in PHASES)
            if abs(online_delta - phase_delta) > 0.02:
                raise ValueError(f"exclusive phase difference mismatch {n}/{block['label']}")
            target_rows.append(row)
        if len(target_rows) != 12:
            raise ValueError(f"expected 12 frozen targets on N{n}")

        timings = {}
        for key in (ONLINE,) + PHASES + NESTED:
            rows = [row["timing_ms"][key] for row in target_rows]
            timings[key] = {
                "mean_baseline_ms": statistics.mean(row["baseline"] for row in rows),
                "mean_candidate_ms": statistics.mean(row["candidate"] for row in rows),
                **paired_mean_interval([row["baseline_minus_candidate"] for row in rows]),
            }
        counters = {}
        for key in COUNTERS:
            rows = [row["counters"][key] for row in target_rows]
            counters[key] = {
                "mean_baseline": statistics.mean(row["baseline"] for row in rows),
                "mean_candidate": statistics.mean(row["candidate"] for row in rows),
                "mean_candidate_minus_baseline": statistics.mean(
                    row["candidate"] - row["baseline"] for row in rows),
            }
        online = timings[ONLINE]["mean_baseline_minus_candidate_ms"]
        pdp = timings["target_pdp_charged"]["mean_baseline_minus_candidate_ms"]
        curves[str(n)] = {
            "target_count": len(target_rows), "phase_timing_ms": timings,
            "counters": counters,
            "pdp_share_of_mean_online_difference": pdp / online if online else None,
            "target_rows": target_rows,
        }

    result = {
        "kind": "modal_vm_s3_mechanism_diagnostic_v1",
        "scope": "Same 24 targets as the VM panel; descriptive paired phase attribution only",
        "candidate_id": None, "workload_id": None,
        "raw_sha256": sha(compressed), "custody_audit_sha256": sha(audit_bytes),
        "source_sha256": panel["host"]["source_sha256"],
        "strict_host_isolation_receipt": None, "controlled_speedup": None,
        "confirmatory_targets_per_curve": None,
        "uncertainty_status": "descriptive t intervals; N53 A/A and physical-host gates failed",
        "nested_timing_note": "rank_* costs are nested within online phases; index_build is reusable setup outside the online interval; do not add these diagnostics to the exclusive phases",
        "s3_counter_note": "rank_target_s3_calls counts logical S3 queries; paired roots can share a generic inversion, but fallback frequency was not measured",
        "curves": curves,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(output), "raw_sha256": result["raw_sha256"],
                      "pdp_share": {n: curves[n]["pdp_share_of_mean_online_difference"]
                                    for n in ("41", "53")}}, sort_keys=True))


if __name__ == "__main__":
    main()
