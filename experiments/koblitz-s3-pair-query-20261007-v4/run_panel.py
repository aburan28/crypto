#!/usr/bin/env python3
"""Run the frozen, one-target S3 baseline/candidate panel in alternating order."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import time


HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs/panel"
ROWS = HERE / "measurement_rows_raw.jsonl"
SUMMARY = HERE / "paired_summary_exploratory.json"
VARIANTS = ("baseline", "candidate")
REPETITIONS = 5


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def run_number(repetition: int, variant: str) -> str:
    manifest = load(HERE / f"{variant}_manifest.json")
    workload = load(HERE / "workload.json")
    return f"{manifest['candidate_id']}W{workload['workload_id']}R{repetition}"


def semantic_digest(record: dict) -> str:
    """Retain exact solution and witness order, ignoring clocks and S3 counters."""
    ignored = {"timing_ms", "peak_rss_bytes", "target_s3_calls", "rank_target_s3_calls"}
    value = {key: item for key, item in record.items() if key not in ignored}
    value["rank_relation_witnesses"] = [
        {key: item for key, item in witness.items() if key != "target_s3_calls"}
        for witness in record["rank_relation_witnesses"]
    ]
    return hashlib.sha256(canonical(value)).hexdigest()


def validate_freeze() -> tuple[dict, dict]:
    freeze = load(HERE / "freeze_receipt.json")
    workload = load(HERE / "workload.json")
    assert freeze["workload_sha256"] == sha(HERE / "workload.json")
    assert freeze["workload_id"] == workload["workload_id"]
    assert workload["identity_sha256"] == hashlib.sha256(
        canonical(workload["record"])
    ).hexdigest()
    assert load(HERE / "sage_runtime_info.json")
    manifests = {}
    for variant in VARIANTS:
        path = HERE / f"{variant}_manifest.json"
        manifest = load(path)
        assert sha(path) == freeze["variants"][variant]["manifest_sha256"]
        assert manifest["candidate_id"] == freeze["variants"][variant]["candidate_id"]
        candidate = manifest["candidate_freeze"]
        assert sha(Path(candidate["source_path"])) == candidate["source_sha256"]
        assert sha(Path(candidate["binary_path"])) == candidate["binary_sha256"]
        assert sha(Path(candidate["target_path"])) == candidate["target_sha256"]
        assert candidate["target_sha256"] == freeze["target_sha256"]
        assert candidate["workload_id"] == workload["workload_id"]
        manifests[variant] = manifest
    return workload, manifests


def execute(variant: str, repetition: int, order: int, workload: dict,
            manifest: dict) -> dict:
    frozen = manifest["candidate_freeze"]
    run_id = run_number(repetition, variant)
    base = RUNS / f"R{repetition}_{variant}"
    output = base.with_suffix(".jsonl")
    stdout = base.with_suffix(".stdout")
    stderr = base.with_suffix(".stderr")
    for path in (output, stdout, stderr):
        if path.exists():
            raise FileExistsError(f"refusing to replace frozen run artifact: {path}")
    command = [
        frozen["binary_path"], "53", "0", "244", "20260928",
        frozen["target_path"], str(output), "14",
    ]
    started = time.time_ns()
    with stdout.open("wb") as out, stderr.open("wb") as err:
        process = subprocess.run(command, stdout=out, stderr=err, check=False)
    ended = time.time_ns()
    result = None
    parse_error = None
    if output.is_file():
        try:
            result = json.loads(output.read_text().splitlines()[-1])
        except (IndexError, ValueError) as exc:
            parse_error = str(exc)
    target = workload["record"]["targets"][0]
    reported_success = bool(
        process.returncode == 0 and result is not None
        and result.get("group_verified") is True
        and result.get("target") == target
        and result.get("recovered_scalar") is not None
    )
    row = {
        "candidate_id": manifest["candidate_id"],
        "workload_id": workload["workload_id"],
        "run_id": run_id,
        "variant": variant,
        "repetition": repetition,
        "execution_order": order,
        "status": "reported_success_pending_sage" if reported_success else "failed",
        "process_exit_code": process.returncode,
        "parse_error": parse_error,
        "execution_start_unix_ns": started,
        "execution_end_unix_ns": ended,
        "process_wall_ms_including_launch_and_output": (ended - started) / 1e6,
        "host": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "cpu_count_visible": os.cpu_count(),
            "isolation_receipt": None,
        },
        "resource_envelope": workload["record"]["resource_envelope"],
        "target": target,
        "target_sha256": frozen["target_sha256"],
        "manifest_sha256": sha(HERE / f"{variant}_manifest.json"),
        "binary_sha256": frozen["binary_sha256"],
        "artifacts": {
            "raw": str(output),
            "raw_sha256": sha(output) if output.exists() else None,
            "stdout": str(stdout),
            "stdout_sha256": sha(stdout),
            "stderr": str(stderr),
            "stderr_sha256": sha(stderr),
        },
        "sage_verified": False,
        "semantic_digest_excluding_clocks_memory_and_s3_counters": (
            semantic_digest(result) if reported_success else None
        ),
    }
    if result is not None:
        timings = result.get("timing_ms", {})
        row.update({
            "online_interval": (
                "after reusable setup through recovered-scalar replay; "
                "fixture, launch, input, base, and index setup excluded"
            ),
            "online_ms": timings.get("target_online_after_reusable_setup"),
            "online_exclusive_phases_ms": {
                key: timings.get(key) for key in (
                    "target_query", "target_pdp_charged", "target_relation_check",
                    "target_descent", "target_recovery_check",
                )
            },
            "setup_ms": timings.get("reusable_setup_total"),
            "cold_after_launch_ms": timings.get("cold_after_launch"),
            "target_pdp_ms": timings.get("target_pdp"),
            "rank_pdp_wall_ms": timings.get("rank_pdp_wall"),
            "rank_pdp_worker_time_sum_ms": timings.get("rank_pdp_worker_time_sum"),
            "base_construction_ms": timings.get("factor_base_construct"),
            "index_build_ms": timings.get("index_build"),
            "actual_factor_base_points": result.get("factor_base_points"),
            "effective_columns": result.get("orbit_columns"),
            "index_entries": result.get("root_table_entries"),
            "attempts": result.get("rank_attempts"),
            "failed_attempts": result.get("rank_failures"),
            "verified_relations": result.get("rank_verified_relations"),
            "novel_rows": result.get("rank_new_rows"),
            "final_rank": result.get("rank"),
            "target_state_probes": result.get("target_state_probes"),
            "rank_state_probes": result.get("rank_state_probes"),
            "target_s3_calls": result.get("target_s3_calls"),
            "rank_s3_calls": result.get("rank_target_s3_calls"),
            "solved_targets": int(reported_success),
            "recovered_scalar_reported": result.get("recovered_scalar"),
            "peak_rss_bytes": result.get("peak_rss_bytes"),
        })
    return row


def main() -> None:
    if ROWS.exists() or SUMMARY.exists():
        raise FileExistsError("panel artifacts already exist; preserve raw runs")
    workload, manifests = validate_freeze()
    RUNS.mkdir(parents=True, exist_ok=True)
    rows = []
    order = 0
    for repetition in range(1, REPETITIONS + 1):
        variants = VARIANTS if repetition % 2 else VARIANTS[::-1]
        for variant in variants:
            order += 1
            row = execute(variant, repetition, order, workload, manifests[variant])
            rows.append(row)
            with ROWS.open("a") as stream:
                stream.write(json.dumps(row, sort_keys=True) + "\n")
                stream.flush()
                os.fsync(stream.fileno())
            print(json.dumps({
                "run_id": row["run_id"], "status": row["status"],
                "online_ms": row.get("online_ms"),
                "rank_pdp_wall_ms": row.get("rank_pdp_wall_ms"),
            }), flush=True)

    pairs = []
    for repetition in range(1, REPETITIONS + 1):
        pair = {row["variant"]: row for row in rows if row["repetition"] == repetition}
        equivalent = (
            all(pair[variant]["status"] == "reported_success_pending_sage"
                for variant in VARIANTS)
            and pair["baseline"]["semantic_digest_excluding_clocks_memory_and_s3_counters"]
            == pair["candidate"]["semantic_digest_excluding_clocks_memory_and_s3_counters"]
        )
        pairs.append({
            "repetition": repetition,
            "baseline_run_id": pair["baseline"]["run_id"],
            "candidate_run_id": pair["candidate"]["run_id"],
            "semantic_equal": equivalent,
            "baseline_online_ms": pair["baseline"].get("online_ms"),
            "candidate_online_ms": pair["candidate"].get("online_ms"),
            "exploratory_online_ratio": (
                pair["baseline"]["online_ms"] / pair["candidate"]["online_ms"]
                if equivalent and pair["candidate"]["online_ms"] else None
            ),
        })
    successful = [pair for pair in pairs if pair["semantic_equal"]]
    summary = {
        "kind": "s3_pair_query_exploratory_panel",
        "workload_id": workload["workload_id"],
        "repetitions": REPETITIONS,
        "run_order": [row["variant"] for row in rows],
        "candidate_ids": {variant: manifests[variant]["candidate_id"] for variant in VARIANTS},
        "single_public_target": workload["record"]["targets"][0],
        "all_pairs_semantically_equal": len(successful) == REPETITIONS,
        "pairs": pairs,
        "online_ms_median": {
            variant: statistics.median(row["online_ms"] for row in rows
                                       if row["variant"] == variant
                                       and row["status"] == "reported_success_pending_sage")
            for variant in VARIANTS
        },
        "exploratory_paired_online_ratio_median": (
            statistics.median(pair["exploratory_online_ratio"] for pair in successful)
            if successful else None
        ),
        "promoted_speedup": None,
        "claim_limit": (
            "Host-wide CPU isolation and noise gates were not recorded; "
            "ratios are exploratory. Independent Sage replay is pending."
        ),
        "measurement_rows_sha256": sha(ROWS),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, sort_keys=True))
    if not summary["all_pairs_semantically_equal"]:
        raise SystemExit("at least one pair failed or changed semantics; inspect raw rows")


if __name__ == "__main__":
    main()
