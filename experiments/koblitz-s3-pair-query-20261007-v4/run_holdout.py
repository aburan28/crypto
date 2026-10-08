#!/usr/bin/env python3
"""Run each fresh public point as its own paired one-target workload."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import time

from run_panel import canonical, semantic_digest


HERE = Path(__file__).resolve().parent
HOLDOUT = HERE / "holdout"
RUNS = HOLDOUT / "runs"
ROWS = HOLDOUT / "measurement_rows_raw.jsonl"
SUMMARY = HOLDOUT / "paired_summary_exploratory.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def main() -> None:
    if ROWS.exists() or SUMMARY.exists():
        raise FileExistsError("holdout measurements already exist")
    fixtures_path = HOLDOUT / "fixtures.json"
    fixtures = load(fixtures_path)
    assert fixtures["count"] == len(fixtures["fixtures"]) == 12
    freeze = load(HERE / "freeze_receipt.json")
    manifests = {}
    for variant in ("baseline", "candidate"):
        path = HERE / f"{variant}_manifest.json"
        manifest = load(path)
        frozen = manifest["candidate_freeze"]
        assert sha(path) == freeze["variants"][variant]["manifest_sha256"]
        assert sha(Path(frozen["source_path"])) == frozen["source_sha256"]
        assert sha(Path(frozen["binary_path"])) == frozen["binary_sha256"]
        manifests[variant] = manifest
    RUNS.mkdir(exist_ok=True)
    rows = []
    for fixture in fixtures["fixtures"]:
        label = fixture["label"]
        workload = load(HOLDOUT / label / "workload.json")
        assert workload["workload_id"] == fixture["workload_id"]
        assert workload["identity_sha256"] == hashlib.sha256(
            canonical(workload["record"])
        ).hexdigest()
        assert workload["record"]["targets"] == [fixture["public_point"]]
        target_path = HOLDOUT / label / "public_target.json"
        assert sha(target_path) == fixture["public_target_sha256"]
        order = ("baseline", "candidate") if int(label[1:]) % 2 else (
            "candidate", "baseline"
        )
        for within_pair, variant in enumerate(order, 1):
            manifest = manifests[variant]
            frozen = manifest["candidate_freeze"]
            run_id = f"{manifest['candidate_id']}W{workload['workload_id']}R1"
            stem = RUNS / f"{label}_{variant}"
            output = stem.with_suffix(".jsonl")
            stdout = stem.with_suffix(".stdout")
            stderr = stem.with_suffix(".stderr")
            for path in (output, stdout, stderr):
                if path.exists():
                    raise FileExistsError(f"refusing to replace raw run: {path}")
            command = [
                frozen["binary_path"], "53", "0", "244", "20260928",
                str(target_path), str(output), "14",
            ]
            started = time.time_ns()
            with stdout.open("wb") as out, stderr.open("wb") as err:
                process = subprocess.run(command, stdout=out, stderr=err, check=False)
            ended = time.time_ns()
            result = None
            parse_error = None
            if output.exists():
                try:
                    result = json.loads(output.read_text().splitlines()[-1])
                except (IndexError, ValueError) as exc:
                    parse_error = str(exc)
            reported_success = bool(
                process.returncode == 0 and result is not None
                and result.get("group_verified") is True
                and result.get("target") == fixture["public_point"]
                and result.get("recovered_scalar") == fixture["fixture_scalar"]
            )
            row = {
                "candidate_id": manifest["candidate_id"],
                "workload_id": workload["workload_id"],
                "run_id": run_id,
                "target_label": label,
                "variant": variant,
                "pair_order": within_pair,
                "status": (
                    "reported_success_pending_sage" if reported_success else "failed"
                ),
                "process_exit_code": process.returncode,
                "parse_error": parse_error,
                "start_unix_ns": started,
                "end_unix_ns": ended,
                "process_wall_ms_including_launch_and_output": (ended-started)/1e6,
                "target": fixture["public_point"],
                "public_target_sha256": fixture["public_target_sha256"],
                "manifest_sha256": sha(HERE / f"{variant}_manifest.json"),
                "binary_sha256": frozen["binary_sha256"],
                "resource_envelope": workload["record"]["resource_envelope"],
                "host": {
                    "platform": platform.platform(),
                    "machine": platform.machine(),
                    "cpu_count_visible": os.cpu_count(),
                    "isolation_receipt": None,
                },
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
                timing = result["timing_ms"]
                row.update({
                    "online_interval": (
                        "after reusable base and index setup through scalar point replay"
                    ),
                    "online_ms": timing["target_online_after_reusable_setup"],
                    "online_exclusive_phases_ms": {
                        phase: timing[phase] for phase in (
                            "target_query", "target_pdp_charged",
                            "target_relation_check", "target_descent",
                            "target_recovery_check",
                        )
                    },
                    "setup_ms": timing["reusable_setup_total"],
                    "cold_after_launch_ms": timing["cold_after_launch"],
                    "target_pdp_ms": timing["target_pdp"],
                    "rank_pdp_wall_ms": timing["rank_pdp_wall"],
                    "rank_pdp_worker_time_sum_ms": timing["rank_pdp_worker_time_sum"],
                    "base_construction_ms": timing["factor_base_construct"],
                    "index_build_ms": timing["index_build"],
                    "actual_factor_base_points": result["factor_base_points"],
                    "effective_columns": result["orbit_columns"],
                    "index_entries": result["root_table_entries"],
                    "attempts": result["rank_attempts"],
                    "failed_attempts": result["rank_failures"],
                    "verified_relations": result["rank_verified_relations"],
                    "novel_rows": result["rank_new_rows"],
                    "final_rank": result["rank"],
                    "target_state_probes": result["target_state_probes"],
                    "rank_state_probes": result["rank_state_probes"],
                    "target_s3_calls": result["target_s3_calls"],
                    "rank_s3_calls": result["rank_target_s3_calls"],
                    "recovered_scalar_reported": result["recovered_scalar"],
                    "solved_targets": int(reported_success),
                    "peak_rss_bytes": result["peak_rss_bytes"],
                })
            rows.append(row)
            with ROWS.open("a") as stream:
                stream.write(json.dumps(row, sort_keys=True) + "\n")
                stream.flush()
                os.fsync(stream.fileno())
            print(json.dumps({"label": label, "variant": variant,
                              "status": row["status"],
                              "online_ms": row.get("online_ms")}), flush=True)

    pairs = []
    for fixture in fixtures["fixtures"]:
        label = fixture["label"]
        pair = {row["variant"]: row for row in rows if row["target_label"] == label}
        equal = (
            all(pair[v]["status"] == "reported_success_pending_sage"
                for v in ("baseline", "candidate"))
            and pair["baseline"]["semantic_digest_excluding_clocks_memory_and_s3_counters"]
            == pair["candidate"]["semantic_digest_excluding_clocks_memory_and_s3_counters"]
        )
        pairs.append({
            "target_label": label,
            "workload_id": fixture["workload_id"],
            "semantic_equal": equal,
            "baseline_run_id": pair["baseline"]["run_id"],
            "candidate_run_id": pair["candidate"]["run_id"],
            "baseline_online_ms": pair["baseline"].get("online_ms"),
            "candidate_online_ms": pair["candidate"].get("online_ms"),
            "exploratory_online_ratio": (
                pair["baseline"]["online_ms"] / pair["candidate"]["online_ms"]
                if equal and pair["candidate"]["online_ms"] else None
            ),
        })
    good = [pair for pair in pairs if pair["semantic_equal"]]
    summary = {
        "kind": "s3_pair_fresh_one_target_panel_exploratory",
        "question": fixtures["question"],
        "target_count": len(fixtures["fixtures"]),
        "run_count": len(rows),
        "all_pairs_semantically_equal": len(good) == len(pairs),
        "verified_by_sage": False,
        "pairs": pairs,
        "observed_solved_target_pairs": len(good),
        "exploratory_paired_ratio_median": (
            statistics.median(pair["exploratory_online_ratio"] for pair in good)
            if good else None
        ),
        "controlled_speedup": None,
        "claim_limit": (
            "No host-wide CPU isolation receipt; these paired ratios are "
            "exploratory. Independent Sage replay is pending."
        ),
        "fixtures_sha256": sha(fixtures_path),
        "measurement_rows_raw_sha256": sha(ROWS),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, sort_keys=True), flush=True)
    if not summary["all_pairs_semantically_equal"]:
        raise SystemExit("holdout pair failed or changed semantics; see raw rows")


if __name__ == "__main__":
    main()
