#!/usr/bin/env python3
"""Derive per-run evidence from the immutable n=73/n=83 receipts.

The output retains the historical wall interval and records its limitations.
It does not turn a shared-host timing into a controlled speedup claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
EXPERIMENTS = {
    73: "experiments/koblitz-single-target-n73-20261003",
    83: "experiments/koblitz-single-target-n83-20261006",
}
IC_PHASES = (
    "target_query_ms",
    "target_pdp_ms",
    "target_relation_check_ms",
    "target_descent_ms",
    "target_recovery_check_ms",
)


def tracked_bytes(path: str) -> bytes:
    return subprocess.check_output(
        ["git", "show", f"HEAD:{path}"], cwd=ROOT, stderr=subprocess.DEVNULL
    )


def tracked_json(path: str) -> dict:
    return json.loads(tracked_bytes(path))


def one_jsonl(path: str) -> dict:
    rows = [json.loads(line) for line in tracked_bytes(path).splitlines() if line.strip()]
    if len(rows) != 1:
        raise ValueError(f"expected one JSON row in {path}, found {len(rows)}")
    return rows[0]


def run_paths(prefix: str) -> list[str]:
    names = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", "HEAD", f"{prefix}/runs"], cwd=ROOT
    ).decode().splitlines()
    paths = sorted(name.removesuffix("/run.json") for name in names if name.endswith("/run.json"))
    if len(paths) != 3:
        raise ValueError(f"expected three retained paired runs in {prefix}, found {len(paths)}")
    return paths


def derive(degree: int) -> dict:
    prefix = EXPERIMENTS[degree]
    fixture = tracked_json(f"{prefix}/frozen/fixture.json")
    candidate = tracked_json(f"{prefix}/candidate-manifest.json") if degree == 73 else None
    workload = tracked_json(f"{prefix}/workload.json") if degree == 73 else None
    target = fixture["public_target"]
    generator = fixture["generator"]
    scalar = int(fixture["fixture_scalar_validation_only"])
    rows = []
    for path in run_paths(prefix):
        names = ("run.json", "ic.jsonl", "rho.jsonl", "ic_execution.json", "rho_execution.json", "independent_replay.json")
        raw_hashes = {name: hashlib.sha256(tracked_bytes(f"{path}/{name}")).hexdigest() for name in names}
        run = tracked_json(f"{path}/run.json")
        ic = one_jsonl(f"{path}/ic.jsonl")
        rho = one_jsonl(f"{path}/rho.jsonl")
        ic_exec = tracked_json(f"{path}/ic_execution.json")
        rho_exec = tracked_json(f"{path}/rho_execution.json")
        replay = tracked_json(f"{path}/independent_replay.json")
        if run["status"] != "PRODUCERS_COMPLETE" or run["target_count"] != 1:
            raise ValueError(f"{path}: incomplete producer run or target count")
        if replay["status"] != "PASS" or not ic["group_verified"] or not rho["verified"]:
            raise ValueError(f"{path}: verification failed")
        if ic["published_q"] != target or rho["published_q"] != target:
            raise ValueError(f"{path}: target point differs from frozen fixture")
        if ic["generator"] != generator or rho["generator"] != generator:
            raise ValueError(f"{path}: generator differs from frozen fixture")
        if int(ic["recovered_scalar"]) != scalar or int(rho["recovered_fixture_scalar"]) != scalar:
            raise ValueError(f"{path}: recovered scalars differ from frozen fixture")
        phases = {name.removesuffix("_ms"): float(ic[name]) for name in IC_PHASES}
        ic_ms = float(ic["target_ms"])
        if not math.isclose(math.fsum(phases.values()), ic_ms, abs_tol=0.02):
            raise ValueError(f"{path}: IC phases do not sum to target interval")
        rho_phases = {name: float(rho[name]) for name in ("setup_ms", "walk_ms", "validation_ms")}
        rho_legacy_ms = math.fsum(rho_phases.values())
        if not math.isclose(rho_legacy_ms, float(rho["total_ms"]), abs_tol=0.02):
            raise ValueError(f"{path}: rho phases do not sum to recorded interval")
        rows.append({
            "legacy_run_id": run["run_id"],
            "candidate_id": candidate["candidate_id"] if candidate else None,
            "workload_id": workload["workload_id"] if workload else None,
            "public_target_q": target,
            "ic_online_ms": ic_ms,
            "ic_online_phase_ms": phases,
            "rho_recorded_online_ms": rho_legacy_ms,
            "rho_recorded_phase_ms": rho_phases,
            "rho_walk_plus_validation_ms": rho_phases["walk_ms"] + rho_phases["validation_ms"],
            "recorded_wall_ratio": rho_legacy_ms / ic_ms,
            "native_counters": {"ic_s3_root_index_probes": int(ic["probes"]), "rho_walk_steps": int(rho["walk_steps"])},
            "ic_binary_sha256": run["ic_binary_sha256"],
            "rho_binary_sha256": run.get("rho_binary_sha256", run.get("rho_reference_binary_sha256")),
            "ic_peak_rss_bytes": ic_exec["peak_rss_bytes"],
            "rho_peak_rss_bytes": rho_exec["peak_rss_bytes"],
            "common_resource_cap_bytes": run["resource_cap_bytes_per_arm"],
            "ic_rank_seed": run["ic_rank_seed"],
            "rho_walk_seed": run["rho_walk_seed"],
            "independent_replay": "PASS",
            "source_sha256": raw_hashes,
        })
    return {
        "schema_version": 2,
        "stage": "vs_rho",
        "field_degree": degree,
        "target_count": 1,
        "public_target_q": target,
        "public_generator": generator,
        "record_class": "verified_answer_exploratory_wall",
        "host_isolation_receipt": None,
        "controlled_online_speedup": None,
        "rho_interval_note": "The legacy interval includes setup_ms; its target dependence is not established by these receipts. A new paired run must freeze the exact online start event.",
        "counter_note": "IC S3 root-index probes and rho walk steps are distinct native counters; their quotient is not an operation speedup.",
        "missing_for_promotion": [
            "qualifying physical-host isolation receipt",
            "verified matching resource envelope and exact rho online boundary",
            *(["canonical IC1 candidate and workload manifests"] if degree == 83 else []),
        ],
        "runs": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("degree", type=int, choices=tuple(EXPERIMENTS))
    args = parser.parse_args()
    report = derive(args.degree)
    output = ROOT / EXPERIMENTS[args.degree] / "derived_exploratory_claim_report.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(output)


if __name__ == "__main__":
    main()
