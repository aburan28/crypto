#!/usr/bin/env python3
"""Fail-closed receipt for the bounded local N83 cold diagnostics."""

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PANEL = ROOT / "pilot-01"
LOGS = ROOT / "verification"
PUBLIC_SEED = 2026100801
ORDERED_BINARY_SHA256 = "38670992b1596edee12d50e3eb79452fdc9a6ee83e5b64e1864dfed5d8aafe92"
UNORDERED_BINARY_SHA256 = "786183f6becf192efa7df25680f469b855be73aabc306308bec3ca6fe63d0819"
CASES = (
    (256, "ordered", "cold-256.log", "cold-a1-hash-256-seed-2026100801-attempt-1", "8a7d14adf96048147ae448e23dacd8b21aac5a9b", ORDERED_BINARY_SHA256),
    (64, "ordered", "cold-64.log", "cold-a1-hash-64-seed-2026100801-attempt-1", "8a7d14adf96048147ae448e23dacd8b21aac5a9b", ORDERED_BINARY_SHA256),
    (64, "unordered", "cold-64-unordered.log", "cold-a1-hash-64-unordered-seed-2026100801-attempt-1", "b3603f9a987674bb8ec77bab6779f7ad23ebf2ae", UNORDERED_BINARY_SHA256),
)


def load_json(path: Path):
    return json.loads(path.read_text())


def elapsed_seconds(log: str) -> float:
    matches = re.findall(r"^real (\d+(?:\.\d+)?)$", log, re.MULTILINE)
    assert len(matches) == 1, "expected one complete /usr/bin/time receipt"
    return float(matches[0])


def summarize():
    manifest = load_json(PANEL / "manifest.json")
    replay = load_json(PANEL / "replay.json")
    probes = load_json(PANEL / "probes.json")
    upload = load_json(PANEL / "upload-receipt.json")
    assert manifest["status"] == "completed_factor_base_panel"
    assert replay["status"] == upload["status"] == "PASS"
    assert probes["status"] == "completed_relation_stage_panel"
    assert len(manifest["bases"]) == len(upload["objects"]) == 54

    stages = {
        "construction_internal_seconds": manifest["elapsed_seconds_informational"],
        "replay_process_wall_seconds": elapsed_seconds((LOGS / "replay-01.log").read_text()),
        "probes_process_wall_seconds": elapsed_seconds((LOGS / "probes-01.log").read_text()),
        "upload_process_wall_seconds": elapsed_seconds((LOGS / "upload-01.log").read_text()),
    }
    cases = []
    for columns, pair_mode, log_name, directory, source_commit, binary_sha256 in CASES:
        run_dir = PANEL / directory
        cap = load_json(run_dir / "cap.json")
        assert cap["status"] == "UNKNOWN_budget" and cap["budget_seconds"] == 600
        records = [json.loads(line) for line in (run_dir / "cold-run.jsonl").read_text().splitlines() if line]
        assert not any(record.get("kind") == "compact_orbit_dlp_summary" for record in records), "cap and complete summary disagree"
        log = (LOGS / log_name).read_text()
        phases = [json.loads(line.removeprefix("icv1_phase ")) for line in log.splitlines() if line.startswith("icv1_phase ")]
        elapsed = elapsed_seconds(log)
        cases.append({
            "columns": columns,
            "pair_mode": pair_mode,
            "source_commit": source_commit,
            "binary_sha256": binary_sha256,
            "run_directory": directory,
            "status": cap["status"],
            "solver_cap_seconds": 600,
            "process_wall_seconds": elapsed,
            "completed_pipeline_record_count": 0,
            "partial_record_count": len(records),
            "phase_checkpoints": phases,
        })

    approximate_total = sum(stages.values()) + sum(case["process_wall_seconds"] for case in cases)
    assert approximate_total <= 3600, "local pilot active execution exceeded the authorized cap"
    return {
        "schema": "n83.cold-diagnostics/v1",
        "study": manifest["study"],
        "panel_manifest_blake3": replay["panel_manifest_blake3"],
        "curve_arm": "diagnostic_a1_53bit",
        "public_seed": PUBLIC_SEED,
        "source_stage_times_seconds": stages,
        "cases": cases,
        "approximate_active_execution_seconds": round(approximate_total, 3),
        "time_basis_note": "construction is the producer's internal elapsed time; later stages are /usr/bin/time process wall; the sum is an approximate pilot-budget audit",
        "local_pilot_cap_seconds": 3600,
        "primary_a0_complete_cold_runs": 0,
        "selected_best_total_runtime": None,
    }


def main():
    path = PANEL / "cold-diagnostics.json"
    encoded = (json.dumps(summarize(), indent=2, sort_keys=True) + "\n").encode()
    if path.exists():
        assert path.read_bytes() == encoded, "existing receipt disagrees with immutable source logs"
    else:
        with path.open("xb") as output:
            output.write(encoded)
    print(path)


if __name__ == "__main__":
    main()
