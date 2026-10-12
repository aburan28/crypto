#!/usr/bin/env python3
"""Compose the Stage 176 schedule-screen result."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


STAGE = Path(__file__).resolve().parent
SCREEN = STAGE / "development" / "screen"
BATCHES = (1, 4, 12, 24, 39, 64, 128, 256)
EXPECTED_EQUATIONS = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"
EXPECTED_OPS = 319_313_687_585
EXPECTED_PERFORMED = 147_794_583_858
BASELINE = {
    "wall_seconds": 34.68066983300014,
    "total_core_seconds": 267.171734,
    "peak_rss_bytes": 4_030_480_384,
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def receipt(path: Path) -> dict:
    return {"path": str(path.relative_to(STAGE)), "bytes": path.stat().st_size, "sha256": sha256(path)}


runs = []
checks = []
for batch in BATCHES:
    root = SCREEN / f"batch-{batch}"
    metrics_path = root / "metrics.json"
    stdout_path = root / "stdout.json"
    stderr_path = root / "stderr.txt"
    process = json.loads(metrics_path.read_text())
    report = json.loads(stdout_path.read_text())
    extra = report["cost"]["extra"]
    ok = (
        process["returncode"] == 0
        and not process["timed_out"]
        and report["status"] == "unsat"
        and report["exhaustive"] is True
        and report["source_instance_verified"] is True
        and report["regenerated_source_exact"] is True
        and report["fixed_x1_masks_visited"] == 512
        and report["fixed_x1_systems_constructed"] == 242
        and report["fixed_x1_systems_completed"] == 242
        and report["solver_equations_blake3"] == EXPECTED_EQUATIONS
        and report["cost"]["ops"] == EXPECTED_OPS
        and extra["word_xors_performed"] == EXPECTED_PERFORMED
        and report["conflicts"] is None
    )
    checks.append({"batch": batch, "pass": ok})
    metrics = process["metrics"]
    runs.append(
        {
            "batch": batch,
            "process": process,
            "report": report,
            "ratios_over_batch512_median": {
                "wall": metrics["wall_seconds"] / BASELINE["wall_seconds"],
                "core": metrics["total_core_seconds"] / BASELINE["total_core_seconds"],
                "rss": metrics["peak_rss_bytes"] / BASELINE["peak_rss_bytes"],
            },
            "artifacts": {
                "metrics": receipt(metrics_path),
                "stdout": receipt(stdout_path),
                "stderr": receipt(stderr_path),
            },
        }
    )
if not all(item["pass"] for item in checks):
    raise SystemExit("screen correctness failure")

joint = [
    run
    for run in runs
    if run["process"]["metrics"]["wall_seconds"] < BASELINE["wall_seconds"]
    and run["process"]["metrics"]["total_core_seconds"] < BASELINE["total_core_seconds"]
]
best_core = min(runs, key=lambda run: (run["process"]["metrics"]["total_core_seconds"], run["process"]["metrics"]["wall_seconds"]))
charge = {
    "components": len(runs),
    "wall_seconds_sum": sum(run["process"]["metrics"]["wall_seconds"] for run in runs),
    "total_core_seconds_sum": sum(run["process"]["metrics"]["total_core_seconds"] for run in runs),
    "peak_rss_bytes_max": max(run["process"]["metrics"]["peak_rss_bytes"] for run in runs),
}
result = {
    "schema": "koblitz_stage176_current_f4_x1_schedule.v1",
    "claim_boundary": "One opened n=59 decomposition target and a scheduling diagnostic only; not a full index-calculus run or SOTA.",
    "source_commit": "c38626a5e2674dbfdff796474ea83af0c863917a",
    "binary_sha256": "b156f3320564c1e9fed2b168aba7a0b4f3c875d677c1a9354751bba527ee5bd6",
    "batch512_stage175_median": BASELINE,
    "runs": runs,
    "correctness_checks": checks,
    "best_core_diagnostic": {
        "batch": best_core["batch"],
        "metrics": best_core["process"]["metrics"],
        "ratios": best_core["ratios_over_batch512_median"],
    },
    "decision": {
        "status": "REJECTED_NO_JOINT_SCREEN_WIN",
        "joint_wall_and_core_winners": [run["batch"] for run in joint],
        "confirmation_run_required": False,
        "reason": "No screen arm improved both wall and total core-seconds over the Stage 175 batch-512 medians.",
    },
    "campaign_charge": charge,
    "single_core_seconds": None,
    "conflicts": None,
    "provenance": receipt(SCREEN / "provenance.json"),
}
(STAGE / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

lines = [
    "# Stage 176: current F4 fixed-X1 schedule screen",
    "",
    "All eight arms used the exact Stage 175 binary and target. Every arm returned exhaustive UNSAT with the same equation fingerprint and exact F4 operation counts.",
    "",
    "| X1 batch | wall (s) | core (s) | peak RSS (bytes) | wall / 512 | core / 512 |",
    "|---:|---:|---:|---:|---:|---:|",
]
for run in runs:
    m = run["process"]["metrics"]
    q = run["ratios_over_batch512_median"]
    lines.append(f"| {run['batch']} | {m['wall_seconds']:.6f} | {m['total_core_seconds']:.6f} | {m['peak_rss_bytes']} | {q['wall']:.3f} | {q['core']:.3f} |")
lines += [
    "",
    f"No arm improved both metrics, so confirmation was not run. Batch {best_core['batch']} had the lowest one-run CPU at {best_core['process']['metrics']['total_core_seconds']:.6f} core-seconds but took {best_core['process']['metrics']['wall_seconds']:.6f} wall seconds.",
    "",
    f"The eight-screen charge is {charge['wall_seconds_sum']:.6f} wall seconds, {charge['total_core_seconds_sum']:.6f} core-seconds, and {charge['peak_rss_bytes_max']} bytes maximum RSS.",
    "",
    "Schedule tuning is rejected. This does not change any SOTA gate.",
    "",
]
(STAGE / "RESULTS.md").write_text("\n".join(lines))
print(json.dumps({"decision": result["decision"], "best_core_diagnostic": result["best_core_diagnostic"], "campaign_charge": charge}, indent=2))
