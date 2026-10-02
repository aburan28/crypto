#!/usr/bin/env python3
"""Verify Stage 182 single-core receipts and rejection."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import statistics


STAGE = Path(__file__).resolve().parent
path = STAGE / "result.json"
result = json.loads(path.read_text())
checks: list[tuple[str, bool]] = []


def check(name: str, ok: bool) -> None:
    checks.append((name, bool(ok)))


def sha256(item: Path) -> str:
    return hashlib.sha256(item.read_bytes()).hexdigest()


check("schema", result["schema"] == "koblitz_stage182_current_f4_full_m4ri_single_core.v1")
check("screen_correct", all(run["correct"] for run in result["screen"]["runs"]))
check("confirmation_correct", all(run["correct"] for run in result["confirmation"]["runs"]))
check("screen_continued", result["screen"]["continued"] is True)
check("wall_regressed", result["confirmation"]["median_paired_ratios"]["wall"] > 1)
check("cpu_regressed", result["confirmation"]["median_paired_ratios"]["core"] > 1)
check("decision_rejected", result["decision"]["status"] == "REJECTED_SINGLE_CORE_TIMING_REGRESSION")
check("no_auto_policy", result["decision"]["automatic_single_core_policy"] is False)
check("default_current", result["decision"]["runtime_default"] == "current_blocktables")
check("conflicts_null", result["conflicts"] is None)
for group in ("screen", "confirmation"):
    for run in result[group]["runs"]:
        metrics = run["process"]["metrics"]
        check(f"single_core:{group}:{run['order']}", metrics["single_core_seconds"] == metrics["total_core_seconds"])
        for receipt in run["artifacts"].values():
            artifact = STAGE / receipt["path"]
            check(f"exists:{receipt['path']}", artifact.is_file())
            check(f"sha:{receipt['path']}", sha256(artifact) == receipt["sha256"])
for field in ("wall", "core", "rss"):
    check(
        f"paired_median:{field}",
        result["confirmation"]["median_paired_ratios"][field]
        == statistics.median(pair[field] for pair in result["confirmation"]["pair_ratios"]),
    )
verification = {
    "schema": "koblitz_stage182_verification.v1",
    "checks": len(checks),
    "passed": sum(ok for _, ok in checks),
    "failed": [name for name, ok in checks if not ok],
    "result_sha256": sha256(path),
}
(STAGE / "verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")
print(json.dumps(verification, indent=2))
raise SystemExit(0 if not verification["failed"] else 1)
