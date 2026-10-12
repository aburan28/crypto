#!/usr/bin/env python3
"""Verify Stage 183 receipts, counters, and accepted stacked decision."""

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


check("schema", result["schema"] == "koblitz_stage183_current_f4_dense_pair_selection.v1")
check("tests_passed", all(item["passed"] for item in result["tests"].values()))
check("screen_correct", all(run["correct"] for run in result["screen"]["runs"]))
check("confirmation_correct", all(run["correct"] for run in result["confirmation"]["runs"]))
check("screen_continued", result["screen"]["continued"] is True)
check("wall_passed", result["confirmation"]["median_paired_ratios"]["wall"] < 0.97)
check("cpu_passed", result["confirmation"]["median_paired_ratios"]["core"] < 0.97)
check("accepted", result["decision"]["status"] == "ACCEPTED_FOR_FULL_M4RI_STACK")
check("default_unchanged", result["decision"]["current_runtime_default"] == "quadratic_with_blocktables")
check("single_core_null", result["single_core_seconds"] is None)
check("conflicts_null", result["conflicts"] is None)
for group in ("screen", "confirmation"):
    for run in result[group]["runs"]:
        for receipt in run["artifacts"].values():
            artifact = STAGE / receipt["path"]
            check(f"exists:{receipt['path']}", artifact.is_file())
            check(f"sha:{receipt['path']}", sha256(artifact) == receipt["sha256"])
for item in [result["build"], *result["tests"].values()]:
    for receipt in item["artifacts"].values():
        artifact = STAGE / receipt["path"]
        check(f"setup_exists:{receipt['path']}", artifact.is_file())
        check(f"setup_sha:{receipt['path']}", sha256(artifact) == receipt["sha256"])
for field in ("wall", "core", "rss"):
    check(
        f"paired_median:{field}",
        result["confirmation"]["median_paired_ratios"][field]
        == statistics.median(pair[field] for pair in result["confirmation"]["pair_ratios"]),
    )
verification = {
    "schema": "koblitz_stage183_verification.v1",
    "checks": len(checks),
    "passed": sum(ok for _, ok in checks),
    "failed": [name for name, ok in checks if not ok],
    "result_sha256": sha256(path),
}
(STAGE / "verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")
print(json.dumps(verification, indent=2))
raise SystemExit(0 if not verification["failed"] else 1)
