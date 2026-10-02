#!/usr/bin/env python3
"""Verify Stage 178 receipts, pair ratios, and frozen decision."""

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


check("schema", result["schema"] == "koblitz_stage178_current_f4_dense_reducer_index.v1")
check("decision_rejected", result["decision"]["status"] == "REJECTED_CPU_THRESHOLD")
check("default_linear", result["decision"]["runtime_default"] == "linear")
check("runs_correct", all(run["correct"] for run in result["runs"]))
check("structural_identity", result["mechanism"]["structural_f4_counts_identical"] is True)
check("wall_passed", result["median_paired_ratios"]["wall"] < 0.97)
check("cpu_failed", result["median_paired_ratios"]["core"] >= 0.97)
check("single_core_null", result["single_core_seconds"] is None)
check("conflicts_null", result["conflicts"] is None)
for run in result["runs"]:
    for receipt in run["artifacts"].values():
        artifact = STAGE / receipt["path"]
        check(f"exists:{receipt['path']}", artifact.is_file())
        check(f"size:{receipt['path']}", artifact.stat().st_size == receipt["bytes"])
        check(f"sha:{receipt['path']}", sha256(artifact) == receipt["sha256"])
for receipt in result["build"]["artifacts"].values():
    artifact = STAGE / receipt["path"]
    check(f"build_sha:{receipt['path']}", sha256(artifact) == receipt["sha256"])
for field in ("wall", "core", "rss"):
    check(
        f"paired_median:{field}",
        result["median_paired_ratios"][field]
        == statistics.median(pair[field] for pair in result["pair_ratios"]),
    )
patch = STAGE / result["rejected_patch"]["path"]
check("patch_sha", sha256(patch) == result["rejected_patch"]["sha256"])
verification = {
    "schema": "koblitz_stage178_verification.v1",
    "checks": len(checks),
    "passed": sum(ok for _, ok in checks),
    "failed": [name for name, ok in checks if not ok],
    "result_sha256": sha256(path),
}
(STAGE / "verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")
print(json.dumps(verification, indent=2))
raise SystemExit(0 if not verification["failed"] else 1)
