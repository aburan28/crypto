#!/usr/bin/env python3
"""Verify Stage 179 receipts, tests, mechanism, and decision."""

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


check("schema", result["schema"] == "koblitz_stage179_current_f4_wide_tables_single_target.v1")
check("decision_rejected", result["decision"]["status"] == "REJECTED_WIDER_TABLE_REGRESSION")
check("five_selected", result["decision"]["selected_table_columns"] == 5)
check("tests_passed", all(item["passed"] for item in result["tests"].values()))
check("runs_correct", all(run["correct"] for run in result["runs"]))
check("performed_regressed", result["mechanism"]["performed_xor_ratio_six_over_five"] > 1)
check("table_memory_regressed", result["mechanism"]["table_memory_ratio_six_over_five"] > 1)
check("single_core_null", result["single_core_seconds"] is None)
check("conflicts_null", result["conflicts"] is None)
for group in ("builds", "tests"):
    for item in result[group].values():
        for receipt in item["artifacts"].values():
            artifact = STAGE / receipt["path"]
            check(f"exists:{receipt['path']}", artifact.is_file())
            check(f"sha:{receipt['path']}", sha256(artifact) == receipt["sha256"])
for run in result["runs"]:
    for receipt in run["artifacts"].values():
        artifact = STAGE / receipt["path"]
        check(f"run_exists:{receipt['path']}", artifact.is_file())
        check(f"run_sha:{receipt['path']}", sha256(artifact) == receipt["sha256"])
for field in ("wall", "core", "rss"):
    check(
        f"paired_median:{field}",
        result["median_paired_ratios"][field]
        == statistics.median(pair[field] for pair in result["pair_ratios"]),
    )
verification = {
    "schema": "koblitz_stage179_verification.v1",
    "checks": len(checks),
    "passed": sum(ok for _, ok in checks),
    "failed": [name for name, ok in checks if not ok],
    "result_sha256": sha256(path),
}
(STAGE / "verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")
print(json.dumps(verification, indent=2))
raise SystemExit(0 if not verification["failed"] else 1)
