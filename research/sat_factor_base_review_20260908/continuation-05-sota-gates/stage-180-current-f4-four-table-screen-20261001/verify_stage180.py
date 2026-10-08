#!/usr/bin/env python3
"""Verify Stage 180 receipts and stop decision."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


STAGE = Path(__file__).resolve().parent
path = STAGE / "result.json"
result = json.loads(path.read_text())
checks: list[tuple[str, bool]] = []


def check(name: str, ok: bool) -> None:
    checks.append((name, bool(ok)))


def sha256(item: Path) -> str:
    return hashlib.sha256(item.read_bytes()).hexdigest()


check("schema", result["schema"] == "koblitz_stage180_current_f4_four_table_screen.v1")
check("decision_rejected", result["decision"]["status"] == "REJECTED_SCREEN_WORK_REGRESSION")
check("no_confirmation", result["decision"]["confirmation_required"] is False)
check("five_selected", result["decision"]["selected_table_columns"] == 5)
check("runs_correct", all(run["correct"] for run in result["runs"]))
check("performed_regressed", result["mechanism"]["performed_xor_ratio_four_over_five"] > 1)
check("table_memory_reduced", result["mechanism"]["table_memory_ratio_four_over_five"] < 1)
check("single_core_null", result["single_core_seconds"] is None)
for run in result["runs"]:
    for receipt in run["artifacts"].values():
        artifact = STAGE / receipt["path"]
        check(f"exists:{receipt['path']}", artifact.is_file())
        check(f"sha:{receipt['path']}", sha256(artifact) == receipt["sha256"])
for receipt in result["builds"]["four"]["artifacts"].values():
    artifact = STAGE / receipt["path"]
    check(f"build_sha:{receipt['path']}", sha256(artifact) == receipt["sha256"])
inherited = STAGE / result["builds"]["inherited_five"]["receipt_path"]
check("inherited_build_sha", sha256(inherited) == result["builds"]["inherited_five"]["receipt_sha256"])
verification = {
    "schema": "koblitz_stage180_verification.v1",
    "checks": len(checks),
    "passed": sum(ok for _, ok in checks),
    "failed": [name for name, ok in checks if not ok],
    "result_sha256": sha256(path),
}
(STAGE / "verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")
print(json.dumps(verification, indent=2))
raise SystemExit(0 if not verification["failed"] else 1)
