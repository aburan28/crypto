#!/usr/bin/env python3
"""Verify Stage 186 screen receipts and stop decision."""

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


check("schema", result["schema"] == "koblitz_stage186_dense_pair_elimination_policy.v1")
check("runs_correct", all(run["correct"] for run in result["screen"]["runs"]))
check("performed_improved", result["screen"]["full_over_blocktables"]["performed_xors"] < 1)
check("core_regressed", result["screen"]["full_over_blocktables"]["core"] > 1)
check("wall_regressed", result["screen"]["full_over_blocktables"]["wall"] > 1)
check("no_confirmation", result["screen"]["continued"] is False)
check("decision_rejected", result["decision"]["status"] == "REJECTED_SCREEN_CPU_REGRESSION")
check("blocktables_selected", result["decision"]["selected_phase_b_elimination"] == "current_five_column_blocktables")
check("dense_selected", result["decision"]["selected_pair_selector"] == "dense_exact")
check("single_core_null", result["single_core_seconds"] is None)
check("conflicts_null", result["conflicts"] is None)
for run in result["screen"]["runs"]:
    for receipt in run["artifacts"].values():
        artifact = STAGE / receipt["path"]
        check(f"exists:{receipt['path']}", artifact.is_file())
        check(f"sha:{receipt['path']}", sha256(artifact) == receipt["sha256"])
verification = {
    "schema": "koblitz_stage186_verification.v1",
    "checks": len(checks),
    "passed": sum(ok for _, ok in checks),
    "failed": [name for name, ok in checks if not ok],
    "result_sha256": sha256(path),
}
(STAGE / "verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")
print(json.dumps(verification, indent=2))
raise SystemExit(0 if not verification["failed"] else 1)
