#!/usr/bin/env python3
"""Verify Stage 176 raw receipts, ratios, and negative decision."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


STAGE = Path(__file__).resolve().parent
result_path = STAGE / "result.json"
result = json.loads(result_path.read_text())
checks: list[tuple[str, bool]] = []


def check(name: str, ok: bool) -> None:
    checks.append((name, bool(ok)))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


check("schema", result["schema"] == "koblitz_stage176_current_f4_x1_schedule.v1")
check("negative_decision", result["decision"]["status"] == "REJECTED_NO_JOINT_SCREEN_WIN")
check("no_joint_winner", result["decision"]["joint_wall_and_core_winners"] == [])
check("no_confirmation", result["decision"]["confirmation_run_required"] is False)
check("single_core_null", result["single_core_seconds"] is None)
check("conflicts_null", result["conflicts"] is None)
check("correctness", all(item["pass"] for item in result["correctness_checks"]))
base = result["batch512_stage175_median"]
for run in result["runs"]:
    metrics = run["process"]["metrics"]
    check(f"not_joint:{run['batch']}", not (metrics["wall_seconds"] < base["wall_seconds"] and metrics["total_core_seconds"] < base["total_core_seconds"]))
    for field, key in (("wall_seconds", "wall"), ("total_core_seconds", "core"), ("peak_rss_bytes", "rss")):
        check(f"ratio:{run['batch']}:{key}", run["ratios_over_batch512_median"][key] == metrics[field] / base[field])
    for receipt in run["artifacts"].values():
        path = STAGE / receipt["path"]
        check(f"exists:{receipt['path']}", path.is_file())
        check(f"size:{receipt['path']}", path.stat().st_size == receipt["bytes"])
        check(f"sha:{receipt['path']}", sha256(path) == receipt["sha256"])
verification = {
    "schema": "koblitz_stage176_verification.v1",
    "checks": len(checks),
    "passed": sum(ok for _, ok in checks),
    "failed": [name for name, ok in checks if not ok],
    "result_sha256": sha256(result_path),
}
(STAGE / "verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")
print(json.dumps(verification, indent=2))
raise SystemExit(0 if not verification["failed"] else 1)
