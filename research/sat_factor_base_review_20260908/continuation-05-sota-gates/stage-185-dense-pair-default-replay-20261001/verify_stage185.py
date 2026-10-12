#!/usr/bin/env python3
"""Verify Stage 185 exact-commit replay and custody."""

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


check("schema", result["schema"] == "koblitz_stage185_dense_pair_default_replay.v1")
check("failed_build_preserved", result["failed_locked_build"]["process"]["returncode"] != 0)
check("build_passed", result["build"]["process"]["returncode"] == 0)
check("tests_passed", all(item["passed"] for item in result["tests"].values()))
check("native_correct", result["replay"]["default_native_correct"] is True)
check("direct_correct", result["replay"]["direct_correct"] is True)
check("single_core_null", result["replay"]["valid_single_core_seconds"] is None)
check("selection_pass", result["decision"]["status"] == "SELECTED_DEFAULT_REPLAY_PASS")
check("lock_restore_required", result["decision"]["lockfile_restore_required"] is True)
check("conflicts_null", result["conflicts"] is None)
for item in [result["failed_locked_build"], result["build"], *result["tests"].values()]:
    for receipt in item["artifacts"].values():
        artifact = STAGE / receipt["path"]
        check(f"exists:{receipt['path']}", artifact.is_file())
        check(f"sha:{receipt['path']}", sha256(artifact) == receipt["sha256"])
for item in result["replay"].values():
    if not isinstance(item, dict) or "artifacts" not in item:
        continue
    for receipt in item["artifacts"].values():
        artifact = STAGE / receipt["path"]
        check(f"replay_exists:{receipt['path']}", artifact.is_file())
        check(f"replay_sha:{receipt['path']}", sha256(artifact) == receipt["sha256"])
verification = {
    "schema": "koblitz_stage185_verification.v1",
    "checks": len(checks),
    "passed": sum(ok for _, ok in checks),
    "failed": [name for name, ok in checks if not ok],
    "result_sha256": sha256(path),
}
(STAGE / "verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")
print(json.dumps(verification, indent=2))
raise SystemExit(0 if not verification["failed"] else 1)
