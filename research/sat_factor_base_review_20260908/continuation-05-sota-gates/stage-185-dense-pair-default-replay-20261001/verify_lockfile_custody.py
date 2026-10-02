#!/usr/bin/env python3
"""Verify the additive Stage-185 lockfile custody correction."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
correction = json.loads((ROOT / "lockfile-custody-correction.json").read_text())
lock = ROOT / "Cargo.lock"
result = ROOT / "result.json"
checks = {
    "schema": correction["schema"] == "koblitz_stage185_lockfile_custody_correction.v1",
    "status": correction["status"] == "evidence_scoped_lock_restored",
    "path": correction["lockfile_path"].endswith("stage-185-dense-pair-default-replay-20261001/Cargo.lock"),
    "bytes": lock.stat().st_size == correction["lockfile_bytes"] == 29115,
    "lock_sha256": hashlib.sha256(lock.read_bytes()).hexdigest() == correction["lockfile_sha256"],
    "result_sha256": hashlib.sha256(result.read_bytes()).hexdigest() == correction["selected_replay_result_sha256"],
    "root_not_selected": correction["root_lock_tracking_selected"] is False,
    "boundaries_unchanged": all(
        correction[key] is False
        for key in ("scientific_result_changed", "timing_changed", "routing_changed", "claim_changed")
    ),
}
value = {
    "schema": "koblitz_stage185_lockfile_custody_verification.v1",
    "checks": checks,
    "passed": all(checks.values()),
}
print(json.dumps(value, indent=2, sort_keys=True))
raise SystemExit(0 if value["passed"] else 1)
