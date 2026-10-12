#!/usr/bin/env python3
"""Read-only source-lock audit; it must never construct a candidate or Q."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PINNED = {
    ".github/workflows/ecc2k130-m3-four-policy.yml",
    ".github/workflows/ecc2k130-m3-base-selector-source.yml",
    "research/notes/ecc2k130/m3_base_selector_20260930/CONFIG.json",
    "research/notes/ecc2k130/m3_base_selector_20260930/PROTOCOL.md",
    "research/notes/ecc2k130/m3_base_selector_20260930/SOURCE_LOCK.md",
    "research/notes/ecc2k130/m3_base_selector_20260930/check_protocol.py",
    "research/notes/ecc2k130/m3_base_selector_20260930/check_source_lock.py",
    "research/notes/ecc2k130/m3_base_selector_20260930/produce.py",
    "research/notes/ecc2k130/m3_base_selector_20260930/score.py",
    "research/notes/ecc2k130/m3_base_selector_20260930/test_score.py",
    "research/notes/ecc2k130/m3_base_selector_20260930/verify.py",
    "research/notes/ecc2k130/m3_four_policy_20260930/CONFIG.json",
    "research/notes/ecc2k130/m3_four_policy_20260930/FROZEN.json",
    "research/notes/ecc2k130/m3_four_policy_20260930/FROZEN_REPLAY.json",
    "research/notes/ecc2k130/m3_four_policy_20260930/PROTOCOL.md",
    "research/notes/ecc2k130/m3_four_policy_20260930/REPLAY_REPAIR_PROTOCOL.md",
    "research/notes/ecc2k130/m3_four_policy_20260930/produce.py",
    "research/notes/ecc2k130/m3_four_policy_20260930/verify.py",
    "research/notes/ecc2k130/m3_four_policy_20260930/test_arithmetic.py",
    "research/notes/ecc2k130/m3_four_policy_20260930/evidence_run_36722040881/result.json",
    "research/ecc2k130_factor_base_pilot_20260924/run.py",
    "research/ecc2k130_factor_base_pilot_20260924/results_final.json",
    "research/ecc2k130_factor_base_replication_20260925/run.py",
    "research/ecc2k130_dual_transport_20260925/dual_transport.py",
    "research/ecc2k130_relations/fastfield.py",
    "research/ecc2k130_relations/relations.py",
    "research/ecc2k130_oriented_transport_20260924/oriented_velu.py",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check() -> dict:
    lock = json.loads((HERE / "FROZEN.json").read_text())
    assert lock["schema"] == "ecc2k130-degree7-m3-base-selector-source-lock-v1"
    assert lock["pre_outcome"] is True
    assert lock["sha256"].keys() == PINNED
    assert lock["protocol_sha256"] == sha(HERE / "PROTOCOL.md")
    for relative, digest in lock["sha256"].items():
        assert sha(ROOT / relative) == digest, relative
    return {"status": "PASS_SOURCE_LOCK", "pinned_files": len(PINNED),
            "source_lock_sha256": sha(HERE / "FROZEN.json"),
            "fixtures_generated_by_check": False}


if __name__ == "__main__":
    print(json.dumps(check(), indent=2, sort_keys=True))
