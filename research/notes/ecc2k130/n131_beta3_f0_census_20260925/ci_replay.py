#!/usr/bin/env python3
"""Frozen-source preflight and archive-only evidence replay."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "rotated_subspace_support_20260925" / "gate.py"
PDP_PARENT = HERE.parent / "rotated_pdp_corpus_20260925" / "verify.py"
LOCAL_FILES = ("INPUT.json", "count.py", "replay.py", "run.py", "ci_replay.py", "PROTOCOL.md")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def preflight() -> dict:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert frozen["domain"] == "ECC2K130-N131-BETA3-F0-CENSUS-20260925-v1"
    for name in LOCAL_FILES:
        assert sha(HERE / name) == frozen["sha256"][name], name
    assert sha(PARENT) == frozen["sha256"]["parent_gate.py"]
    assert sha(PDP_PARENT) == frozen["sha256"]["parent_pdp_verify.py"]
    data = json.loads((HERE / "INPUT.json").read_text())
    assert data["domain"] == frozen["domain"]
    assert data["pilot_masks"] == 32768 and data["full_masks"] == 2097152
    spec = importlib.util.spec_from_file_location("n131_f0_count", HERE / "count.py")
    assert spec is not None and spec.loader is not None
    producer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(producer)
    for ordinal in range(1 << 16):
        mask = ordinal ^ (ordinal >> 1)
        assert producer.gray_ordinal(mask) == ordinal
    # A partner within full V0 but outside the pilot prefix is not a pilot collision.
    partner = 40000
    assert producer.gray_ordinal(partner ^ (partner >> 1)) >= data["pilot_masks"]
    with tempfile.TemporaryDirectory(prefix="n131-f0-selftest-") as tmp:
        report = Path(tmp) / "selftest.json"
        result = subprocess.run([sys.executable, str(HERE / "replay.py"),
                                 "--input", str(HERE / "INPUT.json"),
                                 "--self-test", "--output", str(report)],
                                capture_output=True, text=True, timeout=120)
        assert result.returncode == 0, (result.stdout, result.stderr)
        control = json.loads(report.read_text())
        assert control["status"] == "pass"
        assert [r["n"] for r in control["small_controls"]] == [13, 19]
        assert control["n131_group_law_samples"]["sample_masks"] == 64
    return frozen


def evidence_replay(evidence: Path, frozen: dict) -> None:
    receipt_path = evidence / "receipt.json"
    assert receipt_path.is_file()
    assert (evidence / "receipt.json.sha256").read_text().strip() == sha(receipt_path)
    receipt = json.loads(receipt_path.read_text())
    assert receipt["source_sha256"] == frozen["sha256"]
    assert receipt["frozen_sha256"] == sha(HERE / "FROZEN.json")
    assert receipt["status"] in {"complete", "failed",
                                 "full_censored_by_preregistered_pilot_projection",
                                 "full_censored_by_preregistered_pilot_rss"}
    actual = {str(path.relative_to(evidence)): sha(path)
              for path in sorted(evidence.rglob("*"))
              if path.is_file() and path.name not in ("receipt.json", "receipt.json.sha256")}
    assert receipt["raw_sha256"] == actual
    for stage in ("pilot", "full"):
        summary_path = evidence / stage / "summary.json"
        replay_path = evidence / stage / "replay.json"
        if not summary_path.exists():
            assert not replay_path.exists()
            continue
        summary = json.loads(summary_path.read_text())
        assert summary["stage"] == stage
        assert summary["total_masks"] == (32768 if stage == "pilot" else 2097152)
        assert summary["zero_x_count"] == 1 and summary["one_x_count"] == 0
        assert summary["zero_column_count"] == 1
        assert 4 * summary["nonzero_signed_columns"] >= summary["liftable_nonzero_x"]
        assert summary["nonzero_signed_columns"] <= summary["liftable_nonzero_x"]
        bins = {int(k): v for k, v in summary["multiplicity_histogram"].items()}
        assert set(bins) == {1, 2, 3, 4}
        assert sum(bins.values()) == summary["nonzero_signed_columns"]
        assert sum(k * v for k, v in bins.items()) == summary["liftable_nonzero_x"]
        if replay_path.exists():
            replay = json.loads(replay_path.read_text())
            assert replay["status"] == "pass"
            assert replay["measured"]["row_sha256"] == summary["row_sha256"]
            assert replay["measured"]["liftable_nonzero_x"] == summary["liftable_nonzero_x"]
            assert replay["measured"]["nonzero_signed_columns"] == summary["nonzero_signed_columns"]
        else:
            assert receipt["status"] == "failed"
    if receipt["status"] == "complete":
        assert (evidence / "pilot" / "replay.json").exists()
        assert (evidence / "full" / "replay.json").exists()
    print("n131 beta3 F0 source, controls and archive replay PASS")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    frozen = preflight()
    if args.evidence is None:
        print("n131 beta3 F0 frozen source and controls PASS")
    else:
        evidence_replay(args.evidence, frozen)


if __name__ == "__main__":
    main()
