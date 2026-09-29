#!/usr/bin/env python3
"""Hash-only preflight; optional committed-evidence receipt check."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert sha(Path(__file__)) == frozen["ci_replay_sha256"]
    pins = {
        "protocol_sha256": HERE / "PROTOCOL.md",
        "diagnostic_sha256": HERE / "diagnostic.py",
        "verify_sha256": HERE / "verify.py",
        "runner_sha256": HERE / "run.py",
        "parent_producer_sha256": NOTES / "rotated_subspace_support_20260925/gate.py",
        "corpus_verify_sha256": NOTES / "rotated_pdp_corpus_20260925/verify.py",
        "corpus_archive_sha256": NOTES / "rotated_pdp_corpus_20260925/evidence/raw.tar.gz",
        "gate_inputs_sha256": NOTES / "rotated_m56_export_gate_20260925/INPUTS.json",
        "gate_frozen_sha256": NOTES / "rotated_m56_export_gate_20260925/FROZEN.json",
    }
    for key, path in pins.items():
        assert sha(path) == frozen[key], key
    for arm in ("n13-m5", "n19-m6"):
        path = NOTES / f"rotated_m56_export_gate_20260925/evidence/{arm}.json"
        assert sha(path) == frozen["gate_evidence_sha256"][arm]
    for name in ("diagnostic.py", "verify.py"):
        proc = subprocess.run([sys.executable, str(HERE / name), "--self-test"],
                              capture_output=True, text=True, check=True, timeout=20)
        assert "PASS" in proc.stdout
    if args.evidence is not None:
        receipt = json.loads(args.evidence.read_text())
        assert receipt["protocol"] == frozen["domain"]
        assert receipt["freeze_sha256"] == sha(HERE / "FROZEN.json")
        assert receipt["decision"] in ("PASS", "CENSORED_OR_FAILED")
        for attempt in receipt["attempts"]:
            assert attempt["arm"] in ("n13-m5", "n19-m6")
            assert attempt["phase"] in ("producer", "verifier")
    print(json.dumps({"decision": "PASS", "freeze_sha256": sha(HERE / "FROZEN.json"),
                      "evidence_checked": args.evidence is not None}, sort_keys=True))


if __name__ == "__main__":
    main()
