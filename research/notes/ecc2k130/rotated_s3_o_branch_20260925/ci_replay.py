#!/usr/bin/env python3
"""Hash-only preflight and archive-only receipt audit for branch CNF gate."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evidence(path: Path, frozen: dict):
    receipt = json.loads(path.read_text())
    assert receipt["protocol"] == frozen["domain"]
    assert receipt["freeze_sha256"] == sha(HERE / "FROZEN.json")
    assert receipt["decision"] == "PASS"
    assert [row["phase"] for row in receipt["attempts"]] == ["producer", "verifier"]
    root = path.parent
    for row in receipt["attempts"]:
        assert row["exit_code"] == 0 and not row["external_timeout"]
        assert row["wall_seconds"] <= 195
        output = root / ("producer/result.json" if row["phase"] == "producer" else "verify.json")
        assert sha(output) == row["expected_sha256"]
        for stream in ("stdout", "stderr"):
            assert sha(root / f"{row['phase']}.{stream}.txt") == row[f"{stream}_sha256"]
    producer = json.loads((root / "producer/result.json").read_text())
    verifier = json.loads((root / "verify.json").read_text())
    assert producer["domain"] == verifier["domain"] == frozen["domain"]
    assert verifier["decision"] == "PASS"
    assert verifier["producer_sha256"] == sha(root / "producer/result.json")
    assert [row["panel"] for row in producer["panels"]] == ["n2-m4", "n3-m5", "n4-m4", "n13-m5"]
    assert [row["panel"] for row in verifier["panels"]] == ["n2-m4", "n3-m5", "n4-m4", "n13-m5"]
    for row in producer["panels"]:
        panel = root / "producer" / row["panel"]
        for name, key in (("base.cnf", "base_sha256"),
                          ("schema.json", "schema_sha256"),
                          ("paths.jsonl.gz", "paths_sha256")):
            assert sha(panel / name) == row[key]
    for measured in (producer, verifier):
        assert measured["wall_seconds"] <= 180
        assert measured["peak_rss_bytes"] <= 512 * 1024 * 1024
    assert set(verifier["negative_controls"]) == {
        "nonlift", "sign_incomplete", "mutated_O_clause", "wrong_O_target_literal"}
    manifest = json.loads((root / "MANIFEST.json").read_text())
    assert manifest["freeze_sha256"] == sha(HERE / "FROZEN.json")
    for entry in manifest["files"]:
        file = root / entry["path"]
        assert file.stat().st_size == entry["bytes"] and sha(file) == entry["sha256"]
    return {"model_paths": sum(row["candidate_paths"] for row in verifier["panels"]),
            "targets": sum(row["target_labels"] for row in verifier["panels"])}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    pins = {
        "protocol_sha256": HERE / "PROTOCOL.md",
        "proof_sha256": HERE / "PROOF.md",
        "export_sha256": HERE / "export.py",
        "verify_sha256": HERE / "verify.py",
        "runner_sha256": HERE / "run.py",
        "ci_replay_sha256": HERE / "ci_replay.py",
        "parent_sha256": NOTES / "rotated_subspace_support_20260925/gate.py",
        "corpus_sha256": NOTES / "rotated_pdp_corpus_20260925/evidence/raw.tar.gz",
        "prior_sha256": NOTES / "rotated_s3_candidate_20260925/evidence/n13-m5/producer/result.json",
        "gate_inputs_sha256": NOTES / "rotated_m56_export_gate_20260925/INPUTS.json",
        "gate_frozen_sha256": NOTES / "rotated_m56_export_gate_20260925/FROZEN.json",
        "gate_n13_sha256": NOTES / "rotated_m56_export_gate_20260925/evidence/n13-m5.json",
    }
    for key, path in pins.items():
        assert sha(path) == frozen[key], key
    assert frozen["domain"] == "rational-s3-o-branch-cnf-n2n3n4-n13-v1"
    assert frozen["panels"] == [["n2-m4", 2, 4, 0x7], ["n3-m5", 3, 5, 0xb],
                                 ["n4-m4", 4, 4, 0x13], ["n13-m5", 13, 5, 0x201b]]
    assert frozen["child_wall_cap_seconds"] == 180
    assert frozen["child_rss_cap_bytes"] == 512 * 1024 * 1024
    checked = evidence(args.evidence, frozen) if args.evidence else None
    print(json.dumps({"decision": "PASS", "freeze_sha256": sha(HERE / "FROZEN.json"),
                      "evidence": checked}, sort_keys=True))


if __name__ == "__main__":
    main()
