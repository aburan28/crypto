#!/usr/bin/env python3
"""Hash-only frozen-source preflight and optional committed evidence audit."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "rotated_subspace_support_20260925/gate.py"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_evidence(path: Path, frozen: dict):
    receipt = json.loads(path.read_text())
    assert receipt["protocol"] == frozen["domain"]
    assert receipt["freeze_sha256"] == sha(HERE / "FROZEN.json")
    assert receipt["decision"] == "PASS"
    assert [item["phase"] for item in receipt["attempts"]] == ["producer", "verifier"]
    out = path.parent
    for item in receipt["attempts"]:
        assert item["exit_code"] == 0 and not item["external_timeout"]
        phase = item["phase"]
        expected = out / ("producer/result.json" if phase == "producer" else "verify.json")
        assert sha(expected) == item["expected_sha256"]
        assert sha(out / f"{phase}.stdout.txt") == item["stdout_sha256"]
        assert sha(out / f"{phase}.stderr.txt") == item["stderr_sha256"]
        assert item["wall_seconds"] <= 195
    producer = json.loads((out / "producer/result.json").read_text())
    verifier = json.loads((out / "verify.json").read_text())
    assert producer["domain"] == frozen["domain"]
    assert verifier["domain"] == frozen["domain"] and verifier["decision"] == "PASS"
    assert verifier["producer_sha256"] == sha(out / "producer/result.json")
    assert verifier["pair_rows_checked"] == sum(row["ordered_rational_pairs"]
                                                 for row in producer["fields"])
    assert verifier["chain_rows_checked"] == sum(row["distinct_affine_candidate_paths"]
                                                  for row in producer["chain_panels"])
    for name, key in (("pair_rows.jsonl.gz", "pair_rows_sha256"),
                      ("chain_rows.jsonl.gz", "chain_rows_sha256")):
        assert sha(out / "producer" / name) == producer[key]
    for measured in (producer, verifier):
        assert measured["wall_seconds"] <= 180
        assert measured["peak_rss_bytes"] <= 512 * 1024 * 1024
    assert producer["first_nonlift_formal_root"] is not None
    assert producer["first_sign_restriction_loss"] is not None
    assert all(row["candidate_only_paths"] == 0 and
               row["candidate_paths_missing_target_sign"] == 0
               for row in producer["chain_panels"])
    return {"pair_rows": verifier["pair_rows_checked"],
            "chain_rows": verifier["chain_rows_checked"]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    sources = {
        "protocol_sha256": HERE / "PROTOCOL.md",
        "proof_sha256": HERE / "PROOF.md",
        "producer_sha256": HERE / "theorem_check.py",
        "verify_sha256": HERE / "verify.py",
        "runner_sha256": HERE / "run.py",
        "ci_replay_sha256": HERE / "ci_replay.py",
        "parent_sha256": PARENT,
    }
    for key, path in sources.items():
        assert sha(path) == frozen[key], key
    assert frozen["domain"] == "finite-fibre-s3-n1to7-chain-n3m4-n4m4-n5m3-v1"
    assert frozen["fields"] == {str(n): poly for n, poly in
                                 ((1, 0x3), (2, 0x7), (3, 0xb), (4, 0x13),
                                  (5, 0x25), (6, 0x43), (7, 0x83))}
    assert frozen["chain_panels"] == [[3, 4], [4, 4], [5, 3]]
    assert frozen["child_wall_cap_seconds"] == 180
    assert frozen["child_rss_cap_bytes"] == 512 * 1024 * 1024
    evidence = check_evidence(args.evidence, frozen) if args.evidence else None
    print(json.dumps({"decision": "PASS", "freeze_sha256": sha(HERE / "FROZEN.json"),
                      "evidence": evidence}, sort_keys=True))


if __name__ == "__main__":
    main()
