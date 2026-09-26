#!/usr/bin/env python3
"""Hash-only preregistration and archive-only sparse-CNF evidence audit."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
DENSE = HERE.parent / "rotated_s3_o_branch_20260925"
PANELS = [["n2-m4", 2, 4, 7], ["n3-m5", 3, 5, 11],
          ["n4-m4", 4, 4, 19], ["n13-m5", 13, 5, 8219]]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_freeze(frozen):
    files = {"protocol_sha256": HERE / "PROTOCOL.md",
             "export_sha256": HERE / "export.py",
             "verify_sha256": HERE / "verify.py",
             "runner_sha256": HERE / "run.py",
             "ci_replay_sha256": HERE / "ci_replay.py",
             "dense_export_sha256": DENSE / "export.py",
             "dense_verify_sha256": DENSE / "verify.py",
             "dense_frozen_sha256": DENSE / "FROZEN.json",
             "dense_receipt_sha256": DENSE / "evidence/receipt.json",
             "dense_result_sha256": DENSE / "evidence/producer/result.json"}
    for key, path in files.items():
        assert sha(path) == frozen[key], key
    assert sha(HERE.parent.parent.parent.parent / ".github/workflows/ecc2k130-s3-sparse-cnf.yml") == frozen["workflow_sha256"]
    for relative, digest in frozen["dense_dependency_sha256"].items():
        assert sha(HERE.parent / relative) == digest, relative
    for name, *_ in PANELS:
        for file in ("base.cnf", "schema.json", "paths.jsonl.gz"):
            path = DENSE / "evidence/producer" / name / file
            item = frozen["panels"][name][file]
            assert path.stat().st_size == item["bytes"] and sha(path) == item["sha256"]
    assert frozen["domain"] == "rational-s3-o-sparse-cnf-n2n3n4-n13-v2"
    assert frozen["panels_spec"] == PANELS
    assert frozen["dense_pr"] == 781
    assert frozen["source_main_commit"] == "8c5178b00b2af91548a5b4f88558723b3ad656c9"
    assert frozen["child_wall_cap_seconds"] == 180
    assert frozen["child_rss_cap_bytes"] == 512 * 1024 * 1024
    assert frozen["external_child_cap_seconds"] == 195


def check_first_failure(frozen):
    root = HERE / "evidence/first_failure"
    if not root.exists():
        return None
    receipt = json.loads((root / "receipt.json").read_text())
    assert receipt["decision"] == "CENSORED_OR_FAILED"
    assert receipt["freeze_sha256"] == frozen["first_failure_freeze_sha256"]
    assert [row["phase"] for row in receipt["attempts"]] == [
        "dense_export", "dense_replay", "sparse_export"]
    assert [row["exit_code"] for row in receipt["attempts"]] == [0, 0, 1]
    assert json.loads((root / "dense/verify.json").read_text())["decision"] == "PASS"
    assert json.loads((root / "sparse/producer/failure.json").read_text())["error"] == "AssertionError()"
    manifest = json.loads((root / "MANIFEST.json").read_text())
    assert manifest["freeze_sha256"] == receipt["freeze_sha256"]
    expected_files = sorted(str(path.relative_to(root)) for path in root.rglob("*")
                            if path.is_file() and path.name != "MANIFEST.json")
    assert [row["path"] for row in manifest["files"]] == expected_files
    for row in manifest["files"]:
        path = root / row["path"]
        assert path.stat().st_size == row["bytes"] and sha(path) == row["sha256"]
    return {"phases": 3, "failed_phase": "sparse_export"}


def check_evidence(receipt_path, frozen, expected_freeze_sha):
    receipt = json.loads(receipt_path.read_text())
    root = receipt_path.parent
    assert receipt["protocol"] == frozen["domain"]
    assert receipt["freeze_sha256"] == expected_freeze_sha
    assert receipt["decision"] == "PASS"
    phases = ["dense_export", "dense_replay", "sparse_export", "sparse_replay"]
    assert [row["phase"] for row in receipt["attempts"]] == phases
    outputs = [root / "dense/producer/result.json", root / "dense/verify.json",
               root / "sparse/producer/result.json", root / "sparse/verify.json"]
    for row, path in zip(receipt["attempts"], outputs):
        assert row["exit_code"] == 0 and not row["external_timeout"]
        assert row["wall_seconds"] <= 195
        assert sha(path) == row["expected_sha256"]
        for stream in ("stdout", "stderr"):
            assert sha(root / f"{row['phase']}.{stream}.txt") == row[f"{stream}_sha256"]
    dense_export, dense_verify, sparse_export, sparse_verify = [json.loads(p.read_text()) for p in outputs]
    assert dense_verify["decision"] == sparse_verify["decision"] == "PASS"
    assert sparse_export["domain"] == sparse_verify["domain"] == frozen["domain"]
    assert sparse_verify["producer_sha256"] == sha(outputs[2])
    for measured in (dense_export, dense_verify, sparse_export, sparse_verify):
        assert measured["wall_seconds"] <= 180
        assert measured["peak_rss_bytes"] <= 512 * 1024 * 1024
    assert set(sparse_verify["negative_controls"]) == {
        "forbid_zero_zero_O", "replace_O_output",
        "drop_sequential_clause", "wrong_O_target", "n13_empty_pair_clause"}
    assert [row["panel"] for row in sparse_export["panels"]] == [row[0] for row in PANELS]
    assert [row["panel"] for row in sparse_verify["panels"]] == [row[0] for row in PANELS]
    for row in dense_export["panels"]:
        name = row["panel"]
        for file, key in (("base.cnf", "base_sha256"),
                          ("schema.json", "schema_sha256"),
                          ("paths.jsonl.gz", "paths_sha256")):
            assert sha(root / "dense/producer" / name / file) == row[key] == frozen["panels"][name][file]["sha256"]
    for row in sparse_export["panels"]:
        name = row["panel"]
        for file, key in (("base.cnf", "base_sha256"),
                          ("schema.json", "schema_sha256"),
                          ("paths.jsonl.gz", "paths_sha256")):
            assert sha(root / "sparse/producer" / name / file) == row[key]
            if file == "paths.jsonl.gz":
                assert row[key] == frozen["panels"][name][file]["sha256"]
    control = sparse_verify["negative_controls"]["n13_empty_pair_clause"]
    assert control["rejection"].startswith("point-law output mismatch")
    assert sha(root / "sparse/producer/negative_empty_pair.cnf") == control["mutated_cnf_sha256"]
    summary = json.loads((root / "summary.json").read_text())
    assert summary["decision"] == "PASS" and len(summary["panels"]) == 4
    assert sum(row["candidate_paths"] for row in summary["panels"]) == 24955
    assert sum(row["target_labels"] for row in summary["panels"]) == 61
    assert [summary["child_sha256"][phase] for phase in phases] == [sha(path) for path in outputs]
    manifest = json.loads((root / "MANIFEST.json").read_text())
    assert manifest["freeze_sha256"] == expected_freeze_sha
    expected_files = sorted(str(path.relative_to(root)) for path in root.rglob("*")
                            if path.is_file() and path.name != "MANIFEST.json")
    assert [row["path"] for row in manifest["files"]] == expected_files
    for row in manifest["files"]:
        path = root / row["path"]
        assert path.stat().st_size == row["bytes"] and sha(path) == row["sha256"]
    return {"model_paths": 24955, "targets": 61,
            "n13_variables": summary["panels"][3]["sparse"]["variables"],
            "n13_clauses": summary["panels"][3]["sparse"]["clauses"]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    check_freeze(frozen)
    first_failure = check_first_failure(frozen)
    previous = HERE / "evidence/second_attempt/receipt.json"
    second_attempt = (check_evidence(previous, frozen, frozen["second_attempt_freeze_sha256"])
                      if previous.exists() else None)
    evidence = (check_evidence(args.evidence, frozen, sha(HERE / "FROZEN.json"))
                if args.evidence else None)
    print(json.dumps({"decision": "PASS", "freeze_sha256": sha(HERE / "FROZEN.json"),
                      "first_failure": first_failure, "second_attempt": second_attempt,
                      "evidence": evidence}, sort_keys=True))


if __name__ == "__main__":
    main()
