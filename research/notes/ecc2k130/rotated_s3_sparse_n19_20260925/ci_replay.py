#!/usr/bin/env python3
"""Hash-only preregistration; later verify untouched n19 archive manifest."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
DEPENDENCIES = {
    "protocol_sha256": HERE / "PROTOCOL.md",
    "workflow_sha256": HERE.parents[3] / ".github/workflows/ecc2k130-s3-sparse-n19-growth.yml",
    "export_sha256": HERE / "export.py",
    "verify_sha256": HERE / "verify.py",
    "run_sha256": HERE / "run.py",
    "resource_test_sha256": HERE / "test_resource_cap.py",
    "release_gate_sha256": HERE / "release_gate.py",
    "ci_replay_sha256": HERE / "ci_replay.py",
    "dense_export_sha256": NOTES / "rotated_s3_o_branch_20260925/export.py",
    "sparse_export_sha256": NOTES / "rotated_s3_sparse_cnf_20260925/export.py",
    "sparse_verify_sha256": NOTES / "rotated_s3_sparse_cnf_20260925/verify.py",
    "sparse_final_receipt_sha256": NOTES / "rotated_s3_sparse_cnf_20260925/evidence/final/receipt.json",
    "sparse_final_summary_sha256": NOTES / "rotated_s3_sparse_cnf_20260925/evidence/final/summary.json",
    "sparse_final_n13_base_sha256": NOTES / "rotated_s3_sparse_cnf_20260925/evidence/final/sparse/producer/n13-m5/base.cnf",
    "sparse_final_n13_schema_sha256": NOTES / "rotated_s3_sparse_cnf_20260925/evidence/final/sparse/producer/n13-m5/schema.json",
    "corpus_sha256": NOTES / "rotated_pdp_corpus_20260925/evidence/raw.tar.gz",
    "corpus_frozen_sha256": NOTES / "rotated_pdp_corpus_20260925/FROZEN.json",
    "target_panel_sha256": NOTES / "rotated_s3_candidate_20260925/evidence/n19-m6/producer/result.json",
    "target_frozen_sha256": NOTES / "rotated_s3_candidate_20260925/FROZEN.json",
    "fibre_proof_sha256": NOTES / "rotated_s3_fibre_theorem_20260925/PROOF.md",
    "independent_sha256": NOTES / "rotated_pdp_corpus_20260925/verify.py",
    "parent_verify_sha256": NOTES / "rotated_subspace_support_20260925/verify.py",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_freeze():
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    for key, path in DEPENDENCIES.items():
        if sha(path) != frozen[key]:
            raise AssertionError(f"frozen input/source mismatch: {key}")
    assert frozen["domain"] == "rational-s3-o-sparse-n19-growth-v1"
    assert frozen["required_parent_pr"] == 786
    assert frozen["preregistered_parent_head"] == "1b67641f482370f6a9257b8d3427a73b836f5929"
    assert frozen["required_parent_head"] == "d16a740e8fd5fff8df6fbd4901b26bf6ff4d4df8"
    assert frozen["required_parent_merge_commit"] == "d1c320a1cb176e416f5aab22c625ce65184fea49"
    assert frozen["release_pr_number"] == 795
    assert frozen["release_label"] == "ecc2k130-n19-sparse-measure-once"
    assert frozen["caps"] == {"state_group": 100000, "transition_pairs": 500000,
                              "primary_paths": 250000, "cnf_bytes": 20000000,
                              "export_wall_seconds": 180, "verify_wall_seconds": 600,
                              "child_rss_bytes": 536870912,
                              "external_export_seconds": 195,
                              "external_verify_seconds": 630}
    return frozen


def audit_attempts(root: Path, receipt: dict, frozen: dict) -> list[bool]:
    attempts = receipt["attempts"]
    assert [row["phase"] for row in attempts] == ["export", "verify"][:len(attempts)]
    assert len(attempts) <= 2
    source = Path(receipt["source_root"])
    original_out = Path(receipt["run_out"])
    interpreter = receipt["python_executable"]
    assert source.is_absolute() and original_out.is_absolute()
    assert source.parts[-4:] == ("research", "notes", "ecc2k130",
                                 "rotated_s3_sparse_n19_20260925")
    assert Path(interpreter).is_absolute() and Path(interpreter).name in {"python", "python3", "python3.12"}
    assert receipt["python_version"].startswith("3.12.")
    passed = []
    for row in attempts:
        phase = row["phase"]
        expected_command = ([interpreter, str(source / "export.py"),
                             "--out", str(original_out / "producer")]
                            if phase == "export" else
                            [interpreter, str(source / "verify.py"), "--produced",
                             str(original_out / "producer"), "--out",
                             str(original_out / "verify.json")])
        assert row["command"] == expected_command
        assert row["rss_cap_bytes"] == frozen["caps"]["child_rss_bytes"]
        assert row["resource_stop"] in (
            None, "PROCESS_GROUP_RSS_CAP", "PEAK_RSS_CAP_AT_EXIT",
            "EXTERNAL_WALL_CAP", "PROCESS_GROUP_LEAK", "MONITOR_ERROR")
        assert row["external_timeout"] == (row["resource_stop"] == "EXTERNAL_WALL_CAP")
        assert 0 <= row["wall_seconds"]
        assert 0 <= row["sampled_group_peak_rss_bytes"]
        assert 0 <= row["direct_child_peak_rss_bytes"]
        assert 0 <= row["group_quiescence_seconds"] <= row["wall_seconds"]
        assert row["group_quiesced"] == (row["resource_stop"] != "PROCESS_GROUP_LEAK")
        assert 0 <= row["user_cpu_seconds"] and 0 <= row["system_cpu_seconds"]
        assert row["host"].startswith("Linux-")
        assert row["stdout_sha256"] == sha(root / f"{phase}.stdout.txt")
        assert row["stderr_sha256"] == sha(root / f"{phase}.stderr.txt")
        expected = root / ("producer/result.json" if phase == "export" else "verify.json")
        assert row["expected_sha256"] == (sha(expected) if expected.is_file() else None)
        wall_cap = frozen["caps"][f"external_{phase}_seconds"]
        if row["resource_stop"] == "EXTERNAL_WALL_CAP":
            assert row["wall_seconds"] >= wall_cap
        elif row["resource_stop"] in ("PROCESS_GROUP_RSS_CAP", "PEAK_RSS_CAP_AT_EXIT"):
            assert row["sampled_group_peak_rss_bytes"] >= row["rss_cap_bytes"]
        ok = (row["exit_code"] == 0 and row["resource_stop"] is None
              and row["wall_seconds"] <= wall_cap
              and row["sampled_group_peak_rss_bytes"] < row["rss_cap_bytes"]
              and row["expected_sha256"] is not None)
        passed.append(ok)
    assert all(passed[:-1]), "a failed child was followed by another phase"
    return passed


def audit_archive(root: Path, frozen):
    manifest = json.loads((root / "MANIFEST.json").read_text())
    assert manifest["freeze_sha256"] == sha(HERE / "FROZEN.json")
    actual = sorted(str(path.relative_to(root)) for path in root.rglob("*")
                    if path.is_file() and path.name != "MANIFEST.json")
    listed = [row["path"] for row in manifest["files"]]
    assert listed == actual and len(listed) == len(set(listed))
    for row in manifest["files"]:
        path = root / row["path"]
        assert path.stat().st_size == row["bytes"] and sha(path) == row["sha256"]
    receipt = json.loads((root / "receipt.json").read_text())
    assert receipt["domain"] == frozen["domain"]
    assert receipt["freeze_sha256"] == sha(HERE / "FROZEN.json")
    passed = audit_attempts(root, receipt, frozen)
    if receipt["decision"] == "PASS":
        assert passed == [True, True]
        summary = json.loads((root / "summary.json").read_text())
        produced = json.loads((root / "producer/result.json").read_text())
        verified = json.loads((root / "verify.json").read_text())
        assert summary["decision"] == verified["decision"] == "PASS"
        assert summary["artifact_sha256"] == {"producer": sha(root / "producer/result.json"),
                                               "verifier": sha(root / "verify.json")}
        expected_summary = {"decision": "PASS", "domain": frozen["domain"],
                            "growth": produced["growth"],
                            "transition_cases": produced["transition_cases"],
                            "primary_paths": produced["primary_paths"],
                            "signed_point_tuples": verified["signed_point_tuples"],
                            "target_labels": len(verified["target_rows"]),
                            "variables": produced["variables"],
                            "clauses": produced["clauses"],
                            "bytes": produced["bytes"],
                            "cold_children": {
                                "export": {key: produced[key] for key in
                                           ("wall_seconds", "cpu_seconds", "peak_rss_bytes")},
                                "verify": {key: verified[key] for key in
                                           ("wall_seconds", "cpu_seconds", "peak_rss_bytes")}},
                            "artifact_sha256": {
                                "producer": sha(root / "producer/result.json"),
                                "verifier": sha(root / "verify.json")}}
        assert summary == expected_summary
        assert produced["growth"] == json.loads(
            (root / "producer/growth.json").read_text())["completed"]
        assert produced["transition_cases"] == verified["transition_cases"]
        for phase, raw in zip(("export", "verify"), (produced, verified)):
            attempt = next(row for row in receipt["attempts"] if row["phase"] == phase)
            assert 0 <= raw["wall_seconds"] <= attempt["wall_seconds"] + 1
            assert 0 <= raw["cpu_seconds"] <= (
                attempt["user_cpu_seconds"] + attempt["system_cpu_seconds"] + 1)
            assert 0 <= raw["peak_rss_bytes"] <= attempt["direct_child_peak_rss_bytes"]
        assert len(verified["target_rows"]) == 33
        with tempfile.TemporaryDirectory(prefix="n19-sparse-independent-replay-") as scratch:
            replay_path = Path(scratch) / "verify.json"
            replay = subprocess.run(
                [sys.executable, str(HERE / "verify.py"), "--produced",
                 str((root / "producer").resolve()), "--out", str(replay_path)],
                cwd=HERE.parents[3], capture_output=True, text=True,
                check=False, timeout=frozen["caps"]["external_verify_seconds"])
            if replay.returncode != 0 or not replay_path.is_file():
                raise AssertionError(f"fresh independent archive replay failed: "
                                     f"exit={replay.returncode}, stderr={replay.stderr[-1000:]}")
            fresh = json.loads(replay_path.read_text())
        resource_keys = {"wall_seconds", "cpu_seconds", "peak_rss_bytes"}
        old_semantics = {key: value for key, value in verified.items()
                         if key not in resource_keys}
        new_semantics = {key: value for key, value in fresh.items()
                         if key not in resource_keys}
        assert old_semantics == new_semantics and fresh["decision"] == "PASS"
    else:
        assert receipt["decision"] in {"CENSORED", "FAILED"}
        assert passed != [True, True]
        assert "summary.json" not in actual
        if any(row["resource_stop"] is not None for row in receipt["attempts"]):
            assert receipt["decision"] == "CENSORED"
    return {"decision": receipt["decision"], "files": len(listed),
            "freeze_sha256": receipt["freeze_sha256"]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    frozen = check_freeze()
    if args.evidence is None:
        print(json.dumps({"decision": "FROZEN_HASH_ONLY",
                          "freeze_sha256": sha(HERE / "FROZEN.json")}, sort_keys=True))
    else:
        print(json.dumps(audit_archive(args.evidence, frozen), sort_keys=True))


if __name__ == "__main__":
    main()
