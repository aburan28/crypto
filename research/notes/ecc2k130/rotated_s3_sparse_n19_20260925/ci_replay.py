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
    assert frozen["caps"] == {"state_group": 100000, "transition_pairs": 500000,
                              "primary_paths": 250000, "cnf_bytes": 20000000,
                              "export_wall_seconds": 180, "verify_wall_seconds": 600,
                              "child_rss_bytes": 536870912,
                              "external_export_seconds": 195,
                              "external_verify_seconds": 630}
    return frozen


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
    if receipt["decision"] == "PASS":
        assert len(receipt["attempts"]) == 2
        assert all(row["exit_code"] == 0 and not row["external_timeout"]
                   for row in receipt["attempts"])
        summary = json.loads((root / "summary.json").read_text())
        produced = json.loads((root / "producer/result.json").read_text())
        verified = json.loads((root / "verify.json").read_text())
        assert summary["decision"] == verified["decision"] == "PASS"
        assert summary["artifact_sha256"] == {"producer": sha(root / "producer/result.json"),
                                               "verifier": sha(root / "verify.json")}
        assert (summary["primary_paths"], summary["variables"], summary["clauses"],
                summary["bytes"]) == (produced["primary_paths"], produced["variables"],
                                      produced["clauses"], produced["bytes"])
        assert (summary["signed_point_tuples"], summary["target_labels"]) == (
            verified["signed_point_tuples"], len(verified["target_rows"]))
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
        assert len(receipt["attempts"]) <= 2
        assert "summary.json" not in actual
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
