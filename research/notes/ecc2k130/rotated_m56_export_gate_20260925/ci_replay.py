#!/usr/bin/env python3
"""Hash-pinned preflight and exact deterministic replay of committed gate data."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
CORPUS = HERE.parent / "rotated_pdp_corpus_20260925"
PARENT = HERE.parent / "rotated_subspace_support_20260925"
ARCHIVE = CORPUS / "evidence/raw.tar.gz"
ARMS = ("n13-m5", "n19-m6")
NONDETERMINISTIC = {"wall_seconds", "cpu_seconds", "peak_rss_bytes"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def preflight() -> dict:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    paths = {
        "inputs_sha256": HERE / "INPUTS.json",
        "protocol_sha256": HERE / "PROTOCOL.md",
        "gate_sha256": HERE / "gate.py",
        "runner_sha256": HERE / "run.py",
        "replay_sha256": Path(__file__),
        "archive_sha256": ARCHIVE,
        "corpus_frozen_sha256": CORPUS / "FROZEN.json",
        "corpus_verify_sha256": CORPUS / "verify.py",
        "parent_verify_sha256": PARENT / "verify.py",
    }
    for key, path in paths.items():
        assert sha(path) == frozen[key], (key, path)
    manifest = json.loads((HERE / "INPUTS.json").read_text())
    assert manifest["domain"] == frozen["domain"]
    assert manifest["corpus_merged_commit"] == frozen["corpus_commit"]
    assert manifest["corpus_archive_sha256"] == frozen["archive_sha256"]
    with tarfile.open(ARCHIVE, "r:gz") as tar:
        for arm in ARMS:
            entry = manifest["arms"][arm]
            rows = entry["targets"]
            assert len(rows) == 8
            assert [r["class"] for r in rows] == ["planted"] * 4 + ["negative"] * 4
            for name, digest in (("targets.json", entry["target_file_sha256"]),
                                 ("factors.json", entry["factor_file_sha256"])):
                member = tar.getmember(f"raw/{arm}/{name}")
                raw = tar.extractfile(member).read()
                assert hashlib.sha256(raw).hexdigest() == digest
                if name == "targets.json":
                    source_rows = json.loads(raw)
                    assert rows == [{k: row[k] for k in ("class", "Q", "R", "coset_multiplicities")}
                                    for row in source_rows]
    proc = subprocess.run([sys.executable, str(HERE / "gate.py"), "--self-test"],
                          capture_output=True, text=True, check=False)
    assert proc.returncode == 0 and proc.stdout.strip() == "self-test PASS", proc.stderr
    return frozen


def replay(evidence: Path, frozen: dict) -> None:
    receipt = json.loads((evidence / "receipt.json").read_text())
    assert receipt["protocol"] == frozen["domain"]
    for arm in ARMS:
        row = receipt["arms"][arm]
        result_path = evidence / f"{arm}.json"
        stdout_path = evidence / f"{arm}.stdout.txt"
        stderr_path = evidence / f"{arm}.stderr.txt"
        assert row["exit_code"] == 0
        assert row["result_sha256"] == sha(result_path)
        assert row["stdout_sha256"] == sha(stdout_path)
        assert row["stderr_sha256"] == sha(stderr_path)
        assert not stderr_path.read_text()
        result = json.loads(result_path.read_text())
        assert result["arm"] == arm and result["archive_sha256"] == frozen["archive_sha256"]
        assert result["source_sha256"] == frozen["gate_sha256"]
        assert result["wall_seconds"] <= 600 and result["peak_rss_bytes"] <= 512 * 1024 * 1024
        with tempfile.TemporaryDirectory() as temp:
            fresh_path = Path(temp) / f"{arm}.json"
            proc = subprocess.run([sys.executable, str(HERE / "gate.py"), "--arm", arm,
                                   "--out", str(fresh_path)],
                                  capture_output=True, text=True, check=False)
            assert proc.returncode == 0, proc.stderr
            fresh = json.loads(fresh_path.read_text())
            for key in NONDETERMINISTIC:
                result.pop(key)
                fresh.pop(key)
            assert fresh == result, arm


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    frozen = preflight()
    if args.evidence:
        replay(args.evidence, frozen)
    print("rotated m5/m6 exporter semantic gate replay PASS")


if __name__ == "__main__":
    main()
