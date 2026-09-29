#!/usr/bin/env python3
"""Fail-closed source/input freeze and optional fresh raw-row replay."""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
FROZEN_PATHS = {
    ".github/workflows/ecc2k130-native-fullpoint-edge.yml",
    "research/notes/ecc2k130/native_fullpoint_edge_20260929/PROTOCOL.md",
    *(
        f"research/notes/ecc2k130/native_fullpoint_edge_20260929/{name}.py"
        for name in ("native_relation", "produce", "verify", "run", "ci_replay")
    ),
    "research/notes/ecc2k130/symbolic_dag_fullpoint_20260925/dag.py",
    "research/ecc2k130_factor_base_replication_20260925/exact_smoke.json",
    "research/ecc2k130_factor_base_replication_20260925/exact_replay.json",
    "research/ecc2k130_relations/relations.py",
    "research/ecc2k130_relations/fastfield.py",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_freeze() -> dict:
    frozen_path = HERE / "FROZEN.json"
    frozen = json.loads(frozen_path.read_text())
    assert frozen["schema"] == "ecc2k130-native-fullpoint-edge-freeze-v1"
    assert frozen["source_main_parent"] == "c237ef9bbfcabb286fd37f1448ebe76bd79a4f99"
    assert frozen["fields"] == [[2, 7], [3, 11]]
    assert frozen["toy_a"] == [0, 1]
    assert frozen["toy_b"] == "all_nonzero"
    assert frozen["leaf_lines"] == [[1, 0], [1, 4]]
    assert frozen["exact_modulus_hex"] == "0x800000000000000000000000000002007"
    assert frozen["child_wall_cap_seconds"] == 600
    assert frozen["external_child_cap_seconds"] == 630
    assert frozen["child_rss_cap_bytes"] == 1 << 30
    assert set(frozen["sha256"]) == FROZEN_PATHS
    commit = frozen["source_commit"]
    assert len(commit) == 40 and all(char in "0123456789abcdef" for char in commit)
    subprocess.run(["git", "merge-base", "--is-ancestor", commit, "HEAD"],
                   cwd=ROOT, check=True, capture_output=True)
    for relative, expected in sorted(frozen["sha256"].items()):
        path = (ROOT / relative).resolve()
        assert path.is_relative_to(ROOT) and path.is_file(), relative
        assert sha(path) == expected, relative
        committed = subprocess.run(["git", "show", f"{commit}:{relative}"], cwd=ROOT,
                                   check=True, capture_output=True).stdout
        assert hashlib.sha256(committed).hexdigest() == expected, relative
        if path.suffix == ".py":
            ast.parse(path.read_text(), filename=str(path))
    return frozen


def replay_archive(directory: Path, frozen: dict) -> None:
    from verify import replay

    assert directory.is_dir()
    archived = json.loads((directory / "verify.json").read_text())
    receipt = json.loads((directory / "receipt.json").read_text())
    assert receipt["decision"] == "PASS"
    assert receipt["freeze_sha256"] == sha(HERE / "FROZEN.json")
    assert len(receipt["attempts"]) == 2
    producer = directory / "producer"
    assert sha(directory / "host.json") == receipt["host_sha256"]
    for phase, attempt in zip(("producer", "verifier"), receipt["attempts"], strict=True):
        assert attempt["phase"] == phase
        assert attempt["exit_code"] == 0 and attempt["external_timeout"] is False
        assert attempt["child_wall_seconds"] <= frozen["child_wall_cap_seconds"]
        assert attempt["child_peak_rss_bytes"] <= frozen["child_rss_cap_bytes"]
        assert attempt["wall_seconds"] <= frozen["external_child_cap_seconds"]
        assert sha(directory / f"{phase}.stdout.txt") == attempt["stdout_sha256"]
        assert sha(directory / f"{phase}.stderr.txt") == attempt["stderr_sha256"]
    assert sha(producer / "result.json") == receipt["producer_sha256"]
    assert sha(directory / "verify.json") == receipt["verifier_sha256"]
    assert receipt["attempts"][0]["result_sha256"] == receipt["producer_sha256"]
    assert receipt["attempts"][1]["result_sha256"] == receipt["verifier_sha256"]
    with tempfile.TemporaryDirectory(prefix="native-edge-replay-") as temp:
        fresh = replay(producer, Path(temp) / "verify.json")
    for key in ("domain", "decision", "toy", "leaves", "raw_rows",
                "rows_sha256", "producer_sha256"):
        assert fresh[key] == archived[key], key
    assert fresh["decision"] == "PASS"
    print(json.dumps({"status": "ARCHIVE_REPLAY_PASS", "rows": fresh["raw_rows"],
                      "rows_sha256": fresh["rows_sha256"]}, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    frozen = check_freeze()
    if args.evidence:
        replay_archive(args.evidence, frozen)
    else:
        print(json.dumps({"status": "SOURCE_FREEZE_PASS",
                          "source_commit": frozen["source_commit"],
                          "files": len(frozen["sha256"])}, sort_keys=True))


if __name__ == "__main__":
    main()
