#!/usr/bin/env python3
"""Fail-closed source freeze and fresh independent replay of every archived arm."""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DOMAIN = "ecc2k130-native-m3-chain-20260929-v1"
FROZEN_PATHS = {
    ".github/workflows/ecc2k130-native-m3-chain.yml",
    *(
        f"research/notes/ecc2k130/native_m3_chain_20260929/{name}"
        for name in ("PROTOCOL.md", "chain.py", "reference.py", "panel.py",
                     "verify_panel.py", "run.py", "ci_replay.py")
    ),
    "research/notes/ecc2k130/native_fullpoint_edge_20260929/native_relation.py",
    "research/notes/ecc2k130/native_fullpoint_edge_20260929/produce.py",
    "research/notes/ecc2k130/native_fullpoint_edge_20260929/verify.py",
    "research/notes/ecc2k130/symbolic_dag_fullpoint_20260925/dag.py",
    "research/notes/ecc2k130/m10_export_capacity_20260925/capacity.py",
    "research/notes/ecc2k130/rotated_subspace_support_20260925/gate.py",
    "research/notes/ecc2k130/symbolic_dag_dimacs_gate_20260925/export.py",
    "research/notes/ecc2k130/symbolic_dag_dimacs_gate_20260925/verify.py",
    "research/ecc2k130_factor_base_replication_20260925/exact_smoke.json",
    "research/ecc2k130_factor_base_replication_20260925/exact_replay.json",
    "research/ecc2k130_relations/relations.py",
    "research/ecc2k130_relations/fastfield.py",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_freeze() -> dict:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert frozen["schema"] == "ecc2k130-native-m3-chain-freeze-v1"
    assert frozen["source_main_parent"] == "a218e41182cda97f4be6f636333202882142eec9"
    assert frozen["toy"] == {"n": 5, "modulus": "0x25", "a": 0, "b": 1,
                             "bases": [[1], [2], [4]]}
    assert frozen["leaf"] == {"n": 131, "modulus": "0x800000000000000000000000000002007",
                              "lines": [[1, 0], [1, 4]], "beta": 3,
                              "dimensions": [44, 44, 43]}
    assert frozen["toy_child_wall_seconds"] == 300
    assert frozen["toy_external_wall_seconds"] == 330
    assert frozen["leaf_child_wall_seconds"] == 600
    assert frozen["leaf_external_wall_seconds"] == 630
    assert frozen["toy_rss_cap_bytes"] == 1 << 30
    assert frozen["leaf_rss_cap_bytes"] == 2 << 30
    assert frozen["dag_node_cap"] == 2_000_000
    assert frozen["cnf_byte_cap"] == 256 << 20
    assert set(frozen["sha256"]) == FROZEN_PATHS
    commit = frozen["source_commit"]
    assert len(commit) == 40 and all(char in "0123456789abcdef" for char in commit)
    subprocess.run(["git", "merge-base", "--is-ancestor", commit, "HEAD"],
                   cwd=ROOT, check=True, capture_output=True)
    for relative, expected in sorted(frozen["sha256"].items()):
        path = (ROOT / relative).resolve()
        assert path.is_relative_to(ROOT) and path.is_file(), relative
        assert sha(path) == expected, relative
        committed = subprocess.run(["git", "show", f"{commit}:{relative}"],
                                   cwd=ROOT, check=True, capture_output=True).stdout
        assert hashlib.sha256(committed).hexdigest() == expected, relative
        if path.suffix == ".py":
            ast.parse(path.read_text(), filename=str(path))
    return frozen


def replay_archive(directory: Path, frozen: dict) -> dict:
    assert directory.is_dir()
    receipt = json.loads((directory / "receipt.json").read_text())
    assert receipt["domain"] == DOMAIN
    assert receipt["decision"] in ("PASS", "CAPACITY_CENSORED")
    assert receipt["freeze_sha256"] == sha(HERE / "FROZEN.json")
    assert sha(directory / "host.json") == receipt["host_sha256"]
    assert len(receipt["attempts"]) == 6
    outcomes = []
    for arm in ("toy", "leaf10", "leaf14"):
        arm_dir = directory / arm
        producer = arm_dir / "producer"
        archived = json.loads((arm_dir / "verify.json").read_text())
        result = json.loads((producer / "result.json").read_text())
        assert result["domain"] == archived["domain"] == DOMAIN
        assert result["arm"] == archived["arm"] == arm
        assert result["decision"] == archived["decision"]
        assert sha(producer / "result.json") == receipt["results"][arm]["producer"]
        assert sha(arm_dir / "verify.json") == receipt["results"][arm]["verifier"]
        external_cap = frozen["toy_external_wall_seconds"] if arm == "toy" else frozen["leaf_external_wall_seconds"]
        child_cap = frozen["toy_child_wall_seconds"] if arm == "toy" else frozen["leaf_child_wall_seconds"]
        rss_cap = frozen["toy_rss_cap_bytes"] if arm == "toy" else frozen["leaf_rss_cap_bytes"]
        for phase in ("producer", "verifier"):
            attempt = next(item for item in receipt["attempts"]
                           if item["arm"] == arm and item["phase"] == phase)
            assert attempt["exit_code"] == 0 and attempt["external_timeout"] is False
            assert attempt["result_decision"] == result["decision"]
            assert attempt["external_wall_seconds"] <= external_cap
            assert attempt["child_wall_seconds"] <= child_cap
            assert attempt["child_peak_rss_bytes"] <= rss_cap
            assert sha(directory / f"{arm}.{phase}.stdout.txt") == attempt["stdout_sha256"]
            assert sha(directory / f"{arm}.{phase}.stderr.txt") == attempt["stderr_sha256"]
            assert receipt["results"][arm][phase] == attempt["result_sha256"]
        with tempfile.TemporaryDirectory(prefix=f"native-m3-{arm}-") as temp:
            fresh_path = Path(temp) / "verify.json"
            child = subprocess.run([sys.executable, str(HERE / "verify_panel.py"),
                                    "--arm", arm, "--producer", str(producer),
                                    "--out", str(fresh_path)],
                                   cwd=HERE, capture_output=True, text=True, timeout=external_cap)
            assert child.returncode == 0, (arm, child.stderr)
            fresh = json.loads(fresh_path.read_text())
        stable = ("domain", "arm", "decision", "producer_sha256")
        for key in stable:
            assert fresh[key] == archived[key], (arm, key)
        if arm == "toy":
            fields = ("checks", "factor_triples", "supported_targets")
        elif result["censor"] == "DAG_NODE_CAP":
            fields = ("censor", "partial_counts")
        else:
            fields = ("line", "dag", "cnf_target_O")
        for key in fields:
            assert fresh[key] == archived[key], (arm, key)
        outcomes.append({"arm": arm, "decision": result["decision"],
                         "producer_sha256": receipt["results"][arm]["producer"]})
    assert receipt["decision"] == ("PASS" if all(row["decision"] == "PASS" for row in outcomes)
                                    else "CAPACITY_CENSORED")
    return {"status": "ARCHIVE_REPLAY_PASS", "outcomes": outcomes}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    frozen = check_freeze()
    result = (replay_archive(args.evidence.resolve(), frozen) if args.evidence else
              {"status": "SOURCE_FREEZE_PASS", "source_commit": frozen["source_commit"],
               "files": len(frozen["sha256"])})
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
