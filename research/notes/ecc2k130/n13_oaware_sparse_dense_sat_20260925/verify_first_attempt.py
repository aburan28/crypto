#!/usr/bin/env python3
"""Independently replay the first local preflight and six tiny SAT smokes."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
import subprocess

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
ARCHIVE = HERE / "evidence/first_local_attempt_20260929"
SOURCE_HEAD = "3fa70fc1b30425a78654c20f46342faca7f01772"
FREEZE_SHA256 = "ada0ecb94236cd3430653cb6203eb80782a731b4b508f5d427d879a43b00f717"
SOLVERS = ("cryptominisat5", "kissat", "cadical")
SAT_CNF = b"p cnf 1 1\n1 0\n"
UNSAT_CNF = b"p cnf 1 2\n1 0\n-1 0\n"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git(*args: str) -> bytes:
    return subprocess.check_output(["git", *args], cwd=REPO)


def main() -> None:
    manifest = json.loads((ARCHIVE / "MANIFEST.json").read_text())
    assert manifest["schema"] == "n13-oaware-paired-first-local-attempt-v1"
    assert manifest["source_head"] == SOURCE_HEAD
    assert manifest["freeze_sha256"] == FREEZE_SHA256
    assert manifest["classification"] == "PREFLIGHT_AND_SMOKE_COMPLETE_PANEL_HELD"
    assert sha((HERE / "host_gate.py").read_bytes()) == manifest["host_gate_sha256"]
    assert sha(Path(__file__).read_bytes()) == manifest["verifier_sha256"]
    assert subprocess.run(["git", "merge-base", "--is-ancestor", SOURCE_HEAD,
                           "HEAD"], cwd=REPO).returncode == 0
    frozen = git("show", SOURCE_HEAD + ":research/notes/ecc2k130/n13_oaware_sparse_dense_sat_20260925/FROZEN.json")
    assert sha(frozen) == FREEZE_SHA256
    spec = json.loads(frozen)
    assert spec["status"] == "HELD_NO_MEASURED_PANEL"
    for field, name in (("run_sha256", "run.py"), ("verify_sha256", "verify.py"),
                        ("ci_replay_sha256", "ci_replay.py")):
        raw = git("show", SOURCE_HEAD + ":research/notes/ecc2k130/n13_oaware_sparse_dense_sat_20260925/" + name)
        assert sha(raw) == spec[field]
    expected = manifest["files_sha256"]
    assert manifest["file_count"] == len(expected) == 28
    actual = {str(path.relative_to(ARCHIVE)) for path in ARCHIVE.rglob("*")
              if path.is_file() and path.name != "MANIFEST.json"}
    assert actual == set(expected)
    for name, digest in expected.items():
        path = PurePosixPath(name)
        assert not path.is_absolute() and ".." not in path.parts
        assert sha((ARCHIVE / name).read_bytes()) == digest, name
    for mode, count in (("preflight", 0), ("smoke", 6)):
        receipt = json.loads((ARCHIVE / mode / "receipt.json").read_text())
        result = json.loads((ARCHIVE / mode / "result.json").read_text())
        assert receipt["mode"] == result["mode"] == mode
        assert receipt["decision"] == "COMPLETE"
        assert receipt["source_commit"] == SOURCE_HEAD
        assert receipt["freeze_sha256"] == result["freeze_sha256"] == FREEZE_SHA256
        assert result["pass"] is True and len(result["entries"]) == count
        assert receipt["process_wall_seconds"] > 0
    entries = json.loads((ARCHIVE / "smoke/result.json").read_text())["entries"]
    assert {(row["solver"], row["case"]) for row in entries} == {
        (solver, case) for solver in SOLVERS for case in ("sat", "unsat")}
    for row in entries:
        stem = row["solver"] + "-" + row["case"]
        assert row == json.loads((ARCHIVE / "smoke" / (stem + ".json")).read_text())
        cnf = SAT_CNF if row["case"] == "sat" else UNSAT_CNF
        assert (ARCHIVE / "smoke" / (stem + ".cnf")).read_bytes() == cnf
        assert row["input_sha256"] == sha(cnf)
        stdout = (ARCHIVE / "smoke" / (stem + ".stdout")).read_bytes()
        stderr = (ARCHIVE / "smoke" / (stem + ".stderr")).read_bytes()
        assert row["stdout_sha256"] == sha(stdout) and row["stdout_bytes"] == len(stdout)
        assert row["stderr_sha256"] == sha(stderr) and row["stderr_bytes"] == len(stderr)
        expected_exit = 10 if row["case"] == "sat" else 20
        expected_status = "s SATISFIABLE" if row["case"] == "sat" else "s UNSATISFIABLE"
        assert row["exit_code"] == expected_exit and row["stop_reason"] is None and row["pass"] is True
        lines = stdout.decode().splitlines()
        assert [line for line in lines if line.startswith("s ")] == [expected_status]
        if row["case"] == "sat":
            assert "v 1 0" in lines
    gate = json.loads((HERE / "evidence/host_admission_first_20260929.json").read_text())
    assert gate["schema"] == "n13-paired-host-admission-v1"
    assert gate["decision"] == "HELD" and gate["source_commit"] == SOURCE_HEAD
    assert gate["freeze_sha256"] == FREEZE_SHA256
    assert gate["frozen_preflight"] is True and gate["process_tree_monitor"] is True
    assert gate["logical_cpus"] == 14 and gate["disk_free_bytes"] >= 4 * 1024**3
    assert gate["load_1_5_15"][0] > gate["limits"]["load_1_max"]
    assert gate["load_1_5_15"][1] > gate["limits"]["load_5_max"]
    assert any(item["name"] == "cryptominisat5" for item in gate["competing_solvers"])
    print(json.dumps({"verdict": "PASS_PREFLIGHT_AND_SMOKE_NO_PANEL",
                      "source_head": SOURCE_HEAD, "freeze_sha256": FREEZE_SHA256,
                      "smoke_children": 6, "panel_children": 0,
                      "host_gate": "HELD"}, sort_keys=True))


if __name__ == "__main__":
    main()
