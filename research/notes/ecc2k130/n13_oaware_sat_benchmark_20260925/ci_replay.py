#!/usr/bin/env python3
"""Hash-only prereg gate and solver-free raw archive replay."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import verify

HERE = Path(__file__).resolve().parent
SOLVERS = ("cryptominisat5", "kissat", "cadical")


def freeze():
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert frozen["domain"] == "ecc2k130-n13-m5-oaware-cnf-sat-stage-v1"
    assert frozen["source_main_merge"] == "8c5178b00b2af91548a5b4f88558723b3ad656c9"
    assert frozen["solver_order"] == list(SOLVERS)
    assert frozen["wall_cap_seconds"] == 15.0
    assert frozen["rss_cap_bytes"] == 2 * 1024**3
    paths = {
        "protocol_sha256": HERE / "PROTOCOL.md",
        "input_sha256": HERE / "INPUT.json",
        "run_sha256": HERE / "run.py",
        "verify_sha256": HERE / "verify.py",
        "ci_replay_sha256": HERE / "ci_replay.py",
        "analyze_sha256": HERE / "analyze.py",
        "base_sha256": verify.CNF_DIR / "base.cnf",
        "schema_sha256": verify.CNF_DIR / "schema.json",
        "oaware_export_sha256": verify.CNF_DIR.parent.parent.parent / "export.py",
        "oaware_verify_sha256": verify.CNF_DIR.parent.parent.parent / "verify.py",
        "point_verify_sha256": verify.CORPUS_VERIFY,
        "point_corpus_sha256": verify.CORPUS,
        "workflow_sha256": Path(".github/workflows/ecc2k130-n13-oaware-sat-benchmark.yml"),
    }
    for key, path in paths.items():
        assert verify.sha(path) == frozen[key], (key, str(path))
    schema, truth, curve, parent = verify.schema_and_truth()
    assert [row["id"] for row in schema["targets"][:32]] == [
        f"Q{q}T{t}" for q in range(8) for t in range(4)]
    assert sum(truth.values()) == 5
    inputs = json.loads((HERE / "INPUT.json").read_text())
    assert inputs["domain"] == frozen["domain"]
    assert inputs["targets"] == [{"id": row["id"], "literal": row["assumption_literal"],
                                   "point": row["point"], "point_oracle_positive": truth[row["id"]]}
                                  for row in schema["targets"][:32]]
    assert sorted(frozen["binaries"]) == sorted(SOLVERS)
    assert all(len(item["sha256"]) == 64 and item["version"]
               for item in frozen["binaries"].values())
    return frozen, schema, truth, curve, parent


def raw_matches(row, root, stem):
    for stream in ("stdout", "stderr"):
        path = root / f"{stem}.{stream}"
        assert verify.sha(path) == row[f"{stream}_sha256"]
        assert path.stat().st_size == row[f"{stream}_bytes"]
    assert json.loads((root / f"{stem}.json").read_text()) == row


def replay_smoke(path, frozen):
    result = json.loads((path / "result.json").read_text())
    assert result["mode"] == "smoke" and result["freeze_sha256"] == verify.sha(HERE / "FROZEN.json")
    assert len(result["entries"]) == 6
    successes = []
    for i, row in enumerate(result["entries"]):
        name, case = SOLVERS[i // 2], ("sat", "unsat")[i % 2]
        stem = f"{name}-{case}"
        assert row["solver"] == name and row["case"] == case
        raw_matches(row, path, stem)
        assert verify.sha(path / f"{stem}.cnf") == row["input_sha256"]
        assert row["command"][0] == frozen["binaries"][name]["path"]
        assert row["command"][-1] == str(path / f"{stem}.cnf")
        try:
            status, assignment = verify.parse_solver_output(
                (path / f"{stem}.stdout").read_bytes(), row["exit_code"], 1)
            expected = "SAT" if case == "sat" else "UNSAT"
            okay = status == expected and (case != "sat" or assignment.get(1) is True)
        except (ValueError, AssertionError):
            okay = False
        assert row["pass"] == (okay and row["stop_reason"] is None)
        successes.append(row["pass"])
    assert result["pass"] == all(successes)
    return result["pass"]


def replay_panel(path, frozen, schema, truth, curve, parent):
    result = json.loads((path / "result.json").read_text())
    assert result["mode"] == "panel" and result["freeze_sha256"] == verify.sha(HERE / "FROZEN.json")
    assert result["base_sha256"] == frozen["base_sha256"]
    assert len(result["entries"]) == 96
    base = (verify.CNF_DIR / "base.cnf").read_bytes()
    checked, statuses = 0, []
    for i, target in enumerate(schema["targets"][:32]):
        raw = verify.query_bytes(base, target["assumption_literal"],
                                 schema["variables"], schema["clauses"])
        expected_sha = hashlib.sha256(raw).hexdigest()
        rotated = SOLVERS[i % 3:] + SOLVERS[:i % 3]
        for j, name in enumerate(rotated):
            row = result["entries"][3 * i + j]
            stem = f"{target['id']}-{name}"
            assert (row["id"], row["solver"], row["assumption_literal"]) == (
                target["id"], name, target["assumption_literal"])
            assert row["input_sha256"] == expected_sha and row["input_bytes"] == len(raw)
            assert row["command"][:-1] == ([frozen["binaries"][name]["path"],
                   "--verb=0", "--threads=1"] if name == "cryptominisat5" else
                   [frozen["binaries"][name]["path"]])
            assert row["command"][-1] == str(path / "query.cnf")
            assert row["point_oracle_positive"] == truth[target["id"]]
            raw_matches(row, path, stem)
            if row["stop_reason"]:
                verdict, certificate = "CENSORED", None
            else:
                try:
                    status, assignment = verify.parse_solver_output(
                        (path / f"{stem}.stdout").read_bytes(), row["exit_code"],
                        schema["variables"])
                    if status == "SAT":
                        certificate = verify.model_certificate(
                            schema, target, assignment, curve, parent)
                        verdict = "SAT" if truth[target["id"]] else "CONTRADICTION"
                    elif status == "UNSAT":
                        certificate = None
                        verdict = "UNSAT" if not truth[target["id"]] else "CONTRADICTION"
                    else:
                        certificate, verdict = None, "CENSORED"
                except (ValueError, AssertionError):
                    certificate, verdict = None, "INVALID"
            assert row["verdict"] == verdict and row["certificate"] == certificate
            assert row["wall_seconds"] >= 0 and row["query_setup_wall_seconds"] >= 0
            assert row["verifier_wall_seconds"] >= 0 and row["sampled_peak_rss_bytes"] >= 0
            checked += 1
            statuses.append(verdict)
    return {"queries": checked, "status_counts": {k: statuses.count(k) for k in sorted(set(statuses))},
            "proof_checked_unsat": False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", type=Path)
    parser.add_argument("--panel", type=Path)
    args = parser.parse_args()
    frozen, schema, truth, curve, parent = freeze()
    result = {"decision": "HASH_PASS", "freeze_sha256": verify.sha(HERE / "FROZEN.json")}
    if args.smoke:
        result["smoke_pass"] = replay_smoke(args.smoke, frozen)
    if args.panel:
        result["panel"] = replay_panel(args.panel, frozen, schema, truth, curve, parent)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
