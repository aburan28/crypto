#!/usr/bin/env python3
"""Pre-outcome source/target hash gate and archive-only paired SAT replay."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import analyze
import verify

HERE = Path(__file__).resolve().parent
SOLVERS = ("cryptominisat5", "kissat", "cadical")
REPS = ("dense", "sparse")
EXPORT_ORDER = ("dense", "sparse", "sparse", "dense")


def freeze():
    f = json.loads((HERE / "FROZEN.json").read_text())
    assert f["domain"] == "ecc2k130-n13-m5-paired-dense-sparse-cnf-v1"
    assert f["solver_order"] == list(SOLVERS) and f["export_order"] == list(EXPORT_ORDER)
    assert f["export_wall_cap_seconds"] == 180.0 and f["export_rss_cap_bytes"] == 512*1024**2
    assert f["solver_wall_cap_seconds"] == 15.0 and f["solver_rss_cap_bytes"] == 2*1024**3
    paths = {
        "protocol_sha256": HERE / "PROTOCOL.md",
        "run_sha256": HERE / "run.py",
        "verify_sha256": HERE / "verify.py",
        "ci_replay_sha256": HERE / "ci_replay.py",
        "analyze_sha256": HERE / "analyze.py",
        "workflow_sha256": Path(".github/workflows/ecc2k130-oaware-sparse-dense-sat.yml"),
        "dense_export_sha256": verify.DENSE / "export.py",
        "sparse_export_sha256": verify.SPARSE / "export.py",
        "dense_base_sha256": verify.DENSE_DIR / "base.cnf",
        "dense_schema_sha256": verify.DENSE_DIR / "schema.json",
        "sparse_base_sha256": verify.SPARSE_DIR / "base.cnf",
        "sparse_schema_sha256": verify.SPARSE_DIR / "schema.json",
        "sparse_verify_sha256": verify.SPARSE / "verify.py",
        "sparse_frozen_sha256": verify.SPARSE / "FROZEN.json",
        "sparse_evidence_receipt_sha256": verify.SPARSE / "evidence/final/receipt.json",
        "sparse_summary_sha256": verify.SPARSE / "evidence/final/summary.json",
        "prior_verify_sha256": verify.OLD / "verify.py",
        "prior_input_sha256": verify.OLD / "INPUT.json",
        "prior_frozen_sha256": verify.OLD / "FROZEN.json",
    }
    for key, path in paths.items():
        assert verify.sha(path) == f[key], (key, str(path))
    summary = json.loads((verify.SPARSE / "evidence/final/summary.json").read_text())
    assert summary["decision"] == "PASS"
    assert sorted(f["binaries"]) == sorted(SOLVERS)
    assert all(len(row["sha256"]) == 64 and row["version"] for row in f["binaries"].values())
    schemas, truth, curve, point, old = verify.schemas_and_truth()
    assert [r["id"] for r in schemas["dense"]["targets"][:32]] == [
        f"Q{q}T{t}" for q in range(8) for t in range(4)]
    return f, schemas, truth, curve, point, old, summary


def raw(row, root, stem):
    assert json.loads((root / f"{stem}.json").read_text()) == row
    for stream in ("stdout", "stderr"):
        p = root / f"{stem}.{stream}"
        assert verify.sha(p) == row[f"{stream}_sha256"]
        assert p.stat().st_size == row[f"{stream}_bytes"]


def replay_smoke(root, f, old):
    receipt = json.loads((root / "receipt.json").read_text())
    result = json.loads((root / "result.json").read_text())
    assert receipt["decision"] == "COMPLETE" and result["pass"]
    assert receipt["freeze_sha256"] == result["freeze_sha256"] == verify.sha(HERE / "FROZEN.json")
    assert len(result["entries"]) == 6
    for i, row in enumerate(result["entries"]):
        name, case = SOLVERS[i//2], ("sat", "unsat")[i%2]
        stem = f"{name}-{case}"
        raw(row, root, stem)
        assert row["solver"] == name and row["case"] == case
        assert verify.sha(root / f"{stem}.cnf") == row["input_sha256"]
        assert row["command"][0] == f["binaries"][name]["path"]
        try:
            status, assignment = old.parse_solver_output(
                (root / f"{stem}.stdout").read_bytes(), row["exit_code"], 1)
            expected = "SAT" if case == "sat" else "UNSAT"
            okay = status == expected and (case != "sat" or assignment.get(1) is True)
        except (ValueError, AssertionError):
            okay = False
        assert row["pass"] == (okay and row["stop_reason"] is None)
    return result


def replay_panel(root, f, schemas, truth, curve, point, old):
    receipt = json.loads((root / "receipt.json").read_text())
    result = json.loads((root / "result.json").read_text())
    assert receipt["decision"] == "COMPLETE"
    assert receipt["freeze_sha256"] == result["freeze_sha256"] == verify.sha(HERE / "FROZEN.json")
    assert len(result["exports"]) == 4 and len(result["entries"]) == 192
    for index, rep in enumerate(EXPORT_ORDER):
        row = result["exports"][index]
        stem = f"export{index}-{rep}"
        raw(row, root, stem)
        assert row["representation"] == rep and row["pair"] == (0 if index < 2 else 1)
        assert row["exit_code"] == 0 and row["stop_reason"] is None
        assert row["command"][0] == f["python_executable"]
        assert row["sampled_peak_tree_rss_bytes"] <= f["export_rss_cap_bytes"]
        assert row["expected_base_sha256"] == f[f"{rep}_base_sha256"]
        assert row["expected_schema_sha256"] == f[f"{rep}_schema_sha256"]
        producer = root / f"{stem}.result.json"
        assert verify.sha(producer) == row["producer_result_sha256"]
        n13 = next(r for r in json.loads(producer.read_text())["panels"] if r["panel"] == "n13-m5")
        assert n13["base_sha256"] == f[f"{rep}_base_sha256"]
        assert n13["schema_sha256"] == f[f"{rep}_schema_sha256"]
    bases = {rep: (verify.DENSE_DIR if rep == "dense" else verify.SPARSE_DIR) / "base.cnf"
             for rep in REPS}
    count, statuses = 0, []
    for ti, target in enumerate(schemas["dense"]["targets"][:32]):
        query_hash = {}
        query_bytes = {}
        for rep in REPS:
            data = verify.query_bytes(old, bases[rep].read_bytes(), target, schemas[rep])
            query_hash[rep] = hashlib.sha256(data).hexdigest()
            query_bytes[rep] = len(data)
        engines = SOLVERS[ti%3:] + SOLVERS[:ti%3]
        reps = REPS if ti%2 == 0 else tuple(reversed(REPS))
        for j, name in enumerate(engines):
            for k, rep in enumerate(reps):
                row = result["entries"][6*ti + 2*j + k]
                stem = f"{target['id']}-{name}-{rep}"
                raw(row, root, stem)
                assert (row["id"], row["solver"], row["representation"]) == (
                    target["id"], name, rep)
                assert row["assumption_literal"] == target["assumption_literal"]
                assert row["point_oracle_positive"] == truth[target["id"]]
                assert row["input_sha256"] == query_hash[rep]
                assert row["input_bytes"] == query_bytes[rep]
                assert row["command"][0] == f["binaries"][name]["path"]
                assert row["command"][-1] == str(root / f"query-{rep}.cnf")
                if row["stop_reason"]:
                    verdict, certificate = "CENSORED", None
                else:
                    try:
                        verdict, certificate = verify.classify(
                            old, (root / f"{stem}.stdout").read_bytes(),
                            row["exit_code"], schemas[rep], target, truth, curve, point)
                    except (ValueError, AssertionError):
                        verdict, certificate = "INVALID", None
                assert row["verdict"] == verdict and row["certificate"] == certificate
                statuses.append(verdict)
                count += 1
    return receipt, result, {s: statuses.count(s) for s in sorted(set(statuses))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", type=Path)
    parser.add_argument("--panel", type=Path)
    parser.add_argument("--analysis", type=Path)
    args = parser.parse_args()
    f, schemas, truth, curve, point, old, audit = freeze()
    output = {"decision": "HASH_PASS", "freeze_sha256": verify.sha(HERE / "FROZEN.json")}
    if args.smoke:
        smoke = replay_smoke(args.smoke, f, old)
        output["smoke_pass"] = smoke["pass"]
    if args.panel:
        assert args.smoke
        receipt, panel, statuses = replay_panel(args.panel, f, schemas, truth, curve, point, old)
        output["panel_queries"] = len(panel["entries"])
        output["status_counts"] = statuses
        output["proof_checked_unsat"] = False
        if args.analysis:
            saved = json.loads(args.analysis.read_text())
            audit_wall = {rep+"_replay": audit["cold_children"][rep+"_replay"]["wall_seconds"]
                          for rep in REPS}
            fresh = analyze.summarize(smoke, receipt, panel, audit_wall)
            assert saved["decision"] == fresh["decision"]
            assert saved["proof_checked_unsat"] is False
            for engine in SOLVERS:
                assert saved["engines"][engine]["paired_sparse_over_dense_full_wall"] == (
                    fresh["engines"][engine]["paired_sparse_over_dense_full_wall"])
    print(json.dumps(output, sort_keys=True))


if __name__ == "__main__":
    main()
