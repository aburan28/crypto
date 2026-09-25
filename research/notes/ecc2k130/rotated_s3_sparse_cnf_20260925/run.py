#!/usr/bin/env python3
"""Four cold subprocesses: dense reference/replay, then sparse export/replay."""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
DENSE = HERE.parent / "rotated_s3_o_branch_20260925"
EXTERNAL_CAP = 195


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value):
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def manifest(root: Path, freeze_sha: str):
    files = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.name != "MANIFEST.json":
            files.append({"path": str(path.relative_to(root)), "bytes": path.stat().st_size,
                          "sha256": sha(path)})
    save(root / "MANIFEST.json", {"freeze_sha256": freeze_sha, "files": files})


def child(receipt: dict, root: Path, phase: str, command: list[str], expected: Path):
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    wall = time.perf_counter()
    try:
        process = subprocess.run(command, capture_output=True, text=True,
                                 check=False, timeout=EXTERNAL_CAP)
        exit_code, stdout, stderr, timed_out = (process.returncode,
                                                 process.stdout, process.stderr, False)
    except subprocess.TimeoutExpired as error:
        exit_code, timed_out = None, True
        stdout, stderr = error.stdout or "", error.stderr or ""
        if isinstance(stdout, bytes):
            stdout = stdout.decode(errors="replace")
        if isinstance(stderr, bytes):
            stderr = stderr.decode(errors="replace")
        stderr += f"\nexternal {EXTERNAL_CAP}-second child cap reached\n"
    stdout_path = root / f"{phase}.stdout.txt"
    stderr_path = root / f"{phase}.stderr.txt"
    stdout_path.write_text(stdout)
    stderr_path.write_text(stderr)
    receipt["attempts"].append({
        "phase": phase, "command": command, "started_utc": started,
        "ended_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "wall_seconds": time.perf_counter() - wall, "exit_code": exit_code,
        "external_timeout": timed_out,
        "expected_sha256": sha(expected) if expected.exists() else None,
        "stdout_sha256": sha(stdout_path), "stderr_sha256": sha(stderr_path)})
    save(root / "receipt.json", receipt)
    if exit_code != 0 or not expected.exists():
        raise RuntimeError(f"{phase} failed/censored; retained at {root}")


def check_dense(root: Path, frozen: dict):
    result = json.loads((root / "dense/producer/result.json").read_text())
    for row in result["panels"]:
        name = row["panel"]
        assert name in frozen["panels"]
        for file, key in (("base.cnf", "base_sha256"),
                          ("schema.json", "schema_sha256"),
                          ("paths.jsonl.gz", "paths_sha256")):
            path = root / "dense/producer" / name / file
            assert sha(path) == row[key] == frozen["panels"][name][file]["sha256"]
    assert json.loads((root / "dense/verify.json").read_text())["decision"] == "PASS"
    return result


def summarize(root: Path, frozen: dict):
    dense = check_dense(root, frozen)
    sparse = json.loads((root / "sparse/producer/result.json").read_text())
    dense_verify = json.loads((root / "dense/verify.json").read_text())
    sparse_verify = json.loads((root / "sparse/verify.json").read_text())
    assert sparse_verify["decision"] == "PASS"
    assert [p["panel"] for p in dense["panels"]] == [p["panel"] for p in sparse["panels"]]
    panels = []
    for d, s in zip(dense["panels"], sparse["panels"]):
        assert (s["dense_variables"], s["dense_clauses"], s["dense_bytes"]) == (
            d["variables"], d["clauses"], (root / "dense/producer" / d["panel"] / "base.cnf").stat().st_size)
        panels.append({"panel": d["panel"],
                       "dense": {"variables": d["variables"], "clauses": d["clauses"],
                                 "bytes": s["dense_bytes"]},
                       "sparse": {"variables": s["variables"], "clauses": s["clauses"],
                                  "bytes": s["bytes"]},
                       "exact_ratios_sparse_over_dense": {
                           "variables": [s["variables"], d["variables"]],
                           "clauses": [s["clauses"], d["clauses"]],
                           "bytes": [s["bytes"], s["dense_bytes"]]},
                       "candidate_paths": s["model_paths"],
                       "target_labels": s["target_labels"]})
    return {"decision": "PASS", "domain": frozen["domain"], "panels": panels,
            "cold_children": {phase: {"wall_seconds": row["wall_seconds"],
                                      "cpu_seconds": row["cpu_seconds"],
                                      "peak_rss_bytes": row["peak_rss_bytes"]}
                              for phase, row in (("dense_export", dense),
                                                 ("dense_replay", dense_verify),
                                                 ("sparse_export", sparse),
                                                 ("sparse_replay", sparse_verify))},
            "child_sha256": {"dense_export": sha(root / "dense/producer/result.json"),
                             "dense_replay": sha(root / "dense/verify.json"),
                             "sparse_export": sha(root / "sparse/producer/result.json"),
                             "sparse_replay": sha(root / "sparse/verify.json")}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert sha(Path(__file__)) == frozen["runner_sha256"]
    for key, path in (("protocol_sha256", HERE / "PROTOCOL.md"),
                      ("export_sha256", HERE / "export.py"),
                      ("verify_sha256", HERE / "verify.py"),
                      ("ci_replay_sha256", HERE / "ci_replay.py")):
        assert sha(path) == frozen[key]
    args.out.mkdir(parents=True, exist_ok=False)
    receipt = {"protocol": frozen["domain"], "freeze_sha256": sha(HERE / "FROZEN.json"),
               "decision": "INCOMPLETE", "attempts": []}
    save(args.out / "receipt.json", receipt)
    try:
        dense_root = args.out / "dense"
        dense_root.mkdir()
        child(receipt, args.out, "dense_export",
              [sys.executable, str(DENSE / "export.py"), "--out", str(dense_root / "producer")],
              dense_root / "producer/result.json")
        for row in json.loads((dense_root / "producer/result.json").read_text())["panels"]:
            for file in ("base.cnf", "schema.json", "paths.jsonl.gz"):
                assert sha(dense_root / "producer" / row["panel"] / file) == frozen["panels"][row["panel"]][file]["sha256"]
        child(receipt, args.out, "dense_replay",
              [sys.executable, str(DENSE / "verify.py"), "--producer", str(dense_root / "producer"),
               "--out", str(dense_root / "verify.json")], dense_root / "verify.json")
        sparse_root = args.out / "sparse"
        sparse_root.mkdir()
        child(receipt, args.out, "sparse_export",
              [sys.executable, str(HERE / "export.py"), "--out", str(sparse_root / "producer")],
              sparse_root / "producer/result.json")
        child(receipt, args.out, "sparse_replay",
              [sys.executable, str(HERE / "verify.py"), "--sparse", str(sparse_root / "producer"),
               "--out", str(sparse_root / "verify.json")], sparse_root / "verify.json")
        summary = summarize(args.out, frozen)
        save(args.out / "summary.json", summary)
        receipt["decision"] = "PASS"
    except Exception as error:
        receipt["decision"] = "CENSORED_OR_FAILED"
        receipt["error"] = repr(error)
        save(args.out / "receipt.json", receipt)
        manifest(args.out, sha(HERE / "FROZEN.json"))
        raise
    save(args.out / "receipt.json", receipt)
    manifest(args.out, sha(HERE / "FROZEN.json"))
    print(json.dumps({"decision": "PASS", "receipt_sha256": sha(args.out / "receipt.json")},
                     sort_keys=True))


if __name__ == "__main__":
    main()
