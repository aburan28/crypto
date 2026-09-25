#!/usr/bin/env python3
"""Cold staged n19 sparse export and independent replay with first receipt."""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value) -> None:
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def manifest(root: Path, freeze_sha: str) -> None:
    files = [{"path": str(path.relative_to(root)), "bytes": path.stat().st_size,
              "sha256": sha(path)} for path in sorted(root.rglob("*"))
             if path.is_file() and path.name != "MANIFEST.json"]
    save(root / "MANIFEST.json", {"freeze_sha256": freeze_sha, "files": files})


def verify_parent_merged(frozen: dict) -> None:
    result = subprocess.run(["git", "merge-base", "--is-ancestor",
                             frozen["required_parent_head"], "origin/main"],
                            cwd=REPO, capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError("#786 exact parent head is not yet in fetched origin/main")


def child(receipt: dict, out: Path, name: str, command: list[str],
          expected: Path, timeout: int) -> None:
    start = dt.datetime.now(dt.timezone.utc).isoformat()
    wall = time.perf_counter()
    try:
        result = subprocess.run(command, cwd=REPO, capture_output=True, text=True,
                                check=False, timeout=timeout)
        code, stdout, stderr, timed_out = result.returncode, result.stdout, result.stderr, False
    except subprocess.TimeoutExpired as error:
        code, timed_out = None, True
        stdout, stderr = error.stdout or "", error.stderr or ""
        if isinstance(stdout, bytes):
            stdout = stdout.decode(errors="replace")
        if isinstance(stderr, bytes):
            stderr = stderr.decode(errors="replace")
        stderr += f"\nexternal {timeout}s cap reached\n"
    stdout_path, stderr_path = out / f"{name}.stdout.txt", out / f"{name}.stderr.txt"
    stdout_path.write_text(stdout)
    stderr_path.write_text(stderr)
    receipt["attempts"].append({"phase": name, "command": command,
                                 "started_utc": start,
                                 "ended_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                                 "wall_seconds": time.perf_counter() - wall,
                                 "exit_code": code, "external_timeout": timed_out,
                                 "expected_sha256": sha(expected) if expected.exists() else None,
                                 "stdout_sha256": sha(stdout_path),
                                 "stderr_sha256": sha(stderr_path)})
    save(out / "receipt.json", receipt)
    if code != 0 or not expected.exists():
        raise RuntimeError(f"{name} failed or censored; first artifacts retained")


def run(out: Path):
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    for key, path in (("protocol_sha256", HERE / "PROTOCOL.md"),
                      ("export_sha256", HERE / "export.py"),
                      ("verify_sha256", HERE / "verify.py"),
                      ("run_sha256", Path(__file__)),
                      ("ci_replay_sha256", HERE / "ci_replay.py")):
        assert sha(path) == frozen[key], key
    spec = importlib.util.spec_from_file_location("n19_frozen_archive_gate", HERE / "ci_replay.py")
    assert spec is not None and spec.loader is not None
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)
    assert audit.check_freeze() == frozen
    verify_parent_merged(frozen)
    out.mkdir(parents=True, exist_ok=False)
    receipt = {"domain": frozen["domain"], "freeze_sha256": sha(HERE / "FROZEN.json"),
               "decision": "INCOMPLETE", "attempts": []}
    save(out / "receipt.json", receipt)
    try:
        producer = out / "producer"
        child(receipt, out, "export", [sys.executable, str(HERE / "export.py"),
                                        "--out", str(producer)],
              producer / "result.json", 195)
        verifier = out / "verify.json"
        child(receipt, out, "verify", [sys.executable, str(HERE / "verify.py"),
                                        "--produced", str(producer), "--out", str(verifier)],
              verifier, 630)
        produced_row = json.loads((producer / "result.json").read_text())
        verified_row = json.loads(verifier.read_text())
        assert verified_row["decision"] == "PASS"
        summary = {"decision": "PASS", "domain": frozen["domain"],
                   "growth": produced_row["growth"],
                   "transition_cases": produced_row["transition_cases"],
                   "primary_paths": produced_row["primary_paths"],
                   "signed_point_tuples": verified_row["signed_point_tuples"],
                   "target_labels": len(verified_row["target_rows"]),
                   "variables": produced_row["variables"],
                   "clauses": produced_row["clauses"],
                   "bytes": produced_row["bytes"],
                   "cold_children": {"export": {key: produced_row[key] for key in
                                                ("wall_seconds", "cpu_seconds", "peak_rss_bytes")},
                                     "verify": {key: verified_row[key] for key in
                                                ("wall_seconds", "cpu_seconds", "peak_rss_bytes")}},
                   "artifact_sha256": {"producer": sha(producer / "result.json"),
                                       "verifier": sha(verifier)}}
        save(out / "summary.json", summary)
        receipt["decision"] = "PASS"
    except Exception as error:
        failure = out / "producer/failure.json"
        verify_failure = out / "verify.json"
        statuses = []
        for path in (failure, verify_failure):
            if path.exists():
                try:
                    statuses.append(json.loads(path.read_text()).get("decision"))
                except (ValueError, OSError):
                    pass
        receipt["decision"] = "CENSORED" if "CENSORED" in statuses or any(
            row["external_timeout"] for row in receipt["attempts"]) else "FAILED"
        receipt["error"] = repr(error)
        save(out / "receipt.json", receipt)
        manifest(out, sha(HERE / "FROZEN.json"))
        raise
    save(out / "receipt.json", receipt)
    manifest(out, sha(HERE / "FROZEN.json"))
    print(json.dumps({"decision": "PASS", "receipt_sha256": sha(out / "receipt.json")},
                     sort_keys=True))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.out)


if __name__ == "__main__":
    main()
