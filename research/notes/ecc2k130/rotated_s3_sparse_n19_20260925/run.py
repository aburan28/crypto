#!/usr/bin/env python3
"""Cold staged n19 sparse export and independent replay with first receipt."""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import os
import platform
import signal
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
    for key in ("preregistered_parent_head", "required_parent_head",
                "required_parent_merge_commit"):
        result = subprocess.run(["git", "merge-base", "--is-ancestor",
                                 frozen[key], "origin/main"],
                                cwd=REPO, capture_output=True, check=False)
        if result.returncode != 0:
            raise RuntimeError(f"#786 {key} is not in fetched origin/main")


def require_linux_proc() -> None:
    if sys.platform != "linux" or not Path("/proc/self/status").is_file():
        raise RuntimeError("Linux /proc is required before the first attempt")


def process_group_snapshot(pgid: int) -> tuple[int, int]:
    """Return sampled resident bytes and non-zombie processes in a session group."""
    total = live = 0
    for entry in Path("/proc").iterdir():
        if not entry.name.isdecimal():
            continue
        try:
            stat = (entry / "stat").read_text()
            fields = stat[stat.rfind(")") + 2:].split()
            if int(fields[2]) != pgid or fields[0] == "Z":
                continue
            live += 1
            for line in (entry / "status").read_text().splitlines():
                if line.startswith("VmRSS:"):
                    total += int(line.split()[1]) * 1024
                    break
        except (FileNotFoundError, PermissionError, ProcessLookupError, ValueError):
            continue
    return total, live


def child(receipt: dict, out: Path, name: str, command: list[str],
          expected: Path, timeout: float, rss_cap: int) -> None:
    require_linux_proc()
    start = dt.datetime.now(dt.timezone.utc).isoformat()
    wall = time.perf_counter()
    stdout_path, stderr_path = out / f"{name}.stdout.txt", out / f"{name}.stderr.txt"
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = subprocess.Popen(command, cwd=REPO, stdout=stdout, stderr=stderr,
                                   start_new_session=True)
        stop = monitor_error = None
        sampled_peak = 0
        status = usage = None
        try:
            while True:
                try:
                    rss, _ = process_group_snapshot(process.pid)
                    sampled_peak = max(sampled_peak, rss)
                    if sampled_peak >= rss_cap:
                        stop = "PROCESS_GROUP_RSS_CAP"
                    elif time.perf_counter() - wall >= timeout:
                        stop = "EXTERNAL_WALL_CAP"
                    if stop is not None:
                        break
                    pid, status, usage = os.wait4(process.pid, os.WNOHANG)
                    if pid:
                        break
                    status = usage = None
                    time.sleep(0.05)
                except Exception as error:
                    stop = "MONITOR_ERROR"
                    monitor_error = repr(error)
                    break
        finally:
            # Also terminate descendants if the direct child exited first.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            if usage is None:
                _, status, usage = os.wait4(process.pid, 0)
            quiescence_start = time.perf_counter()
            group_quiesced = False
            while time.perf_counter() - quiescence_start < 5:
                _, live = process_group_snapshot(process.pid)
                if live == 0:
                    group_quiesced = True
                    break
                time.sleep(0.05)
            group_quiescence_seconds = time.perf_counter() - quiescence_start
        code = os.waitstatus_to_exitcode(status)
        process.returncode = code
        # wait4 covers a direct child's transient RSS, while descendants have
        # only the 50 ms process-group samples. This is a sampled group cap.
        maxrss = usage.ru_maxrss * 1024
        sampled_peak = max(sampled_peak, maxrss)
        if sampled_peak >= rss_cap and stop is None:
            stop = "PEAK_RSS_CAP_AT_EXIT"
        if not group_quiesced:
            stop = "PROCESS_GROUP_LEAK"
        if stop is not None:
            stderr.write(f"\nexternal resource stop: {stop}\n".encode())
        if monitor_error is not None:
            stderr.write(f"monitor error: {monitor_error}\n".encode())
    receipt["attempts"].append({"phase": name, "command": command,
                                 "started_utc": start,
                                 "ended_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                                 "wall_seconds": time.perf_counter() - wall,
                                 "exit_code": code,
                                 "external_timeout": stop == "EXTERNAL_WALL_CAP",
                                 "resource_stop": stop,
                                 "rss_cap_bytes": rss_cap,
                                 "sampled_group_peak_rss_bytes": sampled_peak,
                                 "direct_child_peak_rss_bytes": maxrss,
                                 "group_quiesced": group_quiesced,
                                 "group_quiescence_seconds": group_quiescence_seconds,
                                 "user_cpu_seconds": usage.ru_utime,
                                 "system_cpu_seconds": usage.ru_stime,
                                 "host": platform.platform(),
                                 "expected_sha256": sha(expected) if expected.exists() else None,
                                 "stdout_sha256": sha(stdout_path),
                                 "stderr_sha256": sha(stderr_path)})
    save(out / "receipt.json", receipt)
    if stop is not None or code != 0 or not expected.exists():
        raise RuntimeError(f"{name} failed or censored; first artifacts retained")


def run(out: Path, release_gate: Path):
    out = out.resolve()
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    for key, path in (("protocol_sha256", HERE / "PROTOCOL.md"),
                      ("export_sha256", HERE / "export.py"),
                      ("verify_sha256", HERE / "verify.py"),
                      ("run_sha256", Path(__file__)),
                      ("ci_replay_sha256", HERE / "ci_replay.py"),
                      ("resource_test_sha256", HERE / "test_resource_cap.py"),
                      ("release_gate_sha256", HERE / "release_gate.py")):
        assert sha(path) == frozen[key], key
    spec = importlib.util.spec_from_file_location("n19_frozen_archive_gate", HERE / "ci_replay.py")
    assert spec is not None and spec.loader is not None
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)
    assert audit.check_freeze() == frozen
    verify_parent_merged(frozen)
    require_linux_proc()
    gate = json.loads(release_gate.read_text())
    head = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                   cwd=REPO, text=True).strip()
    if not (gate["decision"] == "FIRST_RELEASE_APPROVED"
            and gate["pr"] == frozen["release_pr_number"]
            and gate["head_sha"] == head
            and gate["repository"] == "aburan28/crypto"
            and gate["run_attempt"] == 1
            and isinstance(gate["run_id"], int) and gate["run_id"] > 0
            and isinstance(gate["label_event_id"], int) and gate["label_event_id"] > 0
            and gate["freeze_sha256"] == sha(HERE / "FROZEN.json")):
        raise RuntimeError("first-release gate receipt does not bind this source and run")
    out.mkdir(parents=True, exist_ok=False)
    (out / "release_gate.json").write_bytes(release_gate.read_bytes())
    receipt = {"domain": frozen["domain"], "freeze_sha256": sha(HERE / "FROZEN.json"),
               "release_gate_sha256": sha(out / "release_gate.json"),
               "release_head_sha": head, "release_run_id": gate["run_id"],
               "release_label_event_id": gate["label_event_id"],
               "decision": "INCOMPLETE", "attempts": [],
               "run_out": str(out), "source_root": str(HERE),
               "python_executable": sys.executable,
               "python_version": platform.python_version()}
    save(out / "receipt.json", receipt)
    try:
        producer = out / "producer"
        child(receipt, out, "export", [sys.executable, str(HERE / "export.py"),
                                        "--out", str(producer)],
              producer / "result.json", frozen["caps"]["external_export_seconds"],
              frozen["caps"]["child_rss_bytes"])
        verifier = out / "verify.json"
        child(receipt, out, "verify", [sys.executable, str(HERE / "verify.py"),
                                        "--produced", str(producer), "--out", str(verifier)],
              verifier, frozen["caps"]["external_verify_seconds"],
              frozen["caps"]["child_rss_bytes"])
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
            row["resource_stop"] is not None or row["external_timeout"] or
            (row["exit_code"] is not None and row["exit_code"] < 0)
            for row in receipt["attempts"]) else "FAILED"
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
    parser.add_argument("--release-gate", type=Path, required=True)
    args = parser.parse_args()
    run(args.out, args.release_gate)


if __name__ == "__main__":
    main()
