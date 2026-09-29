#!/usr/bin/env python3
"""Frozen cold CNF-solver panel; run only after prereg CI and host release."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import resource
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import psutil
import verify

HERE = Path(__file__).resolve().parent
SOLVERS = ("cryptominisat5", "kissat", "cadical")
CAP_WALL = 15.0
CAP_RSS = 2 * 1024**3
POLL = 0.02


def utc():
    return datetime.now(timezone.utc).isoformat()


def save(path, value):
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def sha_bytes(raw):
    return hashlib.sha256(raw).hexdigest()


def rss(process):
    if process is None:
        return 0
    try:
        family = [process] + process.children(recursive=True)
        return sum(item.memory_info().rss for item in family if item.is_running())
    except (psutil.Error, ProcessLookupError):
        return 0


def ru_cpu():
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    return (r.ru_utime, r.ru_stime, r.ru_maxrss)


def execute(command, stdout_path, stderr_path):
    before = ru_cpu()
    start_utc = utc()
    start = time.perf_counter()
    with stdout_path.open("wb") as out, stderr_path.open("wb") as err:
        child = subprocess.Popen(command, stdout=out, stderr=err, start_new_session=True)
        try:
            proc = psutil.Process(child.pid)
        except psutil.NoSuchProcess:
            proc = None
        peak = 0
        reason = None
        while child.poll() is None:
            peak = max(peak, rss(proc))
            if peak > CAP_RSS:
                reason = "rss_cap"
            elif time.perf_counter() - start > CAP_WALL:
                reason = "wall_cap"
            if reason is not None:
                try:
                    os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                break
            time.sleep(POLL)
        exit_code = child.wait()
        peak = max(peak, rss(proc))
    after = ru_cpu()
    # Darwin ru_maxrss is bytes, Linux is KiB; this high water spans prior children.
    scale = 1 if sys.platform == "darwin" else 1024
    return {"command": command, "start_utc": start_utc, "end_utc": utc(),
            "wall_seconds": time.perf_counter() - start, "exit_code": exit_code,
            "stop_reason": reason, "sampled_peak_rss_bytes": peak,
            "children_ru_maxrss_before_bytes": before[2] * scale,
            "children_ru_maxrss_after_bytes": after[2] * scale,
            "user_cpu_seconds": after[0] - before[0],
            "system_cpu_seconds": after[1] - before[1],
            "stdout_sha256": verify.sha(stdout_path), "stderr_sha256": verify.sha(stderr_path),
            "stdout_bytes": stdout_path.stat().st_size, "stderr_bytes": stderr_path.stat().st_size}


def command_for(name, binary, cnf):
    if name == "cryptominisat5":
        return [str(binary), "--verb=0", "--threads=1", str(cnf)]
    return [str(binary), str(cnf)]


def parse_and_certify(record, target, schema, truth, curve, parent, stdout):
    if record["stop_reason"]:
        return "CENSORED", None, record["stop_reason"]
    try:
        verdict, assignment = verify.parse_solver_output(
            stdout.read_bytes(), record["exit_code"], schema["variables"])
        if verdict == "SAT":
            certificate = verify.model_certificate(schema, target, assignment, curve, parent)
            if not truth[target["id"]]:
                return "CONTRADICTION", certificate, "SAT contradicts complete point oracle"
            return "SAT", certificate, None
        if verdict == "UNSAT":
            if truth[target["id"]]:
                return "CONTRADICTION", None, "UNSAT contradicts complete point oracle"
            return "UNSAT", None, None
        return "CENSORED", None, "solver UNKNOWN"
    except (ValueError, AssertionError) as exc:
        return "INVALID", None, f"{type(exc).__name__}: {exc}"


def preflight():
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    checks = {
        "protocol_sha256": HERE / "PROTOCOL.md",
        "input_sha256": HERE / "INPUT.json",
        "run_sha256": HERE / "run.py",
        "verify_sha256": HERE / "verify.py",
        "ci_replay_sha256": HERE / "ci_replay.py",
        "base_sha256": verify.CNF_DIR / "base.cnf",
        "schema_sha256": verify.CNF_DIR / "schema.json",
        "oaware_export_sha256": verify.CNF_DIR.parent.parent.parent / "export.py",
        "oaware_verify_sha256": verify.CNF_DIR.parent.parent.parent / "verify.py",
        "point_verify_sha256": verify.CORPUS_VERIFY,
        "point_corpus_sha256": verify.CORPUS,
    }
    for key, path in checks.items():
        assert verify.sha(path) == frozen[key], key
    assert frozen["source_main_merge"] == "8c5178b00b2af91548a5b4f88558723b3ad656c9"
    assert frozen["wall_cap_seconds"] == CAP_WALL and frozen["rss_cap_bytes"] == CAP_RSS
    assert frozen["solver_order"] == list(SOLVERS)
    assert platform.python_version() == frozen["python_version"]
    assert psutil.__version__ == frozen["psutil_version"]
    binaries = {}
    for name in SOLVERS:
        binary = Path(frozen["binaries"][name]["path"])
        assert binary.is_file() and verify.sha(binary) == frozen["binaries"][name]["sha256"]
        binaries[name] = binary
    return frozen, binaries


def smoke(outdir, frozen, binaries, schema, truth, curve, parent):
    # SAT/UNSAT checks only parser/exit convention; the 32 corpus queries remain untouched.
    templates = (("sat", b"p cnf 1 1\n1 0\n", "SAT"),
                 ("unsat", b"p cnf 1 2\n1 0\n-1 0\n", "UNSAT"))
    entries = []
    for name in SOLVERS:
        for case, raw, expected in templates:
            stem = f"{name}-{case}"
            cnf = outdir / f"{stem}.cnf"
            cnf.write_bytes(raw)
            stdout, stderr = outdir / f"{stem}.stdout", outdir / f"{stem}.stderr"
            record = execute(command_for(name, binaries[name], cnf), stdout, stderr)
            try:
                status, assignment = verify.parse_solver_output(
                    stdout.read_bytes(), record["exit_code"], 1)
                valid = status == expected and (expected != "SAT" or assignment.get(1) is True)
            except (ValueError, AssertionError) as exc:
                status, valid = f"parse-error: {exc}", False
            record.update({"solver": name, "case": case, "input_sha256": sha_bytes(raw),
                           "parsed_status": status, "pass": valid and not record["stop_reason"]})
            save(outdir / f"{stem}.json", record)
            entries.append(record)
    result = {"mode": "smoke", "freeze_sha256": verify.sha(HERE / "FROZEN.json"),
              "start_utc": entries[0]["start_utc"], "end_utc": utc(),
              "pass": all(row["pass"] for row in entries), "entries": entries}
    save(outdir / "result.json", result)
    return result


def panel(outdir, frozen, binaries, schema, truth, curve, parent):
    base_path = verify.CNF_DIR / "base.cnf"
    start = time.perf_counter()
    start_utc = utc()
    base = base_path.read_bytes()
    assert sha_bytes(base) == frozen["base_sha256"]
    initial_setup = time.perf_counter() - start
    results = []
    targets = schema["targets"][:32]
    for ti, target in enumerate(targets):
        qstart = time.perf_counter()
        raw = verify.query_bytes(base, target["assumption_literal"],
                                 schema["variables"], schema["clauses"])
        cnf = outdir / "query.cnf"
        cnf.write_bytes(raw)
        query_hash = sha_bytes(raw)
        query_setup = time.perf_counter() - qstart
        rotated = SOLVERS[ti % 3:] + SOLVERS[:ti % 3]
        for name in rotated:
            stem = f"{target['id']}-{name}"
            stdout, stderr = outdir / f"{stem}.stdout", outdir / f"{stem}.stderr"
            record = execute(command_for(name, binaries[name], cnf), stdout, stderr)
            verify_start = time.perf_counter()
            verdict, certificate, error = parse_and_certify(
                record, target, schema, truth, curve, parent, stdout)
            verify_wall = time.perf_counter() - verify_start
            record.update({"id": target["id"], "solver": name, "assumption_literal":
                           target["assumption_literal"], "input_sha256": query_hash,
                           "input_bytes": len(raw), "query_setup_wall_seconds": query_setup,
                           "verifier_wall_seconds": verify_wall, "verdict": verdict,
                           "point_oracle_positive": truth[target["id"]],
                           "certificate": certificate, "error": error})
            save(outdir / f"{stem}.json", record)
            results.append(record)
        cnf.unlink()
    result = {"mode": "panel", "freeze_sha256": verify.sha(HERE / "FROZEN.json"),
              "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
              "host": platform.platform(), "python": sys.version, "psutil": psutil.__version__,
              "start_utc": start_utc, "end_utc": utc(),
              "full_runner_wall_seconds": time.perf_counter() - start,
              "base_load_wall_seconds": initial_setup, "base_sha256": frozen["base_sha256"],
              "entries": results}
    save(outdir / "result.json", result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("smoke", "panel"))
    parser.add_argument("outdir", type=Path)
    args = parser.parse_args()
    assert not args.outdir.exists(), "never overwrite an attempt"
    process_start = time.perf_counter()
    frozen, binaries = preflight()
    schema, truth, curve, parent = verify.schema_and_truth()
    preflight_wall = time.perf_counter() - process_start
    args.outdir.mkdir(parents=True)
    result = (smoke if args.mode == "smoke" else panel)(
        args.outdir, frozen, binaries, schema, truth, curve, parent)
    result["preflight_wall_seconds"] = preflight_wall
    result["full_process_wall_seconds"] = time.perf_counter() - process_start
    save(args.outdir / "result.json", result)
    print(json.dumps({"mode": args.mode, "pass": result.get("pass"),
                      "entries": len(result["entries"]), "outdir": str(args.outdir)}, sort_keys=True))
    if args.mode == "smoke" and not result["pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
