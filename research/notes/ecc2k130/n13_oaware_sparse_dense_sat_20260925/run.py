#!/usr/bin/env python3
"""Frozen ABBA-export / paired dense-sparse cold SAT stage (no adaptive queries)."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import resource
import shutil
import signal
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import psutil
import verify

HERE = Path(__file__).resolve().parent
SOLVERS = ("cryptominisat5", "kissat", "cadical")
EXPORT_ORDER = ("dense", "sparse", "sparse", "dense")
EXPORT_WALL, EXPORT_RSS = 180.0, 512 * 1024**2
SOLVER_WALL, SOLVER_RSS = 15.0, 2 * 1024**3
POLL = 0.02


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha_bytes(raw):
    return hashlib.sha256(raw).hexdigest()


def save(path, obj):
    path.write_text(json.dumps(obj, sort_keys=True, separators=(",", ":")) + "\n")


def usage():
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    return r.ru_utime, r.ru_stime, r.ru_maxrss


def tree_rss(proc):
    if proc is None:
        return 0
    try:
        family = [proc] + proc.children(recursive=True)
        return sum(item.memory_info().rss for item in family if item.is_running())
    except (psutil.Error, ProcessLookupError):
        return 0


def execute(command, out, err, wall_cap, rss_cap, deadline=None):
    before = usage()
    started = utc()
    clock = time.perf_counter()
    with out.open("wb") as stdout, err.open("wb") as stderr:
        child = subprocess.Popen(command, stdout=stdout, stderr=stderr, start_new_session=True)
        try:
            proc = psutil.Process(child.pid)
        except psutil.NoSuchProcess:
            proc = None
        peak, stop = 0, None
        while child.poll() is None:
            peak = max(peak, tree_rss(proc))
            elapsed = time.perf_counter() - clock
            if deadline is not None and time.perf_counter() > deadline:
                stop = "portfolio_cap"
            elif peak > rss_cap:
                stop = "rss_cap"
            elif elapsed > wall_cap:
                stop = "wall_cap"
            if stop:
                try:
                    os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                break
            time.sleep(POLL)
        exit_code = child.wait()
        peak = max(peak, tree_rss(proc))
    after = usage()
    scale = 1 if sys.platform == "darwin" else 1024
    return {"command": command, "start_utc": started, "end_utc": utc(),
            "wall_seconds": time.perf_counter() - clock,
            "exit_code": exit_code, "stop_reason": stop,
            "sampled_peak_tree_rss_bytes": peak,
            "children_ru_maxrss_before_bytes": before[2] * scale,
            "children_ru_maxrss_after_bytes": after[2] * scale,
            "user_cpu_seconds": after[0] - before[0],
            "system_cpu_seconds": after[1] - before[1],
            "stdout_sha256": verify.sha(out), "stderr_sha256": verify.sha(err),
            "stdout_bytes": out.stat().st_size, "stderr_bytes": err.stat().st_size}


def solver_command(name, binary, cnf):
    return ([str(binary), "--verb=0", "--threads=1", str(cnf)]
            if name == "cryptominisat5" else [str(binary), str(cnf)])


def preflight():
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert frozen["domain"] == "ecc2k130-n13-m5-paired-dense-sparse-cnf-v1"
    assert frozen["solver_order"] == list(SOLVERS)
    assert frozen["export_order"] == list(EXPORT_ORDER)
    assert frozen["export_wall_cap_seconds"] == EXPORT_WALL
    assert frozen["export_rss_cap_bytes"] == EXPORT_RSS
    assert frozen["solver_wall_cap_seconds"] == SOLVER_WALL
    assert frozen["solver_rss_cap_bytes"] == SOLVER_RSS
    assert platform.python_version() == frozen["python_version"]
    assert str(Path(sys.executable).resolve()) == frozen["python_executable"]
    assert psutil.__version__ == frozen["psutil_version"]
    pinned = {
        "protocol_sha256": HERE / "PROTOCOL.md",
        "run_sha256": HERE / "run.py",
        "verify_sha256": HERE / "verify.py",
        "ci_replay_sha256": HERE / "ci_replay.py",
        "analyze_sha256": HERE / "analyze.py",
        "dense_export_sha256": verify.DENSE / "export.py",
        "sparse_export_sha256": verify.SPARSE / "export.py",
        "dense_base_sha256": verify.DENSE_DIR / "base.cnf",
        "dense_schema_sha256": verify.DENSE_DIR / "schema.json",
        "sparse_base_sha256": verify.SPARSE_DIR / "base.cnf",
        "sparse_schema_sha256": verify.SPARSE_DIR / "schema.json",
        "sparse_verify_sha256": verify.SPARSE / "verify.py",
        "sparse_frozen_sha256": verify.SPARSE / "FROZEN.json",
        "sparse_evidence_receipt_sha256": verify.SPARSE / "evidence/final/receipt.json",
        "prior_verify_sha256": verify.OLD / "verify.py",
        "prior_input_sha256": verify.OLD / "INPUT.json",
        "prior_frozen_sha256": verify.OLD / "FROZEN.json",
    }
    for key, path in pinned.items():
        assert verify.sha(path) == frozen[key], (key, str(path))
    binaries = {}
    for name in SOLVERS:
        binary = Path(frozen["binaries"][name]["path"])
        assert binary.is_file() and verify.sha(binary) == frozen["binaries"][name]["sha256"]
        binaries[name] = binary
    schemas, truth, curve, point, old = verify.schemas_and_truth()
    for rep, source in (("dense", verify.DENSE_DIR), ("sparse", verify.SPARSE_DIR)):
        base = (source / "base.cnf").read_bytes()
        assert sha_bytes(base) == frozen[f"{rep}_base_sha256"]
        assert base.split(b"\n", 2)[1] == (
            f"p cnf {schemas[rep]['variables']} {schemas[rep]['clauses']}".encode())
    return frozen, binaries, schemas, truth, curve, point, old


def smoke(outdir, frozen, binaries, old):
    entries = []
    for name in SOLVERS:
        for label, raw, expected in (("sat", b"p cnf 1 1\n1 0\n", "SAT"),
                                     ("unsat", b"p cnf 1 2\n1 0\n-1 0\n", "UNSAT")):
            stem = f"{name}-{label}"
            cnf = outdir / f"{stem}.cnf"
            cnf.write_bytes(raw)
            stdout, stderr = outdir / f"{stem}.stdout", outdir / f"{stem}.stderr"
            row = execute(solver_command(name, binaries[name], cnf), stdout, stderr,
                          SOLVER_WALL, SOLVER_RSS)
            try:
                status, assignment = old.parse_solver_output(stdout.read_bytes(), row["exit_code"], 1)
                valid = (status == expected and (label != "sat" or assignment.get(1) is True))
            except (ValueError, AssertionError):
                valid = False
            row.update({"solver": name, "case": label, "input_sha256": sha_bytes(raw),
                        "pass": valid and row["stop_reason"] is None})
            save(outdir / f"{stem}.json", row)
            entries.append(row)
    result = {"mode": "smoke", "freeze_sha256": verify.sha(HERE / "FROZEN.json"),
              "entries": entries, "pass": all(row["pass"] for row in entries)}
    save(outdir / "result.json", result)
    return result


def export_abba(outdir, frozen, deadline):
    rows = []
    work = outdir / "work"
    work.mkdir()
    for index, rep in enumerate(EXPORT_ORDER):
        stem = f"export{index}-{rep}"
        output = work / stem
        command = [str(Path(sys.executable).resolve()),
                   str((verify.DENSE if rep == "dense" else verify.SPARSE) / "export.py"),
                   "--out", str(output)]
        stdout, stderr = outdir / f"{stem}.stdout", outdir / f"{stem}.stderr"
        row = execute(command, stdout, stderr, EXPORT_WALL, EXPORT_RSS, deadline)
        row.update({"representation": rep, "pair": 0 if index < 2 else 1,
                    "expected_base_sha256": frozen[f"{rep}_base_sha256"],
                    "expected_schema_sha256": frozen[f"{rep}_schema_sha256"]})
        if row["exit_code"] == 0 and row["stop_reason"] is None:
            source = output / "n13-m5"
            assert verify.sha(source / "base.cnf") == row["expected_base_sha256"]
            assert verify.sha(source / "schema.json") == row["expected_schema_sha256"]
            row["producer_result_sha256"] = verify.sha(output / "result.json")
            shutil.copyfile(output / "result.json", outdir / f"{stem}.result.json")
        save(outdir / f"{stem}.json", row)
        rows.append(row)
        if row["exit_code"] != 0 or row["stop_reason"] is not None:
            raise RuntimeError(f"{stem} failed/censored; partial outputs retained")
    return rows


def paired_panel(outdir, frozen, binaries, schemas, truth, curve, point, old, deadline):
    exports = export_abba(outdir, frozen, deadline)
    clock = time.perf_counter()
    base = {"dense": (verify.DENSE_DIR / "base.cnf").read_bytes(),
            "sparse": (verify.SPARSE_DIR / "base.cnf").read_bytes()}
    base_load_wall = time.perf_counter() - clock
    rows = []
    for ti in range(32):
        target = schemas["dense"]["targets"][ti]
        prepared = {}
        for rep in ("dense", "sparse"):
            qclock = time.perf_counter()
            raw = verify.query_bytes(old, base[rep], target, schemas[rep])
            path = outdir / f"query-{rep}.cnf"
            path.write_bytes(raw)
            prepared[rep] = {"path": path, "sha256": sha_bytes(raw),
                             "bytes": len(raw), "setup_wall_seconds": time.perf_counter() - qclock}
        engines = SOLVERS[ti % 3:] + SOLVERS[:ti % 3]
        reps = ("dense", "sparse") if ti % 2 == 0 else ("sparse", "dense")
        for name in engines:
            for rep in reps:
                stem = f"{target['id']}-{name}-{rep}"
                stdout, stderr = outdir / f"{stem}.stdout", outdir / f"{stem}.stderr"
                row = execute(solver_command(name, binaries[name], prepared[rep]["path"]),
                              stdout, stderr, SOLVER_WALL, SOLVER_RSS, deadline)
                proof_clock = time.perf_counter()
                if row["stop_reason"]:
                    verdict, certificate, error = "CENSORED", None, row["stop_reason"]
                else:
                    try:
                        verdict, certificate = verify.classify(
                            old, stdout.read_bytes(), row["exit_code"], schemas[rep],
                            target, truth, curve, point)
                        error = None
                    except (ValueError, AssertionError) as exc:
                        verdict, certificate, error = "INVALID", None, f"{type(exc).__name__}: {exc}"
                row.update({"id": target["id"], "solver": name, "representation": rep,
                            "assumption_literal": target["assumption_literal"],
                            "input_sha256": prepared[rep]["sha256"],
                            "input_bytes": prepared[rep]["bytes"],
                            "query_setup_wall_seconds": prepared[rep]["setup_wall_seconds"],
                            "verify_wall_seconds": time.perf_counter() - proof_clock,
                            "point_oracle_positive": truth[target["id"]],
                            "verdict": verdict, "certificate": certificate, "error": error})
                save(outdir / f"{stem}.json", row)
                rows.append(row)
                if row["stop_reason"] == "portfolio_cap":
                    raise RuntimeError("3900-second aggregate cap; partial raw panel retained")
        for item in prepared.values():
            item["path"].unlink()
    shutil.rmtree(outdir / "work")
    return exports, rows, base_load_wall


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("smoke", "panel"))
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--smoke-receipt", type=Path)
    args = parser.parse_args()
    assert not args.outdir.exists(), "never overwrite an attempt"
    process_start = time.perf_counter()
    frozen, binaries, schemas, truth, curve, point, old = preflight()
    preflight_wall = time.perf_counter() - process_start
    if args.mode == "panel":
        assert args.smoke_receipt is not None
        prior = json.loads(args.smoke_receipt.read_text())
        assert prior["mode"] == "smoke" and prior["pass"]
        assert prior["freeze_sha256"] == verify.sha(HERE / "FROZEN.json")
    args.outdir.mkdir(parents=True)
    receipt = {"mode": args.mode, "freeze_sha256": verify.sha(HERE / "FROZEN.json"),
               "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
               "python": sys.version, "psutil": psutil.__version__, "host": platform.platform(),
               "start_utc": utc(), "preflight_wall_seconds": preflight_wall,
               "decision": "INCOMPLETE"}
    save(args.outdir / "receipt.json", receipt)
    try:
        if args.mode == "smoke":
            result = smoke(args.outdir, frozen, binaries, old)
            assert result["pass"]
        else:
            exports, entries, base_wall = paired_panel(
                args.outdir, frozen, binaries, schemas, truth, curve, point, old,
                process_start + 3900.0)
            result = {"mode": "panel", "freeze_sha256": receipt["freeze_sha256"],
                      "exports": exports, "entries": entries,
                      "base_load_wall_seconds": base_wall}
            save(args.outdir / "result.json", result)
        receipt["decision"] = "COMPLETE"
        receipt["end_utc"] = utc()
        receipt["process_wall_seconds"] = time.perf_counter() - process_start
        save(args.outdir / "receipt.json", receipt)
        print(json.dumps({"mode": args.mode, "decision": receipt["decision"],
                          "entries": len(result["entries"])}, sort_keys=True))
    except Exception as exc:
        receipt["decision"] = "CENSORED_OR_FAILED"
        receipt["error"] = repr(exc)
        receipt["traceback"] = traceback.format_exc()
        receipt["end_utc"] = utc()
        receipt["process_wall_seconds"] = time.perf_counter() - process_start
        save(args.outdir / "receipt.json", receipt)
        raise


if __name__ == "__main__":
    main()
