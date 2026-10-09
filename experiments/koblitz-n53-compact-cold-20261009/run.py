#!/usr/bin/env python3
"""Freeze executables, then capture one fail-closed n53 IC/rho pair.

`freeze` prints a JSON record to stdout. Commit that record as freeze.json
before using `run`; the latter refuses an uncommitted or modified freeze.
"""

import argparse
import hashlib
import json
import os
import platform
import signal
import subprocess
import sys
import time
from pathlib import Path

from verify_workload import W, scale

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE_FILES = (
    "Cargo.lock",
    "examples/koblitz_orbit_dlp_fast_online.rs",
    "examples/koblitz_rho_fixture.rs",
    "experiments/koblitz-n53-compact-cold-20261009/PROTOCOL.md",
    "experiments/koblitz-n53-compact-cold-20261009/run.py",
    "experiments/koblitz-n53-compact-cold-20261009/target_points.jsonl",
    "experiments/koblitz-n53-compact-cold-20261009/verify_workload.py",
    "experiments/koblitz-n53-compact-cold-20261009/workload.json",
)
CAP_SECONDS = 900.0
CAP_RSS_BYTES = 16 * 1024**3


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git(*args):
    return subprocess.check_output(("git", *args), cwd=ROOT).decode().strip()


def assert_clean_sources():
    changed = git("status", "--porcelain", "--untracked-files=no")
    if changed:
        raise RuntimeError("tracked source is dirty; commit or discard changes first")


def freeze(ic_binary, rho_binary):
    assert_clean_sources()
    for binary in (ic_binary, rho_binary):
        if not Path(binary).is_file() or not os.access(binary, os.X_OK):
            raise RuntimeError(f"executable missing: {binary}")
    record = {
        "schema_version": 1,
        "kind": "n53_compact_cold_source_freeze",
        "source_commit": git("rev-parse", "HEAD"),
        "source_sha256": {name: sha(ROOT / name) for name in SOURCE_FILES},
        "ic_binary_sha256": sha(ic_binary),
        "rho_binary_sha256": sha(rho_binary),
        "ic_argv_suffix": ["construct:53:0:244", "target_points.jsonl", "20261009"],
        "rho_argv_suffix": ["53", "0", "signed_frobenius", "1", "strong", "20261009", "7948768810114"],
        "wall_cap_seconds_per_arm": CAP_SECONDS,
        "rss_cap_bytes_per_arm": CAP_RSS_BYTES,
    }
    print(json.dumps(record, sort_keys=True, indent=2))


def verify_freeze(path, ic_binary, rho_binary):
    path = Path(path).resolve()
    try:
        relative = path.relative_to(ROOT).as_posix()
    except ValueError as exc:
        raise RuntimeError("freeze.json must be inside the source repository") from exc
    committed = subprocess.check_output(("git", "show", f"HEAD:{relative}"), cwd=ROOT)
    if committed != path.read_bytes():
        raise RuntimeError("freeze.json differs from the committed version")
    record = json.loads(committed)
    if record.get("kind") != "n53_compact_cold_source_freeze" or record.get("schema_version") != 1:
        raise RuntimeError("unknown freeze schema")
    assert_clean_sources()
    subprocess.run(("git", "merge-base", "--is-ancestor", record["source_commit"], "HEAD"), cwd=ROOT, check=True)
    if record["source_sha256"] != {name: sha(ROOT / name) for name in SOURCE_FILES}:
        raise RuntimeError("source/input hash differs from freeze")
    if record["ic_binary_sha256"] != sha(ic_binary) or record["rho_binary_sha256"] != sha(rho_binary):
        raise RuntimeError("executable hash differs from freeze")
    if record["wall_cap_seconds_per_arm"] != CAP_SECONDS or record["rss_cap_bytes_per_arm"] != CAP_RSS_BYTES:
        raise RuntimeError("resource cap differs from protocol")
    return record


def rss_bytes(pid):
    # On Darwin and Linux, ps reports resident size in KiB. An unavailable
    # monitor is a terminal attempt failure, never an unbounded run.
    result = subprocess.run(("/bin/ps", "-o", "rss=", "-p", str(pid)), text=True,
                            capture_output=True, timeout=2)
    if result.returncode or not result.stdout.strip():
        raise OSError(f"RSS monitor failed: {result.returncode}: {result.stderr.strip()}")
    return int(result.stdout.strip()) * 1024


def terminate_group(proc):
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def run_arm(argv, stdout_path, stderr_path, env):
    started = time.monotonic_ns()
    peak_sampled = 0
    terminal = "exited"
    monitor_error = None
    with open(stdout_path, "wb") as out, open(stderr_path, "wb") as err:
        try:
            proc = subprocess.Popen(argv, stdout=out, stderr=err, env=env, start_new_session=True)
        except OSError as exc:
            return {"argv": [str(part) for part in argv], "status": "launch_failure",
                    "exit_code": None, "wall_ms": (time.monotonic_ns() - started) / 1e6,
                    "sampled_peak_rss_bytes": 0, "child_peak_rss_bytes": None,
                    "monitor_error": str(exc), "stdout_sha256": sha(stdout_path),
                    "stderr_sha256": sha(stderr_path)}
        while True:
            child, status, usage = os.wait4(proc.pid, os.WNOHANG)
            if child:
                break
            elapsed = (time.monotonic_ns() - started) / 1e9
            if elapsed >= CAP_SECONDS:
                terminal = "timeout"
                terminate_group(proc)
                child, status, usage = os.wait4(proc.pid, 0)
                break
            try:
                peak_sampled = max(peak_sampled, rss_bytes(proc.pid))
            except (OSError, subprocess.TimeoutExpired, ValueError) as exc:
                # A fast child may exit between wait4(WNOHANG) and ps.
                child, status, usage = os.wait4(proc.pid, os.WNOHANG)
                if child:
                    break
                terminal = "monitor_failure"
                monitor_error = str(exc)
                terminate_group(proc)
                child, status, usage = os.wait4(proc.pid, 0)
                break
            if peak_sampled > CAP_RSS_BYTES:
                terminal = "memory_cap"
                terminate_group(proc)
                child, status, usage = os.wait4(proc.pid, 0)
                break
            time.sleep(0.1)
    stopped = time.monotonic_ns()
    exit_code = os.waitstatus_to_exitcode(status)
    proc.returncode = exit_code
    peak_child = int(usage.ru_maxrss) * (1 if sys.platform == "darwin" else 1024)
    if terminal == "exited" and peak_child > CAP_RSS_BYTES:
        terminal = "memory_cap"
    return {
        "argv": [str(part) for part in argv],
        "status": terminal if terminal != "exited" else ("success" if exit_code == 0 else "process_failure"),
        "exit_code": exit_code,
        "wall_ms": (stopped - started) / 1e6,
        "sampled_peak_rss_bytes": peak_sampled,
        "child_peak_rss_bytes": peak_child,
        "monitor_error": monitor_error,
        "stdout_sha256": sha(stdout_path),
        "stderr_sha256": sha(stderr_path),
    }


def one_json_line(path, kind):
    rows = [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]
    matches = [row for row in rows if row.get("kind") == kind]
    if len(matches) != 1:
        raise ValueError(f"expected one {kind} record in {path}, got {len(matches)}")
    return matches[0]


def analyse(ic, rho, run_dir):
    target = tuple(W["primary_target"])
    scalar = W["verification_scalar"]
    result = {"status": "unverified", "ic_online_ms": None, "rho_online_ms": None,
              "online_speedup": None, "independent_scalar_replay": False}
    if ic["status"] != "success" or rho["status"] != "success":
        result["status"] = "arm_failure"
        return result
    try:
        ic_target = one_json_line(run_dir / "ic-target.jsonl", "compact_orbit_dlp_target")
        ic_summary = one_json_line(run_dir / "ic.stdout.jsonl", "compact_orbit_dlp_summary")
        rho_row = one_json_line(run_dir / "rho.stdout.jsonl", "rho_public_fixture")
        phases = [ic_target[name] for name in (
            "target_query_ms", "target_pdp_ms", "target_relation_check_ms",
            "target_descent_ms", "target_recovery_check_ms")]
        online = float(ic_target["online_ms"])
        if abs(sum(phases) - online) > max(1e-6, online * 1e-6):
            raise ValueError("IC exclusive target phases do not sum to online interval")
        if abs(float(ic_target["target_phase_sum_ms"]) - online) > max(1e-6, online * 1e-6):
            raise ValueError("producer target phase sum differs from online interval")
        cold = float(ic_summary["cold_in_process_ms"])
        cold_phases = ic_summary["cold_phase_ms"]
        if cold_phases is None or abs(sum(float(value) for value in cold_phases.values()) - cold) > max(1e-6, cold * 1e-6):
            raise ValueError("IC cold phases do not sum to cold interval")
        if tuple(ic_target["published_q"]) != target or tuple(rho_row["published_q"]) != target:
            raise ValueError("IC/rho point mismatch")
        if int(ic_target["recovered_scalar"]) != scalar or int(rho_row["recovered_fixture_scalar"]) != scalar:
            raise ValueError("IC/rho scalar mismatch")
        if ic_target["group_verified"] is not True or rho_row["verified"] is not True:
            raise ValueError("producer scalar verification failed")
        if scale(scalar, tuple(W["generator"])) != target:
            raise ValueError("independent scalar replay failed")
        if ic_summary["rank"] != 244 or ic_summary["orbit_columns"] != 244:
            raise ValueError("IC relation matrix did not reach full rank")
        if ic_summary["factor_base_points"] != 25864:
            raise ValueError("actual factor base differs from protocol")
        if rho_row["reference_grade"] != "strong" or rho_row["quotient_mode"] != "signed_frobenius":
            raise ValueError("rho reference mismatch")
        result.update({"status": "pending_independent_relation_replay",
                       "ic_online_ms": online,
                       "rho_online_ms": float(rho_row["walk_ms"]) + float(rho_row["validation_ms"]),
                       "independent_scalar_replay": True,
                       "ic_cold_in_process_ms": cold,
                       "ic_cold_phase_ms": cold_phases,
                       "actual_factor_base_points": ic_summary["factor_base_points"],
                       "orbit_columns": ic_summary["orbit_columns"],
                       "rank_attempts": ic_summary["rank_attempts"],
                       "rank_relations": ic_summary["rank_relations"],
                       "rank_failures": ic_summary["rank_failures"],
                       "rho_walk_steps": rho_row["walk_steps"],
                       "rho_automorphism_size": rho_row["automorphism_size"]})
    except (ValueError, TypeError, KeyError, IndexError) as exc:
        result["validation_error"] = str(exc)
    return result


def execute(freeze_path, ic_binary, rho_binary, run_dir):
    record = verify_freeze(freeze_path, ic_binary, rho_binary)
    run_dir = Path(run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=False)
    started = {"kind": "n53_compact_cold_attempt_started", "source_freeze_sha256": sha(freeze_path),
               "source_commit": record["source_commit"], "workload_sha256": record["source_sha256"][SOURCE_FILES[-1]],
               "host": {"platform": platform.platform(), "processor": platform.processor(), "logical_cpus": os.cpu_count()},
               "environment_kic": {k: v for k, v in os.environ.items() if k.startswith("KIC_")}}
    (run_dir / "attempt_started.json").write_text(json.dumps(started, sort_keys=True, indent=2) + "\n")
    preflight_error = None
    if started["environment_kic"]:
        preflight_error = "unexpected ambient KIC_* variables"
    else:
        try:
            rss_bytes(os.getpid())
        except (OSError, subprocess.TimeoutExpired, ValueError) as exc:
            preflight_error = f"RSS monitor unavailable: {exc}"
    if preflight_error is not None:
        failure = {"kind": "n53_compact_cold_preflight_failure",
                   "reason": preflight_error, "attempt_started_sha256": sha(run_dir / "attempt_started.json")}
        (run_dir / "receipt.json").write_text(json.dumps(failure, sort_keys=True, indent=2) + "\n")
        print(json.dumps({"status": "preflight_failure", "receipt": str(run_dir / "receipt.json")}, sort_keys=True))
        raise SystemExit(2)
    ic_args = [str(Path(ic_binary).resolve()), "construct:53:0:244", str(HERE / "target_points.jsonl"),
               "20261009", str(run_dir / "ic-target.jsonl")]
    ic_env = os.environ.copy()
    ic_env["KIC_DUMP_BASE"] = str(run_dir / "base.jsonl")
    ic_env["KIC_DUMP_RANK"] = str(run_dir / "rank.jsonl")
    ic = run_arm(ic_args, run_dir / "ic.stdout.jsonl", run_dir / "ic.stderr.txt", ic_env)
    rho_args = [str(Path(rho_binary).resolve()), "53", "0", "signed_frobenius", "1", "strong",
                "20261009", str(W["verification_scalar"])]
    rho = run_arm(rho_args, run_dir / "rho.stdout.jsonl", run_dir / "rho.stderr.txt", os.environ.copy())
    files = {name: sha(run_dir / name) for name in (
        "attempt_started.json", "ic.stdout.jsonl", "ic.stderr.txt", "rho.stdout.jsonl", "rho.stderr.txt")}
    for optional in ("ic-target.jsonl", "base.jsonl", "rank.jsonl"):
        if (run_dir / optional).exists():
            files[optional] = sha(run_dir / optional)
    receipt = {"kind": "n53_compact_cold_attempt", "ic": ic, "rho": rho,
               "analysis": analyse(ic, rho, run_dir), "files_sha256": files,
               "source_freeze_sha256": sha(freeze_path), "isolation_verified": False}
    (run_dir / "receipt.json").write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": receipt["analysis"]["status"], "receipt": str(run_dir / "receipt.json")}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser()
    modes = parser.add_subparsers(dest="mode", required=True)
    for mode in ("freeze", "run"):
        sub = modes.add_parser(mode)
        sub.add_argument("--ic-binary", required=True, type=Path)
        sub.add_argument("--rho-binary", required=True, type=Path)
        if mode == "run":
            sub.add_argument("--freeze", required=True, type=Path)
            sub.add_argument("--run-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.mode == "freeze":
        freeze(args.ic_binary, args.rho_binary)
    else:
        execute(args.freeze, args.ic_binary, args.rho_binary, args.run_dir)


if __name__ == "__main__":
    main()
