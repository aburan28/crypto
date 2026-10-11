#!/usr/bin/env python3
"""Run the frozen one-target n=73 IC/rho pair and record per-process RSS."""
import hashlib
import json
import os
import platform
import signal
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
BASE = ROOT / "experiments/koblitz-single-target-n73-20261002/base_n73_K600.jsonl"
TARGETS = HERE / "frozen/target_points.jsonl"
IC_BINARY = ROOT / "target/release/examples/koblitz_orbit_dlp_fast"
RHO_BINARY = ROOT / "target/release/examples/koblitz_rho_fixture"


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def run_arm(name, argv, env, run_dir, run_record):
    stdout_path = run_dir / (name + ("_summary.jsonl" if name == "ic" else ".jsonl"))
    stderr_path = run_dir / (name + ".stderr.txt")
    started_utc = utc_now()
    started_monotonic = time.monotonic()
    receipt = {
        "arm": name,
        "status": "RUNNING",
        "argv": [str(value) for value in argv],
        "started_at_utc": started_utc,
        "memory_limit_mechanism": "No finite kernel cap; wait4 records per-process peak RSS",
        "memory_limit_bytes": protocol["resource_cap_bytes_per_arm"],
        "known_answer_environment_variable_present": "KIC_RHO_FIXED_TARGET_SCALAR" in env,
        "pid": None,
    }
    write_json(run_dir / (name + "_execution.json"), receipt)
    run_record["status"] = "RUNNING_" + name.upper()
    run_record["active_arm"] = name
    run_record["active_arm_started_at_utc"] = started_utc
    write_json(run_dir / "run.json", run_record)

    stdout_fd = os.open(stdout_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    stderr_fd = os.open(stderr_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    try:
        pid = os.fork()
        if pid == 0:
            try:
                os.dup2(stdout_fd, 1)
                os.dup2(stderr_fd, 2)
                os.close(stdout_fd)
                os.close(stderr_fd)
                os.chdir(ROOT)
                os.execve(str(argv[0]), [str(value) for value in argv], env)
            except BaseException as exc:
                os.write(2, ("runner child exec failure: " + repr(exc) + "\n").encode())
                os._exit(127)
        receipt["pid"] = pid
        write_json(run_dir / (name + "_execution.json"), receipt)
        last_heartbeat = time.monotonic()
        while True:
            waited_pid, wait_status, usage = os.wait4(pid, os.WNOHANG)
            if waited_pid == pid:
                return_code = os.waitstatus_to_exitcode(wait_status)
                peak_rss_bytes = int(usage.ru_maxrss)
                break
            time.sleep(0.5)
            if time.monotonic() - last_heartbeat >= 30:
                receipt["heartbeat_at_utc"] = utc_now()
                receipt["elapsed_process_wall_s"] = time.monotonic() - started_monotonic
                write_json(run_dir / (name + "_execution.json"), receipt)
                last_heartbeat = time.monotonic()
    except KeyboardInterrupt:
        try:
            os.kill(pid, signal.SIGTERM)
        except (NameError, ProcessLookupError):
            pass
        raise
    finally:
        os.close(stdout_fd)
        os.close(stderr_fd)

    receipt.update({
        "status": "COMPLETED" if return_code == 0 else "FAILED",
        "finished_at_utc": utc_now(),
        "elapsed_process_wall_s": time.monotonic() - started_monotonic,
        "return_code": return_code,
        "stdout_sha256": sha256(stdout_path),
        "stderr_sha256": sha256(stderr_path),
        "peak_rss_bytes": peak_rss_bytes,
        "peak_rss_source": "wait4(2) ru_maxrss for this exact child process (bytes on macOS)",
    })
    write_json(run_dir / (name + "_execution.json"), receipt)
    run_record[name + "_return_code"] = return_code
    run_record[name + "_finished_at_utc"] = receipt["finished_at_utc"]
    run_record["active_arm"] = None
    run_record["status"] = "RHO_COMPLETE_IC_PENDING" if name == "rho" else "BOTH_ARMS_COMPLETE"
    write_json(run_dir / "run.json", run_record)
    return return_code


protocol = json.loads((HERE / "protocol.json").read_text())
run_id = protocol["run_id"]
run_dir = RUNS / run_id
run_record = json.loads((run_dir / "run.json").read_text())
if protocol["target_count"] != 1 or run_record["target_count"] != 1:
    raise SystemExit("refusing to launch anything except the frozen one-target workload")
if protocol["paired_run_order"] != ["rho", "ic"]:
    raise SystemExit("protocol run order differs from frozen pair")
if len([line for line in TARGETS.read_text().splitlines() if line.strip()]) != 1:
    raise SystemExit("frozen input does not contain exactly one public point")
expected_ic = json.loads((HERE / "candidate-manifest.json").read_text())["record"]["implementation"]["component_sha256"]["compact_orbit_ic_binary"]
expected_rho = protocol["rho_reference"]["binary_sha256"]
if sha256(IC_BINARY) != expected_ic or sha256(RHO_BINARY) != expected_rho:
    raise SystemExit("binary hash no longer matches the frozen candidate/workload")
if run_record["status"] != "FROZEN_BEFORE_RUN":
    raise SystemExit("run is not in the pre-run frozen state; refusing duplicate execution")

target = json.loads(TARGETS.read_text().splitlines()[0])
rho_env = os.environ.copy()
for key in ("KIC_RHO_FIXED_TARGET_SCALAR", "KIC_RHO_BATCH_CORPUS", "KIC_RHO_FIXTURE_OFFSET"):
    rho_env.pop(key, None)
rho_env["KIC_RHO_PUBLIC_TARGET_POINT"] = json.dumps(target, separators=(",", ":"))
rho_env["KIC_RHO_WALK_SEED"] = str(protocol["rho_walk_seed"])
ic_env = os.environ.copy()
for key in ("KIC_RHO_FIXED_TARGET_SCALAR", "KIC_RHO_PUBLIC_TARGET_POINT", "KIC_RHO_BATCH_CORPUS"):
    ic_env.pop(key, None)

run_record.update({
    "platform": platform.platform(),
    "python": sys.version.split()[0],
    "launch_started_at_utc": utc_now(),
    "paired_run_order": protocol["paired_run_order"],
    "target_count": 1,
    "known_answer_sent_to_ic": False,
    "known_answer_sent_to_rho": False,
    "sidecar_validation_after_run": True,
    "memory_limit_mechanism": "No finite kernel cap; per-process peak RSS recorded with wait4(2)",
})
write_json(run_dir / "run.json", run_record)

rho_argv = [RHO_BINARY, "73", "0", "signed_frobenius", "1", "packed"]
ic_argv = [IC_BINARY, BASE, TARGETS, str(protocol["ic_rank_seed"]), run_dir / "ic.jsonl"]
rho_status = run_arm("rho", rho_argv, rho_env, run_dir, run_record)
ic_status = run_arm("ic", ic_argv, ic_env, run_dir, run_record)
run_record["launch_finished_at_utc"] = utc_now()
run_record["status"] = "PRODUCERS_COMPLETE" if rho_status == 0 and ic_status == 0 else "PRODUCER_FAILURE"
write_json(run_dir / "run.json", run_record)
print(json.dumps({
    "run_id": run_id,
    "status": run_record["status"],
    "rho_return_code": rho_status,
    "ic_return_code": ic_status,
    "run_dir": str(run_dir),
}, sort_keys=True))
if rho_status != 0 or ic_status != 0:
    raise SystemExit(1)
