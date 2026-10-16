#!/usr/bin/env python3
"""Run the three frozen n=83 a=1 one-target IC/rho pairs sequentially.

Mirrors the n=73 frozen contract: same retained base, same frozen public
target, same pinned binaries; rho then IC per run, per-process peak RSS
via wait4(2).  IC precompute uses the parallel guided rank
(KIC_RANK_THREADS=12); the online target stage is unchanged.
"""
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
BASE = HERE / "base_n83_K600.jsonl"
TARGETS = HERE / "frozen/target_points.jsonl"
IC_BINARY = ROOT / "target/release/examples/koblitz_orbit_dlp_fast"
RHO_BINARY = ROOT / "target/release/examples/koblitz_rho_fixture"
RANK_THREADS = 12

SEEDS = {
    "R1": {"rho_walk_seed": 202610061, "ic_rank_seed": 21},
    "R2": {"rho_walk_seed": 202610062, "ic_rank_seed": 22},
    "R3": {"rho_walk_seed": 202610063, "ic_rank_seed": 23},
}


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


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
        "memory_limit_bytes": None,
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


def main():
    RUNS.mkdir(exist_ok=True)
    target = json.loads(TARGETS.read_text().splitlines()[0])
    if len([line for line in TARGETS.read_text().splitlines() if line.strip()]) != 1:
        raise SystemExit("frozen input does not contain exactly one public point")
    fixture = json.loads((HERE / "frozen/fixture.json").read_text())

    for tag, seeds in SEEDS.items():
        run_id = f"N83A1K600W{RANK_THREADS}We{seeds['rho_walk_seed']}R{seeds['ic_rank_seed']}{tag}"
        run_dir = RUNS / run_id
        if run_dir.exists():
            raise SystemExit(f"run directory {run_dir} already exists; refusing duplicate execution")
        run_dir.mkdir(parents=True)
        run_record = {
            "active_arm": None,
            "candidate_id": "IC1N83A1Ckb1fb99600PDP4rootRCguidedLAgaussTDdirectISO0parallel",
            "curve": {
                "n": 83,
                "a": 1,
                "b": 1,
                "field_modulus": "x^83 + x^7 + x^4 + x^2 + 1",
                "subgroup_order": fixture["subgroup_order"],
                "cofactor": fixture["cofactor"],
            },
            "ic_binary_sha256": sha256(IC_BINARY),
            "ic_rank_seed": seeds["ic_rank_seed"],
            "ic_rank_threads": RANK_THREADS,
            "ic_return_code": None,
            "known_answer_sent_to_ic": False,
            "known_answer_sent_to_rho": False,
            "launch_started_at_utc": utc_now(),
            "memory_limit_mechanism": "No finite kernel cap; per-process peak RSS recorded with wait4(2)",
            "paired_run_order": ["rho", "ic"],
            "platform": platform.platform(),
            "python": sys.version.split()[0],
            "resource_cap_bytes_per_arm": None,
            "rho_binary_sha256": sha256(RHO_BINARY),
            "rho_return_code": None,
            "rho_walk_seed": seeds["rho_walk_seed"],
            "run_id": run_id,
            "sidecar_validation_after_run": True,
            "status": "FROZEN_BEFORE_RUN",
            "target": [str(v) for v in target],
            "target_count": 1,
            "validation_scalar_sidecar_path": "frozen/fixture.json",
        }
        write_json(run_dir / "run.json", run_record)

        rho_env = os.environ.copy()
        for key in ("KIC_RHO_FIXED_TARGET_SCALAR", "KIC_RHO_BATCH_CORPUS", "KIC_RHO_FIXTURE_OFFSET"):
            rho_env.pop(key, None)
        rho_env["KIC_RHO_PUBLIC_TARGET_POINT"] = json.dumps(target, separators=(",", ":"))
        rho_env["KIC_RHO_WALK_SEED"] = str(seeds["rho_walk_seed"])
        ic_env = os.environ.copy()
        for key in ("KIC_RHO_FIXED_TARGET_SCALAR", "KIC_RHO_PUBLIC_TARGET_POINT", "KIC_RHO_BATCH_CORPUS"):
            ic_env.pop(key, None)
        ic_env["KIC_RANK_THREADS"] = str(RANK_THREADS)

        rho_argv = [RHO_BINARY, "83", "1", "signed_frobenius", "1", "packed"]
        ic_argv = [IC_BINARY, BASE, TARGETS, str(seeds["ic_rank_seed"]), run_dir / "ic.jsonl"]
        rho_status = run_arm("rho", rho_argv, rho_env, run_dir, run_record)
        ic_status = run_arm("ic", ic_argv, ic_env, run_dir, run_record)
        run_record["launch_finished_at_utc"] = utc_now()
        run_record["status"] = (
            "PRODUCERS_COMPLETE" if rho_status == 0 and ic_status == 0 else "PRODUCER_FAILURE"
        )
        write_json(run_dir / "run.json", run_record)
        print(json.dumps({
            "run_id": run_id,
            "status": run_record["status"],
            "rho_return_code": rho_status,
            "ic_return_code": ic_status,
        }, sort_keys=True))
        if rho_status != 0 or ic_status != 0:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
