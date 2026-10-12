#!/usr/bin/env python3
"""Run one frozen compact/rho child; retain raw failures and process costs."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import psutil
import subprocess
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SPEC = HERE / "input_spec.json"
GiB = 1024**3
SOURCES = {
    "train": "examples/koblitz_s5_sat_instance.rs",
    "compact": "examples/koblitz_s5_sat_instance.rs",
    "rho": "examples/koblitz_rho_batch_ks.rs",
}
EXECUTABLES = {
    "train": "target/release/examples/koblitz_s5_sat_instance",
    "compact": "target/release/examples/koblitz_s5_sat_instance",
    "rho": "target/release/examples/koblitz_rho_batch_ks",
}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def monitor_process(process: subprocess.Popen, started_ns: int,
                    timeout_s: int, rss_cap_bytes: int):
    """Poll child RSS on macOS/Linux and retain Darwin/Linux wait4 peak."""
    child = psutil.Process(process.pid)
    sampled_peak = 0
    samples = 0
    termination = None
    while True:
        pid, status, usage = os.wait4(process.pid, os.WNOHANG)
        if pid:
            break
        elapsed_s = (time.monotonic_ns()-started_ns)/1e9
        try:
            current_rss = child.memory_info().rss
            sampled_peak = max(sampled_peak, current_rss)
            samples += 1
        except (psutil.NoSuchProcess, psutil.ZombieProcess):
            current_rss = None
        if elapsed_s >= timeout_s or (current_rss is not None and current_rss >= rss_cap_bytes):
            termination = "TIMEOUT" if elapsed_s >= timeout_s else "RSS_CAP"
            process.kill()
            pid, status, usage = os.wait4(process.pid, 0)
            assert pid == process.pid
            break
        time.sleep(0.05)
    wait4_peak = int(usage.ru_maxrss * (1 if platform.system() == "Darwin" else 1024))
    return status, usage, termination, sampled_peak, wait4_peak, samples


def run_one(mode: str, n: int, block: int | None, length: int | None,
            out: Path, timeout_s: int, rss_cap_bytes: int = 2*GiB) -> dict:
    assert mode in SOURCES and n in (37, 41)
    assert timeout_s in (180, 600, 1800)
    if mode == "train":
        assert block is None and length is None and timeout_s == 1800
    else:
        assert block in (0, 1, 2) and length in (8, 32)
        assert timeout_s == (180 if length == 8 else 600)
    assert rss_cap_bytes == 2*GiB
    spec_bytes = SPEC.read_bytes()
    spec = json.loads(spec_bytes)
    arm = spec["arms"][str(n)]
    source_path = SOURCES[mode]
    source = (REPO / source_path).read_bytes()
    assert sha(source) == spec["source_sha256"][source_path]
    assert subprocess.run(["git", "merge-base", "--is-ancestor", spec["base_commit"], "HEAD"],
                          cwd=REPO).returncode == 0
    source_commit = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                            cwd=REPO, text=True).strip()
    assert not subprocess.check_output(["git", "status", "--porcelain"],
                                       cwd=REPO, text=True).strip(), "run from clean checkout"
    executable = REPO / EXECUTABLES[mode]
    assert executable.is_file(), f"build pinned executable first: {executable}"
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    env = {key: value for key, value in os.environ.items() if not key.startswith("KIC_")}
    env["RAYON_NUM_THREADS"] = "1"
    input_name = "target_scalars.txt" if mode == "train" else "target_points.jsonl"
    input_source = (HERE / arm["training_file"] if mode == "train" else
                    HERE / arm["blocks"][block]["files"][str(length)]["name"])
    input_data = input_source.read_bytes()
    expected_sha = (arm["training_sha256"] if mode == "train" else
                    arm["blocks"][block]["files"][str(length)]["sha256"])
    assert sha(input_data) == expected_sha
    input_file = out / input_name
    input_file.write_bytes(input_data)
    if mode == "rho":
        seed = 202609250000 + 100*n + block
        command = [str(executable), str(n), "0", "signed_frobenius",
                   str(length), str(seed)]
        env.update({"KIC_RHO_TARGET_POINTS_JSONL": str(input_file),
                    "KIC_RHO_DP_BITS": "4", "KIC_RHO_PRECOMPUTE_WALKS": "0",
                    "KIC_RHO_BATCH_CORPUS": f"sparse-shared-n{n}-b{block}-L{length}"})
    else:
        command = [str(executable), str(n), "0", str(arm["eta_numerator"]),
                   "1000000", "natural", "1", "2000", "1", "internal"]
        env.update({"KIC_ALGEBRA_ENCODING": "orbit_factorized",
                    "KIC_ORBIT_LAZY_RELATIVE_SUPPORT": "1",
                    "KIC_ORBIT_BRANCH_ORDER": "pair_then_pair",
                    "KIC_ORBIT_REP_ENCODING": "one_hot", "KIC_ORBIT_BATCH_ONLY": "1",
                    "KIC_ORBIT_INCLUDE_BASE_HEADER": "1"})
        env["KIC_ORBIT_TARGET_SCALARS" if mode == "train" else
            "KIC_ORBIT_TARGET_POINTS_JSONL"] = str(input_file)
    manifest = {"schema_version": "1.0", "mode": mode, "n": n, "block": block,
                "length": length, "R": arm["R"], "F": arm["F"],
                "input_spec_sha256": sha(spec_bytes), "input_file": input_name,
                "input_sha256": sha(input_data), "source_commit": source_commit,
                "source_path": source_path, "source_sha256": sha(source),
                "executable_sha256": sha(executable.read_bytes()),
                "command": command,
                "environment": {key: env[key] for key in sorted(env)
                                if key.startswith("KIC_") or key == "RAYON_NUM_THREADS"},
                "timeout_seconds": timeout_s, "rss_cap_bytes": rss_cap_bytes,
                "host": platform.node(), "platform": platform.platform()}
    write_json(out / "manifest.json", manifest)
    stdout, stderr = out / "producer.stdout.jsonl", out / "producer.stderr.txt"
    started = time.monotonic_ns()
    with stdout.open("wb") as out_stream, stderr.open("wb") as err_stream:
        process = subprocess.Popen(command, cwd=REPO, env=env,
                                   stdout=out_stream, stderr=err_stream)
        status, usage, termination, sampled_peak, wait4_peak, samples = monitor_process(
            process, started, timeout_s, rss_cap_bytes)
    wall_ms = (time.monotonic_ns()-started)/1e6
    process.returncode = os.waitstatus_to_exitcode(status)
    peak_rss = max(sampled_peak, wait4_peak)
    receipt = {"schema_version": "1.0", "returncode": process.returncode,
               "termination": termination, "timed_out": termination == "TIMEOUT",
               "rss_capped": termination == "RSS_CAP", "wall_ms": wall_ms,
               "user_cpu_s": usage.ru_utime, "system_cpu_s": usage.ru_stime,
               "peak_rss_bytes": peak_rss,
               "sampled_peak_rss_bytes": sampled_peak,
               "wait4_peak_rss_bytes": wait4_peak,
               "rss_sample_count": samples,
               "rss_monitor": "psutil.Process.memory_info().rss; 50ms poll",
               "psutil_version": psutil.__version__,
               "stdout_sha256": sha(stdout.read_bytes()),
               "stderr_sha256": sha(stderr.read_bytes()),
               "manifest_sha256": sha((out / "manifest.json").read_bytes())}
    write_json(out / "receipt.json", receipt)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=tuple(SOURCES), required=True)
    parser.add_argument("--n", type=int, choices=(37, 41), required=True)
    parser.add_argument("--block", type=int, choices=(0, 1, 2))
    parser.add_argument("--length", type=int, choices=(8, 32))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--timeout", type=int, required=True)
    args = parser.parse_args()
    print(json.dumps(run_one(args.mode, args.n, args.block, args.length,
                             args.out, args.timeout), sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
