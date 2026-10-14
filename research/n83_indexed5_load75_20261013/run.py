#!/usr/bin/env python3
"""Bounded N83 indexed-five-byte dense-request single-target capacity gate."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / "research/n83_indexed5_20261013"
SOURCE = ROOT / "examples/koblitz_orbit_dlp_fast.rs"
RHO_SOURCE = ROOT / "examples/n83_rho_optimization_eval.rs"
LOCK = ROOT / "Cargo.lock"
POINT = STUDY / "primary_target.jsonl"
ANSWER = STUDY / "answer.jsonl"
RHO_CONFIG = ROOT / "research/n83_dense_root_20261013/rho_config.json"
EXPECTED = {
    SOURCE: "80e54cba2327bceece22901163bcf32e84c096acc1af795fd4b5d02f0542e287",
    LOCK: "ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90",
    POINT: "21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8",
    ANSWER: "651eb7d0ef08a1655036746311b19b6989e340ca4752b6d675606055e59c12a1",
    RHO_CONFIG: "cf76ebb5372acfc1355b88374f875980be66cdffba1a82cb6e4c2f28bcfa723b",
    RHO_SOURCE: "fb14072f4fa5d08e18d8b79130b6d58a220adae0f76369cc9be949691c83eaa0",
}
ARMS = ((3114, True, 2400), (3000, True, 2000),
        (2543, True, 1400), (2543, False, 1400))
RSS_CAP_BYTES = 14 * 1024**3


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def read_text(path):
    try:
        return path.read_text()
    except OSError:
        return None


def cgroup_root():
    for line in (read_text(Path("/proc/self/cgroup")) or "").splitlines():
        if line.startswith("0::"):
            relative = line.split("::", 1)[1].lstrip("/")
            path = Path("/sys/fs/cgroup") / relative
            if path.exists():
                return path
    return Path("/sys/fs/cgroup")


def cgroup_snapshot(group):
    return {name: read_text(group / name) for name in
            ("memory.current", "memory.peak", "memory.max", "memory.events",
             "cpu.stat", "cpu.max", "cpu.pressure", "memory.pressure")}


def host_receipt(group):
    cpuinfo = read_text(Path("/proc/cpuinfo")) or ""
    meminfo = read_text(Path("/proc/meminfo")) or ""
    total_kib = next((int(line.split()[1]) for line in meminfo.splitlines()
                      if line.startswith("MemTotal:")), None)
    model = next((line.split(":", 1)[1].strip() for line in cpuinfo.splitlines()
                  if line.startswith("model name")), None)
    return {
        "schema": "n83.indexed5-load75-tier1-host/v1", "recorded_utc": now(),
        "provider": "github_actions", "instance_id": "gh-{}-{}-{}".format(
            os.environ.get("GITHUB_RUN_ID", "local"),
            os.environ.get("GITHUB_RUN_ATTEMPT", "0"),
            os.environ.get("GITHUB_JOB", "local")),
        "run_url": "{}/{}/actions/runs/{}".format(
            os.environ.get("GITHUB_SERVER_URL", "https://github.com"),
            os.environ.get("GITHUB_REPOSITORY", "aburan28/crypto"),
            os.environ.get("GITHUB_RUN_ID", "local")),
        "git_sha": os.environ.get("GITHUB_SHA"),
        "runner_name": os.environ.get("RUNNER_NAME"),
        "runner_os": os.environ.get("RUNNER_OS"),
        "runner_arch": os.environ.get("RUNNER_ARCH"),
        "image_version": os.environ.get("ImageVersion"),
        "platform": platform.platform(), "uname": platform.uname()._asdict(),
        "cpu_model": model,
        "visible_affinity": (sorted(os.sched_getaffinity(0))
                             if hasattr(os, "sched_getaffinity") else None),
        "mem_total_kib": total_kib, "cgroup_path": str(group),
        "cgroup_before": cgroup_snapshot(group),
        "disk_free_bytes": shutil.disk_usage(ROOT).free,
        "python": sys.version,
    }


def check_inputs():
    actual = {str(path.relative_to(ROOT)): digest(path) for path in EXPECTED}
    bad = [str(path.relative_to(ROOT)) for path, expected in EXPECTED.items()
           if actual[str(path.relative_to(ROOT))] != expected]
    if bad:
        raise ValueError("frozen input/source hash mismatch: " + ", ".join(bad))
    return actual


def sampled_rss(pid):
    status = read_text(Path(f"/proc/{pid}/status")) or ""
    for line in status.splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) * 1024
    return None


def reap_bounded(child, cap_seconds):
    started = time.monotonic()
    sampled_peak = 0
    stop_reason = None
    status = None
    usage = None
    while True:
        pid, status, usage = os.wait4(child.pid, os.WNOHANG)
        if pid:
            break
        rss = sampled_rss(child.pid)
        if rss is not None:
            sampled_peak = max(sampled_peak, rss)
            if rss > RSS_CAP_BYTES:
                stop_reason = "STOPPED_RESOURCE_GUARD"
        if time.monotonic() - started >= cap_seconds:
            stop_reason = "TIMEOUT"
        if stop_reason:
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                try:
                    child.kill()
                except ProcessLookupError:
                    pass
            _, status, usage = os.wait4(child.pid, 0)
            break
        time.sleep(0.5)
    code = os.waitstatus_to_exitcode(status)
    child.returncode = code
    return {
        "exit_status": code, "stop_reason": stop_reason,
        "wall_seconds": time.monotonic() - started,
        "user_cpu_seconds": usage.ru_utime,
        "system_cpu_seconds": usage.ru_stime,
        "peak_rss_bytes": usage.ru_maxrss * 1024,
        "sampled_peak_rss_bytes": sampled_peak,
        "involuntary_context_switches": usage.ru_nivcsw,
        "major_faults": usage.ru_majflt,
    }


def one_run(folder, command, environment, cap_seconds, group):
    folder.mkdir(parents=True, exist_ok=False)
    write_json(folder / "launch.json", {
        "command": command, "environment_flags": {
            key: value for key, value in environment.items() if key.startswith("KIC_")},
        "cap_seconds": cap_seconds, "rss_cap_bytes": RSS_CAP_BYTES,
        "started_utc": now(),
    })
    before = cgroup_snapshot(group)
    with (folder / "stdout.jsonl").open("wb") as stdout, (folder / "stderr.log").open("wb") as stderr:
        child = subprocess.Popen(command, stdout=stdout, stderr=stderr,
                                 cwd=ROOT, env=environment, start_new_session=True)
        usage = reap_bounded(child, cap_seconds)
    after = cgroup_snapshot(group)
    usage["finished_utc"] = now()
    usage["cgroup_before"] = before
    usage["cgroup_after"] = after
    usage["status"] = (usage["stop_reason"] or
                       ("PASS_PRODUCER" if usage["exit_status"] == 0 else "FAILED"))
    write_json(folder / "usage.json", usage)
    (folder / "exit_status.txt").write_text(str(usage["exit_status"]) + "\n")
    return usage


def run(binary, rho_binary, output):
    hashes = check_inputs()
    if output.exists():
        raise ValueError("output directory already exists")
    output.mkdir(parents=True)
    group = cgroup_root()
    host = host_receipt(group)
    write_json(output / "host.json", host)
    write_json(output / "source.json", {
        "schema": "n83.indexed5-load75-tier1-source/v1", "hashes": hashes,
        "binary_sha256": digest(binary), "rho_binary_sha256": digest(rho_binary),
    })
    if (host["mem_total_kib"] is None or host["mem_total_kib"] < 15_000_000
            or host["disk_free_bytes"] < 2 * 1024**3
            or len(host["visible_affinity"] or []) < 4):
        write_json(output / "preflight.json", {"status": "REJECTED_HOST_CAPACITY",
                                               "mem_total_kib": host["mem_total_kib"],
                                               "disk_free_bytes": host["disk_free_bytes"],
                                               "visible_affinity": host["visible_affinity"]})
        return 2
    write_json(output / "preflight.json", {"status": "PASS_CAPACITY_PREFLIGHT",
                                           "mem_total_kib": host["mem_total_kib"],
                                           "visible_affinity": host["visible_affinity"]})
    for k, dense_request, cap in ARMS:
        arm = f"k{k}_{'dense' if dense_request else 'regular'}"
        print(f"starting {arm} at {now()}", flush=True)
        folder = output / "runs" / "primary" / arm
        env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin"),
               "KIC_PACKED_ROOT_TABLE": "1", "KIC_PACKED_ROOT_TABLE_DENSE": "0",
               "KIC_PACKED_ROOT_TABLE_LOAD75": "0", "KIC_STATE28": "1",
               "KIC_INDEXED_ROOT5": "1", "KIC_ROOT_BLOOM": "0",
               "KIC_SYMMETRIC_S3": "1",
               "KIC_INDEXED5_LOAD75": "1" if dense_request else "0",
               "KIC_SPARSE_CANONICAL": "0", "KIC_RANK_THREADS": "4",
               "KIC_TARGET_THREADS": "1", "KIC_DUMP_BASE": str(folder / "base_dump.jsonl")}
        command = [str(binary), f"construct:83:1:{k}", str(POINT), "29",
                   str(folder / "targets_result.jsonl")]
        arm_usage = one_run(folder, command, env, cap, group)
        print(f"{arm}: {arm_usage['status']} after "
              f"{arm_usage['wall_seconds']:.3f}s", flush=True)
    rho_folder = output / "runs" / "rho_primary"
    rho_command = [str(rho_binary), str(POINT), "2026101301", str(RHO_CONFIG),
                   str(rho_folder / "receipt.json")]
    print(f"starting same-point strong rho at {now()}", flush=True)
    rho_usage = one_run(rho_folder, rho_command,
                        {"PATH": os.environ.get("PATH", "/usr/bin:/bin")}, 120, group)
    print(f"strong rho: {rho_usage['status']} after "
          f"{rho_usage['wall_seconds']:.3f}s", flush=True)
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--rho-binary", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.preflight:
        print(json.dumps({"hashes": check_inputs(), "host": host_receipt(cgroup_root())},
                         sort_keys=True))
        return 0
    if not args.binary or not args.rho_binary or not args.output:
        parser.error("--binary, --rho-binary and --output are required")
    return run(args.binary.resolve(), args.rho_binary.resolve(), args.output.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
