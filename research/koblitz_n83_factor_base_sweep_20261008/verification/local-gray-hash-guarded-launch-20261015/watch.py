#!/usr/bin/env python3
"""Wait briefly for an isolated Mac, then launch the authorized N83 pair once."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time


HERE = Path(__file__).resolve().parent
PAIR = HERE.parent / "k600-gray-hash-paired-protocol-20261014"
sys.path.insert(0, str(PAIR))
import preflight  # noqa: E402
import supervisor  # noqa: E402


IMAGE = "sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31"
QUIET_SAMPLES = 3
INTERVAL_SECONDS = 45
MAX_WAIT_SECONDS = 900


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def call(argv: list[str], timeout: int = 30) -> dict:
    started = time.monotonic()
    try:
        done = subprocess.run(argv, capture_output=True, text=True,
                              timeout=timeout, check=False)
        return {"argv": argv, "exit_code": done.returncode,
                "stdout": done.stdout, "stderr": done.stderr,
                "wall_seconds_L0": time.monotonic() - started}
    except (OSError, subprocess.SubprocessError) as error:
        return {"argv": argv, "exit_code": None, "stdout": "",
                "stderr": f"{type(error).__name__}: {error}",
                "wall_seconds_L0": time.monotonic() - started}


def sample() -> dict:
    top = call(["top", "-l", "1", "-n", "0"], 30)
    load = call(["sysctl", "-n", "vm.loadavg"])
    cores = call(["sysctl", "-n", "hw.logicalcpu"])
    docker = call(["docker", "ps", "--format", "{{.Names}}"])
    result = {"observed_at_utc": datetime.now(timezone.utc).isoformat(),
              "top": top, "load": load, "cores": cores, "docker": docker,
              "system_free_bytes": shutil.disk_usage("/private/tmp").free}
    if any(row["exit_code"] != 0 for row in (top, load, cores, docker)):
        result["status"] = "PROBE_FAILURE"
        return result
    match = re.search(r"CPU usage: .*?([0-9]+(?:\.[0-9]+)?)% idle", top["stdout"])
    numbers = re.findall(r"[0-9]+(?:\.[0-9]+)?", load["stdout"])
    try:
        result["idle_percent"] = float(match.group(1)) if match else None
        result["load_1m"] = float(numbers[0])
        result["logical_cpus"] = int(cores["stdout"].strip())
    except (ValueError, IndexError, AttributeError):
        result["status"] = "PROBE_PARSE_FAILURE"
        return result
    result["other_containers"] = [name for name in docker["stdout"].splitlines()
                                  if name and name != "conductor-db"]
    result["status"] = "QUIET" if is_quiet(result) else "BUSY"
    return result


def is_quiet(row: dict) -> bool:
    return (row.get("idle_percent") is not None and row["idle_percent"] >= 30.0
            and isinstance(row.get("logical_cpus"), int) and row["logical_cpus"] > 0
            and row.get("load_1m", float("inf")) <= row["logical_cpus"]
            and row.get("other_containers") == []
            and row.get("system_free_bytes", 0) >= 10 * (1 << 30))


def final_admission(auth: Path, archives: Path, output: Path, commit: str,
                    sources: dict[str, str]) -> dict:
    if supervisor.PRIOR.clean_commit(preflight.ROOT) != commit:
        raise ValueError("source commit changed while waiting")
    if supervisor.source_hashes(preflight.ROOT) != sources:
        raise ValueError("frozen source bytes changed while waiting")
    if output.exists():
        raise ValueError("authorized output path already exists")
    plan = preflight.make_plan(archives)
    plan_bytes = preflight.PLAN.read_bytes()
    if plan_bytes != (json.dumps(plan, indent=2, sort_keys=True) + "\n").encode():
        raise ValueError("paired plan replay changed")
    record = json.loads(auth.read_text())
    error = supervisor.authorization_error(record, output, hashlib.sha256(plan_bytes).hexdigest())
    if error:
        raise ValueError(error)
    image = call(["docker", "image", "inspect", IMAGE,
                  "--format", "{{.Id}} {{.Architecture}} {{.Os}}"])
    if image["exit_code"] != 0 or image["stdout"].strip() != f"{IMAGE} arm64 linux":
        raise ValueError("pinned local Linux ARM image unavailable")
    cgroup = call(["docker", "run", "--rm", "--pull", "never",
                   "--platform", "linux/arm64", "--network", "none",
                   "--cpus", "1", "--memory", "4096m", "--memory-swap", "4096m",
                   IMAGE, "sh", "-c",
                   "cat /sys/fs/cgroup/memory.max /sys/fs/cgroup/memory.swap.max"])
    if cgroup["exit_code"] != 0 or cgroup["stdout"].splitlines() != ["4294967296", "0"]:
        raise ValueError("local cgroup failed 4 GiB zero-swap admission")
    return {"plan_sha256": hashlib.sha256(plan_bytes).hexdigest(),
            "authorization_sha256": sha256(auth), "image": image,
            "cgroup": cgroup,
            "archive_sha256": {p.name: sha256(p) for p in sorted(archives.glob("*.jsonl.gz"))}}


def write_json(path: Path, value: object) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def run(auth: Path, archives: Path, output: Path, watch_dir: Path,
        max_wait_seconds: int = MAX_WAIT_SECONDS) -> dict:
    if watch_dir.exists() or output.exists():
        raise ValueError("watch and authorized output paths must be fresh")
    record = json.loads(auth.read_text())
    if Path(record["output_dir"]).resolve() != output.resolve():
        raise ValueError("watch output differs from authorization")
    checkout = preflight.ROOT.resolve(strict=True)
    commit = supervisor.PRIOR.clean_commit(checkout)
    sources = supervisor.source_hashes(checkout)
    watch_dir.mkdir(parents=True, exist_ok=False)
    config = {"schema": "n83.local-paired-quiet-watch-config/v1",
              "source_commit": commit, "source_sha256": sources,
              "watcher_sha256": sha256(HERE / "watch.py"),
              "authorization_sha256": sha256(auth),
              "max_wait_seconds": max_wait_seconds,
              "quiet_samples_required": QUIET_SAMPLES,
              "sample_interval_seconds": INTERVAL_SECONDS,
              "authorized_output": str(output.resolve())}
    write_json(watch_dir / "config.json", config)
    start = time.monotonic()
    deadline = start + max_wait_seconds
    quiet = 0
    launched = False
    launch_error = None
    supervisor_exit = None
    admitted = None
    with (watch_dir / "samples.jsonl").open("x", encoding="utf-8") as samples:
        while time.monotonic() < deadline:
            row = sample()
            samples.write(json.dumps(row, sort_keys=True) + "\n")
            samples.flush()
            os.fsync(samples.fileno())
            quiet = quiet + 1 if row["status"] == "QUIET" else 0
            if quiet >= QUIET_SAMPLES:
                # Recheck the immutable inputs and cgroup at the point of launch.
                try:
                    admitted = final_admission(auth, archives, output, commit, sources)
                    last = sample()
                    samples.write(json.dumps(last, sort_keys=True) + "\n")
                    samples.flush()
                    os.fsync(samples.fileno())
                    if last["status"] != "QUIET":
                        quiet = 0
                        continue
                    write_json(watch_dir / "admission.json", admitted)
                    command = [sys.executable, str(PAIR / "supervisor.py"),
                               "--authorization", str(auth), "--output-dir", str(output),
                               "--archive-dir", str(archives)]
                    with (watch_dir / "supervisor.stdout").open("xb") as stdout, \
                         (watch_dir / "supervisor.stderr").open("xb") as stderr:
                        child = subprocess.Popen(command, stdout=stdout, stderr=stderr,
                                                 start_new_session=True)
                        launched = True
                        try:
                            supervisor_exit = child.wait(timeout=5100)
                        except subprocess.TimeoutExpired as error:
                            os.killpg(child.pid, signal.SIGTERM)
                            try:
                                supervisor_exit = child.wait(timeout=20)
                            except subprocess.TimeoutExpired:
                                os.killpg(child.pid, signal.SIGKILL)
                                supervisor_exit = child.wait()
                            raise ValueError("paired supervisor exceeded 5,100-second outer cap") from error
                    break
                except (OSError, ValueError, subprocess.SubprocessError) as error:
                    launch_error = f"{type(error).__name__}: {error}"
                    break
            remaining = deadline - time.monotonic()
            if remaining > 0:
                time.sleep(min(INTERVAL_SECONDS, remaining))
    status = ("PRODUCER_FAILURE_GUARDED_LAUNCH" if launch_error else
              "PRODUCER_FAILURE_SUPERVISOR_EXIT" if launched and supervisor_exit != 0 else
              "LAUNCHED_PAIRED_SUPERVISOR" if launched else
              "WAIT_EXPIRED_NO_SOLVER_LAUNCH")
    result = {"schema": "n83.local-paired-quiet-watch-result/v1",
              "status": status, "launched": launched, "launch_error": launch_error,
              "supervisor_exit_code": supervisor_exit,
              "admission_sha256": sha256(watch_dir / "admission.json") if admitted else None,
              "output_exists": output.exists(),
              "elapsed_watch_seconds_L0": time.monotonic() - start,
              "conservative_search_charge_seconds": 0 if not launched else None}
    write_json(watch_dir / "result.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authorization", type=Path, required=True)
    parser.add_argument("--archives", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--watch-dir", type=Path, required=True)
    parser.add_argument("--max-wait-seconds", type=int, default=MAX_WAIT_SECONDS)
    args = parser.parse_args()
    result = run(args.authorization, args.archives, args.output, args.watch_dir,
                 args.max_wait_seconds)
    print(json.dumps(result, sort_keys=True))
    return 1 if result["status"].startswith("PRODUCER_FAILURE") else 0


if __name__ == "__main__":
    raise SystemExit(main())
