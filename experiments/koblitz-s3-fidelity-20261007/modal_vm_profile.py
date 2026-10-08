"""Separate whole-process perf diagnostic in a Modal VM Sandbox.

Run: python3 experiments/koblitz-s3-fidelity-20261007/modal_vm_profile.py \
       --output NEW_PROFILE.json
The VM is terminated in a finally block. Counters cover setup and the target;
they are not target-online IPC or host-level isolation evidence.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import modal


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "modal/vm_perf_profile.json"
REMOTE_CODE = r'''
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import tempfile
import time

root = Path("/root/crypto")
base = root / "experiments/koblitz-s3-fidelity-20261007/pilot"
bins = root / "target/release/examples"

def read(path):
    try:
        return Path(path).read_text(errors="replace")[:100000]
    except OSError:
        return None

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

host = {
    "platform": platform.platform(), "uname": list(os.uname()),
    "affinity": sorted(os.sched_getaffinity(0)),
    "cpuinfo": read("/proc/cpuinfo"),
    "cgroup": read("/proc/self/cgroup"),
    "cpuset_effective": read("/sys/fs/cgroup/cpuset.cpus.effective"),
    "cpuset_partition": read("/sys/fs/cgroup/cpuset.cpus.partition"),
    "mems_effective": read("/sys/fs/cgroup/cpuset.mems.effective"),
    "cpu_max": read("/sys/fs/cgroup/cpu.max"),
    "cpu_pressure": read("/proc/pressure/cpu"),
    "memory_pressure": read("/proc/pressure/memory"),
    "perf_event_paranoid": read("/proc/sys/kernel/perf_event_paranoid"),
    "thread_siblings_cpu0": read("/sys/devices/system/cpu/cpu0/topology/thread_siblings_list"),
    "modal_image_id": os.environ.get("MODAL_IMAGE_ID"),
    "modal_cloud_provider": os.environ.get("MODAL_CLOUD_PROVIDER"),
    "modal_region": os.environ.get("MODAL_REGION"),
}
records = []
for n in (41, 53):
    fixture = json.loads((base / f"n{n}/fixtures.json").read_text())["fixtures"][0]
    target = base / f"n{n}/T001/public_target.json"
    for arm in ("baseline", "candidate"):
        binary = bins / f"s3_pair_{arm}_frozen"
        with tempfile.TemporaryDirectory(prefix="s3-vm-perf-") as tmp:
            raw = Path(tmp) / "raw.jsonl"
            argv = ["taskset", "-c", "0-14", "numactl", "--membind=0",
                    "perf", "stat", "-x", ",", "-e", "cycles,instructions,task-clock",
                    "--", str(binary), str(n), "0", "244", "20260928",
                    str(target), str(raw), "14"]
            started = time.time_ns()
            result = subprocess.run(argv, capture_output=True, text=True, timeout=180)
            finished = time.time_ns()
            report = json.loads(raw.read_text()) if raw.exists() else None
            records.append({
                "n": n, "arm": arm, "argv": argv,
                "binary_sha256": sha(binary), "target_sha256": sha(target),
                "exit_code": result.returncode, "stderr_perf_raw": result.stderr[:20000],
                "stdout_sha256": hashlib.sha256(result.stdout.encode()).hexdigest(),
                "outer_wall_ns": finished - started,
                "online_ms_inside_process": (report or {}).get("timing_ms", {}).get(
                    "target_online_after_reusable_setup"),
                "recovered_scalar": (report or {}).get("recovered_scalar"),
                "fixture_scalar": int(fixture["fixture_scalar"]),
                "verified_native_and_fixture": bool(
                    result.returncode == 0 and report and report.get("group_verified") is True
                    and report.get("target") == fixture["public_point"]
                    and report.get("recovered_scalar") == int(fixture["fixture_scalar"])),
                "raw_sha256": sha(raw) if raw.exists() else None,
            })
print(json.dumps({
    "kind": "modal_vm_s3_whole_process_perf_diagnostic_v1",
    "scope": "perf covers process setup plus one target; counters cannot be assigned to target-online interval",
    "strict_host_isolation_receipt": None, "controlled_speedup": None,
    "host": host, "records": records,
}, separators=(",", ":")), flush=True)
'''


def main(output: Path) -> None:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    spec = importlib.util.spec_from_file_location("s3_modal_app", HERE / "modal_app.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    app = modal.App("s3-pair-vm-perf-diagnostic")
    with app.run():
        sandbox = modal.Sandbox.create(
            "python3", "-c", REMOTE_CODE,
            app=app, image=module.image, cpu=16, memory=8192, timeout=240,
            experimental_options={"vm_runtime": True},
        )
        try:
            sandbox.wait()
            stdout = sandbox.stdout.read()
            stderr = sandbox.stderr.read()
            if stderr:
                raise RuntimeError(f"VM diagnostic stderr: {stderr[:2000]}")
            result = json.loads(stdout)
            if len(result["records"]) != 4:
                raise RuntimeError("VM diagnostic did not return four profile records")
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
            print(json.dumps({"output": str(output),
                              "verified": sum(r["verified_native_and_fixture"] for r in result["records"]),
                              "perf_exits": [r["exit_code"] for r in result["records"]]}))
        finally:
            sandbox.terminate()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    main(parser.parse_args().output)
