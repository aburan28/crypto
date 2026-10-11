"""Read-only Linux preflight for the paired-S3 strict CPU measurement gate.

Run on a proposed worker before uploading binaries or starting timed work.
This is a rejection probe, not an isolab host receipt or a grant of isolation.
It makes no host changes and uses only the Python standard library.
"""

from __future__ import annotations

import datetime
import hashlib
import json
import os
import platform
from pathlib import Path
import shutil


def read(path: str | Path) -> str | None:
    try:
        return Path(path).read_text(errors="replace").strip()
    except OSError:
        return None


def main() -> None:
    allowed = sorted(os.sched_getaffinity(0))
    topology = []
    for cpu in allowed:
        root = Path(f"/sys/devices/system/cpu/cpu{cpu}")
        topology.append({
            "cpu": cpu,
            "package_id": read(root / "topology/physical_package_id"),
            "core_id": read(root / "topology/core_id"),
            "thread_siblings": read(root / "topology/thread_siblings_list"),
            "numa_node": next((node.name for node in root.glob("node*")
                               if node.name[4:].isdigit()), None),
        })
    visible_cores = {(row["package_id"], row["core_id"]) for row in topology}
    cpuinfo = read("/proc/cpuinfo") or ""
    model = next((line.partition(":")[2].strip() for line in cpuinfo.splitlines()
                  if line.startswith("model name")), None)
    cgroup = read("/proc/self/cgroup") or ""
    v2_controllers = read("/sys/fs/cgroup/cgroup.controllers")
    partition = read("/sys/fs/cgroup/cpuset.cpus.partition")
    exclusive = read("/sys/fs/cgroup/cpuset.cpus.exclusive.effective")
    perf = bool(shutil.which("perf"))
    reasons = []
    if v2_controllers is None:
        reasons.append("cgroup v2 controllers not exposed")
    if partition != "isolated" or not exclusive:
        reasons.append("no visible isolated and exclusive cpuset partition")
    if len(visible_cores) < 15:
        reasons.append("fewer than 15 visible physical cores for SMT-isolated execution")
    if not perf:
        reasons.append("perf executable unavailable for hardware-counter preflight")
    unverified = [
        "physical-host-wide exclusivity and IRQ routing require an isolab host receipt",
        "hardware PMU access requires an actual perf counter test",
        "quiet-settle and A/A noise checks require timed intervals",
    ]
    result = {
        "kind": "s3_read_only_host_isolation_preflight_v1",
        "captured_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "scope": "visible process namespace only; no timed workload executed",
        "os": platform.system(), "kernel": platform.release(),
        "architecture": platform.machine(), "effective_uid": os.geteuid(),
        "cpu_model_visible": model,
        "dmi_product_name": read("/sys/class/dmi/id/product_name"),
        "affinity_cpus": allowed, "visible_physical_core_count": len(visible_cores),
        "visible_topology": topology,
        "self_cgroup_sha256": hashlib.sha256(cgroup.encode()).hexdigest(),
        "docker_cgroup_visible": "/docker/" in cgroup,
        "cgroup_v2_controllers": v2_controllers,
        "cpuset_partition": partition,
        "cpuset_exclusive_effective": exclusive,
        "cpuset_v1": read("/sys/fs/cgroup/cpuset/cpuset.cpus"),
        "cpu_quota_v1": read("/sys/fs/cgroup/cpu/cpu.cfs_quota_us"),
        "cpu_period_v1": read("/sys/fs/cgroup/cpu/cpu.cfs_period_us"),
        "cpu_max_v2": read("/sys/fs/cgroup/cpu.max"),
        "cpu_pressure": read("/proc/pressure/cpu"),
        "memory_pressure": read("/proc/pressure/memory"),
        "perf_event_paranoid": read("/proc/sys/kernel/perf_event_paranoid"),
        "perf_on_path": perf, "numactl_on_path": bool(shutil.which("numactl")),
        "isolab_on_path": bool(shutil.which("isolab")),
        "strict_preflight_pass": False if reasons else None,
        "strict_rejection_reasons": reasons,
        "unverified_requirements": unverified,
        "strict_host_receipt": None, "controlled_speedup": None,
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
