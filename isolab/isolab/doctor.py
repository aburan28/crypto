"""`isolab doctor`: is this host ready for strict runs, and what fixes it.

Each item is a fact about the host, a verdict, and the command that changes
it. ``--apply`` runs the fixes that take effect at runtime (sysfs and sysctl
writes, irqbalance); the kernel command line is printed, never edited.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Any, Callable

from . import perfstat
from .inventory import Paths, format_cpulist, kernel_version, read_text
from .isolation import max_tier


def _item(name: str, status: str, value: Any = None, fix: str | None = None, apply: Callable[[], None] | None = None,
          reboot: bool = False, note: str | None = None) -> dict[str, Any]:
    return {"name": name, "status": status, "value": value, "fix": fix, "apply": apply, "reboot": reboot, "note": note}


def _write(path: str, value: str) -> Callable[[], None]:
    def go() -> None:
        Path(path).write_text(value)
    return go


def _sysctl(key: str, value: str) -> Callable[[], None]:
    return _write("/proc/sys/" + key.replace(".", "/"), value)


def diagnose(inv: dict[str, Any], caps: dict[str, Any], lab_cpus: list[int], paths: Paths = Paths(),
             launcher: str | None = None, which=shutil.which) -> list[dict[str, Any]]:
    items: list[dict[str, Any]] = []
    kv = kernel_version(inv.get("kernel"))
    items.append(_item("linux", "ok" if caps["linux"] else "fail", inv.get("os")))
    items.append(_item("root", "ok" if caps["root"] else "fail", caps["root"], "run the worker as root (systemd unit does)"))
    items.append(_item("kernel", "ok" if kv >= (6, 1) else ("warn" if kv >= (5, 15) else "fail"), inv.get("kernel"),
                       "6.1+ for isolated cpuset partitions", note="5.15+ gives root partitions only"))
    cg = inv.get("cgroups") or {}
    items.append(_item("cgroup_v2_cpuset", "ok" if cg.get("v2") and "cpuset" in (cg.get("controllers") or []) else "fail",
                       cg.get("controllers"), "boot with systemd.unified_cgroup_hierarchy=1 and cgroup_no_v1=all", reboot=True))
    items.append(_item("cpuset_partition", "ok" if caps["isolated_partition"] else ("warn" if caps["cpu_partition"] else "fail"),
                       "isolated" if caps["isolated_partition"] else ("root" if caps["cpu_partition"] else None)))
    tier = max_tier(caps)
    items.append(_item("max_isolation_tier", "ok" if tier == "A" else ("warn" if tier == "B" else "fail"), tier))
    cmd = inv.get("cmdline") or {}
    lab = set(lab_cpus)
    hk = sorted(set(inv["topology"]["online"]) - lab)
    covered = {k: lab <= set(cmd.get(k) or []) for k in ("isolcpus", "nohz_full", "rcu_nocbs")}
    grub = (f"isolcpus=managed_irq,domain,{format_cpulist(lab)} nohz_full={format_cpulist(lab)} "
            f"rcu_nocbs={format_cpulist(lab)} irqaffinity={format_cpulist(hk)}")
    items.append(_item("kernel_cmdline_isolation", "ok" if all(covered.values()) else "warn", covered,
                       f"add to GRUB_CMDLINE_LINUX: {grub}  (then update-grub and reboot)", reboot=True,
                       note="optional: the cpuset partition already excludes the lab cpus; this also stops the tick and RCU callbacks there"))
    if which("systemctl"):
        rc = subprocess.run(["systemctl", "is-active", "irqbalance"], capture_output=True, text=True).stdout.strip()
        items.append(_item("irqbalance", "ok" if rc != "active" else "warn", rc,
                           "systemctl disable --now irqbalance   (or IRQBALANCE_BANNED_CPULIST in /etc/default/irqbalance)",
                           apply=(lambda: subprocess.run(["systemctl", "disable", "--now", "irqbalance"], check=False)) if rc == "active" else None))
    cf = inv.get("cpufreq") or {}
    if cf.get("available"):
        govs = set((cf.get("governors") or {}).values())
        base = Path(paths.sys) / "devices/system/cpu"

        def set_perf() -> None:
            for c in lab_cpus:
                _write(str(base / f"cpu{c}/cpufreq/scaling_governor"), "performance")()
        items.append(_item("governor", "ok" if govs == {"performance"} else "warn", sorted(govs),
                           "cpupower frequency-set -g performance", apply=set_perf if govs != {"performance"} else None))
        t = cf.get("turbo") or {}
        if t.get("control"):
            want = "1" if t["control"].endswith("no_turbo") else "0"
            items.append(_item("turbo", "ok" if not t.get("enabled") else "warn", "on" if t.get("enabled") else "off",
                               f"echo {want} > {t['control']}", apply=_write(t["control"], want) if t.get("enabled") else None))
        else:
            items.append(_item("turbo", "info", "no control exposed", note="on AMD/Arm hosts use the firmware setting"))
    else:
        items.append(_item("cpufreq", "info", "not exposed", note="a VM or a fixed-frequency platform; frequency cannot be pinned from here"))
    kn = inv.get("knobs") or {}
    items.append(_item("nmi_watchdog", "ok" if kn.get("nmi_watchdog") in (0, None) else "warn", kn.get("nmi_watchdog"),
                       "sysctl -w kernel.nmi_watchdog=0", apply=_sysctl("kernel.nmi_watchdog", "0") if kn.get("nmi_watchdog") else None))
    items.append(_item("numa_balancing", "ok" if kn.get("numa_balancing") in (0, None) else "warn", kn.get("numa_balancing"),
                       "sysctl -w kernel.numa_balancing=0", apply=_sysctl("kernel.numa_balancing", "0") if kn.get("numa_balancing") else None))
    items.append(_item("ksm", "ok" if kn.get("ksm_run") in (0, None) else "warn", kn.get("ksm_run"),
                       "echo 0 > /sys/kernel/mm/ksm/run", apply=_write("/sys/kernel/mm/ksm/run", "0") if kn.get("ksm_run") else None))
    items.append(_item("thp", "info", kn.get("thp_enabled"), "echo never > /sys/kernel/mm/transparent_hugepage/enabled",
                       note="recorded per run; set to never or madvise for memory-heavy workloads with low variance"))
    items.append(_item("timer_migration", "ok" if kn.get("timer_migration") in (0, None) else "warn", kn.get("timer_migration"),
                       "sysctl -w kernel.timer_migration=0", apply=_sysctl("kernel.timer_migration", "0") if kn.get("timer_migration") else None))
    mem = inv.get("memory") or {}
    items.append(_item("swap", "ok" if not mem.get("swap_total_kb") else "warn", mem.get("swap_total_kb"), "swapoff -a",
                       apply=(lambda: subprocess.run(["swapoff", "-a"], check=False)) if mem.get("swap_total_kb") else None,
                       note="jobs run with swap.max=0 regardless"))
    items.append(_item("psi", "ok" if caps["psi"] else "fail", caps["psi"], "boot with psi=1", reboot=True))
    pe = kn.get("perf_event_paranoid")
    items.append(_item("perf", "ok" if caps["perf_hw"] else ("warn" if caps["perf_sw"] else "fail"),
                       {"software": caps["perf_sw"], "hardware": caps["perf_hw"], "paranoid": pe},
                       "apt install linux-tools-$(uname -r)  (hardware counters need bare metal or PMU passthrough)"))
    tools = inv.get("tools") or {}
    items.append(_item("podman", "ok" if tools.get("podman") else "warn", (tools.get("podman") or {}).get("version"), "apt install podman crun"))
    items.append(_item("crun", "ok" if tools.get("crun") or tools.get("runc") else "warn",
                       (tools.get("crun") or tools.get("runc") or {}).get("version"), "apt install crun"))
    items.append(_item("runsc", "ok" if tools.get("runsc") else "info", (tools.get("runsc") or {}).get("version"),
                       "deploy/install.sh --gvisor", note="only needed for runtime.oci_runtime=runsc"))
    items.append(_item("launcher", "ok" if launcher and Path(launcher).is_file() else "fail", launcher, "isolab launcher-build"))
    items.append(_item("numactl", "ok" if tools.get("numactl") else "info", bool(tools.get("numactl")), "apt install numactl",
                       note="only the direct backend without cgroups needs it"))
    virt = inv.get("virtualization") or {}
    items.append(_item("bare_metal", "ok" if not virt.get("is_vm") else "warn", virt.get("detect") or virt.get("hypervisor") or "bare-metal",
                       note="strict runs require bare metal unless the spec drops that requirement"))
    imgs = [n for img in inv.get("images") or [] for n in (img.get("names") or [])]
    items.append(_item("images", "ok" if any("isolab-base" in n for n in imgs) else "warn",
                       [n for n in imgs if "isolab" in n][:6], "isolab images build base"))
    if inv.get("gpus"):
        items.append(_item("gpus", "ok", [g["name"] for g in inv["gpus"]]))
        items.append(_item("cdi", "ok" if inv.get("cdi") else "warn", inv.get("cdi"), "nvidia-ctk cdi generate --output=/etc/cdi/nvidia.yaml"))
        pm = {g.get("persistence_mode") for g in inv["gpus"]}
        items.append(_item("gpu_persistence", "ok" if pm == {"Enabled"} else "warn", sorted(pm), "nvidia-smi -pm 1"))
    return items


def apply(items: list[dict[str, Any]]) -> list[str]:
    done = []
    for it in items:
        if it.get("apply") and it["status"] in ("warn", "fail"):
            try:
                it["apply"]()
                done.append(it["name"])
            except Exception as err:  # noqa: BLE001
                done.append(f"{it['name']} FAILED: {err}")
    return done


def render(items: list[dict[str, Any]]) -> str:
    mark = {"ok": "ok  ", "warn": "WARN", "fail": "FAIL", "info": "info"}
    lines = []
    for it in items:
        val = it["value"] if not isinstance(it["value"], (dict, list)) else str(it["value"])
        line = f"{mark[it['status']]}  {it['name']:<26} {str(val)[:60]}"
        if it["status"] in ("warn", "fail") and it.get("fix"):
            line += f"\n        fix: {it['fix']}" + ("  [reboot]" if it.get("reboot") else "")
        if it.get("note") and it["status"] != "ok":
            line += f"\n        note: {it['note']}"
        lines.append(line)
    return "\n".join(lines)


def as_json(items: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{k: v for k, v in it.items() if k != "apply"} | {"applyable": bool(it.get("apply"))} for it in items]
