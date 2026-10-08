"""The isolation gate: which level a measured run earned.

Levels are cumulative.  A run earns the highest level whose every check, and
every lower level's checks, passed.  A check whose input the host does not
expose is ``unknown``, and ``unknown`` never passes: the level is withheld and
the reason recorded.  This is the fail-closed rule the boundary ledger
already applies to missing fields.

L0  recorded   the command ran under the session's capsule.  Operation counts
               need nothing more, because they do not depend on contention
               (AGENTS.md section 10).
L1  pinned     the child's CPU mask, read back after exec, equals the reserved
               set; the set excludes CPU 0 (housekeeping interrupts); the
               measured tree used at most the declared number of CPUs' worth
               of time; the session's reservation evicted other threads.
L2  quiet      L1, and nothing competed during the window: the child's own
               run-queue delay, hypervisor steal on the pinned CPUs, foreign
               runnable tasks on the pinned CPUs, memory stall and the
               child's involuntary context switches are all under threshold,
               and the session preflight found the machine quiet.
L3  isolated   L2, and the host was configured for measurement: the pinned
               CPUs are in isolcpus or nohz_full, the frequency is fixed
               (performance governor, no turbo/boost), SMT is off or the
               siblings are reserved, and the host is not virtualised.

A wall-clock comparison requires both sides at or above the level the spec
declares (``measurement.isolation_required``) and the same env_class_id.
"""
from __future__ import annotations

from typing import Any

from .environment import parse_cpu_list

LEVELS = ("L0", "L1", "L2", "L3")

# The isolated_bench preflight's own defaults.  A session may run its preflight
# stricter (longer settle, lower limits) and still count as quiet; one that
# loosened any of them did not show a quiet machine by this standard's test.
PREFLIGHT_LIMITS = {"settle_s": 2.0, "max_other_cpu": 0.10, "max_psi": 5.0}

DEFAULT_THRESHOLDS: dict[str, float] = {
    "max_run_delay_fraction": 0.005,       # child runnable-but-waiting / (on CPU + waiting)
    "max_steal_fraction": 0.0,             # steal jiffies on pinned CPUs / all their jiffies
    "max_foreign_sample_fraction": 0.0,    # sampler ticks that saw a foreign runnable task
    "max_memory_full_stall_us": 0.0,       # memory PSI 'full' growth over the window
    "max_involuntary_switches_per_s": 50.0,
    "parallelism_slack": 0.05,             # allowed CPU-time/wall above the declared threads
}

# Thresholds may only be tightened by a spec.  Loosening one is refused, so a
# spec cannot buy itself a level the default would not grant.
_TIGHTER_IS_LOWER = set(DEFAULT_THRESHOLDS)


def effective_thresholds(overrides: dict[str, float] | None) -> dict[str, float]:
    th = dict(DEFAULT_THRESHOLDS)
    for k, v in (overrides or {}).items():
        if k not in DEFAULT_THRESHOLDS:
            raise ValueError(f"unknown isolation threshold {k!r}")
        if k in _TIGHTER_IS_LOWER and v > DEFAULT_THRESHOLDS[k]:
            raise ValueError(f"threshold {k} may only be tightened (default {DEFAULT_THRESHOLDS[k]}, got {v})")
        th[k] = float(v)
    return th


def _c(cid: str, level: str, status: str, observed: Any, threshold: Any, note: str) -> dict[str, Any]:
    return {"id": cid, "level": level, "status": status, "observed": observed, "threshold": threshold, "note": note}


def _bool_status(ok: bool | None) -> str:
    return "unknown" if ok is None else ("pass" if ok else "fail")


def evaluate(execution: dict[str, Any], session: dict[str, Any], capsule_stable: dict[str, Any],
             declared_threads: int, overrides: dict[str, float] | None = None) -> dict[str, Any]:
    th = effective_thresholds(overrides)
    checks: list[dict[str, Any]] = []
    pinned = set(execution.get("pinned_cpus") or [])

    # ---- L1 ------------------------------------------------------------------
    aff = execution.get("child_affinity_observed")
    tree = execution.get("tree_affinity") or {}
    outside = tree.get("tasks_outside")
    aff_ok = None if aff is None else (bool(pinned) and set(aff) == pinned and not outside)
    checks.append(_c("affinity_observed", "L1", _bool_status(aff_ok),
                     {"root_after_exec": aff, "tree_tasks_outside": outside, "tree_cpus_seen": tree.get("cpus_seen")},
                     sorted(pinned),
                     "child CPU mask read back after exec, and every task of the tree on every sample: mask within "
                     "the pinned CPUs and last ran on one of them"))
    checks.append(_c("excludes_cpu0", "L1", _bool_status(bool(pinned) and 0 not in pinned), sorted(pinned), "no CPU 0",
                     "CPU 0 takes most housekeeping interrupts"))
    reservation = session.get("reservation") or {}
    left_user = ((reservation.get("left_on_reserved") or {}).get("user_threads") if reservation else None) or []
    checks.append(_c("reservation", "L1",
                     _bool_status(bool(reservation.get("evicted")) and not left_user if reservation else None),
                     {k: reservation.get(k) for k in ("evicted", "threads_moved", "left_on_reserved")} if reservation else None,
                     "other threads evicted, and no user thread left on the reserved CPUs that the session could not move",
                     "isolated_bench.evict over the whole session; a thread confined to the reservation, or one the "
                     "session lacks permission to move (another user's, when not run as root), fails this check, so "
                     "a session not run as root earns at most L0"))
    sched = execution.get("schedstat") or {}
    wall = execution.get("wall_ns") or 0
    ru = execution.get("rusage") or {}
    ru_cpu_ns = None if ru.get("user_s") is None else (ru["user_s"] + (ru.get("sys_s") or 0)) * 1e9
    # schedstat is sampled, so a process that lived between two samples is
    # missed; rusage of the reaped root includes every waited-for descendant.
    # Parallelism takes the larger, and run delay is unknown when the samples
    # saw too little of the tree's CPU time to speak for it.
    covered = sched.get("run_ns") or 0
    coverage_ok = ru_cpu_ns is None or ru_cpu_ns < 20e6 or covered >= 0.9 * ru_cpu_ns
    if (not sched.get("threads_seen") and ru_cpu_ns is None) or wall <= 0:
        checks.append(_c("parallelism", "L1", "unknown", None, declared_threads, "no schedstat or rusage for the measured tree"))
    else:
        par = max(covered, ru_cpu_ns or 0) / wall
        ok = par <= declared_threads + th["parallelism_slack"] and declared_threads <= max(len(pinned), 1)
        checks.append(_c("parallelism", "L1", _bool_status(ok), round(par, 4),
                         {"declared_threads": declared_threads, "pinned_cpus": len(pinned), "slack": th["parallelism_slack"]},
                         f"CPU time of the measured tree / wall; peak thread count {execution.get('max_threads_observed')} "
                         "is recorded but not gated, because an idle pool thread costs nothing"))

    # ---- L2 ------------------------------------------------------------------
    run_ns, wait_ns = sched.get("run_ns"), sched.get("wait_ns")
    if not sched.get("threads_seen") or (run_ns or 0) + (wait_ns or 0) == 0:
        checks.append(_c("run_delay", "L2", "unknown", None, th["max_run_delay_fraction"], "no schedstat"))
    elif not coverage_ok:
        checks.append(_c("run_delay", "L2", "unknown", {"schedstat_run_ns": run_ns, "rusage_cpu_ns": round(ru_cpu_ns)},
                         th["max_run_delay_fraction"],
                         "the 50 ms samples saw under 90 % of the tree's CPU time (short-lived processes), so its run "
                         "delay is not known"))
    else:
        frac = wait_ns / (run_ns + wait_ns)
        checks.append(_c("run_delay", "L2", _bool_status(frac <= th["max_run_delay_fraction"]), round(frac, 6),
                         th["max_run_delay_fraction"], "share of the child's runnable time spent waiting for a CPU"))
    cont = execution.get("contention") or {}
    jif = cont.get("pinned_cpu_jiffies_delta") or {}
    total = sum(jif.values()) if jif else 0
    if total == 0:
        checks.append(_c("steal", "L2", "unknown", None, th["max_steal_fraction"], "no /proc/stat delta on the pinned CPUs"))
    else:
        frac = jif.get("steal", 0) / total
        checks.append(_c("steal", "L2", _bool_status(frac <= th["max_steal_fraction"]), round(frac, 6),
                         th["max_steal_fraction"], "hypervisor steal on the pinned CPUs; a hypervisor that hides steal reads 0"))
    n = cont.get("samples") or 0
    if n == 0:
        checks.append(_c("foreign_runnable", "L2", "unknown", None, th["max_foreign_sample_fraction"],
                         "run too short for one sample"))
    else:
        frac = (cont.get("samples_with_foreign_runnable") or 0) / n
        checks.append(_c("foreign_runnable", "L2", _bool_status(frac <= th["max_foreign_sample_fraction"]), round(frac, 4),
                         th["max_foreign_sample_fraction"], f"sampler ticks with another runnable task on a pinned CPU ({n} ticks)"))
    mem = ((cont.get("psi_total_delta_us") or {}).get("memory") or {}).get("full")
    checks.append(_c("memory_stall", "L2", _bool_status(None if mem is None else mem <= th["max_memory_full_stall_us"]),
                     mem, th["max_memory_full_stall_us"], "memory PSI 'full' stall during the window"))
    ru = execution.get("rusage") or {}
    if wall > 0 and "nivcsw" in ru:
        rate = ru["nivcsw"] / (wall / 1e9)
        checks.append(_c("involuntary_switches", "L2", _bool_status(rate <= th["max_involuntary_switches_per_s"]),
                         round(rate, 2), th["max_involuntary_switches_per_s"], "preemptions of the child per second"))
    else:
        checks.append(_c("involuntary_switches", "L2", "unknown", None, th["max_involuntary_switches_per_s"], "no rusage"))
    pre = session.get("preflight") or {}
    lim = pre.get("limits") or {}
    loosened = sorted(k for k, v in PREFLIGHT_LIMITS.items()
                      if lim.get(k) is None or (lim[k] < v if k == "settle_s" else lim[k] > v))
    quiet_ok = None if not pre else (bool(pre.get("quiet")) and not loosened)
    checks.append(_c("preflight_quiet", "L2", _bool_status(quiet_ok),
                     {**({k: pre.get(k) for k in ("other_cpu_seconds", "psi_some_avg10_max")} if pre else {}),
                      "limits": lim or None, "loosened": loosened},
                     {"isolated_bench preflight at least as strict as": PREFLIGHT_LIMITS},
                     "the machine was quiet before the session started, by a preflight no looser than the default"))

    # ---- L3 ------------------------------------------------------------------
    topo = capsule_stable.get("topology", {})
    iso = set(parse_cpu_list(topo.get("isolated"))) | set(parse_cpu_list(topo.get("nohz_full")))
    checks.append(_c("kernel_isolation", "L3", _bool_status(bool(pinned) and pinned <= iso),
                     {"isolated": topo.get("isolated"), "nohz_full": topo.get("nohz_full")}, sorted(pinned),
                     "pinned CPUs must be in isolcpus or nohz_full"))
    cpus_info = topo.get("cpus") or {}
    govs = {cpus_info.get(str(c), {}).get("governor") for c in pinned}
    gov_ok = None if govs in (set(), {None}) else govs == {"performance"}
    no_turbo, boost = topo.get("intel_pstate_no_turbo"), topo.get("cpufreq_boost")
    turbo_ok = None if (no_turbo, boost) == (None, None) else (no_turbo == "1" or boost == "0")
    freq_ok = None if gov_ok is None or turbo_ok is None else (gov_ok and turbo_ok)
    checks.append(_c("fixed_frequency", "L3", _bool_status(freq_ok),
                     {"governors": sorted(map(str, govs)), "no_turbo": no_turbo, "boost": boost},
                     "performance governor and turbo/boost off", "frequency scaling moves wall time between runs"))
    smt_active = topo.get("smt_active")
    siblings_ok = None
    if smt_active == "1":
        siblings = set()
        for c in pinned:
            siblings |= set(parse_cpu_list(cpus_info.get(str(c), {}).get("thread_siblings")))
        # A session runs one execution at a time, so a reserved CPU other than
        # the pinned ones is idle: a sibling inside the reservation is quiet.
        siblings_ok = siblings <= (pinned | set(reservation.get("cpus") or []))
    elif smt_active == "0":
        siblings_ok = True
    checks.append(_c("smt", "L3", _bool_status(siblings_ok), {"smt_active": smt_active, "control": topo.get("smt_control")},
                     "SMT off, or every sibling of a pinned CPU reserved", "an SMT sibling shares the core's execution units"))
    # A container reports itself before the VM it runs in ("docker" inside a
    # KVM guest), so the hypervisor is read with --vm and from the CPU flag.
    vz = capsule_stable.get("virtualization") or {}
    vm = vz.get("detect_vm", vz.get("detect_virt"))
    cpu = capsule_stable.get("cpu") or {}
    # Only x86 /proc/cpuinfo has a hypervisor flag; elsewhere its absence says nothing.
    hyper = (cpu.get("features") or {}).get("hypervisor") if cpu.get("vendor_id") else None
    if vm is None and hyper is None:
        metal = None
    else:
        metal = vm in (None, "none") and not hyper
    checks.append(_c("bare_metal", "L3", _bool_status(metal), {"detect_vm": vm, "hypervisor_flag": hyper},
                     {"detect_vm": "none", "hypervisor_flag": False},
                     "a virtual CPU can be descheduled by its host without a trace inside the guest"))

    earned = "L0"
    for lvl in LEVELS[1:]:
        upto = LEVELS.index(lvl)
        if all(c["status"] == "pass" for c in checks if LEVELS.index(c["level"]) <= upto):
            earned = lvl
        else:
            break
    nxt = LEVELS.index(earned) + 1
    blocking = [c["id"] for c in checks if nxt < len(LEVELS) and LEVELS.index(c["level"]) <= nxt and c["status"] != "pass"]
    return {"earned_level": earned, "blocking_next_level": blocking, "checks": checks, "thresholds": th}


def level_at_least(level: str, required: str) -> bool:
    return LEVELS.index(level) >= LEVELS.index(required)
