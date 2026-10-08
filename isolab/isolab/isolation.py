"""Reserve hardware for one job and prove it was quiet.

A :class:`Reservation` is a context manager around a job: it takes the host
lock, builds a cpuset partition for the reserved CPUs, evicts every other
thread, moves interrupts away, pins the frequency where it can, and on exit
restores all of it. Between enter and exit the job runs; ``quiesce`` shows
the machine was quiet before it started and ``post_checks`` says what the
samplers saw while it ran. Every mechanism records whether it worked, and
:func:`max_tier` says in advance how far a host can go.
"""
from __future__ import annotations

import contextlib
import errno
import fcntl
import logging
import os
import platform
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable

from . import protocol
from .inventory import Paths, format_cpulist, kernel_version, parse_cpulist, read_int, read_text
from .samplers import Sampler

log = logging.getLogger("isolab.isolation")
#: interrupts per second per reserved CPU above which something besides the tick is landing there
IRQ_RATE_LIMIT = 1500.0
DEFAULT_LOCK = "/run/isolab.lock" if os.access("/run", os.W_OK) else "/tmp/isolab.lock"


@dataclass
class Check:
    name: str
    phase: str            # pre | run | post
    status: str           # pass | fail | unavailable | info
    value: Any = None
    threshold: Any = None
    detail: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {k: v for k, v in asdict(self).items() if v is not None or k in ("value",)}


class QuiesceError(RuntimeError):
    """The host would not go quiet; an infrastructure condition, retried elsewhere."""


class FidelityUnmet(RuntimeError):
    """A required mechanism could not be provided on this host."""


# -- capabilities ------------------------------------------------------------

def probe_capabilities(paths: Paths = Paths(), cgroup_root: str = "/sys/fs/cgroup",
                       perf: dict[str, Any] | None = None, numactl: bool = False,
                       container_backend: bool = False) -> dict[str, Any]:
    linux = platform.system() == "Linux"
    root = hasattr(os, "geteuid") and os.geteuid() == 0
    cg = Path(cgroup_root)
    controllers = (read_text(cg / "cgroup.controllers") or "").split() if linux else []
    writable = False
    if linux and root and controllers:
        probe = cg / "isolab.probe"
        try:
            probe.mkdir(exist_ok=True)
            writable = True
        except OSError:
            writable = False
        finally:
            with contextlib.suppress(OSError):
                probe.rmdir()
    kv = kernel_version()
    cpuset = writable and "cpuset" in controllers
    partition = cpuset and kv >= (5, 15)
    isolated = cpuset and kv >= (6, 1)
    irq_dir = Path(paths.proc) / "irq"
    irq_ok = linux and root and irq_dir.is_dir() and os.access(irq_dir / "default_smp_affinity", os.W_OK)
    cf = Path(paths.sys) / "devices/system/cpu/cpu0/cpufreq/scaling_governor"
    cpufreq_writable = linux and root and cf.exists() and os.access(cf, os.W_OK)
    turbo = None
    for cand in ("intel_pstate/no_turbo", "cpufreq/boost"):
        p = Path(paths.sys) / "devices/system/cpu" / cand
        if p.exists():
            turbo = str(p)
            break
    perf = perf or {"available": False, "hw_events": False}
    return {
        "linux": linux, "root": root, "kernel": platform.release(),
        "cgroup_v2": bool(controllers), "cgroup_writable": writable, "cpuset_cgroup": cpuset,
        "cpu_partition": partition, "isolated_partition": isolated,
        "memory_cgroup": writable and "memory" in controllers,
        "pids_cgroup": writable and "pids" in controllers,
        "evict": linux and root,
        "irq_affinity": irq_ok,
        "numa_bind": cpuset or numactl or container_backend,
        "numa_bind_method": "cpuset.mems" if cpuset else ("numactl" if numactl else ("container" if container_backend else None)),
        "cpufreq_writable": cpufreq_writable,
        "turbo_control": turbo,
        "turbo_writable": bool(turbo) and root and os.access(turbo, os.W_OK),
        "perf_sw": bool(perf.get("available")), "perf_hw": bool(perf.get("hw_events")),
        "psi": (Path(paths.proc) / "pressure/cpu").exists(),
        "drop_caches": linux and root,
        "affinity": hasattr(os, "sched_setaffinity"),
        "cgroup_kill": linux and kv >= (5, 14),
    }


def max_tier(caps: dict[str, Any]) -> str:
    if not caps.get("linux"):
        return "D"
    if caps["isolated_partition"] and caps["evict"] and caps["irq_affinity"] and caps["numa_bind"] and caps["perf_hw"]:
        return "A"
    if caps["cpuset_cgroup"] and caps["evict"] and caps["irq_affinity"] and caps["numa_bind"]:
        return "B"
    if caps["linux"] and caps["affinity"]:
        return "C"
    return "D"


def tier_from_mechanisms(m: dict[str, Any], caps: dict[str, Any]) -> str:
    """The tier a reservation actually achieved."""
    if m.get("partition_state") in ("isolated",) and m.get("evicted") and m.get("irqs_ok") \
            and m.get("numa_bound") and caps.get("perf_hw"):
        return "A"
    if m.get("cgroup") and m.get("evicted") and m.get("irqs_ok") and m.get("numa_bound"):
        return "B"
    if m.get("pinned"):
        return "C"
    return "D"


# -- cgroup helpers ----------------------------------------------------------

def _write(path: Path, value: str) -> str | None:
    """Write a sysfs/cgroupfs file; return an error string or None."""
    try:
        path.write_text(value)
        return None
    except OSError as err:
        return f"{path}: {err.strerror or err}"


def _enable_controllers(cg: Path, wanted: list[str]) -> list[str]:
    have = (read_text(cg / "cgroup.controllers") or "").split()
    enabled = (read_text(cg / "cgroup.subtree_control") or "").split()
    errs = []
    for c in wanted:
        if c in have and c not in enabled:
            err = _write(cg / "cgroup.subtree_control", f"+{c}")
            if err:
                errs.append(err)
    return errs


def kill_cgroup(cg: Path, grace: float = 3.0) -> None:
    if not cg.exists():
        return
    if (cg / "cgroup.kill").exists():
        _write(cg / "cgroup.kill", "1")
    else:
        import signal
        for pid in _procs(cg):
            with contextlib.suppress(ProcessLookupError, PermissionError):
                os.kill(pid, signal.SIGKILL)
    deadline = time.monotonic() + grace
    while _procs(cg) and time.monotonic() < deadline:
        time.sleep(0.05)


def _procs(cg: Path) -> set[int]:
    out: set[int] = set()
    for f in [cg / "cgroup.procs", *cg.glob("*/cgroup.procs"), *cg.glob("*/*/cgroup.procs")]:
        try:
            out |= {int(x) for x in f.read_text().split()}
        except (OSError, ValueError):
            pass
    return out


def remove_cgroup(cg: Path) -> None:
    if not cg.exists():
        return
    kill_cgroup(cg)
    for child in sorted(cg.glob("*/"), reverse=True):
        if child.is_dir():
            remove_cgroup(child)
    with contextlib.suppress(OSError):
        if (cg / "cpuset.cpus.partition").exists():
            _write(cg / "cpuset.cpus.partition", "member")
        cg.rmdir()


# -- thread eviction ---------------------------------------------------------

def evict(cpus: set[int], keep: set[int], proc: str = "/proc") -> dict[int, set[int]]:
    """Move every other thread off ``cpus``; return the affinities to restore."""
    moved: dict[int, set[int]] = {}
    for task_dir in Path(proc).glob("[0-9]*/task/[0-9]*"):
        try:
            tid = int(task_dir.name)
            pid = int(task_dir.parent.parent.name)
        except ValueError:
            continue
        if pid in keep:
            continue
        try:
            current = os.sched_getaffinity(tid)
        except OSError:
            continue
        if not current & cpus:
            continue
        target = current - cpus
        if not target:
            continue
        try:
            os.sched_setaffinity(tid, target)
            moved[tid] = current
        except OSError:
            pass
    return moved


def restore_affinity(moved: dict[int, set[int]], cpus: set[int] | None = None) -> None:
    """Give back the CPUs :func:`evict` took.

    With ``cpus`` (the reserved set) only those CPUs are added back to each
    thread's current affinity, which commutes with other slots evicting and
    restoring their own disjoint CPUs at the same time; without it the saved
    affinity is written back whole.
    """
    for tid, affinity in moved.items():
        with contextlib.suppress(OSError):
            if cpus is None:
                os.sched_setaffinity(tid, affinity)
            else:
                os.sched_setaffinity(tid, os.sched_getaffinity(tid) | (affinity & cpus))


def remaining_on(cpus: set[int], keep: set[int], proc: str = "/proc") -> dict[str, Any]:
    user, kernel = [], 0
    for task_dir in Path(proc).glob("[0-9]*/task/[0-9]*"):
        try:
            pid = int(task_dir.parent.parent.name)
            if pid in keep:
                continue
            if not os.sched_getaffinity(int(task_dir.name)) & cpus:
                continue
            if (task_dir.parent.parent / "cmdline").read_bytes():
                user.append(f"{(task_dir / 'comm').read_text().strip()}[{task_dir.name}]")
            else:
                kernel += 1
        except (OSError, ValueError):
            continue
    return {"user_threads": user[:20], "user_thread_count": len(user), "kernel_threads": kernel}


# -- interrupts --------------------------------------------------------------

def move_irqs(reserved: set[int], housekeeping: set[int], proc: str = "/proc") -> dict[str, Any]:
    """Point every IRQ that targets a reserved CPU at the housekeeping CPUs.

    Per-CPU and kernel-managed interrupts (the local timer, IPIs, managed
    device queues) refuse the write with EIO or EINVAL on every host; they
    are listed as ``unmovable`` rather than counted as a failed mechanism.
    Whether interrupts then actually land on the reserved CPUs is measured
    during the run (``irqs_on_reserved``); boot-time ``isolcpus=managed_irq``
    and ``irqaffinity`` are what remove the managed ones.
    """
    saved: dict[str, str] = {}
    added: dict[str, list[int]] = {}
    failed: list[str] = []
    unmovable: list[str] = []
    moved = 0
    irq_root = Path(proc) / "irq"
    for d in irq_root.glob("[0-9]*"):
        f = d / "smp_affinity_list"
        cur_text = read_text(f)
        if cur_text is None:
            continue
        try:
            cur = set(parse_cpulist(cur_text))
        except ValueError:
            continue
        if not cur & reserved:
            continue
        target = cur - reserved or housekeeping
        try:
            f.write_text(format_cpulist(target))
        except OSError as err:
            name = read_text(d / "actions") or d.name
            if err.errno in (errno.EIO, errno.EINVAL, errno.ENOSPC):
                unmovable.append(f"irq{d.name}({name})")
            else:
                failed.append(f"irq{d.name}({name}): {err.strerror}")
        else:
            saved[d.name] = cur_text
            added[d.name] = sorted(target - cur)
            moved += 1
    return {"moved": moved, "failed": failed, "unmovable": unmovable, "saved": saved, "added": added,
            "reserved": sorted(reserved)}


def housekeeping_mask(housekeeping: set[int]) -> str:
    mask = 0
    for c in housekeeping:
        mask |= 1 << c
    return f"{mask:x}"


def restore_irqs(state: dict[str, Any], proc: str = "/proc") -> None:
    """Undo :func:`move_irqs` for this reservation's CPUs only.

    Each IRQ gets back the reserved CPUs it had and loses the housekeeping
    CPUs this reservation had to add, starting from whatever its affinity is
    now; another slot doing the same for its own CPUs in either order ends
    at the original mask.
    """
    reserved = set(state.get("reserved") or [])
    for irq, text in state.get("saved", {}).items():
        f = Path(proc) / "irq" / irq / "smp_affinity_list"
        try:
            saved = set(parse_cpulist(text))
            cur = set(parse_cpulist(read_text(f) or ""))
        except ValueError:
            _write(f, text)
            continue
        want = (cur - set(state.get("added", {}).get(irq, []))) | (saved & reserved)
        _write(f, format_cpulist(want or saved))


class HostSettings:
    """Host-wide knobs several concurrent reservations may want changed at once.

    The first holder saves the original value and writes the wanted one; the
    last holder to leave writes the original back. State lives in a JSON
    file guarded by its own flock, so slots in one worker and workers in
    separate processes agree. Holders whose process is gone (a crashed
    worker) are dropped, and the saved original survives them.
    """

    def __init__(self, path: str):
        self.path = Path(path)

    @contextlib.contextmanager
    def _locked(self):
        with open(str(self.path) + ".lock", "a+") as lk:
            fcntl.flock(lk, fcntl.LOCK_EX)
            try:
                import json
                try:
                    state = json.loads(self.path.read_text())
                except (OSError, ValueError):
                    state = {}
                yield state
                tmp = self.path.with_suffix(".tmp")
                tmp.write_text(json.dumps(state, sort_keys=True))
                tmp.replace(self.path)
            finally:
                fcntl.flock(lk, fcntl.LOCK_UN)

    @staticmethod
    def _alive(holder: str) -> bool:
        try:
            os.kill(int(holder.split(":", 1)[0]), 0)
            return True
        except (ValueError, ProcessLookupError):
            return False
        except PermissionError:
            return True

    def acquire(self, key: str, file: str, value: str, holder: str) -> str | None:
        """Hold ``file`` at ``value``; return an error string or None."""
        with self._locked() as state:
            ent = state.get(key)
            if ent is not None:
                ent["holders"] = [h for h in ent["holders"] if self._alive(h)]
            if ent is None or (not ent["holders"] and ent.get("file") != file):
                cur = read_text(Path(file))
                if cur is None:
                    return f"{file}: unreadable"
                ent = {"file": file, "saved": cur, "holders": []}
            err = _write(Path(file), value)
            if err:
                if not ent["holders"]:
                    state.pop(key, None)
                return err
            ent["value"] = value
            ent["holders"].append(holder)
            state[key] = ent
            return None

    def release(self, key: str, holder: str) -> None:
        with self._locked() as state:
            ent = state.get(key)
            if ent is None:
                return
            ent["holders"] = [h for h in ent["holders"] if h != holder and self._alive(h)]
            if not ent["holders"]:
                _write(Path(ent["file"]), ent["saved"])
                state.pop(key, None)


# -- locks -------------------------------------------------------------------

def _flock(path: str, mode: int, deadline: float, what: str):
    handle = open(path, "a+")
    while True:
        try:
            fcntl.flock(handle, mode | fcntl.LOCK_NB)
            return handle
        except BlockingIOError:
            if time.monotonic() >= deadline:
                handle.close()
                raise QuiesceError(f"{what} {path}")
            time.sleep(0.5)


def _drain_pending(drain: Path, max_age_s: float = 6 * 3600.0) -> bool:
    """A job that needs the whole host is waiting; a marker whose waiter died, or a very old one, is ignored."""
    try:
        holder, ts = drain.read_text().split()[:2]
    except (OSError, ValueError):
        return False
    return time.time() - float(ts) < max_age_s and HostSettings._alive(holder)


# -- the reservation ---------------------------------------------------------

class Reservation:
    def __init__(self, job_id: str, placement, fidelity: dict[str, Any], caps: dict[str, Any],
                 housekeeping: list[int], memory_mb: int | None, pids: int,
                 paths: Paths = Paths(), cgroup_root: str = "/sys/fs/cgroup",
                 lab_cgroup: str = "isolab.lab", lock_path: str = DEFAULT_LOCK,
                 lock_wait_s: float = 600.0, keep_pids: set[int] | None = None,
                 slot_lock_path: str | None = None, exclusive: bool = True,
                 co_tenants: Callable[[], dict[str, set[int]]] | None = None):
        self.job_id, self.placement, self.fidelity, self.caps = job_id, placement, fidelity, caps
        self.housekeeping = set(housekeeping)
        self.memory_mb, self.pids = memory_mb, pids
        self.paths, self.cgroup_root = paths, Path(cgroup_root)
        self.lab = self.cgroup_root / lab_cgroup
        self.job_cgroup: Path | None = None
        self.lock_path, self.lock_wait_s = lock_path, lock_wait_s
        # slots: the host lock is shared between slots and exclusive for a job that
        # must have the host to itself; each slot also holds its own lock
        self.slot_lock_path = slot_lock_path
        self.exclusive = exclusive or slot_lock_path is None
        self.co_tenants = co_tenants or (lambda: {})
        self.settings = HostSettings(lock_path + ".state")
        self._holder = f"{os.getpid()}:{job_id}"
        self._held: list[str] = []
        self._slot_lock = None
        self.keep = set(keep_pids or set()) | {os.getpid()}
        self.mechanisms: dict[str, Any] = {}
        self._lock = None
        self._moved: dict[int, set[int]] = {}
        self._irqs: dict[str, Any] = {}
        self._gov_saved: dict[int, str] = {}
        self.checks: list[Check] = []

    # ---- enter / exit
    def __enter__(self) -> "Reservation":
        self._acquire_lock()
        m = self.mechanisms
        reserved = set(self.placement.reserved)
        m["reserved_cpus"] = sorted(reserved)
        m["pinned"] = self.caps["affinity"] or self.caps["cpuset_cgroup"]
        try:
            self._cgroups()
            if self.caps["evict"]:
                self._moved = evict(reserved, self.keep, self.paths.proc)
                left = remaining_on(reserved, self.keep, self.paths.proc)
                m["evicted"] = left["user_thread_count"] == 0
                m["threads_moved"] = len(self._moved)
                m["left_on_reserved"] = left
            else:
                m["evicted"] = False
                m["threads_moved"] = 0
            if self.caps["irq_affinity"] and self.housekeeping:
                self._irqs = move_irqs(reserved, self.housekeeping, self.paths.proc)
                dflt = str(Path(self.paths.proc) / "irq/default_smp_affinity")
                if Path(dflt).exists() and self.settings.acquire(
                        "default_smp_affinity", dflt, housekeeping_mask(self.housekeeping), self._holder) is None:
                    self._held.append("default_smp_affinity")
                m["irqs_moved"] = self._irqs["moved"]
                m["irqs_failed"] = self._irqs["failed"]
                m["irqs_unmovable"] = self._irqs["unmovable"]
                m["irqs_ok"] = not self._irqs["failed"]
            else:
                m["irqs_ok"] = False
                m["irqs_moved"] = 0
            self._cpufreq()
        except Exception:
            self.__exit__(None, None, None)
            raise
        m["tier"] = tier_from_mechanisms(m, self.caps)
        return self

    def __exit__(self, *exc) -> None:
        with contextlib.suppress(Exception):
            self._restore_cpufreq()
        with contextlib.suppress(Exception):
            if self._irqs:
                restore_irqs(self._irqs, self.paths.proc)
        for key in self._held:
            with contextlib.suppress(Exception):
                self.settings.release(key, self._holder)
        self._held = []
        with contextlib.suppress(Exception):
            restore_affinity(self._moved, set(self.placement.reserved))
        with contextlib.suppress(Exception):
            if self.job_cgroup is not None:
                remove_cgroup(self.job_cgroup)
            if self.mechanisms.get("cgroup"):
                remove_cgroup(self.lab)
        for attr in ("_slot_lock", "_lock"):
            h = getattr(self, attr)
            if h is not None:
                with contextlib.suppress(OSError):
                    fcntl.flock(h, fcntl.LOCK_UN)
                    h.close()
                setattr(self, attr, None)

    # ---- pieces
    def _acquire_lock(self) -> None:
        deadline = time.monotonic() + self.lock_wait_s
        drain = Path(self.lock_path + ".drain")
        if self.exclusive:
            if self.slot_lock_path is not None:
                # tell slots not to start new jobs, so shared holders drain instead of starving us
                with contextlib.suppress(OSError):
                    drain.write_text(f"{self._holder} {time.time()}")
            try:
                self._lock = _flock(self.lock_path, fcntl.LOCK_EX, deadline, "another job holds")
            finally:
                if self.slot_lock_path is not None:
                    with contextlib.suppress(OSError):
                        if drain.read_text().split()[0] == self._holder:
                            drain.unlink()
        else:
            while _drain_pending(drain):
                if time.monotonic() >= deadline:
                    raise QuiesceError(f"a job needing the whole host is waiting for {self.lock_path}")
                time.sleep(0.5)
            self._lock = _flock(self.lock_path, fcntl.LOCK_SH, deadline, "a job with the whole host holds")
        if self.slot_lock_path is not None:
            try:
                self._slot_lock = _flock(self.slot_lock_path, fcntl.LOCK_EX, deadline, "another job holds")
            except QuiesceError:
                fcntl.flock(self._lock, fcntl.LOCK_UN)
                self._lock.close()
                self._lock = None
                raise
        self.mechanisms["lock"] = self.lock_path
        self.mechanisms["host_exclusive"] = self.exclusive
        if self.slot_lock_path is not None:
            self.mechanisms["slot_lock"] = self.slot_lock_path

    def _cgroups(self) -> None:
        m = self.mechanisms
        m["cgroup"] = False
        m["partition_state"] = None
        m["numa_bound"] = False
        if not self.caps["cpuset_cgroup"]:
            return
        errs: list[str] = []
        errs += _enable_controllers(self.cgroup_root, ["cpuset", "memory", "pids", "cpu"])
        remove_cgroup(self.lab)  # a crashed worker may have left one
        try:
            self.lab.mkdir(exist_ok=True)
        except OSError as err:
            m["cgroup_error"] = f"mkdir {self.lab}: {err}"
            return
        reserved = format_cpulist(self.placement.reserved)
        mems = format_cpulist(self.placement.mems)
        for f, v in (("cpuset.cpus", reserved), ("cpuset.mems", mems)):
            err = _write(self.lab / f, v)
            if err:
                errs.append(err)
        state = "member"
        if self.caps["cpu_partition"] and not errs:
            for want in (("isolated", "root") if self.caps["isolated_partition"] else ("root",)):
                err = _write(self.lab / "cpuset.cpus.partition", want)
                got = read_text(self.lab / "cpuset.cpus.partition") or ""
                if not err and got == want:
                    state = got
                    break
                m["partition_error"] = err or got
        m["partition_state"] = state
        errs += _enable_controllers(self.lab, ["cpuset", "memory", "pids", "cpu"])
        job = self.lab / f"job-{self.job_id}"
        try:
            job.mkdir(exist_ok=True)
        except OSError as err:
            errs.append(f"mkdir {job}: {err}")
            m["cgroup_errors"] = errs
            return
        self.job_cgroup = job
        for f, v in (("cpuset.cpus", format_cpulist(self.placement.cpus)), ("cpuset.mems", mems)):
            err = _write(job / f, v)
            if err:
                errs.append(err)
        m["numa_bound"] = (read_text(job / "cpuset.mems.effective") or "") == mems or \
                          (read_text(job / "cpuset.mems") or "") == mems
        if self.caps.get("memory_cgroup"):
            _write(job / "memory.max", f"{self.memory_mb * 1024 * 1024}" if self.memory_mb else "max")
            _write(job / "memory.swap.max", "0")
            _write(job / "memory.oom.group", "1")
        if self.caps.get("pids_cgroup"):
            _write(job / "pids.max", str(self.pids))
        m["cgroup"] = True
        m["job_cgroup"] = str(job)
        m["cgroup_parent"] = "/" + str(job.relative_to(self.cgroup_root))
        if errs:
            m["cgroup_errors"] = errs

    def _cpufreq(self) -> None:
        m = self.mechanisms
        base = Path(self.paths.sys) / "devices/system/cpu"
        m["governor_set"] = False
        m["turbo_set"] = False
        if self.fidelity.get("governor") == "performance" and self.caps["cpufreq_writable"]:
            for c in self.placement.reserved:
                f = base / f"cpu{c}/cpufreq/scaling_governor"
                cur = read_text(f)
                if cur and cur != "performance":
                    if _write(f, "performance") is None:
                        self._gov_saved[c] = cur
            m["governor_set"] = True
        if self.fidelity.get("turbo") == "off" and self.caps.get("turbo_writable"):
            p = Path(self.caps["turbo_control"])
            want = "1" if p.name == "no_turbo" else "0"
            if self.settings.acquire("turbo", str(p), want, self._holder) is None:
                self._held.append("turbo")
            m["turbo_set"] = read_text(p) == want

    def _restore_cpufreq(self) -> None:
        base = Path(self.paths.sys) / "devices/system/cpu"
        for c, gov in self._gov_saved.items():
            _write(base / f"cpu{c}/cpufreq/scaling_governor", gov)

    # ---- checks
    def pre_checks(self, host: dict[str, Any]) -> list[Check]:
        m, f, caps = self.mechanisms, self.fidelity, self.caps
        req = set(f.get("require") or [])
        out: list[Check] = []

        def need(name: str) -> bool:
            return name in req

        state = m.get("partition_state")
        if caps["cpuset_cgroup"]:
            ok = state in ("isolated", "root")
            out.append(Check("cpu_partition", "pre", "pass" if ok else ("fail" if need("cpu_partition") else "info"),
                             state, detail=m.get("partition_error")))
        else:
            out.append(Check("cpu_partition", "pre", "fail" if need("cpu_partition") else "unavailable",
                             None, detail="no writable cgroup v2 cpuset controller"))
        if caps["evict"]:
            left = m.get("left_on_reserved", {})
            ok = m.get("evicted", False)
            out.append(Check("evicted", "pre", "pass" if ok else "fail", m.get("threads_moved"),
                             detail=None if ok else f"user threads still allowed on reserved cpus: {left.get('user_threads')}"))
        else:
            out.append(Check("evicted", "pre", "fail" if need("evicted") else "unavailable", detail="not root"))
        if caps["irq_affinity"]:
            ok = m.get("irqs_ok", False)
            unmovable = m.get("irqs_unmovable") or []
            detail = None if ok else f"could not move: {m.get('irqs_failed')}"
            if unmovable:
                detail = (detail + "; " if detail else "") + f"{len(unmovable)} per-cpu/managed irqs stay (kernel refuses): {unmovable[:8]}"
            out.append(Check("irq_moved", "pre", "pass" if ok else "fail", m.get("irqs_moved"), detail=detail))
        else:
            out.append(Check("irq_moved", "pre", "fail" if need("irq_moved") else "unavailable"))
        out.append(Check("numa_bind", "pre", "pass" if m.get("numa_bound") else
                         ("fail" if need("numa_bind") else "unavailable"),
                         {"nodes": self.placement.nodes, "method": caps.get("numa_bind_method")}))
        cf = host.get("cpufreq") or {}
        base = Path(self.paths.sys) / "devices/system/cpu"
        if cf.get("available"):
            govs = {read_text(base / f"cpu{c}/cpufreq/scaling_governor") for c in self.placement.reserved}
            if f.get("governor") == "performance":
                out.append(Check("governor", "pre", "pass" if govs == {"performance"} else "fail", sorted(govs), "performance"))
            else:
                out.append(Check("governor", "pre", "info", sorted(govs)))
            t = cf.get("turbo") or {}
            if t.get("control"):
                p = Path(t["control"])
                raw = read_int(p)
                enabled = (raw == 0) if p.name == "no_turbo" else (raw == 1)
                if f.get("turbo") == "off":
                    out.append(Check("turbo", "pre", "pass" if not enabled else "fail", "on" if enabled else "off", "off"))
                else:
                    out.append(Check("turbo", "pre", "info", "on" if enabled else "off"))
            else:
                out.append(Check("turbo", "pre", "fail" if need("frequency_pinned") else "unavailable", detail="no turbo control exposed"))
        else:
            out.append(Check("governor", "pre", "fail" if need("frequency_pinned") else "unavailable", detail="no cpufreq"))
            out.append(Check("turbo", "pre", "unavailable"))
        out.append(Check("perf_counters", "pre", "pass" if caps["perf_hw"] else
                         ("fail" if need("perf_counters") else ("info" if caps["perf_sw"] else "unavailable")),
                         {"software": caps["perf_sw"], "hardware": caps["perf_hw"]}))
        virt = host.get("virtualization") or {}
        out.append(Check("bare_metal", "pre", ("fail" if virt.get("is_vm") else "pass") if need("bare_metal") else "info",
                         "vm" if virt.get("is_vm") else "bare-metal", detail=virt.get("detect") or virt.get("hypervisor")))
        knobs = host.get("knobs") or {}
        out.append(Check("thp", "pre", ("pass" if knobs.get("thp_enabled") == "never" else "fail") if need("thp_never") else "info",
                         knobs.get("thp_enabled")))
        out.append(Check("aslr", "pre", "info", knobs.get("aslr")))
        out.append(Check("numa_balancing", "pre", "info", knobs.get("numa_balancing")))
        out.append(Check("nmi_watchdog", "pre", "info", knobs.get("nmi_watchdog")))
        topo = host.get("topology") or {}
        out.append(Check("smt", "pre", ("pass" if topo.get("smt_control") in ("off", "forceoff")
                                        or topo.get("summary", {}).get("threads_per_core", 1) == 1 else "fail") if need("smt_off") else "info",
                         {"mode": self.placement.smt, "idle_siblings": self.placement.idle_siblings,
                          "control": topo.get("smt_control")}))
        cmd = host.get("cmdline") or {}
        reserved = set(self.placement.reserved)
        out.append(Check("kernel_isolation", "pre", "info",
                         {"isolcpus": reserved <= set(cmd.get("isolcpus") or []),
                          "nohz_full": reserved <= set(cmd.get("nohz_full") or []),
                          "rcu_nocbs": reserved <= set(cmd.get("rcu_nocbs") or [])}))
        if self.memory_mb:
            # read live: the inventory snapshot is from worker start and memory moves constantly
            from .inventory import cpu_topology
            live_nodes = cpu_topology(self.paths).get("nodes") or {}
            avail = sum((live_nodes.get(n) or {}).get("mem_available_kb") or 0 for n in self.placement.nodes)
            if avail:
                out.append(Check("memory_available", "pre", "pass" if avail >= self.memory_mb * 1024 else "fail",
                                 avail // 1024, self.memory_mb, detail="MiB available on the job's NUMA node(s): free + reclaimable cache"))
        out.append(Check("swap", "pre", "pass" if m.get("cgroup") and caps.get("memory_cgroup") else "info",
                         "job swap.max=0" if m.get("cgroup") else (host.get("memory") or {}).get("swap_total_kb")))
        out.append(Check("exclusive_lock", "pre", "pass" if self.exclusive else "info", self.lock_path,
                         detail=None if self.exclusive else f"host shared with other slots; this slot holds {self.slot_lock_path}"))
        out.append(Check("isolation_tier", "pre",
                         "pass" if protocol.tier_at_least(m.get("tier", "D"), f.get("min_isolation_tier")) else "fail",
                         m.get("tier"), f.get("min_isolation_tier")))
        self.checks = out
        return out

    def enforce_requirements(self, checks: list[Check]) -> None:
        """Refuse to run when a required mechanism failed or the tier is too low."""
        req = set(self.fidelity.get("require") or [])
        # checks named after their requirement, plus the ones whose name differs
        alias = {"thp": "thp_never", "smt": "smt_off", "governor": "frequency_pinned", "turbo": "frequency_pinned"}
        failed = []
        for c in checks:
            if c.status != "fail":
                continue
            if c.name in req or alias.get(c.name) in req:
                failed.append(f"{c.name}={c.value} ({c.detail})" if c.detail else f"{c.name}={c.value}")
            elif c.name == "isolation_tier":
                failed.append(f"isolation tier {c.value} is below the required {c.threshold}")
            elif c.name == "memory_available":
                failed.append(f"node has {c.value} MiB free, job wants {c.threshold} MiB")
        if failed:
            raise FidelityUnmet("; ".join(failed))

    def quiesce(self, settle_s: float, retries: int = 5, sample_period: float = 0.5) -> list[Check]:
        """Show the host is quiet before the clock starts. Raises QuiesceError, except under best_effort."""
        f = self.fidelity
        last: list[Check] = []
        if f.get("policy") == "best_effort":
            retries = 1
        for attempt in range(1, retries + 1):
            s = Sampler(self.placement.reserved, period=min(sample_period, max(settle_s, 0.2)),
                        paths=self.paths, exclude=lambda: set(), other_cpu_threshold=f["max_other_cpu"],
                        co_tenants=self.co_tenants)
            s.start()
            time.sleep(settle_s)
            s.stop()
            summ = s.summary()
            checks = self._settle_checks(summ, settle_s)
            last = checks
            if not any(c.status == "fail" for c in checks):
                return checks
            if attempt < retries:
                log.info("%s: host busy during settle (attempt %d): %s", self.job_id, attempt,
                         [(c.name, c.value) for c in checks if c.status == "fail"])
                time.sleep(min(settle_s * attempt, 30.0))
        bad = ", ".join(f"{c.name}={c.value} (limit {c.threshold})" for c in last if c.status == "fail")
        if f.get("policy") == "best_effort":
            return last  # recorded as failures, never refused
        raise QuiesceError(f"host did not go quiet in {retries} settle windows of {settle_s}s: {bad}")

    def _settle_checks(self, summ: dict[str, Any], settle_s: float) -> list[Check]:
        f = self.fidelity
        out = []
        if summ.get("n_samples", 0) == 0:
            return [Check("settle", "pre", "unavailable", detail="no samples")]
        out.append(Check("settle_other_cpu", "pre", "pass" if summ["other_cpu_ratio"] <= f["max_other_cpu"] else "fail",
                         round(summ["other_cpu_ratio"], 4), f["max_other_cpu"],
                         detail=f"top: {summ['other_top']}" if summ["other_top"] else None))
        if self.caps["psi"]:
            for kind in ("cpu", "memory", "io"):
                v = summ["psi_some_pct"].get(kind)
                if v is None:
                    out.append(Check(f"settle_psi_{kind}", "pre", "unavailable"))
                    continue
                status = "pass" if v <= f["max_psi_some_pct"] else ("fail" if kind != "io" else "info")
                detail = f"% of the settle window with {kind} pressure; avg10 was {summ['psi_some_avg10_max'].get(kind)}"
                if kind == "cpu" and not self.exclusive and status == "fail":
                    # host-wide cpu pressure on a shared host is mostly the other slots and their workers
                    # queueing on the housekeeping cpus; the reserved cpus are judged by the idle, other-cpu
                    # and irq checks here and by the job cgroup's own pressure during the run
                    status = "info"
                    detail += "; host shared with other slots, so host-wide cpu pressure is not this slot's gate"
                out.append(Check(f"settle_psi_{kind}", "pre", status, round(v, 3), f["max_psi_some_pct"], detail=detail))
        else:
            out.append(Check("settle_psi_cpu", "pre", "unavailable", detail="no PSI"))
        if summ.get("co_tenant_jobs"):
            out.append(Check("settle_co_tenants", "pre", "info", round(summ["co_tenant_cpu_s"], 4),
                             detail=f"jobs running in other slots of this host: {summ['co_tenant_jobs']}; "
                                    "their CPU is not counted as other_cpu"))
        busy = summ.get("job_cpu_busy_pct")
        if busy is None:
            out.append(Check("settle_reserved_cpu_idle", "pre", "unavailable", detail="no per-cpu accounting"))
        else:
            # only a failure when the reserved cpus were actually cleared; otherwise it is what it is
            status = "pass" if busy <= 2.0 else ("fail" if self.mechanisms.get("evicted") else "info")
            out.append(Check("settle_reserved_cpu_idle", "pre", status, round(busy, 3), 2.0,
                             detail="busy% of the reserved cpus while nothing should run there"))
        steal = summ.get("job_cpu_steal_pct")
        if steal is None:
            out.append(Check("settle_steal", "pre", "unavailable"))
        else:
            out.append(Check("settle_steal", "pre", "pass" if steal <= f["max_job_cpu_steal_pct"] else "fail",
                             round(steal, 3), f["max_job_cpu_steal_pct"]))
        irqs = summ.get("irqs_on_job_cpus")
        if irqs is not None:
            out.append(Check("settle_irqs_on_reserved", "pre", "info", round(irqs / max(settle_s, 1e-6) / max(len(self.placement.reserved), 1), 1),
                             detail="interrupts per second per reserved cpu"))
        return out


def post_checks(summ: dict[str, Any], counters: dict[str, Any] | None, cg: dict[str, Any] | None,
                launcher: dict[str, Any] | None, fidelity: dict[str, Any], mechanisms: dict[str, Any],
                wall_host_s: float | None) -> list[Check]:
    f = fidelity
    out: list[Check] = []
    if summ.get("n_samples", 0):
        out.append(Check("other_cpu", "post", "pass" if summ["other_cpu_ratio"] <= f["max_other_cpu"] else "fail",
                         round(summ["other_cpu_ratio"], 4), f["max_other_cpu"],
                         detail=f"top: {summ['other_top']}" if summ["other_top"] else None))
        cg_psi = summ.get("cg_psi_some_pct") or {}
        for kind in ("cpu", "memory"):
            host_v = summ["psi_some_pct"].get(kind)
            cg_v = cg_psi.get(kind)
            if cg_v is not None:
                # the job's own cgroup: did its tasks stall waiting for cpu or memory
                out.append(Check(f"psi_{kind}", "post", "pass" if cg_v <= f["max_psi_some_pct"] else "fail", round(cg_v, 3), f["max_psi_some_pct"],
                                 detail=f"% of the run the job's own tasks were stalled on {kind}; host-wide share {None if host_v is None else round(host_v, 3)}"))
                out.append(Check(f"host_psi_{kind}", "post", "info", None if host_v is None else round(host_v, 3)))
            elif host_v is None:
                out.append(Check(f"psi_{kind}", "post", "unavailable"))
            else:
                out.append(Check(f"psi_{kind}", "post", "pass" if host_v <= f["max_psi_some_pct"] else "fail", round(host_v, 3), f["max_psi_some_pct"],
                                 detail=f"% of the run with {kind} pressure somewhere on the host (no job cgroup to read)"))
        out.append(Check("psi_io", "post", "info", summ["psi_some_pct"].get("io")))
        if summ.get("job_cpu_steal_pct") is None:
            out.append(Check("steal", "post", "unavailable"))
            out.append(Check("reserved_cpu_busy", "post", "unavailable"))
        else:
            out.append(Check("steal", "post", "pass" if summ["job_cpu_steal_pct"] <= f["max_job_cpu_steal_pct"] else "fail",
                             round(summ["job_cpu_steal_pct"], 3), f["max_job_cpu_steal_pct"]))
            out.append(Check("reserved_cpu_busy", "post", "info", round(summ["job_cpu_busy_pct"], 2),
                             detail="busy% of the reserved cpus; low means the job did not use what it asked for"))
        if summ.get("irqs_on_job_cpus") is not None:
            ncpu = max(len(mechanisms.get("reserved_cpus") or []), 1)
            rate = summ["irqs_on_job_cpus"] / max(summ["wall_s"], 1e-6) / ncpu
            out.append(Check("irqs_on_reserved", "post", "pass" if rate <= IRQ_RATE_LIMIT else "fail", round(rate, 1), IRQ_RATE_LIMIT,
                             detail="interrupts per second per reserved cpu; the local timer tick alone is HZ (100-1000) on a busy cpu, "
                                    "near zero under nohz_full; more means device interrupts are landing on the job"))
        if summ.get("freq_cv_max") is not None:
            status = "pass" if summ["freq_cv_max"] <= 0.02 else ("fail" if f.get("governor") == "performance" else "info")
            out.append(Check("frequency_stable", "post", status, round(summ["freq_cv_max"], 4), 0.02,
                             detail="max coefficient of variation of sampled frequency over the job cpus"))
        else:
            out.append(Check("frequency_stable", "post", "unavailable"))
        if summ.get("throttle_events") is not None and summ.get("thermal_max_c") is not None:
            out.append(Check("thermal_throttle", "post", "pass" if summ["throttle_events"] == 0 else "fail",
                             summ["throttle_events"], 0, detail=f"max temperature {summ['thermal_max_c']} C"))
        out.append(Check("network_bytes", "post", "info", summ.get("net_bytes")))
        out.append(Check("contended_samples", "post", "info", f"{summ['contended_samples']}/{summ['n_samples']}"))
        if summ.get("co_tenant_jobs"):
            out.append(Check("co_tenants", "post", "info", round(summ["co_tenant_cpu_s"], 4),
                             detail=f"jobs running in other slots of this host during the run: {summ['co_tenant_jobs']}; "
                                    "they share last-level cache, memory bandwidth and the power budget, so the grade is capped at B"))
    else:
        out.append(Check("samples", "post", "unavailable", detail="sampler produced nothing"))
    if counters:
        mig = counters.get("cpu_migrations")
        if mig is None:
            out.append(Check("cpu_migrations", "post", "unavailable"))
        elif mechanisms.get("partition_state") in ("isolated", "root"):
            out.append(Check("cpu_migrations", "post", "pass" if mig == 0 else "fail", mig, 0))
        else:
            out.append(Check("cpu_migrations", "post", "info", mig))
        out.append(Check("hardware_counters", "post", "pass" if counters.get("hardware_counters") else
                         ("fail" if "perf_counters" in (f.get("require") or []) else "info"),
                         counters.get("hardware_counters")))
    else:
        out.append(Check("cpu_migrations", "post", "unavailable", detail="perf not run"))
    if cg:
        thr = cg.get("nr_throttled")
        if thr is not None:
            out.append(Check("cpu_throttled", "post", "pass" if thr == 0 else "fail", thr, 0))
        oom = cg.get("oom_kill")
        if oom is not None:
            out.append(Check("oom_kill", "post", "pass" if oom == 0 else "fail", oom, 0))
    if launcher:
        out.append(Check("involuntary_switches", "post", "info", launcher.get("involuntary_switches")))
        if wall_host_s is not None and launcher.get("wall_s") is not None:
            out.append(Check("launch_overhead_s", "post", "info", round(wall_host_s - launcher["wall_s"], 4),
                             detail="host wall minus in-container wall: exec and runtime overhead, not charged to the job"))
    return out


def grade(tier: str, any_contended_in_summary: bool, shared_host: bool = False) -> str:
    """The verdict on a run: its tier, one step lower if the summary used contended repeats.

    ``shared_host`` (another slot's job ran during a summarised repeat) caps
    the grade at B: disjoint cores still share cache, memory bandwidth and
    the package power budget, which is exactly what A promises is absent.
    """
    if tier not in protocol.TIER_ORDER or tier == "none":
        return "none"
    order = ["A", "B", "C", "D"]
    i = order.index(tier)
    if any_contended_in_summary:
        i = min(i + 1, len(order) - 1)
    if shared_host:
        i = max(i, 1)
    return order[i]
