"""One pinned, sampled execution of a measured command.

The child is pinned with ``sched_setaffinity`` between fork and exec, so the
pin covers its whole lifetime and every thread it creates.  The fork happens
before the sampler thread starts: ``preexec_fn`` is only safe in a
single-threaded parent.

While the child runs, a sampler records every ``interval`` seconds:

* the peak number of threads in the measured process tree;
* each thread's ``/proc/<pid>/task/<tid>/schedstat`` (time on CPU, time
  runnable but waiting for a CPU, timeslices);
* the tasks that are runnable on the pinned CPUs but are not ours.

When the child exits it is left unreaped (``waitid(..., WNOWAIT)``) long enough
to read its own schedstat, then reaped with ``wait4`` for its rusage.  Around
the run, ``/proc/stat`` gives the pinned CPUs' jiffies including hypervisor
steal, and ``/proc/pressure`` gives the PSI totals.

The environment the child sees is built, not inherited: a fixed allowlist of
variables needed to find programs, plus exactly what the spec declares.  An
engine knob exported in the operator's shell (KIC_*, IC_*, GAUDRY_*, ...)
therefore cannot change a measured run without appearing in its spec.
"""
from __future__ import annotations

import os
import resource
import signal
import subprocess
import threading
import time
from dataclasses import dataclass, field
from typing import Any

from .canonical import sha256_file
from .environment import proc_stat, psi

# Variables passed through from the parent.  Everything else must be declared.
PASSTHROUGH = ("PATH", "HOME", "LANG", "LC_ALL", "TZ", "TMPDIR", "USER",
               "CARGO_HOME", "RUSTUP_HOME", "SAGE_ROOT", "LD_LIBRARY_PATH")

# Pinned defaults for every measured child.  A spec may override them only by
# declaring the same key, which then appears in its identity.
DEFAULT_ENV = {
    "RAYON_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "PYTHONHASHSEED": "0",
    "LC_ALL": "C",
}

# Prefixes of variables known to change engine behaviour in these repositories.
# They are reported if present in the parent so an operator can see what was
# dropped.
KNOB_PREFIXES = ("KIC_", "IC_", "GAUDRY_", "SOLVER_", "RAYON_", "OMP_", "WDSAT_", "ICMS_")


def build_env(declared: dict[str, str] | None, extra: dict[str, str] | None = None) -> tuple[dict[str, str], dict[str, Any]]:
    env = {k: os.environ[k] for k in PASSTHROUGH if k in os.environ}
    env.update(DEFAULT_ENV)
    env.update(declared or {})
    env.update(extra or {})
    dropped = sorted(k for k in os.environ if k.startswith(KNOB_PREFIXES) and k not in env)
    return env, {"declared": dict(sorted((declared or {}).items())), "defaults": DEFAULT_ENV,
                 "passthrough": [k for k in PASSTHROUGH if k in os.environ],
                 "dropped_engine_knobs": dropped}


def _tids(pid: int) -> list[int]:
    try:
        return [int(t) for t in os.listdir(f"/proc/{pid}/task")]
    except OSError:
        return []


def _schedstat(pid: int, tid: int) -> tuple[int, int, int] | None:
    try:
        with open(f"/proc/{pid}/task/{tid}/schedstat") as fh:
            a, b, c = fh.read().split()[:3]
        return int(a), int(b), int(c)
    except (OSError, ValueError):
        return None


def _proc_table() -> dict[int, tuple[int, str, int, str]]:
    """pid -> (ppid, state, last_cpu, comm) for every process."""
    out = {}
    for d in os.listdir("/proc"):
        if not d.isdigit():
            continue
        try:
            with open(f"/proc/{d}/stat") as fh:
                text = fh.read()
        except OSError:
            continue
        r = text.rfind(")")
        fields = text[r + 2:].split()
        out[int(d)] = (int(fields[1]), fields[0], int(fields[36]), text[text.find("(") + 1:r])
    return out


def _tree(root: int, table: dict[int, tuple[int, str, int, str]]) -> set[int]:
    kids: dict[int, list[int]] = {}
    for pid, (ppid, *_rest) in table.items():
        kids.setdefault(ppid, []).append(pid)
    out, stack = {root}, [root]
    while stack:
        for c in kids.get(stack.pop(), []):
            if c not in out:
                out.add(c)
                stack.append(c)
    return out


def _foreign_runnable(cpus: set[int], exclude: set[int]) -> list[dict[str, Any]]:
    found = []
    for task in _iter_tasks():
        pid, tid, state, cpu, comm = task
        if state == "R" and cpu in cpus and pid not in exclude:
            found.append({"pid": pid, "tid": tid, "comm": comm, "cpu": cpu})
    return found


def _iter_tasks():
    for d in os.listdir("/proc"):
        if not d.isdigit():
            continue
        pid = int(d)
        for tid in _tids(pid):
            try:
                with open(f"/proc/{pid}/task/{tid}/stat") as fh:
                    text = fh.read()
            except OSError:
                continue
            r = text.rfind(")")
            fields = text[r + 2:].split()
            yield pid, tid, fields[0], int(fields[36]), text[text.find("(") + 1:r]


@dataclass(eq=False)  # a Thread must stay hashable: threading keys _limbo by it
class _Sampler(threading.Thread):
    cpus: set[int]
    own: set[int]
    interval: float
    child_pid: int = 0
    samples: int = 0
    samples_foreign: int = 0
    foreign_examples: list = field(default_factory=list)
    max_threads: int = 0
    schedstat: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        threading.Thread.__init__(self, daemon=True)
        self._halt = threading.Event()

    def poll_tree(self) -> set[int]:
        tree = _tree(self.child_pid, _proc_table()) if self.child_pid else set()
        threads = 0
        for pid in tree:
            for tid in _tids(pid):
                threads += 1
                st = _schedstat(pid, tid)
                if st:
                    self.schedstat[tid] = st
        self.max_threads = max(self.max_threads, threads)
        return tree

    def run(self) -> None:
        while not self._halt.is_set():
            tree = self.poll_tree()
            foreign = _foreign_runnable(self.cpus, self.own | tree)
            self.samples += 1
            if foreign:
                self.samples_foreign += 1
                if len(self.foreign_examples) < 8:
                    self.foreign_examples.append(foreign[0])
            self._halt.wait(self.interval)

    def stop(self) -> None:
        self._halt.set()


def _jiffies(stat: dict[str, Any], cpus: set[int]) -> dict[str, int]:
    keys = ("user", "nice", "system", "idle", "iowait", "irq", "softirq", "steal")
    tot = dict.fromkeys(keys, 0)
    for c in cpus:
        row = stat.get(f"cpu{c}") or {}
        for k in keys:
            tot[k] += row.get(k, 0)
    return tot


def _psi_totals(p: dict[str, Any]) -> dict[str, Any]:
    return {res: (None if kinds is None else {k: v.get("total") for k, v in kinds.items()}) for res, kinds in p.items()}


def run_measured(argv: list[str], cpus: set[int], env: dict[str, str], cwd: str | None,
                 stdout_path: str, stderr_path: str, timeout_s: float | None = None,
                 memory_limit_mib: int | None = None, interval: float = 0.05) -> dict[str, Any]:
    """Run ``argv`` pinned to ``cpus``; return the execution block of a record."""
    psi0, stat0 = psi(), proc_stat()
    sampler = _Sampler(cpus=set(cpus), own={os.getpid()}, interval=interval)

    def preexec() -> None:
        os.sched_setaffinity(0, cpus)
        if memory_limit_mib:
            lim = memory_limit_mib * 1024 * 1024
            resource.setrlimit(resource.RLIMIT_AS, (lim, lim))

    with open(stdout_path, "wb") as out, open(stderr_path, "wb") as err:
        t0 = time.monotonic_ns()
        try:
            proc = subprocess.Popen(argv, cwd=cwd, env=env, stdout=out, stderr=err, stdin=subprocess.DEVNULL,
                                    preexec_fn=preexec, start_new_session=True)
        except OSError as exc:
            return {"argv": argv, "launch_error": str(exc), "wall_ns": 0, "exit": {"returncode": 127,
                    "signal": None, "timed_out": False}}
        sampler.child_pid = proc.pid
        sampler.start()
        try:
            observed_affinity = sorted(os.sched_getaffinity(proc.pid))
        except OSError:
            observed_affinity = None
        timed_out = False
        deadline = None if timeout_s is None else time.monotonic() + timeout_s
        while True:
            try:
                info = os.waitid(os.P_PID, proc.pid, os.WEXITED | os.WNOHANG | os.WNOWAIT)
            except ChildProcessError:
                info = None
            if info is not None:
                t1 = time.monotonic_ns()
                sampler.poll_tree()          # the exited root's schedstat is still readable
                _, status, ru = os.wait4(proc.pid, 0)
                break
            if deadline is not None and time.monotonic() > deadline:
                timed_out = True
                os.killpg(proc.pid, signal.SIGKILL)
                _, status, ru = os.wait4(proc.pid, 0)
                t1 = time.monotonic_ns()
                break
            time.sleep(0.002)
        proc.returncode = 0  # reaped by wait4
    sampler.stop()
    sampler.join()
    psi1, stat1 = psi(), proc_stat()
    code = os.waitstatus_to_exitcode(status)
    j0, j1 = _jiffies(stat0, cpus), _jiffies(stat1, cpus)
    p0, p1 = _psi_totals(psi0), _psi_totals(psi1)
    psi_delta: dict[str, Any] = {}
    for res, kinds in p1.items():
        if kinds is None or p0.get(res) is None:
            psi_delta[res] = None
        else:
            psi_delta[res] = {k: (kinds[k] - p0[res].get(k, 0.0)) for k in kinds}
    wall = t1 - t0
    run_ns = sum(v[0] for v in sampler.schedstat.values())
    wait_ns = sum(v[1] for v in sampler.schedstat.values())
    return {
        "argv": argv,
        "cwd": cwd,
        "pinned_cpus": sorted(cpus),
        "child_affinity_observed": observed_affinity,
        "exit": {"returncode": code if code >= 0 else None, "signal": -code if code < 0 else None,
                 "timed_out": timed_out},
        "wall_ns": wall,
        "rusage": {"user_s": ru.ru_utime, "sys_s": ru.ru_stime, "maxrss_kib": ru.ru_maxrss,
                   "minflt": ru.ru_minflt, "majflt": ru.ru_majflt, "nvcsw": ru.ru_nvcsw,
                   "nivcsw": ru.ru_nivcsw},
        "schedstat": {"run_ns": run_ns, "wait_ns": wait_ns,
                      "slices": sum(v[2] for v in sampler.schedstat.values()),
                      "threads_seen": len(sampler.schedstat),
                      "method": "per-thread /proc schedstat, sampled every interval and read from the unreaped root at exit; "
                                "a thread that exited between samples contributes its last sampled value"},
        "max_threads_observed": sampler.max_threads,
        "contention": {
            "interval_s": interval,
            "samples": sampler.samples,
            "samples_with_foreign_runnable": sampler.samples_foreign,
            "foreign_examples": sampler.foreign_examples,
            "pinned_cpu_jiffies_delta": {k: j1[k] - j0[k] for k in j0},
            "psi_total_delta_us": psi_delta,
            "procs_running_after": stat1.get("procs_running"),
        },
        "outputs": {
            "stdout": {"path": os.path.basename(stdout_path), "sha256": sha256_file(stdout_path),
                       "bytes": os.path.getsize(stdout_path)},
            "stderr": {"path": os.path.basename(stderr_path), "sha256": sha256_file(stderr_path),
                       "bytes": os.path.getsize(stderr_path)},
        },
    }
