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

Nothing the child starts outlives its record.  The runner makes itself a
child subreaper, so a descendant that daemonises or calls ``setsid`` is
reparented to it rather than to init; when the root exits (or times out, or
the runner is interrupted) every remaining descendant is killed and reaped
before the output files are hashed, and the record counts them.

Every sample also reads the CPU mask and the last CPU of every task in the
tree, so a producer that re-pins itself (``taskset``, an OpenMP
``proc_bind``) is caught, not just a child whose first mask was wrong.

The environment the child sees is built, not inherited: a fixed allowlist of
variables needed to find programs, plus exactly what the spec declares.  An
engine knob exported in the operator's shell (KIC_*, IC_*, GAUDRY_*, ...)
therefore cannot change a measured run without appearing in its spec.
"""
from __future__ import annotations

import ctypes
import hashlib
import os
import resource
import shutil
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
    passed = [k for k in PASSTHROUGH if k in os.environ]
    # Values are hashed, not stored: PATH and LD_LIBRARY_PATH decide which
    # binary and libraries run, and a hash shows when two runs differed
    # without writing a home directory into the evidence.
    return env, {"declared": dict(sorted((declared or {}).items())), "defaults": DEFAULT_ENV,
                 "passthrough": passed,
                 "passthrough_sha256": {k: hashlib.sha256(os.environ[k].encode()).hexdigest() for k in passed},
                 "dropped_engine_knobs": dropped}


def resolve_argv0(argv: list[str], env: dict[str, str], cwd: str | None) -> dict[str, Any]:
    """The file argv[0] names, as the child's PATH resolves it, and its hash."""
    prog = argv[0] if argv else ""
    path = prog if os.sep in prog else shutil.which(prog, path=env.get("PATH"))
    if path and not os.path.isabs(path) and cwd:
        path = os.path.join(cwd, path)
    real = os.path.realpath(path) if path else None
    return {"name": prog, "resolved": real, "sha256": sha256_file(real) if real else None}


_PR_SET_CHILD_SUBREAPER = 36


def _become_subreaper() -> bool:
    try:
        libc = ctypes.CDLL(None, use_errno=True)
        return libc.prctl(_PR_SET_CHILD_SUBREAPER, 1, 0, 0, 0) == 0
    except (OSError, AttributeError):
        return False


def _stat_fields(path: str) -> tuple[list[str], str] | None:
    try:
        with open(path) as fh:
            text = fh.read()
    except OSError:
        return None
    r = text.rfind(")")
    return text[r + 2:].split(), text[text.find("(") + 1:r]


def _descendants(root: int, me: int) -> set[int]:
    """Live processes started by the measured child: its session (the child is
    a session leader), anything still parented below it, and anything
    reparented to this runner as subreaper other than the root itself."""
    table: dict[int, tuple[int, int]] = {}
    for d in os.listdir("/proc"):
        if not d.isdigit():
            continue
        got = _stat_fields(f"/proc/{d}/stat")
        if got:
            f = got[0]
            table[int(d)] = (int(f[1]), int(f[3]))  # ppid, session
    out = {pid for pid, (ppid, sid) in table.items() if sid == root or ppid == me}
    out |= _tree(root, {pid: (ppid, "", 0, "") for pid, (ppid, _sid) in table.items()})
    out.discard(root)
    out.discard(me)
    return {p for p in out if p in table}


def _kill_descendants(root: int, me: int) -> int:
    """SIGKILL and reap every descendant of ``root``; return how many there were."""
    killed: set[int] = set()
    for _ in range(5):
        left = _descendants(root, me) - killed
        if not left:
            break
        for pid in left:
            try:
                os.kill(pid, signal.SIGKILL)
                killed.add(pid)
            except ProcessLookupError:
                pass
        for pid in left:
            try:
                os.waitpid(pid, 0)  # ours only when reparented to this subreaper
            except ChildProcessError:
                pass
        time.sleep(0.01)
    return len(killed)


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


def _foreign_runnable(cpus: set[int], exclude: set[int], session: int) -> list[dict[str, Any]]:
    # Membership is decided in this same scan by session id, so a child the
    # measured tree forked a moment ago is ours, not foreign.
    found = []
    for pid, tid, state, cpu, comm, sid in _iter_tasks():
        if state == "R" and cpu in cpus and pid not in exclude and sid != session:
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
            yield pid, tid, fields[0], int(fields[36]), text[text.find("(") + 1:r], int(fields[3])


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
    cpus_seen: set = field(default_factory=set)
    outside: int = 0
    outside_examples: list = field(default_factory=list)

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
                # Where the task may run and where it last ran: a task that
                # re-pinned itself shows up here even if the root's mask held.
                try:
                    mask = os.sched_getaffinity(tid)
                except OSError:
                    mask = None
                got = _stat_fields(f"/proc/{pid}/task/{tid}/stat")
                last = int(got[0][36]) if got else None
                if last is not None:
                    self.cpus_seen.add(last)
                if (mask is not None and not mask <= self.cpus) or (last is not None and last not in self.cpus):
                    self.outside += 1
                    if len(self.outside_examples) < 8:
                        self.outside_examples.append({"pid": pid, "tid": tid, "mask": sorted(mask) if mask else None,
                                                      "last_cpu": last, "comm": got[1] if got else None})
        self.max_threads = max(self.max_threads, threads)
        return tree

    def run(self) -> None:
        while not self._halt.is_set():
            tree = self.poll_tree()
            foreign = _foreign_runnable(self.cpus, self.own | tree, self.child_pid)
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


def _empty_execution(argv: list[str], cwd: str | None, cpus: set[int], interval: float,
                     stdout_path: str, stderr_path: str) -> dict[str, Any]:
    """The full execution shape for a run that never started: every
    observation unknown, the (empty) output files still hashed."""
    return {
        "argv": argv, "cwd": cwd, "pinned_cpus": sorted(cpus), "child_affinity_observed": None,
        "tree_affinity": {"tasks_outside": 0, "cpus_seen": [], "examples": []},
        "descendants_killed": 0,
        "exit": {"returncode": 127, "signal": None, "timed_out": False},
        "wall_ns": 0, "rusage": None,
        "schedstat": {"run_ns": None, "wait_ns": None, "slices": None, "threads_seen": 0,
                      "method": "not measured: the command did not start"},
        "max_threads_observed": None,
        "contention": {"interval_s": interval, "samples": 0, "samples_with_foreign_runnable": 0,
                       "foreign_examples": [], "pinned_cpu_jiffies_delta": None, "psi_total_delta_us": None,
                       "procs_running_after": None},
        "outputs": _outputs(stdout_path, stderr_path),
    }


def _outputs(stdout_path: str, stderr_path: str) -> dict[str, Any]:
    out: dict[str, Any] = {name: {"path": os.path.basename(p), "sha256": sha256_file(p), "bytes": os.path.getsize(p)}
                           for name, p in (("stdout", stdout_path), ("stderr", stderr_path))}
    out["files"] = {}  # the session lists the producer's other files here
    return out


def run_measured(argv: list[str], cpus: set[int], env: dict[str, str], cwd: str | None,
                 stdout_path: str, stderr_path: str, timeout_s: float | None = None,
                 memory_limit_mib: int | None = None, interval: float = 0.05) -> dict[str, Any]:
    """Run ``argv`` pinned to ``cpus``; return the execution block of a record."""
    me = os.getpid()
    subreaper = _become_subreaper()
    psi0, stat0 = psi(), proc_stat()
    sampler = _Sampler(cpus=set(cpus), own={me}, interval=interval)

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
            out.close()
            err.close()
            rec = _empty_execution(argv, cwd, cpus, interval, stdout_path, stderr_path)
            rec["launch_error"] = str(exc)
            return rec
        sampler.child_pid = proc.pid
        sampler.start()
        try:
            observed_affinity = sorted(os.sched_getaffinity(proc.pid))
        except OSError:
            observed_affinity = None
        timed_out = False
        reaped = False
        deadline = None if timeout_s is None else time.monotonic() + timeout_s
        try:
            while True:
                try:
                    info = os.waitid(os.P_PID, proc.pid, os.WEXITED | os.WNOHANG | os.WNOWAIT)
                except ChildProcessError:
                    info = None
                if info is not None:
                    t1 = time.monotonic_ns()
                    sampler.poll_tree()          # the exited root's schedstat is still readable
                    _, status, ru = os.wait4(proc.pid, 0)
                    reaped = True
                    break
                if deadline is not None and time.monotonic() > deadline:
                    timed_out = True
                    os.killpg(proc.pid, signal.SIGKILL)
                    _, status, ru = os.wait4(proc.pid, 0)
                    reaped = True
                    t1 = time.monotonic_ns()
                    break
                time.sleep(0.002)
        finally:
            # Interrupted (^C, an exception): take the whole tree down with us.
            if not reaped:
                sampler.stop()
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                _kill_descendants(proc.pid, me)
                try:
                    os.waitpid(proc.pid, 0)
                except ChildProcessError:
                    pass
        proc.returncode = 0  # reaped by wait4
        sampler.stop()
        sampler.join()
        # Whatever the root left behind dies before its output is hashed.
        leftover = _kill_descendants(proc.pid, me)
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
        "tree_affinity": {"tasks_outside": sampler.outside, "cpus_seen": sorted(sampler.cpus_seen),
                          "examples": sampler.outside_examples},
        "descendants_killed": leftover,
        "subreaper": subreaper,
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
                                "a thread that exited between samples contributes its last sampled value, and one that "
                                "lived between two samples is missed (the gate compares run_ns with rusage CPU)"},
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
        "outputs": _outputs(stdout_path, stderr_path),
    }
