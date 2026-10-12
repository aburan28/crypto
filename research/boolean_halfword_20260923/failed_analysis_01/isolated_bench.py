#!/usr/bin/env python3
"""Run a benchmark on a reserved, pinned CPU and record what else was running.

Wall time is only evidence when the run had the core to itself.  Pinning the
benchmark (``taskset``) is not enough on its own: it keeps the benchmark on one
core but does not keep anything else off that core, and it does nothing about
builds or other jobs on the neighbouring cores that share the cache and the
memory bus.  This tool does all three, and records the evidence:

1. **One benchmark at a time.**  It takes an exclusive ``flock`` on
   ``--lock`` (default ``/tmp/crypto-bench.lock``).  Heavy non-benchmark work
   (``cargo build``, test suites, audit cells) should run through
   ``isolated_bench.py busy -- CMD``, which takes the same lock, so a build
   cannot start in the middle of a timed stage.
2. **A quiet machine before starting.**  It samples every other process's CPU
   use for ``--settle`` seconds and refuses to start if they use more than
   ``--max-other-cpu`` CPUs' worth between them, or if the 10-second CPU or
   memory pressure (PSI) is above ``--max-psi``.
3. **A reserved core.**  It moves every other thread it is allowed to move off
   the benchmark CPUs (``--cpus``), pins the benchmark there, and restores the
   other threads' affinity afterwards.  On a host with SMT it refuses a CPU
   whose sibling is not reserved too.  A process started while a run had
   moved its parent off those CPUs inherits the narrowed mask, and no restore
   reaches it; so when this process's own mask lacks the CPUs asked for, it
   widens the mask to them where the system allows, and the record says so
   (``affinity_widened_from``).
4. **A record of the conditions.**  Per run: wall time, user and system CPU,
   voluntary and involuntary context switches, page faults, load average and
   PSI before and after, and the CPU time every other process used while the
   benchmark ran.  A run where other processes used more than
   ``--max-other-cpu`` is marked ``contended``.

What it cannot do: it cannot see or stop other tenants of a virtual machine's
host, and it cannot fix the CPU frequency.  That residual noise is what an A/A
run measures, and the A/A spread is reported beside every A/B difference
(AGENTS.md section 10).

Usage::

    tools/isolated_bench.py run --cpus 3 --out rec.jsonl -- ./worker < job.json
    tools/isolated_bench.py reserve --cpus 3 --out cond.json -- python3 tournament.py run ... --cpu 3
    tools/isolated_bench.py busy -- cargo build --release

``run`` pins the command itself.  ``reserve`` is for a harness that pins its
own children (the tournament evaluator): it takes the lock, reserves the CPUs
and monitors the machine for the whole command, but leaves pinning to it.
"""

from __future__ import annotations

import argparse
import contextlib
import fcntl
import json
import os
import platform
import subprocess
import sys
import threading
import time
from pathlib import Path

DEFAULT_LOCK = '/tmp/crypto-bench.lock'
TICK = os.sysconf('SC_CLK_TCK')


def parse_cpus(text: str) -> set[int]:
    cpus: set[int] = set()
    for part in text.split(','):
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            lo, hi = part.split('-')
            cpus.update(range(int(lo), int(hi) + 1))
        else:
            cpus.add(int(part))
    if not cpus:
        raise ValueError('no CPUs given')
    return cpus


def smt_siblings(cpu: int) -> set[int]:
    path = Path(f'/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list')
    try:
        return parse_cpus(path.read_text())
    except (OSError, ValueError):
        return {cpu}


def check_cpus(cpus: set[int]) -> list[int] | None:
    """Refuse a reservation that cannot work.  Return the inherited mask if
    it lacked some of ``cpus`` and was widened to them, else ``None``.

    An isolated run moves every other thread off its CPUs and restores the
    threads it moved; a process forked meanwhile inherits the narrowed mask
    and keeps it.  A harness started that way would be refused here on every
    later run, so the mask is widened instead, where the system allows it.
    """
    allowed = os.sched_getaffinity(0)
    inherited = None
    if cpus - allowed:
        with contextlib.suppress(OSError):
            os.sched_setaffinity(0, allowed | cpus)
        now = os.sched_getaffinity(0)
        missing = cpus - now
        if missing:
            raise SystemExit(f'CPUs {sorted(missing)} are outside this process affinity {sorted(allowed)}')
        inherited, allowed = sorted(allowed), now
    if cpus == allowed:
        raise SystemExit('reserving every CPU leaves nowhere to move other work; leave at least one free')
    for cpu in cpus:
        stray = smt_siblings(cpu) - cpus
        if stray:
            raise SystemExit(f'CPU {cpu} shares a core with {sorted(stray)}; reserve the siblings too')
    return inherited


def psi(kind: str) -> dict | None:
    try:
        lines = Path(f'/proc/pressure/{kind}').read_text().split('\n')
    except OSError:
        return None
    out = {}
    for line in lines:
        if not line:
            continue
        name, *fields = line.split()
        out[name] = {k: float(v) for k, v in (f.split('=') for f in fields)}
    return out


def loadavg() -> list[float]:
    return [float(x) for x in Path('/proc/loadavg').read_text().split()[:3]]


def conditions() -> dict:
    return {'unix_time': time.time(), 'loadavg': loadavg(),
            'psi_cpu': psi('cpu'), 'psi_memory': psi('memory')}


def descendants(root: int, snapshot: dict[int, tuple[str, int, int]]) -> set[int]:
    """Find children in the same /proc scan used for CPU accounting.

    A short-lived benchmark can exit between a CPU scan and a second process
    tree scan. In that case its charged ticks look like unrelated load.
    """
    children: dict[int, list[int]] = {}
    for pid, (_, _, ppid) in snapshot.items():
        children.setdefault(ppid, []).append(pid)
    out, stack = set(), [root]
    while stack:
        pid = stack.pop()
        for child in children.get(pid, []):
            if child not in out:
                out.add(child)
                stack.append(child)
    return out


def cpu_snapshot() -> dict[int, tuple[str, int, int]]:
    """Read each process's command, CPU ticks and parent from one stat entry."""
    out = {}
    for entry in Path('/proc').iterdir():
        if not entry.name.isdigit():
            continue
        try:
            stat = (entry / 'stat').read_text()
        except OSError:
            continue
        name = stat[stat.index('(') + 1:stat.rindex(')')]
        fields = stat[stat.rindex(')') + 2:].split()
        out[int(entry.name)] = (name, int(fields[11]) + int(fields[12]), int(fields[1]))
    return out


def cpu_ticks() -> dict[int, tuple[str, int]]:
    """utime+stime of every process, in clock ticks, with its command name."""
    return {pid: (name, ticks) for pid, (name, ticks, _) in cpu_snapshot().items()}


def other_use(before: dict, after: dict, exclude: set[int]) -> dict:
    """CPU seconds used by processes outside ``exclude`` between two samples."""
    used = {}
    for pid, (name, ticks) in after.items():
        if pid in exclude:
            continue
        delta = ticks - before.get(pid, (name, 0))[1]
        if delta > 0:
            used[f'{name}[{pid}]'] = delta / TICK
    return dict(sorted(used.items(), key=lambda kv: -kv[1]))


def settle(seconds: float, exclude: set[int]) -> dict:
    before = cpu_ticks()
    time.sleep(seconds)
    used = other_use(before, cpu_ticks(), exclude)
    return {'seconds': seconds, 'other_cpu_seconds': sum(used.values()), 'top': dict(list(used.items())[:8])}


def evict(cpus: set[int], keep: set[int]) -> dict[int, set[int]]:
    """Move every other thread off ``cpus``; return the affinities to restore."""
    moved = {}
    for task_dir in Path('/proc').glob('[0-9]*/task/[0-9]*'):
        tid = int(task_dir.name)
        pid = int(task_dir.parent.parent.name)
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


def restore(moved: dict[int, set[int]]) -> None:
    for tid, affinity in moved.items():
        with contextlib.suppress(OSError):
            os.sched_setaffinity(tid, affinity)


def unmovable_on(cpus: set[int], keep: set[int]) -> dict:
    """Threads still allowed on the reserved CPUs after eviction.

    Per-CPU kernel threads cannot be moved and are only counted; a user thread
    left here refused the move and is named, since it can still take the core.
    """
    user, kernel = [], 0
    for task_dir in Path('/proc').glob('[0-9]*/task/[0-9]*'):
        pid = int(task_dir.parent.parent.name)
        if pid in keep:
            continue
        try:
            if not os.sched_getaffinity(int(task_dir.name)) & cpus:
                continue
            if (task_dir.parent.parent / 'cmdline').read_bytes():
                user.append(f"{(task_dir / 'comm').read_text().strip()}[{task_dir.name}]")
            else:
                kernel += 1
        except OSError:
            continue
    return {'user_threads': user, 'kernel_threads': kernel}


@contextlib.contextmanager
def locked(path: str, wait: bool):
    handle = open(path, 'a+')
    try:
        flags = fcntl.LOCK_EX if wait else fcntl.LOCK_EX | fcntl.LOCK_NB
        try:
            fcntl.flock(handle, flags)
        except BlockingIOError:
            raise SystemExit(f'another benchmark or busy job holds {path}')
        yield
    finally:
        handle.close()


def host() -> dict:
    model = ''
    with contextlib.suppress(OSError):
        for line in Path('/proc/cpuinfo').read_text().split('\n'):
            if line.startswith('model name'):
                model = line.split(':', 1)[1].strip()
                break
    return {'cpu_model': model, 'logical_cpus': os.cpu_count(), 'kernel': platform.release(),
            'machine': platform.machine(), 'python': platform.python_version()}


def preflight(args, exclude: set[int]) -> dict:
    quiet = settle(args.settle, exclude)
    now = conditions()
    worst = max((p['some']['avg10'] for p in (now['psi_cpu'], now['psi_memory']) if p), default=0.0)
    budget = args.max_other_cpu * args.settle
    if quiet['other_cpu_seconds'] > budget or worst > args.max_psi:
        raise SystemExit('machine is busy: other processes used '
                         f"{quiet['other_cpu_seconds']:.2f} CPU s in {args.settle} s (limit {budget:.2f}), "
                         f"PSI some avg10 {worst} (limit {args.max_psi}); top: {quiet['top']}")
    return {'settle': quiet, 'conditions': now}


def run_pinned(args, cpus: set[int], command: list[str], stdin, stdout=None) -> dict:
    exclude = {os.getpid()}
    before_ticks = cpu_ticks()
    before = conditions()
    inherited = os.sched_getaffinity(0)
    os.sched_setaffinity(0, cpus)
    start = time.monotonic()
    try:
        child = subprocess.Popen(command, stdin=stdin, stdout=stdout)
    finally:
        os.sched_setaffinity(0, inherited)
    _, status, usage = os.wait4(child.pid, 0)
    wall = time.monotonic() - start
    after = conditions()
    others = other_use(before_ticks, cpu_ticks(), exclude | {child.pid})
    other_seconds = sum(others.values())
    return {'command': command, 'cpus': sorted(cpus), 'exit_status': os.waitstatus_to_exitcode(status),
            'wall_seconds': wall, 'user_seconds': usage.ru_utime, 'system_seconds': usage.ru_stime,
            'voluntary_switches': usage.ru_nvcsw, 'involuntary_switches': usage.ru_nivcsw,
            'minor_faults': usage.ru_minflt, 'major_faults': usage.ru_majflt,
            'max_rss_kib': usage.ru_maxrss, 'before': before, 'after': after,
            'other_cpu_seconds': other_seconds, 'other_top': dict(list(others.items())[:8]),
            'contended': other_seconds > args.max_other_cpu * wall}


class Monitor(threading.Thread):
    """Sample other processes' CPU use every ``period`` seconds while a harness runs."""

    def __init__(self, period: float, exclude_root: int, threshold: float):
        super().__init__(daemon=True)
        self.period, self.root, self.threshold = period, exclude_root, threshold
        self.samples: list[dict] = []
        self.stop = threading.Event()

    def run(self) -> None:
        mine = {os.getpid()}
        previous = cpu_ticks()
        while not self.stop.wait(self.period):
            snapshot = cpu_snapshot()
            now = {pid: (name, ticks) for pid, (name, ticks, _) in snapshot.items()}
            used = other_use(previous, now, mine | descendants(self.root, snapshot) | {self.root})
            total = sum(used.values())
            self.samples.append({'unix_time': time.time(), 'other_cpu_seconds': total,
                                 'contended': total > self.threshold * self.period,
                                 'top': dict(list(used.items())[:5]), 'loadavg': loadavg()})
            previous = now


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='mode', required=True)
    for name in ('run', 'reserve', 'busy'):
        p = sub.add_parser(name)
        p.add_argument('--lock', default=DEFAULT_LOCK)
        p.add_argument('--wait', action='store_true', help='wait for the lock instead of failing')
        if name != 'busy':
            p.add_argument('--cpus', required=True, help='benchmark CPUs, e.g. 3 or 2-3')
            p.add_argument('--out', required=True, help='JSON lines file the record is appended to')
            p.add_argument('--settle', type=float, default=2.0)
            p.add_argument('--max-other-cpu', type=float, default=0.10,
                           help='other processes may use at most this many CPUs on average')
            p.add_argument('--max-psi', type=float, default=5.0)
            p.add_argument('--label', default='')
        if name == 'run':
            p.add_argument('--stdin', help='file fed to the command on standard input')
        if name == 'reserve':
            p.add_argument('--period', type=float, default=5.0)
        p.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        parser.error('no command after --')

    if args.mode == 'busy':
        with locked(args.lock, wait=True):
            return subprocess.call(command)

    cpus = parse_cpus(args.cpus)
    inherited = check_cpus(cpus)
    with locked(args.lock, wait=args.wait):
        # Only this process is exempt.  The shell and agent that launched it
        # are moved off the reserved CPUs and charged as contention like
        # anything else: an agent busy during a timed run is still noise.
        mine = {os.getpid()}
        pre = preflight(args, mine)
        moved = evict(cpus, mine)
        try:
            record = {'schema': 'isolated-bench/1', 'mode': args.mode, 'label': args.label,
                      'host': host(), 'reserved_cpus': sorted(cpus), 'threads_moved': len(moved),
                      'left_on_reserved': unmovable_on(cpus, mine), 'preflight': pre}
            if inherited is not None:
                record['affinity_widened_from'] = inherited
            if args.mode == 'run':
                stdin = open(args.stdin, 'rb') if args.stdin else None
                try:
                    record['run'] = run_pinned(args, cpus, command, stdin)
                finally:
                    if stdin:
                        stdin.close()
                code = record['run']['exit_status']
            else:
                child = subprocess.Popen(command)
                monitor = Monitor(args.period, child.pid, args.max_other_cpu)
                monitor.start()
                code = child.wait()
                monitor.stop.set()
                monitor.join()
                record.update(command=command, exit_status=code, samples=monitor.samples,
                              contended_samples=sum(s['contended'] for s in monitor.samples))
        finally:
            restore(moved)
    with open(args.out, 'a') as handle:
        handle.write(json.dumps(record, sort_keys=True) + '\n')
    return code


if __name__ == '__main__':
    sys.exit(main())
