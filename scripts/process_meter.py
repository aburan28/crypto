#!/usr/bin/env python3
"""Execute one command and write exact child CPU/max-RSS metrics as JSON."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import resource
import signal
import subprocess
import threading
import time


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cwd", type=Path, required=True)
    parser.add_argument("--timeout", type=float, required=True)
    parser.add_argument("--stdout", type=Path, required=True)
    parser.add_argument("--stderr", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--exclusive-create", action="store_true")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command and args.command[0] == "--" else args.command
    if not command:
        parser.error("a command is required after --")

    outputs = (args.stdout, args.stderr, args.metrics)
    if len({path.resolve() for path in outputs}) != len(outputs):
        parser.error("stdout, stderr, and metrics paths must be distinct")
    reserved: list[tuple[Path, int]] = []
    try:
        if args.exclusive_create:
            flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
            if hasattr(os, "O_NOFOLLOW"):
                flags |= os.O_NOFOLLOW
            for path in outputs:
                reserved.append((path, os.open(path, flags, 0o644)))
        else:
            for path in outputs:
                reserved.append(
                    (path, os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644))
                )
    except OSError as error:
        for created_path, descriptor in reserved:
            os.close(descriptor)
            if args.exclusive_create:
                try:
                    created_path.unlink()
                except OSError:
                    pass
        parser.error(f"cannot reserve output artifacts: {error}")

    descriptors = {path: descriptor for path, descriptor in reserved}
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.perf_counter()
    timed_out = threading.Event()
    with os.fdopen(descriptors[args.stdout], "w") as stdout, os.fdopen(descriptors[args.stderr], "w") as stderr:
        process = subprocess.Popen(
            command,
            cwd=args.cwd,
            text=True,
            stdout=stdout,
            stderr=stderr,
            start_new_session=True,
        )
        def kill_if_running(sig: int) -> None:
            if process.poll() is None:
                os.killpg(process.pid, sig)

        def expire() -> None:
            timed_out.set()
            kill_if_running(signal.SIGTERM)
            hard_kill = threading.Timer(5, kill_if_running, args=(signal.SIGKILL,))
            hard_kill.daemon = True
            hard_kill.start()

        watchdog = threading.Timer(args.timeout, expire)
        watchdog.daemon = True
        watchdog.start()
        process.wait()
        watchdog.cancel()
        orphan_group_terminated = False
        try:
            # The leader may already be reaped while descendants remain in
            # its process group.  poll() cannot detect those descendants.
            os.killpg(process.pid, 0)
        except ProcessLookupError:
            pass
        else:
            orphan_group_terminated = True
            os.killpg(process.pid, signal.SIGTERM)
            time.sleep(0.05)
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    wall = time.perf_counter() - started
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    user = after.ru_utime - before.ru_utime
    system = after.ru_stime - before.ru_stime
    # Darwin reports bytes; Linux and the BSDs report KiB.
    rss_scale = 1 if platform.system() == "Darwin" else 1024
    record = {
        "command": command,
        "returncode": process.returncode,
        "watchdog_seconds": args.timeout,
        "timed_out": timed_out.is_set(),
        "orphan_group_terminated": orphan_group_terminated,
        "metrics": {
            "wall_seconds": wall,
            "user_seconds": user,
            "system_seconds": system,
            "total_core_seconds": user + system,
            "single_core_seconds": user + system,
            "peak_rss_bytes": int(after.ru_maxrss * rss_scale),
            "meter": "fresh-process getrusage(RUSAGE_CHILDREN)",
        },
    }
    with os.fdopen(descriptors[args.metrics], "w") as metrics:
        metrics.write(json.dumps(record, indent=2) + "\n")
        metrics.flush()
        os.fsync(metrics.fileno())


if __name__ == "__main__":
    main()
