#!/usr/bin/env python3
"""Harmless runner controls: no SAT input, producer, or measured panel."""
from __future__ import annotations

import json
import signal
import sys
import tempfile
import time
from pathlib import Path

import run


def sleeper(root: Path, label: str, wall: float, rss: int, deadline=None):
    return run.execute([sys.executable, "-c", "import time; time.sleep(2)"],
                       root / f"{label}.stdout", root / f"{label}.stderr",
                       wall, rss, deadline)


def main():
    with tempfile.TemporaryDirectory(prefix="oaware-runner-control-") as tmp:
        root = Path(tmp)
        previous_here, previous_argv = run.HERE, sys.argv
        try:
            run.HERE = root / "deliberately-missing-freeze"
            sys.argv = [str(previous_here / "run.py"), "smoke", str(root / "failed-preflight")]
            try:
                run.main()
            except FileNotFoundError:
                pass
            else:
                raise AssertionError("missing freeze unexpectedly passed preflight")
            receipt = json.loads((root / "failed-preflight/receipt.json").read_text())
            assert receipt["decision"] == "CENSORED_OR_FAILED"
            assert receipt["freeze_sha256"] is None and receipt["preflight_wall_seconds"] >= 0
        finally:
            run.HERE, sys.argv = previous_here, previous_argv

        # Exercise child-kill accounting with deterministic empty-child enumeration.
        # This does not certify that the real host permits process-tree monitoring.
        real_children = run.psutil.Process.children
        run.psutil.Process.children = lambda self, recursive=False: []
        try:
            assert sleeper(root, "wall", .05, 1 << 30)["stop_reason"] == "wall_cap"
            assert sleeper(root, "rss", 1, 1)["stop_reason"] == "sampled_rss_cap"
            previous_handler = signal.signal(signal.SIGALRM, run.deadline_signal)
            try:
                signal.setitimer(signal.ITIMER_REAL, .05)
                row = sleeper(root, "portfolio", 1, 1 << 30, time.perf_counter() + 1)
                assert row["stop_reason"] == "portfolio_cap"
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
                signal.signal(signal.SIGALRM, previous_handler)
        finally:
            run.psutil.Process.children = real_children
    print("runner failure archive and mocked-process-tree cap controls PASS; no solver ran")


if __name__ == "__main__":
    main()
