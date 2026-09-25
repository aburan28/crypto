#!/usr/bin/env python3
"""Cold child supervision with enforced address-space/file caps and watchdogs."""
from __future__ import annotations

import hashlib
import os
import resource
import signal
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def _limits(rss_cap: int, file_cap: int | None):
    def apply() -> None:
        resource.setrlimit(resource.RLIMIT_AS, (rss_cap, rss_cap))
        if file_cap is not None:
            resource.setrlimit(resource.RLIMIT_FSIZE, (file_cap, file_cap))
    return apply


def run_child(command: list[str], *, cwd: Path, stdout: Path, stderr: Path,
              wall_cap: float, rss_cap: int, watched_file: Path | None = None,
              watched_file_cap: int | None = None, file_cap: int | None = None) -> dict:
    """Return a raw receipt; a nonzero status or cap is never silently accepted."""
    stdout.parent.mkdir(parents=True, exist_ok=True)
    stderr.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    utc = datetime.now(timezone.utc).isoformat()
    peak = 0
    stop_reason = None
    with stdout.open('wb') as out_stream, stderr.open('wb') as err_stream:
        child = subprocess.Popen(command, cwd=cwd, stdout=out_stream, stderr=err_stream,
                                 start_new_session=True,
                                 preexec_fn=_limits(rss_cap, file_cap))
        while child.poll() is None:
            # POSIX process-group sample; RLIMIT_AS remains the hard per-process
            # address-space cap even if a short process escapes a sample.
            sampled = subprocess.run(['ps', '-axo', 'pgid=,rss='],
                                     capture_output=True, text=True, check=False)
            if sampled.returncode == 0:
                rss = 0
                for row in sampled.stdout.splitlines():
                    pieces = row.split()
                    if len(pieces) == 2 and pieces[0].isdigit() and pieces[1].isdigit():
                        if int(pieces[0]) == child.pid:
                            rss += int(pieces[1]) * 1024
                peak = max(peak, rss)
            elapsed = time.monotonic() - started
            if elapsed > wall_cap:
                stop_reason = 'WALL_CAP'
            elif peak > rss_cap:
                stop_reason = 'RSS_CAP'
            elif (watched_file is not None and watched_file_cap is not None
                  and watched_file.exists() and watched_file.stat().st_size > watched_file_cap):
                stop_reason = 'FILE_CAP'
            if stop_reason:
                try:
                    os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                break
            time.sleep(0.05)
        code = child.wait()
    wall = time.monotonic() - started
    return {'argv': command, 'cwd': str(cwd), 'utc_start': utc,
            'exit_code': code, 'stop_reason': stop_reason,
            'wall_seconds': wall, 'sampled_peak_rss_bytes': peak,
            'rss_address_space_cap_bytes': rss_cap,
            'file_size_cap_bytes': file_cap,
            'stdout_sha256': sha(stdout), 'stderr_sha256': sha(stderr),
            'watched_file_sha256': sha(watched_file) if watched_file else None,
            'watched_file_bytes': watched_file.stat().st_size if watched_file and watched_file.is_file() else None}
