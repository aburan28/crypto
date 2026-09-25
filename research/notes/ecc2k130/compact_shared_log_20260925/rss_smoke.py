#!/usr/bin/env python3
"""Tiny non-producer check of live psutil RSS-cap and wait4 peak handling."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

from run import monitor_process, write_json


def one(name: str, code: str, cap_bytes: int) -> dict:
    started = time.monotonic_ns()
    process = subprocess.Popen([sys.executable, "-c", code],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    status, usage, termination, sampled, wait4_peak, samples = monitor_process(
        process, started, 5, cap_bytes)
    return {"name": name, "exit_code": os.waitstatus_to_exitcode(status),
            "termination": termination, "sampled_peak_rss_bytes": sampled,
            "wait4_peak_rss_bytes": wait4_peak, "rss_sample_count": samples,
            "child_cpu_s": usage.ru_utime+usage.ru_stime,
            "wall_ms": (time.monotonic_ns()-started)/1e6,
            "cap_bytes": cap_bytes}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    normal = one("below_cap", "import time; time.sleep(0.2)", 64*1024*1024)
    capped = one("cross_cap", "import time; x=bytearray(64*1024*1024); time.sleep(1)",
                 32*1024*1024)
    assert normal["exit_code"] == 0 and normal["termination"] is None
    assert normal["rss_sample_count"] >= 1
    assert capped["termination"] == "RSS_CAP" and capped["exit_code"] != 0
    assert capped["sampled_peak_rss_bytes"] >= capped["cap_bytes"]
    report = {"classification": "RSS_CAP_SMOKE_PASS", "platform": platform.platform(),
              "python": sys.version, "normal": normal, "capped": capped}
    write_json(args.out, report)
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
