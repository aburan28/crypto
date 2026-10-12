#!/usr/bin/env python3
"""Record a quiet-host admission check before the frozen paired SAT panel."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import shutil
import sys

import psutil
import run

LOAD_1_MAX = 2.0
LOAD_5_MAX = 3.0
MIN_FREE_BYTES = 4 * 1024**3
SOLVER_NAMES = {"cryptominisat5", "kissat", "cadical"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit("refusing to overwrite a host admission receipt")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    receipt = {"schema": "n13-paired-host-admission-v1",
               "utc": datetime.now(timezone.utc).isoformat(),
               "decision": "HELD", "runner": str(Path(sys.executable).resolve()),
               "python": platform.python_version(),
               "psutil": psutil.__version__,
               "limits": {"load_1_max": LOAD_1_MAX, "load_5_max": LOAD_5_MAX,
                          "min_free_bytes": MIN_FREE_BYTES},
               "frozen_preflight": False, "process_tree_monitor": False}
    try:
        frozen, *_ = run.preflight()
        run.assert_monitor_available()
        receipt["frozen_preflight"] = True
        receipt["process_tree_monitor"] = True
        receipt["freeze_sha256"] = run.verify.sha(run.HERE / "FROZEN.json")
        receipt["source_commit"] = run.verify.source_ancestry(frozen)
        receipt["logical_cpus"] = os.cpu_count()
        receipt["load_1_5_15"] = list(os.getloadavg())
        receipt["disk_free_bytes"] = shutil.disk_usage("/private/tmp").free
        competitors = []
        for process in psutil.process_iter(["pid", "name"]):
            name = process.info["name"]
            if name in SOLVER_NAMES:
                competitors.append({"pid": process.info["pid"], "name": name})
        receipt["competing_solvers"] = sorted(competitors, key=lambda item: item["pid"])
        receipt["decision"] = ("ADMITTED" if
            receipt["load_1_5_15"][0] <= LOAD_1_MAX and
            receipt["load_1_5_15"][1] <= LOAD_5_MAX and
            receipt["disk_free_bytes"] >= MIN_FREE_BYTES and
            not receipt["competing_solvers"] else "HELD")
    except (OSError, psutil.Error, AssertionError, RuntimeError) as exc:
        receipt["error"] = f"{type(exc).__name__}: {exc}"
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"decision": receipt["decision"], "out": str(args.out)}, sort_keys=True))
    return 0 if receipt["decision"] == "ADMITTED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
