#!/usr/bin/env python3
"""Select a reservable CPU set and the CPUs the F5 harness will pin."""

import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from isolated_bench import smt_siblings  # noqa: E402


def choose(threads):
    allowed = set(os.sched_getaffinity(0))
    groups = []
    for cpu in sorted(allowed):
        siblings = smt_siblings(cpu)
        if siblings <= allowed and siblings not in groups:
            groups.append(siblings)
    for start in range(len(groups)):
        reserved = set()
        for group in groups[start:]:
            reserved |= group
            if len(reserved) >= threads:
                if len(reserved) < len(allowed):
                    return {
                        "status": "ok",
                        "allowed": sorted(allowed),
                        "reserved": sorted(reserved),
                        "pin": sorted(reserved)[:threads],
                        "threads": threads,
                    }
                break
    return {"status": "no_reservable_cpu_set", "allowed": sorted(allowed), "threads": threads}


def main():
    threads = int(sys.argv[1])
    path = Path(sys.argv[2])
    result = choose(threads)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if result["status"] != "ok":
        print(json.dumps(result), file=sys.stderr)
        return 2
    print("reserve=" + ",".join(map(str, result["reserved"])))
    print("pin=" + ",".join(map(str, result["pin"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
