#!/usr/bin/env python3
"""Nonpromoting, same-binary local matrix-F5 screen; retain every process."""

import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BIN = ROOT / "target/release/examples/f4_f2_bench"
SEEDS = ("0", "badc0de1")
BASE = {
    "KIC_F5_ECHELON": "2",
    "KIC_F5_FUSED_BUILD": "1",
    "KIC_F5_DIRECT_PACK": "1",
    "KIC_GF2_FORCE_AVX2": "1",
    "KIC_GF2_REUSE_TABLE": "1",
    "KIC_GF2_SIMD": "1",
    "KIC_GF2_DEFER_ABOVE": "0",
    "KIC_GF2_WORD_BATCH": "0",
    "KIC_F5_AVX512_UNPACK": "0",
    "RAYON_NUM_THREADS": "1",
}
ARMS = {
    "tables": (("t4", {"KIC_GF2_TABLES": "4"}),
               ("t6", {"KIC_GF2_TABLES": "6"}),
               ("t8", {"KIC_GF2_TABLES": "8"})),
    "pivot": (("p1", {"KIC_GF2_PIVOT_CANDIDATES": "1"}),
              ("p4", {"KIC_GF2_PIVOT_CANDIDATES": "4"})),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(seed, name, options, phase, pair):
    env = os.environ.copy()
    env.update(BASE)
    env.update(options)
    cmd = [str(BIN), "1", "24", "f5", seed]
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=300)
    cases = {}
    for line in proc.stdout.splitlines():
        try:
            case = json.loads(line)
            cases[case["case"]] = case
        except (KeyError, json.JSONDecodeError):
            pass
    return {"seed": seed, "arm": name, "options": options, "phase": phase,
            "pair": pair, "cmd": cmd, "returncode": proc.returncode,
            "stdout": proc.stdout, "stderr": proc.stderr, "cases": cases}


def main():
    group = sys.argv[1]
    if group not in ARMS:
        raise SystemExit("usage: screen.py tables|pivot")
    arms = ARMS[group]
    report = {
        "group": group, "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": {"platform": platform.platform(), "machine": platform.machine(),
                 "processor": platform.processor(), "cpu_count": os.cpu_count()},
        "git_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "binary_sha256": sha(BIN),
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in (
            ROOT / "src/cryptanalysis/gf2_elim.rs",
            ROOT / "src/cryptanalysis/matrix_f5_f2.rs",
            ROOT / "examples/f4_f2_bench.rs")},
        "runs": [],
    }
    for seed in SEEDS:
        for name, options in arms:
            report["runs"].append(run(seed, name, options, "warmup", 0))
        for pair in range(3):
            for name, options in (arms if pair % 2 == 0 else tuple(reversed(arms))):
                report["runs"].append(run(seed, name, options, "paired", pair))
    output = HERE / f"screen_{group}_{report['timestamp_utc'].replace(':', '')}.json"
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(output)
    for seed in SEEDS:
        for name, _ in arms:
            calls = [r for r in report["runs"] if r["seed"] == seed and
                     r["arm"] == name and r["phase"] == "paired"]
            data = [r["cases"].get("f5_n24_m24_d4", {}) for r in calls]
            print(seed, name, [(x.get("wall_ms"), x.get("reduce_ms"),
                                 x.get("unpack_ms"), x.get("output_terms"),
                                 x.get("row_space_fp")) for x in data])


if __name__ == "__main__":
    main()
