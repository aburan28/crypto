#!/usr/bin/env python3
"""Frozen local flat-row F5 screen; not an eligible x86 claim."""

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
ARMS = (
    ("prior", "0"),
    ("flat", "1"),
)
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


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(seed, arm, phase, pair, position):
    name, flat = arm
    env = os.environ.copy()
    env.update(BASE)
    env["KIC_F5_UNPACK_DIRECT"] = "1"
    env["KIC_F5_FLAT_ROWS"] = flat
    cmd = [str(BIN), "1", "24", "f5", seed]
    start = time.monotonic()
    try:
        proc = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=120)
        cases = [json.loads(line) for line in proc.stdout.splitlines() if line.strip()]
        status = "ok" if proc.returncode == 0 and len(cases) == 7 else "failure"
        return {"seed": seed, "arm": name, "phase": phase, "pair": pair,
                "position": position, "options": {"KIC_F5_UNPACK_DIRECT": "1",
                "KIC_F5_FLAT_ROWS": flat}, "cmd": cmd,
                "status": status, "returncode": proc.returncode,
                "process_wall_ms": (time.monotonic() - start) * 1000,
                "stdout": proc.stdout, "stderr": proc.stderr,
                "cases": {case["case"]: case for case in cases}}
    except subprocess.TimeoutExpired as exc:
        return {"seed": seed, "arm": name, "phase": phase, "pair": pair,
                "position": position, "cmd": cmd, "status": "timeout",
                "process_wall_ms": (time.monotonic() - start) * 1000,
                "stdout": exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else exc.stdout,
                "stderr": exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else exc.stderr}


def save(path, report):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def main():
    output = HERE / ("screen_" + time.strftime("%Y-%m-%dT%H%M%SZ", time.gmtime()) + ".json")
    report = {
        "status": "running", "protocol": "PROTOCOL.md", "seeds": SEEDS,
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
    save(output, report)
    for seed in SEEDS:
        sequence = [("warmup", 0, i, arm) for i, arm in enumerate(ARMS)]
        for pair in range(5):
            sequence.extend(("paired", pair, position, ARMS[(pair + position) % 2])
                            for position in range(2))
        for phase, pair, position, arm in sequence:
            record = run(seed, arm, phase, pair, position)
            report["runs"].append(record)
            save(output, report)
            if record["status"] != "ok":
                report["status"] = record["status"]
                save(output, report)
                print(output)
                return 1
    report["status"] = "complete"
    save(output, report)
    print(output)
    for seed in SEEDS:
        for name, _ in ARMS:
            primary = [r["cases"]["f5_n24_m24_d4"] for r in report["runs"]
                       if r["seed"] == seed and r["phase"] == "paired" and r["arm"] == name]
            print(seed, name, [(p["wall_ms"], p["reduce_ms"], p["unpack_ms"],
                                p["output_terms"], p["row_space_fp"]) for p in primary])
    return 0


if __name__ == "__main__":
    sys.exit(main())
