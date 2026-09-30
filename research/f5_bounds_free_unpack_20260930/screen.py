#!/usr/bin/env python3
"""Frozen local same-binary F5 bounds-free unpack screen; exploratory only."""

import hashlib
import itertools
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BIN = Path(os.environ.get("F5_SCREEN_BINARY", ROOT / "target/release/examples/f4_f2_bench"))
SEEDS = ("0", "badc0de1")
PRIMARY = "f5_n24_m24_d4"
CASES = {
    "f5_n12_m12_d4", "f5_n16_m16_d3", "f5_n16_m16_d4",
    "f5_n20_m20_d3", "f5_n20_m20_d4", "f5_n24_m24_d3", PRIMARY,
}
SIGNATURE = (
    "row_space_fp", "rows_fp", "output_terms", "rank", "rows_f4",
    "rows_built", "cols", "rows_pruned", "criterion_word_ops", "reduce_word_ops",
)
BASE = {
    "KIC_F5_ECHELON": "2",
    "KIC_F5_FUSED_BUILD": "1",
    "KIC_F5_DIRECT_PACK": "1",
    "KIC_F5_UNPACK_DIRECT": "1",
    "KIC_GF2_BRANCHLESS_STRIP": "1",
    "KIC_GF2_REUSE_TABLE": "1",
    "KIC_GF2_TABLES": "4",
    "KIC_GF2_SIMD": "1",
    "KIC_GF2_DEFER_ABOVE": "0",
    "KIC_GF2_WORD_BATCH": "0",
    "KIC_F5_AVX512_UNPACK": "0",
    "RAYON_NUM_THREADS": "1",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, report):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    temp.replace(path)


def run(seed, mode, phase, pair, position):
    env = os.environ.copy()
    env.update(BASE)
    env["KIC_F5_UNPACK_UNCHECKED"] = str(mode)
    cmd = [str(BIN), "1", "24", "f5", seed]
    record = {
        "seed": seed, "mode": mode, "phase": phase, "pair": pair,
        "position": position, "cmd": cmd, "options": dict(BASE),
        "load_before": os.getloadavg(),
    }
    record["options"]["KIC_F5_UNPACK_UNCHECKED"] = str(mode)
    start = time.monotonic()
    try:
        proc = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=120)
        record.update(
            process_wall_ms=(time.monotonic() - start) * 1000,
            returncode=proc.returncode, stdout=proc.stdout, stderr=proc.stderr,
            load_after=os.getloadavg(),
        )
        if proc.returncode != 0:
            record["status"] = "failure"
            return record
        cases = [json.loads(line) for line in proc.stdout.splitlines() if line.strip()]
        record["cases"] = {case["case"]: case for case in cases}
        record["status"] = (
            "ok" if len(cases) == len(CASES) and set(record["cases"]) == CASES
            else "output_error"
        )
    except subprocess.TimeoutExpired as exc:
        record.update(
            status="timeout", process_wall_ms=(time.monotonic() - start) * 1000,
            stdout=exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else exc.stdout,
            stderr=exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else exc.stderr,
            load_after=os.getloadavg(),
        )
    except (json.JSONDecodeError, KeyError, TypeError) as exc:
        record.update(status="output_error", error=str(exc), load_after=os.getloadavg())
    return record


def fingerprint(record):
    return {
        name: {field: case[field] for field in SIGNATURE}
        for name, case in record["cases"].items()
    }


def bootstrap(values):
    ordered = sorted(
        statistics.median(values[i] for i in indices)
        for indices in itertools.product(range(len(values)), repeat=len(values))
    )
    return [ordered[int((len(ordered) - 1) * p)] for p in (0.025, 0.975)]


def pairs(records, seed, case, phase, field):
    values = []
    for pair in range(5):
        group = [r for r in records if r["seed"] == seed and r["phase"] == phase and r["pair"] == pair]
        if len(group) != 2:
            raise ValueError(f"incomplete {phase} pair {pair}")
        if phase == "aa":
            left, right = sorted(group, key=lambda r: r["position"])
        else:
            left = next(r for r in group if r["mode"] == 0)
            right = next(r for r in group if r["mode"] == 1)
        values.append(left["cases"][case][field] / right["cases"][case][field])
    return values


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
    expected = {}
    for seed in SEEDS:
        sequence = [("warmup", 0, 0, 0), ("warmup", 0, 1, 1)]
        for pair in range(5):
            sequence.extend([("aa", pair, 0, 0), ("aa", pair, 1, 0)])
        for pair in range(5):
            modes = [0, 1] if pair % 2 == 0 else [1, 0]
            sequence.extend(("paired", pair, position, mode) for position, mode in enumerate(modes))
        for phase, pair, position, mode in sequence:
            record = run(seed, mode, phase, pair, position)
            if record["status"] == "ok":
                try:
                    actual = fingerprint(record)
                    if seed not in expected:
                        expected[seed] = actual
                    elif actual != expected[seed]:
                        record["status"] = "output_mismatch"
                    for name, case in record["cases"].items():
                        selected = bool(mode == 1 and name == PRIMARY)
                        if case["unchecked_unpack_used"] != selected:
                            record["status"] = "route_mismatch"
                except (KeyError, TypeError) as exc:
                    record.update(status="output_error", error=str(exc))
            report["runs"].append(record)
            save(output, report)
            if record["status"] != "ok":
                report["status"] = record["status"]
                save(output, report)
                print(output)
                return 1
    report["summary"] = {}
    for seed in SEEDS:
        report["summary"][seed] = {}
        for case in sorted(CASES):
            report["summary"][seed][case] = {}
            for field in ("wall_ms", "reduce_ms", "unpack_ms"):
                aa = pairs(report["runs"], seed, case, "aa", field)
                ab = pairs(report["runs"], seed, case, "paired", field)
                report["summary"][seed][case][field] = {
                    "aa_ratios": aa, "aa_min": min(aa), "aa_max": max(aa),
                    "paired_ratios": ab, "paired_median": statistics.median(ab),
                    "bootstrap_95pct": bootstrap(ab),
                }
    report["status"] = "complete"
    save(output, report)
    print(output)
    for seed in SEEDS:
        print(seed, json.dumps(report["summary"][seed][PRIMARY], sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
