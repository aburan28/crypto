#!/usr/bin/env python3
"""Frozen local F5 stratified-pivot screen; not an eligible x86 gain claim."""

import hashlib
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
BIN = ROOT / "target/release/examples/f4_f2_bench"
SEEDS = ("0", "badc0de1")
CASES = {
    "f5_n12_m12_d4", "f5_n16_m16_d3", "f5_n16_m16_d4",
    "f5_n20_m20_d3", "f5_n20_m20_d4", "f5_n24_m24_d3",
    "f5_n24_m24_d4",
}
PRIMARY = "f5_n24_m24_d4"
SIGNATURE = (
    "row_space_fp", "rank", "criterion_word_ops", "rows_f4",
    "rows_built", "cols", "rows_pruned",
)
BASE = {
    "KIC_F5_ECHELON": "2", "KIC_F5_FUSED_BUILD": "1",
    "KIC_F5_DIRECT_PACK": "1", "KIC_F5_UNPACK_DIRECT": "1",
    "KIC_GF2_FORCE_AVX2": "1", "KIC_GF2_REUSE_TABLE": "1",
    "KIC_GF2_SIMD": "1", "KIC_GF2_DEFER_ABOVE": "0",
    "KIC_GF2_WORD_BATCH": "0", "KIC_F5_AVX512_UNPACK": "0",
    "KIC_GF2_AVX2_TABLE_BUILD": "0", "RAYON_NUM_THREADS": "1",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, report):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def run(seed, mode, phase, pair, position):
    env = os.environ.copy()
    env.update(BASE)
    env["KIC_GF2_STRATIFIED_WEIGHT_PIVOT"] = str(mode)
    command = [str(BIN), "1", "24", "f5", seed]
    record = {
        "seed": seed, "mode": mode, "phase": phase,
        "pair": pair, "position": position, "command": command,
        "load_before": os.getloadavg(),
    }
    start = time.monotonic()
    try:
        proc = subprocess.run(command, env=env, capture_output=True,
                              text=True, timeout=120)
        record.update(
            process_wall_ms=(time.monotonic() - start) * 1000,
            returncode=proc.returncode, stdout=proc.stdout,
            stderr=proc.stderr, load_after=os.getloadavg(),
        )
        if proc.returncode != 0:
            record["status"] = "failure"
            return record
        rows = [json.loads(line) for line in proc.stdout.splitlines()
                if line.strip()]
        record["cases"] = {row["case"]: row for row in rows}
        record["status"] = "ok" if len(rows) == len(CASES) and set(record["cases"]) == CASES else "output_error"
    except subprocess.TimeoutExpired as exc:
        record.update(
            status="timeout", process_wall_ms=(time.monotonic() - start) * 1000,
            stdout=exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else exc.stdout,
            stderr=exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else exc.stderr,
            load_after=os.getloadavg(),
        )
    except (json.JSONDecodeError, KeyError, TypeError) as exc:
        record.update(status="output_error", error=str(exc),
                      load_after=os.getloadavg())
    return record


def signature(record):
    return {name: {field: case[field] for field in SIGNATURE}
            for name, case in record["cases"].items()}


def pair_ratios(runs, seed, phase, case, field):
    ratios = []
    for pair in range(5):
        rows = [r for r in runs if r["seed"] == seed and
                r["phase"] == phase and r["pair"] == pair]
        assert len(rows) == 2
        if phase == "aa":
            rows.sort(key=lambda r: r["position"])
        else:
            rows.sort(key=lambda r: r["mode"])
        ratios.append(rows[0]["cases"][case][field] /
                      rows[1]["cases"][case][field])
    return ratios


def main():
    output = HERE / ("screen_" + time.strftime("%Y-%m-%dT%H%M%SZ", time.gmtime()) + ".json")
    report = {
        "status": "running", "protocol": "PROTOCOL.md", "seeds": SEEDS,
        "host": {"platform": platform.platform(), "machine": platform.machine(),
                 "processor": platform.processor(), "cpu_count": os.cpu_count(),
                 "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
                 "load_start": os.getloadavg()},
        "git_sha": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                            cwd=ROOT, text=True).strip(),
        "binary_sha256": sha(BIN),
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in (
            ROOT / "src/cryptanalysis/matrix_f5_f2.rs",
            ROOT / "src/cryptanalysis/gf2_elim.rs",
            ROOT / "examples/f4_f2_bench.rs")},
        "runs": [],
    }
    save(output, report)
    expected = {}
    for seed in SEEDS:
        sequence = [("warmup", 0, mode, mode) for mode in range(2)]
        for pair in range(5):
            sequence.extend(("aa", pair, pos, 0) for pos in range(2))
        for pair in range(5):
            modes = (0, 1) if pair % 2 == 0 else (1, 0)
            sequence.extend(("paired", pair, pos, mode)
                            for pos, mode in enumerate(modes))
        for phase, pair, position, mode in sequence:
            record = run(seed, mode, phase, pair, position)
            if record["status"] == "ok":
                try:
                    sig = signature(record)
                    if seed not in expected:
                        expected[seed] = sig
                    elif sig != expected[seed]:
                        record["status"] = "output_mismatch"
                        record["expected"] = expected[seed]
                        record["actual"] = sig
                    if record["cases"][PRIMARY]["stratified_pivot_used"] != bool(mode):
                        record["status"] = "route_mismatch"
                except (KeyError, TypeError) as exc:
                    record["status"] = "output_error"
                    record["error"] = str(exc)
            report["runs"].append(record)
            save(output, report)
            if record["status"] != "ok":
                report["status"] = record["status"]
                save(output, report)
                print(output)
                return 1
    report["signatures"] = expected
    report["summary"] = {}
    for seed in SEEDS:
        report["summary"][seed] = {}
        for case in sorted(CASES):
            aa = pair_ratios(report["runs"], seed, "aa", case, "wall_ms")
            ab = pair_ratios(report["runs"], seed, "paired", case, "wall_ms")
            report["summary"][seed][case] = {
                "aa_ratios": aa, "aa_min": min(aa), "aa_max": max(aa),
                "prior_new_ratios": ab, "prior_new_median": statistics.median(ab),
            }
    report["host"]["load_end"] = os.getloadavg()
    report["status"] = "complete"
    save(output, report)
    print(output)
    for seed in SEEDS:
        print(seed, PRIMARY, report["summary"][seed][PRIMARY])
    return 0


if __name__ == "__main__":
    sys.exit(main())
