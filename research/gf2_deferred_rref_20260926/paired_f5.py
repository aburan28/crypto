#!/usr/bin/env python3
"""Frozen same-binary matrix-F5 A/A and A/B stage comparison.

The example reports timings inside the F5 call; process launch and input
construction are outside that interval. All calls, including failures, are
retained in the output file before the next call starts.
"""

import argparse
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


PRIMARY = "f5_n24_m24_d4"
EXPECTED_CASES = {
    "f5_n12_m12_d4",
    "f5_n16_m16_d3",
    "f5_n16_m16_d4",
    "f5_n20_m20_d3",
    "f5_n20_m20_d4",
    "f5_n24_m24_d3",
    PRIMARY,
}
SIGNATURE_FIELDS = ("rows_fp", "rank", "rows_pruned", "criterion_word_ops")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def percentiles(values):
    ordered = sorted(values)
    return [ordered[int((len(ordered) - 1) * p)] for p in (0.025, 0.975)]


def median_bootstrap_interval(values):
    # Five pairs give exactly 5^5 = 3125 bootstrap resamples.
    samples = (
        statistics.median(values[i] for i in indices)
        for indices in itertools.product(range(len(values)), repeat=len(values))
    )
    return percentiles(samples)


def write_receipt(path, receipt):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def run_one(binary, mode, phase, pair, position):
    env = os.environ.copy()
    env["RAYON_NUM_THREADS"] = "1"
    env["KIC_GF2_DEFER_ABOVE"] = str(mode)
    record = {
        "phase": phase,
        "pair": pair,
        "position": position,
        "mode": mode,
        "load_before": os.getloadavg(),
    }
    t0 = time.monotonic()
    try:
        proc = subprocess.run(
            [str(binary), "1", "24", "f5"],
            capture_output=True,
            text=True,
            env=env,
            timeout=120,
            check=False,
        )
        record.update(
            process_wall_ms=(time.monotonic() - t0) * 1000,
            exit_code=proc.returncode,
            stdout=proc.stdout,
            stderr=proc.stderr,
            load_after=os.getloadavg(),
        )
        if proc.returncode != 0:
            record["status"] = "failure"
            return record
        cases = [json.loads(line) for line in proc.stdout.splitlines() if line.strip()]
        record["cases"] = {case["case"]: case for case in cases}
        record["status"] = (
            "ok"
            if len(record["cases"]) == len(cases)
            and set(record["cases"]) == EXPECTED_CASES
            else "output_error"
        )
    except subprocess.TimeoutExpired as exc:
        record.update(
            status="timeout",
            process_wall_ms=(time.monotonic() - t0) * 1000,
            stdout=exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else exc.stdout,
            stderr=exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else exc.stderr,
            load_after=os.getloadavg(),
        )
    except (json.JSONDecodeError, KeyError, TypeError) as exc:
        record.update(status="output_error", error=str(exc), load_after=os.getloadavg())
    return record


def signature(record):
    return {
        name: {field: case[field] for field in SIGNATURE_FIELDS}
        for name, case in record["cases"].items()
    }


def pair_ratios(records, phase, field):
    pairs = sorted({r["pair"] for r in records if r["phase"] == phase})
    ratios = []
    for pair in pairs:
        group = [r for r in records if r["phase"] == phase and r["pair"] == pair]
        if len(group) != 2:
            raise ValueError(f"incomplete {phase} pair {pair}")
        if phase == "aa":
            reference, candidate = sorted(group, key=lambda r: r["position"])
        else:
            reference = next(r for r in group if r["mode"] == 0)
            candidate = next(r for r in group if r["mode"] == 1)
        a = reference["cases"][PRIMARY][field]
        b = candidate["cases"][PRIMARY][field]
        ratios.append(a / b)
    return ratios


def summarize(records):
    summary = {}
    for field in ("reduce_ms", "wall_ms"):
        aa = pair_ratios(records, "aa", field)
        ab = pair_ratios(records, "ab", field)
        summary[field] = {
            "aa_ratios": aa,
            "aa_min": min(aa),
            "aa_max": max(aa),
            "ab_ratios": ab,
            "ab_median": statistics.median(ab),
            "ab_min": min(ab),
            "ab_bootstrap_95pct": median_bootstrap_interval(ab),
        }
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--pairs", type=int, default=5)
    args = parser.parse_args()
    if args.pairs < 1:
        parser.error("--pairs must be positive")
    binary = args.binary.resolve()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    affinity = None
    if hasattr(os, "sched_getaffinity"):
        allowed = os.sched_getaffinity(0)
        affinity = min(allowed)
        os.sched_setaffinity(0, {affinity})
    receipt = {
        "status": "running",
        "primary_case": PRIMARY,
        "pairs_per_phase": args.pairs,
        "binary_sha256": sha256(binary),
        "benchmark_source_sha256": sha256(Path("examples/f4_f2_bench.rs")),
        "kernel_source_sha256": sha256(Path("src/cryptanalysis/gf2_elim.rs")),
        "git_sha": os.environ.get("GITHUB_SHA"),
        "host": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cpu_count": os.cpu_count(),
            "pinned_cpu": affinity,
            "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
            "load_start": os.getloadavg(),
        },
        "runs": [],
    }
    write_receipt(output, receipt)
    expected = None
    sequence = [("warmup", 0, 0, 0), ("warmup", 0, 1, 1)]
    for pair in range(args.pairs):
        sequence.extend((("aa", pair, 0, 0), ("aa", pair, 1, 0)))
    for pair in range(args.pairs):
        arms = (0, 1) if pair % 2 == 0 else (1, 0)
        sequence.extend(("ab", pair, position, mode) for position, mode in enumerate(arms))

    for phase, pair, position, mode in sequence:
        record = run_one(binary, mode, phase, pair, position)
        if record["status"] == "ok":
            try:
                actual = signature(record)
            except (KeyError, TypeError) as exc:
                record["status"] = "output_error"
                record["error"] = str(exc)
            else:
                if expected is None:
                    expected = actual
                    receipt["output_signature"] = expected
                elif actual != expected:
                    record["status"] = "output_mismatch"
                    record["expected_signature"] = expected
                    record["actual_signature"] = actual
        receipt["runs"].append(record)
        write_receipt(output, receipt)
        if record["status"] != "ok":
            receipt["status"] = record["status"]
            write_receipt(output, receipt)
            return 1

    receipt["summary"] = summarize(receipt["runs"])
    receipt["host"]["load_end"] = os.getloadavg()
    receipt["status"] = "complete"
    write_receipt(output, receipt)
    print(json.dumps(receipt["summary"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
