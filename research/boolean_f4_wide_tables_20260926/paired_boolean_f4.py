#!/usr/bin/env python3
"""Frozen same-binary Boolean F4 table-width A/A and A/B stage comparison.

The example reports timings inside the F4 call; process launch and input
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


PRIMARY = "n20_m30"
WORKLOADS = {
    "frozen": 0,
    "holdout_a": 0x1AC0FFEE,
    "holdout_b": 0x2468ACE0,
}
EXPECTED_CASES = {
    "n12_m12",
    "n14_m14",
    "n16_m16",
    "n16_m24",
    "n18_m18",
    "n20_m20",
    PRIMARY,
}
SIGNATURE_FIELDS = (
    "basis_fp", "basis_len", "steps", "divisor_tests", "word_xors",
    "matrix_rows_max", "matrix_cols_max",
)


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


def run_one(binaries, workload, seed_xor, threads, mode, phase, pair, position):
    env = os.environ.copy()
    env["RAYON_NUM_THREADS"] = str(threads)
    env["F4_F2_BITMAP_SEEN"] = "1"
    record = {
        "phase": phase,
        "workload": workload,
        "seed_xor_hex": f"{seed_xor:016x}",
        "rayon_threads": threads,
        "pair": pair,
        "position": position,
        "mode": mode,
        "load_before": os.getloadavg(),
    }
    binary = binaries[mode]
    t0 = time.monotonic()
    try:
        proc = subprocess.run(
            [str(binary), "1", "20", "f4", f"{seed_xor:016x}"],
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


def pair_ratios(records, workload, case, phase, field):
    pairs = sorted({r["pair"] for r in records if r["workload"] == workload and r["phase"] == phase})
    ratios = []
    for pair in pairs:
        group = [r for r in records if r["workload"] == workload and r["phase"] == phase and r["pair"] == pair]
        if len(group) != 2:
            raise ValueError(f"incomplete {phase} pair {pair}")
        if phase == "aa":
            reference, candidate = sorted(group, key=lambda r: r["position"])
        else:
            reference = next(r for r in group if r["mode"] == 0)
            candidate = next(r for r in group if r["mode"] == 1)
        a = reference["cases"][case][field]
        b = candidate["cases"][case][field]
        ratios.append(a / b)
    return ratios


def summarize(records):
    summary = {}
    for workload in WORKLOADS:
        summary[workload] = {}
        for case in sorted(EXPECTED_CASES):
            summary[workload][case] = {}
            for field in ("build_ms", "eliminate_ms", "wall_ms"):
                aa = pair_ratios(records, workload, case, "aa", field)
                ab = pair_ratios(records, workload, case, "ab", field)
                summary[workload][case][field] = {
                    "aa_ratios": aa,
                    "aa_min": min(aa),
                    "aa_max": max(aa),
                    "ab_ratios": ab,
                    "ab_median": statistics.median(ab),
                    "ab_min": min(ab),
                    "ab_bootstrap_95pct": median_bootstrap_interval(ab),
                }
    return summary


def linux_host_details():
    details = {}
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists():
        fields = {
            key.strip(): value.strip()
            for key, value in (line.split(":", 1) for line in cpuinfo.read_text().splitlines() if ":" in line)
        }
        details["cpu_model"] = fields.get("model name", "").strip()
        flags = set(fields.get("flags", "").split())
        details["cpu_features"] = sorted(flags & {"popcnt", "avx2", "avx512f", "pclmulqdq", "bmi2"})
    meminfo = Path("/proc/meminfo")
    if meminfo.exists():
        details["mem_total"] = next((line for line in meminfo.read_text().splitlines() if line.startswith("MemTotal:")), "")
    return details


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference_binary", type=Path)
    parser.add_argument("candidate_binary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--pairs", type=int, default=5)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--candidate-feature", choices=("f4-wide-tables", "f4-five-tables", "default-five"), default="f4-wide-tables")
    args = parser.parse_args()
    if args.pairs < 1:
        parser.error("--pairs must be positive")
    if args.threads < 1:
        parser.error("--threads must be positive")
    binaries = (args.reference_binary.resolve(), args.candidate_binary.resolve())
    binary_hashes = [sha256(binary) for binary in binaries]
    if binary_hashes[0] == binary_hashes[1]:
        parser.error("reference and candidate binaries have identical hashes")
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    affinity = None
    if hasattr(os, "sched_getaffinity"):
        allowed = sorted(os.sched_getaffinity(0))
        if len(allowed) < args.threads:
            parser.error(f"{args.threads} threads requested but only {len(allowed)} CPUs allowed")
        affinity = allowed[:args.threads]
        os.sched_setaffinity(0, set(affinity))
    receipt = {
        "status": "running",
        "primary_case": PRIMARY,
        "workloads": {name: f"{seed:016x}" for name, seed in WORKLOADS.items()},
        "pairs_per_phase": args.pairs,
        "rayon_threads": args.threads,
        "binary_sha256": binary_hashes,
        "candidate_feature": args.candidate_feature,
        "benchmark_source_sha256": sha256(Path("examples/f4_f2_bench.rs")),
        "kernel_source_sha256": sha256(Path("src/cryptanalysis/pq_f4_f2.rs")),
        "git_sha": os.environ.get("GITHUB_SHA"),
        "host": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cpu_count": os.cpu_count(),
            "pinned_cpus": affinity,
            "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
            "load_start": os.getloadavg(),
            **linux_host_details(),
        },
        "runs": [],
    }
    write_receipt(output, receipt)
    expected = {}
    sequence = []
    for workload, seed_xor in WORKLOADS.items():
        sequence.extend(((workload, seed_xor, "warmup", 0, 0, 0), (workload, seed_xor, "warmup", 0, 1, 1)))
        for pair in range(args.pairs):
            sequence.extend(((workload, seed_xor, "aa", pair, 0, 0), (workload, seed_xor, "aa", pair, 1, 0)))
        for pair in range(args.pairs):
            arms = (0, 1) if pair % 2 == 0 else (1, 0)
            sequence.extend((workload, seed_xor, "ab", pair, position, mode) for position, mode in enumerate(arms))

    for workload, seed_xor, phase, pair, position, mode in sequence:
        record = run_one(binaries, workload, seed_xor, args.threads, mode, phase, pair, position)
        if record["status"] == "ok":
            try:
                actual = signature(record)
            except (KeyError, TypeError) as exc:
                record["status"] = "output_error"
                record["error"] = str(exc)
            else:
                if workload not in expected:
                    expected[workload] = actual
                    receipt["output_signatures"] = expected
                elif actual != expected[workload]:
                    record["status"] = "output_mismatch"
                    record["expected_signature"] = expected[workload]
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
    print(json.dumps({name: cases[PRIMARY] for name, cases in receipt["summary"].items()}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
