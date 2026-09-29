#!/usr/bin/env python3
"""Frozen same-binary matrix-F5 default/prior/AVX-512-unpack comparison.

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
WORKLOADS = {
    "frozen": 0,
    "holdout_a": 0xBADC0DE1,
    "holdout_b": 0x5EED2026,
    "holdout_c": 0xF5C02A28,
}
EXPECTED_CASES = {
    "f5_n12_m12_d4",
    "f5_n16_m16_d3",
    "f5_n16_m16_d4",
    "f5_n20_m20_d3",
    "f5_n20_m20_d4",
    "f5_n24_m24_d3",
    PRIMARY,
}
SIGNATURE_FIELDS = ("row_space_fp", "rank", "rows_pruned", "criterion_word_ops")


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


def run_one(binary, workload, seed_xor, threads, mode, phase, pair, position):
    env = os.environ.copy()
    env["RAYON_NUM_THREADS"] = str(threads)
    env["KIC_GF2_DEFER_ABOVE"] = "0"
    env["KIC_GF2_WORD_BATCH"] = "0"
    env["KIC_F5_ECHELON"] = "2" if mode else "0"
    env["KIC_F5_FUSED_BUILD"] = "1" if mode else "0"
    env["KIC_GF2_FORCE_AVX2"] = "1" if mode else "0"
    env["KIC_F5_AVX512_UNPACK"] = "1" if mode == 2 else "0"
    env["KIC_GF2_SIMD"] = "1"
    record = {
        "phase": phase,
        "workload": workload,
        "seed_xor_hex": f"{seed_xor:016x}",
        "pair": pair,
        "position": position,
        "mode": mode,
        "rayon_threads": threads,
        "load_before": os.getloadavg(),
    }
    t0 = time.monotonic()
    try:
        proc = subprocess.run(
            [str(binary), "1", "24", "f5", f"{seed_xor:016x}"],
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


def pair_ratios(records, workload, case, phase, field, left=0, right=0):
    pairs = sorted({r["pair"] for r in records if r["workload"] == workload and r["phase"] == phase})
    ratios = []
    for pair in pairs:
        group = [r for r in records if r["workload"] == workload and r["phase"] == phase and r["pair"] == pair]
        if len(group) != (2 if phase == "aa" else 3):
            raise ValueError(f"incomplete {phase} pair {pair}")
        if phase == "aa":
            reference, candidate = sorted(group, key=lambda r: r["position"])
        else:
            reference = next(r for r in group if r["mode"] == left)
            candidate = next(r for r in group if r["mode"] == right)
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
            for field in ("wall_ms", "f5_build_ms", "reduce_ms", "unpack_ms"):
                aa = pair_ratios(records, workload, case, "aa", field)
                comparisons = {
                    name: pair_ratios(records, workload, case, "triad", field, left, right)
                    for name, left, right in (
                        ("default_prior", 0, 1),
                        ("default_new", 0, 2),
                        ("prior_new", 1, 2),
                    )
                }
                summary[workload][case][field] = {
                    "aa_ratios": aa,
                    "aa_min": min(aa),
                    "aa_max": max(aa),
                    **{
                        name: {
                            "ratios": ratios,
                            "median": statistics.median(ratios),
                            "bootstrap_95pct": median_bootstrap_interval(ratios),
                        }
                        for name, ratios in comparisons.items()
                    },
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
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--pairs", type=int, default=5)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if args.pairs < 1:
        parser.error("--pairs must be positive")
    if args.threads < 1:
        parser.error("--threads must be positive")
    binary = args.binary.resolve()
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
        "binary_sha256": sha256(binary),
        "benchmark_source_sha256": sha256(Path("examples/f4_f2_bench.rs")),
        "f5_source_sha256": sha256(Path("src/cryptanalysis/matrix_f5_f2.rs")),
        "gf2_source_sha256": sha256(Path("src/cryptanalysis/gf2_elim.rs")),
        "mono_source_sha256": sha256(Path("src/cryptanalysis/pq_groebner_f2.rs")),
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
    if not {"avx2", "avx512f", "bmi2"}.issubset(receipt["host"].get("cpu_features", [])):
        receipt["status"] = "unsupported_host"
        write_receipt(output, receipt)
        print("AVX-512F, AVX2 and BMI2 required for this experiment", file=sys.stderr)
        return 2
    write_receipt(output, receipt)
    expected = {}
    raw_expected = {}
    sequence = []
    for workload, seed_xor in WORKLOADS.items():
        sequence.extend((workload, seed_xor, "warmup", 0, mode, mode) for mode in range(3))
        for pair in range(args.pairs):
            sequence.extend(((workload, seed_xor, "aa", pair, 0, 0), (workload, seed_xor, "aa", pair, 1, 0)))
        for pair in range(args.pairs):
            arms = [(pair + position) % 3 for position in range(3)]
            sequence.extend((workload, seed_xor, "triad", pair, position, mode) for position, mode in enumerate(arms))

    for workload, seed_xor, phase, pair, position, mode in sequence:
        record = run_one(binary, workload, seed_xor, args.threads, mode, phase, pair, position)
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
                if record["status"] == "ok":
                    raw = {name: case["rows_fp"] for name, case in record["cases"].items()}
                    raw_key = f"{workload}:{mode}"
                    if raw_key not in raw_expected:
                        raw_expected[raw_key] = raw
                        receipt["raw_output_signatures"] = raw_expected
                    elif raw != raw_expected[raw_key]:
                        record["status"] = "raw_output_mismatch"
                        record["expected_raw_signature"] = raw_expected[raw_key]
                        record["actual_raw_signature"] = raw
                    if record["status"] == "ok" and mode == 2:
                        prior = raw_expected.get(f"{workload}:1")
                        if prior != raw:
                            record["status"] = "prior_new_raw_mismatch"
                    if record["status"] == "ok" and mode in (1, 2):
                        counts = {name: case["output_terms"] for name, case in record["cases"].items()}
                        key = f"{workload}:terms"
                        if key not in raw_expected:
                            raw_expected[key] = counts
                        elif raw_expected[key] != counts:
                            record["status"] = "prior_new_terms_mismatch"
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
