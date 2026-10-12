#!/usr/bin/env python3
"""One-seed same-binary local row-basis F5 comparison.

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
SIGNATURE_FIELDS = (
    "row_space_fp", "rows_fp", "output_terms", "rank", "rows_f4",
    "rows_built", "cols", "rows_pruned", "criterion_word_ops", "reduce_word_ops",
    "direct_pack_used", "direct_unpack_used",
)


def compare_structural(reference, candidate, slack):
    """Check selective echelon against the exact row-basis certificate."""
    if reference["status"] != "ok" or candidate["status"] != "ok":
        return ["incomplete process"]
    errors = []
    invariants = (
        "row_space_fp", "rank", "criterion_word_ops", "rows_f4",
        "rows_built", "cols", "rows_pruned", "direct_pack_used",
        "direct_unpack_used",
    )
    for name in sorted(EXPECTED_CASES):
        a, b = reference["cases"][name], candidate["cases"][name]
        for field in invariants:
            if a[field] != b[field]:
                errors.append(f"{name}: {field} mismatch")
        if a["column_cert_attempted"] or a["selected_original_used"]:
            errors.append(f"{name}: reference certificate route active")
        if name == PRIMARY:
            if not b["column_cert_attempted"] or not b["selected_original_used"]:
                errors.append(f"{name}: certificate route inactive or failed")
            if b["selected_cols"] != min(a["cols"], a["rows_built"] + slack):
                errors.append(f"{name}: selected width mismatch")
            if b["reduce_word_ops"] > a["reduce_word_ops"]:
                errors.append(f"{name}: counted work exceeded reference")
        else:
            if b["column_cert_attempted"] or b["selected_original_used"]:
                errors.append(f"{name}: inactive certificate route active")
            for field in ("rows_fp", "output_terms", "reduce_word_ops"):
                if a[field] != b[field]:
                    errors.append(f"{name}: inactive {field} mismatch")
    return errors


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
    env["KIC_F5_ECHELON"] = "4" if mode else "2"
    env["KIC_F5_COLUMN_SLACK"] = "512"
    env["KIC_F5_FUSED_BUILD"] = "1"
    env["KIC_GF2_FORCE_AVX2"] = "1"
    env["KIC_F5_AVX512_UNPACK"] = "0"
    env["KIC_F5_DIRECT_PACK"] = "1"
    env["KIC_GF2_REUSE_TABLE"] = "1"
    env["KIC_GF2_TABLES"] = "4"
    env["KIC_GF2_SIMD"] = "1"
    env["KIC_F5_UNPACK_DIRECT"] = "1"
    env["KIC_GF2_AVX2_TABLE_BUILD"] = "0"
    env["KIC_GF2_BRANCHLESS_STRIP"] = "1"
    env["KIC_GF2_ROW_BASIS_RANK"] = "1" if mode else "0"
    rank_tables = "8" if mode else "4"
    env["KIC_GF2_ROW_BASIS_TABLES"] = rank_tables
    record = {
        "phase": phase,
        "workload": workload,
        "seed_xor_hex": f"{seed_xor:016x}",
        "pair": pair,
        "position": position,
        "mode": mode,
        "column_certificate_mode": str(mode),
        "row_basis_rank_tables": rank_tables,
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
        if len(group) != 2:
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


def summarize(records, active_workloads):
    summary = {}
    for workload in active_workloads:
        summary[workload] = {}
        for case in sorted(EXPECTED_CASES):
            summary[workload][case] = {}
            for field in ("wall_ms", "f5_build_ms", "reduce_ms", "unpack_ms"):
                aa = pair_ratios(records, workload, case, "aa", field)
                ratios = pair_ratios(records, workload, case, "paired", field, 0, 1)
                summary[workload][case][field] = {
                    "aa_ratios": aa,
                    "aa_min": min(aa),
                    "aa_max": max(aa),
                    "prior_new": {
                        "ratios": ratios,
                        "median": statistics.median(ratios),
                        "bootstrap_95pct": median_bootstrap_interval(ratios),
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
    parser.add_argument("--workload", required=True, choices=WORKLOADS)
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
    active_workloads = {args.workload: WORKLOADS[args.workload]}
    receipt = {
        "status": "running",
        "primary_case": PRIMARY,
        "workloads": {name: f"{seed:016x}" for name, seed in active_workloads.items()},
        "pairs_per_phase": args.pairs,
        "rayon_threads": args.threads,
        "binary_sha256": sha256(binary),
        "benchmark_source_sha256": sha256(Path("examples/f4_f2_bench.rs")),
        "f5_source_sha256": sha256(Path("src/cryptanalysis/matrix_f5_f2.rs")),
        "row_builder_source_sha256": sha256(Path("src/cryptanalysis/koblitz_groebner.rs")),
        "gf2_source_sha256": sha256(Path("src/cryptanalysis/gf2_elim.rs")),
        "mono_source_sha256": sha256(Path("src/cryptanalysis/pq_groebner_f2.rs")),
        "git_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "head_ref": os.environ.get("GITHUB_HEAD_REF"),
        "host": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cpu_count": os.cpu_count(),
            "pinned_cpus": affinity,
            "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
            "load_start": os.getloadavg(),
            **linux_host_details(),
        },
        "qualified_isolation": False,
        "isolation_reason": "Apple ARM64 shared host: no Linux affinity/PSI reservation",
        "runs": [],
    }
    write_receipt(output, receipt)
    expected = {}
    reference_record = {}
    sequence = []
    for workload, seed_xor in active_workloads.items():
        sequence.extend((workload, seed_xor, "warmup", 0, mode, mode) for mode in range(2))
        for pair in range(args.pairs):
            sequence.extend(((workload, seed_xor, "aa", pair, 0, 0), (workload, seed_xor, "aa", pair, 1, 0)))
        for pair in range(args.pairs):
            arms = [0, 1] if pair % 2 == 0 else [1, 0]
            sequence.extend((workload, seed_xor, "paired", pair, position, mode) for position, mode in enumerate(arms))

    for workload, seed_xor, phase, pair, position, mode in sequence:
        record = run_one(binary, workload, seed_xor, args.threads, mode, phase, pair, position)
        if record["status"] == "ok":
            try:
                actual = signature(record)
            except (KeyError, TypeError) as exc:
                record["status"] = "output_error"
                record["error"] = str(exc)
            else:
                signatures = expected.setdefault(workload, {})
                if mode not in signatures:
                    signatures[mode] = actual
                    receipt["output_signatures"] = expected
                elif actual != signatures[mode]:
                    record["status"] = "output_mismatch"
                    record["expected_signature"] = signatures[mode]
                    record["actual_signature"] = actual
                if record["status"] == "ok":
                    primary = record["cases"][PRIMARY]
                    if not primary["direct_pack_used"]:
                        record["status"] = "primary_direct_pack_fallback"
                    elif not primary["direct_unpack_used"]:
                        record["status"] = "direct_unpack_route_mismatch"
                    elif mode == 0:
                        reference_record[workload] = record
                        if primary["column_cert_attempted"] or primary["selected_original_used"]:
                            record["status"] = "reference_route_mismatch"
                    else:
                        errors = compare_structural(reference_record[workload], record, 512)
                        if primary["output_terms"] > reference_record[workload]["cases"][PRIMARY]["output_terms"] * 0.25:
                            errors.append("primary output terms exceeded structural gate")
                        if errors:
                            record["status"] = "candidate_structural_mismatch"
                            record["errors"] = errors
        receipt["runs"].append(record)
        write_receipt(output, receipt)
        if record["status"] != "ok":
            receipt["status"] = record["status"]
            write_receipt(output, receipt)
            return 1

    receipt["summary"] = summarize(receipt["runs"], active_workloads)
    receipt["host"]["load_end"] = os.getloadavg()
    receipt["status"] = "complete"
    write_receipt(output, receipt)
    print(json.dumps({name: cases[PRIMARY] for name, cases in receipt["summary"].items()}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
