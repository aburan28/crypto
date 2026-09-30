#!/usr/bin/env python3
"""Two-seed exact-output gate for opt-in F5 byte-colex row packing."""

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SEEDS = ("0", "badc0de1")
CASES = {
    "f5_n12_m12_d4", "f5_n16_m16_d3", "f5_n16_m16_d4",
    "f5_n20_m20_d3", "f5_n20_m20_d4", "f5_n24_m24_d3",
    "f5_n24_m24_d4",
}
EXACT = (
    "rows_fp", "row_space_fp", "output_terms", "rank", "cols",
    "rows_built", "rows_f4", "rows_pruned", "criterion_word_ops",
    "reduce_word_ops", "direct_pack_used", "direct_unpack_used",
)
BASE = {
    "KIC_F5_ECHELON": "2", "KIC_F5_FUSED_BUILD": "1",
    "KIC_F5_DIRECT_PACK": "1", "KIC_F5_UNPACK_DIRECT": "1",
    "KIC_GF2_FORCE_AVX2": "1", "KIC_GF2_REUSE_TABLE": "1",
    "KIC_GF2_TABLES": "4", "KIC_GF2_SIMD": "1",
    "KIC_GF2_BRANCHLESS_STRIP": "1", "KIC_GF2_DEFER_ABOVE": "0",
    "KIC_GF2_WORD_BATCH": "0", "KIC_F5_AVX512_UNPACK": "0",
    "KIC_GF2_AVX2_TABLE_BUILD": "0", "RAYON_NUM_THREADS": "1",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, sort_keys=True, indent=2) + "\n")
    temporary.replace(path)


def run(binary, seed, byte_colex):
    env = os.environ.copy()
    config = {**BASE, "KIC_F5_COLEX_BYTES": str(byte_colex)}
    env.update(config)
    command = [str(binary), "1", "24", "f5", seed]
    record = {
        "seed": seed, "byte_colex": byte_colex, "command": command,
        "environment": config, "load_before": os.getloadavg(),
    }
    started = time.monotonic()
    try:
        proc = subprocess.run(command, env=env, text=True, capture_output=True, timeout=120)
        record.update(
            process_wall_ms=(time.monotonic() - started) * 1000,
            exit_code=proc.returncode, stdout=proc.stdout, stderr=proc.stderr,
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
            status="timeout", process_wall_ms=(time.monotonic() - started) * 1000,
            stdout=exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else exc.stdout,
            stderr=exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else exc.stderr,
            load_after=os.getloadavg(),
        )
    except (json.JSONDecodeError, KeyError, TypeError) as exc:
        record.update(status="output_error", error=str(exc), load_after=os.getloadavg())
    return record


def compare(reference, candidate):
    errors = []
    if reference["status"] != "ok" or candidate["status"] != "ok":
        return ["incomplete process"]
    for name in sorted(CASES):
        a, b = reference["cases"][name], candidate["cases"][name]
        for field in EXACT:
            if a[field] != b[field]:
                errors.append(f"{name}: {field} mismatch")
        if a["byte_colex_used"]:
            errors.append(f"{name}: reference selected byte colex")
        if b["byte_colex_used"] != b["direct_pack_used"]:
            errors.append(f"{name}: byte colex route mismatch")
    return errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    report = {
        "status": "running", "protocol": "PROTOCOL.md", "seeds": SEEDS,
        "host": {
            "platform": platform.platform(), "machine": platform.machine(),
            "processor": platform.processor(), "cpu_count": os.cpu_count(),
            "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
            "load_start": os.getloadavg(),
        },
        "git_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "binary_sha256": sha(binary),
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in (
            ROOT / "src/cryptanalysis/koblitz_groebner.rs",
            ROOT / "src/cryptanalysis/matrix_f5_f2.rs",
            ROOT / "examples/f4_f2_bench.rs")},
        "runs": [],
    }
    save(output, report)
    for seed in SEEDS:
        for byte_colex in (0, 1):
            report["runs"].append(run(binary, seed, byte_colex))
            save(output, report)
    results = {}
    for seed in SEEDS:
        reference, candidate = [r for r in report["runs"] if r["seed"] == seed]
        results[seed] = compare(reference, candidate)
    report["exactness_errors"] = results
    report["host"]["load_end"] = os.getloadavg()
    report["status"] = (
        "advance_to_local_pairs" if all(not errors for errors in results.values())
        else "reject_exactness"
    )
    save(output, report)
    print(output)
    print(json.dumps({"status": report["status"], "errors": results}, sort_keys=True))
    return 0 if report["status"] == "advance_to_local_pairs" else 1


if __name__ == "__main__":
    sys.exit(main())
