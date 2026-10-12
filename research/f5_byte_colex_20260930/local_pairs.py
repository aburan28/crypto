#!/usr/bin/env python3
"""Unqualified Apple ARM64 A/A and A/B screen for F5 byte-colex packing."""

import argparse
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

from structural_screen import BASE, CASES, EXACT, SEEDS, sha, save

ROOT = Path(__file__).resolve().parents[2]
PRIMARY = "f5_n24_m24_d4"


def run(binary, seed, arm, label):
    enabled = arm == "B"
    config = {**BASE, "KIC_F5_COLEX_BYTES": str(int(enabled))}
    env = os.environ.copy()
    env.update(config)
    command = [str(binary), "1", "24", "f5", seed]
    record = {
        "seed": seed, "arm": arm, "label": label, "command": command,
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


def ratio(reference, other, case, field):
    return reference["cases"][case][field] / other["cases"][case][field]


def pair_summary(pairs):
    return {
        case: {
            "full_call_ratios": [ratio(a, b, case, "wall_ms") for a, b in pairs],
            "build_ratios": [ratio(a, b, case, "f5_build_ms") for a, b in pairs],
            "full_call_median": statistics.median(
                ratio(a, b, case, "wall_ms") for a, b in pairs
            ),
            "build_median": statistics.median(
                ratio(a, b, case, "f5_build_ms") for a, b in pairs
            ),
        }
        for case in sorted(CASES)
    }


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
        "qualified_isolation": False,
        "isolation_reason": "Apple ARM64 lacks Linux affinity/PSI reservation in tools/isolated_bench.py; shared busy lock only",
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
        "runs": [], "pairs": [],
    }
    save(output, report)
    for seed in SEEDS:
        for arm in ("A", "B"):
            report["runs"].append(run(binary, seed, arm, f"warmup_{arm}"))
            save(output, report)
        for round_number in range(5):
            indices = []
            for position in range(2):
                indices.append(len(report["runs"]))
                report["runs"].append(run(binary, seed, "A", f"aa_{round_number}_{position}"))
                save(output, report)
            report["pairs"].append({"seed": seed, "kind": "AA", "round": round_number,
                                    "reference": indices[0], "other": indices[1]})
            save(output, report)
        for round_number in range(5):
            order = ("A", "B") if round_number % 2 == 0 else ("B", "A")
            indices = {}
            for position, arm in enumerate(order):
                indices[arm] = len(report["runs"])
                report["runs"].append(run(binary, seed, arm, f"ab_{round_number}_{position}"))
                save(output, report)
            report["pairs"].append({"seed": seed, "kind": "AB", "round": round_number,
                                    "reference": indices["A"], "other": indices["B"]})
            save(output, report)
    errors = []
    summaries = {}
    for seed in SEEDS:
        baseline = next(r for r in report["runs"] if r["seed"] == seed and r["arm"] == "A")
        if baseline["status"] != "ok":
            errors.append(f"{seed}/warmup_A: {baseline['status']}")
            continue
        for record in (r for r in report["runs"] if r["seed"] == seed):
            if record["status"] != "ok":
                errors.append(f"{seed}/{record['label']}: {record['status']}")
                continue
            for case in CASES:
                expected, actual = baseline["cases"][case], record["cases"][case]
                for field in EXACT:
                    if expected[field] != actual[field]:
                        errors.append(f"{seed}/{record['label']}/{case}: {field} mismatch")
                if actual["byte_colex_used"] != (record["arm"] == "B" and actual["direct_pack_used"]):
                    errors.append(f"{seed}/{record['label']}/{case}: route mismatch")
        if errors:
            continue
        per_kind = {}
        for kind in ("AA", "AB"):
            pairs = [(report["runs"][p["reference"]], report["runs"][p["other"]])
                     for p in report["pairs"] if p["seed"] == seed and p["kind"] == kind]
            per_kind[kind] = pair_summary(pairs)
        summaries[seed] = per_kind
    report["exactness_errors"] = errors
    report["summaries"] = summaries
    report["host"]["load_end"] = os.getloadavg()
    if errors:
        report["status"] = "reject_exactness"
    elif all(summaries[seed]["AB"][PRIMARY]["full_call_median"] < 0.95 for seed in SEEDS):
        report["status"] = "reject_local_screen"
    else:
        report["status"] = "advance_to_linux_qualified_pairs"
    save(output, report)
    print(output)
    print(json.dumps({"status": report["status"], "errors": errors,
                      "primary": {seed: summaries.get(seed, {}).get("AB", {}).get(PRIMARY)
                                  for seed in SEEDS}}, sort_keys=True, indent=2))
    return 0 if report["status"] == "advance_to_linux_qualified_pairs" else 1


if __name__ == "__main__":
    sys.exit(main())
