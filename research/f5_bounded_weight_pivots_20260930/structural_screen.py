#!/usr/bin/env python3
"""Frozen two-seed term-count gate for bounded F5 pivot choices."""

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SEEDS = ("0", "badc0de1")
CASES = {
    "f5_n12_m12_d4", "f5_n16_m16_d3", "f5_n16_m16_d4",
    "f5_n20_m20_d3", "f5_n20_m20_d4", "f5_n24_m24_d3",
    "f5_n24_m24_d4",
}
PRIMARY = "f5_n24_m24_d4"
INVARIANTS = (
    "row_space_fp", "rank", "criterion_word_ops", "rows_f4",
    "rows_built", "cols", "rows_pruned",
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
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    temp.replace(path)


def run(binary, seed, choices, phase):
    env = os.environ.copy()
    env.update(BASE)
    env["KIC_GF2_PIVOT_CHOICES"] = str(choices)
    command = [str(binary), "1", "24", "f5", seed]
    record = {
        "seed": seed, "choices": choices, "phase": phase,
        "command": command, "environment": {**BASE, "KIC_GF2_PIVOT_CHOICES": str(choices)},
        "load_before": os.getloadavg(),
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


def exact_signature(record):
    fields = (*INVARIANTS, "rows_fp", "output_terms", "reduce_word_ops",
              "bounded_weight_pivot_used", "direct_pack_used", "direct_unpack_used")
    return {name: {field: case[field] for field in fields}
            for name, case in record["cases"].items()}


def matching(reference, candidate):
    errors = []
    if reference["status"] != "ok" or candidate["status"] != "ok":
        return ["incomplete run"]
    for name in CASES:
        a, b = reference["cases"][name], candidate["cases"][name]
        for field in INVARIANTS:
            if a[field] != b[field]:
                errors.append(f"{name}: {field} mismatch")
        expected_route = name == PRIMARY
        if a["bounded_weight_pivot_used"]:
            errors.append(f"{name}: reference selected bounded route")
        if b["bounded_weight_pivot_used"] != expected_route:
            errors.append(f"{name}: bounded route mismatch")
        if not b["bounded_weight_pivot_used"]:
            for field in ("rows_fp", "output_terms", "reduce_word_ops"):
                if a[field] != b[field]:
                    errors.append(f"{name}: inactive {field} changed")
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
        "host": {"platform": platform.platform(), "machine": platform.machine(),
                 "processor": platform.processor(), "cpu_count": os.cpu_count(),
                 "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
                 "load_start": os.getloadavg()},
        "git_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "binary_sha256": sha(binary),
        "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in (
            ROOT / "src/cryptanalysis/matrix_f5_f2.rs",
            ROOT / "src/cryptanalysis/gf2_elim.rs",
            ROOT / "examples/f4_f2_bench.rs")},
        "runs": [],
    }
    save(output, report)
    for seed in SEEDS:
        for choices, phase in ((0, "reference_k8"), (8, "candidate_k8"),
                               (0, "reference_k16"), (16, "candidate_k16")):
            record = run(binary, seed, choices, phase)
            report["runs"].append(record)
            save(output, report)
    if any(record["status"] != "ok" for record in report["runs"]):
        report["status"] = "incomplete"
        save(output, report)
        print(output)
        return 1
    results = {}
    for seed in SEEDS:
        records = [record for record in report["runs"] if record["seed"] == seed]
        baseline = records[0]
        if exact_signature(baseline) != exact_signature(records[2]):
            report["status"] = "reference_mismatch"
            save(output, report)
            print(output)
            return 1
        results[seed] = {}
        for choices, reference, candidate in ((8, records[0], records[1]),
                                              (16, records[2], records[3])):
            errors = matching(reference, candidate)
            terms = candidate["cases"][PRIMARY]["output_terms"]
            reference_terms = reference["cases"][PRIMARY]["output_terms"]
            results[seed][str(choices)] = {
                "errors": errors, "reference_terms": reference_terms,
                "candidate_terms": terms, "term_ratio": terms / reference_terms,
                "reference_rank": reference["cases"][PRIMARY]["rank"],
                "candidate_rank": candidate["cases"][PRIMARY]["rank"],
                "reference_row_space_fp": reference["cases"][PRIMARY]["row_space_fp"],
                "candidate_row_space_fp": candidate["cases"][PRIMARY]["row_space_fp"],
            }
    report["structural_results"] = results
    selected = None
    for choices in (8, 16):
        if all(not results[seed][str(choices)]["errors"] and
               results[seed][str(choices)]["term_ratio"] <= 0.75 for seed in SEEDS):
            selected = choices
            break
    report["selected_choices"] = selected
    report["host"]["load_end"] = os.getloadavg()
    report["status"] = "advance_to_local_pairs" if selected else "reject_structural_gate"
    save(output, report)
    print(output)
    print(json.dumps({"status": report["status"], "selected_choices": selected,
                      "structural_results": results}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
