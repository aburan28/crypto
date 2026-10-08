#!/usr/bin/env python3
"""Bounded numerical positive controls; no curve, key, or DLP solver interfaces."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import statistics
import struct
import subprocess
import sys
import time
import tracemalloc

SCHEMA = 1
MODULUS = 257
VARIANTS = ("power_sum", "horner")
PHASES = ("initialize", "encode_input", "decode_input", "warmup", "compute",
          "encode_output", "verify")


def seed_for(*parts):
    return int.from_bytes(hashlib.sha256("/".join(map(str, parts)).encode()).digest()[:8], "big")


def cover_trial(variant, states, seed, budget):
    """Cover a ring using nearest-neighbor steps or independent uniform refreshes."""
    if variant not in ("nearest", "refresh") or states not in (32, 128) or budget < 1:
        raise ValueError("unsupported toy coverage configuration")
    rng = random.Random(seed)
    position = rng.randrange(states)
    seen = {position}
    for operations in range(1, budget + 1):
        if variant == "nearest":
            position = (position + (1 if rng.getrandbits(1) else -1)) % states
        else:
            position = rng.randrange(states)
        seen.add(position)
        if len(seen) == states:
            return {"status": "completed", "operations": operations, "covered": len(seen)}
    return {"status": "censored", "operations": budget, "covered": len(seen)}


def power_sum(coefficients, x):
    return sum(a * pow(x, i, MODULUS) for i, a in enumerate(coefficients)) % MODULUS


def horner(coefficients, x):
    result = 0
    for a in reversed(coefficients):
        result = (result * x + a) % MODULUS
    return result


def reference(coefficients, x):
    # Independent unbounded-integer expression, reduced only at the end.
    return sum(a * x ** i for i, a in enumerate(coefficients)) % MODULUS


def pack(values):
    return struct.pack(f"<{len(values)}H", *values)


def unpack(data):
    return struct.unpack(f"<{len(data) // 2}H", data)


def worker(config):
    variant, seed, count = config["variant"], config["seed"], config["count"]
    if variant not in VARIANTS or not 1 <= count <= 16384:
        raise ValueError("unsupported arithmetic configuration")
    memory_pass = config.get("memory_pass", False)
    if memory_pass:
        tracemalloc.start()
    started = time.perf_counter_ns()
    timings = {}

    def phase(name, fn):
        before = time.perf_counter_ns()
        result = fn()
        timings[name] = time.perf_counter_ns() - before
        return result

    def initialize():
        rng = random.Random(seed)
        return [rng.randrange(MODULUS) for _ in range(25)], [rng.randrange(MODULUS) for _ in range(count)]

    coefficients, xs = phase("initialize", initialize)
    encoded = phase("encode_input", lambda: pack(coefficients + xs))
    decoded = phase("decode_input", lambda: unpack(encoded))
    coefficients, xs = decoded[:25], decoded[25:]
    fn = power_sum if variant == "power_sum" else horner
    phase("warmup", lambda: [fn(coefficients, x) for x in xs[:64]])
    values = phase("compute", lambda: [fn(coefficients, x) for x in xs])
    output = phase("encode_output", lambda: pack(values))

    def verify():
        oracle = [reference(coefficients, x) for x in range(MODULUS)]
        return (all(fn(coefficients, x) == oracle[x] for x in range(MODULUS))
                and all(value == oracle[x] for x, value in zip(xs, values))
                and tuple(values) == unpack(output))

    verified = phase("verify", verify)
    elapsed = time.perf_counter_ns() - started
    memory = {"python_traced_peak_bytes": None, "process_peak_rss_bytes": None}
    if memory_pass:
        memory["python_traced_peak_bytes"] = tracemalloc.get_traced_memory()[1]
        tracemalloc.stop()
        try:
            import resource
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            memory["process_peak_rss_bytes"] = rss if sys.platform == "darwin" else rss * 1024
        except ImportError:
            pass
    return {"status": "completed" if verified else "incorrect", "verified": verified,
            "variant": variant, "seed": seed, "count": count, "memory_pass": memory_pass,
            "phases_ns": timings, "worker_ns": elapsed,
            "unattributed_worker_ns": elapsed - sum(timings.values()),
            "input_bytes": len(encoded), "output_bytes": len(output),
            "input_sha256": hashlib.sha256(encoded).hexdigest(),
            "output_sha256": hashlib.sha256(output).hexdigest(), "memory": memory}


def launch(config, timeout):
    started = time.perf_counter_ns()
    try:
        proc = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker"],
                              input=json.dumps(config), text=True, capture_output=True,
                              timeout=timeout, check=False)
    except subprocess.TimeoutExpired:
        return {**config, "status": "timeout", "verified": False,
                "parent_wall_ns": time.perf_counter_ns() - started}
    wall = time.perf_counter_ns() - started
    if proc.returncode:
        return {**config, "status": "failed", "verified": False, "parent_wall_ns": wall,
                "returncode": proc.returncode, "stderr": proc.stderr[-2000:]}
    try:
        result = json.loads(proc.stdout)
        if not isinstance(result, dict) or result.get("variant") != config["variant"]:
            raise ValueError("invalid worker receipt")
    except (ValueError, TypeError):
        return {**config, "status": "invalid_receipt", "verified": False, "parent_wall_ns": wall}
    wall = time.perf_counter_ns() - started
    result["parent_wall_ns"] = wall
    # Includes interpreter startup/teardown, IPC, receipt parsing, and scheduling.
    result["process_overhead_ns"] = wall - result["worker_ns"]
    return result


def bootstrap_ratio(pairs, seed, samples=1000):
    """Paired bootstrap of ratio of sums; positive values required."""
    if not pairs:
        return None
    if any(a <= 0 or b <= 0 for a, b in pairs):
        raise ValueError("ratios require positive observations")
    rng = random.Random(seed)
    ratios = []
    for _ in range(samples):
        draws = [pairs[rng.randrange(len(pairs))] for _ in pairs]
        ratios.append(sum(a for a, _ in draws) / sum(b for _, b in draws))
    ratios.sort()
    return {"ratio": sum(a for a, _ in pairs) / sum(b for _, b in pairs),
            "ci95": [ratios[int(.025 * samples)], ratios[min(samples - 1, int(.975 * samples))]],
            "pairs": len(pairs), "bootstrap_samples": samples}


def coverage_summary(rows, budget):
    completed = sorted(r["operations"] for r in rows if r["status"] == "completed")
    # Quantiles of the full empirical completion CDF; unreached quantiles remain unknown.
    def quantile(q):
        index = math.ceil(q * len(rows)) - 1
        return completed[index] if index < len(completed) else None
    return {"trials": len(rows), "completed": len(completed),
            "censored": len(rows) - len(completed), "budget": budget,
            "restricted_mean_operations": statistics.mean(r["operations"] for r in rows),
            "completion_median": quantile(.5), "completion_p90": quantile(.9)}


def summarize(rows, budget):
    summary = {"coverage": {}, "arithmetic": {}}
    for split in ("discovery", "holdout"):
        for states in (32, 128):
            groups = {v: [r for r in rows if r["study"] == "coverage" and r["split"] == split
                          and r["states"] == states and r["variant"] == v]
                      for v in ("nearest", "refresh")}
            pairs = [(a["operations"], b["operations"]) for a, b in zip(groups["nearest"], groups["refresh"])]
            summary["coverage"][f"{split}/{states}"] = {
                **{v: coverage_summary(rs, budget) for v, rs in groups.items()},
                "restricted_mean_ratio": bootstrap_ratio(pairs, seed_for(split, states, "bootstrap")),
                "claim": "finite coverage positive control; no DLP inference"}
        selected = [r for r in rows if r["study"] == "arithmetic" and r["split"] == split]
        timings = [r for r in selected if not r["memory_pass"]]
        groups = {v: [r for r in timings if r["variant"] == v] for v in VARIANTS}
        for rs in groups.values():
            rs.sort(key=lambda r: r["trial"])
        valid = bool(timings) and all(r["status"] == "completed" and r["verified"] for r in selected)
        valid = valid and len(groups["power_sum"]) == len(groups["horner"])
        paired = list(zip(groups["power_sum"], groups["horner"]))
        if valid:
            valid = all(a["trial"] == b["trial"] and a["input_sha256"] == b["input_sha256"]
                        and a["output_sha256"] == b["output_sha256"] for a, b in paired)
        result = {"timing_trials": len(timings), "all_verified_and_paired": valid,
                  "statuses": {status: sum(r["status"] == status for r in selected)
                               for status in sorted({r["status"] for r in selected})},
                  "kernel_ratio": None, "full_wall_ratio": None,
                  "decision": "inconclusive", "memory_passes": [r for r in selected if r["memory_pass"]]}
        if valid:
            result["kernel_ratio"] = bootstrap_ratio(
                [(a["phases_ns"]["compute"], b["phases_ns"]["compute"]) for a, b in paired], seed_for(split, "kernel"))
            result["full_wall_ratio"] = bootstrap_ratio(
                [(a["parent_wall_ns"], b["parent_wall_ns"]) for a, b in paired], seed_for(split, "wall"))
            result["median_phases_ms"] = {v: {p: statistics.median(r["phases_ns"][p] for r in rs) / 1e6
                                                   for p in PHASES} for v, rs in groups.items()}
            result["median_parent_wall_ms"] = {v: statistics.median(r["parent_wall_ns"] for r in rs) / 1e6
                                                for v, rs in groups.items()}
            result["decision"] = ("finite_end_to_end_improvement" if result["full_wall_ratio"]["ci95"][0] > 1
                                  else "no_end_to_end_improvement_established")
        summary["arithmetic"][split] = result
    return summary


def run_studies(output, trials=128, arithmetic_trials=24, count=4096, budget=8192, timeout=15):
    if not (2 <= trials <= 1024 and 2 <= arithmetic_trials <= 128 and 1 <= count <= 16384
            and 1 <= budget <= 65536 and 0 < timeout <= 60):
        raise ValueError("study parameters outside bounded limits")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter_ns()
    manifest = {"schema": SCHEMA, "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "python": sys.version, "platform": platform.platform(), "machine": platform.machine(),
                "cpu_count": os.cpu_count(), "trials_per_split_and_size": trials,
                "arithmetic_pairs_per_split": arithmetic_trials, "arithmetic_count": count,
                "coverage_budget": budget, "child_timeout_s": timeout,
                "scope": "32/128-state ring coverage and degree-24 polynomials modulo 257",
                "transfer_scope": "host encode/decode and subprocess IPC only; no GPU transfers measured",
                "selection": "fixed variants; discovery never changes holdout configuration",
                "censoring": "fixed operation budget; report E[min(T,budget)], not mean completion time",
                "timing": "fresh child per arithmetic variant; randomized within-pair order",
                "memory": "separate traced passes; process RSS includes interpreter; no timing claims from traced passes",
                "novelty": "known positive controls; no new cryptanalytic algorithm or solver speedup"}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    rows = []
    with (output / "trials.jsonl").open("x") as ledger:
        def record(row):
            rows.append(row)
            ledger.write(json.dumps(row, sort_keys=True) + "\n")
            ledger.flush()

        for split in ("discovery", "holdout"):
            print(f"{split}: coverage", flush=True)
            for states in (32, 128):
                for trial in range(trials):
                    seed = seed_for("coverage", split, states, trial)
                    for variant in ("nearest", "refresh"):
                        record({"study": "coverage", "split": split, "states": states,
                                "trial": trial, "seed": seed, "variant": variant,
                                **cover_trial(variant, states, seed, budget)})
            print(f"{split}: arithmetic and full-cost measurements", flush=True)
            for trial in range(arithmetic_trials):
                order = list(VARIANTS)
                random.Random(seed_for(split, trial, "order")).shuffle(order)
                for position, variant in enumerate(order):
                    config = {"variant": variant, "seed": seed_for("arithmetic", split, trial),
                              "count": count, "memory_pass": False}
                    record({"study": "arithmetic", "split": split, "trial": trial, "order": position,
                            **launch(config, timeout)})
            for variant in VARIANTS:
                config = {"variant": variant, "seed": seed_for("memory", split),
                          "count": count, "memory_pass": True}
                record({"study": "arithmetic", "split": split, "trial": "memory", **launch(config, timeout)})
    summary = summarize(rows, budget)
    summary["suite_before_report_write_ns"] = time.perf_counter_ns() - started
    summary["trials_sha256"] = hashlib.sha256((output / "trials.jsonl").read_bytes()).hexdigest()
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--trials", type=int, default=128)
    parser.add_argument("--arithmetic-trials", type=int, default=24)
    parser.add_argument("--count", type=int, default=4096)
    parser.add_argument("--budget", type=int, default=8192)
    parser.add_argument("--timeout", type=float, default=15)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(json.load(sys.stdin))))
    elif args.output:
        result = run_studies(args.output, args.trials, args.arithmetic_trials, args.count, args.budget, args.timeout)
        if not all(r["all_verified_and_paired"] for r in result["arithmetic"].values()):
            sys.exit(1)
    else:
        parser.error("--output is required")


if __name__ == "__main__":
    main()
