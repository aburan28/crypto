#!/usr/bin/env python3
"""Verify saved result accounting and source/input/output provenance, without Sage."""
import hashlib
import json
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parent


def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def main():
    r = json.loads((ROOT / "results/run-02.json").read_text())
    assert r["status"] == "passed"
    for key, name in (("source_sha256", "benchmark.py"), ("protocol_sha256", "PROTOCOL.md")):
        assert r[key] == hashlib.sha256((ROOT / name).read_bytes()).hexdigest(), name
    assert r["exhaustive"]["cases"] == 2310
    assert r["exhaustive"]["candidate_equalities"] == 9240
    assert {(p["m"], p["seed"]) for p in r["panels"]} == {
        (m, seed) for m in (31, 83, 131) for seed in (20260929, 20260930)}
    methods = {"binary", "binary_naf", "tau_naf", "reduced_tau_naf", "sage_native"}
    for p in r["panels"]:
        assert len(p["cases"]) == 24 and len(p["rounds"]) == 7 and len(p["aa"]) == 5
        assert p["input_sha256"] == digest(p["cases"])
        for pair in p["aa"]:
            assert pair["a"]["output_sha256"] == pair["b"]["output_sha256"] == p["output_sha256"]
            assert pair["ratio_a_over_b"] == pair["a"]["ns"] / pair["b"]["ns"]
        for round_ in p["rounds"]:
            assert set(round_["order"]) == set(round_["measurements"]) == methods
            for v in round_["measurements"].values():
                assert v["output_sha256"] == p["output_sha256"]
                assert v["ns_per_scalar"] == v["ns"] / 24
        for method, s in p["summary"].items():
            samples = [q["measurements"][method]["ns_per_scalar"] for q in p["rounds"]]
            assert s["median_ns"] == statistics.median(samples)
            assert s["min_ns"] == min(samples)
            ratios = [q["measurements"][method]["ns"] / q["measurements"]["binary_naf"]["ns"] for q in p["rounds"]]
            assert s["median_paired_cost_ratio_to_binary_naf"] == statistics.median(ratios)
            if method != "sage_native":
                assert len(s["counted_per_case"]) == len(s["digits_per_case"]) == 24
                for k, v in s["mean_counts"].items():
                    assert v == statistics.mean(row[k] for row in s["counted_per_case"])
                for k, v in s["mean_digits"].items():
                    assert v == statistics.mean(row[k] for row in s["digits_per_case"])
        noise = max(abs(a["ratio_a_over_b"]-1) for a in p["aa"])
        ratio = p["summary"]["reduced_tau_naf"]["median_paired_cost_ratio_to_binary_naf"]
        assert p["aa_max_deviation"] == noise
        assert p["arithmetic_screen_passed"] == (ratio <= .9 and 1-ratio > noise)
    for key in ("end_to_end_speedup", "gpu_throughput", "whole_walk_cost", "ecdlp_operations", "eligible_runtime_fraction"):
        assert r[key] is None
    failed = json.loads((ROOT / "results/run-01.json").read_text())
    assert failed["status"] == "failed" and "integer_representation" in failed["error"]
    print("PASS: six panels, all paired digests/counters, source hashes, retained failure; end-to-end fields null")


if __name__ == "__main__":
    main()
