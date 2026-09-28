#!/usr/bin/env python3
"""Verify completeness and summarize the preregistered bounded experiment."""
from __future__ import annotations
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import random
import statistics
import sys


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bootstrap(values):
    rng = random.Random(20260922)
    medians = sorted(statistics.median(rng.choices(values, k=len(values))) for _ in range(4000))
    return [medians[100], medians[3899]]


def main(root: Path):
    protocol, metadata = read(root / "protocol.json"), read(root / "metadata.json")
    assert metadata["complete"], "incomplete campaign cannot be promoted"
    for name, expected in metadata["source_hashes"].items():
        assert sha(root / name) == expected, f"source mismatch: {name}"
    test = read(root / "test_receipt.json")
    assert test["exit_code"] == 0 and not test["timed_out"]
    for stream in ("stdout", "stderr"):
        assert hashlib.sha256(test[stream].encode()).hexdigest() == test[f"{stream}_sha256"]
    expected_cells = {
        f"n{n}-{split}-{seed}-{family}-b{batch}"
        for n in protocol["variables"]
        for split in ("discovery", "holdout")
        for seed in protocol[f"{split}_seeds"]
        for family in protocol["families"]
        for batch in protocol["batches"]
    }
    receipts = read(root / "receipts.json")
    assert len(receipts) == len(expected_cells)
    receipts = {r["cell"]: r for r in receipts}
    assert set(receipts) == expected_cells
    raw = defaultdict(list)
    for n in protocol["variables"]:
        for line in (root / f"raw-n{n}.jsonl").read_text().splitlines():
            row = json.loads(line)
            decoded = json.loads(row["raw_line"])
            assert decoded == {k: v for k, v in row.items() if k not in ("cell", "split", "raw_line")}
            raw[row["cell"]].append(row)
    assert set(raw) == expected_cells
    summaries, by_group = [], defaultdict(list)
    for cell, rows in sorted(raw.items()):
        receipt = receipts[cell]
        assert receipt["exit_code"] == 0 and not receipt["timed_out"]
        reconstructed = "".join(r["raw_line"] for r in rows).encode()
        assert hashlib.sha256(reconstructed).hexdigest() == receipt["stdout_sha256"], cell
        assert receipt["stderr_sha256"] == hashlib.sha256(b"").hexdigest(), cell
        assert rows[0]["type"] == "fixture"
        fixture = rows[0]
        for field in ("n", "seed", "family", "batch", "split"):
            assert fixture[field] == receipt[field]
        samples = rows[1:]
        assert len(samples) == protocol["repetitions"] * len(protocol["variants"])
        pairs = {(s["rep"], s["variant"]): s for s in samples}
        assert len(pairs) == len(samples)
        assert set(pairs) == {(rep, arm) for rep in range(protocol["repetitions"]) for arm in protocol["variants"]}
        for rep in range(protocol["repetitions"]):
            ordered = [s["variant"] for s in samples if s["rep"] == rep]
            names = protocol["variants"]
            assert ordered == [names[(rep + i) % len(names)] for i in range(len(names))]
        for sample in samples:
            assert sample["verified_outputs"] == fixture["batch"]
            assert sample["total_ns"] >= sample["setup_ns"] + sample["apply_ns"] + sample["validation_ns"]
            if sample["variant"] in ("schedule", "matrix_cache"):
                expected_hits = fixture["batch"] if fixture["family"] == "repeat" else (fixture["batch"] + 1) // 2
                assert sample["hits"] == expected_hits
                assert sample["fallbacks"] == fixture["batch"] - expected_hits
                assert sample["changed_hits"] == 0
            by_group[(fixture["split"], fixture["family"], fixture["batch"], sample["variant"])].append(sample)
        outputs = {s["output_bytes"] for s in samples}
        assert len(outputs) == 1
        arms = {}
        for variant in protocol["variants"]:
            arm = [s for s in samples if s["variant"] == variant]
            arms[variant] = {f"median_{field}": statistics.median(s[field] for s in arm)
                             for field in ("setup_ns", "apply_ns", "validation_ns", "total_ns", "retained_bytes")}
            arms[variant].update(hits=arm[0]["hits"], fallbacks=arm[0]["fallbacks"], changed_hits=arm[0]["changed_hits"])
        summaries.append({"cell": cell, "split": fixture["split"], "family": fixture["family"],
                          "n": fixture["n"], "batch": fixture["batch"], "seed": fixture["seed"],
                          "fixture_sha256": hashlib.sha256(fixture["raw_line"].encode()).hexdigest(),
                          "arms": arms, "peak_rss_bytes": receipt["peak_rss_bytes"]})

    gates = []
    for n in protocol["variables"]:
        for batch in (16, 64):
            group = [cell for cell in summaries if cell["split"] == "holdout" and cell["family"] == "repeat"
                     and cell["n"] == n and cell["batch"] == batch]
            for control in ("layout", "matrix_cache"):
                ratios = []
                for cell in group:
                    pairs = {(s["rep"], s["variant"]): s for s in raw[cell["cell"]][1:]}
                    ratios.extend(pairs[rep, control]["total_ns"] / pairs[rep, "schedule"]["total_ns"]
                                  for rep in range(protocol["repetitions"]))
                interval = bootstrap(ratios)
                gates.append({"n": n, "batch": batch, "control": control,
                              "paired_ratio_median": statistics.median(ratios),
                              "ci95_paired_median": interval, "threshold": 1.05, "pass": interval[0] > 1.05})
    result = {
        "schema_version": 1, "scope": protocol["scope"],
        "mathematical_finding": protocol["mathematical_boundary"],
        "correctness": "PASS", "cells": len(summaries),
        "paired_batch_samples": sum(len(rows) - 1 for rows in raw.values()),
        "verified_outputs": sum(s["verified_outputs"] for rows in raw.values() for s in rows[1:]),
        "exact_schedule_changed_input_hits": sum(s["changed_hits"] for rows in raw.values() for s in rows[1:] if s["variant"] == "schedule"),
        "cross_instance_reuse": "REJECTED: exact full-support signatures cannot hit changed canonical systems",
        "performance_gate": "PASS" if all(g["pass"] for g in gates) else "REJECTED",
        "gate_details": gates,
        "ci_scope": "paired bootstrap over measured repetitions and two fixed holdout seeds; not a population claim",
        "memory_scope": protocol["memory"],
        "classification": protocol["classification"],
        "comparisons": summaries,
    }
    (root / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    lines = ["# Boolean symbolic-product schedule experiment", "", result["mathematical_finding"], "",
             f"Correctness: **PASS** across {result['cells']} fixed cells, {result['paired_batch_samples']} batch-arm samples, and {result['verified_outputs']} verified matrix outputs.", "",
             f"Cross-instance schedule hits: **{result['exact_schedule_changed_input_hits']}**. Performance promotion against both retained controls: **{result['performance_gate']}**.", "",
             "The table below reports cold batch construction plus output validation in milliseconds, at batch 64 on holdout fixtures. Values are medians over four variable sizes, two seeds and eight repetitions; the per-size acceptance intervals remain in results.json. These are standalone construction diagnostics, not solver or cryptanalytic performance.", "",
             "| Family | Variant | Cold batch median (ms) | Retained structure median (bytes) | Median hits / 64 |", "|---|---|---:|---:|---:|"]
    for family in protocol["families"]:
        for variant in protocol["variants"]:
            samples = by_group["holdout", family, 64, variant]
            lines.append(f"| {family} | {variant} | {statistics.median(s['total_ns'] for s in samples)/1e6:.6f} | {statistics.median(s['retained_bytes'] for s in samples):.0f} | {statistics.median(s['hits'] for s in samples):.0f} |")
    lines += ["", "Cache setup, exact-key checks, fallbacks, output allocation, output validation and destruction are charged. Common fixture and independent reference generation are outside timing; fresh-process receipts cover the whole worker, including those costs. Retained bytes exclude allocator metadata; RSS includes all variants and the common reference corpus.", "",
              "Changing a coefficient from 1 to 0 over F2 changes support. Reusing a schedule across such changes would require a different, parameterized support-envelope contract with explicit cancellation and degree-drop rules. This experiment does not implement or validate that different contract.", "",
              "No production solver path is changed. The full-support cache key is evaluated exactly as proposed, including the stronger packed-matrix-cache control. No curve targets, key work, relation collection, or rho comparison are part of this study."]
    (root / "RESULT.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({key: result[key] for key in ("correctness", "cells", "paired_batch_samples", "verified_outputs", "performance_gate", "exact_schedule_changed_input_hits")}, sort_keys=True))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
