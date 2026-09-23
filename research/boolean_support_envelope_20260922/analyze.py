#!/usr/bin/env python3
"""Verify completeness and summarize the preregistered bounded experiment."""
from __future__ import annotations
from collections import defaultdict
import hashlib
import json
import math
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
        base, envelope, inputs = fixture["base"], fixture["envelope"], fixture["inputs"]
        n, degree, batch = fixture["n"], fixture["degree"], fixture["batch"]
        assert degree == protocol["matrix_degree"] and fixture["active"] == (1 << n) - 1
        assert len(base) == len(envelope) == n and len(inputs) == batch
        assert inputs[0] == base
        assert envelope == [sorted(set(p) | {0} | {1 << i for i in range(n)}) for p in base]
        assert all(len(e) <= protocol["limits"]["terms_per_generator"] for e in envelope)
        same = [value == base for value in inputs]
        inside = [len(value) == n and all(set(p) <= set(e) for p, e in zip(value, envelope)) for value in inputs]
        if fixture["family"] == "repeat":
            assert all(same)
        else:
            assert sum(same) == 1
            assert len({json.dumps(value) for value in inputs}) == batch
        assert inside == [not (fixture["family"] == "escape" and i % 4 == 3) for i in range(batch)]
        def eligible_count(poly):
            d = max((m.bit_count() for m in poly), default=degree + 1)
            return sum(math.comb(n, k) for k in range(min(n, degree - d) + 1))
        new_multipliers = sum(max(0, eligible_count(p) - eligible_count(b))
                              for value in inputs for b, p in zip(base, value))
        assert fixture["newly_required_multipliers"] == new_multipliers
        if fixture["family"] == "degree_cycle" and batch > 1:
            assert new_multipliers > 0
        else:
            assert new_multipliers == 0
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
            if sample["variant"] in ("schedule", "matrix_cache", "envelope", "direct"):
                hits = inside if sample["variant"] == "envelope" else same
                if sample["variant"] == "direct":
                    hits = [False] * batch
                assert sample["hits"] == sum(hits)
                assert sample["fallbacks"] == (0 if sample["variant"] == "direct" else batch - sum(hits))
                assert sample["changed_hits"] == sum(h and not s for h, s in zip(hits, same))
            if sample["variant"] == "envelope":
                assert sample["retained_bytes"] <= protocol["limits"]["plan_retained_bytes"]
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
                          "newly_required_multipliers": new_multipliers,
                          "arms": arms, "peak_rss_bytes": receipt["peak_rss_bytes"]})

    gates = []
    for n in protocol["variables"]:
        for batch in (16, 64):
            for family in ("coefficients", "degree_cycle"):
                group = [cell for cell in summaries if cell["split"] == "holdout" and cell["family"] == family
                         and cell["n"] == n and cell["batch"] == batch]
                for control in ("direct", "layout", "matrix_cache"):
                    ratios = []
                    for cell in group:
                        pairs = {(s["rep"], s["variant"]): s for s in raw[cell["cell"]][1:]}
                        ratios.extend(pairs[rep, control]["total_ns"] / pairs[rep, "envelope"]["total_ns"]
                                      for rep in range(protocol["repetitions"]))
                    interval = bootstrap(ratios)
                    gates.append({"n": n, "batch": batch, "family": family, "control": control,
                                  "paired_ratio_median": statistics.median(ratios),
                                  "ci95_paired_median": interval, "threshold": 1.05, "pass": interval[0] > 1.05})
    result = {
        "schema_version": 1, "scope": protocol["scope"],
        "mathematical_finding": protocol["mathematical_boundary"],
        "correctness": "PASS", "cells": len(summaries),
        "paired_batch_samples": sum(len(rows) - 1 for rows in raw.values()),
        "verified_outputs": sum(s["verified_outputs"] for rows in raw.values() for s in rows[1:]),
        "exact_schedule_changed_input_hits": sum(s["changed_hits"] for rows in raw.values() for s in rows[1:] if s["variant"] == "schedule"),
        "envelope_changed_input_hits": sum(s["changed_hits"] for rows in raw.values() for s in rows[1:] if s["variant"] == "envelope"),
        "envelope_fallbacks": sum(s["fallbacks"] for rows in raw.values() for s in rows[1:] if s["variant"] == "envelope"),
        "newly_required_multipliers_per_grid": sum(c["newly_required_multipliers"] for c in summaries),
        "cross_instance_reuse": "PASS: changing coefficients, cancellations and degree drops independently checked",
        "performance_gate": "PASS" if all(g["pass"] for g in gates) else "REJECTED",
        "gates_passed": sum(g["pass"] for g in gates),
        "gate_details": gates,
        "ci_scope": "paired bootstrap over measured repetitions and two fixed holdout seeds; not a population claim or independent reproduction",
        "memory_scope": protocol["memory"],
        "classification": protocol["classification"],
        "full_pipeline_cost": None, "rho_ratio": None,
        "comparisons": summaries,
    }
    (root / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    lines = ["# Parameterized Boolean support-envelope experiment", "", result["mathematical_finding"], "",
             f"Correctness: **PASS** across {result['cells']} fixed cells, {result['paired_batch_samples']} batch-arm samples, and {result['verified_outputs']} verified matrix outputs.", "",
             f"Changed-input envelope hits: **{result['envelope_changed_input_hits']}**; support-escape fallbacks: **{result['envelope_fallbacks']}**. Exact-support schedule changed hits: **{result['exact_schedule_changed_input_hits']}**.", "",
             f"Cold performance promotion: **{result['performance_gate']}**, {result['gates_passed']} / {len(gates)} gates passed. Newly eligible multipliers before cancellation, over the fixture grid without arm/repetition multiplication: **{result['newly_required_multipliers_per_grid']}**.", "",
             "The table reports cold batch construction plus output validation in milliseconds, at batch 64 on holdouts. Medians pool four sizes, two seeds and ten repetitions; paired intervals remain separate by size and family in results.json. Ratios below are descriptive ratios of pooled medians against direct construction and packed-matrix caching. These are construction-stage diagnostics.", "",
             "| Family | Variant | Cold batch (ms) | Direct / arm | Matrix cache / arm | Retained bytes | Hits / 64 | Correctness |",
             "|---|---|---:|---:|---:|---:|---:|---|"]
    for family in protocol["families"]:
        baseline = statistics.median(s['total_ns'] for s in by_group['holdout', family, 64, 'direct'])
        matrix = statistics.median(s['total_ns'] for s in by_group['holdout', family, 64, 'matrix_cache'])
        for variant in protocol["variants"]:
            samples = by_group["holdout", family, 64, variant]
            total = statistics.median(s['total_ns'] for s in samples)
            lines.append(f"| {family} | {variant} | {total/1e6:.6f} | {baseline/total:.3f} | {matrix/total:.3f} | {statistics.median(s['retained_bytes'] for s in samples):.0f} | {statistics.median(s['hits'] for s in samples):.0f} | PASS |")
    lines += ["", "Cold totals charge compilation, guards, applications, fallbacks, fresh output, validation and destruction. Common fixture/envelope generation and independent reference generation are outside timing; process receipts include them. Retained bytes exclude allocator metadata. Worker RSS covers all variants plus the reference corpus, not candidate-specific peak usage.", "",
              "The retained envelope explicitly represents varying coefficients. It is not a cached completed matrix. Product parity, actual generator degree, newly eligible multipliers and actual output columns are evaluated at each application. Support escapes leave the retained plan unchanged and rebuild directly.", "",
              "No production solver path changes. This is finite generic Boolean matrix construction, with no curve inputs, solving, relation collection, target import, scalar recovery or rho comparison. Full-pipeline costs and normalized cryptanalytic ratios remain null."]
    (root / "RESULT.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({key: result[key] for key in ("correctness", "cells", "paired_batch_samples", "verified_outputs", "performance_gate", "gates_passed", "envelope_changed_input_hits")}, sort_keys=True))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
