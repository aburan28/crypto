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
        n, family, batch = fixture["n"], fixture["family"], fixture["batch"]
        inputs = fixture["inputs"]
        expected_active = sum(1 << (i*n//8) for i in range(8)) if family == "restricted_cycle" else (1 << n)-1
        assert fixture["degree"] == 3 and fixture["active"] == expected_active
        assert len(inputs) == batch and len({json.dumps(x) for x in inputs}) == batch
        for i, polys in enumerate(inputs):
            assert len(polys) == protocol["fixture_generators"]
            for j, poly in enumerate(polys):
                assert poly == sorted(set(poly)) and len(poly) <= protocol["limits"]["terms_per_generator"]
                assert all(0 <= term < (1 << n) for term in poly)
                d = max((m.bit_count() for m in poly), default=-1)
                expected_degree = 2
                if j == 0 and family != "quadratic":
                    if i % 4 == 1: expected_degree = 1
                    if i % 4 == 3: expected_degree = -1
                    if family == "restricted_cycle" and i % 4 == 2: expected_degree = 0
                assert d == expected_degree
        names = protocol["bridge_variants"] if n == protocol["bridge_n"] else protocol["variants"]
        samples = rows[1:]
        assert len(samples) == protocol["repetitions"] * len(names)
        pairs = {(s["rep"],s["variant"]):s for s in samples}
        assert len(pairs) == len(samples)
        assert set(pairs) == {(rep,arm) for rep in range(protocol["repetitions"]) for arm in names}
        for rep in range(protocol["repetitions"]):
            selected = [s for s in samples if s["rep"] == rep]
            assert [s["order"] for s in selected] == list(range(len(names)))
            assert [s["variant"] for s in selected] == [names[(rep+i)%len(names)] for i in range(len(names))]
        for sample in samples:
            arm = sample["variant"]
            assert sample["verified_outputs"] == batch
            assert sample["total_ns"] >= sample["setup_ns"] + sample["apply_ns"] + sample["validation_ns"]
            assert sample["hits"] == (0 if arm == "sorted" else batch)
            expected_entries = 4*(n+1) if arm in ("ranked","sparse_rank") else (1<<n) if arm == "dense" else 0
            assert sample["lookup_entries"] == expected_entries
            assert sample["total_rows"] <= batch * protocol["limits"]["rows"]
            assert sample["total_columns"] <= batch * protocol["limits"]["columns"]
            by_group[(fixture["split"],n,family,batch,arm)].append(sample)
        for field in ("output_bytes","total_rows","total_columns"):
            assert len({s[field] for s in samples}) == 1
        arms = {}
        for arm in names:
            values = [s for s in samples if s["variant"] == arm]
            arms[arm] = {f"median_{f}":statistics.median(s[f] for s in values)
                         for f in ("setup_ns","apply_ns","validation_ns","total_ns","retained_bytes")}
            arms[arm].update(lookup_entries=values[0]["lookup_entries"], hits=values[0]["hits"])
        summaries.append({"cell":cell,"n":n,"family":family,"batch":batch,"seed":fixture["seed"],"split":fixture["split"],
                          "arms":arms,"total_rows":samples[0]["total_rows"],"total_columns":samples[0]["total_columns"],
                          "output_bytes":samples[0]["output_bytes"],"peak_rss_bytes":receipt["peak_rss_bytes"],
                          "fixture_sha256":hashlib.sha256(fixture["raw_line"].encode()).hexdigest()})
    candidate = protocol.get("candidate_arm", "ranked")
    gates, bridge = [], []
    for n in protocol["variables"]:
        for family in protocol["families"]:
            cells = [c for c in summaries if c["split"] == "holdout" and c["n"] == n and c["family"] == family and c["batch"] == 32]
            controls = ["dense"] if n == protocol["bridge_n"] else protocol.get("comparison_controls", ["sorted","binary"])
            for control in controls:
                ratios = []
                for cell in cells:
                    pairs = {(s["rep"],s["variant"]):s for s in raw[cell["cell"]][1:]}
                    ratios.extend(pairs[rep,control]["total_ns"] / pairs[rep,candidate]["total_ns"] for rep in range(protocol["repetitions"]))
                interval = bootstrap(ratios)
                entry = {"n":n,"family":family,"batch":32,"control":control,"paired_ratio_median":statistics.median(ratios),"ci95_paired_median":interval}
                if control == "dense":
                    bridge.append(entry)
                else:
                    threshold = 2.0 if control == "sorted" else 1.05
                    entry.update(threshold=threshold, **{"pass":interval[0]>threshold})
                    gates.append(entry)
    result = {"schema_version":1,"scope":protocol["scope"],"correctness":"PASS","cells":len(summaries),
              "paired_batch_samples":sum(len(rows)-1 for rows in raw.values()),
              "verified_outputs":sum(s["verified_outputs"] for rows in raw.values() for s in rows[1:]),
              "dramatic_scaling_gate":"PASS" if all(g["pass"] for g in gates) else "REJECTED",
              "gates_passed":sum(g["pass"] for g in gates),"gate_details":gates,"bridge_comparisons":bridge,
              "dense_not_executed":[{"n":n,"status":"NOT_EXECUTED","total_ns":None,"reason":"predeclared dense lookup dimension cap"}
                                    for n in protocol["variables"] if n>protocol["bridge_n"]],
              "classification":protocol["classification"],"full_solver_cost":None,"rho_ratio":None,
              "ci_scope":"paired bootstrap over repetitions on two fixed holdout seeds; not a population result or independent reproduction",
              "memory_scope":protocol["memory"],"comparisons":summaries}
    if protocol.get("enable_sparse"):
        result["candidate_arm"] = candidate
    (root/"results.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    repetition_label = "twelve" if protocol["repetitions"] == 12 else str(protocol["repetitions"])
    lines = ["# Combinatorial monomial coordinates: scaling experiment","",
             f"Correctness **PASS**: {result['cells']} cells, {result['paired_batch_samples']} batch-arm samples, {result['verified_outputs']} oracle-verified outputs.","",
             f"Dramatic scalable-construction gate: **{result['dramatic_scaling_gate']}**, {result['gates_passed']} / {len(gates)} comparisons passed.","",
             f"Cold batch32 milliseconds below include setup, applications, allocations, exact equality validation and destruction. Values are medians over two holdout seeds and {repetition_label} balanced repetitions. Ratios here compare pooled medians; acceptance uses the per-size paired intervals in results.json.","",
             "| Variables | Family | Variant | Cold batch (ms) | Sorted / arm | Retained bytes | Lookup entries | Correctness |",
             "|---:|---|---|---:|---:|---:|---:|---|"]
    for n in protocol["variables"]:
        names = protocol["bridge_variants"] if n == protocol["bridge_n"] else protocol["variants"]
        for family in protocol["families"]:
            base = statistics.median(s["total_ns"] for s in by_group['holdout',n,family,32,'sorted'])
            for arm in names:
                samples = by_group['holdout',n,family,32,arm]
                total = statistics.median(s['total_ns'] for s in samples)
                lines.append(f"| {n} | {family} | {arm} | {total/1e6:.6f} | {base/total:.3f} | {statistics.median(s['retained_bytes'] for s in samples):.0f} | {samples[0]['lookup_entries']} | PASS |")
    lines += ["","The dense lookup executes only at n=12. Larger dense cells are NOT_EXECUTED by design and retain null costs; no failure or hypothetical runtime is inferred.","",
              "Ranked coordinates use O(nD) prefix-count entries, but the ambient basis, multiplier lists and dense row widths remain charged. Whole-worker RSS includes common references and all arms. Fixtures and independent oracle construction are outside arm timing and inside process receipts.","",
              "These are generic Boolean matrix-construction diagnostics. No curve inputs, scalar recovery, relation collection, production solver integration, full-solver timing or rho comparison occurs."]
    if protocol.get("enable_sparse"):
        lines += ["", "This additive run evaluates sparse_rank against the fresh sorted, binary and ranked controls. Only touched nonzero intermediate words are retained; the final dense matrix and all validation costs remain charged. Its fresh protocol and holdouts do not replace the initial failed scaling run."]
    (root/"RESULT.md").write_text("\n".join(lines)+"\n")
    print(json.dumps({k:result[k] for k in ('correctness','cells','paired_batch_samples','verified_outputs','dramatic_scaling_gate','gates_passed')},sort_keys=True))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
