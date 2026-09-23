#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import platform
import socket
import statistics
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
FROZEN = HERE / "frozen"


def load(path: Path):
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_lines(path: Path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.startswith("{")]


def median(values):
    return statistics.median(values)


def kernel_summary(run: str):
    rows = json_lines(RUNS / run / "stdout.txt")
    groups = defaultdict(list)
    for row in rows:
        assert row["correct"] is True
        groups[(row["case"], row["reduce_above"])].append(row)
    cells = []
    for (case, reduce_above), samples in sorted(groups.items()):
        assert len({(row["rank"], row["word_ops"], row["input_hash"]) for row in samples}) == 1
        allocating_ns = median(row["allocating_ns"] for row in samples)
        arena_ns = median(row["arena_ns"] for row in samples)
        cells.append({
            "case": case,
            "reduce_above": reduce_above,
            "samples": len(samples),
            "rank": samples[0]["rank"],
            "word_ops": samples[0]["word_ops"],
            "input_hash": samples[0]["input_hash"],
            "allocating_median_ns": allocating_ns,
            "arena_median_ns": arena_ns,
            "median_speedup": allocating_ns / arena_ns,
            "paired_ratio_median": median(
                row["allocating_ns"] / row["arena_ns"] for row in samples
            ),
        })
    return {
        "receipt": load(RUNS / run / "receipt.json"),
        "samples": len(rows),
        "cells": cells,
        "aggregate_paired_ratio_median": median(
            row["allocating_ns"] / row["arena_ns"] for row in rows
        ),
    }


EXACT_STAGE_FIELDS = [
    "curve", "a", "n", "m", "ell", "targets", "decomposed",
    "verdict_digest", "word_ops", "f4_calls", "f4_oversize",
    "reductions", "infeasible_branches", "propagations", "splits",
    "matrix_rows", "matrix_cols", "rows_pruned", "criterion_word_ops",
    "specialise_word_ops", "exhausted", "oversize",
]


def stage_summary(group: str, repetitions: int, raw: dict):
    observations = {"arena": [], "allocating": []}
    for policy in observations:
        for rep in range(repetitions):
            run = raw["runs"][f"{policy}_{group}_r{rep}"]
            receipt, stage = run["receipt"], run["stage"]
            observations[policy].append({
                "process_wall_ms": receipt["wall_ms"],
                "child_cpu_ms": receipt["child_user_ms"] + receipt["child_system_ms"],
                "peak_rss_bytes": receipt["peak_rss_bytes"],
                "stage_wall_ms": sum(row["wall_ns"] for row in stage["rows"]) / 1e6,
                "build_ms": sum(row["build_ns"] for row in stage["rows"]) / 1e6,
                "reduce_ms": sum(row["reduce_ns"] for row in stage["rows"]) / 1e6,
                "readback_ms": sum(row["readback_ns"] for row in stage["rows"]) / 1e6,
                "receipt": receipt,
                "stage": stage,
            })
    for rep in range(repetitions):
        arena = observations["arena"][rep]["stage"]["rows"]
        allocating = observations["allocating"][rep]["stage"]["rows"]
        assert len(arena) == len(allocating)
        for selected, control in zip(arena, allocating):
            for field in EXACT_STAGE_FIELDS:
                assert selected.get(field) == control.get(field), (group, rep, field)

    metrics = [
        "process_wall_ms", "child_cpu_ms", "peak_rss_bytes",
        "stage_wall_ms", "build_ms", "reduce_ms", "readback_ms",
    ]
    medians = {
        policy: {
            metric: median(row[metric] for row in observations[policy])
            for metric in metrics
        }
        for policy in observations
    }
    def ratio(control, selected):
        if selected == 0:
            return 1.0 if control == 0 else None
        return control / selected

    speedups = {
        metric: ratio(medians["allocating"][metric], medians["arena"][metric])
        for metric in metrics
        if metric != "peak_rss_bytes"
    }
    speedups["allocating_over_arena_peak_rss"] = (
        medians["allocating"]["peak_rss_bytes"] / medians["arena"]["peak_rss_bytes"]
    )
    return {
        "repetitions_per_policy": repetitions,
        "medians": medians,
        "speedups": speedups,
        "exact_scientific_fields": EXACT_STAGE_FIELDS,
        "exact_gate": "PASS",
    }


protocol = load(HERE / "protocol.json")
stage_raw = load(HERE / "stage_panel_raw.json")
synthetic = kernel_summary("kernel_synthetic")
semaev = kernel_summary("kernel_semaev")
frozen = stage_summary("frozen", 13, stage_raw)
holdout = stage_summary("holdout", 15, stage_raw)
exact_receipt = load(RUNS / "exact_equivalence" / "receipt.json")
assert exact_receipt["exit_code"] == 0
assert all(cell["median_speedup"] > 1 for cell in synthetic["cells"])
assert all(cell["median_speedup"] > 1 for cell in semaev["cells"])

result = {
    "schema_version": "1.0",
    "task_id": protocol["task_id"],
    "classification": "INHERITED_F4_M4RI_ARENA_ROOT_KERNEL_PASS_STAGE_TIMING_INDETERMINATE",
    "selected_policy": protocol["selected_policy"],
    "exact_equivalence": {
        "status": "PASS",
        "matrix_cases": 96,
        "receipt": exact_receipt,
    },
    "synthetic_kernel": synthetic,
    "semaev_root_kernel": semaev,
    "stage_panels": {"frozen": frozen, "holdout": holdout},
    "interpretation": {
        "root_kernel": "PASS: every retained cell improved with identical algebraic output and counted work",
        "whole_stage": "INDETERMINATE: inherited F4 leaves root elimination as a small fraction of total cost and host-sensitive CPU and wall metrics disagree",
        "scientific_promotion": "engineering only",
    },
    "claim_boundary": protocol["claim_boundary"],
}
(HERE / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

largest = next(cell for cell in semaev["cells"] if cell["case"].startswith("binary-n31"))
source_hashes = {
    path.name: sha256(path) for path in sorted(FROZEN.iterdir()) if path.is_file()
}
stage_hashes = {
    run["receipt"]["executable_sha256"]
    for run in stage_raw["runs"].values()
}
assert len(stage_hashes) == 1
claim = {
    "schema_version": 2,
    "task_id": protocol["task_id"],
    "stage": "decomposition",
    "n_or_bits": 31,
    "m_summands": 2,
    "dimension_l_or_dim": 16,
    "unknowns": 32,
    "system_degree": "quadratic generators, degree-3 Macaulay root",
    "eq_var_ratio": 31 / 32,
    "ffd_or_degree_of_regularity": 3,
    "oracle_class": "inherited matrix-F4 with reusable block-4 M4RI scratch",
    "median_ms_per_target": {"n31_root_arena_ms": largest["arena_median_ns"] / 1e6},
    "largest_solvable": {"n": 31, "m": 2, "unknowns": 32, "root_rank": largest["rank"]},
    "fixture_hash": largest["input_hash"],
    "executable_or_source_hash": {
        "sources": source_hashes,
        "stage_executable": next(iter(stage_hashes)),
    },
    "host_id": {
        "hostname": socket.gethostname(),
        "system": platform.system(),
        "machine": platform.machine(),
    },
    "resource_caps": {
        "threads": 1,
        "block_width": 4,
        "synthetic_repetitions": 100,
        "semaev_root_repetitions": 100,
        "frozen_stage_repetitions_per_policy": 13,
        "holdout_stage_repetitions_per_policy": 15,
    },
    "seeds": [17, 937, 20260914, 66142],
    "claim_boundary": protocol["claim_boundary"],
    "claim_boundary_non_claims": [
        "not relation-yield evidence",
        "not an end-to-end index-calculus crossover",
        "not an asymptotic or sub-rho result",
        "not external or private targets",
        "not a key-recovery claim"
    ],
}
(HERE / "claim_report_decomposition.json").write_text(
    json.dumps(claim, indent=2, sort_keys=True) + "\n"
)
required_global = [
    "fixture_hash", "executable_or_source_hash", "host_id", "resource_caps",
    "seeds", "claim_boundary_non_claims",
]
required_stage = [
    "n_or_bits", "m_summands", "dimension_l_or_dim", "unknowns",
    "system_degree", "eq_var_ratio", "ffd_or_degree_of_regularity",
    "oracle_class", "median_ms_per_target", "largest_solvable",
]
missing_global = [field for field in required_global if field not in claim]
missing_stage = [field for field in required_stage if field not in claim]
(HERE / "claim_check_decomposition.json").write_text(json.dumps({
    "schema_version": 2,
    "stage": "decomposition",
    "fail_closed": True,
    "required_global_provenance": required_global,
    "required_stage_fields": required_stage,
    "missing_global_provenance": missing_global,
    "missing_stage_fields": missing_stage,
    "status": "PASS" if not missing_global and not missing_stage else "SCHEMA_INCOMPLETE",
}, indent=2, sort_keys=True) + "\n")

files = {
    str(path.relative_to(HERE)): sha256(path)
    for path in sorted(HERE.rglob("*"))
    if path.is_file() and path.name != "manifest.json"
}
(HERE / "manifest.json").write_text(json.dumps({
    "schema_version": "1.0",
    "task_id": protocol["task_id"],
    "files": files,
}, indent=2, sort_keys=True) + "\n")
