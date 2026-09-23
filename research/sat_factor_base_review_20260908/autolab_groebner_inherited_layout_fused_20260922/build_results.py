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
POLICIES = ("fused", "materialized", "uncached")


def load(path: Path):
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def median(values):
    return statistics.median(values)


def json_lines(path: Path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.startswith("{")]


def fused_kernel():
    rows = json_lines(RUNS / "fused_root_kernel" / "stdout.txt")
    groups = defaultdict(list)
    for row in rows:
        assert row["correct"] is True
        groups[row["case"]].append(row)
    cells = []
    for case, samples in sorted(groups.items()):
        assert len({(row["rows"], row["cols"], row["fixture_hash"]) for row in samples}) == 1
        materialized = median(row["materialized_ns"] for row in samples)
        fused = median(row["fused_ns"] for row in samples)
        cells.append({
            "case": case,
            "samples": len(samples),
            "rows": samples[0]["rows"],
            "cols": samples[0]["cols"],
            "fixture_hash": samples[0]["fixture_hash"],
            "materialized_median_ns": materialized,
            "fused_median_ns": fused,
            "median_speedup": materialized / fused,
            "paired_ratio_median": median(
                row["materialized_ns"] / row["fused_ns"] for row in samples
            ),
        })
    return {
        "receipt": load(RUNS / "fused_root_kernel" / "receipt.json"),
        "samples": len(rows),
        "cells": cells,
        "aggregate_paired_ratio_median": median(
            row["materialized_ns"] / row["fused_ns"] for row in rows
        ),
    }


EXACT_FIELDS = [
    "curve", "a", "n", "m", "ell", "targets", "decomposed",
    "verdict_digest", "word_ops", "f4_calls", "f4_oversize",
    "reductions", "infeasible_branches", "propagations", "splits",
    "matrix_rows", "matrix_cols", "rows_pruned", "criterion_word_ops",
    "specialise_word_ops", "exhausted", "oversize",
]


def stage_panel(group: str, raw: dict):
    observations = {policy: [] for policy in POLICIES}
    rung_samples = defaultdict(lambda: {policy: [] for policy in POLICIES})
    for policy in POLICIES:
        for rep in range(10):
            run = raw["runs"][f"{policy}_{group}_r{rep}"]
            receipt, stage = run["receipt"], run["stage"]
            observation = {
                "process_wall_ms": receipt["wall_ms"],
                "child_cpu_ms": receipt["child_user_ms"] + receipt["child_system_ms"],
                "peak_rss_bytes": receipt["peak_rss_bytes"],
                "stage_wall_ms": sum(row["wall_ns"] for row in stage["rows"]) / 1e6,
                "build_ms": sum(row["build_ns"] for row in stage["rows"]) / 1e6,
                "reduce_ms": sum(row["reduce_ns"] for row in stage["rows"]) / 1e6,
            }
            observations[policy].append(observation)
            for row in stage["rows"]:
                rung_samples[(row["curve"], row["m"], row["ell"])][policy].append(row)

    for rep in range(10):
        selected = raw["runs"][f"fused_{group}_r{rep}"]["stage"]["rows"]
        for policy in ("materialized", "uncached"):
            control = raw["runs"][f"{policy}_{group}_r{rep}"]["stage"]["rows"]
            assert len(selected) == len(control)
            for a, b in zip(selected, control):
                for field in EXACT_FIELDS:
                    assert a.get(field) == b.get(field), (group, rep, policy, field)

    metrics = tuple(observations["fused"][0])
    medians = {
        policy: {
            metric: median(row[metric] for row in observations[policy])
            for metric in metrics
        }
        for policy in POLICIES
    }
    speedups = {
        control: {
            metric: medians[control][metric] / medians["fused"][metric]
            for metric in metrics
        }
        for control in ("materialized", "uncached")
    }
    rungs = []
    for (curve, m, ell), policies in rung_samples.items():
        row = {
            "curve": curve,
            "m": m,
            "ell": ell,
            "layout_hits_median": median(sample["layout_hits"] for sample in policies["fused"]),
            "layout_misses_median": median(
                sample["layout_misses"] for sample in policies["fused"]
            ),
            "medians": {},
        }
        for policy in POLICIES:
            row["medians"][policy] = {
                "stage_wall_ms": median(sample["wall_ns"] for sample in policies[policy]) / 1e6,
                "build_ms": median(sample["build_ns"] for sample in policies[policy]) / 1e6,
            }
        row["speedups"] = {
            control: {
                metric: row["medians"][control][metric] / row["medians"]["fused"][metric]
                for metric in ("stage_wall_ms", "build_ms")
            }
            for control in ("materialized", "uncached")
        }
        if m == 2:
            assert row["speedups"]["uncached"]["stage_wall_ms"] > 1
        else:
            assert row["layout_hits_median"] == 0 and row["layout_misses_median"] == 0
        rungs.append(row)
    return {
        "repetitions_per_policy": 10,
        "medians": medians,
        "speedups": speedups,
        "rungs": rungs,
        "exact_scientific_fields": EXACT_FIELDS,
        "exact_gate": "PASS",
    }


protocol = load(HERE / "protocol.json")
raw = load(HERE / "stage_panel_raw.json")
kernel = fused_kernel()
frozen = stage_panel("frozen", raw)
holdout = stage_panel("holdout", raw)
exact_receipt = load(RUNS / "exact_tests" / "receipt.json")
assert exact_receipt["exit_code"] == 0
assert all(cell["median_speedup"] > 1 for cell in kernel["cells"])

result = {
    "schema_version": "1.0",
    "task_id": protocol["task_id"],
    "classification": "INHERITED_F4_QUADRATIC_LAYOUT_FUSED_ROOT_AND_STAGE_PASS",
    "selected_policy": protocol["selected_policy"],
    "exact_tests": {"status": "PASS", "receipt": exact_receipt},
    "fused_root_kernel": kernel,
    "stage_panels": {"frozen": frozen, "holdout": holdout},
    "interpretation": {
        "root_packing": "PASS: every real Semaev root cell improved with byte-identical matrices",
        "quadratic_stage": "PASS: every quadratic frozen and holdout rung improved against uncached roots",
        "cubic_policy": "BYPASS: the automatic policy leaves the unstable cubic support uncached",
        "scientific_promotion": "engineering only",
    },
    "claim_boundary": protocol["claim_boundary"],
}
(HERE / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

largest = next(cell for cell in kernel["cells"] if cell["case"].startswith("binary-n31"))
source_hashes = {path.name: sha256(path) for path in sorted(FROZEN.iterdir()) if path.is_file()}
stage_hashes = {run["receipt"]["executable_sha256"] for run in raw["runs"].values()}
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
    "oracle_class": "inherited matrix-F4 with exact quadratic root-layout reuse and fused nested packing",
    "median_ms_per_target": {"n31_fused_root_pack_ms": largest["fused_median_ns"] / 1e6},
    "largest_solvable": {"n": 31, "m": 2, "unknowns": 32, "root_rows": largest["rows"]},
    "fixture_hash": largest["fixture_hash"],
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
        "root_repetitions": 100,
        "stage_repetitions_per_policy": 10,
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
