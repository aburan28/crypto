#!/usr/bin/env python3
"""Summarize the frozen matched suite without mutating raw evidence."""

from __future__ import annotations

import argparse
import collections
import json
import math
from pathlib import Path
import random
import statistics


def percentile(values, fraction):
    ordered = sorted(values)
    if not ordered:
        return None
    position = fraction * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def paired_interval(ratios, *, seed=20260921, repetitions=4000):
    logs = [math.log(value) for value in ratios]
    rng = random.Random(seed)
    samples = []
    for _ in range(repetitions):
        draw = [logs[rng.randrange(len(logs))] for _ in logs]
        samples.append(math.exp(statistics.mean(draw)))
    return {
        "geometric_mean": math.exp(statistics.mean(logs)),
        "ci95": [percentile(samples, 0.025), percentile(samples, 0.975)],
        "pairs": len(ratios),
        "bootstrap_repetitions": repetitions,
        "seed": seed,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run")
    args = parser.parse_args()
    run = Path(args.run).resolve()
    processes = [json.loads(line) for line in (run / "processes.jsonl").read_text().splitlines()]
    certificates = json.loads((run / "certificates.json").read_text())
    cells = {row["label"]: row for row in json.loads((run / "cells.json").read_text())}
    assert all(row["status"] == "COMPLETED" for row in processes)

    grouped = collections.defaultdict(list)
    for row in processes:
        grouped[(row["cell"], row["target_index"], row["variant"])].append(row)
    target_rows = []
    for (cell, target_index, variant), rows in sorted(grouped.items()):
        certificate = certificates[f"{cell}/target-{target_index}"]
        target_rows.append({
            "cell": cell,
            "target_index": target_index,
            "target": rows[0]["target"],
            "variant": variant,
            "producing": certificate["symmetric_roots"] > 0,
            "roots": certificate["symmetric_roots"],
            "wall_median_seconds": statistics.median(row["wall_seconds"] for row in rows),
            "max_f4_degree": max(row["max_f4_degree"] for row in rows),
            "peak_round_matrix_area": max(row["peak_round_matrix_area"] for row in rows),
            "max_matrix_rows": max(row["max_matrix_rows"] for row in rows),
            "max_matrix_columns": max(row["max_matrix_columns"] for row in rows),
            "system_sha256": rows[0]["system_sha256"],
            "repetitions": len(rows),
        })

    summaries = []
    for cell in cells:
        producing_targets = sorted({row["target_index"] for row in target_rows
                                    if row["cell"] == cell and row["producing"]})
        variants = cells[cell]["variants"]
        direct_by_pair = {(row["target_index"], row["repetition"]): row
                          for row in processes if row["cell"] == cell
                          and row["variant"] == "direct"
                          and row["target_index"] in producing_targets}
        radical_by_pair = {(row["target_index"], row["repetition"]): row
                           for row in processes if row["cell"] == cell
                           and row["variant"] == "radical_symmetric"
                           and row["target_index"] in producing_targets}
        keys = sorted(set(direct_by_pair) & set(radical_by_pair))
        ratios = [radical_by_pair[key]["wall_seconds"] / direct_by_pair[key]["wall_seconds"]
                  for key in keys]
        stage_interval = paired_interval(ratios) if ratios else None
        variant_stats = {}
        for variant in variants:
            selected = [row for row in target_rows if row["cell"] == cell
                        and row["variant"] == variant and row["producing"]]
            variant_stats[variant] = {
                "producing_targets": len(selected),
                "wall_median_seconds": statistics.median(
                    row["wall_median_seconds"] for row in selected) if selected else None,
                "f4_degrees": sorted({row["max_f4_degree"] for row in selected}),
                "peak_matrix_area_median": statistics.median(
                    row["peak_round_matrix_area"] for row in selected) if selected else None,
                "peak_matrix_area_max": max(
                    (row["peak_round_matrix_area"] for row in selected), default=None),
            }
        direct_all = [row for row in target_rows if row["cell"] == cell
                      and row["variant"] == "direct"]
        radical_all = [row for row in target_rows if row["cell"] == cell
                       and row["variant"] == "radical_symmetric"]
        direct_batch = sum(row["wall_median_seconds"] for row in direct_all)
        radical_solver_batch = sum(row["wall_median_seconds"] for row in radical_all)
        preprocessing = cells[cell]["radical_preprocessing"].get(
            "total_seconds",
            cells[cell]["radical_preprocessing"].get("total_preprocessing_seconds", 0.0))
        radical_cold_batch = preprocessing + radical_solver_batch
        summaries.append({
            "cell": cell,
            "factor_base_size": cells[cell]["factor_base_size"],
            "targets": len(cells[cell]["targets"]),
            "producing_targets": len(producing_targets),
            "variants": variant_stats,
            "radical_over_direct_paired_wall": stage_interval,
            "preprocessing_seconds": preprocessing,
            "batch_direct_solver_seconds": direct_batch,
            "batch_radical_solver_seconds": radical_solver_batch,
            "batch_radical_cold_seconds": radical_cold_batch,
            "batch_cold_radical_over_direct": radical_cold_batch / direct_batch,
            "classification": "engineering solver-stage diagnostic",
            "total_common_operations": None,
            "S": None,
            "rho_ratio": None,
            "generic_floor_ratio": None,
        })

    output = {
        "status": "complete",
        "processes": len(processes),
        "cells": len(cells),
        "targets": len(certificates),
        "all_processes_completed": True,
        "all_certificates_pass": all(
            cert["presentation_equivalence"] and cert["nondegenerate_group_equivalence"]
            for cert in certificates.values()),
        "primary_metrics": {
            "total_common_operations": None,
            "S": None,
            "rho_ratio": None,
            "generic_floor_ratio": None,
            "reason": "solver-stage suite only; no complete ECDLP operation conversion",
        },
        "target_results": target_rows,
        "cell_results": summaries,
    }
    (run / "summary.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": output["status"], "processes": len(processes),
                      "targets": len(certificates), "cells": len(cells)}, indent=2))


if __name__ == "__main__":
    main()
