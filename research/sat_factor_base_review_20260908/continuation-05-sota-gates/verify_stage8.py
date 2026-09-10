#!/usr/bin/env python3
"""Fail-closed verifier and summary for the frozen degree-15 panel."""

from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def distribution(values: list[float | int]) -> dict:
    return {
        "min": min(values),
        "median": statistics.median(values),
        "mean": statistics.fmean(values),
        "max": max(values),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    expected = [(17, 20260981), (53, 20260982), (89, 20260983), (137, 20260984), (199, 20260985)]
    rows = []
    for secret, seed in expected:
        ic = load(args.directory / f"ic-{secret}.json")
        ic_meter = load(args.directory / f"ic-{secret}.metrics.json")
        rho = load(args.directory / f"rho-auto-{secret}.json")
        rho_meter = load(args.directory / f"rho-auto-{secret}.metrics.json")

        assert ic_meter["returncode"] == 0 and not ic_meter["timed_out"]
        assert rho_meter["returncode"] == 0 and not rho_meter["timed_out"]
        assert ic["n"] == rho["n"] == 15 and ic["a"] == rho["a"] == 1
        assert ic["seed"] == rho["seed"] == seed
        assert ic["target"] == rho["target"]
        assert ic["factor_base_predicate"]["factor_spec"] == "divisor:2,4"
        assert ic["factor_base_predicate"]["divisor_indices"] == [2, 4]
        assert ic["factor_base_predicate"]["enumerates_target_subgroup"] is False
        assert ic["factor_base_predicate"]["uses_discrete_log_labels"] is False
        assert ic["factor_base_predicate"]["factor_base_logs_constructed"] is False
        assert ic["report"]["collapse_projected_orbits"] is True
        assert ic["report"]["m_cofactor_admissible"] is True
        assert ic["report"]["verified_unknown_scalar_recovery"] is True
        assert ic["report"]["recovered_scalar"] == str(secret)
        assert ic["report"]["direct_relation"] is False
        assert ic["report"]["sat_invalid_models"] == 0
        assert rho["automorphism_optimized"] is True
        assert rho["verified_unknown_scalar_recovery"] is True

        ic_metrics = ic_meter["metrics"]
        rho_metrics = rho_meter["metrics"]
        rows.append(
            {
                "secret": secret,
                "seed": seed,
                "index_calculus": {
                    "relations": ic["report"]["relations"],
                    "trials": ic["report"]["trials"],
                    "conflicts": ic["report"]["sat_conflicts"],
                    "linear_solve_attempts": ic["report"]["linear_solve_attempts"],
                    "linear_algebra_ns": ic["timing_ns"]["linear_algebra"],
                    "wall_seconds": ic_metrics["wall_seconds"],
                    "single_core_seconds": ic_metrics["single_core_seconds"],
                    "total_core_seconds": ic_metrics["total_core_seconds"],
                    "peak_rss_bytes": ic_metrics["peak_rss_bytes"],
                },
                "automorphism_rho": {
                    "iterations": rho["iterations"],
                    "reported_group_additions": rho["reported_group_additions"],
                    "wall_seconds": rho_metrics["wall_seconds"],
                    "single_core_seconds": rho_metrics["single_core_seconds"],
                    "total_core_seconds": rho_metrics["total_core_seconds"],
                    "peak_rss_bytes": rho_metrics["peak_rss_bytes"],
                },
                "ic_over_rho_core_ratio": ic_metrics["total_core_seconds"]
                / rho_metrics["total_core_seconds"],
            }
        )

    summary = {
        "schema": "koblitz_degree15_replication_panel_result.v1",
        "rows": rows,
        "all_checks_pass": True,
        "successes": len(rows),
        "failures": 0,
        "distribution": {
            "index_calculus_relations": distribution([row["index_calculus"]["relations"] for row in rows]),
            "index_calculus_trials": distribution([row["index_calculus"]["trials"] for row in rows]),
            "index_calculus_conflicts": distribution([row["index_calculus"]["conflicts"] for row in rows]),
            "index_calculus_wall_seconds": distribution([row["index_calculus"]["wall_seconds"] for row in rows]),
            "index_calculus_core_seconds": distribution([row["index_calculus"]["total_core_seconds"] for row in rows]),
            "rho_wall_seconds": distribution([row["automorphism_rho"]["wall_seconds"] for row in rows]),
            "rho_core_seconds": distribution([row["automorphism_rho"]["total_core_seconds"] for row in rows]),
            "ic_over_rho_core_ratio": distribution([row["ic_over_rho_core_ratio"] for row in rows]),
        },
        "aggregate": {
            "index_calculus_core_seconds": sum(row["index_calculus"]["total_core_seconds"] for row in rows),
            "rho_core_seconds": sum(row["automorphism_rho"]["total_core_seconds"] for row in rows),
            "maximum_peak_rss_bytes": max(
                max(row["index_calculus"]["peak_rss_bytes"], row["automorphism_rho"]["peak_rss_bytes"])
                for row in rows
            ),
        },
        "claim_boundary": "A five-target toy distribution with rho faster in every pair is not an asymptotic crossover or SOTA result.",
    }
    args.output.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
