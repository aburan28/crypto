#!/usr/bin/env python3
"""Fail-closed verifier for GitHub-hosted Koblitz toy reproduction artifacts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    degree9_rows = []
    for secret in [53, 101]:
        ic = json.loads((args.directory / f"ic-{secret}.json").read_text())
        ic_meter = json.loads((args.directory / f"ic-{secret}.metrics.json").read_text())
        rho = json.loads((args.directory / f"rho-auto-{secret}.json").read_text())
        rho_meter = json.loads((args.directory / f"rho-auto-{secret}.metrics.json").read_text())
        assert ic_meter["returncode"] == 0 and not ic_meter["timed_out"]
        assert rho_meter["returncode"] == 0 and not rho_meter["timed_out"]
        assert ic["target"] == rho["target"]
        assert ic["report"]["verified_unknown_scalar_recovery"] is True
        assert ic["report"]["recovered_scalar"] == str(secret)
        assert ic["report"]["direct_relation"] is False
        assert ic["report"]["sat_invalid_models"] == 0
        assert ic["factor_base_predicate"]["factor_base_logs_constructed"] is False
        assert ic["factor_base_predicate"]["uses_discrete_log_labels"] is False
        assert ic["factor_base_predicate"]["enumerates_target_subgroup"] is False
        assert rho["automorphism_optimized"] is True
        assert rho["verified_unknown_scalar_recovery"] is True
        degree9_rows.append(
            {
                "secret": secret,
                "index_calculus": {
                    "relations": ic["report"]["relations"],
                    "trials": ic["report"]["trials"],
                    "conflicts": ic["report"]["sat_conflicts"],
                    "wall_seconds": ic_meter["metrics"]["wall_seconds"],
                    "core_seconds": ic_meter["metrics"]["total_core_seconds"],
                    "peak_rss_bytes": ic_meter["metrics"]["peak_rss_bytes"],
                },
                "automorphism_rho": {
                    "iterations": rho["iterations"],
                    "reported_group_additions": rho["reported_group_additions"],
                    "wall_seconds": rho_meter["metrics"]["wall_seconds"],
                    "core_seconds": rho_meter["metrics"]["total_core_seconds"],
                    "peak_rss_bytes": rho_meter["metrics"]["peak_rss_bytes"],
                },
            }
        )
    discovery = json.loads((args.directory / "n15-discovery.json").read_text())
    discovery_meter = json.loads((args.directory / "n15-discovery.metrics.json").read_text())
    assert discovery_meter["returncode"] == 0 and not discovery_meter["timed_out"]
    assert discovery["selected"]["divisor_indices"] == [2, 4]
    assert discovery["selected"]["rational_points"] == 281
    assert discovery["selected"]["projected_signed_frobenius_orbits"] == 6
    assert discovery["forbidden_inputs"] == {
        "target_constructed": False,
        "target_subgroup_enumerated": False,
        "discrete_log_labels_constructed": False,
        "relation_yield_used": False,
        "solver_timing_used": False,
    }

    degree15_rows = []
    for secret, seed in [(17, 20260981), (53, 20260982), (89, 20260983), (137, 20260984), (199, 20260985)]:
        ic = json.loads((args.directory / f"n15-ic-{secret}.json").read_text())
        ic_meter = json.loads((args.directory / f"n15-ic-{secret}.metrics.json").read_text())
        rho = json.loads((args.directory / f"n15-rho-auto-{secret}.json").read_text())
        rho_meter = json.loads((args.directory / f"n15-rho-auto-{secret}.metrics.json").read_text())
        assert ic_meter["returncode"] == 0 and not ic_meter["timed_out"]
        assert rho_meter["returncode"] == 0 and not rho_meter["timed_out"]
        assert ic["seed"] == rho["seed"] == seed and ic["target"] == rho["target"]
        assert ic["factor_base_predicate"]["divisor_indices"] == [2, 4]
        assert ic["factor_base_predicate"]["factor_base_logs_constructed"] is False
        assert ic["factor_base_predicate"]["uses_discrete_log_labels"] is False
        assert ic["factor_base_predicate"]["enumerates_target_subgroup"] is False
        assert ic["report"]["collapse_projected_orbits"] is True
        assert ic["report"]["m_cofactor_admissible"] is True
        assert ic["report"]["verified_unknown_scalar_recovery"] is True
        assert ic["report"]["recovered_scalar"] == str(secret)
        assert ic["report"]["direct_relation"] is False
        assert ic["report"]["sat_invalid_models"] == 0
        assert rho["automorphism_optimized"] is True
        assert rho["verified_unknown_scalar_recovery"] is True
        degree15_rows.append(
            {
                "secret": secret,
                "seed": seed,
                "index_calculus": {
                    "relations": ic["report"]["relations"],
                    "trials": ic["report"]["trials"],
                    "conflicts": ic["report"]["sat_conflicts"],
                    "wall_seconds": ic_meter["metrics"]["wall_seconds"],
                    "core_seconds": ic_meter["metrics"]["total_core_seconds"],
                    "peak_rss_bytes": ic_meter["metrics"]["peak_rss_bytes"],
                },
                "automorphism_rho": {
                    "iterations": rho["iterations"],
                    "reported_group_additions": rho["reported_group_additions"],
                    "wall_seconds": rho_meter["metrics"]["wall_seconds"],
                    "core_seconds": rho_meter["metrics"]["total_core_seconds"],
                    "peak_rss_bytes": rho_meter["metrics"]["peak_rss_bytes"],
                },
            }
        )

    result = {
        "schema": "koblitz_external_host_reproduction.v2",
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
        },
        "degree9_rows": degree9_rows,
        "degree15_public_discovery": {
            "selected_divisor_indices": discovery["selected"]["divisor_indices"],
            "rational_points": discovery["selected"]["rational_points"],
            "projected_signed_frobenius_orbits": discovery["selected"]["projected_signed_frobenius_orbits"],
            "wall_seconds": discovery_meter["metrics"]["wall_seconds"],
            "core_seconds": discovery_meter["metrics"]["total_core_seconds"],
            "peak_rss_bytes": discovery_meter["metrics"]["peak_rss_bytes"],
        },
        "degree15_rows": degree15_rows,
        "all_checks_pass": True,
        "scope": "GitHub-hosted fresh-checkout reproduction of public factor-base discovery plus degree-9 and degree-15 toy correctness and accounting; workflow is project-authored and is not an external novelty review",
    }
    (args.directory / "external-reproduction-summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
