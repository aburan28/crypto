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
    degree15_discovery = json.loads((args.directory / "n15-discovery.json").read_text())
    discovery_meter = json.loads((args.directory / "n15-discovery.metrics.json").read_text())
    assert discovery_meter["returncode"] == 0 and not discovery_meter["timed_out"]
    assert degree15_discovery["selected"]["divisor_indices"] == [2, 4]
    assert degree15_discovery["selected"]["rational_points"] == 281
    assert degree15_discovery["selected"]["projected_signed_frobenius_orbits"] == 6
    assert degree15_discovery["forbidden_inputs"] == {
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

    degree23_discoveries = []
    for curve_a in [0, 1]:
        discovery = json.loads((args.directory / f"n23-discovery-a{curve_a}.json").read_text())
        meter = json.loads((args.directory / f"n23-discovery-a{curve_a}.metrics.json").read_text())
        assert meter["returncode"] == 0 and not meter["timed_out"]
        assert discovery["n"] == 23 and discovery["a"] == curve_a
        assert discovery["requested_dimension"] == 12
        assert discovery["sizing_gate"]["passes"] is True
        assert discovery["forbidden_inputs"] == {
            "target_constructed": False,
            "target_subgroup_enumerated": False,
            "discrete_log_labels_constructed": False,
            "relation_yield_used": False,
            "solver_timing_used": False,
        }
        degree23_discoveries.append(
            {
                "curve_a": curve_a,
                "selected_divisor_indices": discovery["selected"]["divisor_indices"],
                "rational_points": discovery["selected"]["rational_points"],
                "projected_signed_frobenius_orbits": discovery["selected"]["projected_signed_frobenius_orbits"],
                "wall_seconds": meter["metrics"]["wall_seconds"],
                "core_seconds": meter["metrics"]["total_core_seconds"],
                "peak_rss_bytes": meter["metrics"]["peak_rss_bytes"],
            }
        )
    assert degree23_discoveries[0]["selected_divisor_indices"] == [0, 2]
    assert degree23_discoveries[0]["rational_points"] == 4281
    assert degree23_discoveries[0]["projected_signed_frobenius_orbits"] == 93

    degree23_ic = json.loads((args.directory / "n23-ic-101.json").read_text())
    degree23_ic_meter = json.loads((args.directory / "n23-ic-101.metrics.json").read_text())
    degree23_rho = json.loads((args.directory / "n23-rho-auto-101.json").read_text())
    degree23_rho_meter = json.loads((args.directory / "n23-rho-auto-101.metrics.json").read_text())
    assert degree23_ic_meter["returncode"] == 0 and not degree23_ic_meter["timed_out"]
    assert degree23_rho_meter["returncode"] == 0 and not degree23_rho_meter["timed_out"]
    assert degree23_ic["target"] == degree23_rho["target"]
    assert degree23_ic["factor_base_predicate"]["divisor_indices"] == [0, 2]
    assert degree23_ic["factor_base_predicate"]["factor_base_logs_constructed"] is False
    assert degree23_ic["factor_base_predicate"]["uses_discrete_log_labels"] is False
    assert degree23_ic["factor_base_predicate"]["enumerates_target_subgroup"] is False
    assert degree23_ic["report"]["verified_unknown_scalar_recovery"] is True
    assert degree23_ic["report"]["recovered_scalar"] == "101"
    assert degree23_ic["report"]["direct_relation"] is False
    assert degree23_ic["report"]["sat_invalid_models"] == 0
    assert degree23_rho["automorphism_optimized"] is True
    assert degree23_rho["verified_unknown_scalar_recovery"] is True

    result = {
        "schema": "koblitz_external_host_reproduction.v3",
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
        },
        "degree9_rows": degree9_rows,
        "degree15_public_discovery": {
            "selected_divisor_indices": degree15_discovery["selected"]["divisor_indices"],
            "rational_points": degree15_discovery["selected"]["rational_points"],
            "projected_signed_frobenius_orbits": degree15_discovery["selected"]["projected_signed_frobenius_orbits"],
            "wall_seconds": discovery_meter["metrics"]["wall_seconds"],
            "core_seconds": discovery_meter["metrics"]["total_core_seconds"],
            "peak_rss_bytes": discovery_meter["metrics"]["peak_rss_bytes"],
        },
        "degree15_rows": degree15_rows,
        "degree23_public_discoveries": degree23_discoveries,
        "degree23_pair": {
            "index_calculus": {
                "relations": degree23_ic["report"]["relations"],
                "trials": degree23_ic["report"]["trials"],
                "conflicts": degree23_ic["report"]["sat_conflicts"],
                "unknown_sat_targets": degree23_ic["report"]["sat_unknowns"],
                "wall_seconds": degree23_ic_meter["metrics"]["wall_seconds"],
                "core_seconds": degree23_ic_meter["metrics"]["total_core_seconds"],
                "peak_rss_bytes": degree23_ic_meter["metrics"]["peak_rss_bytes"],
            },
            "automorphism_rho": {
                "iterations": degree23_rho["iterations"],
                "restarts": degree23_rho["restarts"],
                "reported_group_additions": degree23_rho["reported_group_additions"],
                "wall_seconds": degree23_rho_meter["metrics"]["wall_seconds"],
                "core_seconds": degree23_rho_meter["metrics"]["total_core_seconds"],
                "peak_rss_bytes": degree23_rho_meter["metrics"]["peak_rss_bytes"],
            },
        },
        "all_checks_pass": True,
        "scope": "GitHub-hosted fresh-checkout reproduction of public factor-base discovery plus degree-9, degree-15 and degree-23 toy correctness and accounting; workflow is project-authored and is not an external novelty review",
    }
    (args.directory / "external-reproduction-summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
