#!/usr/bin/env python3
"""Verify and summarize scalar-blind end-to-end runs and rho controls."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RUN = ROOT / "stage-2-unknown-scalar-20260909"


def load(name: str) -> tuple[dict, dict]:
    result = json.loads((RUN / f"{name}.json").read_text())
    metrics = json.loads((RUN / f"{name}.metrics.json").read_text())
    return result, metrics


def main() -> None:
    rows = []
    for secret in [53, 101]:
        ic, ic_meter = load(f"ic-secret-{secret}")
        rho, rho_meter = load(f"rho-secret-{secret}")
        assert ic_meter["returncode"] == 0 and not ic_meter["timed_out"]
        assert rho_meter["returncode"] == 0 and not rho_meter["timed_out"]
        assert ic["target"] == rho["target"]
        assert ic["factor_base_predicate"]["factor_base_logs_constructed"] is False
        assert ic["factor_base_predicate"]["uses_discrete_log_labels"] is False
        assert ic["factor_base_predicate"]["enumerates_target_subgroup"] is False
        assert ic["report"]["direct_relation"] is False
        assert ic["report"]["sat_invalid_models"] == 0
        assert ic["report"]["verified_unknown_scalar_recovery"] is True
        assert rho["verified_unknown_scalar_recovery"] is True
        assert rho["automorphism_optimized"] is False
        rows.append(
            {
                "secret": secret,
                "target": ic["target"],
                "index_calculus": {
                    "relations": ic["report"]["relations"],
                    "trials": ic["report"]["trials"],
                    "sat_conflicts": ic["report"]["sat_conflicts"],
                    "linear_solve_attempts": ic["report"]["linear_solve_attempts"],
                    "wall_seconds": ic_meter["metrics"]["wall_seconds"],
                    "total_core_seconds": ic_meter["metrics"]["total_core_seconds"],
                    "peak_rss_bytes": ic_meter["metrics"]["peak_rss_bytes"],
                    "internal_end_to_end_seconds": ic["timing_ns"]["end_to_end"] / 1e9,
                    "verified": True,
                },
                "generic_rho": {
                    "iterations": rho["iterations"],
                    "wall_seconds": rho_meter["metrics"]["wall_seconds"],
                    "total_core_seconds": rho_meter["metrics"]["total_core_seconds"],
                    "peak_rss_bytes": rho_meter["metrics"]["peak_rss_bytes"],
                    "internal_end_to_end_seconds": rho["timing_ns"]["end_to_end"] / 1e9,
                    "automorphism_optimized": False,
                    "verified": True,
                },
            }
        )
    summary = {
        "schema": "koblitz_unknown_scalar_stage2_summary.v1",
        "curve": {"n": 9, "a": 0, "subgroup_order": 127},
        "factor_base": {
            "kind": "GGMP linearised-polynomial kernel",
            "ell": 6,
            "points": 55,
            "signed_frobenius_orbits": 4,
            "factor_base_logs_constructed": False,
            "target_subgroup_enumerated_for_factor_base_selection": False,
        },
        "rows": rows,
        "checks": {
            "both_unknown_scalars_recovered_and_verified": True,
            "no_direct_relation_shortcut": True,
            "no_invalid_sat_model": True,
            "same_targets_used_for_ic_and_rho": True,
        },
        "claim": "two complete toy unknown-scalar runs over an algebraically defined factor base; generic rho controls only; no automorphism-optimized comparison or SOTA claim",
    }
    (ROOT / "stage-2-summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    lines = [
        "# Stage 2: unknown-scalar runs over an algebraic factor base",
        "",
        "Both runs use the GGMP linearised-polynomial kernel selected from public field parameters. The implementation constructs no factor-base logarithms and no subgroup log table. Each recovered scalar is accepted only after recomputing the public target point.",
        "",
        "| secret | IC relations/trials | SAT conflicts | IC wall / core-s / MiB | generic rho iterations | rho wall / core-s / MiB |",
        "|--:|--:|--:|--:|--:|--:|",
    ]
    for row in rows:
        ic, rho = row["index_calculus"], row["generic_rho"]
        lines.append(
            f"| {row['secret']} | {ic['relations']}/{ic['trials']} | {ic['sat_conflicts']} | "
            f"{ic['wall_seconds']:.6f} / {ic['total_core_seconds']:.6f} / {ic['peak_rss_bytes']/1048576:.2f} | "
            f"{rho['iterations']} | {rho['wall_seconds']:.6f} / {rho['total_core_seconds']:.6f} / {rho['peak_rss_bytes']/1048576:.2f} |"
        )
    lines += [
        "",
        "Each index-calculus run collected five relations in five trials, made one modular linear-algebra attempt, avoided the direct-relation shortcut, and had zero invalid SAT models. The two runs used 122 and 125 SAT conflicts.",
        "",
        "The first IC process has a large cold-start wall outlier relative to its 0.004329 core-seconds and 0.002152-second internal total. With only two runs, no stable timing ratio is inferred. Both generic-rho controls are faster on internal time and total core-seconds.",
        "",
        "These rho controls do not quotient by signed Frobenius orbits. They therefore do not satisfy the requested automorphism-optimized Pollard-rho gate. This stage establishes scalar-blind toy completion, not a competitive result.",
        "",
    ]
    (ROOT / "STAGE2_RESULTS.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
