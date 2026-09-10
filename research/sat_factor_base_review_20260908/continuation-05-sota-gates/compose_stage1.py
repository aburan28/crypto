#!/usr/bin/env python3
"""Compose v2 CryptoMiniSat data and v3 clean-build WDSat data."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
V2 = ROOT / "stage-1-seed-20260909-v2" / "result.json"
V3 = ROOT / "stage-1-wdsat-rebuild-20260909-v3" / "result.json"


def keyed(report: dict) -> dict:
    return {
        (row["cell"]["n"], row["cell"]["ell"], row["cell"]["m"], row["cell"]["basis"]): row
        for row in report["instances"]
    }


def solver(row: dict, name: str) -> dict:
    return next(item for item in row.get("solvers", []) if item["solver"] == name)


def metric(record: dict, name: str):
    return record.get("metrics", {}).get(name)


def main() -> None:
    v2 = json.loads(V2.read_text())
    v3 = json.loads(V3.read_text())
    a, b = keyed(v2), keyed(v3)
    assert a.keys() == b.keys()
    rows = []
    for key in a:
        old, repaired = a[key], b[key]
        ma, mb = old.get("manifest"), repaired.get("manifest")
        if ma is None or mb is None:
            rows.append(
                {
                    "cell": old["cell"],
                    "status": "failed_operational",
                    "reason": old.get("artifact_status", repaired.get("artifact_status")),
                }
            )
            continue
        assert ma["exports"] == mb["exports"], f"source export drift at {key}"
        assert ma["factor_base_basis_bitmasks"] == mb["factor_base_basis_bitmasks"]
        native = ma["native_sat"]
        mitm = ma["direct_meet_in_the_middle"]
        cms = solver(old, "cryptominisat")
        wdsat = solver(repaired, "wdsat")
        build = repaired["wdsat_build"]
        rows.append(
            {
                "cell": old["cell"],
                "status": "completed",
                "representation": ma["representation"],
                "factor_base_predicate": ma["factor_base_predicate"],
                "source": {
                    "variables": ma["source_variables"],
                    "equations": ma["source_equations"],
                    "max_degree": ma["source_max_degree"],
                    "anf_blake3": ma["exports"]["wdsat_anf"]["blake3"],
                    "xor_dimacs_blake3": ma["exports"]["cryptominisat_xor_dimacs"]["blake3"],
                    "magma_blake3": ma["exports"]["magma_boolean_f4"]["blake3"],
                    "v2_v3_exports_identical": True,
                },
                "combined_exporter_process": {
                    "wall_seconds": metric(old["generator"], "wall_seconds"),
                    "total_core_seconds": metric(old["generator"], "total_core_seconds"),
                    "peak_rss_bytes": metric(old["generator"], "peak_rss_bytes"),
                    "includes": "predicate, target, source/export, native SAT, direct MITM",
                },
                "native_sat": {
                    "status": native["result"],
                    "source_model_valid": native["source_model_valid"],
                    "conflicts": native["stats"]["conflicts"],
                    "wall_seconds": ma["timing_ns"]["native_encoding_and_solve"] / 1e9,
                    "process_core_and_peak_note": "available only for the combined exporter process in stage 1",
                },
                "direct_mitm": {
                    "status": mitm["status"],
                    "wall_seconds": mitm["wall_ns"] / 1e9,
                    "factor_points": mitm["factor_points"],
                    "group_additions": mitm["group_additions"],
                },
                "wdsat": {
                    "status": wdsat["status"],
                    "source_model_valid": wdsat.get("source_model_valid"),
                    "conflicts": wdsat.get("conflicts"),
                    "wall_seconds": metric(wdsat, "wall_seconds"),
                    "total_core_seconds": metric(wdsat, "total_core_seconds"),
                    "peak_rss_bytes": metric(wdsat, "peak_rss_bytes"),
                    "build_total_core_seconds": metric(build, "total_core_seconds"),
                    "build_peak_rss_bytes": metric(build, "peak_rss_bytes"),
                    "binary_sha256": build["binary_sha256"],
                    "rational_point_decomposition": (
                        False if old["cell"]["basis"] == "ggmp" else (True if wdsat["status"] == "sat" else None)
                    ),
                },
                "cryptominisat": {
                    "status": cms["status"],
                    "source_model_valid": cms.get("source_model_valid"),
                    "conflicts": cms.get("conflicts"),
                    "wall_seconds": metric(cms, "wall_seconds"),
                    "total_core_seconds": metric(cms, "total_core_seconds"),
                    "peak_rss_bytes": metric(cms, "peak_rss_bytes"),
                    "rational_point_decomposition": (
                        False if old["cell"]["basis"] == "ggmp" else (True if cms["status"] == "sat" else None)
                    ),
                },
                "magma_f4": {"status": "unavailable_operational"},
                "independent_geometry_qualification": (
                    {
                        "status": "degenerate_algebraic_model_check",
                        "factor_curve_points": 1,
                        "only_point": "T=(0,1)",
                        "planted_relation": "T+T+T=T",
                        "external_model_x_coordinates": [0, 65536, 65536],
                        "nonlifting_coordinates": [65536, 65536],
                        "meaningful_point_decomposition_model": False,
                    }
                    if old["cell"]["basis"] == "ggmp"
                    else None
                ),
            }
        )

    summary = {
        "schema": "koblitz_pdp_stage1_composition.v1",
        "source_runs": {
            "v2": str(V2.relative_to(ROOT)),
            "v3_wdsat_replacement": str(V3.relative_to(ROOT)),
        },
        "composition_rule": "v3 replaces only v2 WDSat cells; v2 supplies native SAT, direct MITM and CryptoMiniSat",
        "rows": rows,
        "checks": {
            "all_available_sat_models_source_valid": all(
                arm.get("status") != "sat" or arm.get("source_model_valid") is True
                for row in rows
                if row["status"] == "completed"
                for arm in [row["native_sat"], row["wdsat"], row["cryptominisat"]]
            ),
            "all_v2_v3_exports_identical": all(
                row.get("source", {}).get("v2_v3_exports_identical", False)
                for row in rows
                if row["status"] == "completed"
            ),
            "unknown_never_counted_as_unsat": True,
            "factor_base_predicates_deny_scalar_labels": all(
                not row["factor_base_predicate"]["uses_discrete_log_labels"]
                and not row["factor_base_predicate"]["enumerates_target_subgroup"]
                for row in rows
                if row["status"] == "completed"
            ),
            "all_sat_models_are_verified_point_decompositions": False,
            "point_model_shortfall": "external SAT models were checked against source ANF; the n=31 GGMP models independently fail curve lifting",
        },
        "gate_status": {
            "1_full_cost_all_stages": "partial: stage-1 PDP setup/build/solve resources are charged; relation collection and linear algebra await the end-to-end stage",
            "2_matched_solver_matrix": "partial: native XOR, WDSat, CryptoMiniSat and direct MITM ran; Magma F4 is unavailable; the n=31 GGMP row is a degenerate algebraic-model check",
            "3_required_metrics": "partial: external solver process metrics are complete; native SAT and MITM need isolated process receipts instead of the combined exporter receipt",
            "4_n31_n41_and_larger_pdp": "partial: planted targets are point-decomposable by direct MITM at n=31, n=41 and n=59, but external SAT outputs are only source-ANF validated and the GGMP n=31 model is nonlifting",
            "5_unknown_scalar_end_to_end": "pending",
            "6_automorphism_optimized_pollard_rho": "pending",
            "7_external_reproduction_and_novelty": "pending",
        },
        "claim": "Strong internal engineering and planted-PDP evidence only; not a Koblitz index-calculus SOTA result",
    }
    (ROOT / "stage-1-composed-summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    completed = [row for row in rows if row["status"] == "completed"]
    lines = [
        "# Stage 1: matched planted-PDP solver matrix",
        "",
        "This stage defines every factor base from public field algebra. No run enumerates the target subgroup or constructs factor-base discrete-log labels. Each cell exports one source system to WDSat ANF, CryptoMiniSat CNF-XOR, and Magma Boolean F4; native SAT and direct meet-in-the-middle use the same target and factor-base predicate.",
        "",
        "The table reports solver-kernel wall time. WDSat and CryptoMiniSat also have isolated total core-seconds and peak RSS in `stage-1-composed-summary.json`. Native SAT and MITM have exact internal wall clocks, while their process CPU/RSS are combined in the exporter receipt; that accounting split is a remaining gate, not silently estimated here.",
        "",
        "| n | ell | base | vars/eqs | native result, conflicts, s | WDSat result, conflicts, s | CMS result, conflicts, s | MITM result, s |",
        "|--:|--:|:--|--:|:--|:--|:--|:--|",
    ]
    for row in completed:
        cell, src = row["cell"], row["source"]
        native, wdsat, cms, mitm = row["native_sat"], row["wdsat"], row["cryptominisat"], row["direct_mitm"]
        def arm(value: dict) -> str:
            conflicts = "-" if value.get("conflicts") is None else str(value["conflicts"])
            wall = "-" if value.get("wall_seconds") is None else f"{value['wall_seconds']:.6f}"
            return f"{value['status']}, {conflicts}, {wall}"
        lines.append(
            f"| {cell['n']} | {cell['ell']} | {cell['basis']} | {src['variables']}/{src['equations']} | "
            f"{arm(native)} | {arm(wdsat)} | {arm(cms)} | {mitm['status']}, {mitm['wall_seconds']:.6f} |"
        )
    lines += [
        "",
        "| n | base | combined exporter core-s / peak MiB | WDSat core-s / peak MiB | WDSat build core-s / peak MiB | CMS core-s / peak MiB |",
        "|--:|:--|--:|--:|--:|--:|",
    ]
    for row in completed:
        cell = row["cell"]
        combined, wdsat, cms = row["combined_exporter_process"], row["wdsat"], row["cryptominisat"]
        def resource_cell(core, rss) -> str:
            if core is None or rss is None:
                return "-"
            return f"{core:.6f} / {rss / 1048576:.2f}"
        lines.append(
            f"| {cell['n']} | {cell['basis']} | "
            f"{resource_cell(combined['total_core_seconds'], combined['peak_rss_bytes'])} | "
            f"{resource_cell(wdsat['total_core_seconds'], wdsat['peak_rss_bytes'])} | "
            f"{resource_cell(wdsat['build_total_core_seconds'], wdsat['build_peak_rss_bytes'])} | "
            f"{resource_cell(cms['total_core_seconds'], cms['peak_rss_bytes'])} |"
        )
    lines += [
        "",
        "At n=31 and n=41, every available solver returned a source-ANF-validated SAT model. Source validation alone does not prove that the model coordinates lift to curve points. At n=59, direct MITM found the planted decomposition in about 3.13 seconds; native SAT stopped at 100,000 conflicts, WDSat reached the 120-second watchdog, and CryptoMiniSat reached the same watchdog. Both capped solver outcomes are inconclusive.",
        "",
        "The n=31 GGMP cell derives its five-dimensional x-domain as the kernel of a linearised polynomial obtained from a factor of T^31-1. Independent review found that only T=(0,1) lifts from this domain, the planted relation is T+T+T=T, and both external models contain the nonlifting coordinate x=65536. This row is a degenerate algebraic-model/exporter check, not a meaningful GGMP point-decomposition benchmark. The generator now rejects any such domain with fewer than three distinct curve points.",
        "",
        "Magma is not installed on this host, so the generated `.magma` inputs are retained but F4 is an operationally missing baseline. The n=67 cell exceeds the current n<64 field-bitmask implementation and asserts nothing about PDP hardness.",
        "",
        "Direct MITM confirms planted point decompositions through n=59, but the SAT outputs are not yet lift-validated end to end. Fully isolated process accounting for every PDP arm, Magma F4, repeated scaling, external reproduction, and novelty review remain open. The result is not a Koblitz index-calculus SOTA claim.",
        "",
        "Primary comparisons: Trimoska-Ionica-Dequen, *A SAT-Based Approach for Index Calculus on Binary Elliptic Curves* (ePrint 2019/313); Galbraith-Granger-Merz-Petit, *On Index Calculus Algorithms for Subfield Curves* (ePrint 2020/1315).",
        "",
    ]
    (ROOT / "STAGE1_RESULTS.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
