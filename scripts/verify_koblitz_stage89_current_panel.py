#!/usr/bin/env python3
"""Compose and verify the current selected n=53 panel and unknown-scalar run."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import statistics
import tempfile
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage33_n41_unknown_scalar as stage33
import run_koblitz_stage42_n53_same_target as stage42
import run_koblitz_stage57_n53_unknown as stage57
import run_koblitz_stage79_specialized_pair_batch as stage79
import run_koblitz_stage81_specialized_pair_sums as stage81
import verify_koblitz_stage49_width_results as common


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
PREDECESSOR = STAGE / "stage-73-current-gate-audit-20260912/audit.json"
DEFAULT_OUTPUT = STAGE / "stage-89-current-selected-panel-20260912"
SCHEMA = "koblitz_stage89_current_selected_panel.v1"
SEAL_SCHEMA = "koblitz_stage89_current_selected_panel_seal.v1"

CONFIG: dict[int, dict[str, Any]] = {
    83: {
        "run_id": 34724087748,
        "commit": "22a3c7ac54265b0557dc3b763dfbb9896648119b",
        "artifact_id": 10307027649,
        "artifact_digest": "sha256:681b5b9c9a43e45e774e68f49e37aa0256c2344e37e73378a8885e21abab38c4",
        "archive": "koblitz-stage83-selected-run1-evidence-34724087748.tar.gz",
        "archive_root": "koblitz-stage83-selected-run1-34724087748",
        "prefix": "stage79",
        "job": "n53-specialized-pair-batch",
        "runner": stage79,
        "selected_role": "candidate",
    },
    84: {
        "run_id": 34724270013,
        "commit": "22a3c7ac54265b0557dc3b763dfbb9896648119b",
        "artifact_id": 10306169449,
        "artifact_digest": "sha256:ba27cd6549c6159c1d0f358ab51b1c750899a1cc7ceba5b3b1a2ff7071d7c6c2",
        "archive": "koblitz-stage84-selected-run2-evidence-34724270013.tar.gz",
        "archive_root": "koblitz-stage84-selected-run2-34724270013",
        "prefix": "stage79",
        "job": "n53-specialized-pair-batch",
        "runner": stage79,
        "selected_role": "candidate",
    },
    85: {
        "run_id": 34725484816,
        "commit": "6feb8a5807fc6d3bfacc4bbe0437a5890ac357fa",
        "artifact_id": 10307647532,
        "artifact_digest": "sha256:2cca9135b6899bf708be3bc13ef6884da6f202b185c562214a1ca51e230bb86f",
        "archive": "koblitz-stage85-selected-run3-evidence-34725484816.tar.gz",
        "archive_root": "koblitz-stage85-selected-run3-34725484816",
        "prefix": "stage81",
        "job": "n53-specialized-pair-sums",
        "runner": stage81,
        "selected_role": "baseline",
    },
    86: {
        "run_id": 34725643637,
        "commit": "6feb8a5807fc6d3bfacc4bbe0437a5890ac357fa",
        "artifact_id": 10308085686,
        "artifact_digest": "sha256:b8d0ce2bed0a2a3bf5ab6d9a876f8af504fd8762485463bcbe1b506403d0ff04",
        "archive": "koblitz-stage86-selected-run4-evidence-34725643637.tar.gz",
        "archive_root": "koblitz-stage86-selected-run4-34725643637",
        "prefix": "stage81",
        "job": "n53-specialized-pair-sums",
        "runner": stage81,
        "selected_role": "baseline",
    },
    87: {
        "run_id": 34725960358,
        "commit": "8190f5bee6fdcee2f989125083eb2999229f4142",
        "artifact_id": 10307144266,
        "artifact_digest": "sha256:50ad876d57d6feb86266e77322d71a311f0f9793dc271bc898e3542ad5738d06",
        "archive": "koblitz-stage87-selected-current-evidence-34725960358.tar.gz",
        "archive_root": "koblitz-stage87-selected-current-34725960358",
        "prefix": "stage79",
        "job": "n53-specialized-pair-batch",
        "runner": stage79,
        "selected_role": "candidate",
    },
}

UNKNOWN_CONFIG: dict[str, Any] = {
    "run_id": 34724295607,
    "commit": "22a3c7ac54265b0557dc3b763dfbb9896648119b",
    "artifact_id": 10308010013,
    "artifact_digest": "sha256:2ce9f565805ca0ff8571ccc2f7895d871185053704acf06f276283bb3cb1b53b",
    "archive": "koblitz-stage88-unknown-current-evidence-34724295607.tar.gz",
    "archive_root": "koblitz-stage88-unknown-current-34724295607",
    "prefix": "stage69",
    "job": "n53-unknown-optimized",
}


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def selected_result(stage: int, config: dict[str, Any]) -> dict[str, Any]:
    with tempfile.TemporaryDirectory() as directory:
        root = common.extract(stage, config, Path(directory))
        prefix = config["prefix"]
        verification = config["runner"].verify(
            root / f"{prefix}-build", root / f"{prefix}-run"
        )
        embedded = load(
            root / f"{prefix}-verification.json", f"embedded Stage-{stage} verification"
        )
        require(verification == embedded, f"Stage-{stage} embedded verification changed")
        run = load(root / f"{prefix}-run/result.json", f"Stage-{stage} run")
        build_seal = custody.sha256_file(
            root / f"{prefix}-build/result-seal.json", f"Stage-{stage} build seal"
        )
        run_seal = custody.sha256_file(
            root / f"{prefix}-run/result-seal.json", f"Stage-{stage} run seal"
        )
    role = config["selected_role"]
    observation = run[role]
    full_base = observation["base"]
    full_summary = observation["summary"]
    base_keys = (
        "base_hash",
        "factor_base_points",
        "orbit_columns",
        "field_x_values_scanned",
        "support_index_entries",
        "support_table_allocated_bytes",
        "support_x_prefilter_bits",
        "support_x_prefilter_kind",
        "support_x_prefilter_hash_strategy",
        "support_x_prefilter_insert_hash_reuse",
        "parallel_support_expansion",
        "pipelined_support_expansion",
        "specialized_n53_support_expansion",
        "specialized_n53_pair_sums",
        "base_construction_ms",
        "support_index_ms",
        "total_setup_ms",
        "field_mul_backend",
        "field_square_backend",
        "field_reduction_backend",
        "field_product_pipeline",
        "field_product_dispatch",
        "selection_uses_scalar_labels",
        "subgroup_membership_verified",
        "frobenius_closed",
        "negation_closed",
    )
    summary_keys = (
        "published_q",
        "published_q_point_key",
        "recovered_fixture_scalar",
        "factor_base_logs_known_by_construction",
        "linear_solution_verified",
        "all_relations_group_verified",
        "query_mode",
        "query_parallel_threads",
        "query_parallel_waves",
        "query_parallel_chunks",
        "query_specialized_n53_pair_batch",
        "query_compact_pair_scratch",
        "query_dual_sign_denominator_sharing",
        "query_shared_sign_denominator_inputs",
        "factor_base_points",
        "orbit_columns",
        "matrix_columns",
        "admitted_relations",
        "target_trials",
        "full_rank_at_relation",
        "setup_ms",
        "collection_ms",
        "charged_total_ms",
        "support_queries",
        "query_batch_inversions",
        "query_x_filter_rejections",
        "query_exact_table_misses",
        "rank_target_slot_index_allocated_bytes",
        "linear_solve_ms",
        "solution_validation_ms",
        "reference_validation_ms",
        "field_mul_backend",
        "field_square_backend",
        "field_reduction_backend",
        "field_product_pipeline",
        "field_product_dispatch",
    )
    base = {key: full_base.get(key) for key in base_keys}
    summary = {key: full_summary.get(key) for key in summary_keys}
    process = run[f"{role}_process"]["metrics"]
    rho_process = run["rho_process"]["metrics"]
    ratio = process["wall_seconds"] / rho_process["wall_seconds"]
    require(summary.get("admitted_relations") == 95, f"Stage-{stage} relation count changed")
    require(summary.get("full_rank_at_relation") == 95, f"Stage-{stage} rank endpoint changed")
    require(summary.get("support_queries") == 48_531_878, f"Stage-{stage} support queries changed")
    require(summary.get("query_exact_table_misses") == 377_803, f"Stage-{stage} exact misses changed")
    require(summary.get("query_specialized_n53_pair_batch") is True, f"Stage-{stage} selected query batch changed")
    require(summary.get("query_compact_pair_scratch") is True, f"Stage-{stage} compact scratch changed")
    require(summary.get("query_dual_sign_denominator_sharing") is True, f"Stage-{stage} dual-sign sharing changed")
    require(summary.get("linear_solution_verified") is True, f"Stage-{stage} solution failed")
    require(summary.get("all_relations_group_verified") is True, f"Stage-{stage} relation verification failed")
    return {
        "source_commit": config["commit"],
        "selected_role": role,
        "archive": common.archive_identity(stage, config),
        "workflow": common.workflow_accounting(stage, config),
        "artifact": common.artifact_accounting(stage, config),
        "build_seal_sha256": build_seal,
        "run_seal_sha256": run_seal,
        "verification": verification,
        "base": base,
        "summary": summary,
        "process": process,
        "rho_process": rho_process,
        "direct_over_rho_wall_ratio": ratio,
        "whole_process_crossover": ratio < 1,
        "build_outer": run["build_binding"]["outer_resources"],
        "science_outer": run["outer_resources"],
        "charged_total_core_seconds_available": run[
            "charged_total_core_seconds_available"
        ],
        "maximum_sampled_process_tree_rss_bytes": run[
            "maximum_sampled_process_tree_rss_bytes"
        ],
    }


def unknown_result() -> dict[str, Any]:
    stage = 88
    config = UNKNOWN_CONFIG
    with tempfile.TemporaryDirectory() as directory:
        root = common.extract(stage, config, Path(directory))
        # Stage 88 predates host-identity and sharded-selection fields now
        # required by the live Stage 57 runner. Its archive, build, and run
        # seals remain authoritative for this historical replay.
        stage42.validate_build(root / "stage69-build")
        stage33.verify_seal(
            root / "stage69-run", stage57.RUN_SEAL_SCHEMA, "run_frozen"
        )
        embedded = load(root / "stage69-verification.json", "embedded Stage-88 verification")
        verification = embedded
        run = load(root / "stage69-run/result.json", "Stage-88 run")
        build_seal = custody.sha256_file(
            root / "stage69-build/result-seal.json", "Stage-88 build seal"
        )
        run_seal = custody.sha256_file(
            root / "stage69-run/result-seal.json", "Stage-88 run seal"
        )
    summary = run["direct_summary"]
    require(run.get("optimized_stack") is True, "Stage-88 optimized stack changed")
    require(verification.get("target_scalar_constructed_or_supplied") is False, "Stage-88 target scalar was supplied")
    require(verification.get("factor_base_logs_known_by_construction") is False, "Stage-88 factor-base logs were supplied")
    require(summary.get("query_specialized_n53_pair_batch") is True, "Stage-88 specialized batch changed")
    require(summary.get("admitted_relations") == summary.get("full_rank_at_relation") == 95, "Stage-88 rank endpoint changed")
    require(run.get("recovered_scalar") == 7_892_094_459_170, "Stage-88 recovered scalar changed")
    return {
        "source_commit": config["commit"],
        "archive": common.archive_identity(stage, config),
        "workflow": common.workflow_accounting(stage, config),
        "artifact": common.artifact_accounting(stage, config),
        "build_seal_sha256": build_seal,
        "run_seal_sha256": run_seal,
        "verification": verification,
        "summary": summary,
        "rho": run["rho"],
        "direct_process": run["direct_process"]["metrics"],
        "rho_process": run["rho_process"]["metrics"],
        "comparison": run["comparison"],
        "build_outer": run["build_binding"]["outer_resources"],
        "science_outer": run["outer_resources"],
        "charged_total_core_seconds_available": run[
            "charged_total_core_seconds_available"
        ],
        "maximum_sampled_process_tree_rss_bytes": run[
            "maximum_sampled_process_tree_rss_bytes"
        ],
    }


def result() -> dict[str, Any]:
    predecessor = load(PREDECESSOR, "Stage-73 predecessor")
    require(
        predecessor.get("schema") == "koblitz_stage73_current_gate_audit.v1",
        "Stage-73 predecessor changed",
    )
    selected = {str(stage): selected_result(stage, config) for stage, config in CONFIG.items()}
    ratios = [selected[str(stage)]["direct_over_rho_wall_ratio"] for stage in CONFIG]
    expected = [
        0.9464834379305384,
        1.019315178402027,
        1.0949740741801435,
        0.932363424377663,
        1.0187756063132052,
    ]
    require(ratios == expected, "selected-stack ratio panel changed")
    median = statistics.median(ratios)
    wins = sum(ratio < 1 for ratio in ratios)
    require(wins == 2 and median == 1.0187756063132052, "selected-stack panel conclusion changed")
    current = selected["87"]
    base = current["base"]
    summary = current["summary"]
    require(current["source_commit"] == "8190f5bee6fdcee2f989125083eb2999229f4142", "selected current commit changed")
    require(base.get("specialized_n53_support_expansion") is False, "variable support expansion became selected")
    require(base.get("specialized_n53_pair_sums") is False, "variable pair sums became selected")
    require(base.get("support_x_prefilter_hash_strategy") == "single_mix_split_29_bit_indices", "selected Bloom strategy changed")
    require(base.get("support_x_prefilter_insert_hash_reuse") is True, "selected insert-hash reuse changed")
    unknown = unknown_result()
    return {
        "schema": SCHEMA,
        "status": "current_selected_panel_and_unknown_verified",
        "predecessor_stage73": {
            "path": str(PREDECESSOR.relative_to(REPO)),
            "sha256": custody.sha256_file(PREDECESSOR, "Stage-73 predecessor"),
        },
        "selected_runs": selected,
        "selected_panel": {
            "run_count": len(ratios),
            "ratios": ratios,
            "wins": wins,
            "losses": len(ratios) - wins,
            "median_direct_over_rho_wall_ratio": median,
            "repeatable_whole_process_crossover": False,
        },
        "current_selected": current,
        "unknown_scalar": unknown,
        "factor_base_algebraically_point_defined": True,
        "factor_base_selection_uses_scalar_labels": False,
        "target_subgroup_enumerated_for_factor_base": False,
        "factor_base_logs_known_by_construction": False,
        "unknown_scalar_recovery_verified": True,
        "same_public_target_panel": True,
        "exact_relation_hash_equality": True,
        "sat_conflicts": None,
        "licensed_magma_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "remaining_blockers": [
            "selected n=53 stack wins two of five matched runs but has a 1.018776 median direct/rho wall ratio",
            "fresh-build cost remains far above one rho target",
            "current selected-stack forced-single-core evidence remains absent",
            "licensed same-instance Magma and unaffiliated reproduction/novelty review remain absent",
            "the complete growing-n WDSat/CryptoMiniSat/Magma/MITM/GGMP panel remains incomplete",
        ],
    }


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-89 output must be new")
    output.mkdir()
    custody.write_json_new(output / "verification.json", result())
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    payload = {
        "schema": SEAL_SCHEMA,
        "status": "current_selected_panel_frozen",
        "inventory": inventory,
        "inventory_sha256": custody.canonical_sha256(inventory),
    }
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "result-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "verification.json", "committed Stage-89 result")
    require(committed == result(), "committed Stage-89 result differs from source evidence")
    seal = load(output / "result-seal.json", "Stage-89 result seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(
        seal.get("schema") == SEAL_SCHEMA
        and claimed == custody.canonical_sha256(payload),
        "Stage-89 result seal changed",
    )
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    require(
        inventory == seal["inventory"]
        and custody.canonical_sha256(inventory) == seal["inventory_sha256"],
        "Stage-89 result inventory changed",
    )
    return committed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    compose = sub.add_parser("compose")
    compose.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    check = sub.add_parser("verify")
    check.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    try:
        value = (
            freeze(args.output.resolve())
            if args.command == "compose"
            else verify(args.output.resolve())
        )
        print(json.dumps(value, indent=2, sort_keys=True))
    except (
        OSError,
        ValueError,
        KeyError,
        VerificationError,
        common.VerificationError,
        stage57.Stage57Error,
        stage79.Stage79Error,
        stage81.Stage81Error,
        custody.PhaseBError,
    ) as error:
        raise SystemExit(f"stage89-current-panel: {error}")


if __name__ == "__main__":
    main()
