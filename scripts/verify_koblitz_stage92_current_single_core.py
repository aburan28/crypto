#!/usr/bin/env python3
"""Compose and verify the current selected-stack n=53 CPU-0 result."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage33_n41_unknown_scalar as stage33
import run_koblitz_stage42_n53_same_target as stage42
import run_koblitz_stage71_single_core as stage71
import verify_koblitz_stage49_width_results as common


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
PREDECESSOR = STAGE / "stage-90-current-gate-audit-20260912/audit.json"
DEFAULT_OUTPUT = STAGE / "stage-92-current-single-core-20260913"
SCHEMA = "koblitz_stage92_current_single_core.v1"
SEAL_SCHEMA = "koblitz_stage92_current_single_core_seal.v1"
CONFIG: dict[str, Any] = {
    "run_id": 34726819616,
    "commit": "8190f5bee6fdcee2f989125083eb2999229f4142",
    "artifact_id": 10307264889,
    "artifact_digest": "sha256:0dfb12d5a4d6cf42c46f4fbfcfaa1bc3a2dfe1e956d1dd4caecbee54ff17d49e",
    "archive": "koblitz-stage91-current-single-core-evidence-34726819616.tar.gz",
    "archive_root": "koblitz-stage91-current-single-core-34726819616",
    "job": "n53-single-core",
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


def result() -> dict[str, Any]:
    predecessor = load(PREDECESSOR, "Stage-90 predecessor")
    require(
        predecessor.get("schema") == "koblitz_stage90_current_gate_audit.v1",
        "Stage-90 predecessor changed",
    )
    with tempfile.TemporaryDirectory() as directory:
        root = common.extract(91, CONFIG, Path(directory))
        # Stage 91 is a frozen pre-sharding CPU-0 receipt. Validate its
        # original build/run seals and embedded result without requiring host
        # fields added later to the live Stage 71 runner.
        stage42.validate_build(root / "stage71-build")
        stage33.verify_seal(
            root / "stage71-run", stage71.RUN_SEAL_SCHEMA, "run_frozen"
        )
        embedded = load(root / "stage71-verification.json", "embedded Stage-91 verification")
        verification = embedded
        run = load(root / "stage71-run/result.json", "Stage-91 run")
        build_seal = custody.sha256_file(
            root / "stage71-build/result-seal.json", "Stage-91 build seal"
        )
        run_seal = custody.sha256_file(
            root / "stage71-run/result-seal.json", "Stage-91 run seal"
        )
    affinity = run["affinity"]
    direct = run["direct"]
    base = direct["base"]
    summary = direct["summary"]
    require(affinity.get("inheritance_verified") is True, "CPU affinity inheritance failed")
    require(
        affinity.get("parent_effective_cpus")
        == affinity.get("child_effective_cpus")
        == [0],
        "effective CPU set changed",
    )
    require(summary.get("query_parallel_threads") == 1, "query thread count changed")
    require(
        summary.get("query_specialized_n53_pair_batch") is True,
        "selected specialized query batch changed",
    )
    require(summary.get("query_compact_pair_scratch") is True, "compact scratch changed")
    require(
        summary.get("query_dual_sign_denominator_sharing") is True,
        "dual-sign sharing changed",
    )
    require(
        base.get("specialized_n53_support_expansion") is False
        and base.get("specialized_n53_pair_sums") is False,
        "variable support specialization became selected",
    )
    require(
        summary.get("admitted_relations") == summary.get("full_rank_at_relation") == 95,
        "rank endpoint changed",
    )
    require(summary.get("linear_solution_verified") is True, "solution failed")
    base_keys = (
        "base_hash",
        "factor_base_points",
        "orbit_columns",
        "support_index_entries",
        "support_table_allocated_bytes",
        "support_x_prefilter_hash_strategy",
        "support_x_prefilter_insert_hash_reuse",
        "specialized_n53_support_expansion",
        "specialized_n53_pair_sums",
        "field_product_pipeline",
        "field_product_dispatch",
    )
    summary_keys = (
        "published_q",
        "recovered_fixture_scalar",
        "query_mode",
        "query_parallel_threads",
        "query_specialized_n53_pair_batch",
        "query_compact_pair_scratch",
        "query_dual_sign_denominator_sharing",
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
        "query_exact_table_misses",
        "linear_solve_ms",
        "solution_validation_ms",
        "reference_validation_ms",
        "linear_solution_verified",
        "all_relations_group_verified",
    )
    compact_base = {key: base.get(key) for key in base_keys}
    compact_summary = {key: summary.get(key) for key in summary_keys}
    comparison = run["comparison"]
    require(
        comparison.get("direct_over_rho_wall_ratio") == 1.3543032952267562,
        "single-core wall ratio changed",
    )
    require(
        comparison.get("direct_core_over_rho_core_ratio") == 1.3401796163922106,
        "single-core CPU ratio changed",
    )
    return {
        "schema": SCHEMA,
        "status": "current_selected_n53_cpu0_result_verified",
        "source_commit": CONFIG["commit"],
        "predecessor_stage90": {
            "path": str(PREDECESSOR.relative_to(REPO)),
            "sha256": custody.sha256_file(PREDECESSOR, "Stage-90 audit"),
        },
        "archive": common.archive_identity(91, CONFIG),
        "workflow": common.workflow_accounting(91, CONFIG),
        "artifact": common.artifact_accounting(91, CONFIG),
        "build_seal_sha256": build_seal,
        "run_seal_sha256": run_seal,
        "verification": verification,
        "affinity": affinity,
        "base": compact_base,
        "summary": compact_summary,
        "rho": run["rho"],
        "direct_process": run["direct_process"]["metrics"],
        "rho_process": run["rho_process"]["metrics"],
        "comparison": comparison,
        "build_outer": run["build_binding"]["outer_resources"],
        "science_outer": run["outer_resources"],
        "charged_total_core_seconds_available": run[
            "charged_total_core_seconds_available"
        ],
        "charged_sequential_wall_seconds_available": run[
            "charged_sequential_wall_seconds_available"
        ],
        "maximum_sampled_process_tree_rss_bytes": run[
            "maximum_sampled_process_tree_rss_bytes"
        ],
        "factor_base_selection_uses_scalar_labels": False,
        "target_subgroup_enumerated_for_factor_base": False,
        "factor_base_logs_known_by_construction": False,
        "licensed_magma_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-92 output must be new")
    output.mkdir()
    custody.write_json_new(output / "verification.json", result())
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    payload = {
        "schema": SEAL_SCHEMA,
        "status": "current_single_core_frozen",
        "inventory": inventory,
        "inventory_sha256": custody.canonical_sha256(inventory),
    }
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "result-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "verification.json", "committed Stage-92 result")
    require(committed == result(), "committed Stage-92 result differs from source evidence")
    seal = load(output / "result-seal.json", "Stage-92 result seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(
        seal.get("schema") == SEAL_SCHEMA
        and claimed == custody.canonical_sha256(payload),
        "Stage-92 result seal changed",
    )
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    require(
        inventory == seal["inventory"]
        and custody.canonical_sha256(inventory) == seal["inventory_sha256"],
        "Stage-92 result inventory changed",
    )
    return committed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    create = sub.add_parser("compose")
    create.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
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
        stage71.Stage71Error,
        custody.PhaseBError,
    ) as error:
        raise SystemExit(f"stage92-single-core: {error}")


if __name__ == "__main__":
    main()
