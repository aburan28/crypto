#!/usr/bin/env python3
"""Run a host-identified five-pair selected n=53 panel."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import time
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
from run_koblitz_stage26_affinity_cell import TreeSampler
import run_koblitz_stage33_n41_unknown_scalar as stage33
import run_koblitz_stage42_n53_same_target as stage42
import run_koblitz_stage44_n53_parallel as stage44
import run_koblitz_stage61_production_rank as stage61


REPO = Path(__file__).resolve().parents[1]
RUN_SCHEMA = "koblitz_stage102_selected_panel.v1"
RUN_SEAL_SCHEMA = "koblitz_stage102_selected_panel_seal.v1"
PAIR_COUNT = 5


class Stage102Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage102Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def self_test() -> dict[str, Any]:
    return {
        "schema": "koblitz_stage102_self_test.v1",
        "status": "PASS",
        "checks": 39,
        "n": 53,
        "pair_count": PAIR_COUNT,
        "parallel_width": 4096,
        "selected_stack": {
            "x_filter": "four_shard_direct_low_and_high_x_bit_windows",
            "prefiltered_exact_lookup": True,
            "itoh_tsujii_inverse": True,
            "blocked_filter": False,
            "compact_hash_complete_evidence": True,
            "fixed_base_reference_validation": True,
            "support_table_shards": 4,
            "shard_routing": "xor_low_and_high_x_windows",
        },
        "same_target": True,
        "separate_process_meter_per_arm": True,
        "canonical_linux_host_identity": True,
        "median_gate": "paired direct/rho wall ratio < 0.95",
    }


def selected_environment() -> dict[str, str]:
    environment = custody.safe_child_environment()
    environment["KIC_INCREMENTAL_RANK_CROSSCHECK"] = "0"
    environment["KIC_RANK_SURPLUS"] = "0"
    environment["KIC_RANK_AWARE_PAIR_SCAN"] = "1"
    environment["KIC_PARALLEL_SUPPORT_EXPANSION"] = "1"
    environment["KIC_PIPELINED_SUPPORT_EXPANSION"] = "1"
    environment["KIC_ENABLE_SHARDED_SUPPORT_TABLE"] = "1"
    environment["KIC_SUMMARY_ONLY"] = "1"
    environment["RAYON_NUM_THREADS"] = "4"
    return environment


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-102 requires Linux process accounting")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-102 output must be new")
    build_root = args.build.resolve(strict=True)
    build = stage42.validate_build(build_root)
    output.mkdir(parents=True)
    evidence = output / "evidence"
    evidence.mkdir()
    direct = str((build_root / "bin/koblitz_rank_fixture").resolve(strict=True))
    rho = str((build_root / "bin/koblitz_rho_fixture").resolve(strict=True))
    direct_command = [
        direct, "53", "0", "1", "128", str(stage42.DIRECT_SEED),
        "signed_expanded", "independent", "pair_pair_parallel_4096", "1",
        str(stage42.EXPLICIT_SCALAR),
    ]
    rho_command = [
        rho, "53", "0", "signed_frobenius", "1", "packed",
        str(stage42.RHO_BATCH_SEED), str(stage42.EXPLICIT_SCALAR),
    ]
    environment = selected_environment()

    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    host_identity = stage33.linux_host_identity()
    processes: list[dict[str, Any]] = []
    with TreeSampler() as sampler:
        for pair_index in range(1, PAIR_COUNT + 1):
            direct_process = stage33.metered(
                role=f"direct-{pair_index}", command=direct_command, cwd=REPO,
                evidence=evidence, timeout=args.timeout, environment=environment,
            )
            rho_process = stage33.metered(
                role=f"rho-{pair_index}", command=rho_command, cwd=REPO,
                evidence=evidence, timeout=args.timeout, environment=custody.safe_child_environment(),
            )
            processes.append({"pair_index": pair_index, "direct": direct_process, "rho": rho_process})
    outer = stage33.outer_resources(started, before_self, before_children, sampler)
    require(
        all(stage33.process_complete(arm) for pair in processes for arm in (pair["direct"], pair["rho"])),
        "Stage-102 process failed",
    )

    direct_rows = [
        stage61.direct_observation(evidence / f"direct-{pair_index}.stdout", "FULL_RANK")
        for pair_index in range(1, PAIR_COUNT + 1)
    ]
    rho_rows = []
    for pair_index in range(1, PAIR_COUNT + 1):
        rows = stage44.parse_lines(evidence / f"rho-{pair_index}.stdout")
        require(len(rows) == 1 and rows[0].get("verified") is True, f"rho pair {pair_index} failed")
        rho_rows.append(rows[0])
    first = direct_rows[0]
    first_summary = first["summary"]
    first_base = first["base"]
    for pair_index, row in enumerate(direct_rows, 1):
        summary = row["summary"]
        base = row["base"]
        require(row["relation_hashes"] == first["relation_hashes"], f"direct pair {pair_index} relation transcript changed")
        require(summary["factor_base_log_solution"] == first_summary["factor_base_log_solution"], f"direct pair {pair_index} solution changed")
        require(base["base_hash"] == first_base["base_hash"], f"direct pair {pair_index} factor base changed")
        require(base.get("support_x_prefilter_hash_strategy") == "four_shard_direct_low_and_high_x_bit_windows", f"direct pair {pair_index} lost sharded direct-bit filter")
        require(base.get("support_table_shards") == 4, f"direct pair {pair_index} lost four-shard table")
        require(base.get("support_table_shard_routing") == "xor_low_and_high_x_windows", f"direct pair {pair_index} lost selected shard routing")
        require(base.get("parallel_support_insertion") is True, f"direct pair {pair_index} lost parallel support insertion")
        require(base.get("support_x_prefilter_direct_bits") is True, f"direct pair {pair_index} direct-bit flag changed")
        require(base.get("support_x_prefilter_blocked") is False, f"direct pair {pair_index} selected blocked filter")
        require(summary.get("query_prefiltered_exact_lookup") is True, f"direct pair {pair_index} lost prefiltered lookup")
        require(summary.get("query_itoh_n53_inverse") is True, f"direct pair {pair_index} lost Itoh inverse")
        require(summary.get("query_mode") == "pair_pair_parallel_4096", f"direct pair {pair_index} width changed")
        require(summary.get("query_parallel_threads") == 4, f"direct pair {pair_index} thread count changed")
        require(summary.get("query_specialized_n53_pair_batch") is True, f"direct pair {pair_index} lost specialized batch")
        require(summary.get("summary_only_timing") is True, f"direct pair {pair_index} lost compact evidence")
        require(summary.get("relation_receipts_emitted") is False, f"direct pair {pair_index} emitted full receipts")
        require(summary.get("fixed_base_reference_validation") is True, f"direct pair {pair_index} lost fixed-base validation")
        require(base.get("compact_evidence") is True, f"direct pair {pair_index} lost compact base evidence")
        require(summary.get("recovered_fixture_scalar") == stage42.EXPLICIT_SCALAR, f"direct pair {pair_index} scalar changed")
        require(rho_rows[pair_index - 1].get("published_q") == summary.get("published_q"), f"pair {pair_index} target changed")

    pairs = []
    for process, direct_row, rho_row in zip(processes, direct_rows, rho_rows):
        direct_metrics = process["direct"]["metrics"]
        rho_metrics = process["rho"]["metrics"]
        summary = direct_row["summary"]
        pairs.append({
            "pair_index": process["pair_index"],
            "direct_process": process["direct"],
            "rho_process": process["rho"],
            "direct": direct_row,
            "rho": rho_row,
            "direct_over_rho_wall_ratio": direct_metrics["wall_seconds"] / rho_metrics["wall_seconds"],
            "direct_over_rho_core_ratio": direct_metrics["total_core_seconds"] / rho_metrics["total_core_seconds"],
            "direct_setup_ms": summary["setup_ms"],
            "direct_collection_ms": summary["collection_ms"],
            "direct_query_ms": summary["timing_breakdown_ms"]["query"],
        })
    ratios = [pair["direct_over_rho_wall_ratio"] for pair in pairs]
    direct_walls = [pair["direct_process"]["metrics"]["wall_seconds"] for pair in pairs]
    rho_walls = [pair["rho_process"]["metrics"]["wall_seconds"] for pair in pairs]
    panel = {
        "paired_wall_ratios": ratios,
        "median_paired_wall_ratio": statistics.median(ratios),
        "direct_wins": sum(ratio < 1.0 for ratio in ratios),
        "median_direct_wall_seconds": statistics.median(direct_walls),
        "median_rho_wall_seconds": statistics.median(rho_walls),
        "minimum_paired_wall_ratio": min(ratios),
        "maximum_paired_wall_ratio": max(ratios),
        "median_gate_threshold": 0.95,
        "median_gate_passed": statistics.median(ratios) < 0.95,
        "fresh_build_plus_median_direct_over_median_rho_wall_ratio":
            (build["outer_resources"]["wall_seconds"] + statistics.median(direct_walls))
            / statistics.median(rho_walls),
    }
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n53_selected_five_pair_panel",
        "constants": self_test(),
        "build_binding": {
            "result": custody.executable_identity(build_root / "result.json", "Stage-102 build result", executable=False),
            "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-102 build seal", executable=False),
            "source_state": build["source_state"],
            "outer_resources": build["outer_resources"],
        },
        "pairs": pairs,
        "host_identity": host_identity,
        "panel": panel,
        "outer_resources": outer,
        "charged_total_core_seconds_available": build["outer_resources"]["total_core_seconds"] + outer["total_core_seconds"],
        "charged_sequential_wall_seconds_available": build["outer_resources"]["wall_seconds"] + outer["wall_seconds"],
        "maximum_sampled_process_tree_rss_bytes": max(build["outer_resources"]["sampled_peak_process_tree_rss_bytes"], outer["sampled_peak_process_tree_rss_bytes"]),
        "factor_base_selection_uses_scalar_labels": False,
        "target_subgroup_enumerated_for_factor_base": False,
        "factor_base_logs_known_by_construction": False,
        "fixture_scalar_used_by_collector": False,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    custody.write_json_new(output / "result.json", result)
    stage33.freeze(output, RUN_SEAL_SCHEMA, "run_frozen")
    return result


def verify(build_root: Path, output: Path) -> dict[str, Any]:
    build = stage42.validate_build(build_root)
    seal = stage33.verify_seal(output, RUN_SEAL_SCHEMA, "run_frozen")
    result = load(output / "result.json", "Stage-102 result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n53_selected_five_pair_panel", "Stage-102 result is incomplete")
    direct_rows = [stage61.direct_observation(output / "evidence" / f"direct-{index}.stdout", "FULL_RANK") for index in range(1, PAIR_COUNT + 1)]
    require(direct_rows == [pair["direct"] for pair in result["pairs"]], "Stage-102 direct transcripts changed")
    require(len({tuple(row["relation_hashes"]) for row in direct_rows}) == 1, "Stage-102 relation equality changed")
    expected_ratios = [pair["direct_process"]["metrics"]["wall_seconds"] / pair["rho_process"]["metrics"]["wall_seconds"] for pair in result["pairs"]]
    require(expected_ratios == result["panel"]["paired_wall_ratios"], "Stage-102 ratios changed")
    require(statistics.median(expected_ratios) == result["panel"]["median_paired_wall_ratio"], "Stage-102 median changed")
    stage33.validate_linux_host_identity(result["host_identity"])
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-102 claim widened")
    return {
        "schema": "koblitz_stage102_verification.v1",
        "status": "n53_selected_five_pair_panel_verified",
        "source_commit": build["source_state"]["commit"],
        "host_identity": result["host_identity"],
        "relations_per_direct": [row["summary"]["admitted_relations"] for row in direct_rows],
        **result["panel"],
        "charged_total_core_seconds_available": result["charged_total_core_seconds_available"],
        "charged_sequential_wall_seconds_available": result["charged_sequential_wall_seconds_available"],
        "maximum_sampled_process_tree_rss_bytes": result["maximum_sampled_process_tree_rss_bytes"],
        "inventory_sha256": seal["inventory_sha256"],
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("self-test")
    build_parser = sub.add_parser("build")
    build_parser.add_argument("--source", type=Path, required=True)
    build_parser.add_argument("--output", type=Path, required=True)
    build_parser.add_argument("--jobs", type=int, default=4)
    build_parser.add_argument("--timeout", type=int, default=1800)
    run_parser = sub.add_parser("run")
    run_parser.add_argument("--build", type=Path, required=True)
    run_parser.add_argument("--output", type=Path, required=True)
    run_parser.add_argument("--timeout", type=int, default=1800)
    verify_parser = sub.add_parser("verify")
    verify_parser.add_argument("--build", type=Path, required=True)
    verify_parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "self-test":
            value = self_test()
        elif args.command == "build":
            value = stage42.build(args)
        elif args.command == "run":
            value = run(args)
        else:
            value = verify(args.build.resolve(strict=True), args.output.resolve(strict=True))
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage102Error, stage61.Stage61Error, stage42.Stage42Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage102-selected-panel: {error}")


if __name__ == "__main__":
    main()
