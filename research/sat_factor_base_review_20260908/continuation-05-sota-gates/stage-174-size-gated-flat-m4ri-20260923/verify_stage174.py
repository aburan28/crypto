#!/usr/bin/env python3
"""Independently verify the committed Stage 174 evidence."""

from __future__ import annotations

import copy
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import subprocess


HERE = Path(__file__).resolve().parent
GATES = HERE.parent
DEV = HERE / "development"
RESULT = HERE / "result.json"
STAGE171 = GATES / "stage-171-dense-symbolic-sets-20260923/result.json"
STAGE173 = GATES / "stage-173-ggmp-same-cell-availability-20260923/result.json"
COMMIT = "e51efb219edd4511a08d545de21184928e5e771c"
TREE = "321c58f585735c3e2c1d430e4c39a2555b654838"
FINGERPRINT = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_report(report: dict) -> dict:
    value = copy.deepcopy(report)
    for key in (
        "timing_ns",
        "solver_description",
        "solver_echelon_policy",
        "solver_flat_m4ri_control",
        "solver_flat_m4ri_min_words",
    ):
        value.pop(key, None)
    cost = value.get("cost", {})
    cost.pop("wall_ns", None)
    extra = cost.get("extra", {})
    for key in (
        "build_ns",
        "eliminate_ns",
        "pair_update_ns",
        "m4ri_scratch_bytes_max",
        "flat_m4ri_matrices",
        "flat_m4ri_order_bytes_max",
    ):
        extra.pop(key, None)
    return value


def metric(item: dict, key: str) -> float:
    return item["process"]["metrics"][key]


def scan_metrics() -> dict:
    rows = []
    for path in sorted(DEV.rglob("*")):
        if not path.is_file() or not (
            path.name == "metrics.json"
            or path.name.endswith(".metrics.json")
            or path.name.endswith(".metrics")
        ):
            continue
        try:
            record = load(path)
        except (OSError, json.JSONDecodeError):
            continue
        metrics = record.get("metrics", record)
        if all(
            key in metrics
            for key in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")
        ):
            rows.append(metrics)
    return {
        "resource_components": len(rows),
        "summed_wall_seconds": math.fsum(row["wall_seconds"] for row in rows),
        "total_core_seconds": math.fsum(row["total_core_seconds"] for row in rows),
        "peak_rss_bytes": max(row["peak_rss_bytes"] for row in rows),
    }


def main() -> None:
    value = load(RESULT)
    checks = 0

    def require(condition: bool, message: str) -> None:
        nonlocal checks
        if not condition:
            raise AssertionError(message)
        checks += 1

    require(value["schema"] == "koblitz_stage174_size_gated_flat_m4ri.v1", "schema")
    require(value["status"] == "complete_same_target_parallel_f4_improvement", "status")
    require(value["koblitz_index_calculus_sota"] is False, "SOTA boundary")
    require(value["fresh_target_holdout_complete"] is False, "fresh-target boundary")
    require(value["licensed_magma_f4_complete"] is False, "Magma boundary")
    require(value["full_cost_gate_passed"] is False, "full-cost boundary")
    require(value["independent_external_reproduction_satisfied"] is False, "external boundary")
    require(value["gates"]["all_seven_gates_passed"] is False, "seven gates")

    source = value["source_revision"]
    require(source["commit"] == COMMIT and source["tree"] == TREE, "source revision")
    repository = HERE.parents[3]
    subprocess.run(["git", "cat-file", "-e", f"{COMMIT}^{{commit}}"], cwd=repository, check=True)
    actual_tree = subprocess.run(
        ["git", "rev-parse", f"{COMMIT}^{{tree}}"],
        cwd=repository,
        text=True,
        capture_output=True,
        check=True,
    ).stdout.strip()
    require(actual_tree == TREE, "source tree")

    artifacts = []

    def collect(node: object) -> None:
        if isinstance(node, dict):
            if set(node) == {"path", "bytes", "sha256"}:
                artifacts.append(node)
            else:
                for child in node.values():
                    collect(child)
        elif isinstance(node, list):
            for child in node:
                collect(child)

    collect(value)
    require(bool(artifacts), "artifact references")
    for receipt in artifacts:
        path = (HERE / receipt["path"]).resolve()
        require(path.is_relative_to(HERE.resolve()), f"artifact escapes stage: {path}")
        require(path.is_file(), f"artifact missing: {path}")
        require(path.stat().st_size == receipt["bytes"], f"artifact size: {path}")
        require(sha256(path) == receipt["sha256"], f"artifact hash: {path}")

    pairs = value["clean_paired_repeats"]
    candidates = pairs["candidates"]
    controls = pairs["controls"]
    require(len(candidates) == len(controls) == 3, "paired run count")
    baseline = stable_report(controls[0]["report"])
    require(
        all(stable_report(item["report"]) == baseline for item in candidates + controls),
        "scientific reports differ",
    )
    for item in candidates:
        report = item["report"]
        require(report["status"] == "unsat" and report["exhaustive"] is True, "candidate terminal")
        require(report["solver_equations_blake3"] == FINGERPRINT, "candidate equations")
        require(report["cost"]["ops"] == 99_199_976_264, "candidate work")
        require(report["cost"]["extra"]["flat_m4ri_matrices"] == 241, "candidate flat count")
        require(report["solver_flat_m4ri_min_words"] == 1_048_576, "candidate threshold")
    for item in controls:
        report = item["report"]
        require(report["status"] == "unsat" and report["exhaustive"] is True, "control terminal")
        require(report["solver_equations_blake3"] == FINGERPRINT, "control equations")
        require(report["cost"]["ops"] == 99_199_976_264, "control work")
        require(report["cost"]["extra"]["flat_m4ri_matrices"] == 0, "control flat count")

    recomputed_pairs = []
    for index, (candidate, control) in enumerate(zip(candidates, controls), 1):
        cm = candidate["process"]["metrics"]
        rm = control["process"]["metrics"]
        cx = candidate["report"]["cost"]["extra"]
        rx = control["report"]["cost"]["extra"]
        recomputed_pairs.append(
            {
                "pair": index,
                "candidate_over_control_wall": cm["wall_seconds"] / rm["wall_seconds"],
                "candidate_over_control_core": cm["total_core_seconds"] / rm["total_core_seconds"],
                "candidate_over_control_rss": cm["peak_rss_bytes"] / rm["peak_rss_bytes"],
                "candidate_over_control_build_ns": cx["build_ns"] / rx["build_ns"],
                "candidate_over_control_eliminate_ns": cx["eliminate_ns"] / rx["eliminate_ns"],
            }
        )
    require(recomputed_pairs == pairs["pair_ratios"], "pair ratios")
    medians = {
        key: statistics.median(item[key] for item in recomputed_pairs)
        for key in recomputed_pairs[0]
        if key != "pair"
    }
    require(medians == pairs["median_pair_ratios"], "pair medians")
    require(medians["candidate_over_control_wall"] < 1, "wall improvement")
    require(medians["candidate_over_control_core"] < 1, "CPU improvement")

    direct = value["same_binary_direct_mitm"]
    direct_medians = {
        key: statistics.median(metric(item, key) for item in direct["runs"])
        for key in ("wall_seconds", "total_core_seconds", "peak_rss_bytes")
    }
    require(direct_medians == direct["median_metrics"], "direct medians")
    require(
        value["comparisons"]["selected_over_direct_mitm_median"]["wall_ratio"] > 1,
        "direct wall boundary",
    )
    require(
        value["comparisons"]["selected_over_direct_mitm_median"]["core_ratio"] > 1,
        "direct CPU boundary",
    )

    single = value["single_core_ablation"]
    require(single["segmented"]["report"]["single_thread_requested"] is True, "single control")
    require(single["hybrid"]["report"]["single_thread_requested"] is True, "single hybrid")
    require(single["hybrid_over_segmented"]["core_ratio"] >= 1, "single-core non-promotion")
    sat = value["validation"]["sat_witness_control"]
    for item in sat.values():
        require(item["report"]["status"] == "sat", "SAT terminal")
        require(item["report"]["source_model_valid"] is True, "SAT model")
        require(item["report"]["source_witness_valid"] is True, "SAT witness")

    incremental = scan_metrics()
    recorded = value["campaign_accounting"]["stage174_incremental_measured_lower_bound"]
    for key in incremental:
        require(incremental[key] == recorded[key], f"incremental accounting {key}")
    previous = load(STAGE173)["campaign_accounting"]["measured_lower_bound"]
    campaign = value["campaign_accounting"]["measured_lower_bound_through_stage174"]
    require(
        campaign["resource_components"]
        == previous["resource_components"] + incremental["resource_components"],
        "campaign components",
    )
    require(
        math.isclose(
            campaign["summed_wall_seconds"],
            previous["summed_wall_seconds"] + incremental["summed_wall_seconds"],
            abs_tol=1e-9,
            rel_tol=0,
        ),
        "campaign wall",
    )
    require(
        math.isclose(
            campaign["total_core_seconds"],
            previous["total_core_seconds"] + incremental["total_core_seconds"],
            abs_tol=1e-9,
            rel_tol=0,
        ),
        "campaign core",
    )

    # External exact-build artifacts are optional on independent runners. If
    # the local custody paths exist, bind them to the committed hashes.
    optional = {
        Path("/Volumes/SSD990/koblitz-native-f4-build23-e51efb21/source.tar"): source[
            "source_archive_sha256"
        ],
        Path("/Volumes/SSD990/koblitz-native-f4-build23-e51efb21/bin/koblitz_pdp_backend"): source[
            "backend_sha256"
        ],
        Path("/Volumes/SSD990/koblitz-native-f4-build23-e51efb21/bin/koblitz_pdp_export"): source[
            "exporter_sha256"
        ],
    }
    checked_optional = 0
    for path, expected in optional.items():
        if path.is_file():
            require(sha256(path) == expected, f"external custody hash: {path}")
            checked_optional += 1

    verification = {
        "schema": "koblitz_stage174_verification.v1",
        "ok": True,
        "checks": checks,
        "artifact_references": len(artifacts),
        "optional_external_artifacts_checked": checked_optional,
        "result_sha256": sha256(RESULT),
        "source_commit": COMMIT,
        "selected_pair_medians": medians,
        "claim_boundary": "same-target native-F4 engineering; not SOTA",
    }
    if os.environ.get("KIC_VERIFY_NO_WRITE") != "1":
        (HERE / "verification.json").write_text(
            json.dumps(verification, indent=2, sort_keys=True) + "\n"
        )
    print(json.dumps(verification, sort_keys=True))


if __name__ == "__main__":
    main()
