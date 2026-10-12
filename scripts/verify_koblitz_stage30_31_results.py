#!/usr/bin/env python3
"""Verify retained Stage-30 search and Stage-31 unknown-scalar results."""

from __future__ import annotations

from datetime import datetime
import json
import math
from pathlib import Path
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage31_n31_unknown_scalar as stage31_runner


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
ROOT30 = STAGE / "stage-30-n31-factor-search-result-20260911"
ROOT31 = STAGE / "stage-31-n31-unknown-scalar-result-20260911"
RUN30 = 34639879925
COMMIT30 = "c6f93c9494197ba545eae4c4dbdb6f3f4a2e0b09"
RUN31 = 34642615206
COMMIT31 = "bfebb2c743885371ebd5d9364936ccaea61ff075"


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"JSON object required: {path}")
    return value


def identity(path: Path) -> dict[str, Any]:
    return {"bytes": path.stat().st_size, "sha256": custody.sha256_file(path, "retained result")}


def verify_seal(root: Path, number: int, run: int, commit: str) -> dict[str, Any]:
    seal = load(root / "result-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(
        seal.get("schema") == f"koblitz_stage{number}_result_seal.v1"
        and seal.get("status") == "result_frozen"
        and seal.get("workflow_run_id") == run
        and seal.get("workflow_commit") == commit
        and claimed == custody.canonical_sha256(payload),
        f"Stage-{number} result seal is invalid",
    )
    inventory = custody.all_regular_inventory(root, {"result-seal.json"})
    require(inventory == seal.get("inventory") and custody.canonical_sha256(inventory) == seal.get("inventory_sha256"), f"Stage-{number} inventory changed")
    return seal


def verify_workflow(root: Path, run: int, commit: str, artifact_name: str, digest: str) -> dict[str, Any]:
    workflow = load(root / "workflow.json")
    require(
        workflow.get("databaseId") == run
        and workflow.get("headSha") == commit
        and workflow.get("status") == "completed"
        and workflow.get("conclusion") == "success",
        "workflow metadata is incomplete",
    )
    artifacts = load(root / "artifacts.json")
    rows = artifacts.get("artifacts")
    require(artifacts.get("total_count") == 1 and isinstance(rows, list) and len(rows) == 1, "artifact metadata count changed")
    artifact = rows[0]
    require(
        artifact.get("name") == artifact_name
        and artifact.get("digest") == digest
        and artifact.get("expired") is False
        and isinstance(artifact.get("id"), int)
        and isinstance(artifact.get("size_in_bytes"), int)
        and artifact["size_in_bytes"] > 0,
        "artifact metadata changed",
    )
    created = datetime.fromisoformat(workflow["createdAt"].replace("Z", "+00:00"))
    updated = datetime.fromisoformat(workflow["updatedAt"].replace("Z", "+00:00"))
    return {"url": workflow["url"], "wall_seconds": (updated - created).total_seconds(), "artifact": {key: artifact[key] for key in ("id", "name", "size_in_bytes", "digest", "expires_at")}}


def parse_build_time(path: Path) -> dict[str, Any]:
    fields = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith("User time (seconds):"):
            fields["user"] = float(line.rsplit(":", 1)[1])
        elif line.startswith("System time (seconds):"):
            fields["system"] = float(line.rsplit(":", 1)[1])
        elif line.startswith("Elapsed (wall clock) time"):
            fields["wall"] = line.split("): ", 1)[1]
        elif line.startswith("Maximum resident set size (kbytes):"):
            fields["rss"] = int(line.rsplit(":", 1)[1]) * 1024
        elif line.startswith("Exit status:"):
            fields["exit"] = int(line.rsplit(":", 1)[1])
    require(set(fields) == {"user", "system", "wall", "rss", "exit"} and fields["exit"] == 0, "build-time receipt is incomplete")
    pieces = [float(value) for value in fields["wall"].split(":")]
    wall = pieces[0] * 60 + pieces[1] if len(pieces) == 2 else pieces[0] * 3600 + pieces[1] * 60 + pieces[2]
    return {"total_core_seconds": fields["user"] + fields["system"], "wall_seconds": wall, "peak_rss_bytes": fields["rss"], "receipt": identity(path)}


def verify_stage30() -> dict[str, Any]:
    seal = verify_seal(ROOT30, 30, RUN30, COMMIT30)
    workflow = verify_workflow(ROOT30, RUN30, COMMIT30, f"koblitz-stage30-n31-search-{RUN30}", "sha256:0310cb60af5e2ce5847a555f37ba80ae02565b08ffba0ea8c130760dfbcc4d17")
    wrapper = load(ROOT30 / "stage30-search/wrapper-result.json")
    search = load(ROOT30 / "stage30-search/search.json")
    recipe_path = ROOT30 / "stage30-search/factor-base.json"
    recipe = stage31_runner.validate_factor_base(recipe_path)
    require(wrapper.get("schema") == "koblitz_stage30_n31_factor_base_search.v1" and wrapper.get("status") == "complete_validated_factor_base", "Stage-30 wrapper is incomplete")
    require(wrapper.get("information_boundary") == {"factor_base_discrete_log_labels_used": False, "future_unknown_scalar_target_available": False, "relation_yield_and_holdouts": "public search samples and separate generated validation controls only", "target_subgroup_enumerated": False}, "Stage-30 information boundary changed")
    process = wrapper.get("search_process", {})
    require(process.get("returncode") == 0 and process.get("timed_out") is False and process.get("orphan_group_terminated") is False, "Stage-30 process did not complete cleanly")
    require(wrapper.get("artifacts", {}).get("factor_base", {}).get("sha256") == identity(recipe_path)["sha256"], "Stage-30 recipe identity changed")
    require(wrapper.get("artifacts", {}).get("search", {}).get("sha256") == identity(ROOT30 / "stage30-search/search.json")["sha256"], "Stage-30 search identity changed")
    require(search.get("schema_version") == 1 and search.get("operation") == "search" and search.get("status") == "complete", "Stage-30 search status changed")
    require(search.get("degree") == 31 and search.get("curve_a") == 0 and search.get("summands") == 3 and search.get("targets") == 256 and search.get("candidate_count") == 56, "Stage-30 search dimensions changed")
    require(search.get("selected") == recipe and search.get("best_census") == recipe["spec"], "Stage-30 selected recipe changed")
    selected = search.get("selected_summary", {})
    require(selected.get("spec") == recipe["spec"] and selected.get("points") == 3971 and selected.get("projected_columns") == 32 and selected.get("coverage") == 0.796875 and selected.get("covered") == 204, "Stage-30 selected summary changed")
    validation = search.get("validation", {})
    runs = validation.get("runs")
    require(validation.get("holdout_samples") == 2 and validation.get("validated_top") == 1 and isinstance(runs, list) and len(runs) == 1 and runs[0].get("eligible") is True and len(runs[0].get("holdout", [])) == 2, "Stage-30 holdout validation changed")
    require(all(row.get("status") == "complete" and row.get("result", {}).get("verified") is True for row in runs[0]["holdout"]), "Stage-30 holdout failed")
    outer = wrapper.get("outer_resources", {})
    require(all(isinstance(outer.get(name), (int, float)) and math.isfinite(outer[name]) and outer[name] > 0 for name in ("wall_seconds", "total_core_seconds", "sampled_peak_process_tree_rss_bytes")), "Stage-30 outer resources changed")
    return {"seal_inventory_sha256": seal["inventory_sha256"], "workflow": workflow, "recipe": {**identity(recipe_path), "value": recipe}, "selected_summary": selected, "validation": {"holdouts": 2, "median_process_seconds": runs[0]["median_process_seconds"]}, "search_outer_resources": outer, "search_process_resources": process["metrics"], "build_resources": parse_build_time(ROOT30 / "stage30-build/build.time")}


def verify_stage31(stage30: dict[str, Any]) -> dict[str, Any]:
    seal = verify_seal(ROOT31, 31, RUN31, COMMIT31)
    workflow_meta = verify_workflow(ROOT31, RUN31, COMMIT31, f"koblitz-stage31-n31-unknown-scalar-{RUN31}", "sha256:4e20284b44c32c52a4e1ca2d36ac0bd98f8de311b0973517ce981c8394b89a41")
    result = load(ROOT31 / "stage31-run/result.json")
    workflow = load(ROOT31 / "stage31-run/workflow.stdout")
    require(result.get("schema") == "koblitz_stage31_n31_unknown_scalar_result.v1" and result.get("status") == "complete_public_unknown_scalar_panel", "Stage-31 result is incomplete")
    require(result.get("workflow_result") == workflow, "Stage-31 wrapper and workflow output differ")
    require(result.get("factor_base_recipe", {}).get("sha256") == stage30["recipe"]["sha256"], "Stage-31 used a different factor base")
    require(result.get("target_scalars_constructed_or_supplied") is False and result.get("factor_base_logs_known_by_construction") is False, "Stage-31 information boundary changed")
    stage31_runner.validate_workflow(workflow)
    params = load(ROOT31 / "stage31-run/params.json")
    require(params.get("factor_base") == {"mode": "spec", "spec": stage30["recipe"]["value"]["spec"]}, "Stage-31 params changed the recipe")
    require(params.get("targets") == [{"public_hash_seed": seed} for seed in stage31_runner.PUBLIC_TARGET_SEEDS], "Stage-31 public targets changed")
    logs = load(ROOT31 / "stage31-run/workflow/logs.json")
    require(logs.get("spec") == stage30["recipe"]["value"]["spec"] and len(logs.get("columns", [])) == 32, "Stage-31 factor-base log database changed")
    require(all(isinstance(row.get("log"), str) and row.get("x") and row.get("y") for row in logs["columns"]), "Stage-31 log column is incomplete")
    solutions = load(ROOT31 / "stage31-run/workflow/solutions.json")
    require(len(solutions.get("solutions", [])) == 5 and all(row.get("verified") is True and row.get("expected") == "not_constructed" for row in solutions["solutions"]), "Stage-31 solutions changed")
    baseline = load(ROOT31 / "stage31-run/workflow/baseline.json")
    require(baseline.get("claim_boundary") == "public_hash_unknown_scalar" and baseline.get("ic", {}).get("verified") == 5 and baseline.get("rho", {}).get("verified") == 5, "Stage-31 rho comparison changed")
    recovered = [row["recovered"] for row in solutions["solutions"]]
    require(recovered == [row.get("recovered") for row in baseline.get("targets_detail", [])], "Stage-31 IC and rho recovered different scalars")
    process = result.get("process_receipt", {})
    require(process.get("returncode") == 0 and process.get("timed_out") is False and process.get("orphan_group_terminated") is False, "Stage-31 process did not complete cleanly")
    outer = result.get("outer_resources", {})
    require(all(isinstance(outer.get(name), (int, float)) and math.isfinite(outer[name]) and outer[name] > 0 for name in ("wall_seconds", "total_core_seconds", "sampled_peak_process_tree_rss_bytes")), "Stage-31 outer resources changed")
    return {"seal_inventory_sha256": seal["inventory_sha256"], "workflow": workflow_meta, "targets": solutions["solutions"], "factor_base_logs": len(logs["columns"]), "collection": next(stage for stage in workflow["stages"] if stage["stage"] == "collect"), "linear_algebra": next(stage for stage in workflow["stages"] if stage["stage"] == "logs")["linear_algebra"], "baseline": baseline, "workflow_outer_resources": outer, "workflow_process_resources": process["metrics"], "build_resources": parse_build_time(ROOT31 / "stage31-build/build.time")}


def verify() -> dict[str, Any]:
    stage30 = verify_stage30()
    stage31 = verify_stage31(stage30)
    rho = stage31["baseline"]["rho"]
    ic = stage31["baseline"]["ic"]
    return {
        "schema": "koblitz_stage30_31_result_verification.v1",
        "status": "validated_factor_base_and_five_public_unknown_scalars_verified",
        "stage30": stage30,
        "stage31": stage31,
        "comparisons": {
            "online_ic_over_rho_wall": ic["descent_seconds_per_target"] / rho["seconds_per_target"],
            "amortised_ic_over_rho_wall": ic["amortised_seconds_per_target"] / rho["seconds_per_target"],
            "stage31_outer_over_rho_total_wall": stage31["workflow_outer_resources"]["wall_seconds"] / rho["seconds_total"],
        },
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }


if __name__ == "__main__":
    try:
        print(json.dumps(verify(), indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, VerificationError, stage31_runner.Stage31Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage30-31-results: {error}")
