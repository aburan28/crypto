#!/usr/bin/env python3
"""Verify the Stage-37 corroborating algebraic-walk run."""

from __future__ import annotations

import argparse
import copy
from datetime import datetime
import json
from pathlib import Path, PurePosixPath
import subprocess
import tarfile
import tempfile
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
import run_koblitz_stage35_algebraic_walk as runner
import verify_koblitz_stage35_result as primary


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
ARCHIVE = STAGE / "koblitz-stage37-algebraic-walk-repeat-34694400320.tar.gz"
ARCHIVE_SHA = STAGE / "koblitz-stage37-algebraic-walk-repeat-34694400320.tar.gz.sha256"
WORKFLOW = STAGE / "stage-37-repeat-workflow-34694400320.json"
ARTIFACTS = STAGE / "stage-37-repeat-artifacts-34694400320.json"
DEFAULT_OUTPUT = STAGE / "stage-37-algebraic-walk-repeat-result-20260912"

RUN_ID = 34694400320
RUN_COMMIT = "467b6b652c39730c6069442613e2397a6f173c48"
ARTIFACT_ID = 10297933814
ARTIFACT_DIGEST = "sha256:b1eec6f64f5dd70ee01a512892519f237397449a69c6970b1ba6dc96b7fbe744"
ARCHIVE_ROOT = "koblitz-stage35-algebraic-walk-34694400320"
SCHEMA = "koblitz_stage37_repeat_result.v1"
SEAL_SCHEMA = "koblitz_stage37_repeat_result_seal.v1"

EXPECTED_SOURCE_DELTA = [
    "RESEARCH_AUTOLAB_LOG.md",
    "RESEARCH_QUASI_SUBFIELD.md",
    "examples/quasi_subfield_census.rs",
    "examples/quasi_subfield_reach.rs",
    "src/cryptanalysis/mod.rs",
    "src/cryptanalysis/quasi_subfield.rs",
]
RELEVANT_PATHS = [
    "Cargo.toml",
    "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-20-rust-build/Cargo.lock",
    "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-35-n41-algebraic-walk-params.json",
    "scripts/process_meter.py",
    "scripts/run_koblitz_stage33_n41_unknown_scalar.py",
    "scripts/run_koblitz_stage35_algebraic_walk.py",
    "src/binary_ecc",
    "src/bin/ic.rs",
    "src/bin/ic",
    "src/cryptanalysis/koblitz_factor_base_search.rs",
    "src/cryptanalysis/koblitz_fast.rs",
    "src/cryptanalysis/koblitz_groebner.rs",
    "src/cryptanalysis/koblitz_index_calculus.rs",
    "src/cryptanalysis/koblitz_sparse_la.rs",
]


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def timestamp(value: str) -> datetime:
    require(isinstance(value, str) and value.endswith("Z"), "workflow timestamp is invalid")
    return datetime.fromisoformat(value[:-1] + "+00:00")


def git(*arguments: str) -> str:
    return subprocess.run(["git", *arguments], cwd=REPO, text=True, capture_output=True, check=True).stdout.strip()


def source_equivalence() -> dict[str, Any]:
    delta = git("diff", "--name-only", primary.RUN_COMMIT, RUN_COMMIT).splitlines()
    require(delta == EXPECTED_SOURCE_DELTA, "Stage-37 source delta changed")
    objects = {}
    for path in RELEVANT_PATHS:
        left = git("rev-parse", f"{primary.RUN_COMMIT}:{path}")
        right = git("rev-parse", f"{RUN_COMMIT}:{path}")
        require(left == right, f"Stage-37 relevant source changed: {path}")
        objects[path] = left
    return {
        "classification": "koblitz_path_equivalent_not_whole_tree_identical",
        "primary_commit": primary.RUN_COMMIT,
        "repeat_commit": RUN_COMMIT,
        "changed_paths": delta,
        "changed_paths_relevant_to_koblitz_run": False,
        "relevant_git_objects": objects,
    }


def archive_identity() -> dict[str, Any]:
    fields = ARCHIVE_SHA.read_text().split()
    relative = str(ARCHIVE.relative_to(REPO))
    require(len(fields) == 2 and fields[1] == relative, "Stage-37 archive sidecar is invalid")
    digest = custody.sha256_file(ARCHIVE, "Stage-37 archive")
    require(fields[0] == digest, "Stage-37 archive hash changed")
    return {"path": relative, "bytes": ARCHIVE.stat().st_size, "sha256": digest}


def extract(destination: Path) -> Path:
    with tarfile.open(ARCHIVE, "r:gz") as source:
        members = source.getmembers()
        roots: set[str] = set()
        for member in members:
            path = PurePosixPath(member.name)
            require(not path.is_absolute() and ".." not in path.parts, f"unsafe archive path: {member.name}")
            require(member.isfile() or member.isdir(), f"unsupported archive member: {member.name}")
            if path.parts:
                roots.add(path.parts[0])
        require(roots == {ARCHIVE_ROOT}, "Stage-37 archive root changed")
        source.extractall(destination, filter="data")
    return destination / ARCHIVE_ROOT


def workflow_accounting() -> dict[str, Any]:
    value = load(WORKFLOW, "Stage-37 workflow metadata")
    require(value.get("databaseId") == RUN_ID and value.get("headSha") == RUN_COMMIT, "Stage-37 workflow identity changed")
    require(value.get("status") == "completed" and value.get("conclusion") == "success", "Stage-37 workflow is incomplete")
    jobs = value.get("jobs")
    require(isinstance(jobs, list) and {row.get("name") for row in jobs} == {"validate-control", "algebraic-walk"}, "Stage-37 workflow jobs changed")
    for row in jobs:
        require(row.get("status") == "completed" and row.get("conclusion") == "success", f"Stage-37 job failed: {row.get('name')}")
    production = next(row for row in jobs if row["name"] == "algebraic-walk")
    return {
        "run_id": RUN_ID,
        "head_sha": RUN_COMMIT,
        "url": value.get("url"),
        "workflow_wall_seconds": (timestamp(value["updatedAt"]) - timestamp(value["createdAt"])).total_seconds(),
        "production_job_wall_seconds": (timestamp(production["completedAt"]) - timestamp(production["startedAt"])).total_seconds(),
        "production_job_url": production.get("url"),
    }


def artifact_accounting() -> dict[str, Any]:
    value = load(ARTIFACTS, "Stage-37 artifact metadata")
    rows = value.get("artifacts")
    require(value.get("total_count") == 1 and isinstance(rows, list) and len(rows) == 1, "Stage-37 artifact count changed")
    row = rows[0]
    require(row.get("id") == ARTIFACT_ID and row.get("digest") == ARTIFACT_DIGEST and row.get("expired") is False, "Stage-37 artifact identity changed")
    require(row.get("name") == f"koblitz-stage35-algebraic-walk-{RUN_ID}", "Stage-37 artifact name changed")
    source = row.get("workflow_run", {})
    require(source.get("id") == RUN_ID and source.get("head_sha") == RUN_COMMIT, "Stage-37 artifact source changed")
    return {key: row[key] for key in ("id", "name", "size_in_bytes", "digest", "expires_at")}


def strip_solution_times(items: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = copy.deepcopy(items)
    for row in out:
        row.pop("elapsed_seconds", None)
    return out


def strip_baseline_times(value: dict[str, Any]) -> dict[str, Any]:
    out = copy.deepcopy(value)
    out.pop("elapsed_seconds", None)
    out["ic"].pop("descent_seconds_total", None)
    out["ic"].pop("descent_seconds_per_target", None)
    out["ic"].pop("precompute_seconds", None)
    out["ic"].pop("amortised_seconds_per_target", None)
    out["rho"].pop("seconds_total", None)
    out["rho"].pop("seconds_per_target", None)
    out.pop("ratio", None)
    out.pop("verdict", None)
    for row in out["targets_detail"]:
        row.pop("seconds", None)
    return out


def mathematical_equivalence(primary_root: Path, repeat_root: Path) -> dict[str, Any]:
    p_run = load(primary_root / "stage35-run/result.json", "primary Stage-35 run")
    r_run = load(repeat_root / "stage35-run/result.json", "repeat Stage-37 run")
    p_work, r_work = p_run["workflow_result"], r_run["workflow_result"]
    require(p_work["params_digest"] == r_work["params_digest"], "Stage-37 parameter digest changed")
    require(p_work["factor_base"] == r_work["factor_base"], "Stage-37 factor base changed")
    p_items, r_items = p_work["solutions"]["items"], r_work["solutions"]["items"]
    require(strip_solution_times(p_items) == strip_solution_times(r_items), "Stage-37 solutions or walked-probe counts changed")
    p_base = next(row["vs_rho"] for row in p_work["stages"] if row["stage"] == "baseline")
    r_base = next(row["vs_rho"] for row in r_work["stages"] if row["stage"] == "baseline")
    require(strip_baseline_times(p_base) == strip_baseline_times(r_base), "Stage-37 rho trajectory or charge counts changed")
    require(
        (primary_root / "stage35-run/workflow/factor_base.json").read_bytes()
        == (repeat_root / "stage35-run/workflow/factor_base.json").read_bytes(),
        "Stage-37 factor-base artifact changed",
    )
    require(
        (primary_root / "stage35-run/workflow/logs.json").read_bytes()
        == (repeat_root / "stage35-run/workflow/logs.json").read_bytes(),
        "Stage-37 log database changed",
    )
    return {
        "params_digest": p_work["params_digest"],
        "factor_base_exact": True,
        "factor_base_logs_exact": True,
        "public_targets_recovered_logs_and_walked_probe_counts_exact": True,
        "rho_recovered_logs_iterations_restarts_and_group_additions_exact": True,
        "descent_trials": p_base["ic"]["descent_trials_total"],
        "rho_iterations": p_base["rho"]["iterations_total"],
    }


def result() -> dict[str, Any]:
    source = source_equivalence()
    with tempfile.TemporaryDirectory() as directory:
        temp = Path(directory)
        primary_root = primary.extract(temp / "primary")
        repeat_root = extract(temp / "repeat")
        primary_result = runner.verify(primary_root / "stage35-build", primary_root / "stage35-run")
        repeat_result = runner.verify(repeat_root / "stage35-build", repeat_root / "stage35-run")
        embedded = load(repeat_root / "stage35-verification.json", "embedded Stage-37 verification")
        require(embedded == repeat_result, "Stage-37 embedded verification changed")
        mathematics = mathematical_equivalence(primary_root, repeat_root)
    require(primary_result["ic_over_rho_online_wall_ratio"] < 1.0 and repeat_result["ic_over_rho_online_wall_ratio"] < 1.0, "one Stage-35 run lost the online crossover")
    speedups = [1.0 / primary_result["ic_over_rho_online_wall_ratio"], 1.0 / repeat_result["ic_over_rho_online_wall_ratio"]]
    require(primary_result["full_cost_gate_passed"] is False and repeat_result["full_cost_gate_passed"] is False, "Stage-37 full-cost claim widened")
    return {
        "schema": SCHEMA,
        "status": "two_hosted_runs_koblitz_path_equivalent_online_crossover",
        "source_equivalence": source,
        "archive": archive_identity(),
        "artifact": artifact_accounting(),
        "workflow": workflow_accounting(),
        "mathematical_equivalence": mathematics,
        "primary": primary_result,
        "repeat": repeat_result,
        "online_ic_speedup_over_rho": {
            "primary": speedups[0],
            "repeat": speedups[1],
            "minimum": min(speedups),
            "maximum": max(speedups),
        },
        "claim_boundary": "two project-authored hosted runs with Koblitz-path-equivalent source; finite post-precomputation online crossover only",
        "whole_tree_exact_replicate": False,
        "unaffiliated_reproduction": False,
        "amortised_crossover": False,
        "whole_process_crossover": False,
        "full_cost_gate_passed": False,
        "licensed_magma_complete": False,
        "independent_external_reproduction_satisfied": False,
        "koblitz_index_calculus_sota": False,
    }


def freeze(output: Path) -> dict[str, Any]:
    require(not output.exists() and not output.is_symlink(), "Stage-37 output must be new")
    output.mkdir(parents=False)
    custody.write_json_new(output / "verification.json", result())
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    payload = {"schema": SEAL_SCHEMA, "status": "repeat_frozen", "inventory": inventory, "inventory_sha256": custody.canonical_sha256(inventory)}
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(output / "result-seal.json", seal)
    return seal


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    committed = load(output / "verification.json", "committed Stage-37 result")
    require(committed == result(), "committed Stage-37 result differs from source evidence")
    seal = load(output / "result-seal.json", "Stage-37 result seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == custody.canonical_sha256(payload), "Stage-37 seal is invalid")
    inventory = custody.all_regular_inventory(output, {"result-seal.json"})
    require(inventory == seal.get("inventory") and custody.canonical_sha256(inventory) == seal.get("inventory_sha256"), "Stage-37 inventory changed")
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
        value = freeze(args.output.resolve()) if args.command == "compose" else verify(args.output.resolve())
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, VerificationError, runner.Stage35Error, custody.PhaseBError) as error:
        raise SystemExit(f"stage37-repeat: {error}")


if __name__ == "__main__":
    main()
