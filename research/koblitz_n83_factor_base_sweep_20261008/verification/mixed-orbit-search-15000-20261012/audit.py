#!/usr/bin/env python3
"""Independent receipt audit for the first retained N83 mixed-orbit search."""

import hashlib
import json
from pathlib import Path
import subprocess


EVIDENCE = Path(__file__).resolve().parent / "evidence-01"
ROOT = EVIDENCE.parents[4]
COMMIT = "f0a1e55f7cb6d965b0791ef68cb636ec8b8a3a1a"
IMAGE = "sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31"
PRIOR_AUDIT_SHA = "76a8bf74e46feb0ef4fbe650395f39733f4ff95d31c35f8829cb506a9ad04c7d"
SMALL_SHA = "eb342b574e244e1e67e8e78c805f33ee88edcdb4f3a14d4114e6f413beede380"
ENVELOPE_SHA = "aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a"
BINARY_SHA = "264df6474b2acc4865091a3eef8da97424c694affa0dc994971de7d01df4e346"


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read(path: str) -> dict:
    return json.loads((EVIDENCE / path).read_text())


def require(condition: bool, label: str) -> None:
    if not condition:
        raise ValueError(label)


def option(command: list[str], name: str) -> str:
    return command[command.index(name) + 1]


def audit_stage(stage: dict, name: str, wall: int, cpus: str, memory: str,
                location: Path) -> None:
    require(stage["name"] == name and stage["exit_code"] == 0
            and stage["timed_out"] is False and stage["launch_error"] is None
            and stage["cleanup_verified"] is True and stage["cleanup_error"] is None,
            f"{name} outer disposition")
    state = stage["container_state"]
    require(state["ExitCode"] == 0 and state["OOMKilled"] is False
            and state["Status"] == "exited" and state["Running"] is False
            and stage["wall_elapsed_seconds_L0"] <= wall + 30,
            f"{name} cgroup exit or wall cap")
    command = stage["command"]
    require(option(command, "--pull") == "never"
            and option(command, "--platform") == "linux/arm64"
            and option(command, "--network") == "none"
            and option(command, "--cpus") == cpus
            and option(command, "--memory") == memory
            and option(command, "--memory-swap") == memory
            and IMAGE in command, f"{name} isolation")
    if name != "build":
        require("--read-only" in command and option(command, "--pids-limit") == "64"
                and option(command, "--security-opt") == "no-new-privileges"
                and any("dst=/worker,readonly" in arg for arg in command)
                and any("dst=/small.jsonl.gz,readonly" in arg for arg in command)
                and any("dst=/envelope.jsonl.gz,readonly" in arg for arg in command),
                f"{name} mounts")
    for suffix in ("stdout", "stderr"):
        require(sha(location / f"{name}.{suffix}") == stage[f"{suffix}_sha256"],
                f"{name} {suffix} hash")


def check() -> dict:
    auth = read("authorization.json")
    build = read("build/build.json")
    config = read("run/config.json")
    outer = read("run/outer.json")
    budget = read("run/budget.json")
    worker_config = read("run/run/config.json")
    result = read("run/run/result.json")
    require(read("build-launch.stdout") == build
            and read("run-launch.stdout") == outer
            and (EVIDENCE / "build-launch.stderr").stat().st_size == 0
            and (EVIDENCE / "run-launch.stderr").stat().st_size == 0,
            "launch receipts and stderr")
    event_lines = (EVIDENCE / "run/run/events.jsonl").read_text().splitlines()
    require(len(event_lines) == 1, "exactly one retained trial")
    event = json.loads(event_lines[0])
    auth_sha = sha(EVIDENCE / "authorization.json")
    require(auth == {
        "schema": "n83.additional-local-budget/v2",
        "study": "koblitz_n83_factor_base_sweep_20261008",
        "purpose": "mixed_orbit_two_residual_first_gate",
        "prior_total_charged_seconds": "12600",
        "prior_import_audit_sha256": PRIOR_AUDIT_SHA,
        "additional_authorized_seconds": "2400",
        "build_dir": build["binary"].removesuffix("/worker"),
        "output_dir": str(Path(build["binary"]).parent.parent / "run"),
        "user_approval_reference": "2026-10-12 user yes to 2400 additional local seconds, annotation 1",
    }, "exact authorization")
    require(auth["output_dir"] == str(Path(auth["build_dir"]).parent / "run"),
            "authorization output path")
    require(sha(EVIDENCE / "build/worker") == BINARY_SHA, "retained Linux worker hash")
    require((EVIDENCE / "build/worker").read_bytes()[:4] == b"\x7fELF", "Linux worker ELF")
    require(build["schema"] == "n83.mixed-orbit-linux-build/v1"
            and build["status"] == "PASS_build_only" and build["source_commit"] == COMMIT
            and build["authorization_sha256"] == auth_sha
            and build["binary_sha256"] == BINARY_SHA
            and build["container_image_id"] == IMAGE
            and build["conservative_charge_seconds"] == 600
            and 0 <= build["elapsed_seconds_L0"] <= 600
            and build["retained_n83_search_executed"] is False,
            "frozen build receipt")
    audit_stage(build["stage"], "build", 450, "4", "8192m", EVIDENCE / "build")
    require(build["stage"]["command"][-2:] ==
            ["--example", "koblitz_n83_mixed_orbit_worker"], "build target")
    for relative, digest in build["source_sha256"].items():
        blob = subprocess.run(["git", "show", f"{COMMIT}:{relative}"], cwd=ROOT,
                              check=True, capture_output=True).stdout
        require(hashlib.sha256(blob).hexdigest() == digest, f"source blob {relative}")
    inputs = config["input_sha256"]
    objects = read("objects.json")
    require(objects["schema"] == "n83.retained-object-downloads/v1"
            and [(row["label"], row["sha256"], row["size_bytes"])
                 for row in objects["objects"]] == [
                     ("small", SMALL_SHA, 5956114),
                     ("envelope", ENVELOPE_SHA, 10313647),
                 ], "downloaded object receipts")
    object_names = {
        "small": "76e1cbe4ea5a9d26e5b24a2d27f9401271dc471846529190a0416d8e09837a92",
        "envelope": "82990cda5d57b1c4c4b635673a881ad233c24f707e7ac644a11478aeddddf851",
    }
    for row in objects["objects"]:
        require(row["s3_uri"] == (
            "s3://crypto-autoresearcher/factor-bases/icv1/etc/"
            "koblitz_n83_factor_base_sweep_20261008/v2-size-frontier/a0/objects/"
            + object_names[row["label"]] + ".jsonl.gz")
            and row["head_object_2026_10_12"]["content_length"] == row["size_bytes"],
            "S3 object presence or path")
        local = Path(auth["build_dir"]).parent / "objects" / f"{row['label']}.jsonl.gz"
        if local.exists():
            require(sha(local) == row["sha256"] and local.stat().st_size == row["size_bytes"],
                    f"{row['label']} current downloaded bytes")
    expected_inputs = {
        "authorization.json": auth_sha,
        "build.json": sha(EVIDENCE / "build/build.json"),
        "worker": BINARY_SHA,
        "small.jsonl.gz": SMALL_SHA,
        "envelope.jsonl.gz": ENVELOPE_SHA,
    }
    for path, digest in inputs.items():
        match = next((value for suffix, value in expected_inputs.items()
                      if path.endswith("/" + suffix)), None)
        if match is None:
            relative = str(Path(path).relative_to(ROOT))
            match = build["source_sha256"].get(relative)
        require(digest == match, f"frozen input {path}")
    require(len(inputs) == len(expected_inputs) + len(build["source_sha256"]),
            "frozen input set")
    require(config["schema"] == "n83.mixed-orbit-supervisor-config/v1"
            and config["source_commit"] == COMMIT
            and config["authorization_sha256"] == auth_sha
            and config["container_image_id"] == IMAGE
            and config["build_cap_seconds"] == 600
            and config["verify_cap_seconds"] == 360
            and config["search_cap_seconds"] == 1200
            and config["cleanup_reserve_seconds"] == 240
            and config["conservative_total_charge_seconds"] == 2400
            and config["max_trials"] == 1
            and config["conflicts_per_trial"] == 5000
            and config["max_models_per_trial"] == 32
            and config["small_summands"] == 3
            and config["residual_summands"] == 2
            and config["total_index_calculus_runtime_ms"] is None,
            "run configuration")
    require(outer["schema"] == "n83.mixed-orbit-supervisor-outer/v1"
            and outer["source_commit"] == COMMIT
            and outer["config_sha256"] == sha(EVIDENCE / "run/config.json")
            and outer["authorization_sha256"] == auth_sha
            and outer["status"] == "UNKNOWN_search_cap"
            and outer["import_passed"] is True
            and outer["event_replay_error"] is None
            and outer["frozen_input_changed"] is False
            and 0 <= outer["elapsed_supervisor_seconds_L0"] <= 1800
            and outer["total_index_calculus_runtime_ms"] is None
            and outer["selected_best_total_runtime"] is None,
            "outer receipt")
    require([stage["name"] for stage in outer["stages"]] == ["verify", "search"],
            "launched run stages")
    audit_stage(outer["stages"][0], "verify", 360, "1", "4096m", EVIDENCE / "run")
    audit_stage(outer["stages"][1], "search", 1200, "1", "4096m", EVIDENCE / "run")
    verify_args = outer["stages"][0]["command"]
    search_args = outer["stages"][1]["command"]
    require(verify_args[verify_args.index("/worker") + 1:] ==
            ["verify", "/small-manifest.json", "/small.jsonl.gz",
             "/envelope-manifest.json", "/envelope.jsonl.gz"], "verify invocation")
    require(search_args[search_args.index("/worker") + 1:] ==
            ["run", "/small-manifest.json", "/small.jsonl.gz",
             "/envelope-manifest.json", "/envelope.jsonl.gz", "/corpus.json",
             "1150", "1", "5000", "32", "3", "2", "/out/run"], "search invocation")
    require(read("run/verify.stdout") == {
        "schema": "n83.mixed-orbit-import/v1", "status": "PASS",
        "curve_slug": "icv1-f2m83-tm6151469093347-debefd74",
        "small_object_blake3": "76e1cbe4ea5a9d26e5b24a2d27f9401271dc471846529190a0416d8e09837a92",
        "envelope_object_blake3": "82990cda5d57b1c4c4b635673a881ad233c24f707e7ac644a11478aeddddf851",
        "small_points": 196212, "small_signed_orbits": 1182,
        "residual_points": 143756, "residual_signed_orbits": 866,
        "solver_search_executed": False,
    }, "import receipt")
    require(worker_config["schema"] == "n83.mixed-orbit-partial-config/v1"
            and worker_config["budget_seconds"] == 1150
            and worker_config["max_trials"] == 1
            and worker_config["conflict_budget_per_trial"] == 5000
            and worker_config["max_models_per_trial"] == 32
            and worker_config["small_summands"] == 3
            and worker_config["residual_summands"] == 2
            and worker_config["order_first_pair"] is True
            and worker_config["max_variables"] == 150000
            and worker_config["max_domain_clauses"] == 12000000
            and worker_config["memory_cgroup_limit_bytes"] == 4 * 1024**3
            and worker_config["memory_cgroup_swap_limit_bytes"] == 0,
            "worker configuration")
    require(result["schema"] == "n83.mixed-orbit-partial-result/v1"
            and result["status"] == "UNKNOWN_search_cap"
            and result["trials_attempted"] == result["solver_exhausted_trials"] == 1
            and result["solver_calls"] == 1
            and result["solver_conflicts"] == 5000
            and all(result[key] == 0 for key in (
                "solver_models", "invalid_models", "group_rejected_models",
                "graph_accepted_candidates", "graph_pivot_rows", "graph_full_from_partials",
                "graph_capped_rows", "completed_rows"))
            and result["matrix_rank"] is None
            and result["total_index_calculus_runtime_ms"] is None
            and result["selected_best_total_runtime"] is None
            and 0 < result["memory_cgroup_peak_bytes"] <= 4 * 1024**3,
            "worker search disposition")
    require(event["trial"] == 0
            and event["a"] == "1833078656608906180176350"
            and event["b"] == "163089659398503559150141"
            and event["candidate"] is None and event["graph"] is None
            and event["search"]["solver_calls"] == 1
            and event["search"]["conflicts"] == 5000
            and event["search"]["models"] == 0
            and event["search"]["exhausted"] is True
            and event["search"]["refuted"] is False
            and event["search"]["peak_variables"] <= 150000
            and event["search"]["domain_clauses"] <= 12000000
            and event["elapsed_process_wall_ms"] <= 1200 * 1000,
            "event ledger")
    require(budget["schema"] == "n83.mixed-orbit-supervisor-budget/v2"
            and budget["authorization_sha256"] == auth_sha
            and budget["prior_total_charged_seconds"] == "12600"
            and budget["prior_import_audit_sha256"] == PRIOR_AUDIT_SHA
            and budget["additional_authorized_seconds"] == "2400"
            and budget["build_conservative_charge_seconds"] == 600
            and budget["run_conservative_charge_seconds"] == 1800
            and budget["conservative_charge_seconds"] == 2400
            and budget["remaining_additional_seconds"] == "0"
            and budget["round_status"] == outer["status"]
            and budget["selected_best_total_runtime"] is None,
            "conservative budget")
    return {
        "schema": "n83.mixed-orbit-search-15000-audit/v1",
        "status": "PASS",
        "search_status": "UNKNOWN_search_cap",
        "source_commit": COMMIT,
        "small_object_sha256": SMALL_SHA,
        "envelope_object_sha256": ENVELOPE_SHA,
        "binary_sha256": BINARY_SHA,
        "trials_attempted": 1,
        "solver_conflicts": 5000,
        "solver_models": 0,
        "graph_accepted_candidates": 0,
        "completed_rows": 0,
        "conservative_charge_seconds": 2400,
        "cumulative_conservative_charge_seconds": 15000,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }


if __name__ == "__main__":
    print(json.dumps(check(), indent=2, sort_keys=True))
