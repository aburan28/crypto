#!/usr/bin/env python3
"""Audit retained K600 mixed-orbit receipts without rerunning the solver."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / "evidence-01"
REPO = HERE.parents[3]
STUDY = REPO / "research/koblitz_n83_factor_base_sweep_20261008"
sys.path.insert(0, str(STUDY))

from chain_search_supervisor import parse_point  # noqa: E402
from k600_mixed_orbit_supervisor import (  # noqa: E402
    CORPUS, DOMAIN_CLAUSES_3_PLUS_2, ENVELOPE_BLAKE3,
    ENVELOPE_SHA256, FIRST_TRIAL_A, FIRST_TRIAL_B, IMAGE,
    PRIOR_IMPORT_RECEIPT_SHA256, PRIOR_SEARCH_AUDIT_SHA256,
    SMALL_BLAKE3, SMALL_SHA256, SOURCE_PATHS, replay_events,
    search_receipt_error, verify_receipt_error,
)


COMMIT = "f3d27f68b68cda907e0668042ac7fa13ae6b579d"
BINARY_SHA256 = "26dcb2b49388fcb89853eabbda01db71d148620cd465bd5458099e0606dcf606"
EXTERNAL = Path("/Volumes/SSD990/crypto/worktrees/n83-k600-first-run-evidence-20261012")


def require(condition: bool, label: str) -> None:
    if not condition:
        raise ValueError(label)


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record(path: str) -> dict:
    return json.loads((EVIDENCE / path).read_text())


def option(command: list[str], flag: str) -> str:
    return command[command.index(flag) + 1]


def stage_ok(stage: dict, name: str, wall: int, cpus: str, memory: str,
             folder: str) -> None:
    require(stage["name"] == name and stage["exit_code"] == 0
            and stage["timed_out"] is False and stage["launch_error"] is None
            and stage["cleanup_verified"] is True and stage["cleanup_error"] is None,
            f"{name} disposition or cleanup")
    state = stage["container_state"]
    require(state["ExitCode"] == 0 and state["OOMKilled"] is False
            and state["Running"] is False and state["Status"] == "exited"
            and 0 <= stage["wall_elapsed_seconds_L0"] <= wall,
            f"{name} container exit or wall")
    command = stage["command"]
    require(option(command, "--pull") == "never"
            and option(command, "--platform") == "linux/arm64"
            and option(command, "--network") == "none"
            and option(command, "--cpus") == cpus
            and option(command, "--memory") == memory
            and option(command, "--memory-swap") == memory
            and IMAGE in command, f"{name} isolation")
    if name != "build":
        require("--read-only" in command
                and option(command, "--pids-limit") == "64"
                and option(command, "--security-opt") == "no-new-privileges"
                and all(any(f"dst=/{item},readonly" in arg for arg in command)
                        for item in ("worker", "small.jsonl.gz", "envelope.jsonl.gz")),
                f"{name} mounts")
    for suffix in ("stdout", "stderr"):
        require(sha(EVIDENCE / folder / f"{name}.{suffix}")
                == stage[f"{suffix}_sha256"], f"{name} {suffix} hash")


def check(small: Path | None = None, envelope: Path | None = None) -> dict:
    auth, objects = record("authorization.json"), record("objects.json")
    small = small or Path(auth["build_dir"]).parent / "inputs" / (SMALL_BLAKE3 + ".jsonl.gz")
    envelope = envelope or Path(auth["build_dir"]).parent / "inputs" / (ENVELOPE_BLAKE3 + ".jsonl.gz")
    build, config = record("build/build.json"), record("run/config.json")
    outer, budget = record("run/outer.json"), record("run/budget.json")
    worker_config, result = record("run/run/config.json"), record("run/run/result.json")
    require(auth == {
        "schema": "n83.additional-local-budget/v4",
        "study": "koblitz_n83_factor_base_sweep_20261008",
        "purpose": "mixed_orbit_k600_two_residual_first_gate",
        "prior_total_charged_seconds": "15000",
        "prior_search_audit_sha256": PRIOR_SEARCH_AUDIT_SHA256,
        "prior_import_receipt_sha256": PRIOR_IMPORT_RECEIPT_SHA256,
        "additional_authorized_seconds": "2400",
        "build_dir": str(EXTERNAL / "build"),
        "output_dir": str(EXTERNAL / "run"),
        "user_approval_reference": (
            "2026-10-12 user yes to 2400 additional local seconds for frozen first "
            "retained-N83 mixed-orbit build, import, bounded search, audit, and cleanup; annotation 1"),
    }, "authorization record")
    require(objects["schema"] == "n83.k600-retained-object-inputs/v1"
            and objects["checked_utc_date"] == "2026-10-12"
            and len(objects["objects"]) == 2, "object inventory")
    for row, label, digest, blake3, length, suffix, etag, modified in zip(
        objects["objects"], ("small", "envelope"),
        (SMALL_SHA256, ENVELOPE_SHA256), (SMALL_BLAKE3, ENVELOPE_BLAKE3),
        (3025798, 10313647),
        ("a0/objects/", "v2-size-frontier/a0/objects/"),
        ('"359c027f8c0be18a8e41696f6c5c7aed"',
         '"96006317b77bf3a271eb8330ab209ce5-2"'),
        ("2026-10-09T17:03:04+00:00", "2026-10-10T18:17:29+00:00"),
    ):
        require(row["label"] == label and row["local_sha256"] == digest
                and row["compressed_blake3"] == blake3
                and row["content_length"] == length
                and row["s3_uri"] == (
                    "s3://crypto-autoresearcher/factor-bases/icv1/etc/"
                    "koblitz_n83_factor_base_sweep_20261008/" + suffix + blake3 + ".jsonl.gz")
                and row["etag"] == etag and row["last_modified"] == modified,
                f"{label} S3 identity")
    for path, digest in ((small, SMALL_SHA256), (envelope, ENVELOPE_SHA256)):
        require(sha(path) == digest, "local object hash")
    auth_sha = sha(EVIDENCE / "authorization.json")
    require(build["schema"] == "n83.mixed-orbit-k600-linux-build/v1"
            and build["status"] == "PASS_build_only"
            and build["source_commit"] == COMMIT
            and build["authorization_sha256"] == auth_sha
            and build["container_image_id"] == IMAGE
            and build["binary"] == str(Path(auth["build_dir"]) / "worker")
            and build["binary_sha256"] == BINARY_SHA256
            and build["conservative_charge_seconds"] == 600
            and 0 <= build["elapsed_seconds_L0"] <= 600
            and build["retained_n83_import_executed"] is False
            and build["retained_n83_search_executed"] is False,
            "frozen build receipt")
    require(sha(EVIDENCE / "build/worker") == BINARY_SHA256
            and (EVIDENCE / "build/worker").read_bytes()[:4] == b"\x7fELF",
            "archived Linux worker")
    stage_ok(build["stage"], "build", 450, "4", "8192m", "build")
    require(build["stage"]["command"][-2:] ==
            ["--example", "koblitz_n83_mixed_orbit_worker"], "build command")
    require(set(build["source_sha256"]) == set(SOURCE_PATHS), "frozen source set")
    for relative, digest in build["source_sha256"].items():
        blob = subprocess.run(["git", "show", f"{COMMIT}:{relative}"],
                              cwd=REPO, capture_output=True, check=True).stdout
        require(hashlib.sha256(blob).hexdigest() == digest, f"source blob {relative}")
    inputs = config["input_sha256"]
    fixed = {
        str(Path(auth["build_dir"]) / "worker"): BINARY_SHA256,
        str(Path(auth["build_dir"]) / "build.json"): sha(EVIDENCE / "build/build.json"),
        str(Path(auth["build_dir"]).parent / "authorization.json"): auth_sha,
    }
    for path, digest in fixed.items():
        require(inputs.get(path) == digest, f"input {path}")
    for label, digest in (("small", SMALL_SHA256), ("envelope", ENVELOPE_SHA256)):
        matches = [value for path, value in inputs.items()
                   if path.endswith("/" + objects["objects"][int(label == "envelope")]["compressed_blake3"] + ".jsonl.gz")]
        require(matches == [digest], f"{label} input hash")
    for relative, digest in build["source_sha256"].items():
        matches = [value for path, value in inputs.items() if path.endswith("/" + relative)]
        require(matches == [digest], f"source input {relative}")
    require(len(inputs) == len(SOURCE_PATHS) + 5, "input set cardinality")
    require(config["schema"] == "n83.mixed-orbit-k600-supervisor-config/v1"
            and config["source_commit"] == COMMIT
            and config["authorization_sha256"] == auth_sha
            and config["container_image_id"] == IMAGE
            and config["build_cap_seconds"] == 600
            and config["verify_cap_seconds"] == 360
            and config["search_cap_seconds"] == 1200
            and config["worker_wall_seconds"] == 1150
            and config["cleanup_reserve_seconds"] == 240
            and config["conservative_total_charge_seconds"] == 2400
            and config["max_trials"] == 1
            and config["conflicts_per_trial"] == 5000
            and config["max_models_per_trial"] == 32
            and config["small_summands"] == 3
            and config["residual_summands"] == 2
            and config["total_index_calculus_runtime_ms"] is None,
            "run configuration")
    require(outer["schema"] == "n83.mixed-orbit-k600-supervisor-outer/v1"
            and outer["status"] == "UNKNOWN_search_cap"
            and outer["source_commit"] == COMMIT
            and outer["config_sha256"] == sha(EVIDENCE / "run/config.json")
            and outer["authorization_sha256"] == auth_sha
            and outer["import_passed"] is True
            and outer["event_replay_error"] is None
            and outer["frozen_input_changed"] is False
            and 0 <= outer["elapsed_supervisor_seconds_L0"] <= 1800
            and outer["total_index_calculus_runtime_ms"] is None
            and outer["selected_best_total_runtime"] is None,
            "outer disposition")
    require([row["name"] for row in outer["stages"]] == ["verify", "search"],
            "stage sequence")
    for row in outer["stages"]:
        stage_ok(row, row["name"], 360 if row["name"] == "verify" else 1200,
                 "1", "4096m", "run")
        for source, target in ((EXTERNAL / "build/worker", "worker"),
                               (EXTERNAL / "inputs" / (SMALL_BLAKE3 + ".jsonl.gz"),
                                "small.jsonl.gz"),
                               (EXTERNAL / "inputs" / (ENVELOPE_BLAKE3 + ".jsonl.gz"),
                                "envelope.jsonl.gz")):
            require(f"type=bind,src={source},dst=/{target},readonly" in row["command"],
                    f"{row['name']} input mount {target}")
    require(outer["stages"][0]["command"][-6:] ==
            ["/worker", "verify-alt", "/small-manifest.json", "/small.jsonl.gz",
             "/envelope-manifest.json", "/envelope.jsonl.gz"], "verify command")
    require(outer["stages"][1]["command"][-14:] ==
            ["/worker", "run-alt", "/small-manifest.json", "/small.jsonl.gz",
             "/envelope-manifest.json", "/envelope.jsonl.gz", "/corpus.json",
             "1150", "1", "5000", "32", "3", "2", "/out/run"], "search command")
    require(verify_receipt_error(record("run/verify.stdout")) is None,
            "K600 import receipt")
    require(search_receipt_error(worker_config, result) is None,
            "worker search receipt")
    require(result["status"] == "UNKNOWN_search_cap"
            and result["trials_attempted"] == 1
            and result["solver_exhausted_trials"] == 1
            and result["solver_calls"] == 1
            and result["solver_conflicts"] == 5000
            and 0 < result["memory_cgroup_peak_bytes"] <= 4 * 1024**3
            and all(result[key] == 0 for key in (
                "solver_models", "invalid_models", "group_rejected_models",
                "graph_accepted_candidates", "graph_pivot_rows",
                "graph_full_from_partials", "graph_capped_rows", "completed_rows")),
            "bounded natural trial")
    events = (EVIDENCE / "run/run/events.jsonl").read_text().splitlines()
    require(len(events) == 1, "one event")
    event = json.loads(events[0])
    search = event["search"]
    require(event["trial"] == 0 and (event["a"], event["b"]) ==
            (FIRST_TRIAL_A, FIRST_TRIAL_B)
            and event["candidate"] is None and event["graph"] is None
            and search["domain_clauses"] == DOMAIN_CLAUSES_3_PLUS_2
            and search["conflicts"] == 5000 and search["models"] == 0
            and search["exhausted"] is True and search["refuted"] is False,
            "first-trial event")
    corpus = json.loads((REPO / CORPUS).read_text())
    selected = [row for row in corpus["targets"]
                if row.get("a") == 0 and row.get("fixture") == 0]
    require(len(selected) == 1, "public target selection")
    target = parse_point(selected[0]["point"])
    require(target is not None, "public target point")
    require(replay_events(EVIDENCE / "run/run/events.jsonl", small, envelope,
                          target, result) is None, "independent event replay")
    require(budget == {
        "schema": "n83.mixed-orbit-k600-supervisor-budget/v1",
        "authorization_sha256": auth_sha,
        "prior_total_charged_seconds": "15000",
        "prior_search_audit_sha256": PRIOR_SEARCH_AUDIT_SHA256,
        "prior_import_receipt_sha256": PRIOR_IMPORT_RECEIPT_SHA256,
        "additional_authorized_seconds": "2400",
        "build_conservative_charge_seconds": 600,
        "run_conservative_charge_seconds": 1800,
        "conservative_charge_seconds": 2400,
        "authorized_full_gate_seconds": 2400,
        "remaining_additional_seconds": "0",
        "round_status": "UNKNOWN_search_cap",
        "selected_best_total_runtime": None,
    }, "full conservative charge")
    return {
        "schema": "n83.k600-mixed-orbit-first-gate-audit/v1",
        "status": "PASS",
        "search_status": "UNKNOWN_search_cap",
        "source_commit": COMMIT,
        "binary_sha256": BINARY_SHA256,
        "small_object_sha256": SMALL_SHA256,
        "envelope_object_sha256": ENVELOPE_SHA256,
        "trials_attempted": 1,
        "solver_conflicts": 5000,
        "solver_models": 0,
        "graph_accepted_candidates": 0,
        "completed_rows": 0,
        "domain_clauses": DOMAIN_CLAUSES_3_PLUS_2,
        "conservative_charge_seconds": 2400,
        "cumulative_conservative_charge_seconds": 17400,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--small-object", type=Path)
    parser.add_argument("--envelope-object", type=Path)
    arguments = parser.parse_args()
    print(json.dumps(check(arguments.small_object, arguments.envelope_object),
                     indent=2, sort_keys=True))
