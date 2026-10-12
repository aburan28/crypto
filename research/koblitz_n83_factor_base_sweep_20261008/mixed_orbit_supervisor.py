#!/usr/bin/env python3
"""Build and supervise the pinned N83 one/two-residual producer gate.

Preparing or testing this file launches no retained N83 job. ``build`` and
``run`` require a separate, run-specific grant beyond the fully charged
12,600-second local allowance. Each stage is launched once and charged its
entire declared cap even when it fails early.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import uuid

from chain_search_supervisor import (
    PRIMARY_GENERATOR, PRIMARY_ORDER, canonical_decimal, parse_point,
    point_add, point_mul,
)
from compact_probe_supervisor import clean_commit, sha256_file, write_new_json
from structured_orbit_search_supervisor import stage


STUDY = "koblitz_n83_factor_base_sweep_20261008"
ROOT = Path(__file__).resolve().parents[2]
IMAGE = "sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31"
SMALL_BLAKE3 = "76e1cbe4ea5a9d26e5b24a2d27f9401271dc471846529190a0416d8e09837a92"
ENVELOPE_BLAKE3 = "82990cda5d57b1c4c4b635673a881ad233c24f707e7ac644a11478aeddddf851"
SMALL_SHA256 = "eb342b574e244e1e67e8e78c805f33ee88edcdb4f3a14d4114e6f413beede380"
ENVELOPE_SHA256 = "aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a"
CORPUS_BLAKE3 = "831eaa7be2387634fd610f6706103b2af9b4e274f0db59923c0280af3cbed74b"
SMALL_MANIFEST = f"research/{STUDY}/verification/matched-k1182-k2048-900s-20261011/evidence-01/panel-k1182/manifest.json"
ENVELOPE_MANIFEST = f"research/{STUDY}/verification/remaining-budget-k2048-20261010/manifest.json"
CORPUS = f"research/{STUDY}/pilot-01/probe-corpus.json"
PRIOR_IMPORT_AUDIT = f"research/{STUDY}/verification/mixed-orbit-import-900s-20261012/audit.json"
PRIOR_IMPORT_AUDIT_SHA256 = "76a8bf74e46feb0ef4fbe650395f39733f4ff95d31c35f8829cb506a9ad04c7d"
SOURCE_PATHS = (
    "examples/koblitz_n83_mixed_orbit_worker.rs",
    "src/cryptanalysis/koblitz_large_prime.rs",
    "src/cryptanalysis/koblitz_index_calculus.rs",
    "src/cryptanalysis/binary_semaev_chain_sat.rs",
    f"research/{STUDY}/mixed_orbit_supervisor.py",
    f"research/{STUDY}/chain_search_supervisor.py",
    f"research/{STUDY}/structured_orbit_search_supervisor.py",
    f"research/{STUDY}/compact_probe_supervisor.py",
    SMALL_MANIFEST, ENVELOPE_MANIFEST, CORPUS, PRIOR_IMPORT_AUDIT,
)
BUILD_CAP = 600
BUILD_WALL = 450
VERIFY_CAP = 360
SEARCH_CAP = 1200
WORKER_WALL = 1150
CLEANUP_RESERVE = 240
RUN_CAP = VERIFY_CAP + SEARCH_CAP + CLEANUP_RESERVE
TOTAL_CHARGE = BUILD_CAP + RUN_CAP
TRIALS = 1
CONFLICTS = 5000
MODELS = 32
SMALL_SUMMANDS = 3
RESIDUAL_SUMMANDS = 2
# Independently derived from the worker's published BLAKE3 trial schedule:
# BLAKE3("n83.mixed-orbit-partial-trial/v1" || label || 0u64_le) mod r.
FIRST_TRIAL_A = "1833078656608906180176350"
FIRST_TRIAL_B = "163089659398503559150141"


def authorization_error(record: object, build: Path, output: Path) -> str | None:
    if not isinstance(record, dict):
        return "additional-budget record must be an object"
    expected = {
        "schema": "n83.additional-local-budget/v2",
        "study": STUDY,
        "purpose": "mixed_orbit_two_residual_first_gate",
        "prior_total_charged_seconds": "12600",
        "prior_import_audit_sha256": PRIOR_IMPORT_AUDIT_SHA256,
        "build_dir": str(build.resolve()),
        "output_dir": str(output.resolve()),
    }
    for key, value in expected.items():
        if record.get(key) != value:
            return f"additional-budget {key} mismatch"
    amount = record.get("additional_authorized_seconds")
    if (not isinstance(amount, str) or not amount.isascii() or not amount.isdecimal()
            or amount != str(int(amount)) or Decimal(amount) < TOTAL_CHARGE):
        return "additional-budget amount below full conservative charge"
    approval = record.get("user_approval_reference")
    if not isinstance(approval, str) or not approval.strip():
        return "additional-budget approval reference missing"
    return None


def load_authorization(path: Path, build: Path, output: Path) -> tuple[dict, str]:
    raw = path.read_bytes()
    record = json.loads(raw)
    error = authorization_error(record, build, output)
    if error:
        raise ValueError(error)
    return record, hashlib.sha256(raw).hexdigest()


def outside_fresh_checkout(path: Path, checkout: Path) -> Path:
    path = path.resolve()
    if path.is_relative_to(checkout) or path.exists():
        raise ValueError("output must be a fresh path outside the checkout")
    return path


def image_identity(image: str = IMAGE) -> str:
    inspected = subprocess.run(
        ["docker", "image", "inspect", image, "--format", "{{.Id}} {{.Architecture}} {{.Os}}"],
        capture_output=True, text=True, timeout=15, check=True,
    )
    if inspected.stdout.strip() != f"{IMAGE} arm64 linux" or image != IMAGE:
        raise ValueError("pinned Linux arm64 image unavailable")
    return IMAGE


def source_hashes(checkout: Path) -> dict[str, str]:
    return {relative: sha256_file(checkout / relative) for relative in SOURCE_PATHS}


def prior_import_evidence_error(checkout: Path) -> str | None:
    path = checkout / PRIOR_IMPORT_AUDIT
    try:
        if sha256_file(path) != PRIOR_IMPORT_AUDIT_SHA256:
            return "prior import audit bytes changed"
        evidence = json.loads(path.read_text())
    except (OSError, UnicodeError, json.JSONDecodeError):
        return "prior import audit missing or malformed"
    if (not isinstance(evidence, dict)
            or evidence.get("schema") != "n83.mixed-orbit-import-900-audit/v1"
            or evidence.get("status") != "PASS"
            or evidence.get("conservative_charge_seconds") != 900
            or evidence.get("solver_search_executed") is not False
            or evidence.get("total_index_calculus_runtime_ms") is not None):
        return "prior import audit disposition mismatch"
    return None


def is_linux_elf(path: Path) -> bool:
    if not path.is_file():
        return False
    with path.open("rb") as stream:
        return stream.read(4) == b"\x7fELF"


def build(authorization: Path, build_dir: Path, output_dir: Path,
          checkout: Path = ROOT) -> dict:
    checkout = checkout.resolve(strict=True)
    build_dir = outside_fresh_checkout(build_dir, checkout)
    output_dir = output_dir.resolve()
    if output_dir.is_relative_to(checkout) or output_dir.exists():
        raise ValueError("run output must be a fresh path outside the checkout")
    authorization = authorization.resolve(strict=True)
    _, auth_hash = load_authorization(authorization, build_dir, output_dir)
    commit = clean_commit(checkout)
    prior_error = prior_import_evidence_error(checkout)
    if prior_error:
        raise ValueError(prior_error)
    sources = source_hashes(checkout)
    image = image_identity()
    registry = Path(os.environ.get("CARGO_HOME", str(Path.home() / ".cargo"))) / "registry"
    if not registry.is_dir():
        raise ValueError("offline Cargo registry unavailable")
    started = time.monotonic()
    build_dir.mkdir(parents=True, exist_ok=False)
    target = build_dir / "target"
    target.mkdir()
    container_name = f"n83-mixed-build-{uuid.uuid4().hex[:16]}"
    command = [
        "docker", "run", "--pull", "never", "--platform", "linux/arm64",
        "--name", container_name, "--network", "none", "--cpus", "4",
        "--memory", "8192m", "--memory-swap", "8192m", "--pids-limit", "256",
        "--cap-drop", "ALL", "--security-opt", "no-new-privileges",
        "--user", f"{os.getuid()}:{os.getgid()}",
        "--tmpfs", "/cargo:rw,mode=1777",
        "--mount", f"type=bind,src={registry},dst=/cargo/registry,readonly",
        "--mount", f"type=bind,src={checkout},dst=/src,readonly",
        "--mount", f"type=bind,src={target},dst=/target",
        "--env", "CARGO_HOME=/cargo", "--env", "CARGO_TARGET_DIR=/target",
        "--workdir", "/src", image, "cargo", "build", "--release", "--offline",
        "--locked", "--no-default-features", "--example", "koblitz_n83_mixed_orbit_worker",
    ]
    launched = stage(name="build", command=command, wall=BUILD_WALL,
                     output=build_dir, container_name=container_name)
    artifact = target / "release" / "examples" / "koblitz_n83_mixed_orbit_worker"
    binary = build_dir / "worker"
    try:
        frozen = clean_commit(checkout) == commit and source_hashes(checkout) == sources
    except (OSError, ValueError, subprocess.SubprocessError):
        frozen = False
    state = launched.get("container_state")
    status = "PRODUCER_FAILURE_build"
    if not launched.get("cleanup_verified"):
        status = "PRODUCER_FAILURE_cleanup"
    elif launched.get("timed_out"):
        status = "UNKNOWN_build_wall_cap"
    elif not frozen:
        status = "PRODUCER_FAILURE_source_changed"
    elif (launched.get("exit_code") == 0 and not launched.get("launch_error")
          and isinstance(state, dict) and state.get("ExitCode") == 0
          and state.get("OOMKilled") is False and is_linux_elf(artifact)):
        try:
            shutil.copyfile(artifact, binary)
            binary.chmod(0o755)
            status = "PASS_build_only"
        except OSError:
            status = "PRODUCER_FAILURE_binary"
    elapsed = time.monotonic() - started
    if elapsed > BUILD_CAP:
        status = "PRODUCER_FAILURE_budget_overshoot"
    receipt = {
        "schema": "n83.mixed-orbit-linux-build/v1", "status": status,
        "source_commit": commit, "source_sha256": sources,
        "authorization_sha256": auth_hash, "container_image_id": image,
        "binary_sha256": sha256_file(binary) if binary.is_file() else None,
        "binary": str(binary) if binary.is_file() else None,
        "stage": launched, "elapsed_seconds_L0": elapsed,
        "conservative_charge_seconds": BUILD_CAP,
        "retained_n83_import_executed": False,
        "retained_n83_search_executed": False,
    }
    write_new_json(build_dir / "build.json", receipt)
    return receipt


def verify_receipt_error(record: object) -> str | None:
    expected = {
        "schema": "n83.mixed-orbit-import/v1", "status": "PASS",
        "curve_slug": "icv1-f2m83-tm6151469093347-debefd74",
        "small_object_blake3": SMALL_BLAKE3,
        "envelope_object_blake3": ENVELOPE_BLAKE3,
        "small_points": 196212, "small_signed_orbits": 1182,
        "residual_points": 143756, "residual_signed_orbits": 866,
        "solver_search_executed": False,
    }
    if not isinstance(record, dict):
        return "import receipt missing"
    for key, value in expected.items():
        if type(record.get(key)) is not type(value) or record[key] != value:
            return f"import {key} mismatch"
    return None


def finite_nonnegative(value: object) -> bool:
    return (not isinstance(value, bool) and isinstance(value, (int, float))
            and math.isfinite(value) and value >= 0)


def valid_worker_cap(record: object) -> bool:
    return (isinstance(record, dict)
            and record.get("schema") == "n83.mixed-orbit-partial-cap/v1"
            and record.get("status") == "UNKNOWN_budget"
            and type(record.get("budget_seconds")) is int
            and record["budget_seconds"] == WORKER_WALL
            and finite_nonnegative(record.get("elapsed_process_wall_ms"))
            and record["elapsed_process_wall_ms"] <= SEARCH_CAP * 1000)


def search_receipt_error(config: object, result: object) -> str | None:
    if not isinstance(config, dict) or not isinstance(result, dict):
        return "search config or result missing"
    expected_config = {
        "schema": "n83.mixed-orbit-partial-config/v1", "study": STUDY,
        "small_object_blake3": SMALL_BLAKE3,
        "envelope_object_blake3": ENVELOPE_BLAKE3,
        "public_corpus_path": "/corpus.json",
        "public_corpus_canonical_json_blake3": CORPUS_BLAKE3,
        "budget_seconds": WORKER_WALL, "max_trials": TRIALS,
        "conflict_budget_per_trial": CONFLICTS,
        "max_models_per_trial": MODELS,
        "small_summands": SMALL_SUMMANDS,
        "residual_summands": RESIDUAL_SUMMANDS,
        "order_first_pair": True, "max_variables": 150000,
        "max_domain_clauses": 12000000,
        "memory_cgroup_limit_bytes": 4 * 1024**3,
        "memory_cgroup_swap_limit_bytes": 0,
    }
    for key, value in expected_config.items():
        if type(config.get(key)) is not type(value) or config[key] != value:
            return f"config {key} mismatch"
    if config.get("graph_limits") != {
        "max_merge_depth": 32, "max_row_terms": 4096, "max_pivots": 1000000,
    }:
        return "graph limits mismatch"
    if (not isinstance(config.get("source_blake3"), str)
            or re.fullmatch(r"[0-9a-f]{64}", config["source_blake3"]) is None):
        return "worker source hash missing"
    expected_result = {
        "schema": "n83.mixed-orbit-partial-result/v1", "study": STUDY,
        "small_object_blake3": SMALL_BLAKE3,
        "envelope_object_blake3": ENVELOPE_BLAKE3,
        "public_corpus_canonical_json_blake3": CORPUS_BLAKE3,
        "small_signed_orbits": 1182, "residual_signed_orbits": 866,
        "matrix_rank": None, "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    for key, value in expected_result.items():
        if type(result.get(key)) is not type(value) or result[key] != value:
            return f"result {key} mismatch"
    counters = (
        "trials_attempted", "solver_exhausted_trials", "solver_refuted_trials",
        "completed_rows", "solver_calls", "solver_conflicts", "solver_models",
        "invalid_models", "group_rejected_models", "graph_accepted_candidates",
        "graph_pivot_rows", "graph_full_from_partials", "graph_capped_rows",
        "graph_pivot_reductions",
    )
    if any(type(result.get(key)) is not int or result[key] < 0 for key in counters):
        return "malformed search or graph counters"
    if (not 1 <= result["trials_attempted"] <= TRIALS
            or result["solver_exhausted_trials"] + result["solver_refuted_trials"] > result["trials_attempted"]
            or result["completed_rows"] + result["graph_pivot_rows"] + result["graph_capped_rows"]
            != result["graph_accepted_candidates"]
            or result["graph_full_from_partials"] > result["completed_rows"]
            or result["graph_accepted_candidates"] > result["solver_models"]
            or result["solver_conflicts"] > TRIALS * CONFLICTS
            or result["solver_models"] > TRIALS * MODELS
            or result["invalid_models"] + result["group_rejected_models"] > result["solver_models"]):
        return "search or graph counters disagree"
    if (not finite_nonnegative(result.get("base_import_elapsed_ms_L0"))
            or not finite_nonnegative(result.get("worker_wall_ms_L0"))
            or not isinstance(result.get("memory_cgroup_peak_bytes"), int)
            or isinstance(result["memory_cgroup_peak_bytes"], bool)
            or not 0 < result["memory_cgroup_peak_bytes"] <= 4 * 1024**3):
        return "search timing or cgroup peak malformed"
    status = result.get("status")
    if status not in {
        "PASS_bounded_partial_search", "UNKNOWN_preflight_cap", "UNKNOWN_search_cap",
        "UNKNOWN_graph_cap", "PRODUCER_FAILURE_invalid_model",
    }:
        return "search status inadmissible"
    if status == "PASS_bounded_partial_search" and (
        result["solver_exhausted_trials"] or result["graph_capped_rows"]
        or result["invalid_models"] or result["trials_attempted"] != TRIALS
    ):
        return "bounded-search status hides a cap or invalid model"
    if status == "UNKNOWN_preflight_cap" and (
        result["trials_attempted"] != 1 or result["solver_calls"] != 0
        or result["solver_exhausted_trials"] != 1
    ):
        return "preflight status has solver work"
    if status == "UNKNOWN_search_cap" and result["solver_exhausted_trials"] == 0:
        return "search-cap status has no exhausted trial"
    if status == "UNKNOWN_graph_cap" and result["graph_capped_rows"] == 0:
        return "graph-cap status has no capped row"
    if status == "PRODUCER_FAILURE_invalid_model" and result["invalid_models"] == 0:
        return "invalid-model status has no invalid model"
    return None


def _group_sum(points: list[tuple[int, int] | None]) -> tuple[int, int] | None:
    total = None
    for point in points:
        total = point_add(total, point)
    return total


def replay_equation(a: str, b: str, small_points: list[tuple[int, int]],
                    residual_points: list[tuple[int, int]], target: tuple[int, int]) -> bool:
    if not canonical_decimal(a, PRIMARY_ORDER) or not canonical_decimal(b, PRIMARY_ORDER):
        return False
    lhs = point_add(point_mul(PRIMARY_GENERATOR, int(a)), point_mul(target, int(b)))
    rhs = _group_sum([*small_points, *residual_points])
    return lhs == rhs


def replay_events(path: Path, small_object: Path, target: tuple[int, int],
                  result: dict) -> str | None:
    try:
        events = [json.loads(line) for line in path.read_text().splitlines()]
    except (OSError, UnicodeError, json.JSONDecodeError):
        return "event ledger missing or malformed"
    if len(events) != result["trials_attempted"]:
        return "event ledger length mismatch"
    needed_indices: set[int] = set()
    summed = {key: 0 for key in (
        "solver_calls", "solver_conflicts", "solver_models", "invalid_models",
        "group_rejected_models", "solver_exhausted_trials", "solver_refuted_trials",
        "graph_accepted_candidates", "graph_pivot_rows", "graph_capped_rows",
        "completed_rows", "graph_full_from_partials", "graph_pivot_reductions",
    )}
    for trial, event in enumerate(events):
        if not isinstance(event, dict) or event.get("trial") != trial:
            return "event trial order mismatch"
        if (not canonical_decimal(event.get("a"), PRIMARY_ORDER)
                or not canonical_decimal(event.get("b"), PRIMARY_ORDER)
                or not finite_nonnegative(event.get("elapsed_process_wall_ms"))):
            return "event source equation or wall time malformed"
        if trial == 0 and (event["a"], event["b"]) != (FIRST_TRIAL_A, FIRST_TRIAL_B):
            return "event differs from frozen first-trial coefficient schedule"
        search = event.get("search")
        if not isinstance(search, dict):
            return "event search statistics missing"
        for source, destination in (
            ("solver_calls", "solver_calls"), ("conflicts", "solver_conflicts"),
            ("models", "solver_models"), ("invalid_models", "invalid_models"),
            ("group_rejected_models", "group_rejected_models"),
        ):
            value = search.get(source)
            if type(value) is not int or value < 0:
                return "event search counter malformed"
            if source == "conflicts" and value > CONFLICTS:
                return "event conflict cap exceeded"
            if source == "models" and value > MODELS:
                return "event model cap exceeded"
            summed[destination] += value
        for source, destination in (("exhausted", "solver_exhausted_trials"),
                                    ("refuted", "solver_refuted_trials")):
            if type(search.get(source)) is not bool:
                return "event search disposition malformed"
            summed[destination] += int(search[source])
        candidate = event.get("candidate")
        if candidate is not None:
            if (not isinstance(candidate, dict) or candidate.get("source_id") != trial + 1
                    or candidate.get("a") != event.get("a")
                    or candidate.get("b") != event.get("b")
                    or not isinstance(candidate.get("small_point_indices"), list)
                    or len(candidate["small_point_indices"]) != SMALL_SUMMANDS
                    or any(type(i) is not int or not 0 <= i < 196212
                           for i in candidate["small_point_indices"])):
                return "partial source or small indices mismatch"
            needed_indices.update(candidate["small_point_indices"])
            summed["graph_accepted_candidates"] += 1
        graph = event.get("graph")
        if candidate is None and graph is not None:
            return "graph outcome without partial"
        if candidate is not None and (not isinstance(graph, dict)
                                      or graph.get("kind") not in {"pivot", "full", "UNKNOWN_graph_cap"}):
            return "partial has no graph outcome"
        if candidate is not None:
            if type(graph.get("merge_depth")) is not int or graph["merge_depth"] < 0:
                return "graph merge depth malformed"
            summed["graph_pivot_reductions"] += graph["merge_depth"]
            if graph["kind"] == "pivot":
                summed["graph_pivot_rows"] += 1
            elif graph["kind"] == "full":
                summed["completed_rows"] += 1
                summed["graph_full_from_partials"] += int(graph["merge_depth"] > 0)
            else:
                summed["graph_capped_rows"] += 1
    if any(result[key] != value for key, value in summed.items()):
        return "event ledger counters disagree with worker result"
    if not any(event.get("candidate") is not None for event in events):
        return None
    points: dict[int, tuple[int, int]] = {}
    try:
        with gzip.open(small_object, "rt", encoding="utf-8") as archive:
            header = json.loads(next(archive))
            representatives = [parse_point(value) for value in header["representatives"]]
            if len(representatives) != 1182 or any(p is None for p in representatives):
                return "small orbit representatives malformed"
            for index, line in enumerate(archive):
                if index in needed_indices:
                    point = parse_point(json.loads(line)["point"])
                    if point is None:
                        return "small point is infinity"
                    points[index] = point
    except (OSError, UnicodeError, ValueError, KeyError, StopIteration, json.JSONDecodeError):
        return "small archive replay failed"
    if set(points) != needed_indices:
        return "small point rows missing"
    sources: dict[int, tuple[int, int]] = {}
    for event in events:
        candidate = event.get("candidate")
        if candidate is None:
            continue
        try:
            residual = [parse_point(value) for value in candidate["residual_points"]]
        except (KeyError, TypeError, ValueError):
            return "residual point malformed"
        if (len(residual) != RESIDUAL_SUMMANDS or any(p is None for p in residual)
                or not replay_equation(candidate["a"], candidate["b"],
                                       [points[i] for i in candidate["small_point_indices"]],
                                       residual, target)):
            return "partial point equation failed independent replay"
        sources[candidate["source_id"]] = int(candidate["a"]), int(candidate["b"])
        graph = event["graph"]
        if graph["kind"] == "full":
            try:
                source_terms = graph["source_coeffs"]
                if (not isinstance(source_terms, list)
                        or any(not isinstance(item, list) or len(item) != 2
                               or not canonical_decimal(item[0], result["trials_attempted"] + 1)
                               or not canonical_decimal(item[1], PRIMARY_ORDER)
                               or int(item[0]) not in sources for item in source_terms)):
                    return "completed row source provenance malformed"
                if len({item[0] for item in source_terms}) != len(source_terms):
                    return "completed row repeats a source coefficient"
                combined_a = sum(int(c) * sources[int(i)][0] for i, c in source_terms) % PRIMARY_ORDER
                combined_b = sum(int(c) * sources[int(i)][1] for i, c in source_terms) % PRIMARY_ORDER
                if (not canonical_decimal(graph["a"], PRIMARY_ORDER)
                        or not canonical_decimal(graph["b"], PRIMARY_ORDER)
                        or combined_a != int(graph["a"]) or combined_b != int(graph["b"])):
                    return "completed row source combination disagrees"
                if any(not canonical_decimal(k, 1182)
                       or not canonical_decimal(v, PRIMARY_ORDER)
                       for k, v in graph["small"]):
                    return "completed row terms are not canonical"
                small_terms = [(int(k), int(v)) for k, v in graph["small"]]
                if any(not 0 <= i < 1182 or not 0 <= c < PRIMARY_ORDER
                       for i, c in small_terms):
                    return "completed row term outside signed orbit base"
                rhs = _group_sum([point_mul(representatives[i], c)
                                  for i, c in small_terms])
                lhs = point_add(point_mul(PRIMARY_GENERATOR, int(graph["a"])),
                                point_mul(target, int(graph["b"])))
            except (KeyError, TypeError, ValueError):
                return "completed row malformed"
            if lhs != rhs:
                return "completed row failed independent group replay"
    return None


def _inputs(checkout: Path, authorization: Path, binary: Path,
            small_object: Path, envelope_object: Path, build_receipt: Path) -> dict[str, str]:
    paths = [authorization, binary, small_object, envelope_object, build_receipt]
    paths += [checkout / relative for relative in SOURCE_PATHS]
    return {str(path): sha256_file(path) for path in paths}


def run(authorization: Path, build_dir: Path, output_dir: Path,
        small_object: Path, envelope_object: Path, checkout: Path = ROOT) -> dict:
    checkout = checkout.resolve(strict=True)
    build_dir = build_dir.resolve(strict=True)
    output_dir = outside_fresh_checkout(output_dir, checkout)
    authorization = authorization.resolve(strict=True)
    small_object = small_object.resolve(strict=True)
    envelope_object = envelope_object.resolve(strict=True)
    auth, auth_hash = load_authorization(authorization, build_dir, output_dir)
    commit = clean_commit(checkout)
    prior_error = prior_import_evidence_error(checkout)
    if prior_error:
        raise ValueError(prior_error)
    binary = build_dir / "worker"
    build_receipt = build_dir / "build.json"
    build_record = json.loads(build_receipt.read_text())
    build_stage = build_record.get("stage")
    build_state = build_stage.get("container_state") if isinstance(build_stage, dict) else None
    if (build_record.get("schema") != "n83.mixed-orbit-linux-build/v1"
            or build_record.get("status") != "PASS_build_only"
            or build_record.get("source_commit") != commit
            or build_record.get("source_sha256") != source_hashes(checkout)
            or build_record.get("authorization_sha256") != auth_hash
            or build_record.get("container_image_id") != IMAGE
            or build_record.get("binary") != str(binary)
            or build_record.get("binary_sha256") != sha256_file(binary)
            or build_record.get("conservative_charge_seconds") != BUILD_CAP
            or not isinstance(build_stage, dict)
            or build_stage.get("exit_code") != 0
            or build_stage.get("timed_out") is not False
            or build_stage.get("cleanup_verified") is not True
            or not isinstance(build_state, dict)
            or build_state.get("ExitCode") != 0
            or build_state.get("OOMKilled") is not False
            or not finite_nonnegative(build_record.get("elapsed_seconds_L0"))
            or build_record["elapsed_seconds_L0"] > BUILD_CAP
            or not os.access(binary, os.X_OK) or not is_linux_elf(binary)):
        raise ValueError("frozen Linux build receipt mismatch")
    if sha256_file(small_object) != SMALL_SHA256 or sha256_file(envelope_object) != ENVELOPE_SHA256:
        raise ValueError("downloaded S3 object SHA-256 mismatch")
    image = image_identity()
    inputs = _inputs(checkout, authorization, binary, small_object, envelope_object, build_receipt)
    started = time.monotonic()
    deadline = started + RUN_CAP
    output_dir.mkdir(parents=True, exist_ok=False)
    config = {
        "schema": "n83.mixed-orbit-supervisor-config/v1", "study": STUDY,
        "source_commit": commit, "input_sha256": inputs,
        "authorization_sha256": auth_hash, "container_image_id": image,
        "build_cap_seconds": BUILD_CAP, "verify_cap_seconds": VERIFY_CAP,
        "search_cap_seconds": SEARCH_CAP, "worker_wall_seconds": WORKER_WALL,
        "cleanup_reserve_seconds": CLEANUP_RESERVE,
        "conservative_total_charge_seconds": TOTAL_CHARGE,
        "max_trials": TRIALS, "conflicts_per_trial": CONFLICTS,
        "max_models_per_trial": MODELS,
        "small_summands": SMALL_SUMMANDS, "residual_summands": RESIDUAL_SUMMANDS,
        "total_index_calculus_runtime_ms": None,
    }
    write_new_json(output_dir / "config.json", config)
    config_hash = sha256_file(output_dir / "config.json")
    stages = []
    status = "PRODUCER_FAILURE_preflight"
    import_passed = False
    replay_error = None
    token = uuid.uuid4().hex[:16]
    for name, wall in (("verify", VERIFY_CAP), ("search", SEARCH_CAP)):
        if time.monotonic() + wall + CLEANUP_RESERVE > deadline:
            status = "UNKNOWN_budget_before_stage"
            break
        container_name = f"n83-mixed-{name}-{token}"
        command = [
            "docker", "run", "--pull", "never", "--platform", "linux/arm64",
            "--name", container_name, "--network", "none", "--cpus", "1",
            "--memory", "4096m", "--memory-swap", "4096m", "--pids-limit", "64",
            "--cap-drop", "ALL", "--security-opt", "no-new-privileges",
            "--read-only", "--tmpfs", "/tmp:rw,noexec,nosuid,size=64m",
            "--user", f"{os.getuid()}:{os.getgid()}",
            "--mount", f"type=bind,src={binary},dst=/worker,readonly",
            "--mount", f"type=bind,src={small_object},dst=/small.jsonl.gz,readonly",
            "--mount", f"type=bind,src={envelope_object},dst=/envelope.jsonl.gz,readonly",
            "--mount", f"type=bind,src={checkout / SMALL_MANIFEST},dst=/small-manifest.json,readonly",
            "--mount", f"type=bind,src={checkout / ENVELOPE_MANIFEST},dst=/envelope-manifest.json,readonly",
            "--mount", f"type=bind,src={checkout / CORPUS},dst=/corpus.json,readonly",
            "--mount", f"type=bind,src={output_dir},dst=/out",
            image, "/worker", name == "verify" and "verify" or "run",
            "/small-manifest.json", "/small.jsonl.gz",
            "/envelope-manifest.json", "/envelope.jsonl.gz",
        ]
        if name == "search":
            command += ["/corpus.json", str(WORKER_WALL), str(TRIALS), str(CONFLICTS),
                        str(MODELS), str(SMALL_SUMMANDS), str(RESIDUAL_SUMMANDS), "/out/run"]
        receipt = stage(name=name, command=command, wall=wall,
                        output=output_dir, container_name=container_name)
        stages.append(receipt)
        state = receipt.get("container_state")
        if not receipt.get("cleanup_verified"):
            status = "PRODUCER_FAILURE_cleanup"
            break
        if receipt.get("timed_out"):
            status = "UNKNOWN_wall_cap"
            break
        if state and state.get("OOMKilled") is True:
            status = "UNKNOWN_resource_cap"
            break
        if receipt.get("exit_code") == 124:
            try:
                cap = json.loads((output_dir / "run" / "cap.json").read_text())
            except (OSError, UnicodeError, json.JSONDecodeError):
                cap = None
            status = ("UNKNOWN_worker_wall_cap" if name == "search" and valid_worker_cap(cap)
                      else "PRODUCER_FAILURE_worker_exit")
            break
        if receipt.get("exit_code") in (134, 137):
            status = "UNKNOWN_resource_cap"
            break
        if (receipt.get("launch_error") or receipt.get("exit_code") != 0
                or not isinstance(state, dict) or state.get("ExitCode") != 0):
            status = "PRODUCER_FAILURE_worker_exit"
            break
        if name == "verify":
            try:
                record = json.loads((output_dir / "verify.stdout").read_text())
            except (OSError, UnicodeError, json.JSONDecodeError):
                record = None
            error = verify_receipt_error(record)
            if error:
                status = "PRODUCER_FAILURE_import_receipt"
                break
            status = "PASS_import_only"
            import_passed = True
        else:
            try:
                worker_config = json.loads((output_dir / "run" / "config.json").read_text())
                worker_result = json.loads((output_dir / "run" / "result.json").read_text())
            except (OSError, UnicodeError, json.JSONDecodeError):
                worker_config = worker_result = None
            error = search_receipt_error(worker_config, worker_result)
            if error:
                status = "PRODUCER_FAILURE_search_receipt"
                break
            try:
                corpus = json.loads((checkout / CORPUS).read_text())
                selected = [row for row in corpus["targets"]
                            if row.get("a") == 0 and row.get("fixture") == 0]
                target = parse_point(selected[0]["point"]) if len(selected) == 1 else None
            except (OSError, UnicodeError, ValueError, TypeError, KeyError, json.JSONDecodeError):
                target = None
            if target is None:
                status = "PRODUCER_FAILURE_public_target"
                break
            try:
                replay_error = replay_events(output_dir / "run" / "events.jsonl",
                                             small_object, target, worker_result)
            except (OSError, ValueError, TypeError, KeyError, ArithmeticError) as error:
                replay_error = f"independent replay raised {type(error).__name__}: {error}"
            status = "PRODUCER_FAILURE_event_replay" if replay_error else worker_result["status"]
    try:
        frozen_changed = (clean_commit(checkout) != commit
                          or any(sha256_file(Path(path)) != digest for path, digest in inputs.items())
                          or sha256_file(output_dir / "config.json") != config_hash)
    except (OSError, ValueError):
        frozen_changed = True
    if frozen_changed:
        status = "PRODUCER_FAILURE_frozen_input_changed"
    elapsed = time.monotonic() - started
    if elapsed > RUN_CAP:
        status = "PRODUCER_FAILURE_budget_overshoot"
    outer = {
        "schema": "n83.mixed-orbit-supervisor-outer/v1", "status": status,
        "source_commit": commit, "config_sha256": config_hash,
        "authorization_sha256": auth_hash, "stages": stages,
        "import_passed": import_passed, "event_replay_error": replay_error,
        "frozen_input_changed": frozen_changed,
        "elapsed_supervisor_seconds_L0": elapsed,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    write_new_json(output_dir / "outer.json", outer)
    run_charge = (sum(VERIFY_CAP if receipt["name"] == "verify" else SEARCH_CAP
                      for receipt in stages)
                  + (CLEANUP_RESERVE if stages else 0))
    total_charge = BUILD_CAP + run_charge
    write_new_json(output_dir / "budget.json", {
        "schema": "n83.mixed-orbit-supervisor-budget/v2",
        "authorization_sha256": auth_hash,
        "prior_total_charged_seconds": "12600",
        "prior_import_audit_sha256": PRIOR_IMPORT_AUDIT_SHA256,
        "additional_authorized_seconds": auth["additional_authorized_seconds"],
        "build_conservative_charge_seconds": BUILD_CAP,
        "run_conservative_charge_seconds": run_charge,
        "conservative_charge_seconds": total_charge,
        "authorized_full_gate_seconds": TOTAL_CHARGE,
        "remaining_additional_seconds": str(Decimal(auth["additional_authorized_seconds"]) - total_charge),
        "round_status": status, "selected_best_total_runtime": None,
    })
    return outer


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    for mode in ("build", "run"):
        sub = modes.add_parser(mode)
        sub.add_argument("--authorization", type=Path, required=True)
        sub.add_argument("--build-dir", type=Path, required=True)
        sub.add_argument("--output-dir", type=Path, required=True)
        if mode == "run":
            sub.add_argument("--small-object", type=Path, required=True)
            sub.add_argument("--envelope-object", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "build":
        receipt = build(args.authorization, args.build_dir, args.output_dir)
    else:
        receipt = run(args.authorization, args.build_dir, args.output_dir,
                      args.small_object, args.envelope_object)
    print(json.dumps(receipt, sort_keys=True))
    return 1 if receipt["status"].startswith("PRODUCER_FAILURE") else 0


if __name__ == "__main__":
    raise SystemExit(main())
