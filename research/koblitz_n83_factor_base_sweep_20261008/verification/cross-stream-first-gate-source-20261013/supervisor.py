#!/usr/bin/env python3
"""Frozen first cross-stream N83 import and one-trial gate; opt-in only.

Preparing or checking this module launches nothing. ``run`` requires a fresh
2,400-second grant tied to the audited 17,400-second local ledger.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time
import uuid


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
STUDY = "koblitz_n83_factor_base_sweep_20261008"
STUDY_PATH = ROOT / "research" / STUDY
sys.path.insert(0, str(STUDY_PATH))
import k600_mixed_orbit_supervisor as prior  # noqa: E402
sys.path.insert(0, str(STUDY_PATH / "verification/mixed-orbit-control-reuse-17400-20261012"))
import reuse_control_supervisor as reuse  # noqa: E402
sys.path.insert(0, str(STUDY_PATH / "verification/paired-capacity-screen-20261013"))
import screen as capacity  # noqa: E402


REL = f"research/{STUDY}/verification/cross-stream-first-gate-source-20261013/supervisor.py"
PILOT_MANIFEST = f"research/{STUDY}/pilot-01/manifest.json"
CORPUS = f"research/{STUDY}/pilot-01/probe-corpus.json"
PILOT_MANIFEST_SHA256 = "5bd862987844f709e03beeab9fd265e3d2024011bcadcb82f5af9c72018df6a6"
BASE_HASH = "b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e"
POOL_HASH = "e5600aa76a88246a4c3f73fce83ba73479e79c34e88d6205fffdb91ce09e0a99"
BASE_SHA256 = "3e85c9e419fcddd3c9faed3c84470bf96fcad9178e48e1d73629369335841ecf"
POOL_MD5_FROM_HEAD = "9e41c76de859046eb295aec39e635c3c"
POOL_BYTES_FROM_HEAD = 3_020_387
POINTS_PER_ORBIT = 166
ORBIT_COUNT = 600
PRIOR_TOTAL = 17_400
BUILD_CAP, VERIFY_CAP, SEARCH_CAP, CLEANUP_RESERVE = 600, 360, 1200, 240
BUILD_WALL, VERIFY_WALL, SEARCH_WALL, WORKER_WALL = 450, 300, 1150, 1100
TOTAL_CHARGE = BUILD_CAP + VERIFY_CAP + SEARCH_CAP + CLEANUP_RESERVE
CONFLICTS, MODELS, BASE_SUMMANDS, RESIDUAL_SUMMANDS = 5000, 32, 3, 2


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def authorization_error(record: object, output: Path) -> str | None:
    if not isinstance(record, dict):
        return "authorization is not an object"
    expected = {
        "schema": "n83.additional-local-budget/v6",
        "study": STUDY,
        "purpose": "one_k600_hash_vs_sequential_cross_stream_partial",
        "prior_total_charged_seconds": str(PRIOR_TOTAL),
        "prior_k600_audit_sha256": reuse.AUDIT_SHA256,
        "base_object_blake3": BASE_HASH,
        "pool_object_blake3": POOL_HASH,
        "additional_authorized_seconds": str(TOTAL_CHARGE),
        "output_dir": str(output.resolve()),
    }
    for key, value in expected.items():
        if record.get(key) != value:
            return f"authorization {key} mismatch"
    if not isinstance(record.get("user_approval_reference"), str) or not record["user_approval_reference"].strip():
        return "authorization lacks explicit approval reference"
    return None


def source_hashes(checkout: Path) -> dict[str, str]:
    paths = (*prior.SOURCE_PATHS, REL, reuse.AUDIT)
    return {relative: prior.sha256_file(checkout / relative) for relative in paths}


def selected_rows(checkout: Path) -> tuple[dict, dict]:
    manifest = checkout / PILOT_MANIFEST
    require(prior.sha256_file(manifest) == PILOT_MANIFEST_SHA256,
            "pilot manifest bytes changed")
    rows = json.loads(manifest.read_text())["bases"]
    found = []
    for digest, policy in ((BASE_HASH, "public_x_hash"),
                           (POOL_HASH, "public_x_sequential")):
        selected = [row for row in rows if row.get("compressed_blake3") == digest]
        require(len(selected) == 1, "pilot object identity is not unique")
        row = selected[0]
        require(row.get("a") == 0 and row.get("policy") == policy
                and row.get("seed") == 2026100801
                and row.get("columns") == ORBIT_COUNT
                and row.get("points") == ORBIT_COUNT * POINTS_PER_ORBIT
                and row.get("object") == f"objects/{digest}.jsonl.gz"
                and row.get("s3_uri") == (
                    f"s3://crypto-autoresearcher/factor-bases/icv1/etc/{STUDY}/a0/"
                    f"objects/{digest}.jsonl.gz"),
                "pilot row identity differs")
        found.append(row)
    require(found[1]["bytes"] == POOL_BYTES_FROM_HEAD,
            "pool manifest length differs from S3 HEAD")
    return found[0], found[1]


def object_error(base: Path, pool: Path, checkout: Path) -> str | None:
    try:
        base_row, pool_row = selected_rows(checkout)
        if prior.sha256_file(base) != BASE_SHA256 or base.stat().st_size != base_row["bytes"]:
            return "retained K600 base compressed bytes differ"
        if pool.stat().st_size != pool_row["bytes"]:
            return "cross-stream pool compressed length differs"
        md5 = hashlib.md5(pool.read_bytes()).hexdigest()
        if md5 != POOL_MD5_FROM_HEAD:
            return "cross-stream pool compressed MD5 differs from observed S3 ETag"
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError):
        return "cross-stream object preflight failed"
    return None


def orbit_difference(base: list[tuple[int, int]], pool: list[tuple[int, int]],
                     orbit_size: int = POINTS_PER_ORBIT) -> tuple[int, list[tuple[int, int]]]:
    require(orbit_size > 0 and len(base) % orbit_size == 0
            and len(pool) % orbit_size == 0, "archive orbit length differs")
    base_set = set(base)
    base_x = {x for x, _ in base}
    require(len(base_set) == len(base), "base repeats a point")
    require(len(set(pool)) == len(pool), "pool repeats a point")
    shared = 0
    residual = []
    for offset in range(0, len(pool), orbit_size):
        orbit = pool[offset:offset + orbit_size]
        matched = sum(point in base_set for point in orbit)
        require(matched in (0, orbit_size), "pool has partial signed-orbit overlap")
        if matched == orbit_size:
            shared += 1
        else:
            require(not any(x in base_x for x, _ in orbit),
                    "residual x-domain intersects base")
            residual.extend(orbit)
    require(residual, "cross-stream residual is empty")
    return shared, residual


def archive_points(path: Path, policy: str) -> list[tuple[int, int]]:
    try:
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            header = json.loads(next(stream))
            require(header.get("schema") == "n83.factor-base/v1"
                    and header.get("curve_slug") == "icv1-f2m83-tm6151469093347-debefd74"
                    and header.get("closure") == "signed_frobenius"
                    and header.get("policy") == policy
                    and header.get("seed") == 2026100801
                    and header.get("orbit_columns") == ORBIT_COUNT
                    and header.get("point_count") == ORBIT_COUNT * POINTS_PER_ORBIT,
                    "archive header differs")
            points = []
            for index, line in enumerate(stream):
                row = json.loads(line)
                point = prior.parse_point(row["point"])
                require(point is not None and row.get("column") == index // POINTS_PER_ORBIT
                        and row.get("phase") == (index % POINTS_PER_ORBIT) // 2
                        and row.get("negative") is (index % 2 == 1),
                        "archive point or orbit label differs")
                points.append(point)
            require(len(points) == ORBIT_COUNT * POINTS_PER_ORBIT,
                    "archive point count differs")
            return points
    except (OSError, EOFError, UnicodeError, StopIteration, KeyError, TypeError,
            json.JSONDecodeError) as error:
        raise ValueError("archive replay failed") from error


def verify_error(record: object, shared: int, residual_orbits: int) -> str | None:
    if not isinstance(record, dict):
        return "cross-stream verify receipt missing"
    expected = {
        "schema": "n83.cross-stream-pilot-import/v1",
        "status": "PASS", "curve_slug": "icv1-f2m83-tm6151469093347-debefd74",
        "pool_object_blake3": POOL_HASH,
        "small_object_blake3": BASE_HASH,
        "small_policy": "public_x_hash", "small_seed": 2026100801,
        "small_columns": ORBIT_COUNT,
        "small_points": ORBIT_COUNT * POINTS_PER_ORBIT,
        "pool_policy": "public_x_sequential", "pool_seed": 2026100801,
        "pool_columns": ORBIT_COUNT,
        "pool_points": ORBIT_COUNT * POINTS_PER_ORBIT,
        "shared_signed_orbits": shared,
        "small_only_signed_orbits": ORBIT_COUNT - shared,
        "residual_signed_orbits": residual_orbits,
        "residual_points": residual_orbits * POINTS_PER_ORBIT,
        "solver_search_executed": False,
    }
    for key, value in expected.items():
        if type(record.get(key)) is not type(value) or record[key] != value:
            return f"cross-stream verify {key} differs"
    return None


def search_error(config: object, result: object, shared: int,
                 base: list[tuple[int, int]], residual: list[tuple[int, int]],
                 public_target: tuple[int, int]) -> str | None:
    if not isinstance(config, dict) or not isinstance(result, dict):
        return "cross-stream search receipt missing"
    residual_orbits = len(residual) // POINTS_PER_ORBIT
    required = capacity.requirements(ORBIT_COUNT, ORBIT_COUNT + residual_orbits, 5, 2)
    if (config.get("schema") != "n83.cross-stream-pilot-partial-config/v1"
            or config.get("base_object_blake3") != BASE_HASH
            or config.get("pool_object_blake3") != POOL_HASH
            or config.get("public_corpus_canonical_json_blake3") != prior.CORPUS_BLAKE3
            or config.get("conflict_budget") != CONFLICTS
            or config.get("max_models") != MODELS
            or config.get("base_summands") != BASE_SUMMANDS
            or config.get("residual_summands") != RESIDUAL_SUMMANDS
            or config.get("budget_seconds") != WORKER_WALL
            or config.get("order_first_pair") is not True
            or config.get("max_variables") != 150000
            or config.get("max_domain_clauses") != 12000000
            or config.get("memory_cgroup_limit_bytes") != 4294967296
            or config.get("memory_cgroup_swap_limit_bytes") != 0
            or result.get("schema") != "n83.cross-stream-pilot-partial-result/v1"
            or result.get("base_object_blake3") != BASE_HASH
            or result.get("pool_object_blake3") != POOL_HASH
            or result.get("base_policy") != "public_x_hash"
            or result.get("base_seed") != 2026100801
            or result.get("pool_policy") != "public_x_sequential"
            or result.get("pool_seed") != 2026100801
            or result.get("public_corpus_canonical_json_blake3") != prior.CORPUS_BLAKE3
            or result.get("base_signed_orbits") != ORBIT_COUNT
            or result.get("pool_signed_orbits") != ORBIT_COUNT
            or result.get("shared_signed_orbits") != shared
            or result.get("residual_signed_orbits") != residual_orbits
            or result.get("trial") != 0
            or (result.get("a"), result.get("b")) != (prior.FIRST_TRIAL_A, prior.FIRST_TRIAL_B)
            or result.get("graph_rows") is not None
            or result.get("matrix_rank") is not None
            or result.get("total_index_calculus_runtime_ms") is not None
            or result.get("selected_best_total_runtime") is not None):
        return "cross-stream search configuration or identity differs"
    if (not prior.finite_nonnegative(result.get("base_import_elapsed_ms_L0"))
            or not prior.finite_nonnegative(result.get("worker_wall_ms_L0"))
            or result["base_import_elapsed_ms_L0"] > result["worker_wall_ms_L0"]
            or result["worker_wall_ms_L0"] > WORKER_WALL * 1000
            or type(result.get("memory_cgroup_peak_bytes")) is not int
            or not 0 < result["memory_cgroup_peak_bytes"] <= 4294967296):
        return "cross-stream worker resource counters differ"
    stats = result.get("search")
    if not isinstance(stats, dict):
        return "cross-stream search counters missing"
    max_patterns = 1 << (BASE_SUMMANDS + RESIDUAL_SUMMANDS - 2)
    for key, cap in (("conflicts", CONFLICTS), ("models", MODELS * max_patterns),
                     ("invalid_models", MODELS * max_patterns),
                     ("group_rejected_models", MODELS * max_patterns),
                     ("solver_calls", max_patterns),
                     ("patterns_attempted", max_patterns)):
        value = stats.get(key)
        if type(value) is not int or not 0 <= value <= cap:
            return f"cross-stream search {key} exceeds cap or is malformed"
    for key in ("exhausted", "refuted"):
        if type(stats.get(key)) is not bool:
            return f"cross-stream search {key} malformed"
    if (stats.get("domain_clauses") != required["domain_clauses"]
            or type(stats.get("peak_variables")) is not int
            or not 0 <= stats["peak_variables"] <= required["total_variables"]
            or type(stats.get("peak_xors")) is not int
            or not 0 <= stats["peak_xors"] <= required["xor_rows"]
            or type(stats.get("peak_clauses")) is not int
            or stats["peak_clauses"] < 0):
        return "cross-stream source capacity counters differ"
    candidate = result.get("candidate")
    status = result.get("status")
    if status not in {"PASS_one_bounded_partial", "PASS_one_trial_refuted",
                      "UNKNOWN_search_cap", "UNKNOWN_preflight_cap",
                      "PRODUCER_FAILURE_invalid_model"}:
        return "cross-stream search status unknown"
    if bool(candidate) != (status == "PASS_one_bounded_partial"):
        return "cross-stream candidate/status mismatch"
    if status == "UNKNOWN_search_cap" and not stats["exhausted"]:
        return "cross-stream search cap lacks exhausted counter"
    if status == "PASS_one_trial_refuted" and (not stats["refuted"] or stats["exhausted"]):
        return "cross-stream refutation lacks refuted counter"
    if status == "UNKNOWN_preflight_cap" and (stats["patterns_attempted"] != 0
                                               or stats["solver_calls"] != 0):
        return "cross-stream preflight cap has solver patterns"
    if status == "PRODUCER_FAILURE_invalid_model" and stats["invalid_models"] == 0:
        return "cross-stream invalid-model status lacks invalid model"
    if candidate is not None:
        if (not isinstance(candidate, dict) or candidate.get("source_id") != 1
                or (candidate.get("a"), candidate.get("b")) !=
                (prior.FIRST_TRIAL_A, prior.FIRST_TRIAL_B)):
            return "cross-stream candidate source differs"
        indices = candidate.get("base_point_indices")
        try:
            points = [prior.parse_point(row) for row in candidate["residual_points"]]
        except (KeyError, TypeError, ValueError):
            return "cross-stream candidate residual malformed"
        if (not isinstance(indices, list) or len(indices) != BASE_SUMMANDS
                or any(type(index) is not int or not 0 <= index < len(base)
                       for index in indices)
                or len(points) != RESIDUAL_SUMMANDS
                or any(point is None or point not in set(residual) for point in points)):
            return "cross-stream candidate point outside exact domains"
        if not prior.replay_equation(candidate["a"], candidate["b"],
                                     [base[index] for index in indices], points,
                                     public_target):
            return "cross-stream candidate group equation fails independent replay"
    return None


def run(authorization: Path, output: Path, base_object: Path, pool_object: Path,
        checkout: Path = ROOT) -> dict:
    checkout = checkout.resolve(strict=True)
    output = prior.outside_fresh_checkout(output, checkout)
    authorization = authorization.resolve(strict=True)
    base_object = base_object.resolve(strict=True)
    pool_object = pool_object.resolve(strict=True)
    raw_auth = authorization.read_bytes()
    auth = json.loads(raw_auth)
    error = authorization_error(auth, output)
    if error:
        raise ValueError(error)
    commit = prior.clean_commit(checkout)
    evidence_error = reuse.prior_evidence_error(checkout)
    if evidence_error:
        raise ValueError(evidence_error)
    object_issue = object_error(base_object, pool_object, checkout)
    if object_issue:
        raise ValueError(object_issue)
    image = prior.image_identity()
    registry = Path(os.environ.get("CARGO_HOME", str(Path.home() / ".cargo"))) / "registry"
    require(registry.is_dir(), "offline Cargo registry unavailable")
    sources = source_hashes(checkout)
    inputs = {str(path): prior.sha256_file(path) for path in
              (authorization, base_object, pool_object, checkout / PILOT_MANIFEST,
               checkout / CORPUS)}
    started = time.monotonic()
    deadline = started + TOTAL_CHARGE
    output.mkdir(parents=True, exist_ok=False)
    build_dir = output / "build"
    run_dir = output / "run"
    build_dir.mkdir()
    run_dir.mkdir()
    (build_dir / "target").mkdir()
    auth_hash = hashlib.sha256(raw_auth).hexdigest()
    prior.write_new_json(output / "config.json", {
        "schema": "n83.cross-stream-first-gate-config/v1", "study": STUDY,
        "source_commit": commit, "source_sha256": sources,
        "authorization_sha256": auth_hash, "input_sha256": inputs,
        "container_image_id": image,
        "base_object_blake3": BASE_HASH, "pool_object_blake3": POOL_HASH,
        "prior_total_charged_seconds": PRIOR_TOTAL,
        "conservative_full_charge_seconds": TOTAL_CHARGE,
        "selected_best_total_runtime": None,
    })
    config_hash = prior.sha256_file(output / "config.json")
    stages = []
    status = "PRODUCER_FAILURE_preflight"
    replay_error = None
    overlap = None
    binary = build_dir / "worker"
    binary_sha = None
    token = uuid.uuid4().hex[:16]
    for name, wall in (("build", BUILD_WALL),
                       ("verify", VERIFY_WALL),
                       ("search", SEARCH_WALL)):
        if name != "build":
            try:
                unchanged = (prior.clean_commit(checkout) == commit
                             and source_hashes(checkout) == sources
                             and prior.sha256_file(binary) == binary_sha)
            except (OSError, ValueError):
                unchanged = False
            if not unchanged:
                status = "PRODUCER_FAILURE_source_or_binary_changed"
                break
        if time.monotonic() + wall + CLEANUP_RESERVE > deadline:
            status = "UNKNOWN_budget_before_stage"
            break
        container_name = f"n83-cross-{name}-{token}"
        if name == "build":
            command = [
                "docker", "run", "--pull", "never", "--platform", "linux/arm64",
                "--name", container_name, "--network", "none", "--cpus", "4",
                "--memory", "8192m", "--memory-swap", "8192m", "--pids-limit", "256",
                "--cap-drop", "ALL", "--security-opt", "no-new-privileges",
                "--user", f"{os.getuid()}:{os.getgid()}",
                "--tmpfs", "/cargo:rw,mode=1777",
                "--mount", f"type=bind,src={registry},dst=/cargo/registry,readonly",
                "--mount", f"type=bind,src={checkout},dst=/src,readonly",
                "--mount", f"type=bind,src={build_dir / 'target'},dst=/target",
                "--env", "CARGO_HOME=/cargo", "--env", "CARGO_TARGET_DIR=/target",
                "--workdir", "/src", image, "cargo", "build", "--release", "--offline",
                "--locked", "--no-default-features", "--example", "koblitz_n83_mixed_orbit_worker",
            ]
            stage_output = build_dir
        else:
            command = [
                "docker", "run", "--pull", "never", "--platform", "linux/arm64",
                "--name", container_name, "--network", "none", "--cpus", "1",
                "--memory", "4096m", "--memory-swap", "4096m", "--pids-limit", "64",
                "--cap-drop", "ALL", "--security-opt", "no-new-privileges",
                "--read-only", "--tmpfs", "/tmp:rw,noexec,nosuid,size=64m",
                "--user", f"{os.getuid()}:{os.getgid()}",
                "--mount", f"type=bind,src={binary},dst=/worker,readonly",
                "--mount", f"type=bind,src={base_object},dst=/base.jsonl.gz,readonly",
                "--mount", f"type=bind,src={pool_object},dst=/pool.jsonl.gz,readonly",
                "--mount", f"type=bind,src={checkout / PILOT_MANIFEST},dst=/pilot-manifest.json,readonly",
                "--mount", f"type=bind,src={checkout / CORPUS},dst=/corpus.json,readonly",
                "--mount", f"type=bind,src={run_dir},dst=/out",
                image, "/worker", "verify-cross-pilot" if name == "verify" else "run-cross-pilot",
                BASE_HASH, "/pilot-manifest.json", "/base.jsonl.gz",
                POOL_HASH, "/pilot-manifest.json", "/pool.jsonl.gz",
            ]
            if name == "search":
                command += ["/corpus.json", str(WORKER_WALL), str(CONFLICTS),
                            str(MODELS), str(BASE_SUMMANDS), str(RESIDUAL_SUMMANDS),
                            "/out/search"]
            stage_output = run_dir
        receipt = prior.stage(name=name, command=command, wall=wall,
                              output=stage_output, container_name=container_name)
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
                cap = json.loads((run_dir / "search/cap.json").read_text()) if name == "search" else None
            except (OSError, UnicodeError, json.JSONDecodeError):
                cap = None
            status = ("UNKNOWN_worker_wall_cap" if isinstance(cap, dict)
                      and cap.get("schema") == "n83.cross-stream-pilot-partial-cap/v1"
                      and cap.get("status") == "UNKNOWN_budget"
                      and cap.get("budget_seconds") == WORKER_WALL
                      and prior.finite_nonnegative(cap.get("elapsed_process_wall_ms"))
                      else "PRODUCER_FAILURE_worker_exit")
            break
        if receipt.get("exit_code") in (134, 137):
            status = "UNKNOWN_resource_cap"
            break
        if (receipt.get("launch_error") or receipt.get("exit_code") != 0
                or not isinstance(state, dict) or state.get("ExitCode") != 0):
            status = "PRODUCER_FAILURE_worker_exit"
            break
        if name == "build":
            artifact = build_dir / "target/release/examples/koblitz_n83_mixed_orbit_worker"
            if not prior.is_linux_elf(artifact):
                status = "PRODUCER_FAILURE_binary"
                break
            try:
                shutil.copyfile(artifact, binary)
                binary.chmod(0o755)
                binary_sha = prior.sha256_file(binary)
            except OSError:
                status = "PRODUCER_FAILURE_binary"
                break
            status = "PASS_build_only"
        elif name == "verify":
            try:
                record = json.loads((run_dir / "verify.stdout").read_text())
                base_points = archive_points(base_object, "public_x_hash")
                pool_points = archive_points(pool_object, "public_x_sequential")
                shared, residual_points = orbit_difference(base_points, pool_points)
                overlap = {"shared_signed_orbits": shared,
                           "residual_signed_orbits": len(residual_points) // POINTS_PER_ORBIT}
                replay_error = verify_error(record, shared, overlap["residual_signed_orbits"])
            except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as error:
                replay_error = f"verify replay raised {type(error).__name__}: {error}"
            if replay_error:
                status = "PRODUCER_FAILURE_import_replay"
                break
            status = "PASS_import_only"
        else:
            try:
                config = json.loads((run_dir / "search/config.json").read_text())
                result = json.loads((run_dir / "search/result.json").read_text())
                corpus = json.loads((checkout / CORPUS).read_text())
                selected = [row for row in corpus["targets"]
                            if row.get("a") == 0 and row.get("fixture") == 0]
                target = prior.parse_point(selected[0]["point"]) if len(selected) == 1 else None
                require(target is not None, "public target missing")
                replay_error = search_error(config, result, shared, base_points,
                                            residual_points, target)
            except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as error:
                replay_error = f"search replay raised {type(error).__name__}: {error}"
            status = "PRODUCER_FAILURE_search_replay" if replay_error else result["status"]
    try:
        frozen_changed = (prior.clean_commit(checkout) != commit
                          or source_hashes(checkout) != sources
                          or any(prior.sha256_file(Path(path)) != digest for path, digest in inputs.items())
                          or prior.sha256_file(output / "config.json") != config_hash)
    except (OSError, ValueError):
        frozen_changed = True
    if frozen_changed:
        status = "PRODUCER_FAILURE_frozen_input_changed"
    elapsed = time.monotonic() - started
    if elapsed > TOTAL_CHARGE:
        status = "PRODUCER_FAILURE_budget_overshoot"
    charge = sum({"build": BUILD_CAP, "verify": VERIFY_CAP, "search": SEARCH_CAP}[row["name"]]
                 for row in stages) + (CLEANUP_RESERVE if stages else 0)
    outer = {
        "schema": "n83.cross-stream-first-gate-outer/v1", "status": status,
        "source_commit": commit, "source_sha256": sources,
        "binary_sha256": binary_sha,
        "authorization_sha256": auth_hash, "config_sha256": config_hash,
        "input_sha256": inputs, "stages": stages,
        "overlap": overlap, "replay_error": replay_error,
        "frozen_input_changed": frozen_changed,
        "elapsed_supervisor_seconds_L0": elapsed,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    prior.write_new_json(output / "outer.json", outer)
    prior.write_new_json(output / "budget.json", {
        "schema": "n83.cross-stream-first-gate-budget/v1",
        "prior_total_charged_seconds": PRIOR_TOTAL,
        "prior_k600_audit_sha256": reuse.AUDIT_SHA256,
        "authorization_sha256": auth_hash,
        "conservative_charge_seconds": charge,
        "cumulative_conservative_charge_seconds": PRIOR_TOTAL + charge,
        "remaining_additional_seconds": TOTAL_CHARGE - charge,
        "round_status": status,
        "selected_best_total_runtime": None,
    })
    return outer


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authorization", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--base-object", required=True, type=Path)
    parser.add_argument("--pool-object", required=True, type=Path)
    args = parser.parse_args()
    result = run(args.authorization, args.output_dir, args.base_object,
                 args.pool_object)
    print(json.dumps(result, sort_keys=True))
    return 1 if result["status"].startswith("PRODUCER_FAILURE") else 0


if __name__ == "__main__":
    raise SystemExit(main())
