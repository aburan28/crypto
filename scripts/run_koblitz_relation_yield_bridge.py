#!/usr/bin/env python3
"""Plan, execute, seal, and verify the Stage-21 Koblitz relation-yield bridge."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import stat
import subprocess
import sys
from typing import Any


REPO = Path(__file__).resolve().parents[1]
HERE = REPO / "research" / "sat_factor_base_review_20260908" / "continuation-05-sota-gates"
DEFAULT_PROTOCOL = HERE / "stage-21-relation-yield-protocol.json"
DEFAULT_METER = REPO / "scripts" / "process_meter.py"
FROZEN_LOCK = HERE / "stage-20-rust-build" / "Cargo.lock"
WORKSPACE_LOCK = REPO / "Cargo.lock"
DISCOVERY_SOURCE = REPO / "examples" / "koblitz_public_factor_base_discovery.rs"
YIELD_SOURCE = REPO / "examples" / "koblitz_relation_yield_bridge.rs"
DISCOVERY_BINARY_NAME = "koblitz_public_factor_base_discovery"
YIELD_BINARY_NAME = "koblitz_relation_yield_bridge"
PROTOCOL_SCHEMA = "koblitz_relation_yield_bridge_protocol.v1"
RESULT_SCHEMA = "koblitz_relation_yield_bridge.v1"
RUN_SEAL_SCHEMA = "koblitz_relation_yield_run_seal.v1"
VERIFY_SCHEMA = "koblitz_relation_yield_verification.v1"
VERIFY_SEAL_SCHEMA = "koblitz_relation_yield_verification_seal.v1"
CONTROLLED_ENV = {
    "CARGO_BUILD_JOBS": "1",
    "RAYON_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}
HEX40 = re.compile(r"[0-9a-f]{40}")
HEX64 = re.compile(r"[0-9a-f]{64}")
PRODUCTION_COUNTS = {"natural": 256, "planted_sat": 64, "proven_unsat": 64}
SMOKE_COUNTS = {"natural": 8, "planted_sat": 4, "proven_unsat": 4}
TASK_NAMES = {
    True: ("00-build", "01-discovery-a0", "02-discovery-a1", "03-yield"),
    False: ("00-build", "01-yield-smoke"),
}
DIRECT_SOURCE_PATHS = (
    "Cargo.toml",
    "examples/koblitz_public_factor_base_discovery.rs",
    "examples/koblitz_relation_yield_bridge.rs",
    "scripts/process_meter.py",
    "scripts/run_koblitz_relation_yield_bridge.py",
    "src/cryptanalysis/koblitz_index_calculus.rs",
    "src/cryptanalysis/koblitz_pdp_phase_a.rs",
    "src/cryptanalysis/semaev_sat.rs",
    "src/cryptanalysis/sat.rs",
)
ROW_FIELDS = {
    "global_ordinal", "arm", "arm_ordinal", "packed_target", "x_hex", "y_hex",
    "x_hamming_weight", "y_hamming_weight", "frobenius_orbit_length",
    "candidate_attempts", "candidate_kind", "candidate_attempts_scope",
    "selection_counter", "selection_priority_blake3", "construction_witness_indices",
    "construction_witness_verified", "exact_pair_table_hit", "verified_witness_indices",
    "point_witness_verified",
}
GENERATION_FIELDS = {
    "candidates", "hash_candidates", "uniform_affine_decode_successes",
    "cofactor_projection_scalar_multiplications", "rejected_affine_decode",
    "rejected_infinity", "rejected_duplicate", "rejected_pair_table_hit",
    "selection_pair_table_lookups",
}
ARM_TARGET_HASH_DOMAINS = {
    "natural": b"koblitz-stage21-natural-targets-v1\0",
    "planted_sat": b"koblitz-stage21-planted-targets-v1\0",
    "proven_unsat": b"koblitz-stage21-proven-unsat-targets-v1\0",
}
WITNESS_HASH_DOMAIN = b"koblitz-stage21-witness-sequence-v1\0"


class Stage21Error(RuntimeError):
    pass


def now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


# Minimal unkeyed BLAKE3 needed to replay the producer's public transcript
# commitments without adding a platform dependency to the verifier.
_B3_IV = (
    0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
    0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19,
)
_B3_PERMUTATION = (2, 6, 3, 10, 7, 0, 4, 13, 1, 11, 12, 5, 9, 14, 15, 8)
_B3_CHUNK_START = 1
_B3_CHUNK_END = 2
_B3_PARENT = 4
_B3_ROOT = 8


def _b3_rotr(value: int, shift: int) -> int:
    return ((value >> shift) | (value << (32 - shift))) & 0xFFFFFFFF


def _b3_g(state: list[int], a: int, b: int, c: int, d: int, x: int, y: int) -> None:
    state[a] = (state[a] + state[b] + x) & 0xFFFFFFFF
    state[d] = _b3_rotr(state[d] ^ state[a], 16)
    state[c] = (state[c] + state[d]) & 0xFFFFFFFF
    state[b] = _b3_rotr(state[b] ^ state[c], 12)
    state[a] = (state[a] + state[b] + y) & 0xFFFFFFFF
    state[d] = _b3_rotr(state[d] ^ state[a], 8)
    state[c] = (state[c] + state[d]) & 0xFFFFFFFF
    state[b] = _b3_rotr(state[b] ^ state[c], 7)


def _b3_round(state: list[int], message: list[int]) -> None:
    _b3_g(state, 0, 4, 8, 12, message[0], message[1])
    _b3_g(state, 1, 5, 9, 13, message[2], message[3])
    _b3_g(state, 2, 6, 10, 14, message[4], message[5])
    _b3_g(state, 3, 7, 11, 15, message[6], message[7])
    _b3_g(state, 0, 5, 10, 15, message[8], message[9])
    _b3_g(state, 1, 6, 11, 12, message[10], message[11])
    _b3_g(state, 2, 7, 8, 13, message[12], message[13])
    _b3_g(state, 3, 4, 9, 14, message[14], message[15])


def _b3_words(block: bytes) -> list[int]:
    padded = block + bytes(64 - len(block))
    return [int.from_bytes(padded[index:index + 4], "little") for index in range(0, 64, 4)]


def _b3_compress(
    chaining_value: tuple[int, ...] | list[int],
    block_words: list[int],
    counter: int,
    block_len: int,
    flags: int,
) -> list[int]:
    state = list(chaining_value) + list(_B3_IV[:4]) + [
        counter & 0xFFFFFFFF,
        (counter >> 32) & 0xFFFFFFFF,
        block_len,
        flags,
    ]
    message = list(block_words)
    for _ in range(7):
        _b3_round(state, message)
        message = [message[index] for index in _B3_PERMUTATION]
    return [state[index] ^ state[index + 8] for index in range(8)] + [
        state[index + 8] ^ chaining_value[index] for index in range(8)
    ]


def _b3_chunk_output(chunk: bytes, counter: int) -> tuple[list[int], list[int], int, int, int]:
    blocks = [chunk[index:index + 64] for index in range(0, len(chunk), 64)] or [b""]
    chaining_value = list(_B3_IV)
    for index, block in enumerate(blocks[:-1]):
        flags = _B3_CHUNK_START if index == 0 else 0
        chaining_value = _b3_compress(
            chaining_value, _b3_words(block), counter, len(block), flags
        )[:8]
    final_index = len(blocks) - 1
    final_flags = _B3_CHUNK_END | (_B3_CHUNK_START if final_index == 0 else 0)
    final = blocks[-1]
    return chaining_value, _b3_words(final), counter, len(final), final_flags


def _b3_output_cv(output: tuple[list[int], list[int], int, int, int]) -> list[int]:
    return _b3_compress(*output)[:8]


def _b3_parent_output(left: list[int], right: list[int]) -> tuple[list[int], list[int], int, int, int]:
    return list(_B3_IV), left + right, 0, 64, _B3_PARENT


def blake3_bytes(data: bytes) -> bytes:
    chunks = [data[index:index + 1024] for index in range(0, len(data), 1024)] or [b""]
    stack: list[list[int]] = []
    for chunk_index, chunk in enumerate(chunks[:-1]):
        chaining_value = _b3_output_cv(_b3_chunk_output(chunk, chunk_index))
        total_chunks = chunk_index + 1
        while total_chunks & 1 == 0:
            chaining_value = _b3_output_cv(
                _b3_parent_output(stack.pop(), chaining_value)
            )
            total_chunks >>= 1
        stack.append(chaining_value)
    output = _b3_chunk_output(chunks[-1], len(chunks) - 1)
    while stack:
        output = _b3_parent_output(stack.pop(), _b3_output_cv(output))
    words = _b3_compress(output[0], output[1], 0, output[3], output[4] | _B3_ROOT)
    return b"".join(word.to_bytes(4, "little") for word in words)[:32]


def blake3_hex(data: bytes) -> str:
    return blake3_bytes(data).hex()


def regular_bytes(path: Path, context: str) -> bytes:
    try:
        metadata = path.lstat()
    except OSError as error:
        raise Stage21Error(f"cannot stat {context} {path}: {error}") from error
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
        raise Stage21Error(f"{context} must be a regular non-symlink file: {path}")
    if metadata.st_nlink != 1:
        raise Stage21Error(f"{context} must not be hard-linked: {path}")
    return path.read_bytes()


def file_identity(path: Path, context: str = "file") -> dict[str, Any]:
    data = regular_bytes(path, context)
    return {"path": str(path.resolve()), "bytes": len(data), "sha256": sha256_bytes(data)}


def read_json(path: Path, context: str) -> tuple[dict[str, Any], bytes]:
    data = regular_bytes(path, context)

    def reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise Stage21Error(f"duplicate JSON key {key!r} in {context}")
            result[key] = value
        return result

    try:
        value = json.loads(data, object_pairs_hook=reject_duplicate_pairs)
    except Stage21Error:
        raise
    except json.JSONDecodeError as error:
        raise Stage21Error(f"invalid JSON in {context}: {error}") from error
    if not isinstance(value, dict):
        raise Stage21Error(f"{context} must contain one JSON object")
    return value, data


def write_new(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(path, flags, 0o644)
    except FileExistsError as error:
        raise Stage21Error(f"refusing to overwrite {path}") from error
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(data)
            output.flush()
            os.fsync(output.fileno())
    except BaseException:
        try:
            path.unlink()
        except OSError:
            pass
        raise


def write_json_new(path: Path, value: Any) -> None:
    write_new(path, json.dumps(value, indent=2, sort_keys=True).encode() + b"\n")


def safe_new_directory(path: Path, context: str) -> Path:
    resolved = path.resolve()
    if (
        resolved == Path("/").resolve()
        or resolved == REPO.resolve()
        or REPO.resolve() in resolved.parents
    ):
        raise Stage21Error(f"unsafe {context} path: {resolved}")
    if path.exists() or path.is_symlink():
        raise Stage21Error(f"{context} must not already exist: {path}")
    return resolved


def safe_new_file(path: Path, context: str) -> Path:
    resolved = path.resolve()
    if resolved == REPO.resolve() or REPO.resolve() in resolved.parents:
        raise Stage21Error(f"{context} must be outside the implementation checkout: {resolved}")
    if path.exists() or path.is_symlink():
        raise Stage21Error(f"{context} must not already exist: {path}")
    parent = resolved.parent
    if not parent.is_dir() or parent.is_symlink():
        raise Stage21Error(f"{context} parent must be an existing real directory: {parent}")
    return resolved


def command_output(command: list[str]) -> dict[str, Any]:
    completed = subprocess.run(command, cwd=REPO, text=True, capture_output=True, check=False)
    return {
        "command": command,
        "returncode": completed.returncode,
        "stdout": completed.stdout.strip(),
        "stderr": completed.stderr.strip(),
    }


def require_hex(value: Any, context: str, pattern: re.Pattern[str] = HEX64) -> str:
    if not isinstance(value, str) or pattern.fullmatch(value) is None:
        raise Stage21Error(f"{context} must be lowercase hexadecimal")
    return value


def tool_identity(path: Path, context: str) -> dict[str, Any]:
    invocation = path.absolute()
    resolved = invocation.resolve(strict=True)
    metadata = resolved.stat()
    if not stat.S_ISREG(metadata.st_mode) or not os.access(resolved, os.X_OK):
        raise Stage21Error(f"{context} must be a regular executable: {resolved}")
    data = resolved.read_bytes()
    return {
        # Preserve the invocation basename. Rustup and similar multicall
        # shims dispatch from argv[0]; resolving `cargo` to `rustup` before
        # launch changes the requested tool even though the bytes are equal.
        "path": str(invocation),
        "bytes": len(data),
        "sha256": sha256_bytes(data),
    }


def safe_child_environment(*, cargo_home: Path | None = None, rustc: Path | None = None) -> dict[str, str]:
    environment = {
        "PATH": "/usr/bin:/bin:/usr/sbin:/sbin",
        "LANG": "C",
        "LC_ALL": "C",
        "TZ": "UTC",
        "TMPDIR": "/tmp",
        **CONTROLLED_ENV,
    }
    if cargo_home is not None:
        environment["HOME"] = str(Path.home())
        environment["CARGO_HOME"] = str(cargo_home.resolve())
        environment["CARGO_INCREMENTAL"] = "0"
    if rustc is not None:
        # Preserve a rustup multicall shim's `rustc` basename. Resolving the
        # path to the underlying `rustup` executable changes its dispatch.
        rustc.resolve(strict=True)
        environment["RUSTC"] = str(rustc.absolute())
    return environment


def git_state() -> dict[str, Any]:
    revision = command_output(["git", "rev-parse", "HEAD"])
    if revision["returncode"] != 0 or HEX40.fullmatch(revision["stdout"]) is None:
        raise Stage21Error("cannot resolve the implementation revision")
    branch = command_output(["git", "branch", "--show-current"])
    status_result = command_output(
        ["git", "status", "--porcelain=v1", "--untracked-files=all", "--", "."]
    )
    if status_result["returncode"] != 0:
        raise Stage21Error("cannot inspect implementation cleanliness")
    tree = command_output(["git", "ls-tree", "-r", "--full-tree", "HEAD"])
    if tree["returncode"] != 0:
        raise Stage21Error("cannot inventory the committed implementation tree")
    return {
        "commit": revision["stdout"],
        "branch": branch["stdout"] or None,
        "dirty": bool(status_result["stdout"]),
        "porcelain": status_result["stdout"].splitlines() if status_result["stdout"] else [],
        "committed_tree_sha256": sha256_bytes((tree["stdout"] + "\n").encode()),
    }


def host_identity(toolchain: dict[str, dict[str, Any]]) -> dict[str, Any]:
    memory_bytes = None
    try:
        memory_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (AttributeError, OSError, ValueError):
        pass
    cpu_model = platform.processor() or None
    if platform.system() == "Darwin":
        cpu = command_output(["sysctl", "-n", "machdep.cpu.brand_string"])
        if cpu["returncode"] == 0 and cpu["stdout"]:
            cpu_model = cpu["stdout"]
        memory = command_output(["sysctl", "-n", "hw.memsize"])
        if memory["returncode"] == 0 and memory["stdout"].isdigit():
            memory_bytes = int(memory["stdout"])
    return {
        "platform": platform.platform(),
        "uname": list(platform.uname()),
        "cpu_model": cpu_model,
        "logical_cpus": os.cpu_count(),
        "physical_memory_bytes": memory_bytes,
        "python": command_output([str(Path(sys.executable).resolve()), "--version"]),
        "cargo": command_output([toolchain["cargo"]["path"], "-V"]),
        "rustc": command_output([toolchain["rustc"]["path"], "-Vv"]),
        "toolchain": toolchain,
        "scientific_environment": safe_child_environment(),
    }


def require_number(value: Any, context: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise Stage21Error(f"{context} must be numeric")
    converted = float(value)
    if not math.isfinite(converted) or (positive and converted <= 0):
        raise Stage21Error(f"{context} is outside its admitted range")
    return converted


def validate_protocol(value: dict[str, Any], production: bool) -> dict[str, Any]:
    if value.get("schema") != PROTOCOL_SCHEMA or value.get("status") != "frozen_before_execution":
        raise Stage21Error("the Stage-21 protocol is not frozen under its expected schema")
    scope = value.get("scope", {})
    required_false = (
        "private_or_imported_targets",
        "private_key_recovery",
        "factor_base_discrete_log_labels",
        "target_scalar_labels",
        "target_subgroup_enumeration",
    )
    if scope.get("public_synthetic_controls_only") is not True or any(
        scope.get(field) is not False for field in required_false
    ):
        raise Stage21Error("the protocol widened its public-synthetic or no-log scope")
    curve = value.get("curve", {})
    if production and curve != {
        "family": "Koblitz",
        "equation": "y^2 + x*y = x^3 + a*x^2 + 1",
        "n": 23,
        "a": 0,
        "group_order": "8383412",
        "subgroup_order": "2095853",
        "cofactor": "4",
        "summands": 2,
    }:
        raise Stage21Error("the production curve parameters changed")
    base = value.get("factor_base", {})
    if production and (
        base.get("dimension") != 12
        or base.get("divisor_indices") != [0, 2]
        or base.get("divisor_polynomial") != 5279
        or base.get("linearised_exponents") != [0, 1, 2, 3, 4, 7, 10, 12]
        or base.get("expected_abscissae") != 4096
        or base.get("expected_rational_points") != 4281
        or base.get("expected_signed_frobenius_orbits_before_projection") != 95
        or base.get("expected_projected_signed_frobenius_columns") != 93
        or base.get("selection_reused_without_target_dependent_changes") is not True
        or base.get("selection_rule")
        != "cofactor admissible; maximize rational points; minimize projected signed-Frobenius columns; lower curve-a variant; lexicographically lower divisor indices"
    ):
        raise Stage21Error("the frozen production factor base changed")
    arms = value.get("target_arms", {})
    expected = PRODUCTION_COUNTS
    if set(arms) != set(expected):
        raise Stage21Error("the protocol must contain exactly the three registered arms")
    if not all(isinstance(arms[name], dict) for name in expected):
        raise Stage21Error("each registered target arm must be an object")
    if production and any(arms[name].get("count") != count for name, count in expected.items()):
        raise Stage21Error("the production 256/64/64 target mix changed")
    if arms["natural"].get("relation_oracle_may_influence_selection") is not False:
        raise Stage21Error("the natural arm may not be selected by relation outcome")
    if arms["natural"].get("scalar_is_constructed_or_recorded") is not False:
        raise Stage21Error("the natural arm may not construct target scalar labels")
    expected_seeds = {
        "natural": "6b6f626c69747a2d737461676532312d6e61747572616c2d7631",
        "planted_sat": "6b6f626c69747a2d737461676532312d706c616e7465642d7631",
        "proven_unsat": "6b6f626c69747a2d737461676532312d756e7361742d7631",
    }
    if any(arms[name].get("seed_hex") != seed for name, seed in expected_seeds.items()):
        raise Stage21Error("the domain-separated target seeds changed")
    if (
        arms["planted_sat"].get("witness_required") is not True
        or arms["planted_sat"].get("point_sum_verification_required") is not True
        or arms["proven_unsat"].get("exact_absence_proof_required") is not True
        or arms["proven_unsat"].get("candidate_screening_cost_charged") is not True
    ):
        raise Stage21Error("the planted or proven-UNSAT control contract changed")
    oracle = value.get("exact_reference_oracle", {})
    if (
        oracle.get("canonical_pair_count_formula") != "F*(F+1)/2"
        or oracle.get("complete_for")
        != "two-summand decomposition over the fully materialized frozen factor base with repetition"
        or oracle.get("all_natural_queries_answered") is not True
        or oracle.get("all_planted_witnesses_reverified") is not True
        or oracle.get("all_proven_unsat_queries_absent") is not True
    ):
        raise Stage21Error("the exact reference-oracle contract changed")
    statistics = value.get("statistics", {})
    if (
        not isinstance(statistics, dict)
        or statistics.get("confidence_interval") != "two-sided 95 percent Wilson score interval"
        or "pseudorandom Bernoulli" not in statistics.get("sampling_assumption", "")
        or statistics.get("without_replacement")
        != "Accepted natural targets are distinct, so sampling is without replacement"
        or not isinstance(statistics.get("finite_population_correction"), str)
        or "multiplicity-weighted" not in statistics.get("planted_multiplicity_boundary", "")
        or statistics.get("no_scaling_fit_from_one_n") is not True
    ):
        raise Stage21Error("the frozen statistical interpretation changed")
    producer = value.get("producer", {})
    if (
        producer.get("source") != "examples/koblitz_relation_yield_bridge.rs"
        or producer.get("example") != YIELD_BINARY_NAME
        or producer.get("production_arguments") != []
        or producer.get("smoke_arguments") != ["--smoke", "8", "4", "4"]
        or producer.get("production_schema") != RESULT_SCHEMA
    ):
        raise Stage21Error("the frozen producer command or result schema changed")
    execution = value.get("execution", {})
    if (
        execution.get("task_order")
        != "one fresh-target locked offline build; two public discovery processes in curve-a order; one relation-yield producer"
        or execution.get("parallel_workers") != 1
        or execution.get("requested_threads_per_process") != 1
        or execution.get("in_stage_retries") != 0
        or execution.get("build_watchdog_seconds") != 1800
        or execution.get("build_uses_fresh_target_directory") is not True
        or execution.get("build_offline_and_locked") is not True
        or execution.get("cargo_incremental_disabled") is not True
        or execution.get("discovery_watchdog_seconds") != 300
        or execution.get("producer_watchdog_seconds") != 7200
        or execution.get("whole_driver_watchdog_seconds") != 10800
        or execution.get("process_meter") != "scripts/process_meter.py"
        or execution.get("process_meter_exclusive_create") is not True
        or execution.get("whole_driver_outer_meter_required") is not True
        or execution.get("write_once_output") is not True
        or execution.get("failed_and_capped_processes_retained") is not True
    ):
        raise Stage21Error("the single-worker write-once execution policy changed")
    acceptance = value.get("production_acceptance", {})
    if production and (
        acceptance.get("natural_targets") != 256
        or acceptance.get("planted_sat_targets") != 64
        or acceptance.get("proven_unsat_targets") != 64
        or acceptance.get("planted_sat_hits") != 64
        or acceptance.get("proven_unsat_hits") != 0
        or acceptance.get("all_targets_unique_across_arms") is not True
        or acceptance.get("natural_arm_never_selected_by_oracle_result") is not True
        or acceptance.get("factor_base_logs_absent") is not True
        or acceptance.get("target_scalar_labels_absent") is not True
        or acceptance.get("protocol_path_and_hash_bound") is not True
        or acceptance.get("dependency_lock_bound") is not True
        or acceptance.get("complete_process_and_artifact_inventory") is not True
        or acceptance.get("terminal_status") != "relation_yield_outputs_frozen"
        or acceptance.get("scientific_measurement_admitted_by_control_plane") is not False
        or acceptance.get("pending_independent_payload_replay") is not True
        or acceptance.get("external_portable_verification_satisfied") is not False
    ):
        raise Stage21Error("the production acceptance gate changed")
    return value


def load_protocol(path: Path, production: bool) -> tuple[dict[str, Any], bytes, Path]:
    resolved = path.resolve(strict=True)
    if production and resolved != DEFAULT_PROTOCOL.resolve(strict=True):
        raise Stage21Error("production requires the committed default Stage-21 protocol path")
    value, data = read_json(resolved, "Stage-21 protocol")
    return validate_protocol(value, production), data, resolved


def install_frozen_lock() -> dict[str, Any]:
    frozen = regular_bytes(FROZEN_LOCK, "frozen dependency lock")
    if WORKSPACE_LOCK.exists() or WORKSPACE_LOCK.is_symlink():
        current = regular_bytes(WORKSPACE_LOCK, "workspace dependency lock")
        if current != frozen:
            raise Stage21Error("workspace Cargo.lock differs from the frozen lock")
        installed = False
    else:
        write_new(WORKSPACE_LOCK, frozen)
        installed = True
    return {
        "frozen": file_identity(FROZEN_LOCK, "frozen dependency lock"),
        "workspace": file_identity(WORKSPACE_LOCK, "workspace dependency lock"),
        "installed_by_runner": installed,
    }


def process_metrics(path: Path, context: str) -> dict[str, Any]:
    value, _ = read_json(path, context)
    if set(value) != {
        "command",
        "returncode",
        "watchdog_seconds",
        "timed_out",
        "orphan_group_terminated",
        "metrics",
    }:
        raise Stage21Error(f"{context} uses an unexpected schema")
    if not isinstance(value["command"], list) or not all(isinstance(x, str) for x in value["command"]):
        raise Stage21Error(f"{context}.command must be a string list")
    require_number(value["watchdog_seconds"], f"{context}.watchdog_seconds", positive=True)
    if not isinstance(value["returncode"], int) or isinstance(value["returncode"], bool):
        raise Stage21Error(f"{context}.returncode must be an integer")
    if not isinstance(value["timed_out"], bool) or not isinstance(value["orphan_group_terminated"], bool):
        raise Stage21Error(f"{context} termination flags must be Boolean")
    metrics = value["metrics"]
    if not isinstance(metrics, dict) or set(metrics) != {
        "wall_seconds",
        "user_seconds",
        "system_seconds",
        "total_core_seconds",
        "single_core_seconds",
        "peak_rss_bytes",
        "meter",
    }:
        raise Stage21Error(f"{context}.metrics uses an unexpected schema")
    wall = require_number(metrics["wall_seconds"], f"{context}.wall_seconds", positive=True)
    user = require_number(metrics["user_seconds"], f"{context}.user_seconds")
    system = require_number(metrics["system_seconds"], f"{context}.system_seconds")
    core = require_number(metrics["total_core_seconds"], f"{context}.total_core_seconds")
    alias = require_number(metrics["single_core_seconds"], f"{context}.single_core_seconds")
    if user < 0 or system < 0 or core < 0 or not math.isclose(core, user + system, abs_tol=1e-9):
        raise Stage21Error(f"{context} CPU fields are inconsistent")
    if alias != core:
        raise Stage21Error(f"{context} legacy single_core_seconds alias changed")
    rss = metrics["peak_rss_bytes"]
    if isinstance(rss, bool) or not isinstance(rss, int) or rss < 0:
        raise Stage21Error(f"{context}.peak_rss_bytes must be a nonnegative integer")
    if metrics["meter"] != "fresh-process getrusage(RUSAGE_CHILDREN)":
        raise Stage21Error(f"{context} meter identity changed")
    _ = wall
    return value


def validate_file_identity(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != {"path", "bytes", "sha256"}:
        raise Stage21Error(f"{context} uses an unexpected identity schema")
    if not isinstance(value["path"], str) or not Path(value["path"]).is_absolute():
        raise Stage21Error(f"{context}.path must be absolute")
    if isinstance(value["bytes"], bool) or not isinstance(value["bytes"], int) or value["bytes"] < 0:
        raise Stage21Error(f"{context}.bytes must be a nonnegative integer")
    require_hex(value["sha256"], f"{context}.sha256")
    actual = file_identity(Path(value["path"]), context)
    if actual != value:
        raise Stage21Error(f"{context} changed after its identity was recorded")
    return value


def validate_tool_identity(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != {"path", "bytes", "sha256"}:
        raise Stage21Error(f"{context} uses an unexpected tool identity schema")
    if not isinstance(value["path"], str) or not Path(value["path"]).is_absolute():
        raise Stage21Error(f"{context}.path must be absolute")
    if isinstance(value["bytes"], bool) or not isinstance(value["bytes"], int) or value["bytes"] <= 0:
        raise Stage21Error(f"{context}.bytes must be positive")
    require_hex(value["sha256"], f"{context}.sha256")
    actual = tool_identity(Path(value["path"]), context)
    if actual != value:
        raise Stage21Error(f"{context} changed after its identity was recorded")
    return value


def meter_command(
    *, meter: Path, cwd: Path, timeout: float, stdout: Path, stderr: Path,
    metrics: Path, command: list[str]
) -> list[str]:
    return [
        str(Path(sys.executable).resolve()),
        str(meter.resolve()),
        "--cwd",
        str(cwd.resolve()),
        "--timeout",
        str(float(timeout)),
        "--stdout",
        str(stdout.resolve()),
        "--stderr",
        str(stderr.resolve()),
        "--metrics",
        str(metrics.resolve()),
        "--exclusive-create",
        "--",
        *command,
    ]


def run_metered(
    root: Path, name: str, command: list[str], timeout: float, inputs: list[dict[str, Any]],
    meter: Path, environment: dict[str, str], executable: dict[str, Any]
) -> tuple[dict[str, Any], Path]:
    task = root / "tasks" / name
    task.mkdir(parents=True)
    stdout = task / "stdout"
    stderr = task / "stderr"
    metrics = task / "metrics.json"
    if not command or command[0] != executable.get("path"):
        raise Stage21Error(f"{name} command does not use its bound executable")
    validate_tool_identity(executable, f"{name} executable")
    for index, input_identity in enumerate(inputs):
        validate_file_identity(input_identity, f"{name} input {index}")
    meter_identity = file_identity(meter, "process meter")
    intent = {
        "schema": "koblitz_relation_yield_process_intent.v1",
        "name": name,
        "command": command,
        "cwd": str(REPO.resolve()),
        "watchdog_seconds": timeout,
        "environment": environment,
        "inputs": inputs,
        "meter": meter_identity,
        "meter_exclusive_create": True,
        "executable": executable,
        "stdin": "devnull",
        "close_fds": True,
    }
    write_json_new(task / "intent.json", intent)
    wrapper = meter_command(
        meter=meter,
        cwd=REPO,
        timeout=timeout,
        stdout=stdout,
        stderr=stderr,
        metrics=metrics,
        command=command,
    )
    completed = subprocess.run(
        wrapper,
        cwd=REPO,
        env=environment,
        stdin=subprocess.DEVNULL,
        close_fds=True,
        check=False,
    )
    if completed.returncode != 0:
        raise Stage21Error(f"process meter failed for {name} with status {completed.returncode}")
    receipt = process_metrics(metrics, f"{name} metrics")
    if receipt["command"] != command or receipt["watchdog_seconds"] != timeout:
        raise Stage21Error(f"{name} metrics differ from the frozen command or watchdog")
    if receipt["returncode"] != 0 or receipt["timed_out"] or receipt["orphan_group_terminated"]:
        raise Stage21Error(
            f"{name} did not complete cleanly: returncode={receipt['returncode']} "
            f"timed_out={receipt['timed_out']} orphan={receipt['orphan_group_terminated']}"
        )
    validate_tool_identity(executable, f"{name} executable")
    if file_identity(meter, "process meter") != meter_identity:
        raise Stage21Error(f"{name} process meter changed during execution")
    for index, input_identity in enumerate(inputs):
        validate_file_identity(input_identity, f"{name} input {index}")
    task_receipt = {
        "schema": "koblitz_relation_yield_process_receipt.v1",
        "name": name,
        "process": receipt,
        "intent": file_identity(task / "intent.json", f"{name} intent"),
        "stdout": file_identity(stdout, f"{name} stdout"),
        "stderr": file_identity(stderr, f"{name} stderr"),
        "metrics": file_identity(metrics, f"{name} metrics"),
        "executable": executable,
        "inputs": inputs,
        "environment": environment,
    }
    write_json_new(task / "receipt.json", task_receipt)
    return task_receipt, stdout


def parse_stdout_json(path: Path, context: str) -> dict[str, Any]:
    value, _ = read_json(path, context)
    return value


def validate_discovery(result: dict[str, Any], curve_a: int) -> None:
    if (
        result.get("schema") != "koblitz_public_factor_base_discovery.v1"
        or result.get("n") != 23
        or result.get("a") != curve_a
        or result.get("m") != 2
        or result.get("requested_dimension") != 12
    ):
        raise Stage21Error(f"curve-a={curve_a} discovery output changed shape or parameters")
    forbidden = result.get("forbidden_inputs", {})
    if (
        not isinstance(forbidden, dict)
        or set(forbidden)
        != {
            "target_constructed",
            "target_subgroup_enumerated",
            "discrete_log_labels_constructed",
            "relation_yield_used",
            "solver_timing_used",
        }
        or any(value is not False for value in forbidden.values())
    ):
        raise Stage21Error(f"curve-a={curve_a} discovery consumed a forbidden input")
    if not isinstance(result.get("candidates"), list) or not result["candidates"]:
        raise Stage21Error(f"curve-a={curve_a} discovery emitted no candidates")
    for candidate in result["candidates"]:
        if not isinstance(candidate, dict):
            raise Stage21Error(f"curve-a={curve_a} discovery candidate is not an object")
        for field in (
            "dimension",
            "divisor_polynomial",
            "abscissae",
            "rational_points",
            "signed_frobenius_orbits_before_projection",
            "projected_signed_frobenius_orbits",
        ):
            require_nonnegative_int(candidate.get(field), f"curve-a={curve_a}.{field}")
        if (
            not isinstance(candidate.get("divisor_indices"), list)
            or not candidate["divisor_indices"]
            or not all(isinstance(index, int) and not isinstance(index, bool) and index >= 0 for index in candidate["divisor_indices"])
            or not isinstance(candidate.get("linearised_exponents"), list)
            or not isinstance(candidate.get("m_cofactor_admissible"), bool)
        ):
            raise Stage21Error(f"curve-a={curve_a} discovery candidate fields are invalid")


def cross_curve_winner(discoveries: list[dict[str, Any]]) -> tuple[int, dict[str, Any]]:
    candidates: list[tuple[int, dict[str, Any]]] = []
    for discovery in discoveries:
        for candidate in discovery["candidates"]:
            if candidate.get("m_cofactor_admissible") is True:
                candidates.append((discovery["a"], candidate))
    if not candidates:
        raise Stage21Error("public discovery produced no cofactor-admissible candidate")
    candidates.sort(
        key=lambda item: (
            -item[1]["rational_points"],
            item[1]["projected_signed_frobenius_orbits"],
            item[0],
            item[1]["divisor_indices"],
        )
    )
    return candidates[0]


def wilson_interval_95(hits: int, count: int) -> tuple[float, float]:
    z = 1.959963984540054
    proportion = hits / count
    z2 = z * z
    denominator = 1.0 + z2 / count
    center = (proportion + z2 / (2.0 * count)) / denominator
    margin = z * math.sqrt(
        proportion * (1.0 - proportion) / count + z2 / (4.0 * count * count)
    ) / denominator
    return max(0.0, center - margin), min(1.0, center + margin)


def require_nonnegative_int(value: Any, context: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise Stage21Error(f"{context} must be a nonnegative integer")
    return value


def require_exact_object(value: Any, fields: set[str], context: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != fields:
        raise Stage21Error(f"{context} uses an unexpected schema")
    return value


def json_blake3(value: Any) -> str:
    return blake3_hex(canonical_bytes(value))


def validate_yield_result(
    result: dict[str, Any], production: bool, protocol: dict[str, Any]
) -> dict[str, Any]:
    require_exact_object(
        result,
        {
            "schema", "status", "evidence_class", "production_defaults_used",
            "ledger_promotion_eligible_from_this_output_alone", "measurement_schema",
            "frozen_policy", "factor_base", "exact_pair_table", "arms",
            "ordered_covariate_rows", "ordered_covariate_rows_blake3", "timing_ns",
            "high_level_operation_counts", "retained_size_lower_bounds_bytes", "hashes",
            "resource_accounting_boundary", "claim_boundary", "non_claims",
        },
        "yield producer result",
    )
    if result.get("schema") != RESULT_SCHEMA or result.get("status") != "complete":
        raise Stage21Error("yield producer did not emit a complete supported result")
    if result.get("production_defaults_used") is not production:
        raise Stage21Error("yield producer execution mode differs from the requested mode")
    expected_evidence = (
        "finite_public_synthetic_measurement_pending_independent_replay_and_external_resource_receipt"
        if production
        else "operational_smoke_only_ineligible_for_ledger_promotion"
    )
    if (
        result.get("evidence_class") != expected_evidence
        or result.get("ledger_promotion_eligible_from_this_output_alone") is not False
    ):
        raise Stage21Error("yield producer changed its evidence boundary")
    measurement = result.get("measurement_schema", {})
    require_exact_object(
        measurement,
        {
            "stage", "n_or_bits", "base_id_or_hash", "eta_or_coverage_policy",
            "pr_decomposition_or_hit_rate_with_ci", "trials_per_relation", "target_mix",
        },
        "yield producer measurement schema",
    )
    expected_mix = PRODUCTION_COUNTS if production else SMOKE_COUNTS
    expected_n = 23 if production else 7
    if measurement.get("stage") != "relation_yield" or measurement.get("n_or_bits") != expected_n:
        raise Stage21Error("yield producer omitted the relation-yield measurement schema")
    mix = measurement.get("target_mix", {})
    if mix != expected_mix:
        raise Stage21Error("yield producer omitted or changed a registered arm")
    arms = result.get("arms")
    if not isinstance(arms, dict) or set(arms) != {"natural", "planted_sat", "proven_unsat"}:
        raise Stage21Error("yield producer must retain exactly all three arms")
    for name, expected_count in expected_mix.items():
        arm = arms[name]
        require_exact_object(
            arm,
            {
                "count", "hits", "misses", "hit_rate", "wilson_95_ci",
                "targets_blake3", "witness_results_blake3",
                "verified_point_witnesses", "timing_ns", "high_level_operations",
            },
            f"{name} arm",
        )
        count = require_nonnegative_int(arm.get("count"), f"{name}.count")
        hits = require_nonnegative_int(arm.get("hits"), f"{name}.hits")
        misses = require_nonnegative_int(arm.get("misses"), f"{name}.misses")
        witnesses = require_nonnegative_int(
            arm.get("verified_point_witnesses"), f"{name}.verified_point_witnesses"
        )
        if count != expected_count or hits + misses != count or witnesses != hits:
            raise Stage21Error(f"{name} arm counts or point-witness accounting changed")
        rate = require_number(arm.get("hit_rate"), f"{name}.hit_rate")
        if not math.isclose(rate, hits / count, rel_tol=0, abs_tol=1e-15):
            raise Stage21Error(f"{name} arm hit rate differs from its counts")
        interval = arm.get("wilson_95_ci", {})
        require_exact_object(interval, {"low", "high"}, f"{name} arm Wilson interval")
        low = require_number(interval.get("low"), f"{name}.wilson.low")
        high = require_number(interval.get("high"), f"{name}.wilson.high")
        expected_low, expected_high = wilson_interval_95(hits, count)
        if (
            not math.isclose(low, expected_low, rel_tol=0, abs_tol=1e-14)
            or not math.isclose(high, expected_high, rel_tol=0, abs_tol=1e-14)
        ):
            raise Stage21Error(f"{name} arm Wilson interval differs from its counts")
        require_hex(arm.get("targets_blake3"), f"{name}.targets_blake3")
        require_hex(arm.get("witness_results_blake3"), f"{name}.witness_results_blake3")
        arm_timing = require_exact_object(
            arm.get("timing_ns"),
            {"target_generation_and_selection", "exact_pair_table_queries", "point_witness_verification"},
            f"{name} arm timing",
        )
        for timer_name, timer in arm_timing.items():
            require_nonnegative_int(timer, f"{name}.timing_ns.{timer_name}")
        operations = require_exact_object(
            arm.get("high_level_operations"),
            {"generation", "final_pair_table_lookups", "point_witness_verification_additions"},
            f"{name} arm operations",
        )
        generation = require_exact_object(
            operations.get("generation"), GENERATION_FIELDS, f"{name} generation operations"
        )
        generation = {
            field: require_nonnegative_int(value, f"{name}.generation.{field}")
            for field, value in generation.items()
        }
        if (
            operations.get("final_pair_table_lookups") != count
            or operations.get("point_witness_verification_additions") != hits
        ):
            raise Stage21Error(f"{name} query or witness operation counts changed")
        if name == "planted_sat":
            if generation["candidates"] <= 0 or any(
                generation[field] != 0 for field in GENERATION_FIELDS - {"candidates"}
            ):
                raise Stage21Error("planted generation counters changed their pair-scan meaning")
        elif (
            generation["candidates"] != generation["hash_candidates"]
            or generation["candidates"]
            != generation["rejected_affine_decode"] + generation["uniform_affine_decode_successes"]
            or generation["cofactor_projection_scalar_multiplications"]
            != generation["uniform_affine_decode_successes"]
            or generation["uniform_affine_decode_successes"]
            != count + generation["rejected_infinity"] + generation["rejected_duplicate"]
            + generation["rejected_pair_table_hit"]
        ):
            raise Stage21Error(f"{name} target-generation counters do not partition the attempts")
        if name == "natural" and (
            generation["rejected_pair_table_hit"] != 0
            or generation["selection_pair_table_lookups"] != 0
        ):
            raise Stage21Error("natural target generation consulted the pair oracle")
        if name == "proven_unsat" and generation["selection_pair_table_lookups"] != (
            count + generation["rejected_pair_table_hit"]
        ):
            raise Stage21Error("proven-UNSAT screening lookup count changed")
    planted = arms["planted_sat"]
    unsat = arms["proven_unsat"]
    if planted.get("hits") != planted.get("count") or planted.get("verified_point_witnesses") != planted.get("count"):
        raise Stage21Error("the planted arm failed exact witness validation")
    if unsat.get("hits") != 0 or unsat.get("verified_point_witnesses") != 0:
        raise Stage21Error("the proven-UNSAT arm contains a decomposition")
    natural = arms["natural"]
    expected_trials: float | str = (
        natural["count"] / natural["hits"] if natural["hits"] else "infinity"
    )
    observed_trials = measurement.get("trials_per_relation")
    if isinstance(expected_trials, str):
        if observed_trials != expected_trials:
            raise Stage21Error("zero-hit natural arm changed its trials-per-relation label")
    elif isinstance(observed_trials, bool) or not isinstance(observed_trials, (int, float)) or not math.isclose(
        float(observed_trials), expected_trials, rel_tol=0, abs_tol=1e-15
    ):
        raise Stage21Error("natural trials per relation differs from its hit count")
    primary = require_exact_object(
        measurement.get("pr_decomposition_or_hit_rate_with_ci"),
        {"natural_hit_rate", "wilson_95_ci"},
        "natural primary estimate",
    )
    primary_interval = require_exact_object(
        primary.get("wilson_95_ci"), {"low", "high"}, "natural primary interval"
    )
    if (
        primary.get("natural_hit_rate") != natural["hit_rate"]
        or primary_interval != natural["wilson_95_ci"]
    ):
        raise Stage21Error("natural primary estimate differs from the natural arm")
    policy = result.get("frozen_policy", {})
    require_exact_object(
        policy,
        {
            "schema", "n", "a", "m", "irreducible_low_terms", "group_order",
            "subgroup_order", "cofactor", "divisor_indices", "target_mix", "seed_hex",
            "point_encoding", "hash_to_curve", "selection", "pair_priority", "oracle",
        },
        "yield producer frozen policy",
    )
    hash_to_curve = require_exact_object(
        policy.get("hash_to_curve"),
        {"hash", "domain_hex", "message", "x_draw", "lift", "projection", "rejections", "target_scalar_constructed_or_recorded"},
        "hash-to-curve policy",
    )
    selection = require_exact_object(
        policy.get("selection"),
        {"natural", "planted_sat", "proven_unsat", "arms_disjoint", "hash_arm_attempt_cap_formula"},
        "target-selection policy",
    )
    pair_priority = require_exact_object(
        policy.get("pair_priority"), {"hash", "domain_hex", "message"}, "pair-priority policy"
    )
    expected_curve = protocol["curve"]
    if (
        policy.get("schema") != "koblitz_relation_yield_policy.v1"
        or policy.get("n") != expected_n
        or policy.get("a") != (0 if production else 1)
        or policy.get("m") != 2
        or policy.get("divisor_indices") != [0, 2]
        or policy.get("target_mix") != expected_mix
        or policy.get("seed_hex")
        != {name: protocol["target_arms"][name]["seed_hex"] for name in expected_mix}
        or hash_to_curve.get("target_scalar_constructed_or_recorded") is not False
        or hash_to_curve.get("hash") != "BLAKE3"
        or selection.get("hash_arm_attempt_cap_formula") != "10000*requested + 10000"
        or pair_priority.get("hash") != "BLAKE3"
        or selection.get("arms_disjoint") is not True
        or selection.get("natural")
        != "first distinct nonidentity hash-to-curve/cofactor targets; pair oracle unavailable and unused"
    ):
        raise Stage21Error("yield producer changed its frozen target policy or no-scalar boundary")
    if production and (
        policy.get("group_order") != expected_curve["group_order"]
        or policy.get("subgroup_order") != expected_curve["subgroup_order"]
        or policy.get("cofactor") != expected_curve["cofactor"]
    ):
        raise Stage21Error("yield producer changed the production group parameters")
    factor_base = result.get("factor_base")
    require_exact_object(
        factor_base,
        {"predicate", "predicate_blake3", "factor_base_blake3", "point_order", "selection_and_construction_boundaries"},
        "yield producer factor-base record",
    )
    predicate = factor_base.get("predicate", {})
    require_exact_object(
        predicate,
        {
            "curve", "group_order", "subgroup_order", "cofactor", "m",
            "construction_method", "divisor_indices", "divisor_polynomial_bitmask",
            "linearised_exponents", "dimension", "abscissae", "rational_points",
            "signed_frobenius_orbits_before_projection", "projected_signed_frobenius_columns",
        },
        "yield producer factor-base predicate",
    )
    predicate_curve = require_exact_object(
        predicate.get("curve"), {"n", "a", "irreducible_low_terms"}, "factor-base curve"
    )
    if production and (
        predicate_curve.get("n") != 23
        or predicate_curve.get("a") != 0
        or predicate.get("group_order") != expected_curve["group_order"]
        or predicate.get("subgroup_order") != expected_curve["subgroup_order"]
        or predicate.get("cofactor") != expected_curve["cofactor"]
        or predicate.get("m") != 2
        or predicate.get("divisor_indices") != protocol["factor_base"]["divisor_indices"]
        or predicate.get("divisor_polynomial_bitmask") != protocol["factor_base"]["divisor_polynomial"]
        or predicate.get("linearised_exponents") != protocol["factor_base"]["linearised_exponents"]
        or predicate.get("dimension") != protocol["factor_base"]["dimension"]
        or predicate.get("abscissae") != protocol["factor_base"]["expected_abscissae"]
        or predicate.get("rational_points") != protocol["factor_base"]["expected_rational_points"]
        or predicate.get("signed_frobenius_orbits_before_projection")
        != protocol["factor_base"]["expected_signed_frobenius_orbits_before_projection"]
        or predicate.get("projected_signed_frobenius_columns")
        != protocol["factor_base"]["expected_projected_signed_frobenius_columns"]
    ):
        raise Stage21Error("yield producer changed the frozen algebraic factor-base identity")
    boundaries = factor_base.get("selection_and_construction_boundaries", {})
    require_exact_object(
        boundaries,
        {
            "public_field_and_curve_parameters_only", "target_available_during_selection",
            "target_subgroup_enumerated_for_factor_base",
            "factor_base_discrete_log_labels_constructed",
            "target_discrete_log_labels_constructed", "relation_yield_used_for_selection",
            "solver_timing_used_for_selection",
        },
        "yield producer factor-base boundaries",
    )
    if (
        boundaries.get("public_field_and_curve_parameters_only") is not True
        or boundaries.get("target_available_during_selection") is not False
        or boundaries.get("target_subgroup_enumerated_for_factor_base") is not False
        or boundaries.get("factor_base_discrete_log_labels_constructed") is not False
        or boundaries.get("target_discrete_log_labels_constructed") is not False
        or boundaries.get("relation_yield_used_for_selection") is not False
        or boundaries.get("solver_timing_used_for_selection") is not False
    ):
        raise Stage21Error("factor-base construction crossed a forbidden boundary")
    pair_table = result.get("exact_pair_table", {})
    require_exact_object(
        pair_table,
        {
            "summands", "repeated_factor_points_allowed", "canonical_pair_policy",
            "canonical_pairs", "enumerated_pairs", "unique_target_entries",
            "hash_table_capacity", "duplicate_pair_sums", "packed_key",
            "canonical_pair_transcript_blake3", "unsat_proof_boundary",
        },
        "yield producer pair-table record",
    )
    factor_points = predicate.get("rational_points")
    canonical_pairs = require_nonnegative_int(
        pair_table.get("canonical_pairs"), "exact_pair_table.canonical_pairs"
    )
    enumerated_pairs = require_nonnegative_int(
        pair_table.get("enumerated_pairs"), "exact_pair_table.enumerated_pairs"
    )
    unique_targets = require_nonnegative_int(
        pair_table.get("unique_target_entries"), "exact_pair_table.unique_target_entries"
    )
    duplicate_sums = require_nonnegative_int(
        pair_table.get("duplicate_pair_sums"), "exact_pair_table.duplicate_pair_sums"
    )
    table_capacity = require_nonnegative_int(
        pair_table.get("hash_table_capacity"), "exact_pair_table.hash_table_capacity"
    )
    if (
        isinstance(factor_points, bool)
        or not isinstance(factor_points, int)
        or canonical_pairs != factor_points * (factor_points + 1) // 2
        or enumerated_pairs != canonical_pairs
        or unique_targets > canonical_pairs
        or table_capacity < unique_targets
        or duplicate_sums != canonical_pairs - unique_targets
        or pair_table.get("summands") != 2
        or pair_table.get("repeated_factor_points_allowed") is not True
        or pair_table.get("canonical_pair_policy") != "all indices i<=j exactly once"
    ):
        raise Stage21Error("exact canonical pair-table completeness changed")
    rows = result.get("ordered_covariate_rows")
    if not isinstance(rows, list) or len(rows) != sum(expected_mix.values()):
        raise Stage21Error("ordered covariate rows do not cover the registered targets")
    expected_arms = [
        arm for arm in ("natural", "planted_sat", "proven_unsat")
        for _ in range(expected_mix[arm])
    ]
    packed_targets: set[int] = set()
    arm_ordinals = {name: 0 for name in expected_mix}
    row_hits = {name: 0 for name in expected_mix}
    orbit_applications = 0
    coordinate_width = (expected_n + 3) // 4
    for ordinal, (row, expected_arm) in enumerate(zip(rows, expected_arms, strict=True)):
        require_exact_object(row, ROW_FIELDS, f"ordered covariate row {ordinal}")
        if row.get("global_ordinal") != ordinal or row.get("arm") != expected_arm:
            raise Stage21Error("ordered covariate row order or arm changed")
        if row.get("arm_ordinal") != arm_ordinals[expected_arm]:
            raise Stage21Error("ordered covariate arm ordinal changed")
        arm_ordinals[expected_arm] += 1
        packed = row.get("packed_target")
        if isinstance(packed, bool) or not isinstance(packed, int) or packed <= 0 or packed in packed_targets:
            raise Stage21Error("target arms contain an invalid or duplicate exact point")
        packed_targets.add(packed)
        if any(key in row for key in ("target_scalar", "target_discrete_log", "known_log")):
            raise Stage21Error("a covariate row exposed a target scalar label")
        encoded = packed - 1
        coordinate_mask = (1 << expected_n) - 1
        x = encoded >> expected_n
        y = encoded & coordinate_mask
        if (
            x > coordinate_mask
            or row.get("x_hex") != f"{x:0{coordinate_width}x}"
            or row.get("y_hex") != f"{y:0{coordinate_width}x}"
            or row.get("x_hamming_weight") != x.bit_count()
            or row.get("y_hamming_weight") != y.bit_count()
        ):
            raise Stage21Error("row coordinates or Hamming covariates differ from the packed target")
        orbit = require_nonnegative_int(
            row.get("frobenius_orbit_length"), f"row {ordinal} Frobenius orbit length"
        )
        if orbit == 0 or expected_n % orbit != 0:
            raise Stage21Error("row Frobenius orbit length is incompatible with the field degree")
        orbit_applications += orbit
        attempts = require_nonnegative_int(row.get("candidate_attempts"), f"row {ordinal} candidate attempts")
        counter = require_nonnegative_int(row.get("selection_counter"), f"row {ordinal} selection counter")
        if attempts == 0:
            raise Stage21Error("row candidate-attempt count must be positive")
        is_planted = expected_arm == "planted_sat"
        expected_kind = (
            "complete_canonical_pair_priority_scan"
            if is_planted else "domain_separated_blake3_uniform_affine_then_cofactor_projection"
        )
        expected_scope = (
            "one-based canonical-pair scan position of the selected construction witness"
            if is_planted else "hash candidates since the preceding accepted target in this arm"
        )
        construction = row.get("construction_witness_indices")
        verified = row.get("verified_witness_indices")
        hit = row.get("exact_pair_table_hit")
        if not isinstance(hit, bool):
            raise Stage21Error("row exact-pair-table result must be Boolean")
        if hit:
            row_hits[expected_arm] += 1
        for label, witness in (("construction", construction), ("verified", verified)):
            if witness is not None and (
                not isinstance(witness, list)
                or len(witness) != 2
                or not all(isinstance(index, int) and not isinstance(index, bool) and 0 <= index < factor_points for index in witness)
                or witness[0] > witness[1]
            ):
                raise Stage21Error(f"row {label} witness indices are invalid")
        priority = row.get("selection_priority_blake3")
        if (
            row.get("candidate_kind") != expected_kind
            or row.get("candidate_attempts_scope") != expected_scope
            or (is_planted and attempts != counter + 1)
            or (is_planted and construction is None)
            or (not is_planted and construction is not None)
            or row.get("construction_witness_verified") is not (True if is_planted else None)
            or (is_planted and (not isinstance(priority, str) or HEX64.fullmatch(priority) is None))
            or (not is_planted and priority is not None)
            or (hit and verified is None)
            or (not hit and verified is not None)
            or row.get("point_witness_verified") is not (True if hit else None)
        ):
            raise Stage21Error("row construction, selection, or witness covariates are inconsistent")
    if any(row_hits[name] != arms[name]["hits"] for name in expected_mix):
        raise Stage21Error("row-level pair-table outcomes differ from the arm hit counts")
    timing = result.get("timing_ns", {})
    required_timers = {
        "curve_and_subgroup_construction", "factor_base_predicate_and_materialization",
        "factor_base_projected_column_census", "natural_target_generation",
        "cofactor_class_construction", "canonical_pair_table_and_planted_priority_selection",
        "planted_construction_witness_and_subgroup_verification",
        "proven_unsat_target_screening", "covariate_extraction", "end_to_end",
    }
    require_exact_object(
        timing, required_timers | {"clock", "overlap_note"}, "yield producer timing record"
    )
    if (
        timing.get("clock") != "std::time::Instant monotonic elapsed wall time"
        or timing.get("overlap_note")
        != "planted target selection is performed inside the canonical pair-table scan and is not an additive stage"
        or any(
        require_nonnegative_int(timing[name], f"timing_ns.{name}") > timing["end_to_end"]
        for name in required_timers - {"end_to_end"}
        )
    ):
        raise Stage21Error("internal timing records are incomplete or exceed the enclosing timer")
    if require_nonnegative_int(timing["end_to_end"], "timing_ns.end_to_end") <= 0:
        raise Stage21Error("the end-to-end timer must be positive")
    if (
        timing["natural_target_generation"]
        != arms["natural"]["timing_ns"]["target_generation_and_selection"]
        or timing["canonical_pair_table_and_planted_priority_selection"]
        != arms["planted_sat"]["timing_ns"]["target_generation_and_selection"]
        or timing["proven_unsat_target_screening"]
        != arms["proven_unsat"]["timing_ns"]["target_generation_and_selection"]
    ):
        raise Stage21Error("top-level and arm target-generation timers disagree")
    high_ops = require_exact_object(
        result.get("high_level_operation_counts"),
        {
            "pair_table_group_additions", "factor_base_projection_scalar_multiplications",
            "cofactor_class_scalar_multiplications", "cofactor_class_negations",
            "planted_pair_priority_hashes", "planted_final_subgroup_scalar_multiplications",
            "planted_construction_witness_verification_additions",
            "target_frobenius_applications", "natural", "planted_sat", "proven_unsat",
            "scope_note",
        },
        "yield producer high-level operation counts",
    )
    for name in (
        "pair_table_group_additions", "factor_base_projection_scalar_multiplications",
        "cofactor_class_scalar_multiplications", "cofactor_class_negations",
        "planted_pair_priority_hashes", "planted_final_subgroup_scalar_multiplications",
        "planted_construction_witness_verification_additions", "target_frobenius_applications",
    ):
        require_nonnegative_int(high_ops.get(name), f"high_level_operation_counts.{name}")
    expected_arm_ops = {
        "natural": arms["natural"]["high_level_operations"],
        "planted_sat": {
            "canonical_pair_candidates": canonical_pairs,
            "final_pair_table_lookups": arms["planted_sat"]["count"],
            "witness_verification_group_additions": arms["planted_sat"]["hits"],
        },
        "proven_unsat": arms["proven_unsat"]["high_level_operations"],
    }
    expected_arm_ops["natural"] = {
        "generation": expected_arm_ops["natural"]["generation"],
        "final_pair_table_lookups": expected_arm_ops["natural"]["final_pair_table_lookups"],
        "witness_verification_group_additions": expected_arm_ops["natural"]["point_witness_verification_additions"],
    }
    expected_arm_ops["proven_unsat"] = {
        "generation": expected_arm_ops["proven_unsat"]["generation"],
        "final_pair_table_lookups": expected_arm_ops["proven_unsat"]["final_pair_table_lookups"],
        "witness_verification_group_additions": expected_arm_ops["proven_unsat"]["point_witness_verification_additions"],
    }
    if (
        high_ops["pair_table_group_additions"] != canonical_pairs
        or high_ops["factor_base_projection_scalar_multiplications"] != factor_points
        or high_ops["cofactor_class_scalar_multiplications"] != factor_points
        or high_ops["cofactor_class_negations"] != factor_points
        or high_ops["planted_final_subgroup_scalar_multiplications"] != expected_mix["planted_sat"]
        or high_ops["planted_construction_witness_verification_additions"] != expected_mix["planted_sat"]
        or high_ops["target_frobenius_applications"] != orbit_applications
        or high_ops["natural"] != expected_arm_ops["natural"]
        or high_ops["planted_sat"] != expected_arm_ops["planted_sat"]
        or high_ops["proven_unsat"] != expected_arm_ops["proven_unsat"]
    ):
        raise Stage21Error("high-level operation counts disagree with the payload")
    sizes = result.get("retained_size_lower_bounds_bytes", {})
    require_exact_object(
        sizes,
        {
            "scope", "field_element_bytes", "encoded_point_bytes", "factor_base_payload",
            "canonical_factor_point_clone_payload", "pair_table_live_key_and_witness_payload",
            "pair_table_capacity_key_and_witness_payload", "selected_target_coordinate_payload",
            "sum_using_pair_table_capacity_payload",
        },
        "yield producer retained-size record",
    )
    for name in set(sizes) - {"scope"}:
        require_nonnegative_int(sizes.get(name), f"retained_size_lower_bounds_bytes.{name}")
    field_bytes = (expected_n + 7) // 8
    encoded_point_bytes = 1 + 2 * field_bytes
    expected_factor_payload = predicate["abscissae"] * field_bytes + predicate["dimension"] * field_bytes + factor_points * encoded_point_bytes
    expected_clone_payload = factor_points * encoded_point_bytes
    expected_live_payload = unique_targets * 16
    expected_capacity_payload = table_capacity * 16
    expected_target_payload = sum(expected_mix.values()) * encoded_point_bytes
    if (
        sizes["field_element_bytes"] != field_bytes
        or sizes["encoded_point_bytes"] != encoded_point_bytes
        or sizes["factor_base_payload"] != expected_factor_payload
        or sizes["canonical_factor_point_clone_payload"] != expected_clone_payload
        or sizes["pair_table_live_key_and_witness_payload"] != expected_live_payload
        or sizes["pair_table_capacity_key_and_witness_payload"] != expected_capacity_payload
        or sizes["selected_target_coordinate_payload"] != expected_target_payload
        or sizes["sum_using_pair_table_capacity_payload"]
        != expected_factor_payload + expected_clone_payload + expected_capacity_payload + expected_target_payload
    ):
        raise Stage21Error("retained-size lower bounds disagree with the payload dimensions")
    hashes = result.get("hashes", {})
    require_exact_object(
        hashes,
        {
            "policy_blake3", "predicate_blake3", "factor_base_blake3",
            "canonical_pair_transcript_blake3", "ordered_covariate_rows_blake3",
            "arms", "result_binding_blake3",
        },
        "yield producer hashes",
    )
    for name in (
        "policy_blake3",
        "predicate_blake3",
        "factor_base_blake3",
        "canonical_pair_transcript_blake3",
        "ordered_covariate_rows_blake3",
        "result_binding_blake3",
    ):
        require_hex(hashes.get(name), f"hashes.{name}")
    arm_hashes = require_exact_object(
        hashes.get("arms"),
        {
            "natural_targets", "natural_witness_results", "planted_targets",
            "planted_witness_results", "proven_unsat_targets",
            "proven_unsat_witness_results",
        },
        "yield producer arm hashes",
    )
    for name, value in arm_hashes.items():
        require_hex(value, f"hashes.arms.{name}")
    coverage = require_exact_object(
        measurement.get("eta_or_coverage_policy"),
        {"m", "target_distribution", "oracle", "policy_blake3"},
        "yield producer coverage policy",
    )
    if (
        measurement.get("base_id_or_hash") != hashes["factor_base_blake3"]
        or coverage.get("m") != 2
        or coverage.get("policy_blake3")
        != hashes["policy_blake3"]
        or factor_base.get("predicate_blake3") != hashes["predicate_blake3"]
        or factor_base.get("factor_base_blake3") != hashes["factor_base_blake3"]
        or result.get("ordered_covariate_rows_blake3") != hashes["ordered_covariate_rows_blake3"]
        or pair_table.get("canonical_pair_transcript_blake3")
        != hashes["canonical_pair_transcript_blake3"]
    ):
        raise Stage21Error("yield producer hash bindings disagree")
    recomputed = {
        "policy_blake3": json_blake3(policy),
        "predicate_blake3": json_blake3(predicate),
        "ordered_covariate_rows_blake3": json_blake3(rows),
    }
    rows_by_arm = {
        name: [row for row in rows if row["arm"] == name]
        for name in expected_mix
    }
    recomputed_arm_hashes: dict[str, str] = {}
    arm_field_names = {
        "natural": ("natural_targets", "natural_witness_results"),
        "planted_sat": ("planted_targets", "planted_witness_results"),
        "proven_unsat": ("proven_unsat_targets", "proven_unsat_witness_results"),
    }
    for name, arm_rows in rows_by_arm.items():
        target_bytes = bytearray(ARM_TARGET_HASH_DOMAINS[name])
        witness_bytes = bytearray(WITNESS_HASH_DOMAIN)
        for row in arm_rows:
            packed_bytes = row["packed_target"].to_bytes(8, "little")
            target_bytes.extend(packed_bytes)
            witness_bytes.extend(packed_bytes)
            witness = row["verified_witness_indices"]
            if witness is None:
                witness_bytes.extend((2**64 - 1).to_bytes(8, "little"))
            else:
                witness_bytes.extend(witness[0].to_bytes(4, "little"))
                witness_bytes.extend(witness[1].to_bytes(4, "little"))
        target_field, witness_field = arm_field_names[name]
        recomputed_arm_hashes[target_field] = blake3_hex(bytes(target_bytes))
        recomputed_arm_hashes[witness_field] = blake3_hex(bytes(witness_bytes))
        if (
            arms[name]["targets_blake3"] != recomputed_arm_hashes[target_field]
            or arms[name]["witness_results_blake3"] != recomputed_arm_hashes[witness_field]
        ):
            raise Stage21Error(f"{name} target or witness transcript hash differs from its rows")
    result_binding = {
        "policy_blake3": hashes["policy_blake3"],
        "predicate_blake3": hashes["predicate_blake3"],
        "factor_base_blake3": hashes["factor_base_blake3"],
        "canonical_pair_transcript_blake3": hashes["canonical_pair_transcript_blake3"],
        "ordered_covariate_rows_blake3": hashes["ordered_covariate_rows_blake3"],
        "arm_hashes": arm_hashes,
        "natural_hits": arms["natural"]["hits"],
        "planted_hits": arms["planted_sat"]["hits"],
        "proven_unsat_hits": arms["proven_unsat"]["hits"],
    }
    recomputed["result_binding_blake3"] = json_blake3(result_binding)
    if (
        any(hashes[name] != value for name, value in recomputed.items())
        or arm_hashes != recomputed_arm_hashes
    ):
        raise Stage21Error("derivable BLAKE3 commitments differ from the retained payload")
    resources = result.get("resource_accounting_boundary", {})
    require_exact_object(
        resources,
        {
            "single_core_elapsed_seconds", "user_cpu_seconds", "system_cpu_seconds",
            "total_core_seconds", "peak_rss_bytes", "host_identity", "executable_hash",
            "source_revision", "external_process_receipt_required", "reason",
        },
        "yield producer resource boundary",
    )
    for field in (
        "single_core_elapsed_seconds",
        "user_cpu_seconds",
        "system_cpu_seconds",
        "total_core_seconds",
        "peak_rss_bytes",
        "host_identity",
        "executable_hash",
        "source_revision",
    ):
        if resources.get(field) is not None:
            raise Stage21Error(f"portable producer made an unsupported process claim: {field}")
    if resources.get("external_process_receipt_required") is not True:
        raise Stage21Error("yield producer no longer requires its external process receipt")
    non_claims = result.get("non_claims")
    if not isinstance(non_claims, list) or not non_claims or not all(isinstance(item, str) for item in non_claims):
        raise Stage21Error("yield producer non-claims use an unexpected schema")
    if not isinstance(result.get("claim_boundary"), str) or "SOTA" not in " ".join(non_claims):
        raise Stage21Error("yield producer omitted its narrow claim boundary")
    return {
        "schema": "koblitz_relation_yield_payload_replay.v1",
        "status": "derivable_payload_relationships_and_hashes_replayed",
        "recomputed_blake3": {
            **recomputed,
            "arms": recomputed_arm_hashes,
        },
        "factor_base_materialization_replayed": False,
        "canonical_pair_oracle_replayed": False,
        "pending_independent_payload_replay": [
            "factor-base point materialization and ordered factor-base hash",
            "all canonical pair additions and canonical-pair transcript hash",
            "exact curve re-addition of retained witness indices",
        ],
        "scientific_measurement_admitted": False,
    }


def inventory(root: Path, excluded: set[str]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root).as_posix()
        metadata = path.lstat()
        if stat.S_ISLNK(metadata.st_mode):
            raise Stage21Error(f"artifact tree contains symlink {relative}")
        if stat.S_ISDIR(metadata.st_mode):
            continue
        if not stat.S_ISREG(metadata.st_mode):
            raise Stage21Error(f"artifact tree contains non-regular file {relative}")
        if relative in excluded:
            continue
        data = path.read_bytes()
        records.append({"path": relative, "bytes": len(data), "sha256": sha256_bytes(data)})
    return records


def aggregate_resources(processes: list[dict[str, Any]]) -> dict[str, Any]:
    metrics = [process["process"]["metrics"] for process in processes]
    return {
        "processes": len(metrics),
        "summed_process_wall_seconds": sum(item["wall_seconds"] for item in metrics),
        "summed_user_seconds": sum(item["user_seconds"] for item in metrics),
        "summed_system_seconds": sum(item["system_seconds"] for item in metrics),
        "total_core_seconds": sum(item["total_core_seconds"] for item in metrics),
        "single_core_elapsed_seconds": None,
        "single_core_seconds_legacy_alias": "total_core_seconds",
        "peak_process_rss_bytes": max(item["peak_rss_bytes"] for item in metrics),
        "aggregate_parallel_rss_bytes": None,
        "parallel_rss_reason": "scientific producers are sequential and request one worker; the fresh Cargo build may overlap its parent with a compiler child, so simultaneous aggregate build-tree RSS is unavailable",
        "conflicts": None,
        "conflicts_reason": "the exact pair-table reference oracle is not a conflict-driven SAT solver",
    }


def expected_driver_command(
    protocol: Path,
    output: Path,
    meter: Path,
    smoke: bool,
    allow_dirty: bool,
) -> list[str]:
    command = [
        str(Path(sys.executable).resolve()),
        str(Path(__file__).resolve()),
        "run",
        "--protocol",
        str(protocol.resolve()),
        "--output",
        str(output.resolve()),
        "--meter",
        str(meter.resolve()),
    ]
    if smoke:
        command.append("--smoke")
    if allow_dirty:
        command.append("--allow-dirty")
    return command


def plan(args: argparse.Namespace) -> dict[str, Any]:
    production = not args.smoke
    protocol, protocol_bytes, protocol_path = load_protocol(args.protocol, production)
    if production and args.allow_dirty:
        raise Stage21Error("--allow-dirty is available only for operational smoke")
    output = safe_new_directory(args.output, "planned run output")
    meter = args.meter.resolve(strict=True)
    file_identity(meter, "process meter")
    if meter != DEFAULT_METER.resolve(strict=True):
        raise Stage21Error("Stage-21 requires the committed process meter")
    outer_paths = {
        "stdout": safe_new_file(args.outer_stdout, "outer stdout"),
        "stderr": safe_new_file(args.outer_stderr, "outer stderr"),
        "metrics": safe_new_file(args.outer_metrics, "outer metrics"),
    }
    if len(set(outer_paths.values())) != len(outer_paths):
        raise Stage21Error("outer stdout, stderr, and metrics paths must be distinct")
    if any(path == output or output in path.parents for path in outer_paths.values()):
        raise Stage21Error("outer artifacts must be outside the write-once run output")
    execution = protocol["execution"]
    command = expected_driver_command(
        protocol_path, output, meter, args.smoke, args.allow_dirty
    )
    return {
        "schema": "koblitz_relation_yield_execution_plan.v1",
        "mode": "production" if production else "smoke",
        "scientific_evidence": False,
        "scientific_measurement_eligible_after_verified_execution": production,
        "protocol": {"path": str(protocol_path), "bytes": len(protocol_bytes), "sha256": sha256_bytes(protocol_bytes)},
        "output": str(output),
        "driver_command": command,
        "required_outer_meter_command": meter_command(
            meter=meter,
            cwd=REPO,
            timeout=execution["whole_driver_watchdog_seconds"],
            stdout=outer_paths["stdout"],
            stderr=outer_paths["stderr"],
            metrics=outer_paths["metrics"],
            command=command,
        ),
        "claim_boundary": protocol["claim_boundary"],
    }


def execute(args: argparse.Namespace) -> dict[str, Any]:
    production = not args.smoke
    protocol, protocol_bytes, protocol_path = load_protocol(args.protocol, production)
    output = safe_new_directory(args.output, "run output")
    meter = args.meter.resolve(strict=True)
    file_identity(meter, "process meter")
    if meter != DEFAULT_METER.resolve(strict=True):
        raise Stage21Error("Stage-21 requires the committed process meter")
    source_state = git_state()
    if source_state["dirty"] and not args.allow_dirty:
        raise Stage21Error("execution requires a clean committed checkout unless smoke uses --allow-dirty")
    if production and args.allow_dirty:
        raise Stage21Error("--allow-dirty is available only for operational smoke")
    lock = install_frozen_lock()
    cargo = shutil.which("cargo")
    rustc = shutil.which("rustc")
    if cargo is None or rustc is None:
        raise Stage21Error("the Rust build toolchain is unavailable")
    toolchain = {
        "cargo": tool_identity(Path(cargo), "cargo"),
        "rustc": tool_identity(Path(rustc), "rustc"),
    }
    cargo_home = Path(os.environ.get("CARGO_HOME", Path.home() / ".cargo")).resolve()
    build_environment = safe_child_environment(
        cargo_home=cargo_home,
        rustc=Path(toolchain["rustc"]["path"]),
    )
    scientific_environment = safe_child_environment()
    direct_sources = {
        relative: file_identity(REPO / relative, f"source {relative}")
        for relative in DIRECT_SOURCE_PATHS
    }
    output.mkdir(parents=True)
    (output / "inputs").mkdir()
    (output / "tasks").mkdir()
    (output / "binaries").mkdir()
    write_new(output / "inputs" / "protocol.json", protocol_bytes)
    write_json_new(output / "inputs" / "host.json", host_identity(toolchain))
    write_json_new(output / "inputs" / "source.json", {
        "schema": "koblitz_relation_yield_source_binding.v1",
        "git": source_state,
        "direct_sources": direct_sources,
        "dependency_lock": lock,
        "toolchain": toolchain,
    })
    build_examples = [YIELD_BINARY_NAME] if args.smoke else [DISCOVERY_BINARY_NAME, YIELD_BINARY_NAME]
    build_target = output / "build-target"
    build_command = [
        toolchain["cargo"]["path"],
        "build",
        "--release",
        "--locked",
        "--offline",
        "-j",
        "1",
        "--target-dir",
        str(build_target),
    ]
    for example in build_examples:
        build_command.extend(["--example", example])
    build_receipt, _ = run_metered(
        output,
        "00-build",
        build_command,
        protocol["execution"]["build_watchdog_seconds"],
        list(direct_sources.values()) + [lock["frozen"], lock["workspace"]],
        meter,
        build_environment,
        toolchain["cargo"],
    )
    binaries: dict[str, dict[str, Any]] = {}
    for name in build_examples:
        source = build_target / "release" / "examples" / name
        built_identity = tool_identity(source, f"built {name}")
        data = Path(built_identity["path"]).read_bytes()
        destination = output / "binaries" / name
        write_new(destination, data)
        destination.chmod(0o755)
        binaries[name] = file_identity(destination, f"frozen {name}")
    processes = [build_receipt]
    discoveries: list[dict[str, Any]] = []
    if production:
        for curve_a in (0, 1):
            command = [str((output / "binaries" / DISCOVERY_BINARY_NAME).resolve()), "23", str(curve_a), "2", "12"]
            receipt, stdout = run_metered(
                output,
                f"0{curve_a + 1}-discovery-a{curve_a}",
                command,
                protocol["execution"]["discovery_watchdog_seconds"],
                [binaries[DISCOVERY_BINARY_NAME], direct_sources["examples/koblitz_public_factor_base_discovery.rs"]],
                meter,
                scientific_environment,
                binaries[DISCOVERY_BINARY_NAME],
            )
            result = parse_stdout_json(stdout, f"curve-a={curve_a} discovery stdout")
            validate_discovery(result, curve_a)
            write_json_new(stdout.parent / "result.json", result)
            discoveries.append(result)
            processes.append(receipt)
        winner_a, winner = cross_curve_winner(discoveries)
        expected_base = protocol["factor_base"]
        if (
            winner_a != 0
            or winner.get("dimension") != expected_base["dimension"]
            or winner.get("divisor_indices") != expected_base["divisor_indices"]
            or winner.get("divisor_polynomial") != expected_base["divisor_polynomial"]
            or winner.get("linearised_exponents") != expected_base["linearised_exponents"]
            or winner.get("abscissae") != expected_base["expected_abscissae"]
            or winner.get("rational_points") != expected_base["expected_rational_points"]
            or winner.get("signed_frobenius_orbits_before_projection")
            != expected_base["expected_signed_frobenius_orbits_before_projection"]
            or winner.get("projected_signed_frobenius_orbits")
            != expected_base["expected_projected_signed_frobenius_columns"]
        ):
            raise Stage21Error("fresh public discovery did not reproduce the frozen cross-curve selection")
        discovery_decision = {"curve_a": winner_a, "candidate": winner}
    else:
        discovery_decision = None
    yield_command = [str((output / "binaries" / YIELD_BINARY_NAME).resolve())]
    if args.smoke:
        yield_command.extend(protocol["producer"]["smoke_arguments"])
    yield_receipt, yield_stdout = run_metered(
        output,
        "03-yield" if production else "01-yield-smoke",
        yield_command,
        protocol["execution"]["producer_watchdog_seconds"],
        [binaries[YIELD_BINARY_NAME], direct_sources["examples/koblitz_relation_yield_bridge.rs"]],
        meter,
        scientific_environment,
        binaries[YIELD_BINARY_NAME],
    )
    yield_result = parse_stdout_json(yield_stdout, "yield producer stdout")
    validate_yield_result(yield_result, production, protocol)
    write_json_new(yield_stdout.parent / "result.json", yield_result)
    processes.append(yield_receipt)
    for name, identity in toolchain.items():
        validate_tool_identity(identity, name)
    final_state = git_state()
    if final_state != source_state:
        raise Stage21Error("implementation Git state changed during execution")
    resources = aggregate_resources(processes)
    summary = {
        "schema": "koblitz_relation_yield_run_summary.v1",
        "status": "complete",
        "evidence_class": (
            "finite_public_synthetic_candidate_pending_independent_payload_replay"
            if production else "operational_smoke_only"
        ),
        "production": production,
        "protocol_sha256": sha256_bytes(protocol_bytes),
        "source_revision": source_state,
        "binaries": binaries,
        "public_discovery_decision": discovery_decision,
        "yield_result": yield_result,
        "resources": resources,
        "outer_driver_accounting": "required_before_admission",
        "full_cost_gate_passed": False,
        "full_cost_blockers": [
            "measured single-core elapsed time and simultaneous aggregate Cargo build-tree RSS are unavailable",
            "prior dependency source acquisition in the Cargo cache is outside this bridge receipt",
            "this bridge does not collect a relation matrix or perform modular linear algebra",
            "this bridge does not execute a scalar-hidden end-to-end index-calculus run",
            "this bridge has no matched automorphism-optimized Pollard-rho arm",
            "independent external reproduction and novelty review remain absent",
        ],
        "claim_boundary": protocol["claim_boundary"],
    }
    write_json_new(output / "run-summary.json", summary)
    frozen_inventory = inventory(output, {"run-seal.json"})
    seal_payload = {
        "schema": RUN_SEAL_SCHEMA,
        "status": "relation_yield_outputs_frozen",
        "created_at": now(),
        "production": production,
        "protocol_sha256": sha256_bytes(protocol_bytes),
        "source_commit": source_state["commit"],
        "yield_result_binding_blake3": yield_result.get("hashes", {}).get("result_binding_blake3"),
        "inventory": frozen_inventory,
        "inventory_sha256": canonical_sha256(frozen_inventory),
        "outer_expected_command": expected_driver_command(
            protocol_path, output, meter, args.smoke, args.allow_dirty
        ),
        "outer_driver_accounting": "required_before_admission",
    }
    seal = dict(seal_payload)
    seal["seal_payload_sha256"] = canonical_sha256(seal_payload)
    write_json_new(output / "run-seal.json", seal)
    return seal


def validate_source_binding(
    run_root: Path, summary: dict[str, Any], production: bool
) -> dict[str, Any]:
    source, _ = read_json(run_root / "inputs" / "source.json", "source binding")
    if set(source) != {"schema", "git", "direct_sources", "dependency_lock", "toolchain"} or source.get("schema") != "koblitz_relation_yield_source_binding.v1":
        raise Stage21Error("source binding uses an unexpected schema")
    state = source.get("git")
    if not isinstance(state, dict) or set(state) != {
        "commit", "branch", "dirty", "porcelain", "committed_tree_sha256"
    }:
        raise Stage21Error("source binding has invalid Git provenance")
    require_hex(state.get("commit"), "source commit", HEX40)
    require_hex(state.get("committed_tree_sha256"), "committed tree hash")
    if (
        not isinstance(state.get("dirty"), bool)
        or not isinstance(state.get("porcelain"), list)
        or not all(isinstance(line, str) for line in state["porcelain"])
        or state["dirty"] != bool(state["porcelain"])
        or (production and state["dirty"])
        or summary.get("source_revision") != state
    ):
        raise Stage21Error("source binding cleanliness or summary provenance changed")
    tree = command_output(["git", "ls-tree", "-r", "--full-tree", state["commit"]])
    if tree["returncode"] != 0 or sha256_bytes((tree["stdout"] + "\n").encode()) != state["committed_tree_sha256"]:
        raise Stage21Error("committed implementation tree no longer matches its binding")
    direct = source.get("direct_sources")
    if not isinstance(direct, dict) or set(direct) != set(DIRECT_SOURCE_PATHS):
        raise Stage21Error("source binding does not cover the exact direct-source set")
    for relative in DIRECT_SOURCE_PATHS:
        identity = validate_file_identity(direct[relative], f"source {relative}")
        if identity["path"] != str((REPO / relative).resolve()):
            raise Stage21Error(f"source {relative} path differs from the implementation checkout")
    lock = source.get("dependency_lock")
    if not isinstance(lock, dict) or set(lock) != {"frozen", "workspace", "installed_by_runner"}:
        raise Stage21Error("dependency-lock binding uses an unexpected schema")
    frozen = validate_file_identity(lock["frozen"], "frozen dependency lock")
    workspace = validate_file_identity(lock["workspace"], "workspace dependency lock")
    if (
        frozen["path"] != str(FROZEN_LOCK.resolve())
        or workspace["path"] != str(WORKSPACE_LOCK.resolve())
        or any(frozen[field] != workspace[field] for field in ("bytes", "sha256"))
        or not isinstance(lock["installed_by_runner"], bool)
    ):
        raise Stage21Error("dependency-lock identities disagree")
    toolchain = source.get("toolchain")
    if not isinstance(toolchain, dict) or set(toolchain) != {"cargo", "rustc"}:
        raise Stage21Error("source binding has an unexpected toolchain")
    for name in ("cargo", "rustc"):
        validate_tool_identity(toolchain[name], name)
    return source


def expected_process_contracts(
    run_root: Path,
    protocol: dict[str, Any],
    production: bool,
    source: dict[str, Any],
    binaries: dict[str, Any],
) -> list[dict[str, Any]]:
    build_examples = (
        [DISCOVERY_BINARY_NAME, YIELD_BINARY_NAME]
        if production
        else [YIELD_BINARY_NAME]
    )
    build_command = [
        source["toolchain"]["cargo"]["path"],
        "build", "--release", "--locked", "--offline", "-j", "1",
        "--target-dir", str(run_root / "build-target"),
    ]
    for example in build_examples:
        build_command.extend(["--example", example])
    build_inputs = [source["direct_sources"][name] for name in DIRECT_SOURCE_PATHS]
    build_inputs.extend([
        source["dependency_lock"]["frozen"],
        source["dependency_lock"]["workspace"],
    ])
    contracts = [{
        "name": "00-build",
        "command": build_command,
        "timeout": protocol["execution"]["build_watchdog_seconds"],
        "inputs": build_inputs,
        "executable": source["toolchain"]["cargo"],
        "environment_kind": "build",
    }]
    if production:
        for curve_a in (0, 1):
            contracts.append({
                "name": f"0{curve_a + 1}-discovery-a{curve_a}",
                "command": [
                    binaries[DISCOVERY_BINARY_NAME]["path"],
                    "23", str(curve_a), "2", "12",
                ],
                "timeout": protocol["execution"]["discovery_watchdog_seconds"],
                "inputs": [
                    binaries[DISCOVERY_BINARY_NAME],
                    source["direct_sources"]["examples/koblitz_public_factor_base_discovery.rs"],
                ],
                "executable": binaries[DISCOVERY_BINARY_NAME],
                "environment_kind": "scientific",
            })
    yield_name = "03-yield" if production else "01-yield-smoke"
    yield_command = [binaries[YIELD_BINARY_NAME]["path"]]
    if not production:
        yield_command.extend(protocol["producer"]["smoke_arguments"])
    contracts.append({
        "name": yield_name,
        "command": yield_command,
        "timeout": protocol["execution"]["producer_watchdog_seconds"],
        "inputs": [
            binaries[YIELD_BINARY_NAME],
            source["direct_sources"]["examples/koblitz_relation_yield_bridge.rs"],
        ],
        "executable": binaries[YIELD_BINARY_NAME],
        "environment_kind": "scientific",
    })
    return contracts


def validate_process_tree(
    run_root: Path,
    protocol: dict[str, Any],
    production: bool,
    source: dict[str, Any],
    summary: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    binaries = summary.get("binaries")
    expected_binary_names = (
        {DISCOVERY_BINARY_NAME, YIELD_BINARY_NAME}
        if production else {YIELD_BINARY_NAME}
    )
    if not isinstance(binaries, dict) or set(binaries) != expected_binary_names:
        raise Stage21Error("run summary has an unexpected binary inventory")
    for name, identity in binaries.items():
        validate_file_identity(identity, f"frozen binary {name}")
        if identity["path"] != str((run_root / "binaries" / name).resolve()):
            raise Stage21Error(f"frozen binary {name} path differs from the run archive")
    contracts = expected_process_contracts(run_root, protocol, production, source, binaries)
    task_root = run_root / "tasks"
    actual_names = tuple(sorted(path.name for path in task_root.iterdir()))
    if actual_names != TASK_NAMES[production] or tuple(row["name"] for row in contracts) != TASK_NAMES[production]:
        raise Stage21Error("run task inventory or order differs from the frozen execution")
    receipts: list[dict[str, Any]] = []
    discovery_results: list[dict[str, Any]] = []
    yield_result: dict[str, Any] | None = None
    payload_replay: dict[str, Any] | None = None
    meter_identity = source["direct_sources"]["scripts/process_meter.py"]
    for contract in contracts:
        name = contract["name"]
        task = task_root / name
        receipt, _ = read_json(task / "receipt.json", f"{name} receipt")
        if set(receipt) != {
            "schema", "name", "process", "intent", "stdout", "stderr", "metrics",
            "executable", "inputs", "environment",
        } or receipt.get("schema") != "koblitz_relation_yield_process_receipt.v1" or receipt.get("name") != name:
            raise Stage21Error(f"{name} receipt uses an unexpected schema")
        for field in ("intent", "stdout", "stderr", "metrics"):
            validate_file_identity(receipt[field], f"{name} {field}")
            if receipt[field]["path"] != str((task / ("metrics.json" if field == "metrics" else field + (".json" if field == "intent" else ""))).resolve()):
                raise Stage21Error(f"{name} {field} identity uses the wrong artifact path")
        process = process_metrics(task / "metrics.json", f"{name} metrics")
        if receipt["process"] != process:
            raise Stage21Error(f"{name} receipt differs from its raw process metrics")
        if (
            process["command"] != contract["command"]
            or process["watchdog_seconds"] != contract["timeout"]
            or process["returncode"] != 0
            or process["timed_out"]
            or process["orphan_group_terminated"]
            or receipt["executable"] != contract["executable"]
            or receipt["inputs"] != contract["inputs"]
        ):
            raise Stage21Error(f"{name} differs from its frozen process contract")
        expected_environment = safe_child_environment()
        if contract["environment_kind"] == "build":
            environment = receipt.get("environment", {})
            cargo_home = environment.get("CARGO_HOME")
            if not isinstance(cargo_home, str) or not Path(cargo_home).is_absolute():
                raise Stage21Error("build receipt lacks an absolute Cargo home")
            expected_environment = safe_child_environment(
                cargo_home=Path(cargo_home),
                rustc=Path(source["toolchain"]["rustc"]["path"]),
            )
        if receipt.get("environment") != expected_environment:
            raise Stage21Error(f"{name} child environment differs from the sanitized contract")
        intent, _ = read_json(task / "intent.json", f"{name} intent")
        expected_intent = {
            "schema": "koblitz_relation_yield_process_intent.v1",
            "name": name,
            "command": contract["command"],
            "cwd": str(REPO.resolve()),
            "watchdog_seconds": contract["timeout"],
            "environment": expected_environment,
            "inputs": contract["inputs"],
            "meter": meter_identity,
            "meter_exclusive_create": True,
            "executable": contract["executable"],
            "stdin": "devnull",
            "close_fds": True,
        }
        if intent != expected_intent:
            raise Stage21Error(f"{name} raw intent differs from the frozen process contract")
        expected_files = {"intent.json", "stdout", "stderr", "metrics.json", "receipt.json"}
        if name != "00-build":
            expected_files.add("result.json")
            stdout_result = parse_stdout_json(task / "stdout", f"{name} stdout")
            result, _ = read_json(task / "result.json", f"{name} result")
            if result != stdout_result:
                raise Stage21Error(f"{name} result differs from its raw stdout")
            if "discovery" in name:
                curve_a = 0 if name.endswith("a0") else 1
                validate_discovery(result, curve_a)
                discovery_results.append(result)
            else:
                payload_replay = validate_yield_result(result, production, protocol)
                yield_result = result
        actual_files = {path.name for path in task.iterdir()}
        if actual_files != expected_files:
            raise Stage21Error(f"{name} contains an unexpected artifact set")
        receipts.append(receipt)
    if yield_result is None or payload_replay is None:
        raise Stage21Error("run task inventory lacks the yield result")
    recomputed_resources = aggregate_resources(receipts)
    if summary.get("resources") != recomputed_resources:
        raise Stage21Error("run summary resources differ from the raw process receipts")
    return receipts, discovery_results, yield_result, payload_replay


def validate_run_seal(
    run_root: Path,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    seal, _ = read_json(run_root / "run-seal.json", "run seal")
    if seal.get("schema") != RUN_SEAL_SCHEMA or seal.get("status") != "relation_yield_outputs_frozen":
        raise Stage21Error("run seal is not terminal under the expected schema")
    payload = {key: value for key, value in seal.items() if key != "seal_payload_sha256"}
    if seal.get("seal_payload_sha256") != canonical_sha256(payload):
        raise Stage21Error("run-seal self-hash is invalid")
    actual = inventory(run_root, {"run-seal.json"})
    if actual != seal.get("inventory") or canonical_sha256(actual) != seal.get("inventory_sha256"):
        raise Stage21Error("run artifact inventory changed after sealing")
    summary, _ = read_json(run_root / "run-summary.json", "run summary")
    if summary.get("schema") != "koblitz_relation_yield_run_summary.v1" or summary.get("status") != "complete":
        raise Stage21Error("run summary is not complete")
    production = seal.get("production")
    if not isinstance(production, bool) or summary.get("production") is not production:
        raise Stage21Error("run summary mode differs from its seal")
    protocol, protocol_bytes = read_json(run_root / "inputs" / "protocol.json", "archived protocol")
    validate_protocol(protocol, production)
    protocol_hash = sha256_bytes(protocol_bytes)
    if seal.get("protocol_sha256") != protocol_hash or summary.get("protocol_sha256") != protocol_hash:
        raise Stage21Error("archived protocol differs from the run bindings")
    source = validate_source_binding(run_root, summary, production)
    _, discoveries, yield_result, payload_replay = validate_process_tree(
        run_root, protocol, production, source, summary
    )
    if summary.get("yield_result") != yield_result:
        raise Stage21Error("run summary yield result differs from the raw producer output")
    binding = yield_result.get("hashes", {}).get("result_binding_blake3")
    require_hex(binding, "yield result binding")
    if seal.get("yield_result_binding_blake3") != binding or seal.get("source_commit") != source["git"]["commit"]:
        raise Stage21Error("run seal result or source binding changed")
    outer_command = seal.get("outer_expected_command")
    if not isinstance(outer_command, list) or outer_command.count("--protocol") != 1:
        raise Stage21Error("run seal lacks one protocol-bound outer command")
    protocol_index = outer_command.index("--protocol")
    if protocol_index + 1 >= len(outer_command):
        raise Stage21Error("run seal outer command lacks its protocol path")
    recorded_protocol_path = Path(outer_command[protocol_index + 1])
    if production and recorded_protocol_path.resolve() != DEFAULT_PROTOCOL.resolve():
        raise Stage21Error("production outer command did not use the committed default protocol")
    if production and "--allow-dirty" in outer_command:
        raise Stage21Error("production outer command cannot allow a dirty checkout")
    expected_outer_command = expected_driver_command(
        recorded_protocol_path,
        run_root,
        DEFAULT_METER,
        not production,
        "--allow-dirty" in outer_command,
    )
    if outer_command != expected_outer_command:
        raise Stage21Error("run seal outer command differs from the independently reconstructed command")
    if production:
        winner_a, winner = cross_curve_winner(discoveries)
        expected_decision = {"curve_a": winner_a, "candidate": winner}
        if summary.get("public_discovery_decision") != expected_decision:
            raise Stage21Error("run summary discovery decision is not independently reproducible")
    elif summary.get("public_discovery_decision") is not None:
        raise Stage21Error("smoke run unexpectedly claims a production discovery decision")
    expected_class = (
        "finite_public_synthetic_candidate_pending_independent_payload_replay"
        if production else "operational_smoke_only"
    )
    if (
        summary.get("evidence_class") != expected_class
        or summary.get("outer_driver_accounting") != "required_before_admission"
        or summary.get("full_cost_gate_passed") is not False
        or seal.get("outer_driver_accounting") != "required_before_admission"
    ):
        raise Stage21Error("run evidence class or outer-accounting gate changed")
    return seal, summary, protocol, payload_replay


def verify_outer(
    outer_path: Path,
    seal: dict[str, Any],
    summary: dict[str, Any],
    protocol: dict[str, Any],
) -> dict[str, Any]:
    outer = process_metrics(outer_path, "outer-driver metrics")
    if outer["returncode"] != 0 or outer["timed_out"] or outer["orphan_group_terminated"]:
        raise Stage21Error("outer driver did not terminate cleanly")
    if outer["command"] != seal["outer_expected_command"]:
        raise Stage21Error("outer driver command differs from the frozen run command")
    if outer["watchdog_seconds"] != protocol["execution"]["whole_driver_watchdog_seconds"]:
        raise Stage21Error("outer driver watchdog differs from the frozen protocol")
    charged = summary["resources"]
    metrics = outer["metrics"]
    if metrics["total_core_seconds"] + 1e-9 < charged["total_core_seconds"]:
        raise Stage21Error("outer CPU receipt does not enclose charged child CPU")
    if metrics["wall_seconds"] + 1e-9 < charged["summed_process_wall_seconds"]:
        raise Stage21Error("outer wall receipt does not enclose sequential child wall time")
    if metrics["peak_rss_bytes"] < charged["peak_process_rss_bytes"]:
        raise Stage21Error("outer peak RSS does not enclose child peak RSS")
    return outer


def verify(args: argparse.Namespace) -> dict[str, Any]:
    run_root = args.run_root.resolve(strict=True)
    if run_root.is_symlink() or not run_root.is_dir():
        raise Stage21Error("run root must be a real directory")
    output = safe_new_directory(args.output, "verification output")
    seal, summary, protocol, payload_replay = validate_run_seal(run_root)
    outer = verify_outer(args.outer_metrics.resolve(strict=True), seal, summary, protocol)
    output.mkdir(parents=True)
    verification = {
        "schema": VERIFY_SCHEMA,
        "status": "relation_yield_run_verified",
        "created_at": now(),
        "production": seal["production"],
        "scientific_measurement_admitted": False,
        "measurement_admission_status": (
            "pending_independent_payload_replay"
            if seal["production"] else "operational_smoke_only"
        ),
        "payload_replay": payload_replay,
        "external_portable_verification_satisfied": False,
        "external_portable_verification_blocker": "archived executable and source identities retain original absolute checkout paths; no archive-relative source snapshot replay is implemented",
        "run_root": str(run_root),
        "run_seal": file_identity(run_root / "run-seal.json", "run seal"),
        "run_inventory_sha256": seal["inventory_sha256"],
        "outer_driver_metrics": file_identity(args.outer_metrics.resolve(), "outer metrics"),
        "outer_driver_accounting": outer,
        "measurement_schema": summary["yield_result"]["measurement_schema"],
        "resources": {
            "charged_children": summary["resources"],
            "whole_driver": outer["metrics"],
            "single_core_elapsed_seconds": None,
            "single_core_claim": "unavailable; single_core_seconds in process receipts is a legacy alias of aggregate CPU",
        },
        "full_cost_gate_passed": False,
        "independent_external_reproduction_satisfied": False,
        "claim_boundary": summary["claim_boundary"],
    }
    write_json_new(output / "verification.json", verification)
    verify_inventory = inventory(output, {"verification-seal.json"})
    seal_payload = {
        "schema": VERIFY_SEAL_SCHEMA,
        "status": "verification_frozen",
        "created_at": now(),
        "verification_sha256": file_identity(output / "verification.json")["sha256"],
        "inventory": verify_inventory,
        "inventory_sha256": canonical_sha256(verify_inventory),
    }
    verification_seal = dict(seal_payload)
    verification_seal["seal_payload_sha256"] = canonical_sha256(seal_payload)
    write_json_new(output / "verification-seal.json", verification_seal)
    return verification_seal


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command", required=True)
    for name in ("plan", "run"):
        command = subparsers.add_parser(name)
        command.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
        command.add_argument("--output", type=Path, required=True)
        command.add_argument("--meter", type=Path, default=DEFAULT_METER)
        command.add_argument("--smoke", action="store_true")
        command.add_argument("--allow-dirty", action="store_true")
        if name == "plan":
            command.add_argument("--outer-stdout", type=Path, required=True)
            command.add_argument("--outer-stderr", type=Path, required=True)
            command.add_argument("--outer-metrics", type=Path, required=True)
    verification = subparsers.add_parser("verify")
    verification.add_argument("--run-root", type=Path, required=True)
    verification.add_argument("--outer-metrics", type=Path, required=True)
    verification.add_argument("--output", type=Path, required=True)
    return result


def main() -> None:
    args = parser().parse_args()
    try:
        if args.command == "plan":
            value = plan(args)
        elif args.command == "run":
            value = execute(args)
        else:
            value = verify(args)
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, Stage21Error, subprocess.SubprocessError) as error:
        parser().exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
