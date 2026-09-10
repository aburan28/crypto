#!/usr/bin/env python3
"""Fail-closed verifier for the frozen Stage 18 degree-23 replication panel."""

from __future__ import annotations

import argparse
import copy
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess
import tempfile
import tomllib
from typing import Any


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DEFAULT_PROTOCOL = HERE / "stage-18-degree23-panel-protocol.json"
DEFAULT_AMENDMENT = HERE / "stage-18-amendment-01-lock-correction.json"
DEFAULT_PANEL = HERE / "stage-18-degree23-panel-lock-corrected-20260909"
FAILED_V1_PANEL = HERE / "stage-18-degree23-panel-20260909"

PROTOCOL_SCHEMA = "koblitz_degree23_replication_panel_protocol.v1"
RUN_SCHEMA = "koblitz_degree23_replication_panel_run.v1"
RECEIPT_SCHEMA = "koblitz_degree23_task_receipt.v1"
SUMMARY_SCHEMA = "koblitz_degree23_replication_panel_result.v1"
BUILD_WATCHDOG_SECONDS = 1800.0
ALGORITHM_BASE_COMMIT = "754f76b2e313fff98b3fbed11bc246c3744a1591"
V1_LOCK_SHA256 = "4365fcd166506a05c3ef4bc5ad6894293885e37ef799352816e2f75dc3793719"
LOCK_SHA256 = "b1b9362b067675711e6facc11e2176701ce663dedb86d831a4ca6beb3079e64c"
EXPECTED_PROTOCOL_SHA256 = "baafe7a774bff90323d3739739044bb559e1938b76bb4915ef8a33348ca9b6dc"
EXPECTED_PROTOCOL_FILE_SHA256 = "b7d38807a2b027c929412b81639bc531f39f4c49a6cd20a6069b89b8c9afa13d"
EXPECTED_AMENDMENT_SHA256 = "b6e39fefc60981f3ac53fa37bd35bf18a0c25940f998b25d4f025d7f2bbeca9c"
EXPECTED_AMENDMENT_FILE_SHA256 = "f11ae681eece162d3ea83f0379af0d8fa5596d6a8dd1ccfe4bad49f2e750ef2c"
PROTOCOL_RELATIVE = Path("research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-18-degree23-panel-protocol.json")
AMENDMENT_RELATIVE = Path("research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-18-amendment-01-lock-correction.json")
CORRECTED_LOCK_RELATIVE = Path("research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-18-corrected-Cargo.lock")
V1_LOCK_RELATIVE = Path("research/sat_factor_base_review_20260908/continuation-01/source_snapshots/Cargo.lock")
FAILED_V1_PANEL_RELATIVE = Path("research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-18-degree23-panel-20260909")
CORRECTED_PANEL_RELATIVE = Path("research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-18-degree23-panel-lock-corrected-20260909")
RUNNER_RELATIVE = Path("research/sat_factor_base_review_20260908/continuation-05-sota-gates/run_stage18_degree23_panel.py")
VERIFIER_RELATIVE = Path("research/sat_factor_base_review_20260908/continuation-05-sota-gates/verify_stage18_degree23_panel.py")
METER_RELATIVE = Path("scripts/process_meter.py")
LOCK_ARCHIVE_RELATIVE = Path("dependency-lock/Cargo.lock")
IC_BINARY_RELATIVE = Path("target/release/examples/koblitz_algebraic_e2e")
DISCOVERY_BINARY_RELATIVE = Path("target/release/examples/koblitz_public_factor_base_discovery")
EXPECTED_CORRECTION_PATHS = [
    RUNNER_RELATIVE,
    AMENDMENT_RELATIVE,
    CORRECTED_LOCK_RELATIVE,
    FAILED_V1_PANEL_RELATIVE / "build/metrics.json",
    FAILED_V1_PANEL_RELATIVE / "build/process-invocation.json",
    FAILED_V1_PANEL_RELATIVE / "build/receipt.json",
    FAILED_V1_PANEL_RELATIVE / "build/stderr.txt",
    FAILED_V1_PANEL_RELATIVE / "build/stdout.txt",
    FAILED_V1_PANEL_RELATIVE / "dependency-lock/Cargo.lock",
    FAILED_V1_PANEL_RELATIVE / "outer-attempts/0001/invocation.json",
    FAILED_V1_PANEL_RELATIVE / "outer-attempts/0001/metrics.json",
    FAILED_V1_PANEL_RELATIVE / "outer-attempts/0001/receipt.json",
    FAILED_V1_PANEL_RELATIVE / "outer-attempts/0001/stderr.txt",
    FAILED_V1_PANEL_RELATIVE / "outer-attempts/0001/stdout.txt",
    FAILED_V1_PANEL_RELATIVE / "protocol.json",
    FAILED_V1_PANEL_RELATIVE / "run.json",
    PROTOCOL_RELATIVE,
    VERIFIER_RELATIVE,
]
EXPECTED_CONTROL_DELTA = [
    {"status": "A", "path": str(path), "mode": "100644", "type": "blob"}
    for path in EXPECTED_CORRECTION_PATHS
]
V1_EXPECTED_CONTROL_DELTA = [
    {"status": "A", "path": str(path), "mode": "100644", "type": "blob"}
    for path in (RUNNER_RELATIVE, PROTOCOL_RELATIVE, VERIFIER_RELATIVE)
]
SOURCE_SHA256 = {
    "Cargo.toml": "43611c79a8692d99fc3146caf339ec72b4b93bf0ee1e2a8d50656558fa8beaa5",
    "examples/koblitz_algebraic_e2e.rs": "1bc7ba2ca6d3803a8438b1f32c990914c64f2167a17ea519bcae53fa1d4be263",
    "examples/koblitz_public_factor_base_discovery.rs": "416ab6ec893dbff54e65df1da6b01e7e4b883a0cbf75821c27165bdbe785ab9b",
    "scripts/process_meter.py": "8705343c10c129941ff3a1068be8211c0a495d2ea9a5200bc26b64ad26e3065f",
    "src/cryptanalysis/koblitz_index_calculus.rs": "cc0cde413adf1179afd1b8cd3d9133e13f422acdff395e0c953a310069a7a681",
    "src/cryptanalysis/sat.rs": "9290dd9bab2b640b0a4ab06dd79a36e99feefb067b637b1543b8c084291d2129",
    "src/cryptanalysis/semaev_sat.rs": "18ac7704b399c06f6036d80fefa05dae7f8265f22654bf65f0eeaea7f81bc4b7",
}

EXPECTED_RUNS = [
    (1, 1871679, 7827814920183983948),
    (2, 1056924, 1300251816495109512),
    (3, 1966632, 226800821669791328),
    (4, 1956519, 10471725251496847725),
    (5, 1785876, 11656057938264181787),
]
EXPECTED_TARGETS = {
    1: {"x": "3881370", "y": "1511041"},
    2: {"x": "6067514", "y": "2806957"},
    3: {"x": "3329510", "y": "2418086"},
    4: {"x": "335048", "y": "3444171"},
    5: {"x": "5239510", "y": "6376352"},
}
DISCOVERY_TIMING_FIELDS = {"curve_and_subgroup_construction", "end_to_end"}
DISCOVERY_CANDIDATE_TIMING_FIELDS = {
    "cofactor_admission", "factor_base_construction", "projection_census",
}
IC_TIMING_FIELDS = {
    "curve_and_subgroup_construction", "target_construction",
    "factor_base_predicate_and_materialisation", "projected_orbit_construction",
    "cofactor_admission", "relation_collection", "linear_algebra", "driver_solve",
    "end_to_end",
}
RHO_TIMING_FIELDS = {"curve_construction", "target_construction", "rho", "end_to_end"}
SEMANTIC_FIELDS = (
    "divisor_indices",
    "divisor_polynomial",
    "linearised_exponents",
    "dimension",
    "abscissae",
    "rational_points",
    "signed_frobenius_orbits_before_projection",
    "projected_points",
    "projected_signed_frobenius_orbits",
    "m_cofactor_admissible",
)
EXPECTED_DISCOVERY = {
    0: {
        "group_order": "8383412", "subgroup_order": "2095853", "cofactor": "4",
        "candidates": [
            ([0, 1], 7973, [0, 2, 5, 8, 9, 10, 11, 12], 4235, 94, 4233, 92),
            ([0, 2], 5279, [0, 1, 2, 3, 4, 7, 10, 12], 4281, 95, 4279, 93),
        ],
        "selected": [0, 2],
    },
    1: {
        "group_order": "8393806", "subgroup_order": "4196903", "cofactor": "2",
        "candidates": [
            ([0, 1], 7973, [0, 2, 5, 8, 9, 10, 11, 12], 3957, 87, 3957, 86),
            ([0, 2], 5279, [0, 1, 2, 3, 4, 7, 10, 12], 3911, 86, 3911, 85),
        ],
        "selected": [0, 1],
    },
}


class VerificationError(RuntimeError):
    pass


def _object_no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise VerificationError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(), object_pairs_hook=_object_no_duplicates)
    except VerificationError:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise VerificationError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise VerificationError(f"expected JSON object: {path}")
    return value


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def _canonical_absolute(value: Any, label: str) -> Path:
    require(isinstance(value, str) and value, f"missing {label}")
    path = Path(value)
    require(path.is_absolute() and ".." not in path.parts, f"noncanonical {label}")
    require(str(path) == value, f"noncanonical {label}")
    return path


def _execution_roots(run: dict) -> tuple[Path, Path]:
    roots = run.get("execution_roots")
    require(
        isinstance(roots, dict) and set(roots) == {"repository", "panel"},
        "run lacks exact execution roots",
    )
    repository = _canonical_absolute(roots.get("repository"), "execution repository root")
    panel = _canonical_absolute(roots.get("panel"), "execution panel root")
    require(repository != panel, "execution repository and panel roots collide")
    return repository, panel


def _recorded_relative(value: Any, root: Path, expected: Path, label: str) -> Path:
    path = _canonical_absolute(value, label)
    try:
        relative = path.relative_to(root)
    except ValueError as error:
        raise VerificationError(f"{label} escapes its recorded root") from error
    require(relative == expected, f"{label} changed: {relative}")
    return relative


def _identity_fields(identity: Any, label: str, expected_sha256: str | None = None) -> tuple[str, int, str]:
    require(isinstance(identity, dict), f"missing file identity: {label}")
    path = identity.get("path")
    size = identity.get("bytes")
    digest = identity.get("sha256")
    require(isinstance(path, str), f"missing file-identity path: {label}")
    require(isinstance(size, int) and not isinstance(size, bool) and size >= 0, f"invalid file size: {label}")
    require(
        isinstance(digest, str) and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest),
        f"invalid file hash: {label}",
    )
    if expected_sha256 is not None:
        require(digest == expected_sha256, f"file identity differs from frozen hash: {label}")
    return path, size, digest


def _hex_identifier(value: Any, label: str, length: int = 40) -> str:
    require(
        isinstance(value, str) and len(value) == length
        and all(character in "0123456789abcdef" for character in value),
        f"invalid {label}",
    )
    return value


def git_tracked_delta(repo: Path, commit: str) -> dict:
    """Return the exact tracked tree delta from the frozen algorithm base."""
    _hex_identifier(commit, "tracked-delta commit")
    completed = subprocess.run(
        [
            "git", "diff", "--name-status", "--no-renames", "-z",
            ALGORITHM_BASE_COMMIT, commit, "--",
        ],
        cwd=repo, check=False, capture_output=True,
    )
    require(completed.returncode == 0, "cannot inspect tracked source delta")
    fields = completed.stdout.split(b"\0")
    require(fields and fields[-1] == b"", "malformed tracked source delta")
    fields.pop()
    require(len(fields) % 2 == 0, "malformed tracked source delta fields")
    entries = []
    for offset in range(0, len(fields), 2):
        try:
            status = fields[offset].decode("ascii")
            path = fields[offset + 1].decode("utf-8")
        except UnicodeDecodeError as error:
            raise VerificationError("tracked source delta has a non-UTF-8 field") from error
        require(status and "\t" not in status and path and "\x00" not in path, "invalid tracked source delta entry")
        tree = subprocess.run(
            ["git", "ls-tree", "-z", commit, "--", path],
            cwd=repo, check=False, capture_output=True,
        )
        require(tree.returncode == 0, f"cannot inspect tracked tree entry: {path}")
        mode = object_type = object_id = None
        if tree.stdout:
            require(tree.stdout.endswith(b"\0"), f"malformed tracked tree entry: {path}")
            try:
                metadata, tree_path = tree.stdout[:-1].decode("utf-8").split("\t", 1)
                mode, object_type, object_id = metadata.split(" ", 2)
            except (UnicodeDecodeError, ValueError) as error:
                raise VerificationError(f"malformed tracked tree metadata: {path}") from error
            require(tree_path == path, f"tracked tree path mismatch: {path}")
            _hex_identifier(object_id, f"tracked tree object id: {path}")
        entries.append({
            "status": status, "path": path, "mode": mode,
            "type": object_type, "object_id": object_id,
        })
    projection = [
        {key: entry[key] for key in ("status", "path", "mode", "type")}
        for entry in entries
    ]
    return {
        "schema": "koblitz_stage18_tracked_delta.v1",
        "base_commit": ALGORITHM_BASE_COMMIT,
        "head_commit": commit,
        "entries": entries,
        "matches_expected_control_only_delta": projection == EXPECTED_CONTROL_DELTA,
    }


def require_expected_control_delta(delta: Any) -> None:
    require(isinstance(delta, dict), "missing tracked source delta")
    require(delta.get("schema") == "koblitz_stage18_tracked_delta.v1", "wrong tracked-delta schema")
    require(delta.get("base_commit") == ALGORITHM_BASE_COMMIT, "tracked-delta base changed")
    _hex_identifier(delta.get("head_commit"), "tracked-delta head")
    entries = delta.get("entries")
    require(isinstance(entries, list), "tracked-delta entries are missing")
    projection = []
    for entry in entries:
        require(
            isinstance(entry, dict)
            and set(entry) == {"status", "path", "mode", "type", "object_id"},
            "invalid tracked-delta entry",
        )
        _hex_identifier(entry.get("object_id"), f"tracked-delta object: {entry.get('path')}")
        projection.append({key: entry[key] for key in ("status", "path", "mode", "type")})
    require(
        projection == EXPECTED_CONTROL_DELTA,
        "execution commit differs from the exact Stage 18 correction custody bundle",
    )
    require(delta.get("matches_expected_control_only_delta") is True, "tracked-delta admissibility flag changed")


def _validate_historical_identity(
    identity: Any,
    recorded_root: Path,
    expected_relative: Path,
    label: str,
    expected_sha256: str | None = None,
) -> None:
    path, _, _ = _identity_fields(identity, label, expected_sha256)
    _recorded_relative(path, recorded_root, expected_relative, label)


def _validate_mapped_identity(
    identity: Any,
    recorded_root: Path,
    expected_relative: Path,
    live_path: Path,
    label: str,
    expected_sha256: str | None = None,
) -> None:
    path, size, digest = _identity_fields(identity, label, expected_sha256)
    _recorded_relative(path, recorded_root, expected_relative, label)
    require(live_path.is_file() and not live_path.is_symlink(), f"mapped file is not regular: {label}")
    require(size == live_path.stat().st_size, f"mapped file size changed: {label}")
    actual = sha256_file(live_path)
    require(digest == actual, f"mapped file hash changed: {label}")
    if expected_sha256 is not None:
        require(actual == expected_sha256, f"mapped file differs from protocol hash: {label}")


def _timestamp(value: Any, label: str) -> datetime:
    require(isinstance(value, str) and value.endswith("Z"), f"invalid {label}")
    try:
        return datetime.fromisoformat(value[:-1] + "+00:00")
    except ValueError as error:
        raise VerificationError(f"invalid {label}") from error


def _validate_process_invocation(invocation: dict, command: list[str], cwd: Path,
                                 watchdog: float, meter_path: str, stdout: Path,
                                 stderr: Path, metrics: Path, python_path: str) -> None:
    require(invocation.get("cwd") == str(cwd), "process invocation cwd changed")
    require(invocation.get("command") == command, "process invocation command changed")
    require(invocation.get("watchdog_seconds") == watchdog, "process invocation watchdog changed")
    require(invocation.get("meter_launcher_returncode") == 0, "process meter launcher failed")
    observed = invocation.get("driver_observed_wall_seconds")
    require(isinstance(observed, (int, float)) and math.isfinite(observed) and observed >= 0, "invalid driver-observed wall")
    require(_timestamp(invocation.get("finished_at"), "finished_at") >= _timestamp(invocation.get("started_at"), "started_at"), "process timestamps reversed")
    expected_meter = [
        python_path, meter_path, "--cwd", str(cwd), "--timeout", str(float(watchdog)),
        "--stdout", str(stdout), "--stderr", str(stderr),
        "--metrics", str(metrics), "--", *command,
    ]
    require(invocation.get("meter_command") == expected_meter, "process meter invocation changed")


def derived_run(index: int, subgroup_order: int = 2095853) -> tuple[int, int, str, str]:
    secret_digest = hashlib.sha256(f"crypto-koblitz-stage18-secret-v1:{index}".encode()).digest()
    seed_digest = hashlib.sha256(f"crypto-koblitz-stage18-seed-v1:{index}".encode()).digest()
    secret = 1 + int.from_bytes(secret_digest, "big") % (subgroup_order - 1)
    seed = int.from_bytes(seed_digest[:8], "big")
    return secret, seed, secret_digest.hex(), seed_digest.hex()


def validate_protocol(protocol: dict) -> None:
    require(
        canonical_sha256(protocol) == EXPECTED_PROTOCOL_SHA256,
        "Stage 18 protocol differs from the exact frozen canonical digest",
    )
    require(protocol.get("schema") == PROTOCOL_SCHEMA, "wrong Stage 18 protocol schema")
    require(protocol.get("status") == "frozen_before_execution", "protocol is not frozen before execution")
    binding = protocol.get("implementation_binding")
    require(isinstance(binding, dict), "protocol lacks implementation binding")
    require(binding.get("algorithm_base_commit") == ALGORITHM_BASE_COMMIT, "algorithm base commit changed")
    require(binding.get("source_sha256") == SOURCE_SHA256, "protocol source-hash binding changed")
    require(
        binding.get("tracked_lock_snapshot")
        == "research/sat_factor_base_review_20260908/continuation-01/source_snapshots/Cargo.lock",
        "tracked lock path changed",
    )
    require(binding.get("tracked_lock_sha256") == V1_LOCK_SHA256, "v1 tracked lock hash changed")
    require(binding.get("root_lock_copy_required_for_locked_build") is True, "root lock-copy policy changed")
    curve = protocol.get("curve", {})
    require(curve == {
        "n": 23, "a": 0, "group_order": "8383412", "subgroup_order": "2095853",
        "cofactor": "4", "summands": 2,
    }, "frozen curve parameters changed")
    factor = protocol.get("factor_base", {})
    expected_factor = {
        "factor_spec": "divisor:0,2", "divisor_indices": [0, 2],
        "divisor_polynomial": 5279, "dimension": 12,
        "linearised_exponents": [0, 1, 2, 3, 4, 7, 10, 12], "abscissae": 4096,
        "rational_points": 4281, "signed_frobenius_orbits_before_projection": 95,
        "projected_signed_frobenius_columns": 93,
        "enumerates_target_subgroup": False, "uses_discrete_log_labels": False,
        "factor_base_logs_constructed": False,
    }
    for key, value in expected_factor.items():
        require(factor.get(key) == value, f"frozen factor-base field changed: {key}")
    runs = protocol.get("runs")
    require(isinstance(runs, list) and len(runs) == 5, "protocol must contain five runs")
    observed = []
    for item in runs:
        index = item.get("index")
        require(isinstance(index, int), "run index must be an integer")
        secret, seed, secret_hash, seed_hash = derived_run(index)
        require(item.get("secret") == secret, f"run {index} secret derivation mismatch")
        require(item.get("seed") == str(seed), f"run {index} seed derivation mismatch")
        require(item.get("secret_sha256") == secret_hash, f"run {index} secret hash mismatch")
        require(item.get("seed_sha256") == seed_hash, f"run {index} seed hash mismatch")
        require(item.get("public_target") == EXPECTED_TARGETS.get(index), f"run {index} public target changed")
        observed.append((index, secret, seed))
    require(observed == EXPECTED_RUNS, "frozen run ordering changed")
    require(len({secret for _, secret, _ in observed}) == 5, "frozen secrets are not unique")
    execution = protocol.get("execution", {})
    require(execution.get("inner_process_watchdog_seconds") == 600, "inner watchdog changed")
    require(execution.get("whole_driver_watchdog_seconds") == 7200, "outer watchdog changed")
    require(execution.get("in_stage_retries") == 0, "in-stage retry policy changed")
    require(execution.get("expected_inner_processes") == 12, "expected process count changed")
    ic = protocol.get("index_calculus", {})
    for key, value in {
        "conflict_budget_per_target": 100000, "max_trials": 1000,
        "max_models_per_target": 64, "relation_batch_size": 1, "parallel_threads": 1,
        "direct_relation_forbidden": True, "invalid_models_forbidden": True,
        "cofactor_projection_column_merge": True, "stop_on_verified_rank": True,
    }.items():
        require(ic.get(key) == value, f"index-calculus policy changed: {key}")
    rho = protocol.get("rho", {})
    require(rho.get("automorphism_group_bound") == 46, "rho automorphism bound changed")
    require(rho.get("fresh_deterministic_jump_table_per_restart") is True, "rho jump refresh changed")


def validate_lock_correction() -> None:
    old_path = REPO / V1_LOCK_RELATIVE
    corrected_path = REPO / CORRECTED_LOCK_RELATIVE
    require(
        old_path.is_file() and not old_path.is_symlink()
        and sha256_file(old_path) == V1_LOCK_SHA256,
        "v1 dependency lock snapshot changed",
    )
    require(
        corrected_path.is_file() and not corrected_path.is_symlink()
        and sha256_file(corrected_path) == LOCK_SHA256,
        "corrected dependency lock snapshot changed",
    )
    old_bytes = old_path.read_bytes()
    corrected_bytes = corrected_path.read_bytes()
    insertion_point = b' "hex",\n'
    require(old_bytes.count(insertion_point) == 1, "v1 root dependency insertion point changed")
    require(
        corrected_bytes == old_bytes.replace(
            insertion_point, insertion_point + b' "libc",\n', 1,
        ),
        "corrected lock is not the exact one-line libc insertion",
    )
    try:
        old = tomllib.loads(old_bytes.decode("utf-8"))
        corrected = tomllib.loads(corrected_bytes.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise VerificationError(f"cannot parse dependency lock correction: {error}") from error
    require(old.get("version") == corrected.get("version") == 4, "Cargo lock format version changed")
    old_packages = old.get("package")
    corrected_packages = corrected.get("package")
    require(
        isinstance(old_packages, list) and isinstance(corrected_packages, list)
        and len(old_packages) == len(corrected_packages) == 68,
        "Cargo package inventory changed",
    )
    old_root = [package for package in old_packages if package.get("name") == "crypto"]
    corrected_root = [package for package in corrected_packages if package.get("name") == "crypto"]
    require(len(old_root) == len(corrected_root) == 1, "root crypto package record changed")
    expected_root = copy.deepcopy(old_root[0])
    dependencies = list(expected_root.get("dependencies", []))
    require("libc" not in dependencies and dependencies.count("hex") == 1, "v1 root dependency list changed")
    dependencies.insert(dependencies.index("hex") + 1, "libc")
    expected_root["dependencies"] = dependencies
    require(corrected_root[0] == expected_root, "corrected root package dependency list changed")
    require(
        [package for package in corrected_packages if package.get("name") != "crypto"]
        == [package for package in old_packages if package.get("name") != "crypto"],
        "non-root Cargo package record changed",
    )
    libc_records = [package for package in corrected_packages if package.get("name") == "libc"]
    require(libc_records == [{
        "name": "libc", "version": "0.2.186",
        "source": "registry+https://github.com/rust-lang/crates.io-index",
        "checksum": "68ab91017fe16c622486840e4c83c9a37afeff978bd239b5293d61ece587de66",
    }], "existing libc package identity changed")


def validate_amendment(amendment: dict, protocol: dict) -> None:
    validate_protocol(protocol)
    validate_lock_correction()
    require(
        canonical_sha256(amendment) == EXPECTED_AMENDMENT_SHA256,
        "Stage 18 lock amendment differs from the exact frozen canonical digest",
    )
    require(
        sha256_file(DEFAULT_AMENDMENT) == EXPECTED_AMENDMENT_FILE_SHA256,
        "tracked Stage 18 lock amendment bytes changed",
    )
    require(
        amendment.get("schema") == "koblitz_degree23_replication_panel_amendment.v1"
        and amendment.get("amendment_id") == "STAGE18-AMENDMENT-01-LOCK-CORRECTION"
        and amendment.get("status") == "frozen_before_corrected_execution",
        "wrong Stage 18 lock amendment identity",
    )
    require(sha256_file(DEFAULT_PROTOCOL) == EXPECTED_PROTOCOL_FILE_SHA256, "v1 protocol bytes changed")
    base = amendment.get("base_protocol", {})
    require(base == {
        "path": str(PROTOCOL_RELATIVE),
        "file_sha256": EXPECTED_PROTOCOL_FILE_SHA256,
        "canonical_sha256": EXPECTED_PROTOCOL_SHA256,
    }, "amendment base-protocol binding changed")
    fields = amendment.get("supersession", {}).get("fields")
    require(fields == [
        {
            "json_pointer": "/implementation_binding/tracked_lock_snapshot",
            "prior": protocol["implementation_binding"]["tracked_lock_snapshot"],
            "replacement": str(CORRECTED_LOCK_RELATIVE),
        },
        {
            "json_pointer": "/implementation_binding/tracked_lock_sha256",
            "prior": protocol["implementation_binding"]["tracked_lock_sha256"],
            "replacement": LOCK_SHA256,
        },
    ], "amendment lock supersession changed")
    require(
        amendment.get("supersession", {}).get("all_other_protocol_fields_remain_byte_identical") is True,
        "amendment does not preserve the rest of v1",
    )
    require(
        amendment.get("corrected_execution", {}).get("required_output_path")
        == str(CORRECTED_PANEL_RELATIVE),
        "amendment corrected output path changed",
    )
    require(
        amendment.get("corrected_execution", {}).get("resume_or_retry_forbidden") is True,
        "amendment retry boundary changed",
    )
    require(
        amendment.get("failed_execution", {}).get("panel_path") == str(FAILED_V1_PANEL_RELATIVE)
        and amendment.get("failed_execution", {}).get("scientific_tasks_started") == 0
        and amendment.get("failed_execution", {}).get("resume_or_retry_forbidden") is True,
        "amendment failed-run boundary changed",
    )
    require(
        CORRECTED_LOCK_RELATIVE.is_absolute() is False
        and (REPO / CORRECTED_LOCK_RELATIVE).is_file()
        and not (REPO / CORRECTED_LOCK_RELATIVE).is_symlink()
        and sha256_file(REPO / CORRECTED_LOCK_RELATIVE) == LOCK_SHA256,
        "corrected Stage 18 lock snapshot is missing or changed",
    )


def task_plan(protocol: dict, discovery_binary: str, ic_binary: str) -> list[dict]:
    validate_protocol(protocol)
    tasks = []
    for a in (0, 1):
        tasks.append({
            "id": f"discovery-a{a}", "kind": "discovery", "curve_a": a,
            "command": [discovery_binary, "23", str(a), "2", "12"],
        })
    for index, secret, seed in EXPECTED_RUNS:
        for mode in ("ic", "rho-auto"):
            tasks.append({
                "id": f"row-{index:02d}-{mode}", "kind": mode, "index": index,
                "secret": secret, "seed": seed, "public_target": EXPECTED_TARGETS[index],
                "command": [
                    ic_binary, mode, "23", "0", "2", str(secret), str(seed),
                    "100000", "1000", "divisor:0,2",
                ],
            })
    require(len(tasks) == 12, "internal task-plan error")
    return tasks


def _semantic_candidate(candidate: dict) -> dict:
    return {key: candidate.get(key) for key in SEMANTIC_FIELDS}


def _expected_candidate(row: tuple) -> dict:
    indices, polynomial, exponents, points, orbits, projected_points, projected_orbits = row
    return {
        "divisor_indices": indices, "divisor_polynomial": polynomial,
        "linearised_exponents": exponents, "dimension": 12, "abscissae": 4096,
        "rational_points": points, "signed_frobenius_orbits_before_projection": orbits,
        "projected_points": projected_points,
        "projected_signed_frobenius_orbits": projected_orbits,
        "m_cofactor_admissible": True,
    }


def validate_discovery_result(result: dict, curve_a: int) -> None:
    expected = EXPECTED_DISCOVERY[curve_a]
    require(result.get("schema") == "koblitz_public_factor_base_discovery.v1", "wrong discovery schema")
    for key, value in {
        "n": 23, "a": curve_a, "m": 2, "requested_dimension": 12,
        "group_order": expected["group_order"], "subgroup_order": expected["subgroup_order"],
        "cofactor": expected["cofactor"], "factor_degrees": [1, 11, 11],
    }.items():
        require(result.get(key) == value, f"discovery a={curve_a} changed {key}")
    require(result.get("sizing_gate") == {"m_times_dimension": 24, "n": 23, "passes": True}, "bad sizing gate")
    forbidden = result.get("forbidden_inputs")
    require(isinstance(forbidden, dict) and set(forbidden) == {
        "target_constructed", "target_subgroup_enumerated", "discrete_log_labels_constructed",
        "relation_yield_used", "solver_timing_used",
    } and not any(forbidden.values()), "discovery used a forbidden input")
    candidates = result.get("candidates")
    require(isinstance(candidates, list) and len(candidates) == 2, "wrong discovery candidate inventory")
    require([_semantic_candidate(item) for item in candidates] == [
        _expected_candidate(item) for item in expected["candidates"]
    ], f"discovery a={curve_a} algebra changed")
    selected = result.get("selected")
    require(isinstance(selected, dict), "missing selected discovery candidate")
    selected_expected = next(
        _expected_candidate(item) for item in expected["candidates"] if item[0] == expected["selected"]
    )
    require(_semantic_candidate(selected) == selected_expected, f"discovery a={curve_a} selected candidate changed")
    validate_timing_tree(result.get("timing_ns"), "discovery timing", DISCOVERY_TIMING_FIELDS)
    for item in candidates:
        validate_timing_tree(
            item.get("timing_ns"), "candidate timing", DISCOVERY_CANDIDATE_TIMING_FIELDS,
        )


def validate_timing_tree(value: Any, label: str, expected_fields: set[str]) -> None:
    require(isinstance(value, dict) and value, f"missing {label}")
    require(set(value) == expected_fields, f"{label} field inventory changed")
    for key, number in value.items():
        require(isinstance(key, str), f"invalid {label} key")
        require(isinstance(number, int) and not isinstance(number, bool) and number >= 0, f"invalid {label}.{key}")


def validate_target(target: Any, expected: dict, label: str) -> None:
    require(isinstance(target, dict) and set(target) == {"x", "y"}, f"invalid {label} target")
    for coordinate in ("x", "y"):
        value = target.get(coordinate)
        require(
            isinstance(value, str) and value.isascii() and value.isdecimal()
            and str(int(value)) == value and 0 <= int(value) < (1 << 23),
            f"invalid {label} target coordinate: {coordinate}",
        )
    require(target == expected, f"{label} target does not equal the frozen [secret]G point")


def validate_ic_result(result: dict, task: dict) -> None:
    require(result.get("kind") == "koblitz_unknown_scalar_algebraic_factor_base_e2e", "wrong IC kind")
    for key, value in {"mode": "ic", "n": 23, "a": 0, "m": 2, "subgroup_order": "2095853", "seed": task["seed"]}.items():
        require(result.get(key) == value, f"IC changed {key}")
    predicate = result.get("factor_base_predicate", {})
    for key, value in {
        "kind": "ggmp_divisor_kernel", "factor_spec": "divisor:0,2",
        "divisor_indices": [0, 2], "factor_bitmask": 5279,
        "linearised_exponents": [0, 1, 2, 3, 4, 7, 10, 12],
        "enumerates_target_subgroup": False, "uses_discrete_log_labels": False,
        "factor_base_logs_constructed": False,
    }.items():
        require(predicate.get(key) == value, f"IC factor-base field changed: {key}")
    factor = result.get("factor_base", {})
    for key, value in {
        "ell": 12, "points": 4281, "signed_frobenius_orbits_before_projection": 95,
        "projected_signed_frobenius_orbits": 93, "m_cofactor_admissible": True,
    }.items():
        require(factor.get(key) == value, f"IC factor-base result changed: {key}")
    require(result.get("options") == {
        "conflict_budget_per_target": 100000, "max_trials": 1000,
        "max_models": 64, "parallel_threads": 1,
    }, "IC options changed")
    report = result.get("report", {})
    require(report.get("verified_unknown_scalar_recovery") is True, "IC scalar is not verified")
    require(report.get("recovered_scalar") == str(task["secret"]), "IC recovered the wrong scalar")
    require(report.get("direct_relation") is False, "IC used a direct relation")
    require(report.get("sat_invalid_models") == 0, "IC admitted an invalid SAT model")
    require(report.get("collapse_projected_orbits") is True, "IC projection collapse disabled")
    require(report.get("m_cofactor_admissible") is True, "IC cofactor gate failed")
    counters = ("relations", "trials", "relation_batches", "direct_relations_skipped", "linear_solve_attempts",
                "sat_calls", "sat_models", "sat_refutations", "sat_unknowns", "sat_invalid_models", "sat_conflicts")
    for key in counters:
        require(isinstance(report.get(key), int) and report[key] >= 0, f"invalid IC counter: {key}")
    require(0 < report["trials"] <= 1000, "IC trial budget violated")
    require(report["relation_batches"] == report["trials"], "IC one-target batch/trial mismatch")
    require(
        report["sat_calls"] + report["direct_relations_skipped"] == report["trials"],
        "IC SAT/direct-skip trial counts do not balance",
    )
    require(report["sat_models"] + report["sat_refutations"] + report["sat_unknowns"] == report["sat_calls"], "IC SAT terminal counts do not balance")
    require(report["relations"] <= report["sat_models"], "IC relations exceed SAT models")
    require(report["sat_conflicts"] <= 100000 * report["sat_calls"], "IC conflict accounting exceeds cap")
    require(report["linear_solve_attempts"] >= 1, "IC never attempted linear algebra")
    validate_target(result.get("target"), task["public_target"], "IC")
    validate_timing_tree(result.get("timing_ns"), "IC timing", IC_TIMING_FIELDS)


def validate_rho_result(result: dict, task: dict) -> None:
    require(result.get("kind") == "koblitz_unknown_scalar_automorphism_rho_control", "wrong rho kind")
    for key, value in {"mode": "rho-auto", "n": 23, "a": 0, "subgroup_order": "2095853", "seed": task["seed"]}.items():
        require(result.get(key) == value, f"rho changed {key}")
    require(result.get("verified_unknown_scalar_recovery") is True, "rho scalar is not verified")
    require(result.get("automorphism_optimized") is True, "rho is not automorphism optimized")
    require(result.get("automorphism_group_bound") == 46, "rho automorphism bound changed")
    for key in ("iterations", "walk_step_additions", "restarts", "jump_table_additions",
                "initial_state_additions", "reported_group_additions", "rho_setup_scalar_multiplications",
                "target_construction_scalar_multiplications"):
        require(isinstance(result.get(key), int) and result[key] >= 0, f"invalid rho counter: {key}")
    attempts = result["restarts"] + 1
    require(1 <= attempts <= 64, "rho restart budget violated")
    require(result["iterations"] >= attempts, "rho iterations are fewer than completed attempts")
    require(result["iterations"] <= attempts * (1 << 28), "rho iteration budget violated")
    require(result["walk_step_additions"] == 3 * result["iterations"], "rho walk additions do not balance iterations")
    require(result["jump_table_additions"] == 16 * attempts, "rho jump additions do not balance")
    require(result["initial_state_additions"] == attempts, "rho initial additions do not balance")
    require(result["rho_setup_scalar_multiplications"] == 34 * attempts, "rho setup scalar multiplications do not balance")
    require(result["reported_group_additions"] == result["walk_step_additions"] + 17 * attempts, "rho group additions do not balance")
    require(result["target_construction_scalar_multiplications"] == 1, "rho target construction count changed")
    validate_target(result.get("target"), task["public_target"], "rho")
    validate_timing_tree(result.get("timing_ns"), "rho timing", RHO_TIMING_FIELDS)


def validate_meter(meter: dict, expected_command: list[str], watchdog: float) -> tuple[dict, bool]:
    require(meter.get("command") == expected_command, "meter command differs from frozen task")
    require(meter.get("watchdog_seconds") == float(watchdog), "meter watchdog changed")
    require(isinstance(meter.get("returncode"), int), "meter return code is missing")
    require(isinstance(meter.get("timed_out"), bool), "meter timeout flag is missing")
    require(isinstance(meter.get("orphan_group_terminated"), bool), "meter orphan flag is missing")
    metrics = meter.get("metrics")
    require(isinstance(metrics, dict), "meter resource record is missing")
    for key in ("wall_seconds", "user_seconds", "system_seconds", "total_core_seconds", "single_core_seconds"):
        number = metrics.get(key)
        require(isinstance(number, (int, float)) and not isinstance(number, bool) and math.isfinite(number) and number >= 0, f"invalid meter metric: {key}")
    rss = metrics.get("peak_rss_bytes")
    require(isinstance(rss, int) and not isinstance(rss, bool) and rss > 0, "invalid peak RSS")
    require(metrics.get("meter") == "fresh-process getrusage(RUSAGE_CHILDREN)", "unexpected meter implementation")
    require(math.isclose(metrics["total_core_seconds"], metrics["user_seconds"] + metrics["system_seconds"], rel_tol=0, abs_tol=1e-9), "meter CPU fields do not balance")
    require(math.isclose(metrics["single_core_seconds"], metrics["total_core_seconds"], rel_tol=0, abs_tol=1e-12), "single-core field changed")
    success = meter["returncode"] == 0 and not meter["timed_out"] and not meter["orphan_group_terminated"]
    return metrics, success


def validate_failed_v1_custody(amendment: dict, protocol: dict, panel: Path = FAILED_V1_PANEL) -> dict:
    """Validate the immutable zero-scientific-task v1 build failure."""
    validate_amendment(amendment, protocol)
    require(panel.is_dir() and not panel.is_symlink(), "retained v1 failure panel is missing")
    inventory = amendment["failed_execution"].get("artifact_inventory")
    require(isinstance(inventory, list) and len(inventory) == 13, "failed-panel inventory changed")
    listed = set()
    for row in inventory:
        require(
            isinstance(row, dict) and set(row) == {"path", "bytes", "sha256"},
            "invalid failed-panel inventory row",
        )
        relative = Path(row["path"])
        require(not relative.is_absolute() and ".." not in relative.parts, "unsafe failed-panel path")
        artifact = panel / relative
        require(artifact.is_file() and not artifact.is_symlink(), f"failed-panel artifact is not regular: {relative}")
        require(str(relative) not in listed, f"duplicate failed-panel artifact: {relative}")
        listed.add(str(relative))
        require(
            artifact.stat().st_size == row["bytes"] and sha256_file(artifact) == row["sha256"],
            f"failed-panel artifact changed: {relative}",
        )
    observed = set()
    for artifact in panel.rglob("*"):
        if artifact.is_dir():
            require(not artifact.is_symlink(), f"symlink directory in failed panel: {artifact}")
            continue
        require(artifact.is_file() and not artifact.is_symlink(), f"non-regular failed-panel artifact: {artifact}")
        observed.add(str(artifact.relative_to(panel)))
    require(observed == listed, "retained v1 failure root inventory changed")

    run = read_json(panel / "run.json")
    recorded_repo, recorded_panel = _execution_roots(run)
    require(
        recorded_panel == recorded_repo / FAILED_V1_PANEL_RELATIVE,
        "failed v1 recorded output path changed",
    )
    require(run.get("schema") == RUN_SCHEMA, "failed v1 run schema changed")
    require(run.get("protocol_sha256") == EXPECTED_PROTOCOL_SHA256, "failed v1 protocol hash changed")
    require(run.get("protocol_file_sha256") == EXPECTED_PROTOCOL_FILE_SHA256, "failed v1 protocol bytes changed")
    require(read_json(panel / "protocol.json") == protocol, "failed v1 protocol copy changed")
    require(run.get("evidence_class") == "scientific_candidate", "failed v1 evidence class changed")
    require(run.get("status") == "running", "failed v1 pre-finalization status changed")
    require(run.get("tasks") == {}, "failed v1 unexpectedly contains scientific tasks")
    revision = run.get("source_revision", {})
    failed_commit = amendment["failed_execution"]["source_commit"]
    require(revision.get("commit") == failed_commit and revision.get("dirty") is False, "failed v1 source revision changed")
    recomputed = git_tracked_delta(REPO, failed_commit)
    projection = [
        {key: entry[key] for key in ("status", "path", "mode", "type")}
        for entry in recomputed["entries"]
    ]
    require(projection == V1_EXPECTED_CONTROL_DELTA, "failed v1 execution commit delta changed")
    recorded_delta = revision.get("tracked_delta", {})
    require(
        recorded_delta.get("base_commit") == ALGORITHM_BASE_COMMIT
        and recorded_delta.get("head_commit") == failed_commit
        and recorded_delta.get("entries") == recomputed["entries"]
        and recorded_delta.get("matches_expected_control_only_delta") is True,
        "failed v1 recorded source delta changed",
    )
    implementation = run.get("implementation", {})
    require(set(implementation) == {"runner", "verifier", "meter"}, "failed v1 implementation inventory changed")
    _validate_historical_identity(implementation["runner"], recorded_repo, RUNNER_RELATIVE, "failed v1 runner")
    _validate_historical_identity(implementation["verifier"], recorded_repo, VERIFIER_RELATIVE, "failed v1 verifier")
    _validate_historical_identity(
        implementation["meter"], recorded_repo, METER_RELATIVE, "failed v1 meter",
        SOURCE_SHA256["scripts/process_meter.py"],
    )
    require(sha256_file(panel / LOCK_ARCHIVE_RELATIVE) == V1_LOCK_SHA256, "failed v1 archived lock changed")

    build = panel / "build"
    _regular_inventory(build, {"process-invocation.json", "stdout.txt", "stderr.txt", "metrics.json", "receipt.json"})
    build_receipt = read_json(build / "receipt.json")
    build_command = build_receipt.get("command")
    require(
        isinstance(build_command, list) and len(build_command) == 8
        and Path(build_command[0]).name == "cargo"
        and build_command[1:] == [
            "build", "--release", "--locked", "--example", "koblitz_algebraic_e2e",
            "--example", "koblitz_public_factor_base_discovery",
        ],
        "failed v1 build command changed",
    )
    build_meter = read_json(build / "metrics.json")
    build_metrics, build_success = validate_meter(build_meter, build_command, BUILD_WATCHDOG_SECONDS)
    require(not build_success and build_meter["returncode"] == 101, "failed v1 build terminal changed")
    require(build_receipt == {
        "schema": "koblitz_degree23_build_receipt.v1", "status": "failed",
        "command": build_command, "dependency_lock_sha256": V1_LOCK_SHA256,
        "binaries": {},
    }, "failed v1 build receipt changed")
    build_invocation = read_json(build / "process-invocation.json")
    _validate_process_invocation(
        build_invocation, build_command, recorded_repo, BUILD_WATCHDOG_SECONDS,
        implementation["meter"]["path"], recorded_panel / ".staging-build/stdout.txt",
        recorded_panel / ".staging-build/stderr.txt", recorded_panel / ".staging-build/metrics.json",
        run["host"]["python"]["command"][0],
    )

    outer = panel / "outer-attempts/0001"
    _regular_inventory(outer, {"invocation.json", "stdout.txt", "stderr.txt", "metrics.json", "receipt.json"})
    outer_invocation = read_json(outer / "invocation.json")
    outer_command = outer_invocation.get("command")
    require(outer_command == [
        run["host"]["python"]["command"][0], implementation["runner"]["path"],
        "--inner", "--protocol", str(recorded_repo / PROTOCOL_RELATIVE),
        "--output", str(recorded_panel), "--meter", implementation["meter"]["path"],
    ], "failed v1 outer command changed")
    _validate_process_invocation(
        outer_invocation, outer_command, recorded_repo, 7200.0, implementation["meter"]["path"],
        recorded_panel / "outer-attempts/.staging-0001/stdout.txt",
        recorded_panel / "outer-attempts/.staging-0001/stderr.txt",
        recorded_panel / "outer-attempts/.staging-0001/metrics.json",
        run["host"]["python"]["command"][0],
    )
    outer_meter = read_json(outer / "metrics.json")
    outer_metrics, outer_success = validate_meter(outer_meter, outer_command, 7200.0)
    require(not outer_success and outer_meter["returncode"] == 1, "failed v1 outer terminal changed")
    outer_receipt = read_json(outer / "receipt.json")
    require(outer_receipt == {
        "schema": "koblitz_degree23_outer_receipt.v1", "attempt": 1,
        "returncode": 1, "timed_out": False, "orphan_group_terminated": False,
        "metrics": outer_metrics, "command": outer_command,
    }, "failed v1 outer receipt changed")
    require(run.get("outer_attempts") == [{
        "attempt": 1, "receipt": str(recorded_panel / "outer-attempts/0001/receipt.json"),
        "receipt_sha256": sha256_file(outer / "receipt.json"),
    }], "failed v1 outer index changed")

    frozen_build = amendment["failed_execution"]["build"]
    frozen_outer = amendment["failed_execution"]["outer_attempt"]
    for key in ("wall_seconds", "total_core_seconds", "peak_rss_bytes"):
        require(build_metrics[key] == frozen_build[key], f"failed v1 build cost changed: {key}")
        require(outer_metrics[key] == frozen_outer[key], f"failed v1 outer cost changed: {key}")
    return {
        "schema": "koblitz_stage18_v1_failure_custody.v1",
        "status": "verified_operational_build_failure_zero_scientific_tasks",
        "evidence_class": "operational_failure",
        "scientific_tasks_started": 0,
        "artifact_count": len(inventory),
        "build": frozen_build,
        "outer_attempt": frozen_outer,
        "nested_costs_must_not_be_summed": True,
        "excluded_from_ic_rho_ratios": True,
    }


def _regular_inventory(directory: Path, expected: set[str]) -> None:
    require(directory.is_dir() and not directory.is_symlink(), f"task leaf is not a regular directory: {directory}")
    observed = set()
    for path in directory.iterdir():
        require(path.is_file() and not path.is_symlink(), f"non-regular task artifact: {path}")
        observed.add(path.name)
    require(observed == expected, f"task leaf inventory mismatch at {directory}: {sorted(observed)}")


def verify_task_leaf(protocol: dict, panel: Path, task: dict, write_receipt: bool = False) -> dict:
    leaf = panel / "tasks" / task["id"]
    expected_files = {"invocation.json", "result.json", "stderr.txt", "metrics.json"}
    if (leaf / "receipt.json").exists():
        expected_files.add("receipt.json")
    _regular_inventory(leaf, expected_files)
    invocation = read_json(leaf / "invocation.json")
    require(invocation.get("schema") == "koblitz_degree23_task_invocation.v1", "wrong invocation schema")
    require(invocation.get("task") == task, "task invocation differs from frozen plan")
    require(invocation.get("attempt") == 1, "scientific task was retried")
    run = read_json(panel / "run.json")
    recorded_repo, recorded_panel = _execution_roots(run)
    require(
        recorded_panel == recorded_repo / CORRECTED_PANEL_RELATIVE,
        "corrected Stage 18 recorded output path changed",
    )
    recorded_staging = recorded_panel / ".staging" / task["id"]
    python_path = run.get("host", {}).get("python", {}).get("command", [None])[0]
    meter_path = run.get("implementation", {}).get("meter", {}).get("path")
    require(isinstance(python_path, str) and isinstance(meter_path, str), "run lacks meter launcher identity")
    _validate_process_invocation(
        invocation, task["command"], recorded_repo, 600.0, meter_path,
        recorded_staging / "result.json", recorded_staging / "stderr.txt",
        recorded_staging / "metrics.json", python_path,
    )
    meter = read_json(leaf / "metrics.json")
    metrics, process_success = validate_meter(meter, task["command"], 600)
    receipt = {
        "schema": RECEIPT_SCHEMA, "task_id": task["id"], "attempt": 1,
        "command": task["command"], "metrics": metrics,
        "result_sha256": sha256_file(leaf / "result.json"),
        "stderr_sha256": sha256_file(leaf / "stderr.txt"),
        "metrics_sha256": sha256_file(leaf / "metrics.json"),
    }
    if not process_success:
        receipt.update({
            "status": "inconclusive_process_failure", "verified_terminal": False,
            "returncode": meter["returncode"], "timed_out": meter["timed_out"],
            "orphan_group_terminated": meter["orphan_group_terminated"],
        })
    else:
        result = read_json(leaf / "result.json")
        if task["kind"] == "discovery":
            validate_discovery_result(result, task["curve_a"])
        elif task["kind"] == "ic":
            validate_ic_result(result, task)
        else:
            validate_rho_result(result, task)
        internal_ns = result.get("timing_ns", {}).get("end_to_end")
        require(internal_ns is None or internal_ns / 1e9 <= metrics["wall_seconds"] + 0.05, "internal time exceeds metered wall")
        receipt.update({"status": "verified", "verified_terminal": True, "result": result})
    if (leaf / "receipt.json").exists():
        archived = read_json(leaf / "receipt.json")
        require(archived == receipt, f"task receipt mismatch: {task['id']}")
    elif write_receipt:
        temporary = leaf / "receipt.json.tmp"
        temporary.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        temporary.replace(leaf / "receipt.json")
    return receipt


def _distribution(values: list[float | int]) -> dict | None:
    if not values:
        return None
    return {"min": min(values), "median": statistics.median(values), "mean": statistics.fmean(values), "max": max(values)}


def validate_run_custody(protocol: dict, panel: Path, run: dict, require_post: bool,
                         allow_build_failure: bool = False) -> bool:
    require(run.get("schema") == RUN_SCHEMA, "wrong run schema")
    require(run.get("protocol_sha256") == canonical_sha256(protocol), "run protocol hash mismatch")
    recorded_repo, recorded_panel = _execution_roots(run)
    protocol_copy = panel / "protocol.json"
    require(protocol_copy.is_file() and not protocol_copy.is_symlink(), "frozen protocol copy is missing")
    require(read_json(protocol_copy) == protocol, "frozen protocol copy differs")
    require(run.get("protocol_file_sha256") == sha256_file(protocol_copy), "protocol byte hash changed")
    amendment_copy = panel / "amendment.json"
    require(amendment_copy.is_file() and not amendment_copy.is_symlink(), "frozen amendment copy is missing")
    amendment = read_json(amendment_copy)
    validate_amendment(amendment, protocol)
    require(read_json(DEFAULT_AMENDMENT) == amendment, "panel amendment differs from the tracked amendment")
    require(run.get("amendment_sha256") == canonical_sha256(amendment), "run amendment hash mismatch")
    require(run.get("amendment_file_sha256") == sha256_file(amendment_copy), "amendment byte hash changed")
    require(run.get("failed_v1_custody") == validate_failed_v1_custody(amendment, protocol), "failed v1 custody summary changed")
    revision = run.get("source_revision")
    require(isinstance(revision, dict), "run lacks source revision")
    require(revision.get("algorithm_base_commit") == ALGORITHM_BASE_COMMIT, "run base commit changed")
    require(revision.get("algorithm_base_is_ancestor") is True, "run did not establish base ancestry")
    execution_commit = _hex_identifier(revision.get("commit"), "execution commit")
    ancestry = subprocess.run(
        ["git", "merge-base", "--is-ancestor", ALGORITHM_BASE_COMMIT, execution_commit],
        cwd=REPO, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    require(ancestry.returncode == 0, "recorded execution commit does not descend from the frozen base")
    recomputed_delta = git_tracked_delta(REPO, execution_commit)
    require(revision.get("tracked_delta") == recomputed_delta, "recorded tracked source delta differs from Git")
    require(run.get("evidence_class") in {"scientific_candidate", "operational_smoke"}, "invalid evidence class")
    if run.get("evidence_class") == "scientific_candidate":
        require(revision.get("dirty") is False and revision.get("porcelain") == [], "scientific run began dirty")
        require_expected_control_delta(recomputed_delta)
    host = run.get("host")
    require(isinstance(host, dict), "run lacks host identity")
    require(isinstance(host.get("platform"), str) and host["platform"], "run platform is missing")
    require(isinstance(host.get("uname"), list) and len(host["uname"]) == 6, "run uname is missing")
    require(isinstance(host.get("logical_cpus"), int) and host["logical_cpus"] > 0, "run CPU count is invalid")
    memory = host.get("physical_memory_bytes")
    require(memory is None or isinstance(memory, int) and memory > 0, "run memory size is invalid")
    for tool in ("python", "rustc", "cargo"):
        probe = host.get(tool)
        require(
            isinstance(probe, dict) and probe.get("returncode") == 0
            and isinstance(probe.get("command"), list) and isinstance(probe.get("stdout"), str)
            and probe["stdout"],
            f"run {tool} identity is invalid",
        )
    sources = run.get("sources")
    require(isinstance(sources, dict), "run lacks source identities")
    require(set(sources) == set(SOURCE_SHA256) | {"dependency_lock_snapshot"}, "run source inventory changed")
    for relative, expected_hash in SOURCE_SHA256.items():
        relative_path = Path(relative)
        _validate_mapped_identity(
            sources[relative], recorded_repo, relative_path, REPO / relative_path,
            f"frozen source {relative}", expected_hash,
        )
        require(sources[relative].get("expected_sha256") == expected_hash, f"missing expected source hash: {relative}")
    amendment_source = run.get("amendment_source")
    _validate_mapped_identity(
        amendment_source, recorded_repo, AMENDMENT_RELATIVE, DEFAULT_AMENDMENT,
        "Stage 18 lock amendment", EXPECTED_AMENDMENT_FILE_SHA256,
    )
    lock_relative = CORRECTED_LOCK_RELATIVE
    lock = REPO / lock_relative
    _validate_mapped_identity(
        sources["dependency_lock_snapshot"], recorded_repo, lock_relative, lock,
        "tracked dependency lock snapshot", LOCK_SHA256,
    )
    require(sources["dependency_lock_snapshot"].get("expected_sha256") == LOCK_SHA256, "missing expected lock hash")
    implementation = run.get("implementation")
    require(isinstance(implementation, dict) and set(implementation) == {"runner", "verifier", "meter"}, "run implementation inventory changed")
    _validate_mapped_identity(
        implementation["runner"], recorded_repo, RUNNER_RELATIVE, HERE / RUNNER_RELATIVE.name,
        "Stage 18 runner",
    )
    _validate_mapped_identity(
        implementation["verifier"], recorded_repo, VERIFIER_RELATIVE,
        HERE / "verify_stage18_degree23_panel.py",
        "Stage 18 verifier",
    )
    _validate_mapped_identity(
        implementation["meter"], recorded_repo, METER_RELATIVE, REPO / METER_RELATIVE,
        "process meter", SOURCE_SHA256["scripts/process_meter.py"],
    )
    lock_record = run.get("dependency_lock")
    require(
        isinstance(lock_record, dict)
        and set(lock_record) == {"source", "destination", "archive", "installed_by_runner"}
        and isinstance(lock_record.get("installed_by_runner"), bool),
        "run lacks exact installed-lock receipt",
    )
    _validate_mapped_identity(
        lock_record["source"], recorded_repo, lock_relative, lock,
        "installed-lock source", LOCK_SHA256,
    )
    _validate_historical_identity(
        lock_record["destination"], recorded_repo, Path("Cargo.lock"),
        "execution root Cargo.lock", LOCK_SHA256,
    )
    _validate_mapped_identity(
        lock_record["archive"], recorded_panel, LOCK_ARCHIVE_RELATIVE,
        panel / LOCK_ARCHIVE_RELATIVE, "archived Cargo.lock", LOCK_SHA256,
    )
    tools = run.get("tools")
    require(isinstance(tools, dict) and set(tools) == {"ic_binary", "discovery_binary"}, "run binary inventory changed")
    build = panel / "build"
    require(build.is_dir() and not build.is_symlink(), "build archive missing")
    _regular_inventory(build, {"process-invocation.json", "stdout.txt", "stderr.txt", "metrics.json", "receipt.json"})
    build_receipt = read_json(build / "receipt.json")
    command = build_receipt.get("command")
    require(
        isinstance(command, list) and len(command) == 8
        and Path(command[0]).name == "cargo"
        and command[1:] == [
            "build", "--release", "--locked", "--example", "koblitz_algebraic_e2e",
            "--example", "koblitz_public_factor_base_discovery",
        ],
        "locked release build command changed",
    )
    _canonical_absolute(command[0], "historical cargo executable")
    build_metrics = read_json(build / "metrics.json")
    _, build_success = validate_meter(build_metrics, command, BUILD_WATCHDOG_SECONDS)
    build_invocation = read_json(build / "process-invocation.json")
    require(build_invocation.get("schema") == "koblitz_degree23_process_invocation.v1", "wrong build invocation schema")
    python_path = run.get("host", {}).get("python", {}).get("command", [None])[0]
    _canonical_absolute(python_path, "historical Python launcher")
    recorded_build = recorded_panel / ".staging-build"
    _validate_process_invocation(
        build_invocation, command, recorded_repo, BUILD_WATCHDOG_SECONDS,
        implementation["meter"]["path"], recorded_build / "stdout.txt",
        recorded_build / "stderr.txt", recorded_build / "metrics.json", python_path,
    )
    require(build_receipt.get("dependency_lock_sha256") == LOCK_SHA256, "build lock hash changed")
    require(run.get("build_receipt_sha256") == sha256_file(build / "receipt.json"), "build receipt hash changed")
    if not build_success:
        require(allow_build_failure, "locked release build process failed")
        require(
            build_receipt.get("status") == "failed" and build_receipt.get("binaries") == {},
            "failed build receipt changed",
        )
        require(tools == {
            "ic_binary": {"path": str(recorded_repo / IC_BINARY_RELATIVE), "sha256": None},
            "discovery_binary": {"path": str(recorded_repo / DISCOVERY_BINARY_RELATIVE), "sha256": None},
        }, "failed build unexpectedly produced binary identities")
        require(run.get("status") == "inconclusive_build_failure", "failed build run status changed")
        require(not require_post, "failed build cannot have complete post-execution custody")
        return False
    require(build_receipt.get("status") == "verified", "locked release build was not verified")
    _validate_historical_identity(
        tools["ic_binary"], recorded_repo, IC_BINARY_RELATIVE, "historical IC binary",
    )
    _validate_historical_identity(
        tools["discovery_binary"], recorded_repo, DISCOVERY_BINARY_RELATIVE,
        "historical discovery binary",
    )
    require(build_receipt.get("binaries") == tools, "build/run binary identities differ")
    if require_post:
        post = run.get("post_execution")
        require(isinstance(post, dict), "complete run lacks post-execution custody")
        require(post.get("sources") == sources, "source identity changed after execution")
        require(post.get("binaries") == tools, "binary identity changed after execution")
        require(post.get("root_lock_sha256") == LOCK_SHA256, "root lock changed after execution")
    return True


def validate_outer_attempts(protocol: dict, panel: Path, run: dict,
                            require_success: bool) -> list[dict]:
    recorded_repo, recorded_panel = _execution_roots(run)
    root = panel / "outer-attempts"
    if not root.exists():
        require(not require_success, "complete panel lacks outer driver receipt")
        return []
    require(root.is_dir() and not root.is_symlink(), "outer-attempts is not a regular directory")
    records = []
    directories = sorted(root.iterdir())
    require(len(directories) <= 1, "corrected Stage 18 execution was retried")
    for expected_index, directory in enumerate(directories, start=1):
        require(
            directory.is_dir() and not directory.is_symlink()
            and directory.name == f"{expected_index:04d}",
            "invalid or non-sequential outer attempt directory",
        )
        recorded_directory = recorded_panel / "outer-attempts" / f".staging-{expected_index:04d}"
        expected = {"invocation.json", "stdout.txt", "stderr.txt", "metrics.json", "receipt.json"}
        _regular_inventory(directory, expected)
        invocation = read_json(directory / "invocation.json")
        require(invocation.get("schema") == "koblitz_degree23_outer_invocation.v1", "wrong outer invocation schema")
        require(invocation.get("attempt") == int(directory.name), "outer attempt index changed")
        command = invocation.get("command")
        require(isinstance(command, list) and len(command) >= 11, "invalid outer child command")
        require(command[:4] == [
            run["host"]["python"]["command"][0], run["implementation"]["runner"]["path"],
            "--inner", "--protocol",
        ], "outer child prefix changed")
        require(command[4] == str(recorded_repo / PROTOCOL_RELATIVE), "outer protocol path changed")
        require(command[5:11] == [
            "--amendment", str(recorded_repo / AMENDMENT_RELATIVE),
            "--output", str(recorded_panel), "--meter", run["implementation"]["meter"]["path"],
        ], "outer child paths changed")
        allowed_tail = []
        if run.get("evidence_class") == "operational_smoke":
            allowed_tail.append("--allow-dirty")
        require(command[11:] == allowed_tail, "outer child options changed")
        _validate_process_invocation(
            invocation, command, recorded_repo, 7200.0, run["implementation"]["meter"]["path"],
            recorded_directory / "stdout.txt", recorded_directory / "stderr.txt",
            recorded_directory / "metrics.json",
            run["host"]["python"]["command"][0],
        )
        meter = read_json(directory / "metrics.json")
        metrics, success = validate_meter(meter, command, 7200.0)
        receipt = read_json(directory / "receipt.json")
        require(receipt == {
            "schema": "koblitz_degree23_outer_receipt.v1",
            "attempt": int(directory.name), "returncode": meter["returncode"],
            "timed_out": meter["timed_out"],
            "orphan_group_terminated": meter["orphan_group_terminated"],
            "metrics": metrics, "command": command,
        }, "outer receipt changed")
        records.append({"attempt": int(directory.name), "metrics": metrics, "success": success})
    indexed = run.get("outer_attempts")
    require(isinstance(indexed, list), "run outer-attempt index is missing")
    require(len(indexed) == len(records), "run outer-attempt index is incomplete")
    for item, record in zip(indexed, records):
        receipt_path = root / f"{record['attempt']:04d}" / "receipt.json"
        require(item == {
            "attempt": record["attempt"],
            "receipt": str(recorded_panel / "outer-attempts" / f"{record['attempt']:04d}" / "receipt.json"),
            "receipt_sha256": sha256_file(receipt_path),
        }, "run outer-attempt receipt index changed")
    if require_success:
        require(records and records[-1]["success"], "final outer driver attempt did not complete")
    return records


def validate_final_inventory(panel: Path) -> dict:
    manifest = read_json(panel / "artifact-manifest.json")
    require(manifest.get("schema") == "koblitz_degree23_artifact_manifest.v1", "wrong artifact-manifest schema")
    exclusions = {"run.json", "artifact-manifest.json", "verification.json"}
    require(set(manifest.get("excluded_mutable_controls", [])) == exclusions, "manifest exclusion set changed")
    rows = manifest.get("files")
    require(isinstance(rows, list) and manifest.get("file_count") == len(rows), "manifest file count changed")
    listed = set()
    for row in rows:
        require(isinstance(row, dict) and set(row) == {"path", "bytes", "sha256"}, "invalid manifest row")
        relative = Path(row["path"])
        require(not relative.is_absolute() and ".." not in relative.parts, "unsafe manifest path")
        path = panel / relative
        require(path.is_file() and not path.is_symlink(), f"manifest path is not a regular file: {relative}")
        require(row["path"] not in listed, f"duplicate manifest path: {relative}")
        listed.add(row["path"])
        require(row["bytes"] == path.stat().st_size and row["sha256"] == sha256_file(path), f"manifest hash mismatch: {relative}")
    observed = set()
    for path in panel.rglob("*"):
        if path.is_dir():
            require(not path.is_symlink(), f"symlink directory in panel: {path}")
            continue
        require(path.is_file() and not path.is_symlink(), f"non-regular panel artifact: {path}")
        relative = str(path.relative_to(panel))
        require(not relative.endswith(".tmp"), f"temporary file retained: {relative}")
        observed.add(relative)
    require(observed == listed | exclusions, "final root inventory differs from manifest")
    verification = read_json(panel / "verification.json")
    summary = read_json(panel / "summary.json")
    run = read_json(panel / "run.json")
    require(verification.get("schema") == "koblitz_degree23_panel_verification.v1", "wrong verification schema")
    require(verification.get("status") == summary.get("status"), "verification/summary status mismatch")
    require(verification.get("protocol_sha256") == run.get("protocol_sha256"), "verification protocol hash mismatch")
    require(
        verification.get("amendment_sha256") == run.get("amendment_sha256")
        == summary.get("amendment_sha256"),
        "verification amendment hash mismatch",
    )
    require(verification.get("artifact_manifest_sha256") == sha256_file(panel / "artifact-manifest.json"), "verification manifest hash changed")
    require(verification.get("summary_sha256") == sha256_file(panel / "summary.json"), "verification summary hash changed")
    require(verification.get("exact_frozen_task_count") == 12, "verification task count changed")
    require(verification.get("outer_attempt_finalized") == len(run.get("outer_attempts", [])), "verification outer attempt changed")
    require(verification.get("claim_boundary") == summary.get("claim_boundary"), "verification claim boundary changed")
    require(run.get("status") == summary.get("status"), "run/summary status mismatch")
    require(run.get("summary_sha256") == verification.get("summary_sha256"), "run summary hash mismatch")
    require(run.get("artifact_manifest_sha256") == verification.get("artifact_manifest_sha256"), "run manifest hash mismatch")
    require(run.get("verification_sha256") == sha256_file(panel / "verification.json"), "run verification hash mismatch")
    return manifest


def validate_optional_controls(panel: Path, summary: dict, validate_controls: bool) -> None:
    if not validate_controls:
        return
    controls = [panel / name for name in ("summary.json", "artifact-manifest.json", "verification.json")]
    if any(path.exists() for path in controls):
        require(
            all(path.is_file() and not path.is_symlink() for path in controls),
            "final panel controls are incomplete or non-regular",
        )
        require(read_json(panel / "summary.json") == summary, "archived summary differs from recomputation")
        validate_final_inventory(panel)


def _charged_totals(receipts: list[dict]) -> dict:
    metrics = [item["metrics"] for item in receipts]
    return {
        "processes": len(metrics),
        "wall_seconds_sum": sum(item["wall_seconds"] for item in metrics),
        "user_seconds_sum": sum(item["user_seconds"] for item in metrics),
        "system_seconds_sum": sum(item["system_seconds"] for item in metrics),
        "total_core_seconds_sum": sum(item["total_core_seconds"] for item in metrics),
        "maximum_peak_rss_bytes": max((item["peak_rss_bytes"] for item in metrics), default=None),
        "wall_seconds_distribution": _distribution([item["wall_seconds"] for item in metrics]),
        "total_core_seconds_distribution": _distribution([item["total_core_seconds"] for item in metrics]),
        "peak_rss_bytes_distribution": _distribution([item["peak_rss_bytes"] for item in metrics]),
    }


def _attempted_totals(receipts: list[dict]) -> dict:
    charged = [item for item in receipts if isinstance(item.get("metrics"), dict)]
    totals = _charged_totals(charged)
    totals["status_counts"] = {
        status: sum(item.get("status") == status for item in receipts)
        for status in sorted({item.get("status") for item in receipts})
    }
    totals["receipts_without_valid_metrics"] = len(receipts) - len(charged)
    return totals


def _internal_distributions(receipts: list[dict]) -> dict:
    keys = sorted({key for receipt in receipts for key in receipt.get("result", {}).get("timing_ns", {})})
    return {
        key: _distribution([
            receipt["result"]["timing_ns"][key] / 1e9
            for receipt in receipts if key in receipt.get("result", {}).get("timing_ns", {})
        ])
        for key in keys
    }


def summarize(protocol: dict, panel: Path, allow_incomplete: bool = False,
              validate_controls: bool = True) -> dict:
    validate_protocol(protocol)
    run = read_json(panel / "run.json")
    build_success = validate_run_custody(
        protocol, panel, run, require_post=not allow_incomplete,
        allow_build_failure=allow_incomplete,
    )
    _, recorded_panel = _execution_roots(run)
    outer_attempts = validate_outer_attempts(
        protocol, panel, run, require_success=not allow_incomplete
    )
    if not build_success:
        require(allow_incomplete, "failed build cannot form a complete panel")
        require(run.get("tasks") == {}, "failed build panel contains scientific task records")
        build_meter = read_json(panel / "build/metrics.json")
        build_metrics = build_meter["metrics"]
        outer_core = sum(item["metrics"]["total_core_seconds"] for item in outer_attempts)
        outer_wall = sum(item["metrics"]["wall_seconds"] for item in outer_attempts)
        summary = {
            "schema": SUMMARY_SCHEMA,
            "status": "inconclusive_build_failure",
            "task_panel_complete": False,
            "expected_tasks": 12,
            "attempted_tasks": 0,
            "verified_tasks": 0,
            "status_counts": {},
            "rows": [
                {"index": index, "secret": secret, "seed": str(seed), "status": "not_started"}
                for index, secret, seed in EXPECTED_RUNS
            ],
            "comparison": "inconclusive_no_scientific_tasks",
            "accounting": {
                "scientific_processes": 0,
                "build_operational_metrics": build_metrics,
                "build_returncode": build_meter["returncode"],
                "algorithm_ratios": None,
            },
            "outer_driver_receipts": {
                "attempts": len(outer_attempts),
                "all_completed_attempts_core_seconds": outer_core,
                "all_completed_attempts_wall_seconds": outer_wall,
                "maximum_peak_rss_bytes": max(
                    (item["metrics"]["peak_rss_bytes"] for item in outer_attempts),
                    default=None,
                ),
                "nested_build_cost_must_not_be_added_to_outer_envelope": True,
            },
            "failed_v1_custody": run["failed_v1_custody"],
            "amendment_sha256": run["amendment_sha256"],
            "claim_boundary": read_json(panel / "amendment.json")["claim_boundary"],
        }
        validate_optional_controls(panel, summary, validate_controls)
        return summary
    tools = run.get("tools", {})
    discovery_binary = tools.get("discovery_binary", {}).get("path")
    ic_binary = tools.get("ic_binary", {}).get("path")
    require(isinstance(discovery_binary, str) and isinstance(ic_binary, str), "run lacks binary identities")
    tasks = task_plan(protocol, discovery_binary, ic_binary)
    receipts = {}
    for task in tasks:
        leaf = panel / "tasks" / task["id"]
        if not leaf.exists():
            if allow_incomplete:
                continue
            raise VerificationError(f"missing frozen task: {task['id']}")
        try:
            receipts[task["id"]] = verify_task_leaf(protocol, panel, task)
        except VerificationError as error:
            archived = read_json(leaf / "receipt.json") if (leaf / "receipt.json").is_file() else None
            if not isinstance(archived, dict) or archived.get("status") != "invalid_artifact" or archived.get("error") != str(error):
                raise
            receipts[task["id"]] = archived
    indexed_tasks = run.get("tasks")
    require(isinstance(indexed_tasks, dict), "run task index is missing")
    if not allow_incomplete:
        require(set(indexed_tasks) == {task["id"] for task in tasks}, "run task index is incomplete")
    for task_id, receipt in receipts.items():
        index = indexed_tasks.get(task_id)
        receipt_path = panel / "tasks" / task_id / "receipt.json"
        require(isinstance(index, dict), f"run task index missing: {task_id}")
        require(index.get("status") == receipt.get("status"), f"run task status changed: {task_id}")
        require(
            index.get("receipt") == str(recorded_panel / "tasks" / task_id / "receipt.json"),
            f"run receipt path changed: {task_id}",
        )
        require(index.get("receipt_sha256") == sha256_file(receipt_path), f"run receipt hash changed: {task_id}")
    verified = {key: value for key, value in receipts.items() if value.get("status") == "verified"}
    rows = []
    targets = []
    for index, secret, seed in EXPECTED_RUNS:
        ic = verified.get(f"row-{index:02d}-ic")
        rho = verified.get(f"row-{index:02d}-rho-auto")
        row = {"index": index, "secret": secret, "seed": str(seed)}
        if ic and rho:
            require(ic["result"]["target"] == rho["result"]["target"], f"row {index} targets differ")
            targets.append(canonical_sha256(ic["result"]["target"]))
            ic_metrics, rho_metrics = ic["metrics"], rho["metrics"]
            require(rho_metrics["total_core_seconds"] > 0, f"row {index} rho core time is zero")
            ic_report = ic["result"]["report"]
            rho_result = rho["result"]
            row.update({
                "status": "verified_pair", "target": ic["result"]["target"],
                "index_calculus": {
                    "relations": ic_report["relations"], "trials": ic_report["trials"],
                    "relation_batches": ic_report["relation_batches"],
                    "sat_calls": ic_report["sat_calls"], "sat_models": ic_report["sat_models"],
                    "sat_refutations": ic_report["sat_refutations"],
                    "sat_unknowns": ic_report["sat_unknowns"],
                    "sat_conflicts": ic_report["sat_conflicts"],
                    "direct_relations_skipped": ic_report["direct_relations_skipped"],
                    "linear_solve_attempts": ic_report["linear_solve_attempts"],
                    "timing_ns": ic["result"]["timing_ns"], **ic_metrics,
                },
                "automorphism_rho": {
                    "iterations": rho_result["iterations"], "restarts": rho_result["restarts"],
                    "walk_step_additions": rho_result["walk_step_additions"],
                    "jump_table_additions": rho_result["jump_table_additions"],
                    "initial_state_additions": rho_result["initial_state_additions"],
                    "reported_group_additions": rho_result["reported_group_additions"],
                    "rho_setup_scalar_multiplications": rho_result["rho_setup_scalar_multiplications"],
                    "target_construction_scalar_multiplications": rho_result["target_construction_scalar_multiplications"],
                    "timing_ns": rho_result["timing_ns"],
                    **rho_metrics,
                },
                "ic_over_rho_core_ratio": ic_metrics["total_core_seconds"] / rho_metrics["total_core_seconds"],
            })
        else:
            row.update({
                "status": "inconclusive",
                "ic_status": receipts.get(f"row-{index:02d}-ic", {}).get("status", "missing"),
                "rho_status": receipts.get(f"row-{index:02d}-rho-auto", {}).get("status", "missing"),
            })
        rows.append(row)
    if len(targets) == 5:
        require(len(set(targets)) == 5, "frozen public targets are not unique")
    discovery_receipts = [item for item in (verified.get("discovery-a0"), verified.get("discovery-a1")) if item]
    paired = [row for row in rows if row["status"] == "verified_pair"]
    ic_receipts = [verified[f"row-{index:02d}-ic"] for index, _, _ in EXPECTED_RUNS if f"row-{index:02d}-ic" in verified]
    rho_receipts = [verified[f"row-{index:02d}-rho-auto"] for index, _, _ in EXPECTED_RUNS if f"row-{index:02d}-rho-auto" in verified]
    attempted_discovery = [receipts[task["id"]] for task in tasks if task["kind"] == "discovery" and task["id"] in receipts]
    attempted_ic = [receipts[task["id"]] for task in tasks if task["kind"] == "ic" and task["id"] in receipts]
    attempted_rho = [receipts[task["id"]] for task in tasks if task["kind"] == "rho-auto" and task["id"] in receipts]
    discovery_totals = _charged_totals(discovery_receipts)
    ic_totals = _charged_totals(ic_receipts)
    rho_totals = _charged_totals(rho_receipts)
    discovery_totals["internal_timing_seconds_distributions"] = _internal_distributions(discovery_receipts)
    ic_totals["internal_timing_seconds_distributions"] = _internal_distributions(ic_receipts)
    rho_totals["internal_timing_seconds_distributions"] = _internal_distributions(rho_receipts)
    discovery_core = discovery_totals["total_core_seconds_sum"]
    ic_core = ic_totals["total_core_seconds_sum"]
    rho_core = rho_totals["total_core_seconds_sum"]
    ratios = [row["ic_over_rho_core_ratio"] for row in paired]
    all_tasks_verified = len(verified) == 12
    task_complete = all_tasks_verified and len(paired) == 5 and len(targets) == 5
    complete = task_complete and run.get("evidence_class") == "scientific_candidate" and not allow_incomplete
    if task_complete:
        if all(value > 1 for value in ratios):
            comparison = "rho_faster_all_pairs"
        elif all(value < 1 for value in ratios):
            comparison = "ic_faster_all_pairs"
        else:
            comparison = "mixed"
    else:
        comparison = "inconclusive_partial"
    build_metrics = read_json(panel / "build" / "metrics.json")["metrics"]
    outer_wall = sum(item["metrics"]["wall_seconds"] for item in outer_attempts)
    outer_core = sum(item["metrics"]["total_core_seconds"] for item in outer_attempts)
    inner_wall = build_metrics["wall_seconds"] + sum(
        item["metrics"]["wall_seconds"] for item in receipts.values() if isinstance(item.get("metrics"), dict)
    )
    inner_core = build_metrics["total_core_seconds"] + sum(
        item["metrics"]["total_core_seconds"] for item in receipts.values() if isinstance(item.get("metrics"), dict)
    )
    summary = {
        "schema": SUMMARY_SCHEMA,
        "amendment_sha256": run["amendment_sha256"],
        "failed_v1_custody": run["failed_v1_custody"],
        "status": (
            "complete_verified_panel" if complete
            else "provisional_complete_pending_final_controls" if task_complete
            else "inconclusive_partial"
        ),
        "task_panel_complete": task_complete,
        "expected_tasks": 12, "attempted_tasks": len(receipts), "verified_tasks": len(verified),
        "status_counts": {
            status: sum(item.get("status") == status for item in receipts.values())
            for status in sorted({item.get("status") for item in receipts.values()})
        },
        "rows": rows,
        "comparison": comparison,
        "accounting": {
            "public_discovery_one_time_setup": discovery_totals,
            "index_calculus_online": ic_totals,
            "automorphism_rho": rho_totals,
            "online_ic_over_rho_core_ratio": ic_core / rho_core if len(paired) == 5 and rho_core > 0 else None,
            "strict_setup_charged_ic_over_rho_core_ratio": (discovery_core + ic_core) / rho_core if task_complete and rho_core > 0 else None,
            "per_pair_ratio_distribution": _distribution(ratios),
            "aggregate_ic_relations": sum(row["index_calculus"]["relations"] for row in paired),
            "aggregate_ic_trials": sum(row["index_calculus"]["trials"] for row in paired),
            "aggregate_ic_sat_conflicts": sum(row["index_calculus"]["sat_conflicts"] for row in paired),
            "aggregate_rho_iterations": sum(row["automorphism_rho"]["iterations"] for row in paired),
            "aggregate_rho_reported_group_additions": sum(row["automorphism_rho"]["reported_group_additions"] for row in paired),
            "aggregate_rho_setup_scalar_multiplications": sum(row["automorphism_rho"]["rho_setup_scalar_multiplications"] for row in paired),
            "attempted_process_costs_including_failures": {
                "public_discovery": _attempted_totals(attempted_discovery),
                "index_calculus": _attempted_totals(attempted_ic),
                "automorphism_rho": _attempted_totals(attempted_rho),
                "all_scientific_tasks": _attempted_totals(list(receipts.values())),
            },
            "per_target_factor_base_materialisation_is_already_in_ic_process": True,
            "public_discovery_is_charged_once_and_not_double_counted": True,
            "build_and_outer_driver_are_operational_envelopes": True,
        },
        "outer_driver_receipts": {
            "attempts": len(outer_attempts),
            "all_completed_attempts_core_seconds": outer_core,
            "all_completed_attempts_wall_seconds": outer_wall,
            "maximum_peak_rss_bytes": max((item["metrics"]["peak_rss_bytes"] for item in outer_attempts), default=None),
            "inner_build_plus_scientific_process_core_seconds": inner_core,
            "inner_build_plus_scientific_process_wall_seconds": inner_wall,
            "outer_minus_inner_core_seconds": outer_core - inner_core if outer_attempts else None,
            "outer_minus_sequential_inner_wall_seconds": outer_wall - inner_wall if outer_attempts else None,
            "build_operational_metrics": build_metrics,
            "excluded_from_algorithm_ratios": True,
        },
        "claim_boundary": protocol["claim_boundary"],
    }
    if not allow_incomplete:
        require(len(receipts) == 12, "panel archive is incomplete")
    validate_optional_controls(panel, summary, validate_controls)
    return summary


def self_test() -> dict:
    protocol = read_json(DEFAULT_PROTOCOL)
    validate_protocol(protocol)
    amendment = read_json(DEFAULT_AMENDMENT)
    validate_amendment(amendment, protocol)
    failure = validate_failed_v1_custody(amendment, protocol)
    if failure["scientific_tasks_started"] != 0:
        raise AssertionError("failed v1 custody admitted a scientific task")
    checks = 3
    amendment_mutation = copy.deepcopy(amendment)
    amendment_mutation["supersession"]["fields"][1]["replacement"] = V1_LOCK_SHA256
    try:
        validate_amendment(amendment_mutation, protocol)
    except VerificationError:
        checks += 1
    else:
        raise AssertionError("amendment lock-hash mutation was accepted")
    stage14 = HERE / "stage-14-orbit-admission-20260909"
    for a in (0, 1):
        validate_discovery_result(read_json(stage14 / f"a{a}" / "discovery.json"), a)
        checks += 1
    mutated = read_json(stage14 / "a0" / "discovery.json")
    mutated["forbidden_inputs"]["relation_yield_used"] = True
    try:
        validate_discovery_result(mutated, 0)
    except VerificationError:
        checks += 1
    else:
        raise AssertionError("forbidden discovery input mutation was accepted")
    bound_mutation = copy.deepcopy(protocol)
    bound_mutation["implementation_binding"]["source_sha256"]["Cargo.toml"] = "0" * 64
    try:
        validate_protocol(bound_mutation)
    except VerificationError:
        checks += 1
    else:
        raise AssertionError("protocol source-hash mutation was accepted")
    claim_mutation = copy.deepcopy(protocol)
    claim_mutation["claim_boundary"] = "This establishes SOTA."
    try:
        validate_protocol(claim_mutation)
    except VerificationError:
        checks += 1
    else:
        raise AssertionError("non-source frozen-protocol mutation was accepted")
    expected_delta = {
        "schema": "koblitz_stage18_tracked_delta.v1",
        "base_commit": ALGORITHM_BASE_COMMIT,
        "head_commit": "2" * 40,
        "entries": [dict(entry, object_id="1" * 40) for entry in EXPECTED_CONTROL_DELTA],
        "matches_expected_control_only_delta": True,
    }
    require_expected_control_delta(expected_delta)
    checks += 1
    extra_delta = copy.deepcopy(expected_delta)
    extra_delta["entries"].append({
        "status": "M", "path": "src/binary_ecc/f2m.rs", "mode": "100644",
        "type": "blob", "object_id": "3" * 40,
    })
    extra_delta["matches_expected_control_only_delta"] = False
    try:
        require_expected_control_delta(extra_delta)
    except VerificationError:
        checks += 1
    else:
        raise AssertionError("extra transitive source delta was accepted")
    stage10 = HERE / "stage-10-degree23-public-20260909"
    ic = read_json(stage10 / "ic.json")
    rho = read_json(stage10 / "rho-auto-v2.json")
    # Stage 10 predates the two separately reported public setup timings.
    # Supply zero-valued fields only to exercise the unchanged result schema.
    ic["timing_ns"]["projected_orbit_construction"] = 0
    ic["timing_ns"]["cofactor_admission"] = 0
    old_row = {"secret": 101, "seed": 20261001, "public_target": ic["target"]}
    validate_ic_result(ic, old_row)
    validate_rho_result(rho, old_row)
    checks += 2
    ic_oracle = copy.deepcopy(ic)
    ic_oracle["factor_base_predicate"]["uses_discrete_log_labels"] = True
    try:
        validate_ic_result(ic_oracle, old_row)
    except VerificationError:
        checks += 1
    else:
        raise AssertionError("IC discrete-log oracle mutation was accepted")
    ic_counts = copy.deepcopy(ic)
    ic_counts["report"]["sat_unknowns"] += 1
    try:
        validate_ic_result(ic_counts, old_row)
    except VerificationError:
        checks += 1
    else:
        raise AssertionError("IC counter imbalance was accepted")
    ic_batches = copy.deepcopy(ic)
    ic_batches["report"]["relation_batches"] = 0
    try:
        validate_ic_result(ic_batches, old_row)
    except VerificationError:
        checks += 1
    else:
        raise AssertionError("IC one-target batch imbalance was accepted")
    rho_counts = copy.deepcopy(rho)
    rho_counts["jump_table_additions"] += 1
    try:
        validate_rho_result(rho_counts, old_row)
    except VerificationError:
        checks += 1
    else:
        raise AssertionError("rho operation-ledger mutation was accepted")
    rho_walk = copy.deepcopy(rho)
    rho_walk["walk_step_additions"] -= 1
    rho_walk["reported_group_additions"] -= 1
    try:
        validate_rho_result(rho_walk, old_row)
    except VerificationError:
        checks += 1
    else:
        raise AssertionError("rho walk/iteration imbalance was accepted")
    target_mutation = copy.deepcopy(ic)
    target_mutation["target"] = {"x": [], "y": {}}
    try:
        validate_ic_result(target_mutation, old_row)
    except VerificationError:
        checks += 1
    else:
        raise AssertionError("noncanonical public target was accepted")
    meter = read_json(stage10 / "ic.metrics.json")
    validate_meter(meter, meter["command"], 300.0)
    checks += 1
    meter_mutation = copy.deepcopy(meter)
    meter_mutation["command"][-1] = "divisor:0,1"
    try:
        validate_meter(meter_mutation, meter["command"], 300.0)
    except VerificationError:
        checks += 1
    else:
        raise AssertionError("meter command mutation was accepted")
    with tempfile.TemporaryDirectory(prefix="stage18-verifier-self-test-") as temporary:
        root = Path(temporary)
        duplicate = root / "duplicate.json"
        duplicate.write_text('{"x":1,"x":2}\n')
        try:
            read_json(duplicate)
        except VerificationError:
            checks += 1
        else:
            raise AssertionError("duplicate JSON key was accepted")
        recorded_repo = Path("/recorded/stage18/repository")
        recorded_panel = Path("/recorded/stage18/panel")
        relocated_panel = root / "relocated-panel"
        relocated_lock = relocated_panel / LOCK_ARCHIVE_RELATIVE
        relocated_lock.parent.mkdir(parents=True)
        relocated_lock.write_bytes((REPO / CORRECTED_LOCK_RELATIVE).read_bytes())
        lock_identity = {
            "path": str(recorded_panel / LOCK_ARCHIVE_RELATIVE),
            "bytes": relocated_lock.stat().st_size,
            "sha256": sha256_file(relocated_lock),
        }
        _validate_mapped_identity(
            lock_identity, recorded_panel, LOCK_ARCHIVE_RELATIVE, relocated_lock,
            "self-test relocated lock", LOCK_SHA256,
        )
        historical_binary = {
            "path": str(recorded_repo / IC_BINARY_RELATIVE),
            "bytes": 1,
            "sha256": "1" * 64,
        }
        _validate_historical_identity(
            historical_binary, recorded_repo, IC_BINARY_RELATIVE,
            "self-test historical binary",
        )
        checks += 2
        escaped = copy.deepcopy(lock_identity)
        escaped["path"] = str(recorded_panel) + "/../outside/Cargo.lock"
        try:
            _validate_mapped_identity(
                escaped, recorded_panel, LOCK_ARCHIVE_RELATIVE, relocated_lock,
                "self-test escaped lock", LOCK_SHA256,
            )
        except VerificationError:
            checks += 1
        else:
            raise AssertionError("recorded path escaping its execution root was accepted")
        relocated_lock.write_bytes(b"changed")
        try:
            _validate_mapped_identity(
                lock_identity, recorded_panel, LOCK_ARCHIVE_RELATIVE, relocated_lock,
                "self-test changed relocated lock", LOCK_SHA256,
            )
        except VerificationError:
            checks += 1
        else:
            raise AssertionError("changed relocated lock was accepted")
    return {
        "self_test": "pass", "checks": checks, "frozen_tasks": 12,
        "protocol_sha256": canonical_sha256(protocol),
        "amendment_sha256": canonical_sha256(amendment),
        "failed_v1_scientific_tasks": 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--panel", type=Path, default=DEFAULT_PANEL)
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    try:
        if args.self_test:
            print(json.dumps(self_test(), indent=2, sort_keys=True))
            return
        protocol = read_json(args.protocol)
        print(json.dumps(summarize(protocol, args.panel.resolve(), args.allow_incomplete), indent=2, sort_keys=True))
    except VerificationError as error:
        parser.exit(1, f"verification failed: {error}\n")


if __name__ == "__main__":
    main()
