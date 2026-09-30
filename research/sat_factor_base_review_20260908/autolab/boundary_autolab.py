#!/usr/bin/env python3
"""Agent-facing autolab control plane for the IC boundary ledger.

Reads docs/ic/boundary_targets.json (schema_version 2), fail-closed validates
measurement reports, and launches public-synthetic Koblitz producer runs
against the next admitted beat targets.

This runner is local to the crypto repository. It does not impersonate any
external Autoresearcher dispatcher or Coordinator.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import platform
import resource
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
RESEARCH_DIR = HERE.parent
REPO = HERE.parents[2]
PROTOCOL_PATH = HERE / "protocol.json"
RUNS_DIR = HERE / "runs"
CURRENT_PATH = RUNS_DIR / "current.json"
LOCK_PATH = RUNS_DIR / "autolab.lock"
TASK_ID = "TASK-IC-BOUNDARY-AUTOLAB-20260910"

# Abstract measurement_schema keys -> accepted concrete aliases in reports.
FIELD_ALIASES: dict[str, tuple[str, ...]] = {
    "n_or_bits": ("n_or_bits", "n", "bits"),
    "factor_base_size_F": ("factor_base_size_F", "factor_base_size", "F", "|F|"),
    "orbit_count_K": ("orbit_count_K", "orbit_columns", "K", "orbit_columns_K"),
    "dimension_l_or_dim": ("dimension_l_or_dim", "dimension", "l", "ell", "dim"),
    "construction_method": ("construction_method",),
    "materialized": ("materialized",),
    "construction_wall_ms": ("construction_wall_ms",),
    "retained_bytes": ("retained_bytes",),
    "m_summands": ("m_summands", "m", "summands"),
    "unknowns": ("unknowns",),
    "system_degree": ("system_degree",),
    "eq_var_ratio": ("eq_var_ratio",),
    "ffd_or_degree_of_regularity": (
        "ffd_or_degree_of_regularity",
        "ffd",
        "degree_of_regularity",
        "DoR",
        "dor",
    ),
    "oracle_class": ("oracle_class",),
    "median_ms_per_target": ("median_ms_per_target",),
    "largest_solvable": ("largest_solvable",),
    "base_id_or_hash": ("base_id_or_hash", "base_id", "base_hash"),
    "eta_or_coverage_policy": ("eta_or_coverage_policy", "eta", "coverage_policy"),
    "pr_decomposition_or_hit_rate_with_ci": (
        "pr_decomposition_or_hit_rate_with_ci",
        "pr_decomposition",
        "hit_rate",
        "hit_rate_with_ci",
    ),
    "trials_per_relation": ("trials_per_relation",),
    "target_mix": ("target_mix",),
    "orbit_columns_K": ("orbit_columns_K", "orbit_columns", "K"),
    "relations_collected": ("relations_collected",),
    "relations_needed": ("relations_needed",),
    "surplus": ("surplus",),
    "matrix_dims": ("matrix_dims",),
    "sparse_or_dense": ("sparse_or_dense",),
    "rank_accumulation": ("rank_accumulation",),
    "la_wall_ms_or_la_charged_ms": (
        "la_wall_ms_or_la_charged_ms",
        "la_wall_ms",
        "la_charged_ms",
    ),
    "recovered_d_verified": ("recovered_d_verified", "recovered_d", "d_verified"),
    "stage_timers": ("stage_timers",),
    "claim_boundary": ("claim_boundary",),
    "timing_class": ("timing_class",),
    "target_count": ("target_count",),
    "ic_target_hash": ("ic_target_hash",),
    "rho_target_hash": ("rho_target_hash",),
    "ic_online_wall_ms": ("ic_online_wall_ms",),
    "rho_online_wall_ms": ("rho_online_wall_ms",),
    "online_speedup": ("online_speedup",),
    "online_interval": ("online_interval",),
    "same_resource_envelope": ("same_resource_envelope",),
    "ic_scalar_verified": ("ic_scalar_verified",),
    "rho_scalar_verified": ("rho_scalar_verified",),
    "independent_validation": ("independent_validation",),
    "ic_replay_certificate_sha256": ("ic_replay_certificate_sha256",),
    "rho_replay_certificate_sha256": ("rho_replay_certificate_sha256",),
    "ic_resource_envelope": ("ic_resource_envelope",),
    "rho_resource_envelope": ("rho_resource_envelope",),
    "rho_policy": ("rho_policy",),
    "ic_cost": ("ic_cost",),
    "rho_cost": ("rho_cost",),
    "automorphism_discount": ("automorphism_discount",),
    "all_stages_charged_same_series": ("all_stages_charged_same_series",),
    "verdict": ("verdict",),
    "independent_replay_pointer": ("independent_replay_pointer",),
    "fixture_hash": ("fixture_hash",),
    "executable_or_source_hash": (
        "executable_or_source_hash",
        "executable_hash",
        "source_hash",
    ),
    "host_id": ("host_id",),
    "resource_caps": ("resource_caps",),
    "seeds": ("seeds", "seed"),
    "claim_boundary_non_claims": (
        "claim_boundary_non_claims",
        "non_claims",
        "claim_boundary",
    ),
}

ONLINE_REQUIRED_STAGES: dict[str, tuple[str, ...]] = {
    "ic_included_stages": (
        "target_query",
        "target_PDP",
        "target_relation_check",
        "target_descent",
        "target_recovery_check",
    ),
    "rho_included_stages": ("walk", "collision", "recovery_check"),
}


class AutolabError(RuntimeError):
    """Operator-facing failure that must not be treated as a ledger beat."""


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256(path: Path | str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text())


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AutolabError(message)


def load_protocol() -> dict[str, Any]:
    protocol = read_json(PROTOCOL_PATH)
    require(protocol.get("task_id") == TASK_ID, "protocol task_id mismatch")
    return protocol


def ledger_path(protocol: dict[str, Any]) -> Path:
    return REPO / protocol["ledger"]["path"]


def load_ledger(protocol: dict[str, Any]) -> dict[str, Any]:
    path = ledger_path(protocol)
    require(path.is_file(), f"boundary ledger missing: {path}")
    ledger = read_json(path)
    required = protocol["ledger"]["required_schema_version"]
    require(
        ledger.get("schema_version") == required,
        f"ledger schema_version must be {required}, got {ledger.get('schema_version')}",
    )
    require(
        isinstance(ledger.get("measurement_schema"), dict),
        "ledger missing measurement_schema (fail closed)",
    )
    require(
        ledger["measurement_schema"].get("fail_closed") is True,
        "measurement_schema.fail_closed must be true",
    )
    return ledger


def field_present(report: dict[str, Any], abstract_key: str) -> bool:
    aliases = FIELD_ALIASES.get(abstract_key, (abstract_key,))
    for alias in aliases:
        if alias not in report:
            continue
        value = report[alias]
        if value is None:
            continue
        if isinstance(value, str) and not value.strip():
            continue
        return True
    return False


def missing_fields(report: dict[str, Any], keys: list[str]) -> list[str]:
    return [key for key in keys if not field_present(report, key)]


def validate_claim(
    report: dict[str, Any],
    *,
    stage: str,
    ledger: dict[str, Any],
) -> dict[str, Any]:
    schema = ledger["measurement_schema"]
    require(stage in schema, f"unknown stage for measurement schema: {stage}")
    stage_schema = schema[stage]
    required = list(stage_schema.get("required", []))
    global_required = list(schema.get("global_provenance_required", []))
    missing_stage = missing_fields(report, required)
    missing_global = missing_fields(report, global_required)
    validation_errors: list[str] = []

    if stage == "vs_rho":
        for key in ("candidate_id", "workload_id", "run_id"):
            value = report.get(key)
            if not isinstance(value, str) or not value.strip():
                validation_errors.append(f"{key} must be a nonempty string")
        candidate_id = report.get("candidate_id")
        workload_id = report.get("workload_id")
        run_id = report.get("run_id")
        candidate_manifest_hash = report.get("candidate_manifest_sha256")
        workload_manifest_hash = report.get("workload_manifest_sha256")
        if not isinstance(candidate_manifest_hash, str) or len(candidate_manifest_hash) != 64:
            validation_errors.append("candidate_manifest_sha256 must be a full SHA-256 hex digest")
        elif any(ch not in "0123456789abcdef" for ch in candidate_manifest_hash):
            validation_errors.append("candidate_manifest_sha256 must be lowercase hex")
        if not isinstance(workload_manifest_hash, str) or len(workload_manifest_hash) != 64:
            validation_errors.append("workload_manifest_sha256 must be a full SHA-256 hex digest")
        elif any(ch not in "0123456789abcdef" for ch in workload_manifest_hash):
            validation_errors.append("workload_manifest_sha256 must be lowercase hex")
        if isinstance(candidate_id, str) and isinstance(candidate_manifest_hash, str):
            import re
            match = re.search(r"h([0-9a-f]{12,64})$", candidate_id)
            if not candidate_id.startswith("IC1") or not match:
                validation_errors.append("candidate_id must use the IC1 identity format with a digest suffix")
            elif not candidate_manifest_hash.startswith(match.group(1)):
                validation_errors.append("candidate_id digest must match candidate_manifest_sha256")
        if isinstance(workload_id, str) and isinstance(workload_manifest_hash, str):
            if len(workload_id) != 12 or any(ch not in "0123456789abcdef" for ch in workload_id):
                validation_errors.append("workload_id must be 12 lowercase hex digits")
            elif not workload_manifest_hash.startswith(workload_id):
                validation_errors.append("workload_id must match workload_manifest_sha256")
        if all(isinstance(value, str) and value for value in (candidate_id, workload_id, run_id)):
            import re
            if not re.fullmatch(re.escape(candidate_id) + r"W" + re.escape(workload_id) + r"R[1-9][0-9]*", run_id):
                validation_errors.append("run_id must be <candidate_id>W<workload_id>R<run-number>")
        if type(report.get("target_count")) is not int or report["target_count"] != 1:
            validation_errors.append("target_count must equal 1")
        target_hashes = (report.get("ic_target_hash"), report.get("rho_target_hash"))
        if not all(isinstance(value, str) and value.strip() for value in target_hashes):
            validation_errors.append("IC and rho target hashes must be nonempty strings")
        elif target_hashes[0] != target_hashes[1]:
            validation_errors.append("IC and rho target hashes must match")
        if report.get("timing_class") not in stage_schema.get("timing_class_enum", []):
            validation_errors.append("timing_class must be single_target_online_wall")
        if report.get("same_resource_envelope") is not True:
            validation_errors.append("same_resource_envelope must be true")
        for key in ("ic_scalar_verified", "rho_scalar_verified"):
            if report.get(key) is not True:
                validation_errors.append(f"{key} must be true")
        if report.get("independent_validation") is not True:
            validation_errors.append("independent_validation must be true")
        for key in ("ic_replay_certificate_sha256", "rho_replay_certificate_sha256"):
            digest = report.get(key)
            if not isinstance(digest, str) or len(digest) != 64 or any(
                ch not in "0123456789abcdef" for ch in digest
            ):
                validation_errors.append(f"{key} must be a full lowercase SHA-256 digest")
        ic_resources = report.get("ic_resource_envelope")
        rho_resources = report.get("rho_resource_envelope")
        if not isinstance(ic_resources, dict) or not ic_resources:
            validation_errors.append("ic_resource_envelope must be a nonempty object")
        if not isinstance(rho_resources, dict) or not rho_resources:
            validation_errors.append("rho_resource_envelope must be a nonempty object")
        if isinstance(ic_resources, dict) and isinstance(rho_resources, dict):
            if ic_resources != rho_resources:
                validation_errors.append("IC and rho resource envelopes must match exactly")

        ic_ms = positive_cost(report.get("ic_online_wall_ms"))
        rho_ms = positive_cost(report.get("rho_online_wall_ms"))
        speedup = positive_cost(report.get("online_speedup"))
        if ic_ms is None or rho_ms is None or speedup is None:
            validation_errors.append("online times and speedup must be positive finite numbers")
        elif not math.isclose(speedup, rho_ms / ic_ms, rel_tol=1e-9, abs_tol=1e-12):
            validation_errors.append("online_speedup must equal rho_online_wall_ms / ic_online_wall_ms")

        phase_costs = report.get("ic_online_phase_ms")
        phase_fields = stage_schema.get("ic_online_phase_fields", [])
        if not isinstance(phase_costs, dict):
            validation_errors.append("ic_online_phase_ms must be an object")
        else:
            phase_values = []
            for key in phase_fields:
                value = phase_costs.get(key)
                if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
                    validation_errors.append(f"ic_online_phase_ms.{key} must be a finite nonnegative number")
                else:
                    phase_values.append(float(value))
            if len(phase_values) == len(phase_fields) and ic_ms is not None:
                if not math.isclose(math.fsum(phase_values), ic_ms, rel_tol=1e-6, abs_tol=1e-3):
                    validation_errors.append("IC exclusive phase costs must sum to ic_online_wall_ms")

        interval = report.get("online_interval")
        interval_fields = (
            stage_schema.get("online_interval_event_fields", [])
            + stage_schema.get("online_interval_stage_fields", [])
        )
        if not isinstance(interval, dict):
            validation_errors.append("online_interval must be an object")
        else:
            for key in interval_fields:
                value = interval.get(key)
                valid = (
                    isinstance(value, str) and bool(value.strip())
                    if key.endswith("_event")
                    else isinstance(value, list) and bool(value)
                    and all(isinstance(item, str) and item.strip() for item in value)
                )
                if not valid:
                    validation_errors.append(f"online_interval.{key} is missing or invalid")
            for included_key, required_stages in ONLINE_REQUIRED_STAGES.items():
                included = interval.get(included_key)
                if isinstance(included, list):
                    missing = sorted(set(required_stages) - set(included))
                    if missing:
                        validation_errors.append(
                            f"online_interval.{included_key} is missing required stages: {', '.join(missing)}"
                        )

        policy = report.get("rho_policy")
        if not isinstance(policy, dict):
            validation_errors.append("rho_policy must be an object")
        else:
            for key in stage_schema.get("rho_policy_integer_fields", []):
                value = policy.get(key)
                minimum = stage_schema.get("rho_policy_minimums", {}).get(key)
                if type(value) is not int or (minimum is not None and value < minimum):
                    validation_errors.append(f"rho_policy.{key} must be an integer >= {minimum}")
            for key in stage_schema.get("rho_policy_string_fields", []):
                value = policy.get(key)
                if not isinstance(value, str) or not value.strip():
                    validation_errors.append(f"rho_policy.{key} must be a nonempty string")

    ok = not missing_stage and not missing_global and not validation_errors
    return {
        "schema_version": ledger.get("schema_version"),
        "stage": stage,
        "fail_closed": True,
        "status": "PASS" if ok else "FAIL",
        "missing_stage_fields": missing_stage,
        "missing_global_provenance": missing_global,
        "validation_errors": validation_errors,
        "required_stage_fields": required,
        "required_global_provenance": global_required,
    }


class RunnerLock:
    """Exclusive advisory lock using an atomic create + live-pid check."""

    def __enter__(self) -> "RunnerLock":
        RUNS_DIR.mkdir(parents=True, exist_ok=True)
        while True:
            try:
                fd = os.open(LOCK_PATH, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            except FileExistsError:
                try:
                    record = read_json(LOCK_PATH)
                except Exception:
                    record = {}
                pid = int(record.get("pid") or 0)
                if pid and _pid_alive(pid):
                    raise AutolabError(f"autolab lock is held by pid {pid}")
                LOCK_PATH.unlink(missing_ok=True)
                continue
            payload = {"pid": os.getpid(), "created_at": now(), "task_id": TASK_ID}
            os.write(fd, json.dumps(payload, indent=2, sort_keys=True).encode() + b"\n")
            os.close(fd)
            return self

    def __exit__(self, *_: Any) -> None:
        try:
            record = read_json(LOCK_PATH)
        except Exception:
            record = {}
        if record.get("pid") == os.getpid():
            LOCK_PATH.unlink(missing_ok=True)


def _pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def make_run_id(beat_id: str) -> str:
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    digest = hashlib.sha256(f"{TASK_ID}|{beat_id}|{stamp}".encode()).hexdigest()[:10]
    return f"{stamp}-{digest}"


def current_run_id() -> str:
    require(CURRENT_PATH.is_file(), "no current autolab run")
    run_id = read_json(CURRENT_PATH).get("run_id")
    require(bool(run_id), "current.json missing run_id")
    return str(run_id)


def resolve_run(run_id: str | None) -> Path:
    resolved = run_id or current_run_id()
    path = RUNS_DIR / resolved
    require(path.is_dir(), f"run directory missing: {path}")
    return path


def host_record() -> dict[str, Any]:
    return {
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "machine": platform.machine(),
        "node": platform.node(),
        "cpu_count": os.cpu_count(),
    }


def preflight(protocol: dict[str, Any], ledger: dict[str, Any]) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []

    def add(name: str, ok: bool, detail: str) -> None:
        checks.append({"name": name, "ok": ok, "detail": detail})

    add(
        "ledger_schema_v2",
        ledger.get("schema_version") == 2 and "measurement_schema" in ledger,
        f"schema_version={ledger.get('schema_version')}",
    )
    cargo = shutil.which("cargo")
    add("cargo", cargo is not None, cargo or "cargo not on PATH")
    rustc = shutil.which("rustc")
    add("rustc", rustc is not None, rustc or "rustc not on PATH")
    cms = Path(os.environ.get("KIC_AUTOLAB_CMS", "/opt/homebrew/bin/cryptominisat5"))
    add(
        "cryptominisat5_optional",
        cms.is_file() or shutil.which("cryptominisat5") is not None,
        str(cms if cms.is_file() else shutil.which("cryptominisat5") or "absent"),
    )
    for key, producer in protocol["producers"].items():
        source = REPO / producer["source"]
        add(f"producer_source_{key}", source.is_file(), str(source))
    git_head = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )
    add("git_head", git_head.returncode == 0, git_head.stdout.strip())
    dirty = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=False,
    )
    add(
        "git_status_recorded",
        dirty.returncode == 0,
        f"{len(dirty.stdout.splitlines())} dirty paths",
    )
    ok = all(
        check["ok"]
        for check in checks
        if check["name"] != "cryptominisat5_optional"
    )
    return {
        "schema_version": "1.0",
        "task_id": TASK_ID,
        "kind": "local_autolab_preflight",
        "ok": ok,
        "host": host_record(),
        "ledger_sha256": sha256(ledger_path(protocol)),
        "protocol_sha256": sha256(PROTOCOL_PATH),
        "agent_priorities": ledger.get("agent_priorities", []),
        "checks": checks,
        "claim_boundary": (
            "preflight only; not a fixed-arm cost, relation, rank, memory, rho, "
            "or crossover result"
        ),
        "created_at": now(),
    }


def plan(protocol: dict[str, Any], ledger: dict[str, Any]) -> dict[str, Any]:
    contract = protocol.get("primary_speedup_contract", {})
    contract_ready = contract.get("status") == "ready" and contract.get("target_count") == 1
    beats = []
    for beat_id, beat in sorted(
        protocol["beats"].items(),
        key=lambda item: (item[1].get("priority", 99), item[0]),
    ):
        eligible = contract_ready and bool(beat.get("primary_speedup_eligible")) and (
            beat.get("timing_class_goal") == "single_target_online_wall"
        )
        beats.append(
            {
                "beat_id": beat_id,
                "priority": beat.get("priority"),
                "label": beat.get("label"),
                "regime": beat.get("regime"),
                "stage": beat.get("stage"),
                "n": beat.get("n"),
                "timing_class_goal": beat.get("timing_class_goal"),
                "primary_speedup_eligible": eligible,
                "primary_speedup_blocker": None if eligible else (
                    beat.get("primary_speedup_blocker")
                    or protocol.get("primary_speedup_contract", {}).get("blocker")
                    or "beat is diagnostic-only"
                ),
            }
        )
    return {
        "schema_version": "1.0",
        "task_id": TASK_ID,
        "ledger": str(protocol["ledger"]["path"]),
        "ledger_schema_version": ledger.get("schema_version"),
        "primary_speedup_contract": contract,
        "agent_priorities": ledger.get("agent_priorities", []),
        "beats": beats,
        "how_to_beat": ledger.get("how_to_beat", []),
        "commands": {
            "preflight": "python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py preflight",
            "claim_check": "python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py claim-check --report PATH --stage vs_rho",
        },
    }


def build_producers() -> dict[str, str]:
    examples = ["koblitz_rank_fixture", "koblitz_rho_fixture"]
    for example in examples:
        source = REPO / "examples" / f"{example}.rs"
        require(source.is_file(), f"missing producer source: {source}")
    command = [
        "cargo",
        "build",
        "--release",
        "--example",
        "koblitz_rank_fixture",
        "--example",
        "koblitz_rho_fixture",
    ]
    completed = subprocess.run(command, cwd=REPO, capture_output=True, text=True)
    require(
        completed.returncode == 0,
        "cargo build failed:\n" + completed.stderr[-4000:],
    )
    binaries = {
        "direct": str((REPO / "target/release/examples/koblitz_rank_fixture").resolve()),
        "rho": str((REPO / "target/release/examples/koblitz_rho_fixture").resolve()),
    }
    for path in binaries.values():
        require(Path(path).is_file(), f"built binary missing: {path}")
    return binaries


RESOURCE_RECEIPT_FIELDS = (
    "exit_code",
    "whole_process_wall_ms",
    "whole_process_wall_samples_ms",
    "cold_start_wall_ms",
    "timed_repeats",
    "stdout_stable",
    "children_cpu_user_ms",
    "children_cpu_system_ms",
    "command",
    "seed",
)

DEFAULT_TIMED_REPEATS = 3
# Past this, further repeats cost more than the spread they resolve, so stop
# after the first timed execution.
REPEAT_BUDGET_MS = 20_000.0


def _run_once(
    command: list[str], *, env: dict[str, str], cwd: Path
) -> tuple[subprocess.CompletedProcess[str], float, float, float]:
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
    )
    elapsed_ms = (time.perf_counter() - started) * 1000.0
    usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    # Children CPU deltas (seconds -> ms). Wall is the outer process wait.
    cpu_user_ms = (usage_after.ru_utime - usage_before.ru_utime) * 1000.0
    cpu_system_ms = (usage_after.ru_stime - usage_before.ru_stime) * 1000.0
    return completed, elapsed_ms, cpu_user_ms, cpu_system_ms


def _timing_free(text: str) -> str:
    """Producer stdout with its own timing fields dropped.

    Every row carries measured durations, so raw stdout never repeats exactly.
    What should repeat is the computation: base hash, rank, status, solution.
    """
    rows = []
    for row in parse_json_lines(text):
        rows.append(
            {
                k: v
                for k, v in row.items()
                if not k.endswith("_ms") and k != "timing_breakdown_ms"
            }
        )
    return json.dumps(rows, sort_keys=True)


def run_timed(
    command: list[str],
    *,
    env: dict[str, str],
    cwd: Path,
    repeats: int = DEFAULT_TIMED_REPEATS,
) -> dict[str, Any]:
    """Time a producer without charging it for first-execution cost.

    The first execution of a freshly linked binary pays page-in and, on macOS,
    signature validation: measured here at 183-484 ms against 7 ms warm, which
    is 30x the whole n=13 comparison and larger than either arm at n=37. Timing
    a single run therefore adds a roughly constant per-process term to both
    arms, which drags every ratio toward parity and flatters whichever arm is
    slower.

    The warmup runs the binary with no arguments, which both producers reject
    immediately. That still pays the whole load cost -- after it the real run
    drops from 190 ms to 7.7 ms -- so a rung that takes a quarter of an hour is
    not run twice to save 200 ms.

    The producers take their seed in argv and are deterministic, so repeats are
    the same computation; `stdout_stable` records whether that held, comparing
    the rows with their own timing fields dropped.
    """
    require(repeats >= 1, f"repeats must be >= 1, got {repeats}")
    _, cold_wall_ms, _, _ = _run_once([command[0]], env=env, cwd=cwd)

    walls: list[float] = []
    users: list[float] = []
    systems: list[float] = []
    outputs: list[str] = []
    completed: subprocess.CompletedProcess[str] | None = None
    for _ in range(repeats):
        completed, wall_ms, user_ms, system_ms = _run_once(command, env=env, cwd=cwd)
        walls.append(wall_ms)
        users.append(user_ms)
        systems.append(system_ms)
        outputs.append(_timing_free(completed.stdout))
        if completed.returncode != 0 or wall_ms > REPEAT_BUDGET_MS:
            break

    return {
        "command": command,
        "exit_code": completed.returncode,
        "whole_process_wall_ms": statistics.median(walls),
        "whole_process_wall_samples_ms": walls,
        "cold_start_wall_ms": cold_wall_ms,
        "timed_repeats": len(walls),
        "stdout_stable": len(set(outputs)) <= 1,
        "children_cpu_user_ms": statistics.median(users),
        "children_cpu_system_ms": statistics.median(systems),
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def parse_json_lines(text: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            value = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(value, dict):
            rows.append(value)
    return rows


def seed_for(beat_id: str, arm: str, repetition: int) -> int:
    material = f"{TASK_ID}|{beat_id}|arm={arm}|repetition={repetition}".encode()
    return int.from_bytes(hashlib.sha256(material).digest()[:8], "big")


def positive_cost(value):
    return float(value) if type(value) in (int, float) and math.isfinite(value) and value > 0 else None


def extract_ic_cost(rows: list[dict[str, Any]], timing_class: str) -> float | None:
    if not rows:
        return None
    row = rows[-1]
    if timing_class == "whole_process_wall":
        return None  # filled from outer wait4
    # Online-only and projection-matched slices cannot stand in for a full job.
    for key in ("full_algorithm_charged_total_ms",):
        if key in row:
            return positive_cost(row[key])
        timing = row.get("timing_breakdown_ms")
        if isinstance(timing, dict) and key in timing:
                return positive_cost(timing[key])
    return None


def extract_rho_cost(rows: list[dict[str, Any]]) -> float | None:
    """Return a cost only for exactly one explicit target result row."""
    targets = [row for row in rows if row.get("kind") == "rho_public_fixture"]
    if len(targets) != 1:
        return None
    row = targets[0]
    if "total_ms" in row:
        return positive_cost(row["total_ms"])
    timing = row.get("timing_breakdown_ms")
    if isinstance(timing, dict) and "total_ms" in timing:
        return positive_cost(timing["total_ms"])
    return None


def comparison_integrity(direct_rows, rho_rows, expected_count=None):
    """Bind a comparison to actual complete public targets, never source bytes.

    This checks producer consistency only. Independent group/rank replay is
    still required for a scientific result.
    """
    try:
        require(expected_count == 1, "primary comparison requires exactly one target")
    except AutolabError as exc:
        return {"status": "INVALID_COMPARISON", "reason": str(exc), "fixture_hash": None,
                "independent_validation": False}

    def corpus(rows, kind, verified_key):
        selected = [r for r in rows if r.get('kind') == kind]
        require(bool(selected), 'no complete fixture records')
        headers = [r for r in rows if r.get('kind') == 'point_defined_factor_base']
        if kind == 'relation_rank_summary':
            require(len(headers) == 1, 'missing or ambiguous curve header')
        entries = []
        for row in selected:
            metadata = headers[0] if kind == 'relation_rank_summary' else row
            require(row.get(verified_key) is True, 'unverified fixture')
            require(row.get('recovered_fixture_scalar') is not None, 'missing recovered scalar')
            require(int(row['recovered_fixture_scalar']) == int(row['published_fixture_scalar']), 'wrong scalar')
            entries.append({'index': int(row['fixture_index']), 'n': int(row['n']), 'a': int(row['a']),
                'subgroup_order': int(metadata['subgroup_order']),
                'field_modulus_low_terms': sorted(int(t) for t in metadata['field_modulus_low_terms']),
                'generator': [int(x) for x in row['generator_point_key']],
                'target': [int(x) for x in row['published_q_point_key']]})
        entries.sort(key=lambda r: r['index'])
        count = len(entries) if expected_count is None else expected_count
        require([r['index'] for r in entries] == list(range(count)), 'missing or duplicate fixture')
        return entries
    try:
        direct = corpus(direct_rows, 'relation_rank_summary', 'linear_solution_verified')
        rho = corpus(rho_rows, 'rho_public_fixture', 'verified')
        require(direct == rho, 'IC and rho solved different public targets')
        fingerprint = hashlib.sha256(json.dumps(direct, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        return {'status': 'MATCHED', 'fixture_hash': fingerprint, 'fixtures': len(direct),
                'independent_validation': False}
    except (AutolabError, KeyError, TypeError, ValueError) as exc:
        return {'status': 'INVALID_COMPARISON', 'reason': str(exc), 'fixture_hash': None,
                'independent_validation': False}


def automorphism_discount(n: int) -> dict[str, Any]:
    return {
        "formula": "sqrt(2*n)",
        "A": 2 * n,
        "description": "signed Frobenius automorphism discount family used by Koblitz rho control",
        "n": n,
    }


def draft_vs_rho_claim(
    *,
    beat: dict[str, Any],
    beat_id: str,
    run: Path,
    direct_obs: dict[str, Any],
    rho_obs: dict[str, Any],
    binaries: dict[str, str],
) -> dict[str, Any]:
    timing_class = beat["timing_class_goal"]
    direct_rows = parse_json_lines(direct_obs["stdout"])
    rho_rows = parse_json_lines(rho_obs["stdout"])
    producers_ok = direct_obs["exit_code"] == 0 and rho_obs["exit_code"] == 0
    integrity = comparison_integrity(direct_rows, rho_rows, direct_obs.get('fixtures'))
    ic_cost = (
        float(direct_obs["whole_process_wall_ms"])
        if timing_class == "whole_process_wall"
        else extract_ic_cost(direct_rows, timing_class)
    )
    rho_cost = (
        float(rho_obs["whole_process_wall_ms"])
        if timing_class == "whole_process_wall"
        else extract_rho_cost(rho_rows)
    )
    claim = {
        "schema_version": 2,
        "task_id": TASK_ID,
        "beat_id": beat_id,
        "regime": beat["regime"],
        "stage": "vs_rho",
        "n": beat["n"],
        "n_or_bits": beat["n"],
        "timing_class": timing_class,
        "candidate_id": None,
        "candidate_manifest_sha256": None,
        "workload_id": None,
        "workload_manifest_sha256": None,
        "run_id": None,
        "autolab_run_id": run.name,
        "target_count": integrity.get("fixtures"),
        "ic_target_hash": integrity.get("fixture_hash") if integrity["status"] == "MATCHED" else None,
        "rho_target_hash": integrity.get("fixture_hash") if integrity["status"] == "MATCHED" else None,
        "ic_online_wall_ms": None,
        "rho_online_wall_ms": None,
        "ic_online_phase_ms": None,
        "online_speedup": None,
        "online_interval": None,
        "same_resource_envelope": None,
        "ic_scalar_verified": integrity["status"] == "MATCHED" and producers_ok,
        "rho_scalar_verified": integrity["status"] == "MATCHED" and producers_ok,
        "independent_validation": False,
        "ic_replay_certificate_sha256": None,
        "rho_replay_certificate_sha256": None,
        "ic_resource_envelope": None,
        "rho_resource_envelope": None,
        "rho_policy": None,
        "ic_cost": ic_cost,
        "rho_cost": rho_cost,
        "automorphism_discount": automorphism_discount(int(beat["n"])),
        "all_stages_charged_same_series": bool(producers_ok and integrity['status'] == 'MATCHED'
            and timing_class == 'whole_process_wall'
            and all(isinstance(x, (int, float)) and math.isfinite(x) and x > 0 for x in (ic_cost, rho_cost))),
        "comparison_integrity": integrity,
        "verdict": (
            ("DRAFT_PENDING_INDEPENDENT_VALIDATION" if integrity['status'] == 'MATCHED' else 'INVALID_COMPARISON')
            if producers_ok
            else "PRODUCER_FAILURE"
        ),
        "claim_boundary": (
            "Legacy single-fixture producer diagnostic only. Current producers "
            "report whole-process or operation-counted costs, not target-online "
            "wall intervals; this row is not a primary single-target speedup."
        ),
        "claim_boundary_non_claims": [
            "not key recovery",
            "not asymptotic sub-sqrt",
            "not imported/external points",
            "not ledger promotion until independent validation",
        ],
        "independent_replay_pointer": str(
            (run / "artifacts/claim_draft.json").relative_to(REPO)
        ),
        "fixture_hash": integrity['fixture_hash'],
        "executable_or_source_hash": {
            "direct": sha256(binaries["direct"]),
            "rho": sha256(binaries["rho"]),
            "rank_fixture_source": sha256(REPO / "examples/koblitz_rank_fixture.rs"),
            "rho_fixture_source": sha256(REPO / "examples/koblitz_rho_fixture.rs"),
        },
        "host_id": host_record(),
        "resource_caps": {"common_cap_bytes": beat.get("resource_cap_bytes")},
        "seeds": {
            "direct": direct_obs.get("seed"),
            "rho": rho_obs.get("seed"),
        },
        "producer_exit_codes": {
            "direct": direct_obs["exit_code"],
            "rho": rho_obs["exit_code"],
        },
        "whole_process_wall_ms": {
            "direct": direct_obs["whole_process_wall_ms"],
            "rho": rho_obs["whole_process_wall_ms"],
        },
        "whole_process_wall_method": {
            "description": (
                "median of timed executions after one discarded warmup; "
                "first-execution cost is reported separately and not charged"
            ),
            "cold_start_wall_ms": {
                "direct": direct_obs.get("cold_start_wall_ms"),
                "rho": rho_obs.get("cold_start_wall_ms"),
            },
            "samples_ms": {
                "direct": direct_obs.get("whole_process_wall_samples_ms"),
                "rho": rho_obs.get("whole_process_wall_samples_ms"),
            },
            "timed_repeats": {
                "direct": direct_obs.get("timed_repeats"),
                "rho": rho_obs.get("timed_repeats"),
            },
            "stdout_stable": {
                "direct": direct_obs.get("stdout_stable"),
                "rho": rho_obs.get("stdout_stable"),
            },
        },
        "direct_rows": len(direct_rows),
        "rho_rows": len(rho_rows),
    }
    return claim


def launch(arguments: argparse.Namespace) -> dict[str, Any]:
    protocol = load_protocol()
    ledger = load_ledger(protocol)
    beat_id = arguments.beat
    require(beat_id in protocol["beats"], f"unknown beat id: {beat_id}")
    beat = protocol["beats"][beat_id]
    fixtures = arguments.fixtures
    if fixtures is None:
        fixtures = int(beat.get("fixtures", beat.get("fixtures_default", 1)))
    require(fixtures == 1, "run exactly one target per workload; batch fixtures are disabled")
    repeats = getattr(arguments, "repeats", None) or DEFAULT_TIMED_REPEATS
    require(repeats > 0, "repeats must be positive")

    with RunnerLock():
        preflight_receipt = preflight(protocol, ledger)
        require(preflight_receipt["ok"], "preflight failed; see artifacts after launch dir create")
        run_id = arguments.run_id or make_run_id(beat_id)
        run = RUNS_DIR / run_id
        require(not run.exists(), f"run already exists: {run}")
        for name in ("artifacts", "inputs", "logs", "receipts"):
            (run / name).mkdir(parents=True, exist_ok=True)
        write_json(CURRENT_PATH, {"run_id": run_id, "beat_id": beat_id, "updated_at": now()})
        write_json(run / "artifacts/preflight.json", preflight_receipt)
        write_json(run / "inputs/protocol.json", protocol)
        write_json(run / "inputs/boundary_targets.json", ledger)
        write_json(
            run / "inputs/ledger_pin.json",
            {
                "path": protocol["ledger"]["path"],
                "sha256": sha256(ledger_path(protocol)),
                "schema_version": ledger["schema_version"],
            },
        )

        state: dict[str, Any] = {
            "schema_version": "1.0",
            "task_id": TASK_ID,
            "run_id": run_id,
            "beat_id": beat_id,
            "status": "ACTIVE",
            "phase": "build",
            "fixtures": fixtures,
            "created_at": now(),
            "updated_at": now(),
        }
        write_json(run / "state.json", state)

        if arguments.prepare_only:
            state.update(status="PREPARED", phase="prepared", updated_at=now())
            write_json(run / "state.json", state)
            commands = producer_commands(beat, beat_id, fixtures, binaries=None)
            write_json(run / "artifacts/prepared_commands.json", commands)
            return state

        binaries = build_producers()
        write_json(
            run / "artifacts/binaries.json",
            {key: {"path": path, "sha256": sha256(path)} for key, path in binaries.items()},
        )
        state.update(phase="measurement", updated_at=now())
        write_json(run / "state.json", state)

        env = os.environ.copy()
        # Default to incremental crosscheck for small rungs; beats may override
        # (n=53 dense recompute-after-every-relation is hour-class).
        env.setdefault("KIC_INCREMENTAL_RANK_CROSSCHECK", "1")
        for key, value in (beat.get("env") or {}).items():
            env[str(key)] = str(value)
        direct_seed = seed_for(beat_id, "direct", 0)
        rho_seed = seed_for(beat_id, "rho", 0)
        direct_cmd = [
            binaries["direct"],
            str(beat["n"]),
            str(beat["a"]),
            str(beat["eta"][0]),
            str(beat["eta"][1]),
            str(direct_seed),
            beat["direct"]["pair_mode"],
            beat["direct"]["target_mode"],
            beat["direct"]["query_mode"],
            str(fixtures),
        ]
        rho_cmd = [
            binaries["rho"],
            str(beat["n"]),
            str(beat["a"]),
            beat["rho"]["quotient_mode"],
            str(fixtures),
            beat["rho"]["backend"],
            str(rho_seed),
        ]
        write_json(
            run / "artifacts/commands.json",
            {"direct": direct_cmd, "rho": rho_cmd, "fixtures": fixtures},
        )

        direct_obs = run_timed(direct_cmd, env=env, cwd=REPO, repeats=repeats)
        direct_obs["seed"] = direct_seed
        direct_obs["fixtures"] = fixtures
        (run / "logs/direct.stdout.jsonl").write_text(direct_obs["stdout"])
        (run / "logs/direct.stderr.txt").write_text(direct_obs["stderr"])
        write_json(
            run / "receipts/direct.resource.json",
            {k: direct_obs[k] for k in RESOURCE_RECEIPT_FIELDS},
        )

        rho_obs = run_timed(rho_cmd, env=env, cwd=REPO, repeats=repeats)
        rho_obs["seed"] = rho_seed
        (run / "logs/rho.stdout.jsonl").write_text(rho_obs["stdout"])
        (run / "logs/rho.stderr.txt").write_text(rho_obs["stderr"])
        write_json(
            run / "receipts/rho.resource.json",
            {k: rho_obs[k] for k in RESOURCE_RECEIPT_FIELDS},
        )

        claim = draft_vs_rho_claim(
            beat=beat,
            beat_id=beat_id,
            run=run,
            direct_obs=direct_obs,
            rho_obs=rho_obs,
            binaries=binaries,
        )
        write_json(run / "artifacts/claim_draft.json", claim)
        validation = validate_claim(claim, stage="vs_rho", ledger=ledger)
        write_json(run / "artifacts/claim_check.json", validation)

        producers_ok = direct_obs["exit_code"] == 0 and rho_obs["exit_code"] == 0
        status = "PENDING_INDEPENDENT_VALIDATION" if producers_ok else "PRODUCER_FAILURE"
        if validation["status"] != "PASS":
            status = "SCHEMA_INCOMPLETE"
        if producers_ok and claim['comparison_integrity']['status'] != 'MATCHED':
            status = 'INVALID_COMPARISON'
        elif producers_ok and not claim['all_stages_charged_same_series']:
            status = 'ACCOUNTING_INCOMPLETE'
        state.update(
            status=status,
            phase="analysis",
            updated_at=now(),
            direct_exit_code=direct_obs["exit_code"],
            rho_exit_code=rho_obs["exit_code"],
            claim_check=validation["status"],
        )
        write_json(run / "state.json", state)
        write_json(
            run / "artifacts/candidate.json",
            {
                "schema_version": "1.0",
                "task_id": TASK_ID,
                "run_id": run_id,
                "beat_id": beat_id,
                "status": status,
                "claim_draft_sha256": sha256(run / "artifacts/claim_draft.json"),
                "claim_check": validation,
                "ledger_sha256": sha256(ledger_path(protocol)),
                "created_at": now(),
                "note": (
                    "Draft only. Promote the ledger only after independent validation "
                    "and a claim-check PASS with every required measurement field."
                ),
            },
        )
        files = {
            str(path.relative_to(run)): sha256(path)
            for path in sorted(run.rglob("*"))
            if path.is_file()
        }
        write_json(
            run / "artifacts/review_manifest.json",
            {"schema_version": "1.0", "task_id": TASK_ID, "files": files},
        )
        return state


def producer_commands(
    beat: dict[str, Any],
    beat_id: str,
    fixtures: int,
    binaries: dict[str, str] | None,
) -> dict[str, Any]:
    direct_bin = (
        binaries["direct"]
        if binaries
        else "target/release/examples/koblitz_rank_fixture"
    )
    rho_bin = (
        binaries["rho"] if binaries else "target/release/examples/koblitz_rho_fixture"
    )
    direct_seed = seed_for(beat_id, "direct", 0)
    rho_seed = seed_for(beat_id, "rho", 0)
    return {
        "build": (
            "cargo build --release --example koblitz_rank_fixture "
            "--example koblitz_rho_fixture"
        ),
        "direct": " ".join(
            [
                direct_bin,
                str(beat["n"]),
                str(beat["a"]),
                str(beat["eta"][0]),
                str(beat["eta"][1]),
                str(direct_seed),
                beat["direct"]["pair_mode"],
                beat["direct"]["target_mode"],
                beat["direct"]["query_mode"],
                str(fixtures),
            ]
        ),
        "rho": " ".join(
            [
                rho_bin,
                str(beat["n"]),
                str(beat["a"]),
                beat["rho"]["quotient_mode"],
                str(fixtures),
                beat["rho"]["backend"],
                str(rho_seed),
            ]
        ),
    }


def status(arguments: argparse.Namespace) -> None:
    run = resolve_run(arguments.run_id)
    print((run / "state.json").read_text(), end="")


def verify(arguments: argparse.Namespace) -> dict[str, Any]:
    run = resolve_run(arguments.run_id)
    manifest_path = run / "artifacts/review_manifest.json"
    require(manifest_path.is_file(), f"manifest missing in {run}")
    manifest = read_json(manifest_path)
    mismatches = [
        relative
        for relative, digest in manifest["files"].items()
        if not (run / relative).is_file() or sha256(run / relative) != digest
    ]
    result = {
        "schema_version": "1.0",
        "task_id": TASK_ID,
        "run_id": run.name,
        "status": "PASS" if not mismatches else "FAIL",
        "files": len(manifest["files"]),
        "mismatches": mismatches,
    }
    write_json(run / "artifacts/verification.json", result)
    return result


def claim_check(arguments: argparse.Namespace) -> dict[str, Any]:
    protocol = load_protocol()
    ledger = load_ledger(protocol)
    report = read_json(Path(arguments.report))
    stage = arguments.stage or report.get("stage")
    require(bool(stage), "stage required via --stage or report.stage")
    result = validate_claim(report, stage=str(stage), ledger=ledger)
    if arguments.out:
        write_json(Path(arguments.out), result)
    return result


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    sub = result.add_subparsers(dest="command", required=True)

    sub.add_parser("plan", help="Show ledger priorities and beat launch commands")
    sub.add_parser("preflight", help="Check ledger, deps, and producer sources")

    launch_parser = sub.add_parser("launch", help="Create a run and execute producers")
    launch_parser.add_argument(
        "--beat",
        required=True,
        help="Beat id from protocol.json (e.g. koblitz.vs_rho.n37_wall)",
    )
    launch_parser.add_argument("--run-id")
    launch_parser.add_argument(
        "--fixtures",
        type=int,
        help="Compatibility option; the only accepted value is 1 target per workload",
    )
    launch_parser.add_argument(
        "--repeats",
        type=int,
        help=(
            "Timed executions per arm after a discarded warmup exec "
            f"(default {DEFAULT_TIMED_REPEATS}; stops early once a run "
            f"exceeds {int(REPEAT_BUDGET_MS / 1000)}s)"
        ),
    )
    launch_parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Create run scaffolding and print commands without executing producers",
    )

    for name in ("status", "verify"):
        command = sub.add_parser(name)
        command.add_argument("--run-id")

    claim = sub.add_parser(
        "claim-check",
        help="Fail-closed validate a JSON report against measurement_schema",
    )
    claim.add_argument("--report", required=True)
    claim.add_argument("--stage", help="Stage id (default: report.stage)")
    claim.add_argument("--out", help="Optional path to write the check receipt")
    return result


def main() -> int:
    arguments = parser().parse_args()
    try:
        if arguments.command == "plan":
            protocol = load_protocol()
            ledger = load_ledger(protocol)
            print(json.dumps(plan(protocol, ledger), indent=2, sort_keys=True))
            return 0
        if arguments.command == "preflight":
            protocol = load_protocol()
            ledger = load_ledger(protocol)
            receipt = preflight(protocol, ledger)
            print(json.dumps(receipt, indent=2, sort_keys=True))
            return 0 if receipt["ok"] else 2
        if arguments.command == "launch":
            state = launch(arguments)
            print(json.dumps(state, indent=2, sort_keys=True))
            return 0 if state.get("status") != "PRODUCER_FAILURE" else 3
        if arguments.command == "status":
            status(arguments)
            return 0
        if arguments.command == "verify":
            print(json.dumps(verify(arguments), indent=2, sort_keys=True))
            return 0
        if arguments.command == "claim-check":
            result = claim_check(arguments)
            print(json.dumps(result, indent=2, sort_keys=True))
            return 0 if result["status"] == "PASS" else 4
        raise AutolabError(f"unknown command {arguments.command}")
    except AutolabError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
