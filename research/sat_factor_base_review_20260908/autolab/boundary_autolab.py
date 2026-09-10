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
import os
import platform
import resource
import shutil
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
    ok = not missing_stage and not missing_global
    return {
        "schema_version": ledger.get("schema_version"),
        "stage": stage,
        "fail_closed": True,
        "status": "PASS" if ok else "FAIL",
        "missing_stage_fields": missing_stage,
        "missing_global_provenance": missing_global,
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
    beats = []
    for beat_id, beat in sorted(
        protocol["beats"].items(),
        key=lambda item: (item[1].get("priority", 99), item[0]),
    ):
        beats.append(
            {
                "beat_id": beat_id,
                "priority": beat.get("priority"),
                "label": beat.get("label"),
                "regime": beat.get("regime"),
                "stage": beat.get("stage"),
                "n": beat.get("n"),
                "timing_class_goal": beat.get("timing_class_goal"),
            }
        )
    return {
        "schema_version": "1.0",
        "task_id": TASK_ID,
        "ledger": str(protocol["ledger"]["path"]),
        "ledger_schema_version": ledger.get("schema_version"),
        "agent_priorities": ledger.get("agent_priorities", []),
        "beats": beats,
        "how_to_beat": ledger.get("how_to_beat", []),
        "commands": {
            "preflight": "python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py preflight",
            "smoke": "python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py launch --beat smoke.koblitz.vs_rho.n13",
            "n37_wall": "python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py launch --beat koblitz.vs_rho.n37_wall",
            "n41_charged": "python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py launch --beat koblitz.vs_rho.n41_charged",
            "n37_full": "python3 research/sat_factor_base_review_20260908/autolab/boundary_autolab.py launch --beat koblitz.vs_rho.n37_wall --fixtures 1024",
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


def run_timed(command: list[str], *, env: dict[str, str], cwd: Path) -> dict[str, Any]:
    started = time.perf_counter()
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
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
    return {
        "command": command,
        "exit_code": completed.returncode,
        "whole_process_wall_ms": elapsed_ms,
        "children_cpu_user_ms": cpu_user_ms,
        "children_cpu_system_ms": cpu_system_ms,
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


def extract_ic_cost(rows: list[dict[str, Any]], timing_class: str) -> float | None:
    if not rows:
        return None
    row = rows[-1]
    if timing_class == "whole_process_wall":
        return None  # filled from outer wait4
    for key in (
        "full_algorithm_charged_total_ms",
        "charged_total_ms",
        "projection_matched_charged_total_ms",
        "online_charged_ms",
    ):
        if key in row:
            return float(row[key])
        timing = row.get("timing_breakdown_ms")
        if isinstance(timing, dict) and key in timing:
            return float(timing[key])
    if "charged_total_ms" in row:
        return float(row["charged_total_ms"])
    return None


def extract_rho_cost(rows: list[dict[str, Any]]) -> float | None:
    if not rows:
        return None
    totals = []
    for row in rows:
        if "total_ms" in row:
            totals.append(float(row["total_ms"]))
        elif isinstance(row.get("timing_breakdown_ms"), dict):
            timing = row["timing_breakdown_ms"]
            if "total_ms" in timing:
                totals.append(float(timing["total_ms"]))
    if not totals:
        return None
    return sum(totals) / len(totals)


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
        "ic_cost": ic_cost,
        "rho_cost": rho_cost,
        "automorphism_discount": automorphism_discount(int(beat["n"])),
        "all_stages_charged_same_series": True,
        "verdict": "DRAFT_PENDING_INDEPENDENT_VALIDATION",
        "claim_boundary": (
            "Public synthetic Koblitz fixture comparison only. Not key recovery, "
            "not asymptotic sub-rho, not an imported-point attack, and not a "
            "ledger promotion until independent validation and schema PASS."
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
        "fixture_hash": sha256(REPO / "examples/koblitz_rank_fixture.rs"),
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
    require(fixtures > 0, "fixtures must be positive")

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
        env.setdefault("KIC_INCREMENTAL_RANK_CROSSCHECK", "1")
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

        direct_obs = run_timed(direct_cmd, env=env, cwd=REPO)
        direct_obs["seed"] = direct_seed
        (run / "logs/direct.stdout.jsonl").write_text(direct_obs["stdout"])
        (run / "logs/direct.stderr.txt").write_text(direct_obs["stderr"])
        write_json(
            run / "receipts/direct.resource.json",
            {
                k: direct_obs[k]
                for k in (
                    "exit_code",
                    "whole_process_wall_ms",
                    "children_cpu_user_ms",
                    "children_cpu_system_ms",
                    "command",
                    "seed",
                )
            },
        )

        rho_obs = run_timed(rho_cmd, env=env, cwd=REPO)
        rho_obs["seed"] = rho_seed
        (run / "logs/rho.stdout.jsonl").write_text(rho_obs["stdout"])
        (run / "logs/rho.stderr.txt").write_text(rho_obs["stderr"])
        write_json(
            run / "receipts/rho.resource.json",
            {
                k: rho_obs[k]
                for k in (
                    "exit_code",
                    "whole_process_wall_ms",
                    "children_cpu_user_ms",
                    "children_cpu_system_ms",
                    "command",
                    "seed",
                )
            },
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
        help="Override fixture count (default from beat; use 1024 for full n37/n41)",
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
