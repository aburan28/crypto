#!/usr/bin/env python3
"""Run the frozen stage-13 target-matched PDP replication panel."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import verify_stage13_pdp_panel as verify


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DEFAULT_PROTOCOL = HERE / "stage-13-pdp-panel-protocol.json"
DEFAULT_RUNNER = REPO / "scripts" / "run_koblitz_pdp_matrix.py"
DEFAULT_METER = REPO / "scripts" / "process_meter.py"


def now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_state(path: Path = REPO, ignore: Path | None = None) -> dict:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=path, text=True, capture_output=True, check=True
    ).stdout.strip()
    status_command = ["git", "status", "--porcelain=v1", "--", "."]
    if ignore is not None:
        try:
            relative = ignore.resolve().relative_to(path.resolve())
        except ValueError:
            pass
        else:
            status_command.extend([f":(exclude){relative}", f":(exclude){relative}/**"])
    status = subprocess.run(
        status_command, cwd=path, text=True, capture_output=True, check=True
    ).stdout.splitlines()
    return {"commit": commit, "dirty": bool(status), "porcelain": status}


def atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def task_key(seed: int, cell: dict) -> str:
    return f"seed-{seed}/{cell['id']}"


def config(cell: dict) -> str:
    return ":".join(
        str(cell[field]) for field in ("n", "ell", "basis", "curve_a", "factor_index")
    )


def executable_identity(value: str | None) -> dict:
    if value is None:
        return {"requested": False, "available": False, "path": None, "sha256": None}
    resolved_text = shutil.which(value) if "/" not in value else str(Path(value).resolve())
    resolved = Path(resolved_text) if resolved_text else None
    available = resolved is not None and resolved.is_file()
    return {
        "requested": True,
        "available": available,
        "path": str(resolved) if resolved is not None else value,
        "sha256": sha256_file(resolved) if available else None,
    }


def source_identity(value: str | None) -> dict:
    if value is None:
        return {"requested": False, "available": False, "path": None, "git": None}
    resolved = Path(value).resolve()
    available = (resolved / "src" / "makefile").is_file()
    state = None
    if available:
        try:
            state = git_state(resolved)
        except subprocess.CalledProcessError:
            state = {"commit": None, "dirty": None, "porcelain": None}
    return {"requested": True, "available": available, "path": str(resolved), "git": state}


def implementation_identities(args: argparse.Namespace) -> dict:
    return {
        "python": executable_identity(sys.executable),
        "driver": {"path": str(Path(__file__).resolve()), "sha256": sha256_file(Path(__file__))},
        "verifier": {
            "path": str((HERE / "verify_stage13_pdp_panel.py").resolve()),
            "sha256": sha256_file(HERE / "verify_stage13_pdp_panel.py"),
        },
        "exporter": {"path": str(args.exporter.resolve()), "sha256": sha256_file(args.exporter)},
        "backend": {"path": str(args.backend.resolve()), "sha256": sha256_file(args.backend)},
        "runner": {"path": str(args.runner.resolve()), "sha256": sha256_file(args.runner)},
        "meter": {"path": str(args.meter.resolve()), "sha256": sha256_file(args.meter)},
    }


def requested_tool_identities(args: argparse.Namespace) -> dict:
    return {
        "wdsat_binary": executable_identity(args.wdsat),
        "wdsat_source": source_identity(args.wdsat_source),
        "cryptominisat": executable_identity(args.cryptominisat),
        "magma": executable_identity(args.magma),
    }


def assert_frozen_environment(frozen: dict, args: argparse.Namespace, output: Path) -> None:
    if implementation_identities(args) != frozen["implementation"]:
        raise verify.VerificationError("measurement implementation changed during the panel")
    if requested_tool_identities(args) != frozen["requested_tools"]:
        raise verify.VerificationError("solver or WDSat source identity changed during the panel")
    state = git_state(ignore=output)
    if state["commit"] != frozen["source_revision"]["commit"]:
        raise verify.VerificationError("source commit changed during the panel")
    if frozen["evidence_class"] == "scientific_candidate" and state["dirty"]:
        raise verify.VerificationError("scientific-candidate checkout became dirty during the panel")


def make_runner_command(args: argparse.Namespace, protocol: dict, seed: int, cell: dict, matrix: Path) -> list[str]:
    command = [
        str(Path(sys.executable).resolve()),
        str(args.runner.resolve()),
        "--output",
        str(matrix.resolve()),
        "--exporter",
        str(args.exporter.resolve()),
        "--backend",
        str(args.backend.resolve()),
        "--timeout",
        str(protocol["solver_policy"]["per_process_watchdog_seconds"]),
        "--conflicts",
        str(protocol["solver_policy"]["native_conflict_budget"]),
        "--configs",
        config(cell),
        "--seed",
        str(seed),
    ]
    for option, value, is_source in (
        ("--wdsat", args.wdsat, False),
        ("--wdsat-source", args.wdsat_source, True),
        ("--cryptominisat", args.cryptominisat, False),
        ("--magma", args.magma, False),
    ):
        if value is not None:
            normalized = str(Path(value).resolve()) if is_source else executable_identity(value)["path"]
            command.extend([option, normalized])
    return command


def selected_tasks(protocol: dict, only_cells: list[str], only_seeds: list[int], max_tasks: int | None):
    known_cells = {cell["id"] for cell in protocol["cells"]}
    known_seeds = set(protocol["replicate_seeds"])
    if set(only_cells) - known_cells:
        raise verify.VerificationError(f"unknown --only-cell values: {sorted(set(only_cells) - known_cells)}")
    if set(only_seeds) - known_seeds:
        raise verify.VerificationError(f"unknown --only-seed values: {sorted(set(only_seeds) - known_seeds)}")
    chosen = [
        (seed, cell)
        for seed, cell in verify.tasks(protocol)
        if (not only_cells or cell["id"] in only_cells)
        and (not only_seeds or seed in only_seeds)
    ]
    return chosen[:max_tasks] if max_tasks is not None else chosen


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--exporter", type=Path, required=True)
    parser.add_argument("--backend", type=Path, required=True)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--meter", type=Path, default=DEFAULT_METER)
    parser.add_argument("--wdsat")
    parser.add_argument("--wdsat-source")
    parser.add_argument("--cryptominisat")
    parser.add_argument("--magma")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true", help="for operational smoke runs only")
    parser.add_argument("--only-cell", action="append", default=[])
    parser.add_argument("--only-seed", action="append", type=int, default=[])
    parser.add_argument("--max-tasks", type=int)
    parser.add_argument("--plan", action="store_true", help="validate and print the exact task plan without running")
    args = parser.parse_args()

    protocol = verify.read_json(args.protocol)
    verify.validate_protocol(protocol)
    if not all(path.is_file() for path in (args.runner, args.meter, args.exporter, args.backend)):
        parser.error("runner, meter, built exporter, and built isolated backend must all exist")
    output = args.output.resolve()
    state = git_state(ignore=output)
    chosen = selected_tasks(protocol, args.only_cell, args.only_seed, args.max_tasks)
    if not chosen:
        parser.error("task selection is empty")
    if args.wdsat and args.wdsat_source:
        parser.error("choose one of --wdsat or --wdsat-source")
    protocol_hash = verify.canonical_sha256(protocol)
    frozen_tasks = [
        {
            "seed": seed,
            "cell": cell,
            "config": config(cell),
            "relative_output": str(verify.task_relpath(seed, cell) / "matrix"),
        }
        for seed, cell in verify.tasks(protocol)
    ]
    selection = [task_key(seed, cell) for seed, cell in chosen]
    implementations = implementation_identities(args)
    requested_tools = requested_tool_identities(args)
    auxiliary_dirty = any(
        item.get("git", {}).get("dirty") is True
        for item in requested_tools.values()
        if isinstance(item.get("git"), dict)
    )
    full_selection = len(chosen) == len(verify.tasks(protocol))
    evidence_class = (
        "scientific_candidate"
        if not state["dirty"] and not auxiliary_dirty and full_selection
        else "operational_smoke"
    )
    plan = {
        "schema": "koblitz_target_matched_pdp_execution_plan.v1",
        "protocol": str(args.protocol.resolve()),
        "protocol_sha256": protocol_hash,
        "repo": str(REPO),
        "source_revision": state,
        "implementation": implementations,
        "requested_tools": requested_tools,
        "evidence_class": evidence_class,
        "selected_task_count": len(chosen),
        "full_frozen_task_count": len(verify.tasks(protocol)),
        "scientific_panel_complete_if_run": evidence_class == "scientific_candidate",
        "frozen_tasks": frozen_tasks,
        "invocation_selection": selection,
    }
    if args.plan:
        print(json.dumps(plan, indent=2, sort_keys=True))
        return
    if state["dirty"] and not args.allow_dirty:
        parser.error("full evidence requires a clean checkout; use --allow-dirty only for labelled smoke runs")

    run_path = output / "panel-run.json"
    if output.exists() and not args.resume:
        parser.error("output exists; use a new path or --resume")
    if not output.exists():
        output.mkdir(parents=True)
        run = dict(plan)
        run.update(
            {
                "started_at": now(),
                "status": "running",
                "task_results": {},
                "invocation_history": [{"started_at": now(), "selection": selection}],
            }
        )
        atomic_json(run_path, run)
    else:
        if not run_path.is_file():
            parser.error("resume output lacks panel-run.json")
        run = verify.read_json(run_path)
        immutable = (
            "schema",
            "protocol_sha256",
            "repo",
            "implementation",
            "requested_tools",
            "full_frozen_task_count",
            "frozen_tasks",
        )
        changed = [name for name in immutable if run.get(name) != plan.get(name)]
        if run.get("source_revision", {}).get("commit") != state["commit"]:
            changed.append("source_revision.commit")
        if run.get("evidence_class") not in {"scientific_candidate", "operational_smoke"}:
            changed.append("evidence_class")
        if run.get("evidence_class") == "scientific_candidate" and state["dirty"]:
            changed.append("clean_checkout_required_for_scientific_resume")
        if changed:
            parser.error(f"resume changes frozen execution identity: {', '.join(changed)}")
        run.setdefault("invocation_history", []).append({"started_at": now(), "selection": selection})
        atomic_json(run_path, run)

    for seed, cell in chosen:
        try:
            assert_frozen_environment(run, args, output)
        except verify.VerificationError as error:
            run["status"] = "environment_changed"
            run["environment_error"] = str(error)
            atomic_json(run_path, run)
            raise SystemExit(2) from error
        key = task_key(seed, cell)
        task_dir = output / verify.task_relpath(seed, cell)
        matrix = task_dir / "matrix"
        if matrix.exists():
            if not args.resume:
                parser.error(f"unexpected existing matrix directory {matrix}")
            try:
                receipt = verify.verify_cell(protocol, output, seed, cell, write_receipt=True)
            except verify.VerificationError as error:
                parser.error(f"cannot resume over invalid existing task {key}: {error}")
            run["task_results"][key] = {
                "status": "verified_existing",
                "source_instance_sha256": receipt["source_instance_sha256"],
                "finished_at": now(),
            }
            atomic_json(run_path, run)
            continue
        if task_dir.exists():
            parser.error(
                f"refusing to overwrite incomplete prior attempt {task_dir}; "
                "retain it and use a new panel output path"
            )
        task_dir.mkdir(parents=True, exist_ok=True)
        runner_command = make_runner_command(args, protocol, seed, cell, matrix)
        invocation = {
            "seed": seed,
            "cell": cell,
            "runner_command": runner_command,
            "outer_watchdog_seconds": protocol["solver_policy"]["whole_cell_watchdog_seconds"],
            "started_at": now(),
        }
        atomic_json(task_dir / "invocation.json", invocation)
        meter_command = [
            str(Path(sys.executable).resolve()),
            str(args.meter.resolve()),
            "--cwd",
            str(REPO),
            "--timeout",
            str(protocol["solver_policy"]["whole_cell_watchdog_seconds"]),
            "--stdout",
            str(task_dir / "runner.stdout"),
            "--stderr",
            str(task_dir / "runner.stderr"),
            "--metrics",
            str(task_dir / "outer-metrics.json"),
            "--",
            *runner_command,
        ]
        started = time.perf_counter()
        completed = subprocess.run(meter_command, cwd=REPO, check=False)
        invocation.update(
            {
                "meter_command": meter_command,
                "meter_launcher_returncode": completed.returncode,
                "driver_observed_wall_seconds": time.perf_counter() - started,
                "finished_at": now(),
            }
        )
        atomic_json(task_dir / "invocation.json", invocation)
        try:
            assert_frozen_environment(run, args, output)
        except verify.VerificationError as error:
            run["task_results"][key] = {
                "status": "failed_or_invalid",
                "error": str(error),
                "finished_at": now(),
            }
            atomic_json(run_path, run)
            break
        try:
            receipt = verify.verify_cell(protocol, output, seed, cell, write_receipt=True)
        except verify.VerificationError as error:
            run["task_results"][key] = {
                "status": "failed_or_invalid",
                "error": str(error),
                "finished_at": now(),
            }
            atomic_json(run_path, run)
            continue
        run["task_results"][key] = {
            "status": "verified",
            "source_instance_sha256": receipt["source_instance_sha256"],
            "resource_complete": receipt["resource_complete"],
            "finished_at": now(),
        }
        atomic_json(run_path, run)

    summary = verify.summarize(protocol, output, write_receipts=True, allow_incomplete=True)
    atomic_json(output / "panel-summary.json", summary)
    if summary["panel_artifact_complete"]:
        run["status"] = (
            "scientific_candidate_artifact_set_complete"
            if summary["scientific_evidence_admissible"]
            else "operational_smoke_artifact_set_complete"
        )
    else:
        run["status"] = "partial"
    run["finished_at"] = now()
    run["summary"] = {
        "panel_artifact_complete": summary["panel_artifact_complete"],
        "solver_matrix_terminal_complete": summary["solver_matrix_terminal_complete"],
        "verified_tasks": summary["verified_tasks"],
        "per_backend_resources_complete": summary["per_backend_resources_complete"],
        "magma_executed": summary["magma_executed"],
        "evidence_class": summary["evidence_class"],
        "scientific_evidence_admissible": summary["scientific_evidence_admissible"],
        "full_solver_matrix_gate_passed": summary["full_solver_matrix_gate_passed"],
    }
    atomic_json(run_path, run)
    print(json.dumps(summary, indent=2, sort_keys=True))
    failed_selected = [
        key
        for key in selection
        if run["task_results"].get(key, {}).get("status") == "failed_or_invalid"
    ]
    if failed_selected:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
