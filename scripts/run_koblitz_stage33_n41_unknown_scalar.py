#!/usr/bin/env python3
"""Build, run, seal, and verify the Stage-33 n=41 public-target workflow."""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
import os
from pathlib import Path
import platform
import resource
import shutil
import subprocess
import sys
import time
from typing import Any

import run_koblitz_blind_pdp_phase_b as custody
from run_koblitz_stage26_affinity_cell import Stage26Error, TreeSampler, singleton_affinity


REPO = Path(__file__).resolve().parents[1]
STAGE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
PARAMS = STAGE / "stage-33-n41-public-unknown-params.json"
LOCK = STAGE / "stage-20-rust-build/Cargo.lock"
METER = Path(__file__).with_name("process_meter.py")
TARGET_SEEDS = [41001, 41002, 41003, 41004, 41005]
SEED_MASKS = [1, 2, 4, 8, 16, 32]

BUILD_SCHEMA = "koblitz_stage33_ic_build.v1"
BUILD_SEAL_SCHEMA = "koblitz_stage33_ic_build_seal.v1"
RUN_SCHEMA = "koblitz_stage33_n41_public_unknown_scalar.v1"
RUN_SEAL_SCHEMA = "koblitz_stage33_n41_public_unknown_scalar_seal.v1"


class Stage33Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage33Error(message)


def read_json(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def git(source: Path, *arguments: str) -> str:
    process = subprocess.run(
        ["git", *arguments], cwd=source, text=True, capture_output=True, check=True
    )
    return process.stdout.strip()


def source_state(source: Path) -> dict[str, Any]:
    commit = git(source, "rev-parse", "HEAD")
    custody.require_hex40(commit, "Stage-33 source commit")
    porcelain = git(source, "status", "--porcelain=v1", "--untracked-files=all").splitlines()
    return {"commit": commit, "dirty": bool(porcelain), "porcelain": porcelain}


def validate_params(path: Path = PARAMS) -> dict[str, Any]:
    value = read_json(path, "Stage-33 parameters")
    require(value.get("schema_version") == 1, "Stage-33 parameter schema changed")
    require(value.get("curve") == {"degree": 41, "curve_a": 0}, "Stage-33 curve changed")
    require(value.get("summands") == 3 and value.get("solver") == "pair_table", "Stage-33 decomposition changed")
    require(value.get("seed") == 41031 and value.get("max_trials") == 200000, "Stage-33 workflow bounds changed")
    targets = value.get("targets")
    require(targets == [{"public_hash_seed": seed} for seed in TARGET_SEEDS], "Stage-33 public targets changed")
    require(
        all("known_log" not in target and "random_seed" not in target for target in targets),
        "Stage-33 parameters construct a target scalar",
    )
    search = value.get("factor_base")
    require(isinstance(search, dict) and search.get("mode") == "search", "Stage-33 does not discover its factor base")
    expected_search = {
        "family": "union",
        "min_dimension": 6,
        "max_dimension": 6,
        "max_abscissae": 4096,
        "union_min_seed": 6,
        "union_max_seed": 6,
        "union_samples": 0,
        "targets": 1024,
        "exhaustive_cap": 4096,
        "extra_relations": 2,
        "prune": False,
        "saturate": True,
    }
    require({key: search.get(key) for key in expected_search} == expected_search, "Stage-33 factor-base search changed")
    require(value.get("linear_algebra", {}).get("mode") == "sparse", "Stage-33 sparse linear algebra was disabled")
    require(value.get("collection") == {"unit_trials": 2048, "units": 4, "max_units": 64}, "Stage-33 collection bounds changed")
    baseline = value.get("baseline", {})
    require(
        baseline.get("rho") is True
        and baseline.get("rho_seed") == 5929075965221491791
        and baseline.get("rho_max_iterations") == 268435456,
        "Stage-33 rho control changed",
    )
    return value


def self_test() -> dict[str, Any]:
    value = validate_params()
    require(len(set(TARGET_SEEDS)) == 5, "Stage-33 target seeds are not distinct")
    require(value["factor_base"]["exhaustive_cap"] < 549756390943 - 1, "Stage-33 permits subgroup enumeration")
    return {
        "schema": "koblitz_stage33_self_test.v1",
        "status": "PASS",
        "checks": 11,
        "target_scalars_constructed_or_supplied": False,
        "factor_base_logs_known_by_construction": False,
        "target_subgroup_enumerated": False,
    }


def metered(
    *, role: str, command: list[str], cwd: Path, evidence: Path,
    timeout: int, environment: dict[str, str]
) -> dict[str, Any]:
    stdout = evidence / f"{role}.stdout"
    stderr = evidence / f"{role}.stderr"
    metrics = evidence / f"{role}.metrics.json"
    invocation = [
        str(Path(sys.executable).resolve()),
        str(METER.resolve(strict=True)),
        "--cwd", str(cwd),
        "--timeout", str(timeout),
        "--stdout", str(stdout),
        "--stderr", str(stderr),
        "--metrics", str(metrics),
        "--exclusive-create",
        "--",
        *command,
    ]
    wrapper = subprocess.run(
        invocation,
        stdin=subprocess.DEVNULL,
        close_fds=True,
        env=environment,
        check=False,
    )
    require(wrapper.returncode == 0 and metrics.is_file(), f"Stage-33 meter failed for {role}")
    value = read_json(metrics, f"Stage-33 {role} metrics")
    value["role"] = role
    return value


def process_complete(value: dict[str, Any]) -> bool:
    return (
        value.get("returncode") == 0
        and value.get("timed_out") is False
        and value.get("orphan_group_terminated") is False
        and isinstance(value.get("metrics"), dict)
    )


def linux_host_identity() -> dict[str, Any]:
    """Return a canonical performance-host receipt for Linux measurements."""
    require(platform.system() == "Linux", "host identity requires Linux")
    cpuinfo = Path("/proc/cpuinfo").read_text()
    first_block = cpuinfo.split("\n\n", 1)[0]
    parsed: dict[str, str] = {}
    for line in first_block.splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            parsed[key.strip()] = value.strip()
    fields = {
        key: parsed.get(key)
        for key in (
            "vendor_id",
            "cpu family",
            "model",
            "model name",
            "stepping",
            "microcode",
            "cache size",
            "physical id",
            "siblings",
            "cpu cores",
        )
    }
    require(
        all(isinstance(fields[key], str) and fields[key] for key in ("vendor_id", "model name", "cpu family", "model", "stepping")),
        "Linux CPU identity is incomplete",
    )
    cache = []
    cache_root = Path("/sys/devices/system/cpu/cpu0/cache")
    if cache_root.is_dir():
        for index in sorted(cache_root.glob("index*"), key=lambda path: path.name):
            row = {"index": index.name}
            for name in ("level", "type", "size", "coherency_line_size", "number_of_sets", "ways_of_associativity", "shared_cpu_list"):
                path = index / name
                row[name] = path.read_text().strip() if path.is_file() else None
            cache.append(row)
    def optional(path: str) -> str | None:
        item = Path(path)
        return item.read_text().strip() if item.is_file() else None
    payload = {
        "schema": "koblitz_linux_host_identity.v1",
        "kernel": {
            "system": platform.system(),
            "release": platform.release(),
            "version": platform.version(),
            "machine": platform.machine(),
        },
        "cpu": {
            **fields,
            "flags": sorted(parsed.get("flags", "").split()),
            "logical_cpu_count": os.cpu_count(),
            "effective_affinity": sorted(os.sched_getaffinity(0)),
            "scaling_governor_cpu0": optional("/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"),
            "scaling_min_freq_cpu0": optional("/sys/devices/system/cpu/cpu0/cpufreq/scaling_min_freq"),
            "scaling_max_freq_cpu0": optional("/sys/devices/system/cpu/cpu0/cpufreq/scaling_max_freq"),
        },
        "cache": cache,
        "cgroup": {
            "cpu_max": optional("/sys/fs/cgroup/cpu.max"),
            "cpuset_cpus_effective": optional("/sys/fs/cgroup/cpuset.cpus.effective"),
        },
    }
    return {**payload, "identity_sha256": custody.canonical_sha256(payload)}


def validate_linux_host_identity(value: dict[str, Any]) -> dict[str, Any]:
    require(isinstance(value, dict), "host identity must be an object")
    payload = dict(value)
    claimed = payload.pop("identity_sha256", None)
    require(
        payload.get("schema") == "koblitz_linux_host_identity.v1"
        and claimed == custody.canonical_sha256(payload),
        "host identity self-hash changed",
    )
    cpu = payload.get("cpu")
    require(
        isinstance(cpu, dict)
        and isinstance(cpu.get("model name"), str)
        and isinstance(cpu.get("effective_affinity"), list)
        and bool(cpu["effective_affinity"]),
        "host CPU identity is incomplete",
    )
    return value


def outer_resources(
    started: float,
    before_self: resource.struct_rusage,
    before_children: resource.struct_rusage,
    sampler: TreeSampler,
) -> dict[str, Any]:
    elapsed = time.monotonic() - started
    after_self = resource.getrusage(resource.RUSAGE_SELF)
    after_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    core = math.fsum(
        (
            after_self.ru_utime - before_self.ru_utime,
            after_self.ru_stime - before_self.ru_stime,
            after_children.ru_utime - before_children.ru_utime,
            after_children.ru_stime - before_children.ru_stime,
        )
    )
    return {
        "wall_seconds": elapsed,
        "total_core_seconds": core,
        "average_parallelism": core / elapsed if elapsed else None,
        "sampled_peak_process_tree_rss_bytes": sampler.peak,
        "process_tree_rss_samples": sampler.samples,
        "sample_interval_seconds": 0.02,
    }


def freeze(root: Path, seal_schema: str, status: str) -> dict[str, Any]:
    inventory = custody.all_regular_inventory(root, {"result-seal.json"})
    payload = {
        "schema": seal_schema,
        "status": status,
        "inventory": inventory,
        "inventory_sha256": custody.canonical_sha256(inventory),
    }
    seal = {**payload, "seal_payload_sha256": custody.canonical_sha256(payload)}
    custody.write_json_new(root / "result-seal.json", seal)
    return seal


def verify_seal(root: Path, seal_schema: str, status: str) -> dict[str, Any]:
    seal = read_json(root / "result-seal.json", "Stage-33 result seal")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(
        seal.get("schema") == seal_schema
        and seal.get("status") == status
        and claimed == custody.canonical_sha256(payload),
        "Stage-33 seal is invalid",
    )
    inventory = custody.all_regular_inventory(root, {"result-seal.json"})
    require(
        inventory == seal.get("inventory")
        and custody.canonical_sha256(inventory) == seal.get("inventory_sha256"),
        "Stage-33 sealed inventory changed",
    )
    return seal


def build_environment(output: Path, cargo: Path, rustc: Path) -> dict[str, str]:
    home = Path(os.environ.get("HOME", "/home/runner"))
    path_entries = [str(cargo.parent), "/usr/local/bin", "/usr/bin", "/bin"]
    environment = {
        "PATH": os.pathsep.join(dict.fromkeys(path_entries)),
        "HOME": str(home),
        "CARGO_HOME": str(output / "cargo-home"),
        "RUSTC": str(rustc),
        "CARGO_INCREMENTAL": "0",
        "RUSTFLAGS": "-Awarnings",
        "LANG": "C",
        "LC_ALL": "C",
        "TZ": "UTC",
    }
    if os.environ.get("RUSTUP_HOME"):
        environment["RUSTUP_HOME"] = os.environ["RUSTUP_HOME"]
    elif (home / ".rustup").is_dir():
        environment["RUSTUP_HOME"] = str(home / ".rustup")
    return environment


def build(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-33 build accounting requires Linux /proc")
    source = args.source.resolve(strict=True)
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-33 build output must be new")
    require(1 <= args.jobs <= 64 and args.timeout > 0, "Stage-33 build bounds are invalid")
    state = source_state(source)
    require(state["dirty"] is False, "Stage-33 build source must be clean")
    lock_source = source / LOCK.relative_to(REPO)
    require(lock_source.is_file(), "Stage-33 frozen Cargo.lock is missing")
    root_lock = source / "Cargo.lock"
    require(not root_lock.exists() and not root_lock.is_symlink(), "Stage-33 refuses to replace a root Cargo.lock")
    cargo_name, rustc_name = shutil.which("cargo"), shutil.which("rustc")
    require(cargo_name is not None and rustc_name is not None, "Stage-33 Rust toolchain is missing")
    # Keep the rustup shim invocation names. Resolving either symlink to the
    # shared `rustup` binary would change argv[0] and stop selecting cargo or
    # rustc correctly.
    cargo, rustc = Path(cargo_name).absolute(), Path(rustc_name).absolute()
    output.mkdir(parents=True)
    evidence = output / "evidence"
    evidence.mkdir()
    (output / "cargo-home").mkdir()
    environment = build_environment(output, cargo, rustc)
    target = output / "target"
    binary_path = target / "release/ic"
    packaged = output / "bin/ic"
    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    fetch: dict[str, Any] | None = None
    compilation: dict[str, Any] | None = None
    binary: dict[str, Any] | None = None
    with TreeSampler() as sampler:
        try:
            shutil.copyfile(lock_source, root_lock)
            fetch = metered(
                role="cargo-fetch",
                command=[str(cargo), "fetch", "--locked"],
                cwd=source,
                evidence=evidence,
                timeout=args.timeout,
                environment=environment,
            )
            if process_complete(fetch):
                compilation = metered(
                    role="cargo-build",
                    command=[
                        str(cargo), "build", "--release", "--locked", "--offline",
                        "--jobs", str(args.jobs), "--target-dir", str(target), "--bin", "ic",
                    ],
                    cwd=source,
                    evidence=evidence,
                    timeout=args.timeout,
                    environment=environment,
                )
            if compilation is not None and process_complete(compilation) and binary_path.is_file():
                packaged.parent.mkdir()
                shutil.copyfile(binary_path, packaged)
                packaged.chmod(0o755)
                binary = custody.executable_identity(packaged, "Stage-33 ic binary")
        finally:
            root_lock.unlink(missing_ok=True)
            shutil.rmtree(target, ignore_errors=True)
            shutil.rmtree(output / "cargo-home", ignore_errors=True)
    outer = outer_resources(started, before_self, before_children, sampler)
    after_state = source_state(source)
    complete = (
        process_complete(fetch or {})
        and process_complete(compilation or {})
        and binary is not None
        and after_state == state
    )
    processes = [row for row in (fetch, compilation) if row is not None]
    result = {
        "schema": BUILD_SCHEMA,
        "status": "complete_clean_build" if complete else "incomplete_build",
        "source_state": state,
        "source_objects": {
            "Cargo.toml": git(source, "rev-parse", "HEAD:Cargo.toml"),
            str(LOCK.relative_to(REPO)): git(source, "rev-parse", f"HEAD:{LOCK.relative_to(REPO)}"),
            "src": git(source, "rev-parse", "HEAD:src"),
        },
        "frozen_lock": custody.executable_identity(lock_source, "Stage-33 Cargo.lock", executable=False),
        "requested_parallel_jobs": args.jobs,
        "toolchain": {
            "cargo": custody.executable_identity(cargo, "Stage-33 cargo"),
            "rustc": custody.executable_identity(rustc, "Stage-33 rustc"),
        },
        "processes": processes,
        "process_totals": {
            "total_core_seconds": math.fsum(row["metrics"]["total_core_seconds"] for row in processes),
            "summed_process_wall_seconds": math.fsum(row["metrics"]["wall_seconds"] for row in processes),
            "maximum_individual_process_rss_bytes": max((row["metrics"]["peak_rss_bytes"] for row in processes), default=0),
        },
        "outer_resources": outer,
        "binary": binary,
        "accounting_boundary": {
            "included": "fresh Cargo dependency acquisition, four-job clean release build, driver overhead, and sampled aggregate process-tree RSS",
            "excluded": ["preinstalled operating system", "preinstalled Rust toolchain acquisition"],
        },
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    custody.write_json_new(output / "result.json", result)
    freeze(output, BUILD_SEAL_SCHEMA, "build_frozen")
    require(complete, "Stage-33 clean build did not complete; retain the incomplete sealed output")
    return result


def validate_build(root: Path) -> dict[str, Any]:
    verify_seal(root, BUILD_SEAL_SCHEMA, "build_frozen")
    result = read_json(root / "result.json", "Stage-33 build result")
    require(result.get("schema") == BUILD_SCHEMA and result.get("status") == "complete_clean_build", "Stage-33 build is incomplete")
    require(result.get("source_state", {}).get("dirty") is False, "Stage-33 build source was dirty")
    require(result.get("requested_parallel_jobs") == 4, "Stage-33 build parallelism changed")
    require(len(result.get("processes", [])) == 2 and all(process_complete(row) for row in result["processes"]), "Stage-33 build processes are incomplete")
    binary = root / "bin/ic"
    identity = custody.executable_identity(binary, "Stage-33 ic binary")
    require(
        all(result.get("binary", {}).get(key) == identity[key] for key in ("bytes", "sha256")),
        "Stage-33 build binary changed",
    )
    return result


def union_spec(spec: Any) -> bool:
    if not isinstance(spec, dict):
        return False
    if spec.get("kind") == "frobenius_union":
        return spec.get("seed_masks") == SEED_MASKS
    return spec.get("kind") == "two_torsion_saturated" and union_spec(spec.get("parent"))


def validate_workflow(value: dict[str, Any], workflow_dir: Path) -> dict[str, Any]:
    require(value.get("operation") == "workflow" and value.get("status") == "complete", "Stage-33 workflow is incomplete")
    require(value.get("evidence_scope") == "public_hash_unknown_scalar", "Stage-33 workflow evidence label is wrong")
    require(value.get("degree") == 41 and value.get("curve_a") == 0 and value.get("solver") == "pair_table", "Stage-33 workflow curve or solver changed")
    factor = value.get("factor_base", {})
    require(union_spec(factor.get("spec")), "Stage-33 selected a non-algebraic factor base")
    require(all(isinstance(factor.get(name), int) and factor[name] > 0 for name in ("points", "columns", "signed_orbits", "abscissae")), "Stage-33 factor-base census is incomplete")
    stages = value.get("stages")
    require(isinstance(stages, list), "Stage-33 stage reports are missing")
    by_name = {row.get("stage"): row for row in stages}
    require(set(by_name) == {"select", "collect", "logs", "solve", "baseline"}, "Stage-33 stage inventory changed")
    search = by_name["select"].get("search", {})
    require(
        by_name["select"].get("ran") is True
        and search.get("targets") == 1024
        and search.get("exhaustive_targets") is False
        and search.get("candidates_scored") == 2,
        "Stage-33 factor-base discovery did not run within the frozen boundary",
    )
    require(by_name["collect"].get("status") == "complete" and by_name["collect"].get("trials_total", 0) > 0, "Stage-33 relation collection is incomplete")
    logs = by_name["logs"]
    linear = logs.get("linear_algebra", {})
    require(
        logs.get("status") == "complete"
        and logs.get("columns") == factor["columns"]
        and logs.get("relations", 0) > 0
        and logs.get("rejected") == 0
        and linear.get("mode") == "sparse"
        and linear.get("sparse", {}).get("filter", {}).get("uncovered_columns") == 0,
        "Stage-33 factor-base logarithm solve is incomplete",
    )
    solutions = value.get("solutions", {})
    items = solutions.get("items") if isinstance(solutions, dict) else None
    require(solutions.get("count") == 5 and solutions.get("verified") == 5 and isinstance(items, list) and len(items) == 5, "Stage-33 did not verify five target solutions")
    for index, item in enumerate(items):
        target = item.get("target", {})
        require(
            item.get("index") == index
            and item.get("expected") == "not_constructed"
            and item.get("verified") is True
            and isinstance(item.get("recovered"), str)
            and target.get("kind") == "public_hash_to_curve_cofactor"
            and target.get("target_scalar_constructed") is False
            and target.get("public_hash_seed") == TARGET_SEEDS[index],
            f"Stage-33 public target {index} changed",
        )
    baseline = by_name["baseline"].get("vs_rho", {})
    require(
        baseline.get("claim_boundary") == "public_hash_unknown_scalar"
        and baseline.get("n") == 41
        and baseline.get("targets") == 5
        and baseline.get("ic", {}).get("verified") == 5
        and baseline.get("rho", {}).get("verified") == 5
        and baseline.get("verdict", {}).get("charged_crossover") is False
        and baseline.get("verdict", {}).get("amortised_crossover") is False,
        "Stage-33 IC/rho control is incomplete",
    )
    factor_document = read_json(workflow_dir / "factor_base.json", "Stage-33 factor-base document")
    log_document = read_json(workflow_dir / "logs.json", "Stage-33 log database")
    solution_document = read_json(workflow_dir / "solutions.json", "Stage-33 solution database")
    baseline_document = read_json(workflow_dir / "baseline.json", "Stage-33 rho database")
    require(union_spec(factor_document.get("spec")) and factor_document.get("degree") == 41, "Stage-33 factor-base document changed")
    require(log_document.get("spec") == factor_document.get("spec") and len(log_document.get("columns", [])) == factor["columns"], "Stage-33 log database changed")
    require(solution_document.get("solutions") == items, "Stage-33 solution document differs from the report")
    require(baseline_document == baseline, "Stage-33 rho document differs from the report")
    return baseline


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(platform.system() == "Linux", "Stage-33 scientific run requires Linux CPU affinity")
    output = args.output.resolve()
    require(not output.exists() and not output.is_symlink(), "Stage-33 run output must be new")
    params = args.params.resolve(strict=True)
    validate_params(params)
    build_root = args.build.resolve(strict=True)
    build_result = validate_build(build_root)
    binary = build_root / "bin/ic"
    affinity = singleton_affinity(args.cpu)
    output.mkdir(parents=True)
    parameter_copy = output / "params.json"
    custody.write_new(parameter_copy, params.read_bytes())
    evidence = output / "evidence"
    evidence.mkdir()
    workflow_dir = output / "workflow"
    environment = custody.safe_child_environment()
    command = [
        str(binary.resolve(strict=True)), "--json", "workflow",
        "--params", str(parameter_copy), "--dir", str(workflow_dir),
    ]
    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    with TreeSampler() as sampler:
        process = metered(
            role="workflow",
            command=command,
            cwd=REPO,
            evidence=evidence,
            timeout=args.timeout,
            environment=environment,
        )
    outer = outer_resources(started, before_self, before_children, sampler)
    workflow: dict[str, Any] | None = None
    baseline: dict[str, Any] | None = None
    stdout = evidence / "workflow.stdout"
    try:
        workflow = read_json(stdout, "Stage-33 workflow stdout")
        baseline = validate_workflow(workflow, workflow_dir)
    except (OSError, ValueError, json.JSONDecodeError, Stage33Error):
        workflow = None
        baseline = None
    complete = process_complete(process) and workflow is not None and baseline is not None
    build_outer = build_result["outer_resources"]
    result = {
        "schema": RUN_SCHEMA,
        "status": "complete_n41_public_unknown_scalar" if complete else "incomplete_n41_public_unknown_scalar",
        "params": custody.executable_identity(parameter_copy, "Stage-33 parameters", executable=False),
        "build_binding": {
            "result": custody.executable_identity(build_root / "result.json", "Stage-33 build result", executable=False),
            "seal": custody.executable_identity(build_root / "result-seal.json", "Stage-33 build seal", executable=False),
            "binary": build_result["binary"],
            "source_state": build_result["source_state"],
            "requested_parallel_jobs": build_result["requested_parallel_jobs"],
            "process_totals": build_result["process_totals"],
            "outer_resources": build_outer,
        },
        "affinity": affinity,
        "workflow_process": process,
        "outer_resources": {
            **outer,
            "single_core_elapsed_seconds": outer["wall_seconds"],
            "one_cpu_utilization_percent": 100.0 * outer["total_core_seconds"] / outer["wall_seconds"],
        },
        "workflow_result": workflow,
        "target_scalars_constructed_or_supplied": False,
        "factor_base_logs_known_by_construction": False,
        "factor_base_log_source": "verified relation collection plus group-certified sparse linear algebra",
        "target_subgroup_enumerated_for_factor_base": False,
        "decomposition_solver": "direct pair table; no Groebner basis",
        "sat_conflicts": None,
        "sat_conflict_semantics": "not applicable to the Stage-33 pair-table end-to-end arm; same-instance SAT conflicts are charged in Phase B",
        "charged_total_core_seconds_available": build_outer["total_core_seconds"] + outer["total_core_seconds"],
        "charged_sequential_wall_seconds_available": build_outer["wall_seconds"] + outer["wall_seconds"],
        "maximum_sampled_process_tree_rss_bytes": max(
            build_outer["sampled_peak_process_tree_rss_bytes"],
            outer["sampled_peak_process_tree_rss_bytes"],
        ),
        "stage33_bounded_cost_components_complete": complete,
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    custody.write_json_new(output / "result.json", result)
    freeze(output, RUN_SEAL_SCHEMA, "run_frozen")
    require(complete, "Stage-33 n=41 workflow did not complete; retain the incomplete sealed output")
    return result


def verify(build_root: Path, output: Path) -> dict[str, Any]:
    build_result = validate_build(build_root)
    seal = verify_seal(output, RUN_SEAL_SCHEMA, "run_frozen")
    result = read_json(output / "result.json", "Stage-33 run result")
    require(result.get("schema") == RUN_SCHEMA and result.get("status") == "complete_n41_public_unknown_scalar", "Stage-33 result is incomplete")
    require(result.get("affinity", {}).get("parent_effective") == result.get("affinity", {}).get("child_effective") == [result.get("affinity", {}).get("selected_cpu")], "Stage-33 singleton affinity changed")
    require(result.get("target_scalars_constructed_or_supplied") is False and result.get("factor_base_logs_known_by_construction") is False, "Stage-33 label boundary changed")
    require(result.get("target_subgroup_enumerated_for_factor_base") is False, "Stage-33 subgroup-enumeration boundary changed")
    require(result.get("sat_conflicts") is None, "Stage-33 pair-table arm acquired SAT conflicts")
    workflow = result.get("workflow_result")
    require(isinstance(workflow, dict), "Stage-33 workflow report is absent")
    baseline = validate_workflow(workflow, output / "workflow")
    binding = result.get("build_binding", {})
    require(binding.get("binary") == build_result.get("binary") and binding.get("source_state") == build_result.get("source_state"), "Stage-33 build binding changed")
    require(result.get("stage33_bounded_cost_components_complete") is True, "Stage-33 bounded cost components are incomplete")
    require(result.get("full_cost_gate_passed") is False and result.get("koblitz_index_calculus_sota") is False, "Stage-33 claim widened")
    return {
        "schema": "koblitz_stage33_verification.v1",
        "status": "n41_public_unknown_scalar_build_and_run_verified",
        "source_commit": build_result["source_state"]["commit"],
        "factor_base": workflow["factor_base"],
        "factor_base_search": next(row["search"] for row in workflow["stages"] if row["stage"] == "select"),
        "relations": next(row["relations"] for row in workflow["stages"] if row["stage"] == "logs"),
        "linear_algebra": next(row["linear_algebra"] for row in workflow["stages"] if row["stage"] == "logs"),
        "targets_verified": workflow["solutions"]["verified"],
        "descent_trials": baseline["ic"]["descent_trials_total"],
        "rho_iterations": baseline["rho"]["iterations_total"],
        "rho_over_ic_charged_ratio": baseline["ratio"]["charged"],
        "rho_over_ic_amortised_ratio": baseline["ratio"]["amortised"],
        "single_core_elapsed_seconds": result["outer_resources"]["single_core_elapsed_seconds"],
        "workflow_outer_core_seconds": result["outer_resources"]["total_core_seconds"],
        "charged_total_core_seconds_available": result["charged_total_core_seconds_available"],
        "charged_sequential_wall_seconds_available": result["charged_sequential_wall_seconds_available"],
        "maximum_sampled_process_tree_rss_bytes": result["maximum_sampled_process_tree_rss_bytes"],
        "sat_conflicts": None,
        "score_inventory_sha256": seal["inventory_sha256"],
        "independent_external_reproduction_satisfied": False,
        "licensed_magma_complete": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("self-test")
    build_parser = sub.add_parser("build")
    build_parser.add_argument("--source", type=Path, required=True)
    build_parser.add_argument("--output", type=Path, required=True)
    build_parser.add_argument("--jobs", type=int, default=4)
    build_parser.add_argument("--timeout", type=int, default=3600)
    run_parser = sub.add_parser("run")
    run_parser.add_argument("--build", type=Path, required=True)
    run_parser.add_argument("--params", type=Path, default=PARAMS)
    run_parser.add_argument("--output", type=Path, required=True)
    run_parser.add_argument("--timeout", type=int, default=3600)
    run_parser.add_argument("--cpu", type=int)
    verify_parser = sub.add_parser("verify")
    verify_parser.add_argument("--build", type=Path, required=True)
    verify_parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "self-test":
            value = self_test()
        elif args.command == "build":
            value = build(args)
        elif args.command == "run":
            value = run(args)
        else:
            value = verify(args.build.resolve(strict=True), args.output.resolve(strict=True))
        print(json.dumps(value, indent=2, sort_keys=True))
    except (
        OSError,
        ValueError,
        KeyError,
        subprocess.CalledProcessError,
        Stage26Error,
        Stage33Error,
        custody.PhaseBError,
    ) as error:
        raise SystemExit(f"stage33-n41: {error}")


if __name__ == "__main__":
    main()
