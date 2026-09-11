#!/usr/bin/env python3
"""Rerun the two Stage-26 WDSat capacity failures on one CPU."""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import resource
import sys
import time
from typing import Any

import build_koblitz_phase_b_wdsat as wdsat_builder
import koblitz_stage26_affinity_inputs as packet_tool
import run_koblitz_blind_pdp_phase_b as phase_b
import run_koblitz_stage26_affinity_cell as stage26


REPO = Path(__file__).resolve().parents[1]
PROTOCOL = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-32-wdsat-capacity-protocol.json"
SOURCE_PROTOCOL = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-20-balanced-pdp-phase-b-protocol.json"
PACKET_SHA256 = "c937afb9b172d114768b0b96a4b5ccf66e91fbf278b77c58ed49f0fd74af37f7"
BACKEND_SHA256 = "667fc00681e114acaaf662a8c00c214354e1a08e236a6974fb08e52504711ea0"
FAILED_IDS = [
    "b-c7bd575fcc8ef0dabf0d34f4cf886394276c05b423285925b9bfbe209bf26a6f",
    "b-10c023de48c3028911de41f270180f35b151d8b94828f520285a211997f8fe62",
]
REQUIRED_BUFFERS = [32588, 32804]
SCHEMA = "koblitz_stage32_wdsat_capacity_correction.v1"
SEAL_SCHEMA = "koblitz_stage32_wdsat_capacity_correction_seal.v1"


class Stage32Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage32Error(message)


def load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"JSON object required: {path}")
    return value


def self_test() -> dict[str, Any]:
    protocol = load(PROTOCOL)
    policy = protocol["wdsat_build"]
    require(policy["frozen_config_sha256"] == "4a73c3b5a14ded98f749d282b597594ae3ac05355a39cadcb9345a90bf735352", "Stage-32 config hash changed")
    require(policy["limits"]["max_buffer_size"] == max(REQUIRED_BUFFERS), "Stage-32 buffer limit changed")
    require(len(set(FAILED_IDS)) == 2, "Stage-32 failed-ID set changed")
    return {"schema": "koblitz_stage32_self_test.v1", "status": "PASS", "checks": 3}


def run(args: argparse.Namespace) -> dict[str, Any]:
    require(not args.output.exists() and not args.output.is_symlink(), "output must be new")
    affinity = stage26.singleton_affinity(args.cpu)
    before_self = resource.getrusage(resource.RUSAGE_SELF)
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    sampler = stage26.TreeSampler()
    sampler.__enter__()
    packet_check = packet_tool.verify(args.packet)
    require(packet_check["packet_inventory_sha256"] == PACKET_SHA256, "Stage-32 packet identity changed")
    packet = load(args.packet / "packet-manifest.json")
    records = {row["blind_instance_id"]: row for row in packet["instances"]}
    require(all(blind_id in records for blind_id in FAILED_IDS), "Stage-32 input IDs are absent from the packet")
    wdsat = args.wdsat.resolve(strict=True)
    backend = args.backend.resolve(strict=True)
    require(phase_b.sha256_file(backend, "Stage-32 backend") == BACKEND_SHA256, "Stage-32 backend identity changed")
    identities = {
        "wdsat": phase_b.executable_identity(wdsat, "Stage-32 WDSat"),
        "backend": phase_b.executable_identity(backend, "Stage-32 backend"),
    }
    protocol, _ = phase_b.read_json(PROTOCOL, "Stage-32 capacity protocol")
    phase_b.validate_protocol(protocol)
    source_protocol, _ = phase_b.read_json(SOURCE_PROTOCOL, "Stage-26 source protocol")
    phase_b.validate_protocol(source_protocol)
    state = phase_b.git_state()
    build = wdsat_builder.validate_capsule(
        args.wdsat_build_seal.resolve(strict=True), protocol, identities["wdsat"], state,
        require_clean_implementation=True,
    )
    environment = phase_b.safe_child_environment()
    markers = phase_b.forbidden_material(source_protocol)
    timeout = source_protocol["execution"]["per_process_watchdog_seconds"]
    args.output.mkdir(parents=True)
    (args.output / "tasks").mkdir()
    phase_b.write_json_new(args.output / "affinity.json", affinity)
    results = []
    for index, blind_id in enumerate(FAILED_IDS):
        source_record = records[blind_id]
        task_root = args.output / "tasks" / f"{index:02d}-{blind_id}"
        task_root.mkdir()
        instance_root = task_root / "instance"
        stage26.copy_instance(args.packet, source_record, instance_root)
        manifest_path = instance_root / "manifest.json"
        manifest = load(manifest_path)
        before = phase_b.frozen_source_snapshot(instance_root, manifest)
        requirements = phase_b.wdsat_requirements(instance_root, manifest)
        require(requirements["max_buffer_size"] == REQUIRED_BUFFERS[index], "Stage-32 input buffer requirement changed")
        require(all(requirements[name] <= build["limits"][name] for name in build["limits"]), "Stage-32 WDSat build is undersized")
        source_process = phase_b.run_metered(
            role="source-verification",
            command=[str(backend), "verify-source", str(manifest_path)],
            cwd=instance_root, task_root=task_root,
            input_paths=[manifest_path, *phase_b.manifest_source_paths(instance_root, manifest)],
            timeout=timeout, meter=args.meter.resolve(strict=True), environment=environment,
            markers=markers, expected_executable_sha256=BACKEND_SHA256,
        )
        source_verification = phase_b.source_verification_result(source_process, manifest)
        phase_b.require_source_unchanged(instance_root, manifest, before, "Stage-32 source verification")
        branch_variables = ",".join(str(value) for value in range(1, 3 * source_record["ell"] + 1))
        process = phase_b.run_metered(
            role="wdsat",
            command=[str(wdsat), "-i", str(instance_root / "instance.anf"), "-g", branch_variables],
            cwd=instance_root, task_root=task_root, input_paths=[instance_root / "instance.anf"],
            timeout=timeout, meter=args.meter.resolve(strict=True), environment=environment,
            markers=markers, expected_executable_sha256=identities["wdsat"]["sha256"],
        )
        phase_b.require_source_unchanged(instance_root, manifest, before, "Stage-32 WDSat")
        row = phase_b.external_result(
            solver="wdsat", record=process, manifest=manifest, manifest_path=manifest_path,
            backend=backend, task_root=task_root, instance_root=instance_root, timeout=timeout,
            meter=args.meter.resolve(strict=True), environment=environment, markers=markers,
            backend_sha256=BACKEND_SHA256,
        )
        phase_b.require_source_unchanged(instance_root, manifest, before, "Stage-32 witness validation")
        require(row["status"] not in {"solver_error", "backend_error", "backend_contract_error"}, "corrected WDSat still returned a solver error")
        task = {
            "schema": "koblitz_stage32_wdsat_task.v1", "ordinal": index,
            "blind_instance_id": blind_id, "cell_id": source_record["cell_id"],
            "source_instance_id": source_record["source_instance_id"], "target": source_record["target"],
            "requirements": requirements, "source_before": before,
            "source_verification": source_verification, "result": row,
            "source_after": phase_b.frozen_source_snapshot(instance_root, manifest),
        }
        require(task["source_before"] == task["source_after"], "Stage-32 source changed")
        phase_b.write_json_new(task_root / "task-result.json", task)
        results.append(task)
    sampler.__exit__(None, None, None)
    elapsed = time.monotonic() - started
    after_self = resource.getrusage(resource.RUSAGE_SELF)
    after_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    outer_core = (
        after_self.ru_utime - before_self.ru_utime + after_self.ru_stime - before_self.ru_stime
        + after_children.ru_utime - before_children.ru_utime
        + after_children.ru_stime - before_children.ru_stime
    )
    result = {
        "schema": SCHEMA, "status": "complete_two_input_capacity_correction",
        "superseded_stage26_status": "solver_error", "inputs": FAILED_IDS,
        "corrected_statuses": dict(Counter(task["result"]["status"] for task in results)),
        "packet_verification": packet_check, "wdsat_build": build, "affinity": affinity,
        "outer_resources": {
            "single_core_elapsed_seconds": elapsed, "total_core_seconds": outer_core,
            "one_cpu_utilization_percent": 100.0 * outer_core / elapsed,
            "sampled_peak_process_tree_rss_bytes": sampler.peak,
            "process_tree_rss_samples": sampler.samples, "sample_interval_seconds": 0.02,
        },
        "truth_labels_present": False, "known_witnesses_present": False,
        "full_cost_gate_passed": False, "koblitz_index_calculus_sota": False,
    }
    phase_b.write_json_new(args.output / "result.json", result)
    inventory = phase_b.all_regular_inventory(args.output, {"result-seal.json"})
    payload = {"schema": SEAL_SCHEMA, "status": "correction_frozen", "inventory": inventory, "inventory_sha256": phase_b.canonical_sha256(inventory)}
    seal = {**payload, "seal_payload_sha256": phase_b.canonical_sha256(payload)}
    phase_b.write_json_new(args.output / "result-seal.json", seal)
    return seal


def verify(output: Path) -> dict[str, Any]:
    seal = load(output / "result-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256", None)
    require(seal.get("schema") == SEAL_SCHEMA and claimed == phase_b.canonical_sha256(payload), "Stage-32 seal is invalid")
    inventory = phase_b.all_regular_inventory(output, {"result-seal.json"})
    require(inventory == seal.get("inventory") and phase_b.canonical_sha256(inventory) == seal.get("inventory_sha256"), "Stage-32 inventory changed")
    result = load(output / "result.json")
    require(result.get("schema") == SCHEMA and result.get("status") == "complete_two_input_capacity_correction", "Stage-32 result is incomplete")
    require(result.get("inputs") == FAILED_IDS and result.get("truth_labels_present") is False, "Stage-32 input boundary changed")
    require(sum(result.get("corrected_statuses", {}).values()) == 2 and "solver_error" not in result.get("corrected_statuses", {}), "Stage-32 did not correct both errors")
    return {"schema": "koblitz_stage32_wdsat_correction_verification.v1", "status": "verified", "inventory_sha256": seal["inventory_sha256"], "corrected_statuses": result["corrected_statuses"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("self-test")
    run_parser = sub.add_parser("run")
    run_parser.add_argument("--packet", type=Path, required=True)
    run_parser.add_argument("--output", type=Path, required=True)
    run_parser.add_argument("--wdsat", type=Path, required=True)
    run_parser.add_argument("--wdsat-build-seal", type=Path, required=True)
    run_parser.add_argument("--backend", type=Path, required=True)
    run_parser.add_argument("--meter", type=Path, default=phase_b.DEFAULT_METER)
    run_parser.add_argument("--cpu", type=int)
    verify_parser = sub.add_parser("verify")
    verify_parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "self-test":
            value = self_test()
        elif args.command == "verify":
            value = verify(args.output.resolve())
        else:
            args.packet, args.output = args.packet.resolve(), args.output.resolve()
            value = run(args)
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, Stage32Error, phase_b.PhaseBError, packet_tool.Stage26InputError) as error:
        raise SystemExit(f"stage32-wdsat: {error}")


if __name__ == "__main__":
    main()
