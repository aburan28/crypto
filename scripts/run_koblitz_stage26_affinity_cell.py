#!/usr/bin/env python3
"""Run one truth-free Stage-26 Phase-B cell on one Linux CPU."""

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
import threading
import time
from typing import Any

import koblitz_stage26_affinity_inputs as packet_tool
import run_koblitz_blind_pdp_phase_b as phase_b


REPO = Path(__file__).resolve().parents[1]
DEFAULT_PROTOCOL = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-20-balanced-pdp-phase-b-protocol.json"
SCHEMA = "koblitz_stage26_affinity_cell_result.v1"
SEAL_SCHEMA = "koblitz_stage26_affinity_cell_seal.v1"


class Stage26Error(RuntimeError):
    pass


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise Stage26Error(f"JSON object required: {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    phase_b.write_json_new(path, value)


def singleton_affinity(cpu: int | None) -> dict[str, Any]:
    if platform.system() != "Linux" or not hasattr(os, "sched_getaffinity"):
        raise Stage26Error("Stage-26 execution requires Linux CPU affinity")
    initial = sorted(os.sched_getaffinity(0))
    selected = initial[0] if cpu is None else cpu
    if selected not in initial:
        raise Stage26Error("requested CPU is outside the allowed affinity set")
    os.sched_setaffinity(0, {selected})
    parent = sorted(os.sched_getaffinity(0))
    probe = subprocess.run(
        [sys.executable, "-c", "import json,os;print(json.dumps(sorted(os.sched_getaffinity(0))))"],
        text=True,
        capture_output=True,
        check=True,
    )
    child = json.loads(probe.stdout)
    if parent != [selected] or child != [selected]:
        raise Stage26Error("singleton CPU affinity was not inherited")
    return {
        "schema": "koblitz_stage26_affinity_receipt.v1",
        "platform": platform.system(),
        "uname": list(platform.uname()),
        "initial_allowed_cpus": initial,
        "selected_cpu": selected,
        "parent_effective": parent,
        "child_effective": child,
    }


def proc_tree_rss(root_pid: int) -> int:
    parents: dict[int, int] = {}
    rss: dict[int, int] = {}
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            lines = (entry / "status").read_text().splitlines()
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        fields = {line.split(":", 1)[0]: line.split(":", 1)[1].strip() for line in lines if ":" in line}
        try:
            pid = int(entry.name)
            parents[pid] = int(fields["PPid"].split()[0])
            rss[pid] = int(fields.get("VmRSS", "0 kB").split()[0]) * 1024
        except (KeyError, ValueError):
            continue
    descendants = {root_pid}
    changed = True
    while changed:
        changed = False
        for pid, parent in parents.items():
            if parent in descendants and pid not in descendants:
                descendants.add(pid)
                changed = True
    return sum(rss.get(pid, 0) for pid in descendants)


class TreeSampler:
    def __init__(self) -> None:
        self.stop = threading.Event()
        self.peak = 0
        self.samples = 0
        self.thread = threading.Thread(target=self._run, daemon=True)

    def _run(self) -> None:
        while not self.stop.is_set():
            try:
                self.peak = max(self.peak, proc_tree_rss(os.getpid()))
                self.samples += 1
            except OSError:
                pass
            self.stop.wait(0.02)

    def __enter__(self) -> "TreeSampler":
        self.thread.start()
        return self

    def __exit__(self, *_: object) -> None:
        self.stop.set()
        self.thread.join()
        self.peak = max(self.peak, proc_tree_rss(os.getpid()))
        self.samples += 1


def copy_instance(packet: Path, record: dict[str, Any], destination: Path) -> None:
    destination.mkdir()
    for identity in record["files"]:
        source = packet / identity["path"]
        target = destination / Path(identity["path"]).name
        data = source.read_bytes()
        if len(data) != identity["bytes"] or phase_b.sha256_bytes(data) != identity["sha256"]:
            raise Stage26Error("packet file identity changed before cell execution")
        phase_b.write_new(target, data)


def prepare_task(
    *, ordinal: int, record: dict[str, Any], packet: Path, output: Path,
    protocol: dict[str, Any], tools: dict[str, Path], identities: dict[str, dict[str, Any]],
    environment: dict[str, str],
) -> dict[str, Any]:
    task_root = output / "tasks" / f"{ordinal:06d}-{record['blind_instance_id']}"
    task_root.mkdir()
    instance_root = task_root / "instance"
    copy_instance(packet, record, instance_root)
    manifest_path = instance_root / "manifest.json"
    manifest = read_json(manifest_path)
    if manifest.get("blind_instance_id") != record["blind_instance_id"]:
        raise Stage26Error("packet record and manifest blind IDs differ")
    if manifest.get("source_instance", {}).get("id_blake3") != record["source_instance_id"]:
        raise Stage26Error("packet record and manifest source IDs differ")
    predicate = manifest.get("factor_base_predicate", {})
    if predicate.get("enumerates_target_subgroup") is not False or predicate.get("uses_discrete_log_labels") is not False:
        raise Stage26Error("factor-base predicate violates the algebraic contract")
    source_before = phase_b.frozen_source_snapshot(instance_root, manifest)
    markers = phase_b.forbidden_material(protocol)
    timeout = protocol["execution"]["per_process_watchdog_seconds"]
    verification_record = phase_b.run_metered(
        role="source-verification",
        command=[str(tools["backend"]), "verify-source", str(manifest_path)],
        cwd=instance_root,
        task_root=task_root,
        input_paths=[manifest_path, *phase_b.manifest_source_paths(instance_root, manifest)],
        timeout=timeout,
        meter=tools["meter"],
        environment=environment,
        markers=markers,
        expected_executable_sha256=identities["backend"]["sha256"],
    )
    phase_b.require_source_unchanged(instance_root, manifest, source_before, "source verification")
    source_verification = phase_b.source_verification_result(verification_record, manifest)
    export_result = {
        "schema": "koblitz_stage26_frozen_export_task.v1",
        "ordinal": ordinal,
        "blind_instance_id": record["blind_instance_id"],
        "cell_id": record["cell_id"],
        "target": record["target"],
        "source_instance_id": record["source_instance_id"],
        "input_origin": "exact authenticated Stage-20 export; exporter not rerun",
        "source_artifacts_before": source_before,
        "source_verification": source_verification,
    }
    write_json(task_root / "input-result.json", export_result)
    return {
        "index": ordinal,
        "instance": record,
        "task_root": task_root,
        "instance_root": instance_root,
        "manifest_path": manifest_path,
        "manifest": manifest,
        "source_before": source_before,
        "export_result": export_result,
    }


def process_resources(output: Path) -> dict[str, Any]:
    metrics = [read_json(path)["metrics"] for path in sorted(output.glob("tasks/**/*.metrics.json"))]
    return {
        "process_receipts": len(metrics),
        "total_core_seconds": round(math.fsum(row["total_core_seconds"] for row in metrics), 12),
        "summed_process_wall_seconds": round(math.fsum(row["wall_seconds"] for row in metrics), 12),
        "maximum_individual_process_rss_bytes": max(row["peak_rss_bytes"] for row in metrics),
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    if args.output.exists() or args.output.is_symlink():
        raise Stage26Error("output must be new")
    affinity = singleton_affinity(args.cpu)
    self_before = resource.getrusage(resource.RUSAGE_SELF)
    children_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.monotonic()
    sampler = TreeSampler()
    sampler.__enter__()
    packet_check = packet_tool.verify(args.packet)
    manifest = read_json(args.packet / "packet-manifest.json")
    if packet_check["packet_inventory_sha256"] != args.packet_inventory_sha256:
        raise Stage26Error("packet inventory differs from the required workflow identity")
    if args.cell not in manifest["cells"] or len(manifest["cells"][args.cell]) != 40:
        raise Stage26Error("cell is absent or does not contain exactly 40 inputs")
    protocol, _ = phase_b.read_json(args.protocol, "Stage-20 Phase-B protocol")
    phase_b.validate_protocol(protocol)
    tools = {
        "backend": args.backend.resolve(strict=True),
        "wdsat": args.wdsat.resolve(strict=True),
        "cryptominisat": args.cryptominisat.resolve(strict=True),
        "meter": args.meter.resolve(strict=True),
    }
    identities = {name: phase_b.executable_identity(path, name, executable=name != "meter") for name, path in tools.items()}
    environment = phase_b.safe_child_environment()
    records_by_id = {row["blind_instance_id"]: row for row in manifest["instances"]}
    records = [records_by_id[item] for item in manifest["cells"][args.cell]]
    args.output.mkdir(parents=True)
    (args.output / "tasks").mkdir()
    write_json(args.output / "affinity.json", affinity)
    plan = {
        "schema": "koblitz_stage26_affinity_cell_plan.v1",
        "cell_id": args.cell,
        "instance_count": len(records),
        "backend_order": list(phase_b.BACKENDS),
        "packet_verification": packet_check,
        "tool_identities": identities,
        "affinity": affinity,
        "truth_scoring_status": "withheld_until_all_four_cells_are_frozen",
        "single_thread_requested": True,
    }
    write_json(args.output / "execution-plan.json", plan)
    results = []
    for ordinal, record in enumerate(records):
        prepared = prepare_task(
            ordinal=ordinal, record=record, packet=args.packet, output=args.output,
            protocol=protocol, tools=tools, identities=identities, environment=environment,
        )
        results.append(phase_b.solve_one_task(
            prepared=prepared, protocol=protocol, tools=tools,
            identities=identities, environment=environment,
        ))
    sampler.__exit__(None, None, None)
    elapsed = time.monotonic() - started
    self_after = resource.getrusage(resource.RUSAGE_SELF)
    children_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    total_core = (
        self_after.ru_utime - self_before.ru_utime + self_after.ru_stime - self_before.ru_stime
        + children_after.ru_utime - children_before.ru_utime
        + children_after.ru_stime - children_before.ru_stime
    )
    statuses = {name: Counter() for name in phase_b.BACKENDS}
    conflicts = {name: [] for name in phase_b.BACKENDS}
    for task in results:
        for row in task["backends"]:
            statuses[row["solver"]][row["status"]] += 1
            if row.get("conflicts") is not None:
                conflicts[row["solver"]].append(row["conflicts"])
    result = {
        "schema": SCHEMA,
        "status": "complete_truth_free_single_cpu_cell",
        "cell_id": args.cell,
        "instances": 40,
        "backend_outcomes": 120,
        "status_counts": {name: dict(sorted(value.items())) for name, value in statuses.items()},
        "conflicts": {name: {"reported": len(value), "sum": sum(value), "values": value} for name, value in conflicts.items()},
        "affinity": affinity,
        "packet_verification": packet_check,
        "process_resources": process_resources(args.output),
        "outer_resources": {
            "single_core_elapsed_seconds": elapsed,
            "total_core_seconds": total_core,
            "one_cpu_utilization_percent": 100.0 * total_core / elapsed,
            "sampled_peak_process_tree_rss_bytes": sampler.peak,
            "process_tree_rss_samples": sampler.samples,
            "sample_interval_seconds": 0.02,
            "scope": "packet integrity verification, cell setup, source verification, and all solver execution after singleton affinity",
        },
        "truth_labels_present": False,
        "known_witnesses_present": False,
        "factor_base_logs_known_by_construction": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
        "claim_boundary": "public synthetic PDP solver matrix only",
    }
    write_json(args.output / "result.json", result)
    inventory = phase_b.all_regular_inventory(args.output, {"result-seal.json"})
    payload = {
        "schema": SEAL_SCHEMA,
        "status": "cell_frozen",
        "cell_id": args.cell,
        "packet_inventory_sha256": packet_check["packet_inventory_sha256"],
        "inventory": inventory,
        "inventory_sha256": phase_b.canonical_sha256(inventory),
    }
    seal = {**payload, "seal_payload_sha256": phase_b.canonical_sha256(payload)}
    write_json(args.output / "result-seal.json", seal)
    return seal


def verify(output: Path) -> dict[str, Any]:
    seal = read_json(output / "result-seal.json")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256")
    if seal.get("schema") != SEAL_SCHEMA or claimed != phase_b.canonical_sha256(payload):
        raise Stage26Error("cell seal is invalid")
    inventory = phase_b.all_regular_inventory(output, {"result-seal.json"})
    if inventory != seal["inventory"] or phase_b.canonical_sha256(inventory) != seal["inventory_sha256"]:
        raise Stage26Error("cell inventory changed")
    result = read_json(output / "result.json")
    if result.get("schema") != SCHEMA or result.get("status") != "complete_truth_free_single_cpu_cell":
        raise Stage26Error("cell result is incomplete")
    if result.get("instances") != 40 or result.get("backend_outcomes") != 120:
        raise Stage26Error("cell result count changed")
    if result.get("truth_labels_present") is not False or result.get("known_witnesses_present") is not False:
        raise Stage26Error("cell result violates truth-free execution")
    if result.get("koblitz_index_calculus_sota") is not False:
        raise Stage26Error("cell result overclaims SOTA")
    return {"schema": "koblitz_stage26_affinity_cell_verification.v1", "status": "verified", "cell_id": result["cell_id"], "inventory_sha256": seal["inventory_sha256"]}


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    sub = root.add_subparsers(dest="command", required=True)
    run_p = sub.add_parser("run")
    run_p.add_argument("--packet", type=Path, required=True)
    run_p.add_argument("--packet-inventory-sha256", required=True)
    run_p.add_argument("--cell", required=True)
    run_p.add_argument("--output", type=Path, required=True)
    run_p.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    run_p.add_argument("--backend", type=Path, required=True)
    run_p.add_argument("--wdsat", type=Path, required=True)
    run_p.add_argument("--cryptominisat", type=Path, required=True)
    run_p.add_argument("--meter", type=Path, default=phase_b.DEFAULT_METER)
    run_p.add_argument("--cpu", type=int)
    verify_p = sub.add_parser("verify")
    verify_p.add_argument("--output", type=Path, required=True)
    return root


def main() -> None:
    args = parser().parse_args()
    try:
        value = run(args) if args.command == "run" else verify(args.output.resolve())
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, Stage26Error, phase_b.PhaseBError) as error:
        raise SystemExit(f"stage26-cell: {error}")


if __name__ == "__main__":
    main()
