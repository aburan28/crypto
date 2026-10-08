#!/usr/bin/env python3
"""Package and verify the exact truth-free Stage-20 instance exports for Stage 26."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
from typing import Any

import koblitz_stage23_terminal_evidence as custody
import run_koblitz_blind_pdp_phase_b as phaseb


MANIFEST_SCHEMA = "koblitz_stage26_affinity_input_manifest.v1"
SEAL_SCHEMA = "koblitz_stage26_affinity_input_seal.v1"
INSTANCE_FILES = ("manifest.json", "instance.anf", "instance.xor.cnf", "instance.magma")
EXPORT_NAMES = {
    "instance.anf": "wdsat_anf",
    "instance.xor.cnf": "cryptominisat_xor_dimacs",
    "instance.magma": "magma_boolean_f4",
}
EXPECTED_CELLS = {
    "n31-l5-m3-standard-a1-f0",
    "n31-l5-m3-ggmp-a0-f0",
    "n41-l5-m3-standard-a1-f0",
    "n59-l9-m3-standard-a1-f0",
}
FORBIDDEN_KEYS = {
    "class",
    "class_label",
    "decomposable",
    "ground_truth",
    "known_witness",
    "oracle",
    "truth",
}


class Stage26InputError(RuntimeError):
    pass


def read_json(path: Path, context: str) -> tuple[dict[str, Any], bytes]:
    value, data = custody.read_json(path, context)
    if not isinstance(value, dict):
        raise Stage26InputError(f"{context} must be a JSON object")
    return value, data


def identity(path: Path, recorded_path: str | None = None) -> dict[str, Any]:
    return custody.identity(path, recorded_path=recorded_path)


def validate_no_truth(value: Any, context: str = "manifest") -> None:
    if isinstance(value, dict):
        for key, child in value.items():
            if key.lower() in FORBIDDEN_KEYS:
                raise Stage26InputError(f"{context} contains forbidden truth key {key}")
            validate_no_truth(child, f"{context}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            validate_no_truth(child, f"{context}[{index}]")


def validate_source_run(source_run: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    seal, _ = read_json(source_run / "run-seal.json", "source run seal")
    required = {
        "schema",
        "status",
        "created_at",
        "protocol_sha256",
        "source_phase_a_seal_sha256",
        "blind_bundle_sha256",
        "selected_instance_count",
        "full_instance_count",
        "backend_outcomes",
        "full_panel_complete",
        "inventory",
        "inventory_sha256",
        "truth_scoring_status",
        "seal_payload_sha256",
    }
    if set(seal) != required or seal["status"] != "solver_outputs_frozen":
        raise Stage26InputError("source run seal schema or status changed")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256")
    if claimed != phaseb.canonical_sha256(payload):
        raise Stage26InputError("source run seal self-hash is invalid")
    inventory = seal["inventory"]
    if not isinstance(inventory, list) or phaseb.canonical_sha256(inventory) != seal["inventory_sha256"]:
        raise Stage26InputError("source run inventory hash is invalid")
    actual = phaseb.all_regular_inventory(source_run, {"run-seal.json"})
    if actual != inventory:
        raise Stage26InputError("source run files changed after sealing")
    if (
        seal["selected_instance_count"] != 160
        or seal["full_instance_count"] != 160
        or seal["backend_outcomes"] != 480
        or seal["full_panel_complete"] is not True
        or seal["truth_scoring_status"] != "withheld_until_separate_post_run_step"
    ):
        raise Stage26InputError("source run is not the complete truth-blind Stage-20 panel")
    return seal, inventory


def expected_artifact_record(entry: dict[str, Any], filename: str) -> dict[str, Any]:
    artifacts = entry.get("artifacts")
    if not isinstance(artifacts, dict):
        raise Stage26InputError("export inventory entry lacks artifacts")
    if filename == "manifest.json":
        record = artifacts.get("manifest")
    else:
        exports = artifacts.get("exports")
        record = exports.get(EXPORT_NAMES[filename]) if isinstance(exports, dict) else None
    if not isinstance(record, dict):
        raise Stage26InputError(f"export inventory lacks {filename}")
    return record


def validate_instance(
    *,
    instance_root: Path,
    blind: dict[str, Any],
    exported: dict[str, Any],
    packet_prefix: str,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if exported.get("blind_instance_id") != blind.get("blind_instance_id"):
        raise Stage26InputError("blind/export instance order changed")
    actual_names = {path.name for path in instance_root.iterdir()}
    if actual_names != set(INSTANCE_FILES):
        raise Stage26InputError(f"instance file grammar changed: {sorted(actual_names)}")
    manifest, _ = read_json(instance_root / "manifest.json", "instance manifest")
    validate_no_truth(manifest)
    for field in ("blind_instance_id", "n", "ell", "m", "curve_a", "target"):
        if manifest.get(field) != blind.get(field):
            raise Stage26InputError(f"instance manifest differs from blind bundle field {field}")
    if manifest.get("factor_base_predicate", {}).get("enumerates_target_subgroup") is not False:
        raise Stage26InputError("instance factor-base predicate enumerates the target subgroup")
    if manifest.get("factor_base_predicate", {}).get("uses_discrete_log_labels") is not False:
        raise Stage26InputError("instance factor-base predicate uses discrete-log labels")
    records = []
    for filename in INSTANCE_FILES:
        path = instance_root / filename
        data = custody.regular_bytes(path, f"instance file {filename}")
        expected = expected_artifact_record(exported, filename)
        if len(data) != expected.get("bytes") or custody.sha256(data) != expected.get("sha256"):
            raise Stage26InputError(f"{filename} differs from the frozen export inventory")
        if filename != "manifest.json":
            manifest_export = manifest.get("exports", {}).get(EXPORT_NAMES[filename], {})
            if custody.blake3_bytes(data).hex() != expected.get("manifest_blake3"):
                raise Stage26InputError(f"{filename} BLAKE3 differs from export inventory")
            if custody.blake3_bytes(data).hex() != manifest_export.get("blake3"):
                raise Stage26InputError(f"{filename} BLAKE3 differs from instance manifest")
        records.append(
            {
                "path": f"{packet_prefix}/{filename}",
                "bytes": len(data),
                "sha256": custody.sha256(data),
            }
        )
    return manifest, records


def validate_bundle_and_exports(
    blind_bundle: dict[str, Any], export_inventory: dict[str, Any]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if blind_bundle.get("schema") != "koblitz_pdp_blind_bundle.v1":
        raise Stage26InputError("blind bundle schema changed")
    blind = blind_bundle.get("instances")
    exported = export_inventory.get("source_artifacts")
    if (
        blind_bundle.get("instance_count") != 160
        or export_inventory.get("selected_instances") != 160
        or not isinstance(blind, list)
        or not isinstance(exported, list)
        or len(blind) != 160
        or len(exported) != 160
    ):
        raise Stage26InputError("blind/export inventory does not contain 160 instances")
    ids = [item.get("blind_instance_id") for item in blind]
    export_blind_ids = [item.get("blind_instance_id") for item in exported]
    source_ids = [item.get("source_instance_id") for item in exported]
    if (
        len(set(ids)) != 160
        or ids != export_blind_ids
        or source_ids != export_inventory.get("source_instance_ids")
    ):
        raise Stage26InputError("blind/export instance identity order changed")
    cells: dict[str, int] = {}
    for item in blind:
        validate_no_truth(item, "blind instance")
        cell = item.get("cell_id")
        if cell not in EXPECTED_CELLS:
            raise Stage26InputError(f"unexpected Stage-20 cell {cell}")
        cells[cell] = cells.get(cell, 0) + 1
    if cells != {cell: 40 for cell in EXPECTED_CELLS}:
        raise Stage26InputError(f"cell allocation changed: {cells}")
    return blind, exported


def package(source_run: Path, blind_bundle_path: Path, output: Path) -> dict[str, Any]:
    source_run = custody.real_directory(source_run.resolve(), "source Stage-20 run")
    if output.exists() or output.is_symlink():
        raise Stage26InputError("packet output must be new")
    seal, source_inventory = validate_source_run(source_run)
    blind_bundle, blind_bytes = read_json(blind_bundle_path.resolve(), "blind bundle")
    export_inventory, export_bytes = read_json(
        source_run / "export-inventory.json", "source export inventory"
    )
    blind, exported = validate_bundle_and_exports(blind_bundle, export_inventory)
    if custody.sha256(blind_bytes) != seal["blind_bundle_sha256"]:
        raise Stage26InputError("blind bundle hash differs from source run seal")

    custody.create_new_directory(output, "Stage-26 input packet")
    custody.write_new(output / "source/run-seal.json", custody.regular_bytes(source_run / "run-seal.json", "source run seal"))
    custody.write_new(output / "source/export-inventory.json", export_bytes)
    custody.write_new(output / "source/blind-bundle.json", blind_bytes)
    instance_records = []
    cells: dict[str, list[str]] = {cell: [] for cell in sorted(EXPECTED_CELLS)}
    for index, (blind_item, export_item) in enumerate(zip(blind, exported)):
        blind_id = blind_item["blind_instance_id"]
        source_instance = source_run / "tasks" / f"{index:06d}-{blind_id}" / "instance"
        packet_relative = f"instances/{index:06d}-{blind_id}"
        manifest, file_records = validate_instance(
            instance_root=source_instance,
            blind=blind_item,
            exported=export_item,
            packet_prefix=packet_relative,
        )
        for filename in INSTANCE_FILES:
            custody.write_new(
                output / packet_relative / filename,
                custody.regular_bytes(source_instance / filename, f"source {filename}"),
            )
        cell = blind_item["cell_id"]
        cells[cell].append(blind_id)
        instance_records.append(
            {
                "ordinal": index,
                "blind_instance_id": blind_id,
                "cell_id": cell,
                "n": blind_item["n"],
                "ell": blind_item["ell"],
                "m": blind_item["m"],
                "curve_a": blind_item["curve_a"],
                "factor_index": blind_item["factor_index"],
                "basis": blind_item["basis"],
                "target": blind_item["target"],
                "source_instance_id": manifest["source_instance"]["id_blake3"],
                "files": file_records,
            }
        )
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "status": "truth_free_exact_stage20_exports",
        "source_run": {
            "run_seal": identity(output / "source/run-seal.json", "source/run-seal.json"),
            "inventory_sha256": seal["inventory_sha256"],
            "inventory_records": len(source_inventory),
        },
        "blind_bundle": identity(output / "source/blind-bundle.json", "source/blind-bundle.json"),
        "export_inventory": identity(output / "source/export-inventory.json", "source/export-inventory.json"),
        "instances": instance_records,
        "cells": cells,
        "truth_labels_present": False,
        "known_witnesses_present": False,
        "target_subgroup_enumerated": False,
        "discrete_log_labels_present": False,
    }
    custody.write_json_new(output / "packet-manifest.json", manifest)
    inventory = custody.tree_inventory(output, excluded={"packet-seal.json"})
    seal_payload = {
        "schema": SEAL_SCHEMA,
        "status": "packet_frozen",
        "manifest": identity(output / "packet-manifest.json", "packet-manifest.json"),
        "inventory": inventory,
        "inventory_sha256": custody.canonical_sha256(inventory),
    }
    packet_seal = dict(seal_payload)
    packet_seal["seal_payload_sha256"] = custody.canonical_sha256(seal_payload)
    custody.write_json_new(output / "packet-seal.json", packet_seal)
    verify(output)
    return packet_seal


def verify(packet: Path) -> dict[str, Any]:
    packet = custody.real_directory(packet.resolve(), "Stage-26 input packet")
    seal, _ = read_json(packet / "packet-seal.json", "packet seal")
    required = {
        "schema", "status", "manifest", "inventory", "inventory_sha256", "seal_payload_sha256"
    }
    if set(seal) != required or seal["schema"] != SEAL_SCHEMA or seal["status"] != "packet_frozen":
        raise Stage26InputError("packet seal schema or status changed")
    payload = dict(seal)
    claimed = payload.pop("seal_payload_sha256")
    if claimed != custody.canonical_sha256(payload):
        raise Stage26InputError("packet seal self-hash is invalid")
    actual_inventory = custody.tree_inventory(packet, excluded={"packet-seal.json"})
    if actual_inventory != seal["inventory"] or custody.canonical_sha256(actual_inventory) != seal["inventory_sha256"]:
        raise Stage26InputError("packet inventory changed after sealing")
    manifest, _ = read_json(packet / "packet-manifest.json", "packet manifest")
    if identity(packet / "packet-manifest.json", "packet-manifest.json") != seal["manifest"]:
        raise Stage26InputError("packet manifest identity changed")
    expected_manifest_keys = {
        "schema", "status", "source_run", "blind_bundle", "export_inventory", "instances", "cells",
        "truth_labels_present", "known_witnesses_present", "target_subgroup_enumerated", "discrete_log_labels_present",
    }
    if set(manifest) != expected_manifest_keys or manifest["schema"] != MANIFEST_SCHEMA or manifest["status"] != "truth_free_exact_stage20_exports":
        raise Stage26InputError("packet manifest schema or status changed")
    for field in (
        "truth_labels_present", "known_witnesses_present", "target_subgroup_enumerated", "discrete_log_labels_present"
    ):
        if manifest[field] is not False:
            raise Stage26InputError(f"packet widened {field}")
    run_seal, _ = read_json(packet / "source/run-seal.json", "embedded source run seal")
    if identity(packet / "source/run-seal.json", "source/run-seal.json") != manifest["source_run"]["run_seal"]:
        raise Stage26InputError("embedded source run-seal identity changed")
    run_payload = dict(run_seal)
    run_claimed = run_payload.pop("seal_payload_sha256")
    if run_claimed != phaseb.canonical_sha256(run_payload):
        raise Stage26InputError("embedded source run seal hash is invalid")
    if run_seal["inventory_sha256"] != manifest["source_run"]["inventory_sha256"]:
        raise Stage26InputError("embedded source-run inventory binding changed")
    if len(run_seal["inventory"]) != manifest["source_run"]["inventory_records"]:
        raise Stage26InputError("embedded source-run inventory count changed")
    blind_bundle, _ = read_json(packet / "source/blind-bundle.json", "packet blind bundle")
    export_inventory, _ = read_json(packet / "source/export-inventory.json", "packet export inventory")
    if identity(packet / "source/blind-bundle.json", "source/blind-bundle.json") != manifest["blind_bundle"]:
        raise Stage26InputError("packet blind bundle identity changed")
    if identity(packet / "source/export-inventory.json", "source/export-inventory.json") != manifest["export_inventory"]:
        raise Stage26InputError("packet export inventory identity changed")
    blind, exported = validate_bundle_and_exports(blind_bundle, export_inventory)
    if custody.sha256(custody.regular_bytes(packet / "source/blind-bundle.json", "blind bundle")) != run_seal["blind_bundle_sha256"]:
        raise Stage26InputError("packet blind bundle differs from embedded run seal")
    if not isinstance(manifest["instances"], list) or len(manifest["instances"]) != 160:
        raise Stage26InputError("packet manifest instance count changed")
    expected_paths = {
        "packet-manifest.json", "source/run-seal.json", "source/export-inventory.json", "source/blind-bundle.json"
    }
    cell_ids: dict[str, list[str]] = {cell: [] for cell in sorted(EXPECTED_CELLS)}
    for index, (blind_item, export_item, record) in enumerate(zip(blind, exported, manifest["instances"])):
        blind_id = blind_item["blind_instance_id"]
        prefix = f"instances/{index:06d}-{blind_id}"
        instance_root = packet / prefix
        parsed, files = validate_instance(
            instance_root=instance_root,
            blind=blind_item,
            exported=export_item,
            packet_prefix=prefix,
        )
        expected_record = {
            "ordinal": index,
            "blind_instance_id": blind_id,
            "cell_id": blind_item["cell_id"],
            "n": blind_item["n"],
            "ell": blind_item["ell"],
            "m": blind_item["m"],
            "curve_a": blind_item["curve_a"],
            "factor_index": blind_item["factor_index"],
            "basis": blind_item["basis"],
            "target": blind_item["target"],
            "source_instance_id": parsed["source_instance"]["id_blake3"],
            "files": files,
        }
        if record != expected_record:
            raise Stage26InputError(f"packet instance record {index} is not derivable")
        cell_ids[blind_item["cell_id"]].append(blind_id)
        expected_paths.update(f"{prefix}/{filename}" for filename in INSTANCE_FILES)
    if manifest["cells"] != cell_ids:
        raise Stage26InputError("packet cell partition changed")
    actual_paths = {item["path"] for item in actual_inventory}
    if actual_paths != expected_paths:
        raise Stage26InputError("packet contains missing or extra files")
    return {
        "schema": "koblitz_stage26_affinity_input_verification.v1",
        "status": "truth_free_exact_exports_verified",
        "instances": 160,
        "cells": {cell: len(ids) for cell, ids in cell_ids.items()},
        "source_run_inventory_sha256": run_seal["inventory_sha256"],
        "packet_inventory_sha256": seal["inventory_sha256"],
        "truth_labels_present": False,
        "known_witnesses_present": False,
    }


def self_test() -> dict[str, Any]:
    clean = {"target": {"x": "1", "y": "2"}, "factor_base_predicate": {"uses_discrete_log_labels": False}}
    validate_no_truth(clean)
    try:
        validate_no_truth({"ground_truth": "forbidden"})
    except Stage26InputError:
        pass
    else:
        raise AssertionError("truth label was accepted")
    return {"schema": "koblitz_stage26_affinity_input_self_test.v1", "status": "PASS", "checks": 2}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    package_parser = sub.add_parser("package")
    package_parser.add_argument("--source-run", type=Path, required=True)
    package_parser.add_argument("--blind-bundle", type=Path, required=True)
    package_parser.add_argument("--output", type=Path, required=True)
    verify_parser = sub.add_parser("verify")
    verify_parser.add_argument("--packet", type=Path, required=True)
    sub.add_parser("self-test")
    args = parser.parse_args()
    try:
        if args.command == "package":
            result = package(args.source_run, args.blind_bundle, args.output)
        elif args.command == "verify":
            result = verify(args.packet)
        else:
            result = self_test()
        print(json.dumps(result, indent=2, sort_keys=True))
    except (Stage26InputError, custody.EvidenceError, phaseb.PhaseBError, OSError, ValueError, KeyError) as error:
        parser.exit(1, f"stage26-inputs: {error}\n")


if __name__ == "__main__":
    main()
