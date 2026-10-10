"""Freeze the exact source-minimal Modal transfer and bounded pilot resources."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

from prepare_bundle import make_bundle, canonical
from pilot_modal import expected_arguments


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
STUDY = HERE.parent / "FREEZE.json"
BUILD = HERE / "build_receipt.json"
BUNDLE = HERE / "public_bundle.json"
FREEZE = HERE / "PILOT_FREEZE.json"
TRANSFER = HERE / "TRANSFER_MANIFEST.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def make() -> tuple[dict, dict]:
    study_raw = STUDY.read_bytes()
    study = json.loads(study_raw)
    build_raw = BUILD.read_bytes()
    build = json.loads(build_raw)
    bundle_raw = BUNDLE.read_bytes()
    assert bundle_raw == canonical(make_bundle()) + b"\n"
    assert b"fixture_scalar" not in bundle_raw and b"scalar_derivation" not in bundle_raw
    assert build["kind"] == "s3_crosscurve_modal_linux_build_v1"
    assert build["exit_code"] == 0 and build["source_commit"] == study["source_commit"]
    assert sha((HERE / "build.log").read_bytes()) == build["log_sha256"]
    root_oid = subprocess.check_output(["git", "rev-parse", study["source_commit"] + "^{tree}"],
                                       cwd=ROOT, text=True).strip()
    assert root_oid == study["source_tree_oid"]
    src_oid = subprocess.check_output(["git", "rev-parse", study["source_commit"] + ":src"],
                                      cwd=ROOT, text=True).strip()
    assert src_oid == build["source_tree_git_oid"]
    binary = {}
    for arm, name in (("reference", "s3_pair_candidate_frozen"),
                      ("prepared", "s3_partial_prepared_targetspan")):
        item = build["binaries"][name]
        path = Path(item["path"])
        assert path.is_file() and path.stat().st_size == item["bytes"]
        assert sha(path.read_bytes()) == item["sha256"]
        assert "ELF 64-bit LSB" in item["file"] and "x86-64" in item["file"]
        binary[arm] = item
    import modal
    freeze = {
        "kind": "s3_crosscurve_modal_pilot_freeze_v1",
        "extension_freeze_sha256": sha(study_raw),
        "source_commit": study["source_commit"],
        "source_root_tree_git_oid": root_oid,
        "source_src_tree_git_oid": src_oid,
        "build_receipt_sha256": sha(build_raw),
        "build_log_sha256": build["log_sha256"],
        "build_script_sha256": sha((HERE / "build_linux.py").read_bytes()),
        "bundle_builder_sha256": sha((HERE / "prepare_bundle.py").read_bytes()),
        "public_bundle_sha256": sha(bundle_raw),
        "binary_sha256": {arm: binary[arm]["sha256"] for arm in ("reference", "prepared")},
        "binary_bytes": {arm: binary[arm]["bytes"] for arm in ("reference", "prepared")},
        "pilot_runner_sha256": sha((HERE / "pilot_modal.py").read_bytes()),
        "remote_runner_sha256": sha((HERE / "remote_pilot.py").read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze.py").read_bytes()),
        "prepare_script_sha256": sha(Path(__file__).read_bytes()),
        "modal_sdk_version": modal.__version__,
        "transferred_file_count": 3,
        "transferred_bytes": sum(item["bytes"] for item in binary.values()) + len(bundle_raw),
        "remote_runner_sent_as_sandbox_command": True,
        "requires_explicit_transfer_approval": True,
        "cpu_request_physical_cores": 16,
        "cpu_limit_physical_cores": 16,
        "logical_execution_cpus": 16,
        "memory_request_mib": 32768,
        "memory_limit_mib": 32768,
        "sandbox_timeout_seconds": 1200,
        "workload_deadline_seconds": 1080,
        "per_run_timeout_seconds": 120,
        "minimum_run_seconds": 30,
        "sandbox_creation_attempts": 1,
        "block_network": True,
        "vm_runtime": True,
        "region": None,
        "pilot_target_indices": [25, 28],
        "public_target_count": 4,
        "planned_run_count": 16,
        "paired_order": "ABBA_or_BAAB_by_frozen_curve_and_target_index",
        "operation_count_boundary": study["operation_count_boundary"],
        "perf_counter_boundary": "separate inherited process events attached after child spawn; enabled/running times retained",
        "perf_min_running_fraction_for_ipc": "0.99",
        "single_target_wall_gain_threshold": 2,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    files = [
        {"remote_path": f"/opt/s3/{arm}", "sha256": binary[arm]["sha256"],
         "bytes": binary[arm]["bytes"]} for arm in ("reference", "prepared")
    ] + [{"remote_path": "/opt/s3/public_bundle.json",
          "sha256": sha(bundle_raw), "bytes": len(bundle_raw)}]
    argument_text = json.dumps(expected_arguments(freeze), sort_keys=True,
                               separators=(",", ":")).encode()
    runner_bytes = (HERE / "remote_pilot.py").stat().st_size
    transfer = {
        "kind": "s3_crosscurve_modal_transfer_manifest_v1",
        "pilot_freeze_sha256": sha(json.dumps(freeze, sort_keys=True, indent=2).encode() + b"\n"),
        "file_count": 3, "files": files,
        "transferred_bytes": freeze["transferred_bytes"],
        "fixture_scalars_in_bundle": False,
        "source_repository_transferred": False,
        "remote_runner_sent_as_sandbox_command": True,
        "remote_runner_source_sha256": freeze["remote_runner_sha256"],
        "remote_runner_source_bytes": runner_bytes,
        "command_argument_sha256": sha(argument_text),
        "command_argument_bytes": len(argument_text),
        "total_payload_bytes": freeze["transferred_bytes"] + runner_bytes + len(argument_text),
    }
    return freeze, transfer


if __name__ == "__main__":
    if FREEZE.exists() or TRANSFER.exists():
        raise FileExistsError("Modal pilot freeze or transfer manifest already exists")
    freeze, transfer = make()
    FREEZE.write_text(json.dumps(freeze, sort_keys=True, indent=2) + "\n")
    TRANSFER.write_text(json.dumps(transfer, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"transferred_bytes": freeze["transferred_bytes"],
                      "remote_runner_bytes": transfer["remote_runner_source_bytes"],
                      "planned_runs": freeze["planned_run_count"],
                      "freeze_sha256": sha(FREEZE.read_bytes())}, sort_keys=True))
