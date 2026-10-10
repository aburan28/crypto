"""Create the immutable minimal-transfer Modal pilot manifest after cross-build."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from prepare_bundle import make_bundle, canonical


HERE = Path(__file__).resolve().parent
BUILD = Path("/private/tmp/s3-modal-minimal-build/receipt.json")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def main() -> None:
    output = HERE / "PILOT_FREEZE.json"
    if output.exists():
        raise FileExistsError(output)
    bundle = HERE / "public_bundle.json"
    assert bundle.read_bytes() == canonical(make_bundle()) + b"\n"
    assert b"fixture_scalar" not in bundle.read_bytes()
    build = json.loads(BUILD.read_bytes())
    assert build["exit_code"] == 0
    assert build["target"] == "x86_64-unknown-linux-gnu"
    binaries = {
        "baseline": build["binaries"]["s3_partial_prepared_parallel_targetspan"],
        "expanded": build["binaries"]["s3_partial_prepared_fullrank_targetspan"],
    }
    for binary in binaries.values():
        path = Path(binary["path"])
        assert sha(path.read_bytes()) == binary["sha256"]
        assert "ELF 64-bit LSB" in binary["file"] and "x86-64" in binary["file"]
    import modal

    plan = {
        "kind": "s3_modal_minimal_pilot_freeze_v2",
        "build_receipt_sha256": sha(BUILD.read_bytes()),
        "public_bundle_sha256": sha(bundle.read_bytes()),
        "pilot_runner_sha256": sha((HERE / "pilot_modal.py").read_bytes()),
        "remote_runner_sha256": sha((HERE / "remote_pilot.py").read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze.py").read_bytes()),
        "binary_sha256": {arm: item["sha256"] for arm, item in binaries.items()},
        "transferred_file_count": 3,
        "transferred_bytes": bundle.stat().st_size + sum(
            Path(item["path"]).stat().st_size for item in binaries.values()),
        "modal_sdk_version": modal.__version__,
        "image": "modal.Image.debian_slim",
        "cpu_request_physical_cores": 8,
        "cpu_limit_physical_cores": 8,
        "logical_execution_cpus": 16,
        "memory_request_mib": 32768,
        "memory_limit_mib": 32768,
        "sandbox_timeout_seconds": 1200,
        "workload_deadline_seconds": 1080,
        "minimum_run_seconds": 30,
        "per_run_timeout_seconds": 240,
        "sandbox_creation_attempts": 1,
        "block_network": True, "vm_runtime": True,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    output.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"path": str(output), "sha256": sha(output.read_bytes()),
                      "transferred_files": 3,
                      "transferred_bytes": plan["transferred_bytes"]}, sort_keys=True))


if __name__ == "__main__":
    main()
