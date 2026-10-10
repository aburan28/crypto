"""Upload exactly two frozen Linux binaries and one public-point bundle to Modal."""

from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import json
from pathlib import Path
import time

from prepare_bundle import make_bundle, canonical


HERE = Path(__file__).resolve().parent
PLAN = HERE / "PILOT_FREEZE.json"
BUNDLE = HERE / "public_bundle.json"
REMOTE = HERE / "remote_pilot.py"
BUILD_RECEIPT = Path("/private/tmp/s3-modal-minimal-build/receipt.json")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def preflight() -> tuple[dict, dict, dict]:
    plan = json.loads(PLAN.read_bytes())
    assert plan["kind"] == "s3_modal_minimal_pilot_freeze_v2"
    assert plan["pilot_runner_sha256"] == sha(Path(__file__).read_bytes())
    assert plan["remote_runner_sha256"] == sha(REMOTE.read_bytes())
    assert plan["public_bundle_sha256"] == sha(BUNDLE.read_bytes())
    assert BUNDLE.read_bytes() == canonical(make_bundle()) + b"\n"
    assert b"fixture_scalar" not in BUNDLE.read_bytes()
    build = json.loads(BUILD_RECEIPT.read_bytes())
    assert build["exit_code"] == 0
    assert plan["build_receipt_sha256"] == sha(BUILD_RECEIPT.read_bytes())
    assert build["library_source_tree_git_oid"] == "0a574800502e67882d8f8f88a846c1e80ee795fd"
    paths = {}
    for arm, name in (("baseline", "s3_partial_prepared_parallel_targetspan"),
                      ("expanded", "s3_partial_prepared_fullrank_targetspan")):
        info = build["binaries"][name]
        binary = Path(info["path"])
        assert binary.is_file() and sha(binary.read_bytes()) == info["sha256"]
        assert "ELF 64-bit LSB" in info["file"] and "x86-64" in info["file"]
        assert plan["binary_sha256"][arm] == info["sha256"]
        paths[arm] = binary
    assert plan["transferred_file_count"] == 3
    assert plan["transferred_bytes"] == (sum(path.stat().st_size for path in paths.values())
                                          + BUNDLE.stat().st_size)
    assert plan["cpu_request_physical_cores"] == 8
    assert plan["cpu_limit_physical_cores"] == 8
    assert plan["logical_execution_cpus"] == 16
    assert plan["memory_request_mib"] == 32768
    assert plan["memory_limit_mib"] == 32768
    assert plan["sandbox_timeout_seconds"] == 1200
    assert plan["per_run_timeout_seconds"] == 240
    assert plan["workload_deadline_seconds"] == 1080
    assert plan["block_network"] and plan["vm_runtime"]
    assert plan["sandbox_creation_attempts"] == 1
    bundle = json.loads(BUNDLE.read_bytes())
    assert bundle["run_order"] == [[41, "baseline"], [41, "expanded"],
                                   [41, "expanded"], [41, "baseline"],
                                   [53, "expanded"], [53, "baseline"],
                                   [53, "baseline"], [53, "expanded"]]
    return plan, build, paths


def launch(output: Path, region: str | None) -> None:
    plan, build, paths = preflight()
    import modal

    assert modal.__version__ == plan["modal_sdk_version"]
    output.parent.mkdir(parents=True, exist_ok=True)
    attempt = output.with_suffix(".attempt.json")
    if output.exists() or attempt.exists():
        raise FileExistsError("a Modal creation attempt already exists for this output")
    inventory = [{"remote": f"/opt/s3/{arm}", "sha256": plan["binary_sha256"][arm],
                  "bytes": paths[arm].stat().st_size} for arm in ("baseline", "expanded")]
    inventory.append({"remote": "/opt/s3/public_bundle.json",
                      "sha256": plan["public_bundle_sha256"],
                      "bytes": BUNDLE.stat().st_size})
    started = time.time_ns()
    attempt.write_text(json.dumps({"kind": "s3_modal_minimal_attempt_v1",
                                   "freeze_sha256": sha(PLAN.read_bytes()),
                                   "started_unix_ns": started,
                                   "transfer_inventory": inventory,
                                   "sandbox_creation_attempts": 1,
                                   "cpu_request_physical_cores": plan["cpu_request_physical_cores"],
                                   "cpu_limit_physical_cores": plan["cpu_limit_physical_cores"],
                                   "logical_execution_cpus": plan["logical_execution_cpus"],
                                   "memory_request_mib": plan["memory_request_mib"],
                                   "memory_limit_mib": plan["memory_limit_mib"],
                                   "sandbox_timeout_seconds": plan["sandbox_timeout_seconds"]},
                                  sort_keys=True, indent=2) + "\n")
    receipt = {"kind": "s3_modal_minimal_launch_receipt_v1",
               "freeze_sha256": sha(PLAN.read_bytes()),
               "build_receipt_sha256": sha(BUILD_RECEIPT.read_bytes()),
               "transfer_inventory": inventory,
               "modal_sdk_version": modal.__version__, "sandbox_id": None,
               "started_unix_ns": started, "finished_unix_ns": None,
               "remote_exit_code": None, "stderr": None,
               "error": None, "guest": None,
               "host_exclusive_partition_verified": False,
               "controlled_speedup": None}
    sandbox = None
    try:
        image = modal.Image.debian_slim()
        for arm in ("baseline", "expanded"):
            image = image.add_local_file(str(paths[arm]), f"/opt/s3/{arm}", copy=True)
        image = image.add_local_file(str(BUNDLE), "/opt/s3/public_bundle.json", copy=True)
        expected = {"bundle_sha256": plan["public_bundle_sha256"],
                    "binary_sha256": plan["binary_sha256"],
                    "logical_execution_cpus": plan["logical_execution_cpus"],
                    "workload_deadline_seconds": plan["workload_deadline_seconds"],
                    "minimum_run_seconds": plan["minimum_run_seconds"],
                    "per_run_timeout_seconds": plan["per_run_timeout_seconds"]}
        app = modal.App("s3-minimal-binary-n41-n53-pilot")
        with app.run():
            sandbox = modal.Sandbox.create(
                "python3", "-c", REMOTE.read_text(),
                json.dumps(expected, sort_keys=True, separators=(",", ":")),
                app=app, image=image,
                cpu=(plan["cpu_request_physical_cores"],
                     plan["cpu_limit_physical_cores"]),
                memory=(plan["memory_request_mib"], plan["memory_limit_mib"]),
                timeout=plan["sandbox_timeout_seconds"],
                region=region, block_network=plan["block_network"],
                experimental_options={"vm_runtime": plan["vm_runtime"]})
            receipt["sandbox_id"] = sandbox.object_id
            print(json.dumps({"sandbox_id": sandbox.object_id}), flush=True)
            try:
                sandbox.wait()
                stdout = sandbox.stdout.read()
                stderr = sandbox.stderr.read()
                receipt["remote_exit_code"] = sandbox.returncode
                receipt["stderr"] = stderr[-4000:]
                lines = [line for line in stdout.splitlines() if line.strip()]
                if not lines:
                    raise RuntimeError("Sandbox emitted no compressed guest receipt")
                envelope = json.loads(lines[-1])
                assert envelope["kind"] == "s3_modal_minimal_compressed_receipt_v1"
                assert envelope["encoding"] == "gzip_base64"
                assert sha(envelope["payload"].encode()) == envelope["sha256"]
                guest = json.loads(gzip.decompress(base64.b64decode(envelope["payload"])))
                assert guest["kind"] == "s3_modal_minimal_pilot_guest_receipt_v1"
                assert guest["bundle_sha256"] == plan["public_bundle_sha256"]
                assert guest["binary_sha256"] == plan["binary_sha256"]
                receipt["guest"] = guest
                if sandbox.returncode != 0:
                    raise RuntimeError(f"Sandbox exited {sandbox.returncode}")
            finally:
                sandbox.terminate()
    except BaseException as exc:
        receipt["error"] = f"{type(exc).__name__}:{str(exc)[:1000]}"
    receipt["finished_unix_ns"] = time.time_ns()
    output.write_text(json.dumps(receipt, sort_keys=True, separators=(",", ":")) + "\n")
    print(json.dumps({"output": str(output), "sandbox_id": receipt["sandbox_id"],
                      "error": receipt["error"],
                      "remote_exit_code": receipt["remote_exit_code"],
                      "guest_rows": len((receipt["guest"] or {}).get("rows", []))},
                     sort_keys=True), flush=True)
    if receipt["error"]:
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview", action="store_true")
    parser.add_argument("--launch", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--region")
    args = parser.parse_args()
    if args.preview:
        plan, _, paths = preflight()
        print(json.dumps({"files": 3, "transferred_bytes": plan["transferred_bytes"],
                          "binary_sha256": plan["binary_sha256"],
                          "bundle_sha256": plan["public_bundle_sha256"],
                          "cpu_physical_cores_request_and_limit": [
                              plan["cpu_request_physical_cores"],
                              plan["cpu_limit_physical_cores"]],
                          "guest_logical_execution_cpus": plan["logical_execution_cpus"],
                          "memory_mib_request_and_limit": [
                              plan["memory_request_mib"], plan["memory_limit_mib"]],
                          "sandbox_timeout_seconds": plan["sandbox_timeout_seconds"],
                          "paths": {arm: str(path) for arm, path in paths.items()}},
                         sort_keys=True))
    elif args.launch:
        if args.output is None:
            parser.error("--launch requires --output")
        launch(args.output, args.region)
    else:
        parser.error("choose --preview or --launch")
