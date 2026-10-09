"""Run one bounded Modal VM host audit without uploading repository files.

Only the generic text of remote_probe.py is passed as a Python command
argument. No local directory, source tree, fixture, image layer, or volume is
transferred. Sandbox creation is attempted once; an ambiguous failure is
recorded and never retried automatically.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time


HERE = Path(__file__).resolve().parent
PLAN = HERE / "FREEZE.json"
REMOTE = HERE / "remote_probe.py"
ANALYZER = HERE / "analyze_receipt.py"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def preflight() -> tuple[dict, str]:
    plan = json.loads(PLAN.read_bytes())
    assert plan["kind"] == "source_free_modal_vm_cpu_probe_freeze_v1"
    assert plan["local_runner_sha256"] == sha(Path(__file__).read_bytes())
    assert plan["remote_code_sha256"] == sha(REMOTE.read_bytes())
    assert plan["remote_code_bytes"] == REMOTE.stat().st_size
    assert plan["analyzer_sha256"] == sha(ANALYZER.read_bytes())
    assert plan["cpu_physical_cores"] == 16
    assert plan["memory_mib"] == 32768
    assert plan["sandbox_timeout_seconds"] == 120
    assert plan["vm_runtime"] is True
    assert plan["block_network"] is True
    assert plan["sandbox_creation_attempts"] == 1
    assert plan["repository_files_transferred"] == []
    assert plan["cryptanalysis_sources_or_fixtures_transferred"] == []
    assert plan["modal_volume"] is None
    return plan, REMOTE.read_text()


def launch(output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    plan, remote_code = preflight()
    import modal

    assert getattr(modal, "__version__", None) == plan["modal_sdk_version"]
    started = time.time_ns()
    receipt = {"kind": "source_free_modal_vm_cpu_probe_receipt_v1",
               "freeze_sha256": sha(PLAN.read_bytes()),
               "remote_code_sha256": plan["remote_code_sha256"],
               "modal_sdk_version": modal.__version__,
               "requested_physical_cpu_cores": plan["cpu_physical_cores"],
               "requested_memory_mib": plan["memory_mib"],
               "sandbox_timeout_seconds": plan["sandbox_timeout_seconds"],
               "repository_files_transferred": [],
               "sandbox_id": None, "started_unix_ns": started,
               "finished_unix_ns": None, "remote_exit_code": None,
               "stderr": None, "error": None, "remote": None,
               "host_exclusive_partition_verified": False,
               "controlled_speedup": None}
    sandbox = None
    try:
        app = modal.App("s3-source-free-cpu-fidelity-probe")
        with app.run():
            sandbox = modal.Sandbox.create(
                "python3", "-c", remote_code,
                app=app, image=modal.Image.debian_slim(),
                cpu=plan["cpu_physical_cores"],
                memory=plan["memory_mib"],
                timeout=plan["sandbox_timeout_seconds"],
                block_network=plan["block_network"],
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
                    raise RuntimeError("sandbox returned no JSON probe result")
                remote = json.loads(lines[-1])
                assert remote["kind"] == "source_free_modal_vm_cpu_fidelity_probe_v1"
                assert remote["host_exclusive_partition_verified"] is False
                receipt["remote"] = remote
                if sandbox.returncode != 0:
                    raise RuntimeError(f"remote probe exited {sandbox.returncode}")
            finally:
                sandbox.terminate()
    except BaseException as exc:
        receipt["error"] = f"{type(exc).__name__}:{str(exc)[:1000]}"
    receipt["finished_unix_ns"] = time.time_ns()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(output), "sandbox_id": receipt["sandbox_id"],
                      "error": receipt["error"],
                      "remote_exit_code": receipt["remote_exit_code"]},
                     sort_keys=True), flush=True)
    if receipt["error"]:
        raise SystemExit(1)
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview", action="store_true")
    parser.add_argument("--launch", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.preview:
        plan, _ = preflight()
        print(json.dumps({"remote_code_sha256": plan["remote_code_sha256"],
                          "remote_code_bytes": plan["remote_code_bytes"],
                          "repository_files_transferred": [],
                          "cpu_physical_cores": plan["cpu_physical_cores"],
                          "memory_mib": plan["memory_mib"],
                          "sandbox_timeout_seconds": plan["sandbox_timeout_seconds"]},
                         sort_keys=True))
    elif args.launch:
        if args.output is None:
            parser.error("--launch requires --output")
        launch(args.output)
    else:
        parser.error("choose --preview or --launch")
