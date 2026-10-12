#!/usr/bin/env python3
"""One frozen, fail-closed structured-orbit m5 construction attempt."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
from decimal import Decimal

HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[4]
PRIOR = HERE.parent / "structured-d12-orbit-20261010" / "budget.json"
EXPECTED_OBJECT_SHA256 = "c5de0b26e84ad80bfb7be6a214b89c9682e356effb6079c7abf28cd1151ffb03"
SOURCE = Path("examples/koblitz_n83_structured_orbit_capacity.rs")
PROTOCOL = Path("research/koblitz_n83_factor_base_sweep_20261008/verification/structured-orbit-capacity-20261010/PROTOCOL.md")
RUNNER = Path("research/koblitz_n83_factor_base_sweep_20261008/verification/structured-orbit-capacity-20261010/run_capacity.py")
BUILDER = Path("research/koblitz_n83_factor_base_sweep_20261008/verification/structured-orbit-capacity-20261010/build-linux.sh")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def write_new(name: str, data: object) -> None:
    with (HERE / name).open("x") as stream:
        if isinstance(data, bytes):
            raise TypeError("binary data needs write_bytes_new")
        if isinstance(data, str):
            stream.write(data)
        else:
            json.dump(data, stream, indent=2, sort_keys=True)
            stream.write("\n")


def write_bytes_new(name: str, data: bytes) -> None:
    with (HERE / name).open("xb") as stream:
        stream.write(data)


def docker_call(*args: str, timeout: int = 15) -> subprocess.CompletedProcess[bytes] | None:
    try:
        return subprocess.run(["docker", *args], capture_output=True, timeout=timeout, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return None


def main() -> int:
    if len(sys.argv) != 2:
        raise ValueError("usage: run_capacity.py CONFIG.json")
    config_path = Path(sys.argv[1]).resolve()
    config = json.loads(config_path.read_text())
    if config_path != HERE / "config.json" or config["schema"] != "n83.structured-orbit-capacity-config/v1":
        raise ValueError("frozen config path/schema")
    for name in ("outer.json", "budget.json", "worker.stdout", "worker.stderr"):
        if (HERE / name).exists():
            raise FileExistsError(name)
    if git("rev-parse", "HEAD") != config["source_commit"]:
        raise ValueError("source commit changed")
    if git("diff", "--name-only", "HEAD", "--", str(SOURCE), str(PROTOCOL), str(RUNNER), str(BUILDER)):
        raise ValueError("frozen source/protocol/runner dirty")
    binary = Path(config["binary"]).resolve()
    obj = Path(config["object"]).resolve()
    manifest = (ROOT / config["manifest"]).resolve()
    corpus = (ROOT / config["public_corpus"]).resolve()
    if digest(binary) != config["binary_sha256"] or digest(obj) != EXPECTED_OBJECT_SHA256 or config["object_sha256"] != EXPECTED_OBJECT_SHA256:
        raise ValueError("binary or S3 object hash")
    if digest(manifest) != config["manifest_sha256"] or digest(corpus) != config["public_corpus_sha256"]:
        raise ValueError("manifest or public corpus hash")
    if digest(PRIOR) != config["prior_budget_sha256"]:
        raise ValueError("prior budget hash")
    prior = json.loads(PRIOR.read_text())
    charged = Decimal(prior["cumulative_charged_seconds"])
    authorized = Decimal(prior["authorized_seconds"])
    round_charge = Decimal(prior["remaining_unused_seconds"])
    if charged != Decimal("6922.092994294") or authorized != Decimal("7200") or charged + round_charge != authorized:
        raise ValueError("budget chain")
    if config["wall_cap_seconds"] != 240 or config["memory_limit_bytes"] != 4 * 1024**3:
        raise ValueError("resource contract")
    image = config["image"]
    image_check = docker_call("image", "inspect", image, "--format", "{{.Id}} {{.Architecture}} {{.Os}}")
    if not image_check or image_check.returncode != 0 or image_check.stdout.decode().strip() != f'{config["image_id"]} arm64 linux':
        raise ValueError("local image identity/architecture")
    name = f"n83-structured-capacity-{os.getpid()}"
    command = [
        "docker", "run", "--pull", "never", "--platform", "linux/arm64", "--name", name,
        "--network", "none", "--cpus", "1", "--memory", "4096m", "--memory-swap", "4096m",
        "--pids-limit", "64", "--cap-drop", "ALL", "--security-opt", "no-new-privileges",
        "--read-only", "--tmpfs", "/tmp:rw,noexec,nosuid,size=64m",
        "--user", f"{os.getuid()}:{os.getgid()}",
        "--env", f'N83_SOURCE_COMMIT={config["source_commit"]}',
        "--mount", f"type=bind,src={binary},dst=/worker,readonly",
        "--mount", f"type=bind,src={obj},dst=/object.jsonl.gz,readonly",
        "--mount", f"type=bind,src={manifest},dst=/manifest.json,readonly",
        "--mount", f"type=bind,src={corpus},dst=/probe-corpus.json,readonly",
        image, "/worker", "/object.jsonl.gz", "/manifest.json", "/probe-corpus.json",
    ]
    load_before = os.getloadavg()
    start = time.monotonic()
    process = subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, start_new_session=True)
    timed_out = False
    try:
        client_stdout, client_stderr = process.communicate(timeout=240)
    except subprocess.TimeoutExpired:
        timed_out = True
        docker_call("kill", name)
        os.killpg(process.pid, signal.SIGKILL)
        client_stdout, client_stderr = process.communicate()
    wall = time.monotonic() - start
    state_result = docker_call("inspect", "--format", "{{json .State}}", name)
    state = None
    if state_result and state_result.returncode == 0:
        try:
            state = json.loads(state_result.stdout)
        except (ValueError, UnicodeDecodeError):
            pass
    if not isinstance(state, dict):
        state = None
    logs = docker_call("logs", name)
    stdout = logs.stdout if logs and logs.returncode == 0 else client_stdout
    stderr = logs.stderr if logs and logs.returncode == 0 else client_stderr
    removed = docker_call("rm", "-f", name)
    remaining = docker_call("ps", "-aq", "--filter", f"name=^{name}$")
    cleanup_ok = bool(removed and removed.returncode == 0 and remaining and remaining.returncode == 0 and not remaining.stdout.strip())
    write_bytes_new("worker.stdout", stdout)
    write_bytes_new("worker.stderr", stderr)
    worker = None
    if process.returncode == 0 and not timed_out and state and state.get("ExitCode") == 0:
        try:
            worker = json.loads(stdout)
        except (ValueError, UnicodeDecodeError):
            pass
    if not isinstance(worker, dict):
        worker = None
    valid = bool(cleanup_ok and worker and worker.get("schema") == "n83.structured-orbit-capacity-worker/v1"
        and worker.get("status") == "PASS_model_construction_only"
        and worker.get("source_commit") == config["source_commit"]
        and worker.get("object_compressed_blake3") == config["object_blake3"]
        and worker.get("sat_variables") == 111008
        and worker.get("domain_variables") == 26265
        and worker.get("domain_clauses") == 872260
        and worker.get("domain_xors") == 3320
        and worker.get("memory_cgroup_limit_bytes") == 4 * 1024**3
        and worker.get("memory_cgroup_swap_limit_bytes") == 0
        and type(worker.get("memory_cgroup_peak_bytes")) is int
        and 0 < worker["memory_cgroup_peak_bytes"] <= 4 * 1024**3
        and worker.get("solver_search_executed") is False
        and worker.get("relation_stage_executed") is False
        and worker.get("rank_stage_executed") is False
        and worker.get("total_index_calculus_runtime_ms") is None
        and worker.get("selected_best_total_runtime") is None)
    status = ("PASS_model_construction_only" if valid else
              "PRODUCER_FAILURE_cleanup" if not cleanup_ok else
              "UNKNOWN_wall_cap" if timed_out else
              "UNKNOWN_memory_cap" if state and state.get("OOMKilled") else
              "PRODUCER_FAILURE_worker_exit" if process.returncode else
              "PRODUCER_FAILURE_worker_receipt")
    write_new("outer.json", {
        "schema": "n83.structured-orbit-capacity-outer/v1", "status": status,
        "source_commit": config["source_commit"], "config_sha256": digest(config_path),
        "binary_sha256": config["binary_sha256"], "image_id": config["image_id"],
        "object_sha256": EXPECTED_OBJECT_SHA256, "command": command,
        "docker_client_exit_code": process.returncode, "container_state": state,
        "timed_out": timed_out, "cleanup_verified": cleanup_ok,
        "wall_elapsed_seconds_L0": wall, "load_average_before": load_before,
        "load_average_after": os.getloadavg(), "worker_receipt_valid": valid,
        "worker_stdout_sha256": hashlib.sha256(stdout).hexdigest(),
        "worker_stderr_sha256": hashlib.sha256(stderr).hexdigest(),
        "total_index_calculus_runtime_ms": None, "selected_best_total_runtime": None,
    })
    write_new("budget.json", {
        "schema": "n83.structured-orbit-capacity-budget/v1", "status": "PASS_budget_accounting",
        "prior_budget_sha256": config["prior_budget_sha256"],
        "prior_charged_seconds": str(charged), "round_conservative_charge_seconds": str(round_charge),
        "cumulative_charged_seconds": str(authorized), "authorized_seconds": str(authorized),
        "remaining_unused_seconds": "0", "wall_elapsed_seconds_L0": wall,
        "round_status": status, "selected_best_total_runtime": None,
    })
    print(json.dumps({"status": status, "wall_elapsed_seconds_L0": wall}))
    return 0 if valid else 2


if __name__ == "__main__":
    raise SystemExit(main())
