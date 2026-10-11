#!/usr/bin/env python3
"""Verify archived hashes and stage boundaries for the Linux startup smoke."""

import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
CHECKOUT = HERE.parents[3]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read(name: str) -> bytes:
    return (HERE / name).read_bytes()


def git_bytes(commit: str, path: str) -> bytes:
    return subprocess.run(
        ["git", "-C", str(CHECKOUT), "show", f"{commit}:{path}"],
        capture_output=True, check=True,
    ).stdout


def verify() -> dict:
    config = json.loads(read("config.json"))
    outer = json.loads(read("outer.json"))
    worker = json.loads(read("worker-config.json"))
    source = json.loads(read("source-attestation.json"))
    build = json.loads(read("build-receipt.json"))
    require(build["schema"] == "n83.chain-s3-linux-build/v1" and
            build["status"] == "PASS_build_only", "build receipt shape")
    require(build["build_log_sha256"] == sha256(read("build.log")) and
            b"Finished `release` profile" in read("build.log"), "build log binding")
    require(build["executable_sha256"] == config["binary_sha256"] and
            build["source_commit"] == config["source_commit"], "binary/source binding")
    require("ELF 64-bit" in build["executable_file_description"] and
            "statically linked" in build["executable_file_description"], "Linux ELF shape")
    require(outer["schema"] == "n83.chain-s3-search-outer/v1" and
            outer["status"] == "PRODUCER_FAILURE_worker_exit" and
            outer["worker_exit_code"] == 1, "expected invalid-panel producer failure")
    require(outer["config_sha256"] == sha256(read("config.json")) and
            outer["worker_config_sha256"] == sha256(read("worker-config.json")) and
            outer["stderr_sha256"] == sha256(read("stderr.log")) and
            outer["stdout_sha256"] == sha256(read("stdout.log")), "outer file hashes")
    require(b"panel or replay receipt is not complete and bound" in read("stderr.log"),
            "worker did not reach panel validation")
    require(not any(outer[key] for key in
                    ("source_changed", "config_changed", "binary_changed", "input_changed")),
            "frozen input changed")
    require(config["max_trials"] == worker["max_trials"] == 0 and
            worker["strategy"] == "chain-s3" and
            worker["source_attestation_mode"] == "supervised_snapshot" and
            worker["source_commit"] == source["source_commit"] == config["source_commit"],
            "zero-trial source attestation")
    require(source["schema"] == "n83.chain-s3-source-attestation/v1" and
            source["status_clean"] is True and
            config["source_sha256"]["attestation.json"] == sha256(read("source-attestation.json")),
            "source attestation file")
    require(worker["memory_cgroup_limit_bytes"] == config["memory_cgroup_limit_bytes"] ==
            outer["memory_cgroup_limit_bytes"] == 256 * 1024 * 1024 and
            worker["memory_cgroup_swap_limit_bytes"] ==
            config["memory_cgroup_swap_limit_bytes"] ==
            outer["memory_cgroup_swap_limit_bytes"] == 0, "hard cgroup values")
    for path, digest in config["source_sha256"].items():
        if path != "attestation.json":
            require(sha256(git_bytes(config["source_commit"], path)) == digest,
                    f"source commit differs for {path}")
    for name, digest in config["panel_sha256"].items():
        require(sha256(read(f"empty-panel/{name}")) == digest, f"panel input changed: {name}")
    for name, digest in config["fixture_sha256"].items():
        require(sha256(read(f"fixtures/{name}")) == digest, f"fixture input changed: {name}")
    require(outer["worker_sha256"] is None and outer["report"] is None and
            outer["total_index_calculus_runtime_ms"] is None and
            outer["selected_best_total_runtime"] is None,
            "startup receipt implies a solver result")
    binary = Path(build["executable_path"])
    if binary.is_file():
        require(sha256(binary.read_bytes()) == build["executable_sha256"],
                "local Linux worker differs from build receipt")
    return {
        "status": "PASS_startup_gate_only", "source_commit": config["source_commit"],
        "binary_sha256": build["executable_sha256"], "binary_checked_locally": binary.is_file(),
        "solver_search_executed": False,
        "total_index_calculus_runtime_ms": None,
    }


if __name__ == "__main__":
    print(json.dumps(verify(), sort_keys=True))
