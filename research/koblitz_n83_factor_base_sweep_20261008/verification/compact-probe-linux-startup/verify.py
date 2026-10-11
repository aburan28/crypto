#!/usr/bin/env python3
"""Replay the retained build-only compact-probe startup receipts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys


if sys.flags.optimize:
    raise RuntimeError("startup receipt verification requires enabled assertions")


HERE = Path(__file__).resolve().parent
CHECKOUT = HERE.parents[3]
SOURCE_COMMIT = "2742a09899298fc1a11802a3dbc1e83faa01df20"
IMAGE_ID = "sha256:0dd364ba7e10242f07755449e3a3d0e35f9efd987952737b90def6709ab0c5ce"
BINARY_SHA256 = "2faa775afabdce988dbbf179a721ee746aadf92a3b6ecbede56c96dba9a502d8"


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(value, dict), path
    return value


def checked_source_hashes(config: dict) -> None:
    assert config["source_commit"] == SOURCE_COMMIT
    hashes = config["source_sha256"]
    for relative, expected in hashes.items():
        if relative == "attestation.json":
            content = (HERE / "source-attestation.json").read_bytes()
        else:
            content = subprocess.run(
                ["git", "-C", str(CHECKOUT), "show", f"{SOURCE_COMMIT}:{relative}"],
                capture_output=True, check=True, timeout=15,
            ).stdout
        assert digest(content) == expected, relative


def check_case(name: str, expected_stderr: bytes) -> dict:
    root = HERE / name
    config = read_json(root / "config.json")
    outer = read_json(root / "outer.json")
    assert config["schema"] == "n83.primary-compact-probe-config/v1"
    assert outer["schema"] == "n83.primary-compact-probe-outer/v1"
    assert config["container_image_id"] == IMAGE_ID
    assert config["binary_sha256"] == BINARY_SHA256
    assert config["columns"] == 64 and config["pair_mode"] == "unordered"
    assert config["max_candidate_states"] == 83 * 64 * 65 // 2
    assert config["memory_cgroup_limit_bytes"] == 128 * 1024 * 1024
    assert config["memory_cgroup_swap_limit_bytes"] == 0
    assert config["wall_seconds"] == 30
    checked_source_hashes(config)
    assert digest((root / "config.json").read_bytes()) == outer["config_sha256"]
    assert digest((root / "stdout.log").read_bytes()) == outer["stdout_sha256"]
    assert digest((root / "stderr.log").read_bytes()) == outer["stderr_sha256"]
    assert (root / "stderr.log").read_bytes() == expected_stderr
    assert (root / "stdout.log").read_bytes() == b""
    assert outer["status"] == "PRODUCER_FAILURE_worker_exit"
    assert outer["worker_exit_code"] == 1
    assert outer["worker_sha256"] is None
    assert outer["source_changed"] is False
    assert outer["binary_changed"] is False
    assert outer["config_changed"] is False
    assert outer["total_index_calculus_runtime_ms"] is None
    assert outer["selected_best_total_runtime"] is None
    return {"case": name, "status": outer["status"], "process_wall_ms": outer["process_wall_ms"]}


def main() -> None:
    build = (HERE / "build.log").read_text(encoding="utf-8")
    assert "Finished `release` profile [optimized] target(s) in 3m 36s" in build
    test_summaries = {
        "library-tests.log": "2249 passed; 0 failed; 94 ignored",
        "example-tests.log": "18 passed; 0 failed; 0 ignored",
        "study-python.log": "Ran 15 tests",
        "boundary-python.log": "Ran 16 tests",
    }
    for filename, expected in test_summaries.items():
        outcome = (HERE / filename).read_text(encoding="utf-8")
        assert expected in outcome, filename
        if filename.endswith("-python.log"):
            assert "\nOK\n" in outcome, filename
        else:
            assert "test result: ok." in outcome, filename
    attestation = read_json(HERE / "source-attestation.json")
    assert attestation == {
        "schema": "n83.primary-probe-source-attestation/v1",
        "source_commit": SOURCE_COMMIT,
        "status_clean": True,
    }
    assert (HERE / "invalid-panel" / "manifest.json").read_bytes() == b"{}\n"
    assert (HERE / "invalid-panel" / "replay.json").read_bytes() == b"{}\n"
    cases = [
        check_case("empty-panel", b'Error: Os { code: 2, kind: NotFound, message: "No such file or directory" }\n'),
        check_case("invalid-panel", b'Error: "panel or replay receipt is not complete and bound"\n'),
    ]
    print(json.dumps({"schema": "n83.compact-probe-startup-replay/v1", "status": "PASS_receipts",
                      "source_commit": SOURCE_COMMIT, "binary_sha256": BINARY_SHA256,
                      "build_log_sha256": digest((HERE / "build.log").read_bytes()),
                      "verification_logs_sha256": {
                          name: digest((HERE / name).read_bytes()) for name in test_summaries
                      },
                      "cases": cases}, sort_keys=True))


if __name__ == "__main__":
    main()
