#!/usr/bin/env python3
"""Replay retained native import receipts against frozen Git, binary and archives."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parent / "k600-gray-hash-paired-protocol-20261014"))
import preflight  # noqa: E402
import receipt_audit as audit  # noqa: E402


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def git_bytes(commit: str, path: str) -> bytes:
    return subprocess.run(["git", "show", f"{commit}:{path}"], cwd=ROOT,
                          check=True, capture_output=True).stdout


def replay(archives: Path, evidence: Path, binary: Path) -> dict:
    result = json.loads((evidence / "result.json").read_text())
    plan = json.loads(preflight.PLAN.read_text())
    require(result["status"] == "PASS_NATIVE_IMPORT_ONLY"
            and result["solver_search_executed"] is False
            and result["selected_best_total_runtime"] is None,
            "native result status or claim changed")
    require(result["plan_sha256"] == sha256_file(preflight.PLAN), "plan changed")
    commit = result["source_commit"]
    source = (HERE / "native_import.py").relative_to(ROOT).as_posix()
    require(result["native_import_source_sha256"] == sha256_bytes(git_bytes(commit, source)),
            "native runner source changed")
    for path, digest in result["source_sha256"].items():
        if path == "src:git-tree":
            actual = subprocess.run(["git", "rev-parse", f"{commit}:src"],
                                    cwd=ROOT, check=True, capture_output=True,
                                    text=True).stdout.strip()
        else:
            actual = sha256_bytes(git_bytes(commit, path))
        require(actual == digest, f"frozen Git source changed: {path}")
    require(sha256_file(binary) == result["binary_sha256"], "worker binary changed")
    for path, digest in result["input_sha256"].items():
        actual = archives / Path(path).name
        require(sha256_file(actual) == digest, f"archive changed: {actual.name}")
    require([stage["name"] for stage in result["stages"]]
            == ["build", "gray_pool", "hash_pool"], "stage set changed")
    for stage in result["stages"]:
        require(stage["exit_code"] == 0 and stage["timed_out"] is False,
                f"stage incomplete: {stage['name']}")
        for stream in ("stdout", "stderr"):
            path = evidence / f"{stage['name']}.{stream}"
            require(sha256_file(path) == stage[f"{stream}_sha256"],
                    f"raw stream changed: {path.name}")
    arms = plan["arms"]
    require([arm["role"] for arm in arms] == ["gray_pool", "hash_pool"],
            "plan arms changed")
    domains = {}
    for arm in arms:
        role = arm["role"]
        record = json.loads((evidence / f"{role}.stdout").read_text())
        audit.audit_verify(record, arm)
        base = archives / f"{arm['base_object_blake3']}.jsonl.gz"
        pool = archives / f"{arm['pool_object_blake3']}.jsonl.gz"
        base_points, residual = audit.domains(base, pool, arm["shared_signed_orbits"])
        domains[role] = {"base_points": len(base_points),
                         "residual_points": len(residual),
                         "shared_signed_orbits": arm["shared_signed_orbits"]}
    return {"schema": "n83.k600-native-import-independent-replay/v1",
            "status": "PASS_FROZEN_SOURCE_BINARY_ARCHIVES_AND_DOMAINS",
            "source_commit": commit, "plan_sha256": result["plan_sha256"],
            "binary_sha256": result["binary_sha256"], "domains": domains,
            "solver_search_executed": False,
            "selected_best_total_runtime": None}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-dir", type=Path, required=True)
    parser.add_argument("--evidence-dir", type=Path, default=HERE)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), "replay output must be fresh")
    receipt = replay(args.archive_dir, args.evidence_dir, args.binary)
    args.output.write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(receipt["status"])
