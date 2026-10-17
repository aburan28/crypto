#!/usr/bin/env python3
"""Fail-closed replay of one hosted N83 fixed-witness control artifact."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DOMAIN_SHA256 = "63accc5b9dc9de6327c5b07e93b3e62717f8376a96c3f91f7bb70b77095f23bb"
IMAGE_DIGEST = "sha256:93ce27a88655056a51dbdd8f5f2d7ddc071c7b0070fb288a37b5a285fc83971e"
INPUTS = {
    "examples/koblitz_n83_fixed_witness_control.rs",
    "research/koblitz_n83_factor_base_sweep_20261008/verification/known-witness-ordered-sat-20261014/domains.json",
    "research/koblitz_n83_factor_base_sweep_20261008/verification/mixed-orbit-known-witness-source-20261012/fixture.json",
    "Cargo.lock",
    ".github/workflows/n83-hosted-fixed-witness.yml",
}
PHASES = ("build_started", "circuit_built", "domains_built", "fixed_units_installed")


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_file(commit: str, path: str, root: Path) -> bytes:
    return subprocess.check_output(["git", "show", f"{commit}:{path}"], cwd=root)


def read_hash_lines(path: Path) -> dict[str, str]:
    records = {}
    for line in path.read_text().splitlines():
        digest, name = line.split("  ", 1)
        require(len(digest) == 64 and all(c in "0123456789abcdef" for c in digest),
                f"invalid SHA-256 in {path.name}")
        require(name not in records, f"duplicate hash path {name}")
        records[name] = digest
    return records


def replay(artifact: Path, expected_commit: str, root: Path = ROOT) -> dict:
    artifact = artifact.resolve(strict=True)
    source_commit = (artifact / "source-commit.txt").read_text().strip()
    require(source_commit == expected_commit and len(source_commit) == 40,
            "source commit differs from expected dispatch")
    require((artifact / "source-status.txt").read_bytes() == b"", "dirty hosted checkout")
    source_hashes = read_hash_lines(artifact / "inputs.sha256")
    require(set(source_hashes) == INPUTS, "source/input hash path set differs")
    for path, digest in source_hashes.items():
        require(hashlib.sha256(git_file(source_commit, path, root)).hexdigest() == digest,
                f"Git source/input differs: {path}")
    require(source_hashes[next(p for p in INPUTS if p.endswith("/domains.json"))]
            == DOMAIN_SHA256, "fixed witness domain hash differs")
    image_identity = (artifact / "image-identity.txt").read_text().strip()
    require(" arm64 linux " in image_identity and IMAGE_DIGEST in image_identity,
            "image identity differs")
    worker_hashes = read_hash_lines(artifact / "worker.sha256")
    require(len(worker_hashes) == 1 and next(iter(worker_hashes.values()))
            == sha256(artifact / "worker"), "built worker hash differs")
    require((artifact / "worker").read_bytes()[:4] == b"\x7fELF", "worker is not ELF")
    output_hashes = read_hash_lines(artifact / "output.sha256")
    require({Path(p).name for p in output_hashes}
            == {"solve.stdout", "solve.stderr", "container-inspect.json"},
            "raw output hash path set differs")
    for path, digest in output_hashes.items():
        require(sha256(artifact / Path(path).name) == digest,
                f"raw output hash differs: {path}")
    bounds = dict(line.split("=", 1) for line in
                  (artifact / "solve.bounds").read_text().splitlines())
    require(bounds["conflict_cap"] == "5000" and bounds["wall_cap_seconds"] == "1200",
            "bounded solve contract differs")
    elapsed = int(bounds["elapsed_wall_seconds"])
    require(0 <= elapsed <= 1210, "solve wall exceeds external cap")
    exit_code = int((artifact / "solve.exit-code").read_text().strip())
    require(exit_code in (0, 124), "unexpected worker exit")
    inspect = json.loads((artifact / "container-inspect.json").read_text())
    require(isinstance(inspect, list) and len(inspect) == 1,
            "missing or ambiguous container inspection")
    container = inspect[0]
    host = container["HostConfig"]
    require(host["Memory"] == host["MemorySwap"] == 4 * 1024**3
            and host["NanoCpus"] == 1_000_000_000 and host["PidsLimit"] == 64
            and host["NetworkMode"] == "none" and host["ReadonlyRootfs"] is True,
            "resource isolation differs")
    require(container["Config"]["Image"].endswith(IMAGE_DIGEST),
            "runtime image differs")
    require((artifact / "container-cleanup.stdout").read_text().strip()
            == container["Name"].lstrip("/"), "container cleanup unverified")
    markers = []
    for line in (artifact / "solve.stderr").read_text().splitlines():
        if line.startswith("{"):
            marker = json.loads(line)
            if "phase" in marker:
                markers.append(marker)
    require([m["phase"] for m in markers] == list(PHASES[:len(markers)]),
            "phase marker order differs")
    require(all(isinstance(m.get("elapsed_wall_ms_L0"), (int, float))
                and 0 <= m["elapsed_wall_ms_L0"] <= 1_210_000 for m in markers),
            "phase wall differs")
    require(all(a["elapsed_wall_ms_L0"] <= b["elapsed_wall_ms_L0"]
                for a, b in zip(markers, markers[1:])), "phase walls decreased")
    raw_stdout = (artifact / "solve.stdout").read_text().strip()
    if exit_code == 0:
        result = json.loads(raw_stdout)
        require(len(markers) == len(PHASES)
                and result["schema"] == "n83.fixed-witness-sat-control/v1"
                and result["conflict_cap"] == 5000
                and 0 <= result["conflicts"] <= 5000
                and result["fixed_x_blocks"] == 8
                and result["input_blake3"]
                and result["natural_relations"] is None
                and result["total_index_calculus_runtime_ms"] is None
                and result["selected_best_total_runtime"] is None,
                "completed worker receipt differs")
        status = result["status"]
        require(status in {"SAT_verified_witness_x", "UNKNOWN_conflict_cap",
                           "UNSAT_fixed_witness_encoding",
                           "PRODUCER_FAILURE_invalid_model"}, "worker status differs")
        require((status == "SAT_verified_witness_x") ==
                (result["model_verified"] is True), "model verification differs")
    else:
        require(raw_stdout == "", "timed-out worker emitted a completed receipt")
        if container["State"]["OOMKilled"]:
            status = "UNKNOWN_resource_cap"
        elif markers and markers[-1]["phase"] == "fixed_units_installed":
            status = "UNKNOWN_wall_cap_after_fixed_units"
        else:
            status = "UNKNOWN_wall_cap_before_search"
        result = None
    return {
        "schema": "n83.hosted-fixed-witness-audit/v1",
        "status": status,
        "source_commit": source_commit,
        "worker_sha256": sha256(artifact / "worker"),
        "artifact_file_sha256": {p.name: sha256(p) for p in sorted(artifact.iterdir())
                                 if p.is_file()},
        "exit_code": exit_code,
        "elapsed_wall_seconds_L0": elapsed,
        "phases": markers,
        "solver_conflicts": None if result is None else result["conflicts"],
        "model_verified": None if result is None else result["model_verified"],
        "natural_relations": None,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--expected-commit", required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = replay(args.artifact, args.expected_commit)
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        with args.output.open("x") as stream:
            stream.write(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
