#!/usr/bin/env python3
"""Bounded native import validation of the two frozen K600 S3 pairs; no SAT."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parent / "k600-gray-hash-paired-protocol-20261014"))
import preflight  # noqa: E402
import receipt_audit as audit  # noqa: E402
import supervisor  # noqa: E402


TOOLCHAIN = Path("/Users/adamburan/.rustup/toolchains/1.98-aarch64-apple-darwin/bin")
BUILD_CAP = 1200
IMPORT_CAP = 360
WORKER = "koblitz_n83_mixed_orbit_worker"


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fresh(path: Path) -> Path:
    path = path.resolve()
    require(not path.is_relative_to(ROOT) and not path.exists(),
            "output must be fresh and outside checkout")
    return path


def verify_argv(binary: Path, arm: dict, archives: Path) -> list[str]:
    frozen = arm["worker_verify_argv"]
    require(len(frozen) == 8 and frozen[1] == "verify-cross-pilot"
            and frozen[2] == arm["base_object_blake3"]
            and frozen[5] == arm["pool_object_blake3"],
            "frozen verify command changed")
    manifest = str(preflight.PILOT / "manifest.json")
    return [str(binary), frozen[1], frozen[2], manifest,
            str(archives / f"{frozen[2]}.jsonl.gz"), frozen[5], manifest,
            str(archives / f"{frozen[5]}.jsonl.gz")]


def run_stage(name: str, argv: list[str], output: Path, cap_seconds: int,
              env: dict[str, str]) -> dict:
    stdout_path, stderr_path = output / f"{name}.stdout", output / f"{name}.stderr"
    started = time.monotonic()
    timed_out = False
    with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
        child = subprocess.Popen(["/usr/bin/time", "-l", *argv], cwd=ROOT,
                                 env=env, stdout=stdout, stderr=stderr,
                                 start_new_session=True)
        try:
            exit_code = child.wait(timeout=cap_seconds)
        except subprocess.TimeoutExpired:
            timed_out = True
            try:
                os.killpg(child.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                exit_code = child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                exit_code = child.wait()
    stderr_text = stderr_path.read_text(errors="replace")
    rss = re.search(r"^\s*(\d+)\s+maximum resident set size\s*$",
                    stderr_text, re.MULTILINE)
    return {
        "name": name, "argv": argv, "cap_seconds": cap_seconds,
        "exit_code": exit_code, "timed_out": timed_out,
        "elapsed_wall_seconds_L0": time.monotonic() - started,
        "maximum_resident_set_size_bytes_L0": int(rss.group(1)) if rss else None,
        "stdout_sha256": sha256(stdout_path), "stderr_sha256": sha256(stderr_path),
    }


def run(archives: Path, output: Path) -> dict:
    archives = archives.resolve(strict=True)
    output = fresh(output)
    plan = preflight.make_plan(archives)
    commit = audit.PRIOR.clean_commit(ROOT)
    sources = supervisor.source_hashes(ROOT)
    source_sha256 = sha256(Path(__file__))
    cargo, rustc = TOOLCHAIN / "cargo", TOOLCHAIN / "rustc"
    require(cargo.is_file() and rustc.is_file(), "pinned native toolchain unavailable")
    env = dict(os.environ)
    env["PATH"] = str(TOOLCHAIN) + os.pathsep + env.get("PATH", "")
    env["CARGO_TARGET_DIR"] = str(output / "target")
    env["TMPDIR"] = str(output)
    versions = {}
    for label, path in (("cargo", cargo), ("rustc", rustc)):
        versions[label] = subprocess.run([str(path), "--version"], cwd=ROOT,
                                         capture_output=True, text=True, check=True,
                                         timeout=10).stdout.strip()
    require(versions["cargo"].startswith("cargo 1.98.0 ")
            and versions["rustc"].startswith("rustc 1.98.0 "),
            "native toolchain version changed")
    output.mkdir(parents=True, exist_ok=False)
    plan_sha256 = sha256(preflight.PLAN)
    input_sha256 = {str(archives / f"{digest}.jsonl.gz"):
                    sha256(archives / f"{digest}.jsonl.gz")
                    for digest in (preflight.BASE, *preflight.POOLS.values())}
    stages = []
    status = "PRODUCER_FAILURE_build"
    replay_error = None
    binary_sha256 = None
    try:
        build = run_stage("build", [str(cargo), "build", "--release", "--offline",
                                    "--locked", "--no-default-features", "--example", WORKER],
                          output, BUILD_CAP, env)
        stages.append(build)
        require(build["exit_code"] == 0 and not build["timed_out"],
                "native release build failed or timed out")
        binary = output / "target/release/examples" / WORKER
        require(binary.is_file() and os.access(binary, os.X_OK),
                "native worker binary missing")
        binary_sha256 = sha256(binary)
        for arm in plan["arms"]:
            require(audit.PRIOR.clean_commit(ROOT) == commit
                    and supervisor.source_hashes(ROOT) == sources
                    and sha256(Path(__file__)) == source_sha256
                    and sha256(binary) == binary_sha256
                    and all(sha256(Path(path)) == digest
                            for path, digest in input_sha256.items()),
                    "frozen source, archive or binary changed")
            argv = verify_argv(binary, arm, archives)
            receipt = run_stage(arm["role"], argv, output, IMPORT_CAP, env)
            stages.append(receipt)
            require(receipt["exit_code"] == 0 and not receipt["timed_out"],
                    f"{arm['role']} native import failed or timed out")
            record = json.loads((output / f"{arm['role']}.stdout").read_text())
            audit.audit_verify(record, arm)
            audit.domains(archives / f"{arm['base_object_blake3']}.jsonl.gz",
                          archives / f"{arm['pool_object_blake3']}.jsonl.gz",
                          arm["shared_signed_orbits"])
        status = "PASS_NATIVE_IMPORT_ONLY"
    except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as error:
        replay_error = f"{type(error).__name__}: {error}"
        status = ("BOUNDED_UNKNOWN_native_import_wall_cap"
                  if any(stage["timed_out"] for stage in stages)
                  else "PRODUCER_FAILURE_native_import")
    result = {
        "schema": "n83.k600-gray-hash-native-import/v1", "status": status,
        "source_commit": commit, "source_sha256": sources,
        "native_import_source_sha256": source_sha256,
        "plan_sha256": plan_sha256, "input_sha256": input_sha256,
        "toolchain_versions": versions, "binary_sha256": binary_sha256,
        "stages": stages, "replay_error": replay_error,
        "solver_search_executed": False,
        "graph_rows": None, "matrix_rank": None,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }
    (output / "result.json").write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.archive_dir, args.output_dir)
    print(json.dumps({"status": result["status"], "stage_count": len(result["stages"]),
                      "replay_error": result["replay_error"]}, sort_keys=True))
    return 0 if result["status"] == "PASS_NATIVE_IMPORT_ONLY" else 1


if __name__ == "__main__":
    raise SystemExit(main())
