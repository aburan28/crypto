#!/usr/bin/env python3
"""Bounded, import-only native check of the 24 N83 source-Pareto shapes."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time
from urllib.parse import urlparse

import blake3


HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
ROOT = HERE.parents[3]
PARETO = HERE.parent / "pilot-pair-pareto-20261014" / "receipt.json"
MANIFEST = STUDY / "pilot-01" / "manifest.json"
PLAN = HERE / "plan.json"
PINS = {
    PARETO: "b2053f7a8a21cd04831a05c2cc9c4443cc7d582c1b53360bec3a5064fa72e88e",
    MANIFEST: "5bd862987844f709e03beeab9fd265e3d2024011bcadcb82f5af9c72018df6a6",
    ROOT / "examples/koblitz_n83_mixed_orbit_worker.rs":
        "5f310f7495989a94004a6545bed5682d54d3bc4b4f41db543ed24408a1fa1b7b",
}
BINARY_SHA256 = "c5b9af9930f293533de819c0a673a0b7149ad2418800fa3fee7509d772d207e7"
SRC_TREE = "f5dac3ceb6d264227ef74b197e49f3cf01d1164e"
BUCKET = "crypto-autoresearcher"
PREFIX = "factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/a0/objects/"
PANEL_CAP = 1200


class BoundedStop(Exception):
    pass


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def blake3_file(path: Path) -> str:
    digest = blake3.blake3()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def s3_key(uri: str, digest: str) -> str:
    parsed = urlparse(uri)
    key = parsed.path.lstrip("/")
    require(parsed.scheme == "s3" and parsed.netloc == BUCKET
            and not parsed.query and not parsed.fragment
            and key == f"{PREFIX}{digest}.jsonl.gz", "unexpected S3 object URI")
    return key


def make_plan() -> dict:
    for path, expected in PINS.items():
        require(sha256(path) == expected, f"pinned source changed: {path.name}")
    pareto = json.loads(PARETO.read_text())
    manifest = json.loads(MANIFEST.read_text())
    require(pareto["status"] == "PASS_SOURCE_ONLY_NO_RUNTIME_WINNER"
            and pareto["source_pareto_unique_shape_count"] == 24
            and pareto["source_pareto_pair_count"] == 311
            and pareto["selected_best_total_runtime"] is None,
            "source Pareto receipt changed")
    rows = {row["compressed_blake3"]: row for row in manifest["bases"]}
    require(len(rows) == len(manifest["bases"]), "duplicate manifest object")
    arms = []
    for index, shape in enumerate(pareto["frontier_shapes"]):
        base, pool = shape["example_base_object_blake3"], shape["example_pool_object_blake3"]
        require(base != pool and base in rows and pool in rows, "frontier object missing")
        require(rows[base]["columns"] == shape["example_base_orbits"]
                and rows[pool]["columns"] >= shape["example_residual_signed_orbits"] >= 2,
                "frontier shape and manifest disagree")
        arms.append({"index": index, "base_object_blake3": base,
                     "pool_object_blake3": pool,
                     "base_columns": rows[base]["columns"],
                     "pool_columns": rows[pool]["columns"],
                     "residual_signed_orbits": shape["example_residual_signed_orbits"],
                     "conditional_half_K_row_query_floor":
                         shape["conditional_half_K_row_query_floor"],
                     "source_domain_clauses": shape["source_domain_clauses"],
                     "source_total_variables": shape["source_total_variables"],
                     "tied_pair_count": shape["tied_pair_count"]})
    identities = sorted({arm[key] for arm in arms
                         for key in ("base_object_blake3", "pool_object_blake3")})
    require(len(arms) == 24 and len(identities) == 18,
            "representative panel population changed")
    objects = []
    for digest in identities:
        row = rows[digest]
        s3_key(row["s3_uri"], digest)
        objects.append({key: row[key] for key in (
            "compressed_blake3", "s3_uri", "bytes", "columns", "policy", "seed",
            "plain_blake3", "point_set_blake3")})
    require(sum(row["bytes"] for row in objects) == 27822670,
            "compressed panel size changed")
    return {"schema": "n83.pareto-shape-native-import-plan/v1",
            "status": "SOURCE_READY_IMPORT_ONLY", "parent_pareto_sha256": PINS[PARETO],
            "manifest_sha256": PINS[MANIFEST], "worker_source_sha256":
                PINS[ROOT / "examples/koblitz_n83_mixed_orbit_worker.rs"],
            "rust_src_git_tree_oid": SRC_TREE, "binary_sha256": BINARY_SHA256,
            "panel_wall_cap_seconds": PANEL_CAP, "object_count": len(objects),
            "arm_count": len(arms), "objects": objects, "arms": arms,
            "solver_search_executed": False, "selected_best_total_runtime": None}


def run_stage(name: str, argv: list[str], output: Path, cap: int,
              deadline: float, env: dict[str, str]) -> dict:
    remaining = int(deadline - time.monotonic())
    if remaining < 1:
        raise BoundedStop("panel wall cap reached before stage")
    cap = min(cap, remaining)
    stdout_path = output / "logs" / f"{name}.stdout"
    stderr_path = output / "logs" / f"{name}.stderr"
    started = time.monotonic()
    timed_out = False
    with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
        child = subprocess.Popen(["/usr/bin/time", "-l", *argv], cwd=ROOT,
                                 env=env, stdout=stdout, stderr=stderr,
                                 start_new_session=True)
        try:
            exit_code = child.wait(timeout=cap)
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
    return {"name": name, "argv": argv, "cap_seconds": cap,
            "elapsed_wall_seconds_L0": time.monotonic() - started,
            "exit_code": exit_code, "timed_out": timed_out,
            "maximum_resident_set_size_bytes_L0": int(rss.group(1)) if rss else None,
            "stdout_sha256": sha256(stdout_path), "stderr_sha256": sha256(stderr_path)}


def check_stage(stage: dict) -> None:
    if stage["timed_out"]:
        raise BoundedStop(f"{stage['name']} reached its wall cap")
    require(stage["exit_code"] == 0, f"{stage['name']} exited {stage['exit_code']}")


def verify_record(record: dict, arm: dict, rows: dict) -> None:
    base, pool = rows[arm["base_object_blake3"]], rows[arm["pool_object_blake3"]]
    shared = pool["columns"] - arm["residual_signed_orbits"]
    expected = {"schema": "n83.cross-stream-pilot-import/v1", "status": "PASS",
                "curve_slug": "icv1-f2m83-tm6151469093347-debefd74",
                "small_object_blake3": base["compressed_blake3"],
                "small_policy": base["policy"], "small_seed": base["seed"],
                "small_columns": base["columns"], "small_points": base["columns"] * 166,
                "pool_object_blake3": pool["compressed_blake3"],
                "pool_policy": pool["policy"], "pool_seed": pool["seed"],
                "pool_columns": pool["columns"], "pool_points": pool["columns"] * 166,
                "shared_signed_orbits": shared,
                "small_only_signed_orbits": base["columns"] - shared,
                "residual_signed_orbits": arm["residual_signed_orbits"],
                "residual_points": arm["residual_signed_orbits"] * 166,
                "solver_search_executed": False}
    require(shared >= 0 and expected["small_only_signed_orbits"] >= 0,
            "impossible signed-orbit overlap")
    for key, value in expected.items():
        require(type(record.get(key)) is type(value) and record[key] == value,
                f"worker receipt {key} changed")


def run(binary: Path, output: Path) -> dict:
    output = output.resolve()
    require(not output.exists() and not output.is_relative_to(ROOT),
            "output must be fresh and outside checkout")
    plan = make_plan()
    require(PLAN.read_bytes() == (json.dumps(plan, sort_keys=True, indent=2) + "\n").encode(),
            "frozen panel plan changed")
    plan_sha = sha256(PLAN)
    require(sha256(binary) == BINARY_SHA256, "preserved worker binary changed")
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
                            capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=all"],
                           cwd=ROOT, check=True, capture_output=True,
                           text=True).stdout.strip()
    tree = subprocess.run(["git", "rev-parse", "HEAD:src"], cwd=ROOT, check=True,
                          capture_output=True, text=True).stdout.strip()
    require(not dirty and tree == SRC_TREE, "checkout or pinned Rust source changed")
    output.mkdir(parents=True)
    (output / "archives").mkdir()
    (output / "logs").mkdir()
    env = dict(os.environ, AWS_PAGER="")
    deadline = time.monotonic() + PANEL_CAP
    stages, objects, arms = [], [], []
    status, error = "PRODUCER_FAILURE_IMPORT_PANEL", None
    rows = {row["compressed_blake3"]: row for row in plan["objects"]}
    try:
        for row in plan["objects"]:
            digest = row["compressed_blake3"]
            key = s3_key(row["s3_uri"], digest)
            head = run_stage(f"head-{digest[:12]}",
                             ["aws", "s3api", "head-object", "--bucket", BUCKET,
                              "--key", key, "--output", "json"], output, 30, deadline, env)
            stages.append(head)
            check_stage(head)
            metadata = json.loads((output / "logs" / f"{head['name']}.stdout").read_text())
            require(metadata["ContentLength"] == row["bytes"], "S3 object size changed")
            path = output / "archives" / f"{digest}.jsonl.gz"
            got = run_stage(f"get-{digest[:12]}",
                            ["aws", "s3api", "get-object", "--bucket", BUCKET,
                             "--key", key, str(path), "--output", "json"],
                            output, 60, deadline, env)
            stages.append(got)
            check_stage(got)
            require(path.stat().st_size == row["bytes"] and blake3_file(path) == digest,
                    "download differs from content-addressed manifest")
            objects.append({"compressed_blake3": digest, "s3_uri": row["s3_uri"],
                            "bytes": row["bytes"], "compressed_sha256": sha256(path),
                            "s3_head": {key: metadata.get(key) for key in (
                                "ContentLength", "ETag", "LastModified", "VersionId")}})
        for arm in plan["arms"]:
            require(sha256(binary) == BINARY_SHA256 and sha256(PLAN) == plan_sha
                    and all(sha256(path) == digest for path, digest in PINS.items()),
                    "frozen source or binary changed between arms")
            base, pool = arm["base_object_blake3"], arm["pool_object_blake3"]
            argv = [str(binary), "verify-cross-pilot", base, str(MANIFEST),
                    str(output / "archives" / f"{base}.jsonl.gz"), pool,
                    str(MANIFEST), str(output / "archives" / f"{pool}.jsonl.gz")]
            stage = run_stage(f"arm-{arm['index']:02d}", argv, output, 75, deadline, env)
            stages.append(stage)
            check_stage(stage)
            record = json.loads((output / "logs" / f"{stage['name']}.stdout").read_text())
            verify_record(record, arm, rows)
            arms.append({"index": arm["index"], "base_object_blake3": base,
                         "pool_object_blake3": pool,
                         "residual_signed_orbits": record["residual_signed_orbits"],
                         "shared_signed_orbits": record["shared_signed_orbits"],
                         "import_wall_seconds_L0": stage["elapsed_wall_seconds_L0"],
                         "maximum_resident_set_size_bytes_L0":
                             stage["maximum_resident_set_size_bytes_L0"]})
        status = "PASS_NATIVE_IMPORT_PANEL_ONLY"
    except (BoundedStop, OSError, ValueError, KeyError, TypeError,
            json.JSONDecodeError) as exc:
        error = f"{type(exc).__name__}: {exc}"
        status = ("BOUNDED_UNKNOWN_PANEL_WALL_CAP" if isinstance(exc, BoundedStop)
                  else "PRODUCER_FAILURE_IMPORT_PANEL")
    result = {"schema": "n83.pareto-shape-native-import-panel/v1", "status": status,
              "error": error, "source_commit": commit, "source_sha256": sha256(Path(__file__)),
              "plan_sha256": plan_sha, "binary_sha256": BINARY_SHA256,
              "panel_wall_cap_seconds": PANEL_CAP, "stages": stages,
              "objects": objects, "arms": arms,
              "solver_search_executed": False, "relation_rows": None,
              "selected_best_total_runtime": None}
    (output / "result.json").write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-plan", action="store_true")
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    if args.write_plan:
        require(not PLAN.exists(), "plan already exists")
        PLAN.write_text(json.dumps(make_plan(), sort_keys=True, indent=2) + "\n")
        print(sha256(PLAN))
        return 0
    require(args.binary is not None and args.output_dir is not None,
            "binary and output directory required")
    result = run(args.binary.resolve(strict=True), args.output_dir)
    print(json.dumps({key: result[key] for key in ("status", "error")}, sort_keys=True))
    return 0 if result["status"] == "PASS_NATIVE_IMPORT_PANEL_ONLY" else 1


if __name__ == "__main__":
    raise SystemExit(main())
