#!/usr/bin/env python3
"""Independent row-level replay of the retained Pareto-shape import panel."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess

import blake3


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PLAN = HERE / "plan.json"
PANEL = HERE / "panel.py"
POINTS_PER_ORBIT = 166


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


def archive(path: Path, spec: dict) -> tuple[list[frozenset], set, set]:
    all_points: set[tuple[str, str]] = set()
    all_x: set[str] = set()
    orbits = []
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        header = json.loads(next(stream))
        require(header["schema"] == "n83.factor-base/v1"
                and header["orbit_columns"] == spec["columns"]
                and header["point_count"] == spec["columns"] * POINTS_PER_ORBIT
                and header["policy"] == spec["policy"]
                and header["seed"] == spec["seed"],
                "archive header and manifest disagree")
        for column in range(spec["columns"]):
            orbit = set()
            for offset in range(POINTS_PER_ORBIT):
                point_row = json.loads(next(stream))
                pair = point_row["point"]
                require(isinstance(pair, list) and len(pair) == 2
                        and all(isinstance(value, str) and value.isdecimal()
                                for value in pair)
                        and point_row["column"] == column
                        and point_row["phase"] == offset // 2
                        and point_row["negative"] is (offset % 2 == 1),
                        "archive point or signed-orbit label changed")
                orbit.add((pair[0], pair[1]))
            require(len(orbit) == POINTS_PER_ORBIT
                    and not (orbit & all_points),
                    "duplicate or partial signed orbit inside archive")
            orbits.append(frozenset(orbit))
            all_points.update(orbit)
            all_x.update(point[0] for point in orbit)
        require(next(stream, None) is None
                and len(all_points) == spec["columns"] * POINTS_PER_ORBIT,
                "archive has trailing rows or missing points")
    return orbits, all_points, all_x


def replay(run_dir: Path, binary: Path) -> dict:
    plan = json.loads(PLAN.read_text())
    result = json.loads((run_dir / "result.json").read_text())
    require(result["status"] == "PASS_NATIVE_IMPORT_PANEL_ONLY"
            and result["error"] is None
            and result["solver_search_executed"] is False
            and result["selected_best_total_runtime"] is None
            and len(result["objects"]) == 18
            and len(result["arms"]) == 24
            and len(result["stages"]) == 60,
            "panel result is incomplete")
    require(result["plan_sha256"] == sha256(PLAN)
            and result["binary_sha256"] == sha256(binary),
            "plan or executable differs from run")
    source = PANEL.relative_to(ROOT).as_posix()
    committed_source = subprocess.run(
        ["git", "show", f"{result['source_commit']}:{source}"], cwd=ROOT,
        check=True, capture_output=True).stdout
    require(hashlib.sha256(committed_source).hexdigest() == result["source_sha256"],
            "run source differs from committed source")
    stages = {}
    for stage in result["stages"]:
        name = stage["name"]
        require(name not in stages and stage["exit_code"] == 0
                and stage["timed_out"] is False, "stage missing, duplicated or incomplete")
        stages[name] = stage
        for stream in ("stdout", "stderr"):
            path = run_dir / "logs" / f"{name}.{stream}"
            require(sha256(path) == stage[f"{stream}_sha256"],
                    f"raw stage stream changed: {path.name}")
    specs = {row["compressed_blake3"]: row for row in plan["objects"]}
    require(len(specs) == 18, "plan object set changed")
    parsed = {}
    for row in result["objects"]:
        digest = row["compressed_blake3"]
        spec = specs[digest]
        path = run_dir / "archives" / f"{digest}.jsonl.gz"
        require(path.stat().st_size == spec["bytes"] == row["bytes"]
                and blake3_file(path) == digest
                and sha256(path) == row["compressed_sha256"]
                and row["s3_uri"] == spec["s3_uri"]
                and row["s3_head"]["ContentLength"] == spec["bytes"],
                "downloaded S3 object changed")
        parsed[digest] = archive(path, spec)
    require(set(parsed) == set(specs), "not every planned S3 object was parsed")
    outcomes = []
    for arm, saved in zip(plan["arms"], result["arms"]):
        require(saved["index"] == arm["index"]
                and saved["base_object_blake3"] == arm["base_object_blake3"]
                and saved["pool_object_blake3"] == arm["pool_object_blake3"],
                "arm order or identity changed")
        base_orbits, base_points, base_x = parsed[arm["base_object_blake3"]]
        pool_orbits, _, _ = parsed[arm["pool_object_blake3"]]
        shared = 0
        for orbit in pool_orbits:
            overlap = len(orbit & base_points)
            require(overlap in (0, POINTS_PER_ORBIT), "partial signed-orbit overlap")
            if overlap:
                shared += 1
            else:
                require(not any(point[0] in base_x for point in orbit),
                        "residual x-domain overlaps base")
        residual = len(pool_orbits) - shared
        require(residual == arm["residual_signed_orbits"]
                and residual == saved["residual_signed_orbits"]
                and shared == saved["shared_signed_orbits"],
                "worker overlap differs from independent archive replay")
        record = json.loads((run_dir / "logs" /
                             f"arm-{arm['index']:02d}.stdout").read_text())
        require(record["status"] == "PASS"
                and record["solver_search_executed"] is False
                and record["small_object_blake3"] == arm["base_object_blake3"]
                and record["pool_object_blake3"] == arm["pool_object_blake3"]
                and record["shared_signed_orbits"] == shared
                and record["residual_signed_orbits"] == residual,
                "worker raw receipt differs from independently replayed domains")
        outcomes.append({"index": arm["index"], "shared_signed_orbits": shared,
                         "residual_signed_orbits": residual,
                         "base_points": len(base_points),
                         "pool_points": sum(map(len, pool_orbits))})
    return {"schema": "n83.pareto-shape-native-import-audit/v1",
            "status": "PASS_S3_BYTES_RAW_STAGES_AND_24_ORBIT_PARTITIONS",
            "source_commit": result["source_commit"],
            "result_sha256": sha256(run_dir / "result.json"),
            "plan_sha256": result["plan_sha256"],
            "binary_sha256": result["binary_sha256"],
            "object_count": len(parsed), "arm_count": len(outcomes),
            "outcomes": outcomes, "solver_search_executed": False,
            "selected_best_total_runtime": None}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), "audit output must be fresh")
    value = replay(args.run_dir, args.binary)
    args.output.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")
    print(value["status"])
