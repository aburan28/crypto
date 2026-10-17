#!/usr/bin/env python3
"""Replay the retained first hosted N83 Gray/hash artifact independently."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys

import blake3


HERE = Path(__file__).resolve().parent
STUDY = HERE.parent.parent
PAIR = STUDY / "verification/k600-gray-hash-paired-protocol-20261014"
sys.path.insert(0, str(PAIR))
import preflight  # noqa: E402
import receipt_audit as audit  # noqa: E402


RUN_ID = 38983406418
SOURCE_COMMIT = "cb98e85d1c1c1afdd5000bc70eea00969ee86c6c"
WORKFLOW_COMMIT = "8dd40601e80ae8e632229c6b4b633d4f21925513"
ARCHIVES = {
    "3888e2597fd3519807c9d722f49f83329c93ce5cb0929501cf7cb7cbc3cab1a0":
        ("public_x_sequential", 2026100802,
         "8463aaf30365d6f6f5d6b189aff4e8afcb6c8305ec27c216700776cbdccba29d"),
    "966a0493c96674f0af6d184615ed05309d08fba5c7342957b4186cd50366dc47":
        ("public_x_gray_prefix", 2026100801,
         "df02963d1706fa7a64fd7558d88db79b23ea34801ac58062c4c500113e8c57fe"),
    "b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e":
        ("public_x_hash", 2026100801,
         "3e85c9e419fcddd3c9faed3c84470bf96fcad9178e48e1d73629369335841ecf"),
}


def require(condition: bool, reason: str) -> None:
    if not condition:
        raise ValueError(reason)


def read_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"not an object: {path}")
    return value


def digest(path: Path, factory) -> str:
    state = factory()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            state.update(chunk)
    return state.hexdigest()


def plain_blake3(path: Path) -> str:
    state = blake3.blake3()
    with gzip.open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            state.update(chunk)
    return state.hexdigest()


def replay(root: Path) -> dict:
    root = root.resolve(strict=True)
    require(root.name == f"n83-hosted-pair-{RUN_ID}-1", "wrong run artifact")
    pilot = read_json(root / "n83-pilot-parent/panel/manifest.json")
    require(pilot.get("schema") == "n83.factor-base-panel/v1"
            and pilot.get("status") == "completed_factor_base_panel"
            and pilot.get("completed_base_count") == 54
            and pilot.get("expected_base_count") == 54
            and pilot.get("source_commit") == SOURCE_COMMIT,
            "pilot manifest identity or completeness changed")
    rows = {row["compressed_blake3"]: row for row in pilot["bases"]}
    require(len(rows) == len(pilot["bases"]) == 54, "pilot hash collision")
    archive_dir = root / "n83-archives"
    archive_sha = {}
    for object_hash, (policy, seed, expected_sha) in ARCHIVES.items():
        name = f"{object_hash}.jsonl.gz"
        path = archive_dir / name
        row = rows[object_hash]
        require(row.get("a") == 0 and row.get("columns") == 600
                and row.get("policy") == policy and row.get("seed") == seed
                and row.get("object") == f"objects/{name}"
                and path.stat().st_size == row.get("bytes")
                and digest(path, hashlib.sha256) == expected_sha
                and digest(path, blake3.blake3) == object_hash
                and plain_blake3(path) == row.get("plain_blake3"),
                f"frozen archive changed: {name}")
        archive_sha[name] = expected_sha
    prep = read_json(root / "n83-setup/regeneration.json")
    require(prep.get("schema") == "n83.hosted-regenerated-archives/v1"
            and prep.get("status") == "PASS_3_FROZEN_OBJECTS"
            and prep.get("source_commit") == SOURCE_COMMIT
            and prep.get("workflow_dispatch_head_sha") == WORKFLOW_COMMIT
            and prep.get("archive_sha256") == archive_sha
            and prep.get("pilot_manifest_sha256") == digest(
                root / "n83-pilot-parent/panel/manifest.json", hashlib.sha256)
            and prep.get("pilot_completed_base_count") == 54
            and isinstance(prep.get("hosted_preparation_wall_seconds_L0"), (int, float))
            and 0 <= prep["hosted_preparation_wall_seconds_L0"] <= 1800,
            "regeneration receipt does not bind frozen archives")
    outer_dir = root / "n83-hosted-pair"
    config = read_json(outer_dir / "config.json")
    outer = read_json(outer_dir / "outer.json")
    require(config.get("source_commit") == SOURCE_COMMIT
            and config.get("preparation") == prep
            and config.get("selected_best_total_runtime") is None
            and config.get("github_run_id") == str(RUN_ID)
            and config.get("plan_sha256") == digest(preflight.PLAN, hashlib.sha256),
            "hosted config changed")
    input_sha = {Path(name).name: value
                 for name, value in config["input_sha256"].items()}
    require(input_sha == archive_sha, "hosted config input hashes changed")
    require(outer.get("schema") == "n83.hosted-gray-hash-paired-outer/v1"
            and outer.get("status") == "BOUNDED_PARTIAL_hosted_first_trial"
            and outer.get("new_n83_solver_search_executed") is True
            and outer.get("frozen_input_changed") is False
            and outer.get("selected_best_total_runtime") is None
            and outer.get("conservative_charge_seconds") == 4800
            and outer.get("conservative_charge_including_preparation_seconds_L0")
                == 4800 + math.ceil(prep["hosted_preparation_wall_seconds_L0"]),
            "hosted outer disposition or charge changed")
    plan = preflight.make_plan(archive_dir)
    require(len(outer["arms"]) == len(plan["arms"]) == 2, "paired arms missing")
    arm_receipts = []
    for recorded, arm in zip(outer["arms"], plan["arms"], strict=True):
        role = arm["role"]
        arm_dir = outer_dir / role
        require(recorded == read_json(arm_dir / "outer.json")
                and recorded.get("role") == role
                and recorded.get("status") == "UNKNOWN_search_cap"
                and recorded.get("replay_error") is None
                and recorded.get("matrix_rank") is None
                and recorded.get("total_index_calculus_runtime_ms") is None,
                f"outer arm disposition changed: {role}")
        stages = recorded["stages"]
        require([stage["name"] for stage in stages] == ["build", "import", "search"],
                f"stage sequence changed: {role}")
        for stage in stages:
            name = stage["name"]
            require(stage == read_json(arm_dir / f"{name}.stage.json")
                    and stage.get("exit_code") == 0
                    and stage.get("timed_out") is False
                    and stage.get("cleanup_verified") is True
                    and stage.get("launch_error") is None
                    and stage.get("container_state", {}).get("ExitCode") == 0
                    and stage["container_state"].get("OOMKilled") is False,
                    f"stage not clean: {role}/{name}")
            for suffix in ("stdout", "stderr"):
                path = (arm_dir / f"build.{suffix}" if name == "build"
                        else arm_dir / "run" / f"{name}.{suffix}")
                require(digest(path, hashlib.sha256) == stage[f"{suffix}_sha256"],
                        f"raw stage hash changed: {role}/{name}/{suffix}")
        run_dir = arm_dir / "run"
        audit.audit_verify(read_json(run_dir / "import.stdout"), arm)
        base, residual = audit.domains(
            archive_dir / f"{arm['base_object_blake3']}.jsonl.gz",
            archive_dir / f"{arm['pool_object_blake3']}.jsonl.gz",
            arm["shared_signed_orbits"])
        search_config = read_json(run_dir / "search/config.json")
        result = read_json(run_dir / "search/result.json")
        require(audit.audit_search(search_config, result, arm, base, residual)
                == "UNKNOWN_search_cap", f"search replay changed: {role}")
        stats = result["search"]
        require(stats["solver_calls"] == stats["patterns_attempted"] == 1
                and stats["conflicts"] == 5000 and stats["models"] == 0
                and result["candidate"] is None and result["graph_rows"] is None
                and result["matrix_rank"] is None, f"search evidence changed: {role}")
        altered = dict(result)
        altered["status"] = "PASS_one_bounded_partial"
        try:
            audit.audit_search(search_config, altered, arm, base, residual)
        except ValueError:
            pass
        else:
            raise ValueError(f"tampered search disposition accepted: {role}")
        arm_receipts.append({
            "role": role, "status": result["status"],
            "shared_signed_orbits": result["shared_signed_orbits"],
            "residual_signed_orbits": result["residual_signed_orbits"],
            "conflicts": stats["conflicts"], "models": stats["models"],
            "solver_calls": stats["solver_calls"],
            "patterns_attempted": stats["patterns_attempted"],
            "worker_wall_ms_L0": result["worker_wall_ms_L0"],
            "stage_wall_seconds_L0": {s["name"]: s["wall_elapsed_seconds_L0"]
                                      for s in stages},
            "tampered_status_rejected": True,
        })
    files = [{"path": str(path.relative_to(root)),
              "bytes": path.stat().st_size,
              "sha256": digest(path, hashlib.sha256)}
             for path in sorted(root.rglob("*")) if path.is_file()]
    return {
        "schema": "n83.hosted-paired-retained-audit/v1",
        "status": "PASS_FROZEN_BYTES_TWO_CENSORED_SEARCHES",
        "run_id": RUN_ID, "source_commit": SOURCE_COMMIT,
        "workflow_commit": WORKFLOW_COMMIT,
        "artifact_file_count": len(files),
        "artifact_total_bytes": sum(file["bytes"] for file in files),
        "artifact_files": files,
        "archive_sha256": archive_sha,
        "pilot_completed_base_count": 54,
        "preparation_wall_seconds_L0": prep["hosted_preparation_wall_seconds_L0"],
        "conservative_charge_including_preparation_seconds_L0":
            outer["conservative_charge_including_preparation_seconds_L0"],
        "arms": arm_receipts,
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    result = replay(args.artifact)
    with args.receipt.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, sort_keys=True, indent=2)
        stream.write("\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key != "artifact_files"}, sort_keys=True))


if __name__ == "__main__":
    main()
