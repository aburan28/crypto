#!/usr/bin/env python3
"""Independently replay the hosted fixed N83 witness and S3 artifact bytes."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tarfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
STUDY = HERE.parents[1]
SOURCE_COMMIT = "0891c2e062ff95e923302e570a8101f83a59b6a5"
RUN_ID = 39051579293
PREFIX = f"n83-hosted-witness-{RUN_ID}-1/"
DOMAIN_SHA256 = "63accc5b9dc9de6327c5b07e93b3e62717f8376a96c3f91f7bb70b77095f23bb"
FIXTURE_SHA256 = "fc203db28eed0164e47b0c755c3d9a09f80df1acd5e4c7e09dfe172fa5f0a538"
INPUT_BLAKE3 = "22eb1821dd401063417260041017c43cd416a8f4fc1abea2e9946507029c6d4a"
PROTO = STUDY / "verification/hosted-fixed-witness-20261017/audit.py"
sys.path.insert(0, str(STUDY))
from chain_search_supervisor import field_mul, parse_point, point_add  # noqa: E402


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git_blob(path: str) -> bytes:
    return subprocess.check_output(["git", "show", f"{SOURCE_COMMIT}:{path}"], cwd=ROOT)


def algebra_replay() -> dict:
    fixture_path = "research/koblitz_n83_factor_base_sweep_20261008/verification/mixed-orbit-known-witness-source-20261012/fixture.json"
    domain_path = "research/koblitz_n83_factor_base_sweep_20261008/verification/known-witness-ordered-sat-20261014/domains.json"
    reference_path = "research/koblitz_n83_factor_base_sweep_20261008/chain_search_supervisor.py"
    fixture_bytes = git_blob(fixture_path)
    domain_bytes = git_blob(domain_path)
    require(digest(fixture_bytes) == FIXTURE_SHA256
            and digest(domain_bytes) == DOMAIN_SHA256,
            "frozen witness input differs")
    require(digest((ROOT / reference_path).read_bytes()) == digest(git_blob(reference_path)),
            "reference arithmetic source differs from hosted source")
    fixture = json.loads(fixture_bytes)
    domains = json.loads(domain_bytes)
    require(domains["schema"] == "n83.fixed-witness-domains/v1"
            and domains["identity_mask"] == 0
            and domains["fixture_sha256"] == FIXTURE_SHA256
            and domains["ordered_indices"] == [166, 0, 332, 196212, 196378]
            and len(domains["small_reps"]) == 1182
            and len(domains["residual_reps"]) == 866,
            "fixed domain binding differs")
    terms = [fixture["terms"][i] for i in (1, 0, 2, 3, 4)]
    points = [parse_point(term["point"]) for term in terms]
    require(all(point is not None for point in points), "fixed witness has identity term")
    xs = [str(point[0]) for point in points]
    require(xs == domains["ordered_xs"] and int(xs[0]) < int(xs[1]),
            "ordered witness x differs")
    require(all(x in domains["small_reps"] for x in xs[:3])
            and all(x in domains["residual_reps"] for x in xs[3:]),
            "witness x lies outside exact domains")
    partials = []
    acc = None
    for point in points:
        acc = point_add(acc, point)
        require(acc is not None, "fixed witness has intermediate identity")
        partials.append(acc)
    require([str(point[0]) for point in partials[1:-1]] == domains["partial_xs"]
            and [str(v) for v in partials[-1]] == fixture["equation"]["target"]
            and str(partials[-1][0]) == domains["target_x"],
            "fixed witness partial or target differs")
    residues = []
    for step in range(1, len(points)):
        left = points[0] if step == 1 else partials[step - 1]
        x, y, z = left[0], points[step][0], partials[step][0]
        xy = field_mul(x, y)
        pair_sum = xy ^ field_mul(x, z) ^ field_mul(y, z)
        residues.append(field_mul(xy, z) ^ field_mul(pair_sum, pair_sum) ^ 1)
    require(residues == [0, 0, 0, 0], "original chained S3 witness differs")
    return {"ordered_indices": domains["ordered_indices"],
            "ordered_xs": xs, "partial_xs": domains["partial_xs"],
            "s3_residuals": residues, "target": fixture["equation"]["target"]}


def replay(artifact: Path, s3_tar: Path) -> dict:
    artifact = artifact.resolve(strict=True)
    s3_tar = s3_tar.resolve(strict=True)
    spec = importlib.util.spec_from_file_location("n83_protocol_audit", PROTO)
    require(spec is not None and spec.loader is not None, "protocol auditor unavailable")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    protocol = module.replay(artifact, SOURCE_COMMIT)
    saved_protocol = json.loads((HERE / "protocol-audit.json").read_text())
    require(protocol == saved_protocol
            and protocol["status"] == "SAT_verified_witness_x"
            and protocol["solver_conflicts"] == 787
            and protocol["model_verified"] is True
            and [m["phase"] for m in protocol["phases"]]
                == ["build_started", "circuit_built", "domains_built", "fixed_units_installed"],
            "protocol audit or phase replay differs")
    raw = json.loads((artifact / "solve.stdout").read_text())
    require(raw["status"] == "SAT_verified_witness_x"
            and raw["model_verified"] is True
            and raw["conflicts"] == 787 and raw["conflict_cap"] == 5000
            and raw["identity_mask"] == 0 and raw["fixed_x_blocks"] == 8
            and raw["input_blake3"] == INPUT_BLAKE3
            and raw["fixture_sha256"] == FIXTURE_SHA256
            and raw["before_fixed_units"] == {
                "variables": 101636, "clauses": 714072,
                "xors": 4980, "domain_clauses": 465738,
            }
            and raw["variables"] == 101636 and raw["clauses"] == 714736
            and raw["xors"] == 4980
            and raw["natural_relations"] is None
            and raw["total_index_calculus_runtime_ms"] is None
            and raw["selected_best_total_runtime"] is None,
            "solver witness result differs")
    require(0 <= raw["elapsed_wall_ms_L0"] <= 1_200_000
            and abs(raw["elapsed_wall_ms_L0"] / 1000
                    - protocol["elapsed_wall_seconds_L0"]) < 2,
            "worker and external wall differ")
    receipt = json.loads((HERE / "s3-roundtrip.json").read_text())
    require(receipt["status"] == "PASS_CONTENT_ADDRESSED_GET_AND_IDENTICAL_AUDIT"
            and receipt["run_id"] == RUN_ID
            and receipt["source_commit"] == SOURCE_COMMIT
            and receipt["byte_identical_replay"] is True
            and receipt["protocol_audit_sha256"] == digest((HERE / "protocol-audit.json").read_bytes())
            and receipt["s3_replay_audit_sha256"] == receipt["protocol_audit_sha256"]
            and receipt["archive_bytes"] == s3_tar.stat().st_size
            and receipt["archive_sha256"] == digest(s3_tar.read_bytes())
            and receipt["s3_get_sha256"] == receipt["archive_sha256"]
            and receipt["s3_uri"].endswith(receipt["archive_sha256"] + ".tar.gz"),
            "S3 roundtrip receipt differs")
    seen = {}
    with tarfile.open(s3_tar, "r:gz") as tar:
        for member in tar:
            require(member.isfile() and member.name.startswith(PREFIX)
                    and "/" not in member.name[len(PREFIX):],
                    "S3 artifact member path differs")
            name = member.name[len(PREFIX):]
            require(name not in seen, "duplicate S3 artifact member")
            seen[name] = digest(tar.extractfile(member).read())
    require(seen == protocol["artifact_file_sha256"],
            "S3 archived files differ from GitHub artifact")
    return {
        "schema": "n83.hosted-fixed-witness-retained-audit/v1",
        "status": "PASS_FIXED_WITNESS_SAT_AND_S3_ROUNDTRIP",
        "run_id": RUN_ID, "source_commit": SOURCE_COMMIT,
        "github_artifact_id": receipt["github_artifact_id"],
        "artifact_file_count": len(seen),
        "artifact_sha256": seen,
        "s3_uri": receipt["s3_uri"],
        "s3_archive_sha256": receipt["archive_sha256"],
        "protocol_audit_sha256": receipt["protocol_audit_sha256"],
        "algebra": algebra_replay(),
        "solver_conflicts": raw["conflicts"],
        "worker_wall_ms_L0": raw["elapsed_wall_ms_L0"],
        "outer_wall_seconds_L0": protocol["elapsed_wall_seconds_L0"],
        "phases": protocol["phases"],
        "natural_relations": None, "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--s3-tar", type=Path, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    output = json.dumps(replay(args.artifact, args.s3_tar), indent=2, sort_keys=True) + "\n"
    saved = HERE / "audit.json"
    if args.write:
        with saved.open("x") as stream:
            stream.write(output)
    else:
        require(saved.read_text() == output, "saved audit differs")
    print(json.dumps({"status": "PASS", "audit_sha256": digest(output.encode())}, sort_keys=True))


if __name__ == "__main__":
    main()
