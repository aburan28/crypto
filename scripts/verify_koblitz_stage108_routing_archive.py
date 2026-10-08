#!/usr/bin/env python3
"""Replay and summarize the archived Stage 105--107 routing evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import sys
import tarfile
import tempfile
from typing import Any


REPO = Path(__file__).resolve().parents[1]
EVIDENCE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
FROZEN = EVIDENCE / "stage-108-routing-selection-archive-20260913"
RECEIPTS = FROZEN / "hosted-receipts.json"


class Stage108Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage108Error(message)


def load(path: Path, context: str) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{context} must be a JSON object")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_extract(archive: Path, destination: Path) -> Path:
    with tarfile.open(archive, "r:gz") as source:
        members = source.getmembers()
        require(bool(members), f"{archive.name} is empty")
        root = destination.resolve()
        for member in members:
            target = (destination / member.name).resolve()
            require(target == root or root in target.parents, f"{archive.name} escapes extraction root")
            require(not member.issym() and not member.islnk(), f"{archive.name} contains a link")
        source.extractall(destination, filter="data")
    tops = {member.name.split("/", 1)[0] for member in members if member.name}
    require(len(tops) == 1, f"{archive.name} has multiple roots")
    return destination / next(iter(tops))


def replay_one(receipt: dict[str, Any], temporary: Path) -> dict[str, Any]:
    run_id = receipt["run_id"]
    archive = EVIDENCE / receipt["archive"]
    require(archive.is_file(), f"missing {archive.name}")
    require(sha256(archive) == receipt["archive_sha256"], f"{archive.name} digest changed")
    sidecar = archive.with_name(archive.name + ".sha256")
    require(sidecar.read_text() == f"{receipt['archive_sha256']}  {archive.name}\n", f"{archive.name} sidecar changed")
    root = safe_extract(archive, temporary / str(run_id))
    prefix = receipt["prefix"]
    expected = load(root / f"{prefix}-verification.json", f"run {run_id} hosted verification")
    completed = subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts" / receipt["runner"]),
            "verify",
            "--build",
            str(root / f"{prefix}-build"),
            "--output",
            str(root / f"{prefix}-run"),
        ],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=True,
    )
    require(json.loads(completed.stdout) == expected, f"run {run_id} replay differs from hosted verification")
    return expected


def validate_host(host: dict[str, Any]) -> None:
    require(host.get("schema") == "koblitz_linux_host_identity.v1", "host receipt schema changed")
    payload = dict(host)
    claimed = payload.pop("identity_sha256", None)
    actual = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    require(claimed == actual, "host receipt self-hash changed")


def compose() -> dict[str, Any]:
    receipts = load(RECEIPTS, "Stage-108 hosted receipts")
    require(receipts.get("schema") == "koblitz_stage108_hosted_receipts.v1", "receipt schema changed")
    runs = receipts.get("runs")
    require(isinstance(runs, list) and len(runs) == 8, "receipt run count changed")
    with tempfile.TemporaryDirectory(prefix="koblitz-stage108-") as directory:
        temporary = Path(directory)
        replays = {str(row["run_id"]): replay_one(row, temporary) for row in runs}

    stage105 = replays["34733277351"]
    require(stage105["selected_over_baseline_wall_ratio"] < 1.0, "Stage-105 sharding speedup changed")
    require(stage105["selected_setup_speedup"] > 1.0, "Stage-105 setup speedup changed")
    stage106_ids = [34733682678, 34733871575, 34733875815, 34733877739, 34733880006]
    stage106 = [replays[str(run_id)] for run_id in stage106_ids]
    for result in stage106:
        require(result["direct_over_mixed_wall_ratio"] < 1.0, "Stage-106 routing loss appeared")
        require(result["direct_core_seconds_ratio"] < 1.0, "Stage-106 routing core regression appeared")
        require(result["direct_over_rho_wall_ratio"] < 1.0, "Stage-106 online rho crossover changed")
        validate_host(result["host_identity"])
    route_ratios = [result["direct_over_mixed_wall_ratio"] for result in stage106]
    rho_ratios = [result["direct_over_rho_wall_ratio"] for result in stage106]
    core_ratios = [result["direct_core_seconds_ratio"] for result in stage106]

    unknown = replays["34734272952"]
    require(unknown["factor_base_logs_known_by_construction"] is False, "unknown run gained constructed base logs")
    require(unknown["target_scalar_constructed_or_supplied"] is False, "unknown run gained a supplied scalar")
    require(unknown["whole_process_crossover"] is True, "unknown run online crossover changed")
    validate_host(unknown["host_identity"])

    panel = replays["34734274502"]
    require(panel["direct_wins"] == 5, "selected panel wins changed")
    require(panel["median_gate_passed"] is True, "selected panel gate changed")
    require(panel["fresh_build_plus_median_direct_over_median_rho_wall_ratio"] > 1.0, "full-build boundary changed")
    validate_host(panel["host_identity"])

    host_models = sorted({result["host_identity"]["cpu"]["model name"] for result in stage106})
    return {
        "schema": "koblitz_stage108_routing_selection_archive.v1",
        "status": "n53_routing_selection_archive_verified",
        "stage105": {
            "run_id": 34733277351,
            "selected_candidate": stage105["selected_candidate"],
            "wall_ratio": stage105["selected_over_baseline_wall_ratio"],
            "core_ratio": stage105["selected_core_seconds_ratio"],
            "setup_speedup": stage105["selected_setup_speedup"],
            "selected_over_rho_wall_ratio": stage105["selected_over_rho_wall_ratio"],
        },
        "stage106": {
            "run_ids": stage106_ids,
            "host_models": host_models,
            "direct_route_wins": sum(ratio < 1.0 for ratio in route_ratios),
            "rho_online_wins": sum(ratio < 1.0 for ratio in rho_ratios),
            "route_wall_ratios": route_ratios,
            "median_route_wall_ratio": statistics.median(route_ratios),
            "median_route_core_ratio": statistics.median(core_ratios),
            "direct_over_rho_wall_ratios": rho_ratios,
            "median_direct_over_rho_wall_ratio": statistics.median(rho_ratios),
        },
        "selected_unknown_scalar": {
            "run_id": 34734272952,
            "host_identity": unknown["host_identity"],
            "direct_wall_seconds": unknown["direct_wall_seconds"],
            "rho_wall_seconds": unknown["rho_wall_seconds"],
            "direct_over_rho_wall_ratio": unknown["direct_over_rho_wall_ratio"],
            "recovered_scalar": unknown["recovered_scalar"],
            "whole_process_crossover": unknown["whole_process_crossover"],
        },
        "selected_panel": {
            "run_id": 34734274502,
            "host_identity": panel["host_identity"],
            "paired_wall_ratios": panel["paired_wall_ratios"],
            "median_paired_wall_ratio": panel["median_paired_wall_ratio"],
            "median_direct_wall_seconds": panel["median_direct_wall_seconds"],
            "median_rho_wall_seconds": panel["median_rho_wall_seconds"],
            "direct_wins": panel["direct_wins"],
            "fresh_build_ratio": panel["fresh_build_plus_median_direct_over_median_rho_wall_ratio"],
        },
        "hosted_receipts_sha256": sha256(RECEIPTS),
        "factor_base_selection_uses_scalar_labels": False,
        "target_subgroup_enumerated_for_factor_base": False,
        "factor_base_logs_known_by_construction": False,
        "licensed_magma_complete": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }


def write_new(path: Path, value: dict[str, Any]) -> None:
    require(not path.exists(), f"refusing to overwrite {path}")
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def build_frozen(output: Path) -> dict[str, Any]:
    require(not output.exists(), f"refusing to overwrite {output}")
    output.mkdir(parents=True)
    result = compose()
    result_path = output / "verification.json"
    write_new(result_path, result)
    receipt_rows = load(RECEIPTS, "Stage-108 hosted receipts")["runs"]
    seal = {
        "schema": "koblitz_stage108_routing_selection_archive_seal.v1",
        "status": "result_frozen",
        "verification_sha256": sha256(result_path),
        "hosted_receipts_sha256": sha256(RECEIPTS),
        "archives": {str(row["run_id"]): row["archive_sha256"] for row in receipt_rows},
    }
    write_new(output / "result-seal.json", seal)
    return result


def verify_frozen(output: Path) -> dict[str, Any]:
    result_path = output / "verification.json"
    seal = load(output / "result-seal.json", "Stage-108 result seal")
    require(seal.get("schema") == "koblitz_stage108_routing_selection_archive_seal.v1", "seal schema changed")
    require(sha256(result_path) == seal.get("verification_sha256"), "verification seal changed")
    require(sha256(RECEIPTS) == seal.get("hosted_receipts_sha256"), "receipt seal changed")
    current = compose()
    require(current == load(result_path, "Stage-108 frozen verification"), "frozen verification changed")
    return current


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    build = sub.add_parser("build")
    build.add_argument("--output", type=Path, required=True)
    verify = sub.add_parser("verify")
    verify.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        value = build_frozen(args.output.resolve()) if args.command == "build" else verify_frozen(args.output.resolve(strict=True))
        print(json.dumps(value, indent=2, sort_keys=True))
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage108Error) as error:
        raise SystemExit(f"stage108-routing-archive: {error}")


if __name__ == "__main__":
    main()
