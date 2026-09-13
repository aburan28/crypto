#!/usr/bin/env python3
"""Replay archived n=53 Stages 100--102 and the host-class receipt."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
from typing import Any


REPO = Path(__file__).resolve().parents[1]
EVIDENCE = REPO / "research/sat_factor_base_review_20260908/continuation-05-sota-gates"
FROZEN = EVIDENCE / "stage-103-host-panel-archive-20260913"
RECEIPTS = FROZEN / "hosted-receipts.json"
RUNNERS = {
    100: "run_koblitz_stage100_fixed_base_validation.py",
    101: "run_koblitz_stage101_selected_panel.py",
    102: "run_koblitz_stage102_host_class_panel.py",
}


class Stage103Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage103Error(message)


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
    stage = receipt["stage"]
    require(stage in RUNNERS, f"unsupported stage {stage}")
    archive = EVIDENCE / receipt["archive"]
    require(archive.is_file(), f"missing {archive.name}")
    require(sha256(archive) == receipt["archive_sha256"], f"{archive.name} digest changed")
    sidecar = archive.with_name(archive.name + ".sha256")
    require(sidecar.read_text() == f"{receipt['archive_sha256']}  {archive.name}\n", f"{archive.name} sidecar changed")
    root = safe_extract(archive, temporary / f"stage{stage}")
    build = root / f"stage{stage}-build"
    run = root / f"stage{stage}-run"
    expected = load(root / f"stage{stage}-verification.json", f"Stage-{stage} hosted verification")
    completed = subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts" / RUNNERS[stage]),
            "verify",
            "--build",
            str(build),
            "--output",
            str(run),
        ],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=True,
    )
    replay = json.loads(completed.stdout)
    require(replay == expected, f"Stage-{stage} replay differs from hosted verification")
    return expected


def compose() -> dict[str, Any]:
    receipts = load(RECEIPTS, "Stage-103 hosted receipts")
    require(receipts.get("schema") == "koblitz_stage103_hosted_receipts.v1", "receipt schema changed")
    runs = receipts.get("runs")
    require(isinstance(runs, list) and [row.get("stage") for row in runs] == [100, 101, 102], "receipt stage order changed")
    with tempfile.TemporaryDirectory(prefix="koblitz-stage103-") as directory:
        temporary = Path(directory)
        replays = {str(row["stage"]): replay_one(row, temporary) for row in runs}
    stage100 = replays["100"]
    stage101 = replays["101"]
    stage102 = replays["102"]
    require(stage100["candidate_over_baseline_wall_ratio"] < 1.0, "Stage-100 validation candidate did not improve wall time")
    require(stage100["candidate_whole_process_crossover"] is True, "Stage-100 online crossover changed")
    require(stage101["direct_wins"] == 0 and stage101["median_gate_passed"] is False, "Stage-101 host-sensitive loss changed")
    require(stage102["direct_wins"] == 5 and stage102["median_gate_passed"] is True, "Stage-102 fixed-host panel status changed")
    host = stage102.get("host_identity")
    require(isinstance(host, dict) and host.get("schema") == "koblitz_linux_host_identity.v1", "Stage-102 host receipt is missing")
    payload = dict(host)
    claimed = payload.pop("identity_sha256", None)
    require(claimed == hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest(), "Stage-102 host receipt self-hash changed")
    require(host["cpu"].get("model name") == "AMD EPYC 7763 64-Core Processor", "Stage-102 host model changed")
    return {
        "schema": "koblitz_stage103_host_panel_archive.v1",
        "status": "n53_host_panel_archive_verified",
        "stages": replays,
        "fixed_base_validation": {
            "wall_ratio": stage100["candidate_over_baseline_wall_ratio"],
            "core_ratio": stage100["candidate_over_baseline_core_ratio"],
            "solution_validation_speedup": stage100["candidate_solution_validation_speedup"],
            "relation_validation_speedup": stage100["candidate_reference_validation_speedup"],
        },
        "anonymous_host_panel": {
            "paired_wall_ratios": stage101["paired_wall_ratios"],
            "median_paired_wall_ratio": stage101["median_paired_wall_ratio"],
            "median_direct_wall_seconds": stage101["median_direct_wall_seconds"],
            "median_rho_wall_seconds": stage101["median_rho_wall_seconds"],
            "median_gate_passed": stage101["median_gate_passed"],
        },
        "identified_host_panel": {
            "host_identity": host,
            "paired_wall_ratios": stage102["paired_wall_ratios"],
            "direct_wins": stage102["direct_wins"],
            "median_paired_wall_ratio": stage102["median_paired_wall_ratio"],
            "median_gate_threshold": stage102["median_gate_threshold"],
            "median_gate_passed": stage102["median_gate_passed"],
            "median_direct_wall_seconds": stage102["median_direct_wall_seconds"],
            "median_rho_wall_seconds": stage102["median_rho_wall_seconds"],
            "fresh_build_ratio": stage102["fresh_build_plus_median_direct_over_median_rho_wall_ratio"],
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
    seal = {
        "schema": "koblitz_stage103_host_panel_archive_seal.v1",
        "status": "result_frozen",
        "verification_sha256": sha256(result_path),
        "hosted_receipts_sha256": sha256(RECEIPTS),
        "archives": {
            str(row["stage"]): row["archive_sha256"]
            for row in load(RECEIPTS, "Stage-103 hosted receipts")["runs"]
        },
    }
    write_new(output / "result-seal.json", seal)
    return result


def verify_frozen(output: Path) -> dict[str, Any]:
    result_path = output / "verification.json"
    seal = load(output / "result-seal.json", "Stage-103 result seal")
    require(seal.get("schema") == "koblitz_stage103_host_panel_archive_seal.v1", "seal schema changed")
    require(sha256(result_path) == seal.get("verification_sha256"), "verification seal changed")
    require(sha256(RECEIPTS) == seal.get("hosted_receipts_sha256"), "receipt seal changed")
    current = compose()
    require(current == load(result_path, "Stage-103 frozen verification"), "frozen verification changed")
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage103Error) as error:
        raise SystemExit(f"stage103-host-panel: {error}")


if __name__ == "__main__":
    main()
