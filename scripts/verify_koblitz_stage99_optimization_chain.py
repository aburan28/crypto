#!/usr/bin/env python3
"""Replay the archived n=53 Stage 94--98 optimization chain."""

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
FROZEN = EVIDENCE / "stage-99-optimization-chain-20260913"
RECEIPTS = FROZEN / "hosted-receipts.json"
RUNNERS = {
    94: "run_koblitz_stage94_filter_inverse.py",
    95: "run_koblitz_stage95_itoh_width.py",
    96: "run_koblitz_stage96_blocked_filter.py",
    97: "run_koblitz_stage97_selected_panel.py",
    98: "run_koblitz_stage98_compact_evidence.py",
}


class Stage99Error(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise Stage99Error(message)


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
    receipts = load(RECEIPTS, "Stage-99 hosted receipts")
    require(receipts.get("schema") == "koblitz_stage99_hosted_receipts.v1", "receipt schema changed")
    runs = receipts.get("runs")
    require(isinstance(runs, list) and [row.get("stage") for row in runs] == [94, 95, 96, 97, 98], "receipt stage order changed")
    with tempfile.TemporaryDirectory(prefix="koblitz-stage99-") as directory:
        temporary = Path(directory)
        replays = {str(row["stage"]): replay_one(row, temporary) for row in runs}
    stage94 = replays["94"]
    stage95 = replays["95"]
    stage96 = replays["96"]
    stage97 = replays["97"]
    stage98 = replays["98"]
    require(stage94["combined_over_mixed_wall_ratio"] < 1.0, "Stage-94 selected stack did not improve wall time")
    require(stage95["best_width"] == 4096, "Stage-95 width selection changed")
    require(stage96["blocked_over_direct_bits_wall_ratio"] > 1.0, "Stage-96 blocked-filter rejection changed")
    require(stage97["direct_wins"] == 5 and stage97["median_gate_passed"] is False, "Stage-97 panel status changed")
    require(stage98["compact_over_full_wall_ratio"] < 1.0, "Stage-98 compact evidence did not improve wall time")
    return {
        "schema": "koblitz_stage99_optimization_chain.v1",
        "status": "n53_optimization_chain_verified",
        "stages": replays,
        "selection": {
            "x_filter": "direct_low_and_high_x_bit_windows",
            "prefiltered_exact_lookup": True,
            "inverse": "itoh_tsujii_52_squares_7_multiplies",
            "parallel_width": 4096,
            "blocked_filter_selected": False,
            "compact_evidence_available": True,
        },
        "current_five_pair_panel": {
            "paired_wall_ratios": stage97["paired_wall_ratios"],
            "direct_wins": stage97["direct_wins"],
            "median_paired_wall_ratio": stage97["median_paired_wall_ratio"],
            "median_gate_threshold": stage97["median_gate_threshold"],
            "median_gate_passed": stage97["median_gate_passed"],
            "median_direct_wall_seconds": stage97["median_direct_wall_seconds"],
            "median_rho_wall_seconds": stage97["median_rho_wall_seconds"],
            "fresh_build_ratio": stage97["fresh_build_plus_median_direct_over_median_rho_wall_ratio"],
        },
        "compact_evidence": {
            "wall_ratio": stage98["compact_over_full_wall_ratio"],
            "core_ratio": stage98["compact_over_full_core_ratio"],
            "stdout_ratio": stage98["compact_stdout_over_full_ratio"],
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
        "schema": "koblitz_stage99_optimization_chain_seal.v1",
        "status": "result_frozen",
        "verification_sha256": sha256(result_path),
        "hosted_receipts_sha256": sha256(RECEIPTS),
        "archives": {
            str(row["stage"]): row["archive_sha256"]
            for row in load(RECEIPTS, "Stage-99 hosted receipts")["runs"]
        },
    }
    write_new(output / "result-seal.json", seal)
    return result


def verify_frozen(output: Path) -> dict[str, Any]:
    result_path = output / "verification.json"
    seal = load(output / "result-seal.json", "Stage-99 result seal")
    require(seal.get("schema") == "koblitz_stage99_optimization_chain_seal.v1", "seal schema changed")
    require(sha256(result_path) == seal.get("verification_sha256"), "verification seal changed")
    require(sha256(RECEIPTS) == seal.get("hosted_receipts_sha256"), "receipt seal changed")
    current = compose()
    require(current == load(result_path, "Stage-99 frozen verification"), "frozen verification changed")
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
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError, Stage99Error) as error:
        raise SystemExit(f"stage99-optimization-chain: {error}")


if __name__ == "__main__":
    main()
