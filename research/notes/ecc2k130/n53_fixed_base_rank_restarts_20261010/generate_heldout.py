#!/usr/bin/env python3
"""Generate the preregistered public held-out point once, preserving its fixture."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924"))
from independent_replay import Curve, Field  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=False)


def write_new(path: Path, contents: str) -> None:
    with path.open("x") as stream:
        stream.write(contents)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generator-binary", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    manifest_path = HERE / "HELDOUT_GENERATION.json"
    manifest = json.loads(manifest_path.read_text())
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    pilot = HERE / "PILOT_ANALYSIS.json"
    assert sha(pilot) == manifest["pilot_analysis_sha256"]
    assert json.loads(pilot.read_text())["selected_cap_for_heldout"] == "cap400000"
    assert sha(ROOT / manifest["generator_source"]) == manifest["generator_source_sha256"]
    assert sha(ROOT / manifest["rho_online_source"]) == manifest["rho_online_source_sha256"]
    assert sha(args.generator_binary) == manifest["generator_binary_sha256"]
    assert sha(ROOT / "examples/koblitz_orbit_dlp_fast_online.rs") == frozen["source_sha256"]
    committed = git("show", "HEAD:research/notes/ecc2k130/n53_fixed_base_rank_restarts_20261010/HELDOUT_GENERATION.json")
    assert committed.returncode == 0 and committed.stdout == manifest_path.read_text(), "generation rule is not committed"
    head, upstream = git("rev-parse", "HEAD"), git("rev-parse", "@{u}")
    assert head.returncode == upstream.returncode == 0 and head.stdout == upstream.stdout, "generation rule is not pushed"
    command = [str(args.generator_binary.resolve()) if word == "{generator_binary}" else word
               for word in manifest["generator_command"]]
    if args.dry_run:
        print(json.dumps({"status": "DRY_RUN", "command": command,
                          "head": head.stdout.strip()}, sort_keys=True))
        return

    fixture_path = HERE / "inputs/heldout_fixture_verifier_only.jsonl"
    point_path = HERE / "inputs/heldout_q.jsonl"
    stderr_path = HERE / "inputs/heldout_generation.stderr"
    receipt_path = HERE / "inputs/heldout_generation_receipt.json"
    assert not any(path.exists() for path in
                   (fixture_path, point_path, stderr_path, receipt_path)), "held-out point already generated"
    environment = os.environ.copy()
    for key in ("KIC_RHO_WALK_SEED", "KIC_RHO_FIXED_TARGET_SCALAR",
                "KIC_RHO_TARGET_POINT", "KIC_RHO_BATCH_CORPUS", "KIC_RHO_FIXTURE_OFFSET"):
        environment.pop(key, None)
    receipt = {
        "schema": "n53-rank-restart-heldout-generation-receipt-v1",
        "manifest_sha256": sha(manifest_path),
        "head_commit": head.stdout.strip(),
        "generator_binary_sha256": sha(args.generator_binary),
        "command": command,
    }
    try:
        process = subprocess.run(command, env=environment, cwd=ROOT, capture_output=True,
                                 text=True, timeout=manifest["generator_wall_cap_seconds"])
        write_new(fixture_path, process.stdout)
        write_new(stderr_path, process.stderr)
        receipt["exit_code"] = process.returncode
        assert process.returncode == 0, f"generator exit {process.returncode}"
        records = [json.loads(line) for line in process.stdout.splitlines() if line.strip()]
        assert len(records) == 1
        fixture = records[0]
        assert fixture["verified"] is True and fixture["n"] == 53 and fixture["a"] == 0
        assert fixture["target_kind"] == "public_hash_to_curve_cofactor"
        assert fixture["public_hash_seed"] == manifest["public_hash_seed"]
        assert fixture["published_fixture_scalar"] is None
        point = fixture["published_q"]
        assert len(point) == 2 and all(isinstance(x, int) for x in point)
        curve = Curve(Field(53, fixture["field_modulus_low_terms"]), 0)
        generator = tuple(fixture["generator"])
        assert curve.on_curve(tuple(point)) and curve.on_curve(generator)
        assert curve.mul(int(fixture["subgroup_order"]), tuple(point)) is None
        assert curve.mul(int(fixture["recovered_fixture_scalar"]), generator) == tuple(point)
        write_new(point_path, json.dumps(point, separators=(",", ":")) + "\n")
        old_point = json.loads((HERE / "inputs/development_q.jsonl").read_text())
        assert point != old_point, "generated point duplicates development point"
        patterns = ([f"[{point[0]},{point[1]}]", f"[{point[0]}, {point[1]}]"])
        prior = git("grep", "-l", "-F", "-e", patterns[0], "-e", patterns[1],
                    "HEAD", "--", "research", "experiments", "docs")
        assert prior.returncode in (0, 1), f"prior-input search exit {prior.returncode}: {prior.stderr}"
        matches = prior.stdout.splitlines()
        receipt.update({
            "status": "ELIGIBLE" if not matches else "DUPLICATE_PRIOR_INPUT",
            "point": point,
            "public_hash_counter": fixture["public_hash_counter"],
            "fixture_sha256": sha(fixture_path),
            "point_sha256": sha(point_path),
            "prior_committed_matches": matches,
            "independent_scalar_replay": True,
        })
    except BaseException as error:
        receipt["status"] = "INVALID_OR_INCOMPLETE"
        receipt["error_type"] = type(error).__name__
        receipt["error"] = str(error)
        receipt["traceback"] = traceback.format_exc()
    write_new(receipt_path, json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(json.dumps({key: receipt.get(key) for key in
                      ("status", "point", "public_hash_counter", "prior_committed_matches")},
                     sort_keys=True))


if __name__ == "__main__":
    main()
