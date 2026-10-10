#!/usr/bin/env python3
"""Replay one strong-rho public-point scalar with independent group arithmetic."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924"))
from independent_replay import Curve, Field  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rho", type=Path, required=True)
    parser.add_argument("--public-point", type=Path, required=True)
    parser.add_argument("--verifier-fixture", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "refusing to overwrite rho replay"
    rows = [json.loads(line) for line in args.rho.read_text().splitlines() if line.strip()]
    assert len(rows) == 2
    rho, summary = rows
    point = json.loads(args.public_point.read_text())
    fixture = json.loads(args.verifier_fixture.read_text())
    assert fixture["published_q"] == point and fixture["published_fixture_scalar"] is None
    assert (rho["kind"], summary["kind"]) == (
        "rho_ks_batch_fixture", "rho_ks_batch_summary"
    )
    assert rho["n"] == summary["n"] == fixture["n"] == 53
    assert rho["a"] == summary["a"] == fixture["a"] == 0
    assert rho["batch_seed"] == summary["batch_seed"] == args.seed
    assert rho["target_source"] == summary["target_source"] == "public_point"
    assert rho["published_q"] == point and rho["published_fixture_scalar"] is None
    assert rho["verified"] is True and summary["all_verified"] is True
    assert summary["fixtures"] == 1 and summary["rung"] == 3
    assert (summary["lanes"], summary["dp_bits"], summary["jump_count"]) == (32, 4, 32)
    assert summary["corpus"] is None and summary["explicit_scalar"] is None
    assert rho["online_start_event"] == "after_target_built"
    assert rho["online_stop_event"] == "recovery_check_true"
    assert rho["online_stop_ns"] - rho["online_start_ns"] == rho["online_ns"]
    assert abs(rho["online_ms"] - sum(rho[key] for key in
               ("walk_ms", "collision_ms", "recovery_check_ms"))) < 1e-6
    curve = Curve(Field(53, fixture["field_modulus_low_terms"]), 0)
    generator = tuple(fixture["generator"])
    scalar = int(rho["recovered_fixture_scalar"])
    assert curve.on_curve(generator) and curve.on_curve(tuple(point))
    assert curve.mul(int(fixture["subgroup_order"]), generator) is None
    assert curve.mul(scalar, generator) == tuple(point)
    assert scalar == fixture["recovered_fixture_scalar"]
    receipt = {
        "schema": "n53-strong-rho-public-point-replay-v1",
        "status": "PASS",
        "seed": args.seed,
        "point": point,
        "recovered_scalar": scalar,
        "online_ms": rho["online_ms"],
        "walk_steps": rho["walk_steps"],
        "rho_sha256": sha(args.rho),
        "public_point_sha256": sha(args.public_point),
        "fixture_verifier_only_sha256": sha(args.verifier_fixture),
    }
    args.out.write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
