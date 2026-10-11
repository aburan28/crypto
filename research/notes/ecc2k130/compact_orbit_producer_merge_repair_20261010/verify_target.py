#!/usr/bin/env python3
"""Replay the recorded target scalar with independent binary-curve arithmetic."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(
    0, str(ROOT / "research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924")
)
from independent_replay import Curve, Field  # noqa: E402


def one_jsonl(path: Path) -> dict:
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    assert len(rows) == 1, f"expected one record in {path}"
    return rows[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "refusing to overwrite replay receipt"

    base = one_jsonl(args.base)
    target = one_jsonl(args.target)
    n, a, r = (int(base[key]) for key in ("n", "a", "subgroup_order"))
    assert target["n"] == n and target["a"] == a
    curve = Curve(Field(n, base["field_modulus_low_terms"]), a)
    generator = tuple(target["generator"])
    point = tuple(target["target"])
    scalar = int(target["recovered_scalar"])
    assert curve.on_curve(generator) and curve.on_curve(point)
    assert curve.mul(r, generator) is None
    assert 0 <= scalar < r and curve.mul(scalar, generator) == point
    assert scalar == target["published_fixture_scalar"]

    receipt = {
        "status": "PASS",
        "schema": "compact-orbit-target-independent-replay-v1",
        "n": n,
        "a": a,
        "subgroup_order": r,
        "recovered_scalar": scalar,
        "generator": list(generator),
        "target": list(point),
        "target_record_sha256": hashlib.sha256(args.target.read_bytes()).hexdigest(),
    }
    args.out.write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
