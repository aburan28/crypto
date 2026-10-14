#!/usr/bin/env python3
"""Replay the N83 v7 witness through the first all-affine S3 identity pattern.

This checks group and field equations only. It never constructs CNF or calls SAT.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
STUDY = ROOT / "research/koblitz_n83_factor_base_sweep_20261008"
BASELINE = "d42167ff4e80a9507911338df35572b5d24021c5"
FILES = {
    "research/koblitz_n83_factor_base_sweep_20261008/verification/mixed-orbit-known-witness-source-20261012/fixture.json": "fc203db28eed0164e47b0c755c3d9a09f80df1acd5e4c7e09dfe172fa5f0a538",
    "research/koblitz_n83_factor_base_sweep_20261008/verification/mixed-orbit-control-reuse-19800-20261013/evidence-01/run/control/result.json": "e1c93b0e9d3d91ccf7eca304f39723298fc96bed4aba24f97101ab9e74811465",
    "research/koblitz_n83_factor_base_sweep_20261008/verification/mixed-orbit-control-reuse-19800-20261013/evidence-01/run/audit.json": "6c0926ced564fd5c9e0fff229667cfc1d036f8877974b525270602e27208df6e",
    "research/koblitz_n83_factor_base_sweep_20261008/chain_search_supervisor.py": "fa0443e453e5fb7011019a08b61a219834e1594acc99e007f405f27e44b97a66",
    "src/cryptanalysis/koblitz_large_prime.rs": "6ac5580c444cc65ba3a095fe117557c9c8628283a2e725fcbfa117f144d8b987",
    "src/cryptanalysis/binary_semaev_chain_sat.rs": "0031972aab435d963af54a5f28f0510e45e5301f00c49466661191815becc06e",
}
sys.path.insert(0, str(STUDY))
from chain_search_supervisor import field_mul, parse_point, point_add  # noqa: E402


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def s3_residual(x: int, y: int, z: int) -> int:
    """Binary Koblitz S3: xyz + (xy + xz + yz)^2 + b, with b=1."""
    xy = field_mul(x, y)
    pair_sum = xy ^ field_mul(x, z) ^ field_mul(y, z)
    return field_mul(xy, z) ^ field_mul(pair_sum, pair_sum) ^ 1


def replay() -> dict:
    for relative, expected in FILES.items():
        path = ROOT / relative
        raw = path.read_bytes()
        require(digest(raw) == expected, f"pinned file differs: {relative}")
        committed = subprocess.run(
            ["git", "show", f"{BASELINE}:{relative}"],
            cwd=ROOT, check=True, capture_output=True,
        ).stdout
        require(digest(committed) == expected, f"baseline Git blob differs: {relative}")

    verification = STUDY / "verification"
    fixture = json.loads((verification / "mixed-orbit-known-witness-source-20261012/fixture.json").read_text())
    run = verification / "mixed-orbit-control-reuse-19800-20261013/evidence-01/run"
    result = json.loads((run / "control/result.json").read_text())
    audit = json.loads((run / "audit.json").read_text())
    require(fixture.get("schema") == "n83.mixed-orbit-known-witness/v1"
            and fixture.get("scope") == "known_answer_control_only"
            and fixture.get("equation", {}).get("a") == "0"
            and fixture.get("equation", {}).get("b") == "1"
            and fixture.get("equation", {}).get("small_summands") == 3
            and fixture.get("equation", {}).get("residual_summands") == 2,
            "fixture scope or shape differs")
    stats = result.get("search", {})
    require(result.get("status") == "UNKNOWN_search_cap"
            and result.get("candidate") is None
            and result.get("target") == fixture["equation"]["target"]
            and stats.get("patterns_attempted") == 1
            and stats.get("patterns_refuted") == 0
            and stats.get("solver_calls") == 1
            and stats.get("conflicts") == 5000
            and stats.get("models") == 0
            and stats.get("exhausted") is True
            and stats.get("refuted") is False,
            "bounded first-pattern receipt differs")
    require(audit.get("status") == "VALIDATED_CENSORED"
            and audit.get("worker_status") == "UNKNOWN_search_cap"
            and audit.get("raw_sha256", {}).get("control/result.json")
                == FILES["research/koblitz_n83_factor_base_sweep_20261008/verification/mixed-orbit-control-reuse-19800-20261013/evidence-01/run/control/result.json"],
            "independent audit or raw binding differs")
    # The pinned Rust search visits identity_mask=0 first. In that pattern,
    # all intermediate sums are affine and every chained S3 node is active.
    source = (ROOT / "src/cryptanalysis/koblitz_large_prime.rs").read_text()
    require("for identity_mask in 0..(1u32 << (m - 2))" in source,
            "first-pattern search order differs")
    terms = fixture.get("terms")
    require(isinstance(terms, list) and len(terms) == 5,
            "known-witness term count differs")
    expected_indices = [0, 166, 332, 196212, 196378]
    require([term.get("global_point_index") for term in terms] == expected_indices
            and [term.get("domain") for term in terms] == ["small"] * 3 + ["residual"] * 2
            and all(term.get("coefficient") == "1" and term.get("phase") == 0
                    and term.get("negative") is False for term in terms),
            "known-witness selection differs")
    points = [parse_point(term["point"]) for term in terms]
    require(all(point is not None for point in points), "known-witness point is infinity")
    partials = []
    acc = None
    for point in points:
        acc = point_add(acc, point)
        require(acc is not None, "witness has an intermediate identity")
        partials.append(acc)
    require(partials[-1] == parse_point(fixture["equation"]["target"]),
            "witness sum differs from target")
    residuals = []
    for step in range(1, len(points)):
        left = points[0] if step == 1 else partials[step - 1]
        residual = s3_residual(left[0], points[step][0], partials[step][0])
        require(residual == 0, f"chained S3 node {step} differs")
        residuals.append(residual)
    return {
        "schema": "n83.known-witness-first-pattern-algebra-replay/v1",
        "status": "PASS_group_and_S3_witness_for_first_identity_pattern",
        "scope": "known_answer_control_only",
        "baseline_commit": BASELINE,
        "input_sha256": FILES,
        "identity_mask": 0,
        "partial_sums": [[str(x), str(y)] for x, y in partials],
        "s3_residuals": residuals,
        "target_matches": True,
        "capped_worker_conflicts": 5000,
        "capped_worker_models": 0,
        "SAT_encoding_satisfiability": None,
        "natural_relations": None,
        "selected_best_total_runtime": None,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    path = HERE / "replay.json"
    raw = (json.dumps(replay(), sort_keys=True, indent=2) + "\n").encode()
    if args.write:
        with path.open("xb") as stream:
            stream.write(raw)
    else:
        require(path.read_bytes() == raw, "committed replay receipt differs")
    print("PASS: first all-affine pattern contains the pinned group/S3 witness")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
