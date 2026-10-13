#!/usr/bin/env python3
"""Replay the frozen N83 witness through the exact all-affine S3 chain.

This evaluates field and group equations without constructing or solving SAT.
The archived fixture is byte-pinned; the historical archive-import receipts
remain the authority for the selected points' membership in the two domains.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
FIXTURE = HERE.parent / "mixed-orbit-known-witness-source-20261012/fixture.json"
RECEIPT = HERE / "replay.json"
AUDIT = HERE.parent / "k600-mixed-orbit-first-2400-20261012/audit.json"
AUDIT_SHA256 = "387d393b9805f35ee1ca6328e8602133a1b9bc036088b1be1a3ba5a29bb8469c"
FROZEN_SOURCE_COMMIT = "f3d27f68b68cda907e0668042ac7fa13ae6b579d"
FROZEN_SOURCE_SHA256 = {
    "src/cryptanalysis/binary_semaev_chain_sat.rs": "0031972aab435d963af54a5f28f0510e45e5301f00c49466661191815becc06e",
    "src/cryptanalysis/koblitz_large_prime.rs": "accd5cf0ef7c738aedce38a21d172eeda1b309eaf9c8f4848030dac9d4421ddd",
    "examples/koblitz_n83_mixed_orbit_worker.rs": "e8781e8ba5c06528b1a8967312c5621084c0c3f8e3705e0882f849f4b1b9a152",
}
FIXTURE_SHA256 = "fc203db28eed0164e47b0c755c3d9a09f80df1acd5e4c7e09dfe172fa5f0a538"
SMALL_SHA256 = "eb342b574e244e1e67e8e78c805f33ee88edcdb4f3a14d4114e6f413beede380"
ENVELOPE_SHA256 = "aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a"
CURVE = "icv1-f2m83-tm6151469093347-debefd74"
ORDER = 2417851639230796216685689
MODULUS = (1 << 83) | (1 << 45) | 0b111
FIELD_SIZE = 1 << 83
POINT_INDICES = (0, 166, 332, 196212, 196378)


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value: dict) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def multiply(left: int, right: int) -> int:
    """Carryless product followed by polynomial reduction; separate from SAT."""
    raw = 0
    for bit in range(right.bit_length()):
        if (right >> bit) & 1:
            raw ^= left << bit
    for bit in range(raw.bit_length() - 1, 82, -1):
        if (raw >> bit) & 1:
            raw ^= MODULUS << (bit - 83)
    return raw


def inverse(value: int) -> int:
    require(value != 0, "zero field denominator")
    result = 1
    exponent = FIELD_SIZE - 2
    while exponent:
        if exponent & 1:
            result = multiply(result, value)
        value = multiply(value, value)
        exponent >>= 1
    return result


def point_add(left: tuple[int, int] | None,
              right: tuple[int, int] | None) -> tuple[int, int] | None:
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if y1 != y2 or x1 == 0:
            return None
        slope = x1 ^ multiply(y1, inverse(x1))
        x3 = multiply(slope, slope) ^ slope
        return x3, multiply(x1, x1) ^ multiply(slope ^ 1, x3)
    slope = multiply(y1 ^ y2, inverse(x1 ^ x2))
    x3 = multiply(slope, slope) ^ slope ^ x1 ^ x2
    return x3, multiply(slope, x1 ^ x3) ^ x3 ^ y1


def semaev_s3(x: int, y: int, z: int) -> int:
    """S3 = (xy + xz + yz)^2 + xyz + b, with b=1."""
    xy = multiply(x, y)
    pair_sum = xy ^ multiply(x, z) ^ multiply(y, z)
    return multiply(pair_sum, pair_sum) ^ multiply(xy, z) ^ 1


def parse_point(value: object) -> tuple[int, int]:
    require(isinstance(value, list) and len(value) == 2, "point shape")
    require(all(isinstance(word, str) and word.isascii() and word.isdecimal()
                and (word == "0" or not word.startswith("0")) for word in value),
            "point canonical decimal")
    point = tuple(map(int, value))
    require(all(0 <= word < FIELD_SIZE for word in point), "point outside field")
    x, y = point
    require(multiply(y, y) ^ multiply(x, y)
            == multiply(multiply(x, x), x) ^ 1, "point off curve")
    return point


def check_chain(points: list[tuple[int, int]], target: tuple[int, int]) -> dict:
    """Compute every intermediate identity bit and direct S3 residue."""
    require(len(points) == 5, "five selected points required")
    require(points[0][0] <= points[1][0], "first pair violates lex order")
    current = points[0]
    intermediate = []
    residues = []
    identity_mask = 0
    for step, next_point in enumerate(points[1:], 1):
        following = point_add(current, next_point)
        if step < 4 and following is None:
            identity_mask |= 1 << (step - 1)
        require(current is not None and following is not None,
                "fixture leaves the all-affine identity pattern")
        residue = semaev_s3(current[0], next_point[0], following[0])
        require(residue == 0, f"S3 node {step} differs")
        residues.append(str(residue))
        intermediate.append(str(following[0]))
        current = following
    require(current == target, "ordered group sum differs from target")
    require(identity_mask == 0, "all-affine mask differs")
    return {"identity_mask": identity_mask,
            "intermediate_x": intermediate[:-1],
            "target_x": intermediate[-1], "s3_residues": residues}


def replay(fixture_path: Path = FIXTURE) -> dict:
    require(sha256(fixture_path) == FIXTURE_SHA256, "fixture bytes changed")
    require(sha256(AUDIT) == AUDIT_SHA256, "K600 audit bytes changed")
    for relative, expected in FROZEN_SOURCE_SHA256.items():
        require(sha256(ROOT / relative) == expected, f"frozen source changed: {relative}")
    fixture = json.loads(fixture_path.read_text())
    require(fixture["schema"] == "n83.mixed-orbit-known-witness/v1"
            and fixture["scope"] == "known_answer_control_only"
            and fixture["curve_slug"] == CURVE
            and fixture["subgroup_order"] == str(ORDER)
            and fixture["equation"]["a"] == "0"
            and fixture["equation"]["b"] == "1"
            and fixture["equation"]["small_summands"] == 3
            and fixture["equation"]["residual_summands"] == 2
            and fixture["archive"]["small"]["sha256"] == SMALL_SHA256
            and fixture["archive"]["envelope"]["sha256"] == ENVELOPE_SHA256,
            "fixture metadata differs")
    terms = fixture["terms"]
    require(len(terms) == len(POINT_INDICES), "term count differs")
    points = []
    for expected, term in zip(POINT_INDICES, terms):
        require(term["global_point_index"] == expected
                and term["column"] == expected // 166
                and term["phase"] == 0 and term["negative"] is False
                and term["coefficient"] == "1"
                and term["domain"] == ("small" if expected < 196212 else "residual"),
                "archive-selected term differs")
        points.append(parse_point(term["point"]))
    target = parse_point(fixture["equation"]["target"])
    ordered = list(range(len(points)))
    if points[0][0] > points[1][0]:
        points[:2] = points[1::-1]
        ordered[:2] = ordered[1::-1]
    chain = check_chain(points, target)
    return {
        "schema": "n83.mixed-orbit-control-chain-preflight/v1",
        "status": "PASS_all_affine_chain_witness",
        "scope": "known_answer_control_only",
        "fixture_sha256": FIXTURE_SHA256,
        "k600_audit_sha256": AUDIT_SHA256,
        "frozen_source_commit": FROZEN_SOURCE_COMMIT,
        "frozen_source_sha256": FROZEN_SOURCE_SHA256,
        "source_sha256": sha256(Path(__file__)),
        "curve_slug": CURVE,
        "small_archive_sha256": SMALL_SHA256,
        "envelope_archive_sha256": ENVELOPE_SHA256,
        "ordered_global_point_indices": [POINT_INDICES[i] for i in ordered],
        "first_pair_swapped_for_lex_order": ordered[0] != 0,
        "all_selected_terms_phase_zero": True,
        "all_selected_terms_in_archived_domains_by_fixture": True,
        "chain": chain,
        "sat_model_constructed": False,
        "sat_search_executed": False,
        "natural_relation_yield": None,
        "total_index_calculus_runtime_ms": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = canonical(replay())
    if args.write:
        RECEIPT.write_bytes(rendered)
        print(f"wrote {RECEIPT}")
    else:
        require(RECEIPT.read_bytes() == rendered, "saved preflight receipt differs")
        print("PASS: frozen all-affine chained-S3 preflight")


if __name__ == "__main__":
    main()
