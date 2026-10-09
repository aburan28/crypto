#!/usr/bin/env python3
"""Independently audit the immutable large-prime pool count and bindings.

This reads receipts and does not re-open retained point objects. The original
bounded header-prefix check is separate evidence.
"""

import hashlib
import json
import math
from pathlib import Path


PANEL = Path(__file__).resolve().parent / "pilot-01"
SIZES = (64, 256, 600)
PAIRS = ((64, 256), (64, 600), (256, 600))


def multisets(points: int, count: int) -> int:
    return 1 if count == 0 else math.comb(points + count - 1, count)


def main() -> None:
    manifest_bytes = (PANEL / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    replay = json.loads((PANEL / "replay.json").read_bytes())
    upload = json.loads((PANEL / "upload-receipt.json").read_bytes())
    screen_bytes = (PANEL / "large-prime-pool-screen.json").read_bytes()
    screen = json.loads(screen_bytes)
    budget = json.loads((PANEL / "large-prime-pool-budget.json").read_bytes())
    assert screen["schema"] == "n83.large-prime-pool-support-screen/v1"
    assert screen["status"] == budget["status"] == "PASS_prefix_and_exact_support"
    assert screen["manifest_sha256"] == hashlib.sha256(manifest_bytes).hexdigest()
    assert budget["screen_sha256"] == hashlib.sha256(screen_bytes).hexdigest()
    assert screen["screen_source_sha256"] == hashlib.sha256(
        (PANEL.parent / "large_prime_pool_screen.py").read_bytes()
    ).hexdigest()
    assert screen["replay_panel_manifest_blake3"] == replay["panel_manifest_blake3"]
    assert upload["manifest_blake3"] == replay["panel_manifest_blake3"]
    assert float(budget["total_active_seconds_approx"]) <= float(
        budget["authorized_active_seconds"]
    )
    assert float(budget["charged_active_seconds"]) >= float(
        budget["measured_outer_interval_seconds"]
    )

    by_uri = {row["s3_uri"]: row for row in manifest["bases"]}
    assert len(by_uri) == len(upload["objects"]) == len(replay["checks"]) == 54
    for item in upload["objects"]:
        assert item["downloaded_hash_matches"] is True
        assert item["compressed_blake3"] == by_uri[item["s3_uri"]][
            "compressed_blake3"
        ]
    by_object = {row["object"]: row for row in manifest["bases"]}
    assert all(
        check["status"] == "PASS"
        and check["points_checked"] == by_object[check["object"]]["points"]
        and check["representatives_checked"]
        == by_object[check["object"]]["columns"]
        for check in replay["checks"]
    )

    bindings = screen["group_bindings"]
    assert len(bindings) == 18
    seen = set()
    for binding in bindings:
        key = (binding["curve_a"], binding["policy"], binding["seed"])
        assert key not in seen and binding["prefix_check"] == "PASS"
        seen.add(key)
        for columns in SIZES:
            record = binding["objects"][str(columns)]
            row = by_object[record["object"]]
            assert (
                row["a"],
                row["policy"],
                row["seed"],
                row["columns"],
                row["s3_uri"],
                row["compressed_blake3"],
                row["point_set_blake3"],
            ) == (
                *key,
                columns,
                record["s3_uri"],
                record["compressed_blake3_from_manifest"],
                record["point_set_blake3"],
            )

    cases = screen["cases"]
    assert len(cases) == 90
    seen_cases = set()
    for case in cases:
        arm = case["curve_a"]
        pair = (case["base_columns"], case["envelope_columns"])
        m, j = case["summands"], case["exact_residuals"]
        assert arm in (0, 1) and pair in PAIRS and 2 <= m <= 6 and 0 <= j <= 2
        key = (arm, pair, m, j)
        assert key not in seen_cases
        seen_cases.add(key)
        b, ell = case["base_points"], case["residual_points"]
        assert b == 166 * pair[0] and ell == 166 * (pair[1] - pair[0])
        exact = multisets(b, m - j) * multisets(ell, j)
        cumulative = sum(
            multisets(b, m - count) * multisets(ell, count)
            for count in range(j + 1)
        )
        assert int(case["unordered_multisets_exact_residuals"]) == exact
        assert int(case["unordered_multisets_at_most_residuals"]) == cumulative
        order = int(case["subgroup_order"])
        ceiling = case["uniform_nonidentity_target_hit_ceiling"]
        assert int(ceiling["numerator"]) == min(cumulative, order - 1)
        assert int(ceiling["denominator"]) == order - 1
    print("PASS 54 object bindings; 18 prefix receipts; 90 exact support cases; budget")


if __name__ == "__main__":
    main()
