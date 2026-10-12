#!/usr/bin/env python3
"""Post-run diagnostic for the archived n37/L1024 pin-order verifier failure.

This does not alter the frozen verifier or confer timing eligibility. It first
checks the producer's recorded intermediate x values against every pairing of
the same four selected points, then delegates every other check to the exact
frozen target and cold-run verifiers.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930"))
from run_panel import sha  # noqa: E402
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"))
from verify_rank import point  # noqa: E402
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_point_panel_20260929"))
from verify_panel import check_target as frozen_target_check  # noqa: E402
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/disjoint_cold_q_20260930"))
import verify_cold as frozen_cold  # noqa: E402

FROZEN_COLD_SHA = "c12f35fe8dc07861df657858a0874c14e685aa5eb76727c32c34095cbcc49f6e"
FROZEN_TARGET_SHA = "649459a53b979b9eaf150ca1973c80b1bde29baf3f5567443c668978d3c06e6d"
PAIRINGS = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))
MISMATCHES: list[dict] = []


def check_any_pairing(record: dict, fixture: dict, base: dict, logs: list[int],
                      curve, generator, checked_labels: set) -> None:
    indices = record["point_indices"]
    assert isinstance(indices, list) and len(indices) == 4
    chosen = [point(base["factor_base_point_coordinates"][index]) for index in indices]
    assert all(selected is not None and curve.on_curve(selected) for selected in chosen)
    pinned = record["pinned_intermediates"]
    assert isinstance(pinned, list) and len(pinned) == 2
    candidates = []
    for pairing in PAIRINGS:
        sums = [curve.add(chosen[i], chosen[j]) for i, j in pairing]
        if all(value is not None for value in sums):
            x_values = [value[0] for value in sums]
            if set(x_values) == set(pinned):
                candidates.append((pairing, x_values))
    assert candidates, (record["fixture_index"], pinned, indices)

    canonical = [curve.add(chosen[0], chosen[1]), curve.add(chosen[2], chosen[3])]
    assert all(value is not None for value in canonical)
    canonical_x = [value[0] for value in canonical]
    if set(canonical_x) != set(pinned):
        MISMATCHES.append({"fixture_index": record["fixture_index"],
                           "published_q": record["published_q"],
                           "point_indices": indices,
                           "recorded_pinned_intermediates": pinned,
                           "canonical_pair_x": canonical_x,
                           "matching_pairings": [[list(pair) for pair in pairing]
                                                 for pairing, _ in candidates]})
    # The original checker still proves the factor-base labels, group sum,
    # known-answer discrete log, and first-pair intermediate coordinates.
    shadow = dict(record)
    shadow["pinned_intermediates"] = canonical_x
    frozen_target_check(shadow, fixture, base, logs, curve, generator, checked_labels)


def replay(run_dir: Path) -> dict:
    assert sha(ROOT / "research/notes/ecc2k130/disjoint_cold_q_20260930/verify_cold.py") == FROZEN_COLD_SHA
    assert sha(ROOT / "research/notes/ecc2k130/compact_orbit_point_panel_20260929/verify_panel.py") == FROZEN_TARGET_SHA
    MISMATCHES.clear()
    frozen_cold.check_target = check_any_pairing
    receipt = frozen_cold.verify("n37_L1024", run_dir.resolve(), relocated=True)
    assert receipt["status"] == "PASS"
    assert not receipt["host_timing_candidate"]
    assert len(receipt["checks"]) == 15
    assert sum(row["target_logs_verified"] for row in receipt["checks"]) == 15360
    return {"schema": "ecc2k130-disjoint-cold-pairing-addendum-v1",
            "status": "CORRECTED_REPLAY_PASS_NOT_TIMING_ELIGIBLE",
            "cell": "n37_L1024", "frozen_cold_verifier_sha256": FROZEN_COLD_SHA,
            "frozen_target_verifier_sha256": FROZEN_TARGET_SHA,
            "cold_run_sha256": sha(run_dir / "cold_run.json"),
            "original_hosted_failure_sha256": sha(run_dir / "receipt.json"),
            "verified_children": 15, "verified_target_logs": 15360,
            "pairing_mismatch_count": len(MISMATCHES),
            "pairing_mismatches": MISMATCHES,
            "corrected_core_receipt": receipt,
            "timing_eligible": False}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite a pairing addendum"
    try:
        result = replay(args.run_dir)
    except BaseException as error:
        result = {"status": "FAIL", "error_type": type(error).__name__,
                  "error": str(error), "traceback": traceback.format_exc()}
        args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        raise
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: result[key] for key in ("status", "cell", "verified_children",
                                               "verified_target_logs", "pairing_mismatch_count")},
                     sort_keys=True))


if __name__ == "__main__":
    main()
