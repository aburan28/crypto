#!/usr/bin/env python3
"""Independent target replay with separate selected and x-only witnesses.

The x-only S3 search pins two pair-sum x coordinates before the full-point
lift is chosen.  The selected lift proves the recovered log; a possibly
different sign lift of the same four x coordinates proves the pinned roots.
"""
from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"))
from verify_rank import Curve, point  # noqa: E402


def pinned_x_witness(curve: Curve, selected: list[tuple[int, int]],
                     q: tuple[int, int], pinned: list[int]) -> list[int]:
    """Return sign bits for a four-point lift with the pinned pair roots and Q.

    On this binary curve the two points with a given nonzero x are P and -P.
    Checking all 16 sign assignments is therefore an exact x-only witness
    check, independent of the solver's selected full-point lift.
    """
    assert len(selected) == 4
    assert isinstance(pinned, list) and len(pinned) == 2
    assert all(type(x) is int and 0 <= x < (1 << curve.f.n) for x in pinned)
    options = []
    for first, second in ((0, 1), (2, 3)):
        pair_options = []
        for sign_first in (0, 1):
            for sign_second in (0, 1):
                left = curve.neg(selected[first]) if sign_first else selected[first]
                right = curve.neg(selected[second]) if sign_second else selected[second]
                pair = curve.add(left, right)
                if pair is not None:
                    pair_options.append((pair, (sign_first, sign_second)))
        options.append(pair_options)
    for left, left_signs in options[0]:
        for right, right_signs in options[1]:
            if sorted((left[0], right[0])) != sorted(pinned):
                continue
            if curve.add(left, right) == q:
                return [*left_signs, *right_signs]
    raise AssertionError("pinned x roots have no four-point lift summing to Q")


def check_target(record: dict, fixture: dict, base: dict, logs: list[int],
                 curve: Curve, generator: tuple[int, int],
                 checked_labels: set[tuple[str, int]]) -> dict:
    """Replay a target; the caller must first verify the frozen fixture file."""
    q = tuple(fixture["published_q"])
    scalar = fixture["published_fixture_scalar"]
    r = base["subgroup_order"]
    assert record["kind"] == "compact_orbit_dlp_target"
    assert record["fixture_index"] == fixture["fixture_index"]
    assert record["n"] == fixture["n"] and record["a"] == fixture["a"]
    assert record["target"] == record["published_q"] == list(q)
    assert record["generator"] == list(generator)
    assert record["published_fixture_scalar"] is None
    assert record["recovered_matches_published"] is None
    assert record["target_generation_ms_excluded"] == 0
    assert record["exit_code"] == 0 and record["group_verified"] is True
    assert 1 <= scalar < r and record["recovered_scalar"] == scalar
    assert curve.on_curve(q)

    indices, codes = record["point_indices"], record["x_codes"]
    assert isinstance(indices, list) and len(indices) == 4
    assert isinstance(codes, list) and len(codes) == 4
    base_points = base["factor_base_point_coordinates"]
    labels = base["factor_base_point_labels"]
    selected = []
    recomputed_log = 0
    for index, code in zip(indices, codes):
        assert type(index) is int and 0 <= index < len(base_points)
        chosen = point(base_points[index])
        assert chosen is not None and chosen[0] == code
        assert curve.on_curve(chosen)
        column, coefficient = labels[index]
        assert 0 <= column < len(logs)
        label_key = (base["base_hash"], index)
        if label_key not in checked_labels:
            representative = point(base["factor_base_representatives"][column])
            assert curve.mul(coefficient, representative) == chosen
            checked_labels.add(label_key)
        recomputed_log = (recomputed_log + coefficient * logs[column]) % r
        selected.append(chosen)
    left = curve.add(selected[0], selected[1])
    right = curve.add(selected[2], selected[3])
    assert left is not None and right is not None
    assert curve.add(left, right) == q
    assert recomputed_log == scalar

    pinned = record["pinned_intermediates"]
    signs = pinned_x_witness(curve, selected, q, pinned)
    return {"selected_pair_x": [left[0], right[0]],
            "pinned_pair_x": pinned, "pinned_witness_signs": signs,
            "pinned_matches_selected": sorted((left[0], right[0])) == sorted(pinned)}
