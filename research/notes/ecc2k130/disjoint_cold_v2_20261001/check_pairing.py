#!/usr/bin/env python3
"""Check compact target witnesses without assuming the quartet's pair order."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"))
from verify_rank import point  # noqa: E402

PAIRINGS = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))


def check_target(record: dict, fixture: dict, base: dict, logs: list[int],
                 curve, generator: tuple[int, int],
                 checked_labels: set[tuple[str, int]]) -> None:
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
    assert record["recovered_scalar"] == scalar
    indices, codes = record["point_indices"], record["x_codes"]
    assert isinstance(indices, list) and len(indices) == 4
    assert isinstance(codes, list) and len(codes) == 4
    base_points = base["factor_base_point_coordinates"]
    labels = base["factor_base_point_labels"]
    chosen = []
    recomputed_log = 0
    for index, code in zip(indices, codes):
        assert 0 <= index < len(base_points)
        selected = point(base_points[index])
        assert selected is not None and selected[0] == code
        assert curve.on_curve(selected)
        column, coefficient = labels[index]
        assert 0 <= column < len(logs)
        label_key = (base["base_hash"], index)
        if label_key not in checked_labels:
            representative = point(base["factor_base_representatives"][column])
            assert curve.mul(coefficient, representative) == selected
            checked_labels.add(label_key)
        recomputed_log = (recomputed_log + coefficient * logs[column]) % r
        chosen.append(selected)
    pinned = record["pinned_intermediates"]
    assert isinstance(pinned, list) and len(pinned) == 2
    matches = 0
    for pairing in PAIRINGS:
        left = curve.add(chosen[pairing[0][0]], chosen[pairing[0][1]])
        right = curve.add(chosen[pairing[1][0]], chosen[pairing[1][1]])
        if left is not None and right is not None and set(pinned) == {left[0], right[0]}:
            assert curve.add(left, right) == q
            matches += 1
    assert matches >= 1, (record["fixture_index"], indices, pinned)
    assert recomputed_log == scalar
