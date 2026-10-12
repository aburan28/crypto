#!/usr/bin/env python3
"""Read-only independent replay of the six frozen public-point corpora."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"))
from verify_rank import Curve, Field  # noqa: E402
from prepare import CELLS, COMPACT_SHA, RHO_SHA, SEED, SOURCE_COMMIT, corpus  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def verify() -> dict:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    assert frozen["schema"] == "ecc2k130-point-sum-cold-freeze-v1"
    assert frozen["source_commit"] == SOURCE_COMMIT
    assert frozen["seed"] == SEED
    assert frozen["compact_source_sha256"] == COMPACT_SHA
    assert frozen["rho_source_sha256"] == RHO_SHA
    assert sha(ROOT / "examples/koblitz_orbit_dlp_s3_batch.rs") == COMPACT_SHA
    assert sha(ROOT / "examples/koblitz_rho_batch_ks_v3.rs") == RHO_SHA
    assert sha(HERE / "prepare.py") == frozen["prepare_sha256"]
    inventory = frozen["prior_point_inventory"]
    assert len(inventory) == 25
    prior = {37: set(), 41: set(), 53: set()}
    for item in inventory:
        path = ROOT / item["path"]
        assert sha(path) == item["sha256"]
        original = rows(path)
        assert len(original) == item["rows"] > 0
        assert item["n"] in prior
        prior[item["n"]].update(tuple(q) for q in original)
    assert {str(n): len(prior[n]) for n in prior} == frozen["prior_unique_counts"]
    new = {37: set(), 41: set(), 53: set()}
    checked = []
    assert len(frozen["specs"]) == len(CELLS)
    for n, length, k, prefilter, blocks in CELLS:
        name = f"n{n}_L{length}"
        spec = frozen["specs"][name]
        assert (spec["n"], spec["a"], spec["L"], spec["K"],
                spec["prefilter"], spec["blocks"], spec["seed"], spec["corpus"]) == (
            n, 0, length, k, prefilter, blocks, SEED, corpus(n, length))
        fixture_path = HERE / spec["fixture_file"]
        points_path = HERE / spec["points_file"]
        assert sha(fixture_path) == spec["fixture_sha256"]
        assert sha(points_path) == spec["points_sha256"]
        labels = rows(fixture_path)
        points = rows(points_path)
        assert len(labels) == len(points) == length
        curve = Curve(Field(n, spec["field_modulus_low_terms"]), 0)
        generator = tuple(spec["generator"])
        r = spec["subgroup_order"]
        assert curve.on_curve(generator) and curve.mul(r, generator) is None
        assert spec["automorphism_size"] == 2 * n
        for index, (record, point) in enumerate(zip(labels, points)):
            assert (record["kind"], record["fixture_index"], record["n"],
                    record["a"], record["corpus"], record["batch_seed"]) == (
                "rho_ks_public_fixture", index, n, 0, corpus(n, length), SEED)
            assert record["field_modulus_low_terms"] == spec["field_modulus_low_terms"]
            assert record["generator"] == spec["generator"]
            assert record["subgroup_order"] == r
            assert record["automorphism_size"] == 2 * n
            scalar = record["published_fixture_scalar"]
            q = tuple(point)
            assert record["published_q"] == point
            assert 1 <= scalar < r and curve.on_curve(q)
            assert curve.mul(scalar, generator) == q
            assert q not in prior[n] and q not in new[n]
            new[n].add(q)
        checked.append({"cell": name, "points": length,
                        "points_sha256": spec["points_sha256"],
                        "labels_sha256": spec["fixture_sha256"]})
    assert {str(n): len(new[n]) for n in new} == frozen["new_unique_counts"]
    return {"status": "PASS", "schema": "ecc2k130-point-sum-cold-input-replay-v1",
            "frozen_sha256": sha(HERE / "FROZEN.json"),
            "point_equations_verified": sum(item["points"] for item in checked),
            "cells": checked}


if __name__ == "__main__":
    print(json.dumps(verify(), sort_keys=True))
