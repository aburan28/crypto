#!/usr/bin/env python3
"""Independent group-law replay of every frozen disjoint public point."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"))
from verify_rank import Curve, Field  # noqa: E402
from prepare import (BASE_SEED, CELLS, COMPACT_SHA, LOCK_SHA, PREREG_COMMIT,
                     PROTOCOL_SHA, PRIOR_DIGEST,
                     PRIOR_FILES, PRIOR_ROWS, RHO_SHA, SOURCE_FREEZE,
                     SOURCE_FREEZE_SHA, block_corpus, block_seed)  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list:
    return [json.loads(line) for line in path.read_bytes().splitlines() if line.strip()]


def replay() -> dict:
    freeze_path = HERE / "FROZEN.json"
    frozen = json.loads(freeze_path.read_text())
    assert frozen["schema"] == "ecc2k130-disjoint-cold-q-freeze-v1"
    assert frozen["preregistration_commit"] == PREREG_COMMIT
    assert frozen["protocol_sha256"] == PROTOCOL_SHA == sha(HERE / "PROTOCOL.md")
    assert frozen["prepare_sha256"] == sha(HERE / "prepare.py")
    assert sha(SOURCE_FREEZE) == SOURCE_FREEZE_SHA
    source = frozen["source"]
    assert (source["source_freeze_sha256"], source["compact_source_sha256"],
            source["rho_source_sha256"], source["source_lock_sha256"]) == (
        SOURCE_FREEZE_SHA, COMPACT_SHA, RHO_SHA, LOCK_SHA)
    assert len(source["rho_generator_binary_sha256"]) == 64
    assert source["rustc_version_verbose"].startswith("rustc ")
    assert source["cargo_version"].startswith("cargo ")

    inventory = frozen["prior_point_inventory"]
    assert len(inventory) == PRIOR_FILES
    prior = {37: set(), 41: set(), 53: set()}
    canonical = []
    for item in inventory:
        path = ROOT / item["path"]
        assert path.is_file() and HERE not in path.parents
        match = re.match(r"n(37|41|53)_", path.name)
        assert match is not None
        n = int(match.group(1))
        assert item["n"] == n and sha(path) == item["sha256"]
        original = rows(path)
        assert len(original) == item["rows"] > 0
        assert all(isinstance(point, list) and len(point) == 2 for point in original)
        prior[n].update(tuple(point) for point in original)
        canonical.append({key: item[key] for key in ("path", "sha256", "rows")})
    assert canonical == sorted(canonical, key=lambda item: item["path"])
    digest = hashlib.sha256(json.dumps(canonical, sort_keys=True,
                                       separators=(",", ":")).encode()).hexdigest()
    assert digest == PRIOR_DIGEST == frozen["prior_inventory_digest"]
    assert sum(item["rows"] for item in inventory) == PRIOR_ROWS
    assert {str(n): len(prior[n]) for n in prior} == frozen["prior_unique_counts"]
    current = sorted((ROOT / "research/notes/ecc2k130").glob("**/*.points.jsonl"))
    current_prior = [str(path.relative_to(ROOT)) for path in current if HERE not in path.parents]
    assert current_prior == [item["path"] for item in inventory], "prior inventory changed"

    seen = {37: set(), 41: set(), 53: set()}
    checked = []
    assert len(frozen["specs"]) == len(CELLS)
    assert BASE_SEED == 2026093091000
    for cell_index, (n, length, k, prefilter, blocks) in enumerate(CELLS):
        cell = f"n{n}_L{length}"
        spec = frozen["specs"][cell]
        assert (spec["n"], spec["a"], spec["L"], spec["K"],
                spec["prefilter"], spec["blocks"]) == (n, 0, length, k, prefilter, blocks)
        assert spec["automorphism_size"] == 2 * n
        curve = Curve(Field(n, spec["field_modulus_low_terms"]), 0)
        generator = tuple(spec["generator"])
        order = spec["subgroup_order"]
        assert curve.on_curve(generator) and curve.mul(order, generator) is None
        assert len(spec["block_specs"]) == blocks
        for block, record in enumerate(spec["block_specs"]):
            seed = block_seed(cell_index, block)
            corpus = block_corpus(n, length, block)
            assert (record["block"], record["seed"], record["corpus"]) == (
                block, seed, corpus)
            expected_prefix = f"fixtures/{cell}_b{block:02d}"
            assert record["fixture_file"] == expected_prefix + ".fixture.jsonl"
            assert record["points_file"] == expected_prefix + ".points.jsonl"
            fixture_path = HERE / record["fixture_file"]
            points_path = HERE / record["points_file"]
            assert sha(fixture_path) == record["fixture_sha256"]
            assert sha(points_path) == record["points_sha256"]
            labels = rows(fixture_path)
            points = rows(points_path)
            assert len(labels) == len(points) == length
            for index, (label, point) in enumerate(zip(labels, points)):
                assert (label["kind"], label["n"], label["a"],
                        label["fixture_index"], label["batch_seed"], label["corpus"]) == (
                    "rho_ks_public_fixture", n, 0, index, seed, corpus)
                assert (label["field_modulus_low_terms"], label["generator"],
                        label["subgroup_order"], label["automorphism_size"]) == (
                    spec["field_modulus_low_terms"], spec["generator"], order, 2 * n)
                scalar = label["published_fixture_scalar"]
                q = tuple(point)
                assert isinstance(scalar, int) and 1 <= scalar < order
                assert isinstance(point, list) and len(point) == 2
                assert all(isinstance(x, int) and 0 <= x < 1 << n for x in point)
                assert label["published_q"] == point and curve.on_curve(q)
                assert curve.mul(scalar, generator) == q
                assert q not in prior[n] and q not in seen[n], (cell, block, index)
                seen[n].add(q)
            checked.append({"cell": cell, "block": block, "targets": length,
                            "points_sha256": record["points_sha256"],
                            "fixture_sha256": record["fixture_sha256"]})
    assert {str(n): len(seen[n]) for n in seen} == frozen["new_unique_counts"]
    assert sum(item["targets"] for item in checked) == 15390
    return {"status": "PASS", "schema": "ecc2k130-disjoint-cold-q-input-replay-v1",
            "frozen_sha256": sha(freeze_path),
            "point_equations_verified": 15390,
            "blocks": checked}


if __name__ == "__main__":
    print(json.dumps(replay(), sort_keys=True))
