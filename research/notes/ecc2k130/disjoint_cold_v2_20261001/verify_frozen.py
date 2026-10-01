#!/usr/bin/env python3
"""Independently regenerate every selected Q and rejected orbit candidate."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"))
from verify_rank import Curve, Field  # noqa: E402
from prepare import (BASE_SEED, CELLS, COMPACT_SHA, DOMAIN, LOCK_SHA,
                     OLD_INPUT_FREEZE, OLD_INPUT_FREEZE_SHA, PREREG_COMMIT,
                     PRIOR_DIGEST, PROTOCOL_SHA, RHO_SHA, SOURCE_FREEZE,
                     SOURCE_FREEZE_SHA, block_corpus, block_seed, cell_name)  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list:
    return [json.loads(line) for line in path.read_bytes().splitlines() if line.strip()]


def canonical_x(x: int, field: Field, n: int) -> int:
    smallest = x
    current = x
    for _ in range(n):
        smallest = min(smallest, current)
        current = field.mul(current, current)
    assert current == x
    return smallest


def replay() -> dict:
    frozen_path = HERE / "FROZEN.json"
    frozen = json.loads(frozen_path.read_text())
    assert frozen["schema"] == "ecc2k130-disjoint-cold-v2-freeze-v1"
    assert frozen["preregistration_commit"] == PREREG_COMMIT
    assert frozen["protocol_sha256"] == PROTOCOL_SHA == sha(HERE / "PROTOCOL.md")
    assert frozen["prepare_sha256"] == sha(HERE / "prepare.py")
    assert frozen["candidate_hash_domain"] == DOMAIN
    assert sha(SOURCE_FREEZE) == SOURCE_FREEZE_SHA
    assert sha(OLD_INPUT_FREEZE) == OLD_INPUT_FREEZE_SHA
    assert frozen["source"]["source_freeze_sha256"] == SOURCE_FREEZE_SHA
    assert frozen["source"]["compact_source_sha256"] == COMPACT_SHA
    assert frozen["source"]["rho_source_sha256"] == RHO_SHA
    assert frozen["source"]["source_lock_sha256"] == LOCK_SHA
    assert frozen["source"]["rustc_version_verbose"].startswith("rustc ")
    assert frozen["source"]["cargo_version"].startswith("cargo ")

    inventory = frozen["prior_point_inventory"]
    assert len(inventory) == 76
    canonical = [{key: item[key] for key in ("path", "sha256", "rows")}
                 for item in inventory]
    assert canonical == sorted(canonical, key=lambda item: item["path"])
    digest = hashlib.sha256(json.dumps(canonical, sort_keys=True,
                                       separators=(",", ":")).encode()).hexdigest()
    assert digest == PRIOR_DIGEST == frozen["prior_inventory_digest"]
    tracked = [path for path in subprocess.check_output(
        ["git", "ls-files", "research/notes/ecc2k130/**/*.points.jsonl"],
        cwd=ROOT, text=True).splitlines()
        if not path.startswith(str(HERE.relative_to(ROOT)) + "/")]
    canonical_paths = [item["path"] for item in canonical]
    # The equality used to cover every ecc2k130 points file. Later freezes
    # add their own corpora beside this one. A new file inside a directory
    # this inventory already listed is still a hole; a file in a new
    # directory is a later experiment and does not change these hashes.
    historical_dirs = {str(Path(path).parent) for path in canonical_paths}
    historical_tracked = [path for path in tracked if str(Path(path).parent) in historical_dirs]
    assert canonical_paths == historical_tracked

    old = json.loads(OLD_INPUT_FREEZE.read_text())
    curves = {}
    prior = {}
    count_by_n = {"37": 0, "41": 0, "53": 0}
    for n in (37, 41, 53):
        spec = old["specs"][cell_name(n, 1)]
        field = Field(n, spec["field_modulus_low_terms"])
        curve = Curve(field, 0)
        generator = tuple(spec["generator"])
        order = spec["subgroup_order"]
        assert curve.on_curve(generator) and curve.mul(order, generator) is None
        curves[n] = (field, curve, generator, order, spec["field_modulus_low_terms"])
        prior[n] = set()
    for item in inventory:
        path = ROOT / item["path"]
        assert path.is_file() and HERE not in path.parents
        assert sha(path) == item["sha256"]
        n = next((value for value in curves if path.name.startswith(f"n{value}_")), None)
        assert n is not None
        entries = rows(path)
        assert len(entries) == item["rows"] > 0
        field, curve, _generator, _order, _modulus = curves[n]
        for q in entries:
            assert isinstance(q, list) and len(q) == 2
            assert curve.on_curve(tuple(q))
            prior[n].add(canonical_x(q[0], field, n))
        count_by_n[str(n)] += len(entries)
    assert count_by_n == frozen["prior_rows_by_n"] == {
        "37": 8216, "41": 13321, "53": 13321}
    assert {str(n): len(prior[n]) for n in prior} == frozen["prior_orbits_by_n"]

    assert set(frozen["specs"]) == {cell_name(n, length) for n, length, *_ in CELLS}
    assert BASE_SEED == 2026100110000
    seen = {n: set() for n in curves}
    checked = []
    rejects_total = {"zero_scalar": 0, "prior_orbit": 0, "new_orbit": 0}
    for cell_index, (n, length, k, prefilter, blocks) in enumerate(CELLS):
        cell = cell_name(n, length)
        spec = frozen["specs"][cell]
        field, curve, generator, order, modulus = curves[n]
        assert (spec["n"], spec["a"], spec["L"], spec["K"],
                spec["prefilter"], spec["blocks"]) == (n, 0, length, k, prefilter, blocks)
        assert (spec["field_modulus_low_terms"], spec["generator"],
                spec["subgroup_order"], spec["automorphism_size"]) == (
            modulus, list(generator), order, 2 * n)
        assert len(spec["block_specs"]) == blocks
        for block, record in enumerate(spec["block_specs"]):
            assert (record["block"], record["seed"], record["corpus"]) == (
                block, block_seed(cell_index, block), block_corpus(n, length, block))
            prefix = f"fixtures/{cell}_b{block:02d}"
            assert record["fixture_file"] == prefix + ".fixture.jsonl"
            assert record["points_file"] == prefix + ".points.jsonl"
            fixture_path = HERE / record["fixture_file"]
            points_path = HERE / record["points_file"]
            assert sha(fixture_path) == record["fixture_sha256"]
            assert sha(points_path) == record["points_sha256"]
            labels, points = rows(fixture_path), rows(points_path)
            assert len(labels) == len(points) == length
            rejected = {"zero_scalar": 0, "prior_orbit": 0, "new_orbit": 0}
            accepted = 0
            for candidate_index in range(record["candidate_attempts"]):
                text = f"{DOMAIN}|{cell}|{block}|{candidate_index}".encode("utf-8")
                digest = hashlib.sha256(text).digest()
                scalar = int.from_bytes(digest, "big") % order
                if scalar == 0:
                    rejected["zero_scalar"] += 1
                    continue
                q = curve.mul(scalar, generator)
                assert q is not None and curve.on_curve(q)
                key = canonical_x(q[0], field, n)
                if key in prior[n]:
                    rejected["prior_orbit"] += 1
                    continue
                if key in seen[n]:
                    rejected["new_orbit"] += 1
                    continue
                assert accepted < length
                seen[n].add(key)
                label, public = labels[accepted], points[accepted]
                assert label == {
                    "kind": "disjoint_hash_public_fixture", "n": n, "a": 0,
                    "fixture_index": accepted, "batch_seed": record["seed"],
                    "corpus": record["corpus"], "candidate_index": candidate_index,
                    "candidate_sha256": digest.hex(), "selection_domain": DOMAIN,
                    "field_modulus_low_terms": modulus, "generator": list(generator),
                    "subgroup_order": order, "automorphism_size": 2 * n,
                    "published_fixture_scalar": scalar, "published_q": list(q)}
                assert public == list(q)
                accepted += 1
            assert accepted == length and record["candidate_attempts"] > 0
            assert labels[-1]["candidate_index"] == record["candidate_attempts"] - 1
            assert rejected == record["candidate_rejections"]
            for reason in rejects_total:
                rejects_total[reason] += rejected[reason]
            checked.append({"cell": cell, "block": block, "targets": length,
                            "points_sha256": record["points_sha256"],
                            "fixture_sha256": record["fixture_sha256"],
                            "candidate_attempts": record["candidate_attempts"],
                            "candidate_rejections": rejected})
    assert {str(n): len(seen[n]) for n in seen} == frozen["new_orbits_by_n"]
    assert sum(row["targets"] for row in checked) == 15390
    return {"status": "PASS", "schema": "ecc2k130-disjoint-cold-v2-input-replay-v1",
            "frozen_sha256": sha(frozen_path), "point_equations_verified": 15390,
            "prior_orbits_by_n": {str(n): len(prior[n]) for n in prior},
            "new_orbits_by_n": {str(n): len(seen[n]) for n in seen},
            "candidate_rejections": rejects_total, "blocks": checked}


if __name__ == "__main__":
    print(json.dumps(replay(), sort_keys=True))
