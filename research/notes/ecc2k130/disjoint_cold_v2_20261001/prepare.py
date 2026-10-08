#!/usr/bin/env python3
"""Freeze target-blind orbit-disjoint public Q for the v2 cold panel."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"))
from verify_rank import Curve, Field  # noqa: E402

PREREG_COMMIT = "8b4da5f88f49bd47c3b962d6f38a5debf3a0e6ad"
PROTOCOL_SHA = "e8d65dbd1ac2796f3cdfd2f7a090c16788525fb50f2b7e0a91911779bdb1d3fa"
SOURCE_FREEZE = ROOT / "research/notes/ecc2k130/compact_s3_prefilter_20260930/FROZEN.json"
SOURCE_FREEZE_SHA = "3e9f67cc2cd6de5a8458badb3525983d561d9c3449c05118819fb424c5093b2b"
OLD_INPUT_FREEZE = ROOT / "research/notes/ecc2k130/disjoint_cold_q_20260930/FROZEN.json"
OLD_INPUT_FREEZE_SHA = "8a3d4bbcb3bf894f522f1395b99b09eed6a3266af6cf5109b472bb4079963601"
PRIOR_DIGEST = "5dc98a163504b7c24f1ec781c8fbc6fc80a875d300129babff067863a21947fb"
COMPACT_SHA = "702a0a05709bc14bc10bafbb08edb2a2f5f86794e970d84c317da3bc65cdcf38"
RHO_SHA = "98b6e8a27d821ebfc7ae716410184dffe89aadf27dfe36bcac8ce834c39cf04c"
LOCK_SHA = "7d671f48c2da93f133d98802e80f858d1d9ea3b86996f7037f758990e1566627"
RANK_VERIFIER_SHA = "254869e52bdbd417e14d74c87940293490160faaa7fbb45f4edcb6b773a4b0f6"
CELLS = ((37, 1, 7, "off", 20), (37, 1024, 42, "off", 5),
         (41, 1, 85, "off", 5), (41, 1024, 255, "blocked", 5),
         (53, 1, 220, "off", 5), (53, 1024, 440, "blocked", 5))
BASE_SEED = 2026100110000
DOMAIN = "ecc2k130-disjoint-cold-v2-20261001"


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha(path: Path) -> str:
    return sha_bytes(path.read_bytes())


def cell_name(n: int, length: int) -> str:
    return f"n{n}_L{length}"


def block_seed(cell_index: int, block: int) -> int:
    return BASE_SEED + 100 * cell_index + block


def block_corpus(n: int, length: int, block: int) -> str:
    return f"compact-disjoint-cold-v2-n{n}-L{length}-b{block:02d}-20261001"


def candidate(cell: str, block: int, index: int, order: int) -> tuple[int, str]:
    assert index >= 0 and order > 1
    message = f"{DOMAIN}|{cell}|{block}|{index}".encode("utf-8")
    digest = hashlib.sha256(message).digest()
    return int.from_bytes(digest, "big") % order, digest.hex()


def orbit_key(x: int, field: Field, n: int) -> int:
    """Sign shares x; Frobenius squares x, so the minimum is a quotient key."""
    values = []
    current = x
    for _ in range(n):
        values.append(current)
        current = field.mul(current, current)
    assert current == x
    return min(values)


def fixed_curves() -> dict[int, tuple[Curve, Field, tuple[int, int], int, list[int]]]:
    assert sha(OLD_INPUT_FREEZE) == OLD_INPUT_FREEZE_SHA
    old = json.loads(OLD_INPUT_FREEZE.read_text())
    curves = {}
    for n in (37, 41, 53):
        spec = old["specs"][cell_name(n, 1)]
        field = Field(n, spec["field_modulus_low_terms"])
        curve = Curve(field, 0)
        generator = tuple(spec["generator"])
        order = spec["subgroup_order"]
        assert curve.on_curve(generator) and curve.mul(order, generator) is None
        assert spec["automorphism_size"] == 2 * n
        curves[n] = (curve, field, generator, order, spec["field_modulus_low_terms"])
    return curves


def prior_inventory(curves: dict) -> tuple[dict[int, set[int]], list[dict], str, dict[str, int]]:
    old = json.loads(OLD_INPUT_FREEZE.read_text())
    inventory = [{key: row[key] for key in ("path", "sha256", "rows")}
                 for row in old["prior_point_inventory"]]
    for spec in old["specs"].values():
        for block in spec["block_specs"]:
            inventory.append({"path": str(OLD_INPUT_FREEZE.parent.relative_to(ROOT) /
                                           block["points_file"]),
                              "sha256": block["points_sha256"], "rows": spec["L"]})
    inventory.sort(key=lambda row: row["path"])
    tracked = subprocess.check_output(
        ["git", "ls-files", "research/notes/ecc2k130/**/*.points.jsonl"],
        cwd=ROOT, text=True).splitlines()
    assert [row["path"] for row in inventory] == tracked
    digest = sha_bytes(json.dumps(inventory, sort_keys=True,
                                  separators=(",", ":")).encode())
    assert len(inventory) == 76 and sum(row["rows"] for row in inventory) == 34858
    assert digest == PRIOR_DIGEST
    prior = {n: set() for n in curves}
    rows_by_n = {str(n): 0 for n in curves}
    for record in inventory:
        path = ROOT / record["path"]
        assert path.is_file() and sha(path) == record["sha256"], record["path"]
        name = path.name
        n = next((candidate_n for candidate_n in curves
                  if name.startswith(f"n{candidate_n}_")), None)
        assert n is not None, name
        data = [json.loads(line) for line in path.read_bytes().splitlines() if line.strip()]
        assert len(data) == record["rows"]
        curve, field, _generator, _order, _modulus = curves[n]
        for q in data:
            assert isinstance(q, list) and len(q) == 2
            assert curve.on_curve(tuple(q))
            prior[n].add(orbit_key(q[0], field, n))
        rows_by_n[str(n)] += len(data)
    assert rows_by_n == {"37": 8216, "41": 13321, "53": 13321}
    return prior, inventory, digest, rows_by_n


def source_info(source: Path) -> dict:
    assert sha(SOURCE_FREEZE) == SOURCE_FREEZE_SHA
    freeze = json.loads(SOURCE_FREEZE.read_text())
    assert freeze["source_sha256"]["examples/koblitz_orbit_dlp_s3_batch.rs"] == COMPACT_SHA
    assert freeze["source_sha256"]["examples/koblitz_rho_batch_ks_v3.rs"] == RHO_SHA
    assert freeze["source_sha256"]["research/notes/ecc2k130/compact_shared_log_20260925/Cargo.lock"] == LOCK_SHA
    for path, expected in freeze["source_sha256"].items():
        assert sha(source / path) == expected, path
    assert sha(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929/verify_rank.py") == RANK_VERIFIER_SHA
    return {"source_freeze_sha256": SOURCE_FREEZE_SHA,
            "compact_source_sha256": COMPACT_SHA,
            "rho_source_sha256": RHO_SHA,
            "source_lock_sha256": LOCK_SHA,
            "rustc_version_verbose": subprocess.check_output(
                ["rustc", "--version", "--verbose"], text=True).strip(),
            "cargo_version": subprocess.check_output(["cargo", "--version"], text=True).strip()}


def generate(source: Path) -> dict:
    assert sha(HERE / "PROTOCOL.md") == PROTOCOL_SHA
    prereg = subprocess.check_output(
        ["git", "show", f"{PREREG_COMMIT}:research/notes/ecc2k130/disjoint_cold_v2_20261001/PROTOCOL.md"],
        cwd=ROOT)
    assert sha_bytes(prereg) == PROTOCOL_SHA
    assert not (HERE / "FROZEN.json").exists() and not (HERE / "fixtures").exists()
    curves = fixed_curves()
    prior, inventory, prior_digest, rows_by_n = prior_inventory(curves)
    source = source_info(source)
    seen: dict[int, set[int]] = {n: set() for n in curves}
    specs, pending = {}, []
    for cell_index, (n, length, k, prefilter, blocks) in enumerate(CELLS):
        cell = cell_name(n, length)
        curve, field, generator, order, modulus = curves[n]
        block_specs = []
        for block in range(blocks):
            seed = block_seed(cell_index, block)
            corpus = block_corpus(n, length, block)
            rejected = {"zero_scalar": 0, "prior_orbit": 0, "new_orbit": 0}
            labels, points = [], []
            index = 0
            while len(labels) < length:
                scalar, digest = candidate(cell, block, index, order)
                if scalar == 0:
                    rejected["zero_scalar"] += 1
                else:
                    q = curve.mul(scalar, generator)
                    assert q is not None and curve.on_curve(q)
                    key = orbit_key(q[0], field, n)
                    if key in prior[n]:
                        rejected["prior_orbit"] += 1
                    elif key in seen[n]:
                        rejected["new_orbit"] += 1
                    else:
                        seen[n].add(key)
                        record = {"kind": "disjoint_hash_public_fixture",
                                  "n": n, "a": 0, "fixture_index": len(labels),
                                  "batch_seed": seed, "corpus": corpus,
                                  "candidate_index": index, "candidate_sha256": digest,
                                  "selection_domain": DOMAIN,
                                  "field_modulus_low_terms": modulus,
                                  "generator": list(generator),
                                  "subgroup_order": order,
                                  "automorphism_size": 2 * n,
                                  "published_fixture_scalar": scalar,
                                  "published_q": list(q)}
                        labels.append(record)
                        points.append(list(q))
                index += 1
            prefix = f"fixtures/{cell}_b{block:02d}"
            fixture_bytes = "".join(json.dumps(row, sort_keys=True,
                                                separators=(",", ":")) + "\n"
                                    for row in labels).encode()
            points_bytes = "".join(json.dumps(row, separators=(",", ":")) + "\n"
                                   for row in points).encode()
            fixture_file, points_file = prefix + ".fixture.jsonl", prefix + ".points.jsonl"
            pending.extend(((fixture_file, fixture_bytes), (points_file, points_bytes)))
            block_specs.append({"block": block, "seed": seed, "corpus": corpus,
                                "fixture_file": fixture_file,
                                "fixture_sha256": sha_bytes(fixture_bytes),
                                "points_file": points_file,
                                "points_sha256": sha_bytes(points_bytes),
                                "candidate_attempts": index,
                                "candidate_rejections": rejected})
        specs[cell] = {"n": n, "a": 0, "L": length, "K": k,
                       "prefilter": prefilter, "blocks": blocks,
                       "field_modulus_low_terms": modulus,
                       "generator": list(generator), "subgroup_order": order,
                       "automorphism_size": 2 * n, "block_specs": block_specs}
    frozen = {"schema": "ecc2k130-disjoint-cold-v2-freeze-v1",
              "preregistration_commit": PREREG_COMMIT,
              "protocol_sha256": PROTOCOL_SHA,
              "prepare_sha256": sha(Path(__file__)),
              "source": source,
              "prior_point_inventory": inventory,
              "prior_inventory_digest": prior_digest,
              "prior_rows_by_n": rows_by_n,
              "prior_orbits_by_n": {str(n): len(prior[n]) for n in curves},
              "new_orbits_by_n": {str(n): len(seen[n]) for n in curves},
              "candidate_hash_domain": DOMAIN,
              "specs": specs}
    (HERE / "fixtures").mkdir()
    for relative, data in pending:
        (HERE / relative).write_bytes(data)
    (HERE / "FROZEN.json").write_text(json.dumps(frozen, indent=2, sort_keys=True) + "\n")
    return frozen


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    args = parser.parse_args()
    frozen = generate(args.source_root.resolve())
    print(json.dumps({"status": "PASS", "cells": sorted(frozen["specs"]),
                      "new_orbits_by_n": frozen["new_orbits_by_n"]}, sort_keys=True))


if __name__ == "__main__":
    main()
