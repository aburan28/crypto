#!/usr/bin/env python3
"""Freeze disjoint public Q per paired block before any cold measurement."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"))
from verify_rank import Curve, Field  # noqa: E402

SOURCE_FREEZE = ROOT / "research/notes/ecc2k130/compact_s3_prefilter_20260930/FROZEN.json"
SOURCE_FREEZE_SHA = "3e9f67cc2cd6de5a8458badb3525983d561d9c3449c05118819fb424c5093b2b"
PREREG_COMMIT = "758a26ce13c80aa8e43cf2a2834d5428b728efff"
PROTOCOL_SHA = "59b039cc3af63d69f770b39a69ced0df1c8c7440dab36d45d906d64ed5a0d0f1"
COMPACT_SHA = "702a0a05709bc14bc10bafbb08edb2a2f5f86794e970d84c317da3bc65cdcf38"
RHO_SHA = "98b6e8a27d821ebfc7ae716410184dffe89aadf27dfe36bcac8ce834c39cf04c"
LOCK_SHA = "7d671f48c2da93f133d98802e80f858d1d9ea3b86996f7037f758990e1566627"
RANK_VERIFIER_SHA = "254869e52bdbd417e14d74c87940293490160faaa7fbb45f4edcb6b773a4b0f6"
PRIOR_FILES = 31
PRIOR_ROWS = 19468
PRIOR_DIGEST = "30a0b11482bad46c6efe994278001941a2e59c84cf10d366ec433b785637c9cd"
CELLS = (
    (37, 1, 7, "off", 20),
    (37, 1024, 42, "off", 5),
    (41, 1, 85, "off", 5),
    (41, 1024, 255, "blocked", 5),
    (53, 1, 220, "off", 5),
    (53, 1024, 440, "blocked", 5),
)
BASE_SEED = 2026093091000


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha(path: Path) -> str:
    return sha_bytes(path.read_bytes())


def cell_name(n: int, length: int) -> str:
    return f"n{n}_L{length}"


def block_seed(cell_index: int, block: int) -> int:
    return BASE_SEED + 100 * cell_index + block


def block_corpus(n: int, length: int, block: int) -> str:
    return f"compact-disjoint-cold-n{n}-L{length}-b{block:02d}-20260930-v1"


def prior_points() -> tuple[dict[int, set[tuple[int, int]]], list[dict], str]:
    files = sorted((ROOT / "research/notes/ecc2k130").glob("**/*.points.jsonl"))
    files = [path for path in files if HERE not in path.parents]
    inventory = []
    prior: dict[int, set[tuple[int, int]]] = {37: set(), 41: set(), 53: set()}
    digest_rows = []
    for path in files:
        match = re.match(r"n(37|41|53)_", path.name)
        assert match is not None, path
        n = int(match.group(1))
        data = path.read_bytes()
        rows = [json.loads(line) for line in data.splitlines() if line.strip()]
        assert rows and all(isinstance(q, list) and len(q) == 2 and
                            all(isinstance(x, int) and 0 <= x < 1 << n for x in q)
                            for q in rows), path
        relative = str(path.relative_to(ROOT))
        record = {"path": relative, "sha256": sha_bytes(data), "rows": len(rows)}
        digest_rows.append(record)
        inventory.append({**record, "n": n})
        prior[n].update(tuple(q) for q in rows)
    digest = sha_bytes(json.dumps(digest_rows, sort_keys=True,
                                  separators=(",", ":")).encode())
    assert (len(files), sum(row["rows"] for row in inventory), digest) == (
        PRIOR_FILES, PRIOR_ROWS, PRIOR_DIGEST), "published point inventory changed; refreeze protocol"
    return prior, inventory, digest


def check_source(source: Path, rho: Path) -> dict:
    assert sha(SOURCE_FREEZE) == SOURCE_FREEZE_SHA
    freeze = json.loads(SOURCE_FREEZE.read_text())
    assert freeze["source_sha256"]["examples/koblitz_orbit_dlp_s3_batch.rs"] == COMPACT_SHA
    assert freeze["source_sha256"]["examples/koblitz_rho_batch_ks_v3.rs"] == RHO_SHA
    assert freeze["source_sha256"]["research/notes/ecc2k130/compact_shared_log_20260925/Cargo.lock"] == LOCK_SHA
    assert sha(source / "examples/koblitz_orbit_dlp_s3_batch.rs") == COMPACT_SHA
    assert sha(source / "examples/koblitz_rho_batch_ks_v3.rs") == RHO_SHA
    assert sha(source / "research/notes/ecc2k130/compact_shared_log_20260925/Cargo.lock") == LOCK_SHA
    assert sha(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929/verify_rank.py") == RANK_VERIFIER_SHA
    assert rho.is_file()
    return {"source_freeze_sha256": SOURCE_FREEZE_SHA,
            "compact_source_sha256": COMPACT_SHA,
            "rho_source_sha256": RHO_SHA,
            "source_lock_sha256": LOCK_SHA,
            "rho_generator_binary_sha256": sha(rho),
            "rustc_version_verbose": subprocess.check_output(
                ["rustc", "--version", "--verbose"], text=True).strip(),
            "cargo_version": subprocess.check_output(["cargo", "--version"], text=True).strip()}


def generate(source: Path, rho: Path, out: Path) -> dict:
    assert out.resolve() == HERE.resolve(), "this protocol freezes inputs in its versioned directory"
    assert sha(HERE / "PROTOCOL.md") == PROTOCOL_SHA
    prereg = subprocess.check_output(
        ["git", "show", f"{PREREG_COMMIT}:research/notes/ecc2k130/disjoint_cold_q_20260930/PROTOCOL.md"],
        cwd=ROOT)
    assert sha_bytes(prereg) == PROTOCOL_SHA, "preregistration commit missing or changed"
    assert not (out / "FROZEN.json").exists(), "refusing to overwrite a target freeze"
    assert not (out / "fixtures").exists(), "refusing to overwrite target files"
    source_info = check_source(source, rho)
    prior, inventory, inventory_digest = prior_points()
    new: dict[int, set[tuple[int, int]]] = {37: set(), 41: set(), 53: set()}
    expected_by_n: dict[int, tuple] = {}
    pending: list[tuple[str, bytes]] = []
    specs: dict[str, dict] = {}
    for cell_index, (n, length, k, prefilter, blocks) in enumerate(CELLS):
        cell = cell_name(n, length)
        block_specs = []
        for block in range(blocks):
            seed = block_seed(cell_index, block)
            corpus = block_corpus(n, length, block)
            env = {key: value for key, value in os.environ.items() if not key.startswith("KIC_")}
            env.update({"KIC_RHO_GENERATE_ONLY": "1", "KIC_RHO_BATCH_CORPUS": corpus,
                        "RAYON_NUM_THREADS": "1", "LC_ALL": "C"})
            command = [str(rho.resolve()), str(n), "0", "signed_frobenius",
                       str(length), str(seed)]
            result = subprocess.run(command, env=env, capture_output=True, check=True)
            assert not result.stderr.strip(), (cell, block, result.stderr.decode(errors="replace"))
            labels = [json.loads(line) for line in result.stdout.splitlines() if line.strip()]
            assert len(labels) == length, (cell, block)
            first = labels[0]
            identity = (first["field_modulus_low_terms"], first["generator"],
                        first["subgroup_order"], first["automorphism_size"])
            if n in expected_by_n:
                assert identity == expected_by_n[n], (cell, block, "generator changed")
            expected_by_n[n] = identity
            field = Field(n, identity[0])
            curve = Curve(field, 0)
            generator = tuple(identity[1])
            order = identity[2]
            assert curve.on_curve(generator) and curve.mul(order, generator) is None
            assert identity[3] == 2 * n
            points: list[list[int]] = []
            for index, record in enumerate(labels):
                assert (record["kind"], record["n"], record["a"],
                        record["fixture_index"], record["batch_seed"], record["corpus"]) == (
                    "rho_ks_public_fixture", n, 0, index, seed, corpus)
                assert (record["field_modulus_low_terms"], record["generator"],
                        record["subgroup_order"], record["automorphism_size"]) == identity
                scalar = record["published_fixture_scalar"]
                point = record["published_q"]
                q = tuple(point)
                assert isinstance(scalar, int) and 1 <= scalar < order
                assert isinstance(point, list) and len(point) == 2
                assert all(isinstance(x, int) and 0 <= x < 1 << n for x in point)
                assert curve.on_curve(q) and curve.mul(scalar, generator) == q
                assert q not in prior[n] and q not in new[n], (cell, block, index, "Q collision")
                new[n].add(q)
                points.append(point)
            prefix = f"fixtures/{cell}_b{block:02d}"
            fixture_file = prefix + ".fixture.jsonl"
            points_file = prefix + ".points.jsonl"
            fixture_bytes = result.stdout
            points_bytes = "".join(json.dumps(q, separators=(",", ":")) + "\n"
                                   for q in points).encode()
            pending.extend(((fixture_file, fixture_bytes), (points_file, points_bytes)))
            block_specs.append({"block": block, "seed": seed, "corpus": corpus,
                                "fixture_file": fixture_file,
                                "fixture_sha256": sha_bytes(fixture_bytes),
                                "points_file": points_file,
                                "points_sha256": sha_bytes(points_bytes)})
        specs[cell] = {"n": n, "a": 0, "L": length, "K": k,
                       "prefilter": prefilter, "blocks": blocks,
                       "field_modulus_low_terms": expected_by_n[n][0],
                       "generator": expected_by_n[n][1],
                       "subgroup_order": expected_by_n[n][2],
                       "automorphism_size": expected_by_n[n][3],
                       "block_specs": block_specs}
    frozen = {"schema": "ecc2k130-disjoint-cold-q-freeze-v1",
              "preregistration_commit": PREREG_COMMIT,
              "protocol_sha256": PROTOCOL_SHA,
              "prepare_sha256": sha(Path(__file__)),
              "source": source_info,
              "prior_point_inventory": inventory,
              "prior_inventory_digest": inventory_digest,
              "prior_unique_counts": {str(n): len(prior[n]) for n in prior},
              "new_unique_counts": {str(n): len(new[n]) for n in new},
              "specs": specs,
              "generator_command": "KIC_RHO_GENERATE_ONLY=1 KIC_RHO_BATCH_CORPUS=<corpus> rho <n> 0 signed_frobenius <L> <block_seed>"}
    (out / "fixtures").mkdir()
    for relative, data in pending:
        (out / relative).write_bytes(data)
    (out / "FROZEN.json").write_text(json.dumps(frozen, indent=2, sort_keys=True) + "\n")
    return frozen


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--rho", type=Path, required=True)
    args = parser.parse_args()
    frozen = generate(args.source_root.resolve(), args.rho.resolve(), HERE)
    print(json.dumps({"status": "PASS", "cells": sorted(frozen["specs"]),
                      "new_unique_counts": frozen["new_unique_counts"]}, sort_keys=True))


if __name__ == "__main__":
    main()
