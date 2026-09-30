#!/usr/bin/env python3
"""Freeze six disjoint public-Q streams and separate verifier-only labels."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"))
from verify_rank import Curve, Field  # noqa: E402

SOURCE_COMMIT = "d8d67d3b683334db105195b0fc74d8e660dc4b5e"
COMPACT_SHA = "7b2393731cd7da047b9404585c85dc2f384ee857e43778661e806649c6cf552d"
RHO_SHA = "98b6e8a27d821ebfc7ae716410184dffe89aadf27dfe36bcac8ce834c39cf04c"
SEED = 2026093007
CELLS = ((37, 1, 7, "off", 20), (37, 1024, 42, "off", 5),
         (41, 1, 85, "off", 5), (41, 1024, 255, "blocked", 5),
         (53, 1, 220, "off", 5), (53, 1024, 440, "blocked", 5))


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha(path: Path) -> str:
    return sha_bytes(path.read_bytes())


def corpus(n: int, length: int) -> str:
    return f"point-sum-cold-n{n}-L{length}-20260930-v1"


def prior_points() -> tuple[dict[int, set[tuple[int, int]]], list[dict]]:
    candidates = sorted((ROOT / "research/notes/ecc2k130").glob("**/*.points.jsonl"))
    assert len(candidates) == 25, "published point inventory changed; refreeze protocol"
    prior = {37: set(), 41: set(), 53: set()}
    inventory = []
    for path in candidates:
        matching = [n for n in prior if f"n{n}_" in path.name]
        assert len(matching) == 1, path
        n = matching[0]
        rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
        assert rows and all(isinstance(row, list) and len(row) == 2 for row in rows)
        prior[n].update(tuple(row) for row in rows)
        inventory.append({"path": str(path.relative_to(ROOT)),
                          "sha256": sha(path), "n": n, "rows": len(rows)})
    return prior, inventory


def generate(rho: Path, out: Path) -> dict:
    assert not (out / "FROZEN.json").exists(), "refusing to overwrite a frozen target corpus"
    assert not (out / "fixtures").exists(), "refusing to overwrite target files"
    assert sha(ROOT / "examples/koblitz_orbit_dlp_s3_batch.rs") == COMPACT_SHA
    assert sha(ROOT / "examples/koblitz_rho_batch_ks_v3.rs") == RHO_SHA
    assert (ROOT / "research/notes/ecc2k130/point_sum_cold_20260930/PROTOCOL.md").is_file()
    prior, inventory = prior_points()
    all_new = {37: set(), 41: set(), 53: set()}
    pending = []
    specs = {}
    for n, length, k, prefilter, blocks in CELLS:
        name = f"n{n}_L{length}"
        named_corpus = corpus(n, length)
        env = {key: value for key, value in os.environ.items() if not key.startswith("KIC_")}
        env.update({"KIC_RHO_BATCH_CORPUS": named_corpus,
                    "KIC_RHO_GENERATE_ONLY": "1",
                    "RAYON_NUM_THREADS": "1", "LC_ALL": "C"})
        command = [str(rho.resolve()), str(n), "0", "signed_frobenius",
                   str(length), str(SEED)]
        result = subprocess.run(command, env=env, capture_output=True, text=True, check=True)
        assert not result.stderr.strip(), result.stderr
        rows = [json.loads(line) for line in result.stdout.splitlines() if line.strip()]
        assert len(rows) == length
        first = rows[0]
        field = Field(n, first["field_modulus_low_terms"])
        curve = Curve(field, 0)
        generator = tuple(first["generator"])
        r = first["subgroup_order"]
        assert curve.on_curve(generator) and curve.mul(r, generator) is None
        assert first["automorphism_size"] == 2 * n
        points = []
        for index, row in enumerate(rows):
            assert (row["kind"], row["n"], row["a"], row["fixture_index"],
                    row["batch_seed"], row["corpus"]) == (
                "rho_ks_public_fixture", n, 0, index, SEED, named_corpus)
            assert all(row[key] == first[key] for key in (
                "generator", "subgroup_order", "automorphism_size", "field_modulus_low_terms"))
            scalar = row["published_fixture_scalar"]
            q = tuple(row["published_q"])
            assert 1 <= scalar < r and curve.on_curve(q)
            assert curve.mul(scalar, generator) == q
            assert q not in prior[n] and q not in all_new[n], (name, index)
            all_new[n].add(q)
            points.append(q)
        fixture_bytes = result.stdout.encode()
        point_bytes = "".join(json.dumps(q, separators=(",", ":")) + "\n"
                              for q in points).encode()
        fixture_rel = f"fixtures/{name}.fixture.jsonl"
        point_rel = f"fixtures/{name}.points.jsonl"
        pending.extend(((fixture_rel, fixture_bytes), (point_rel, point_bytes)))
        specs[name] = {
            "n": n, "a": 0, "L": length, "K": k, "blocks": blocks,
            "prefilter": prefilter, "corpus": named_corpus, "seed": SEED,
            "fixture_file": fixture_rel, "fixture_sha256": sha_bytes(fixture_bytes),
            "points_file": point_rel, "points_sha256": sha_bytes(point_bytes),
            "subgroup_order": r, "automorphism_size": 2 * n,
            "generator": first["generator"],
            "field_modulus_low_terms": first["field_modulus_low_terms"],
        }
    out.mkdir(parents=True)
    for relative, data in pending:
        path = out / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    frozen = {
        "schema": "ecc2k130-point-sum-cold-freeze-v1",
        "source_commit": SOURCE_COMMIT,
        "compact_source_sha256": COMPACT_SHA, "rho_source_sha256": RHO_SHA,
        "rho_generator_binary_sha256": sha(rho),
        "prepare_sha256": sha(Path(__file__)),
        "seed": SEED, "specs": specs,
        "prior_point_inventory": inventory,
        "prior_unique_counts": {str(n): len(prior[n]) for n in prior},
        "new_unique_counts": {str(n): len(all_new[n]) for n in all_new},
        "generator_command": "KIC_RHO_GENERATE_ONLY=1 KIC_RHO_BATCH_CORPUS=<corpus> rho <n> 0 signed_frobenius <L> 2026093007",
    }
    (out / "FROZEN.json").write_text(json.dumps(frozen, sort_keys=True, indent=2) + "\n")
    return frozen


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rho", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=HERE)
    args = parser.parse_args()
    frozen = generate(args.rho, args.out)
    print(json.dumps({"status": "PASS", "specs": sorted(frozen["specs"]),
                      "new_unique_counts": frozen["new_unique_counts"]}, sort_keys=True))


if __name__ == "__main__":
    main()
