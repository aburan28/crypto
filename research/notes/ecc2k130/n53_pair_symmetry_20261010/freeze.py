#!/usr/bin/env python3
"""Freeze source, binaries, base and candidate identities before Q generation."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PARENT = ROOT / "research/notes/ecc2k130/n53_fixed_base_rank_restarts_20261010"
BASE = PARENT / "inputs/base_n53_k220.jsonl"
SOURCE = ROOT / "examples/koblitz_orbit_dlp_fast_online.rs"
GENERATOR = ROOT / "examples/koblitz_rho_fixture.rs"
RHO_SOURCE = ROOT / "examples/koblitz_rho_batch_ks_strong_online.rs"
MAIN = "9be831f4ec75c9a18334a67033333aeeb4b42fad"
BASE_HASH = "7af2460c8b5a2c29f9d1aa7fefecbcde3a6ce761293d0dab0cecfc3bc980b973"
PREFIX = "IC1N53Ce0fb23320PDP4rootRCguidedLAgaussTDdirectISO0h"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def save_new_or_equal(path: Path, value: object) -> None:
    encoded = json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False) + "\n"
    if path.exists():
        assert path.read_text() == encoded, f"frozen content changed: {path}"
    else:
        path.write_text(encoded)


def replace_source(value: object, old: str, new: str) -> object:
    if isinstance(value, dict):
        return {key: replace_source(item, old, new) for key, item in value.items()}
    if isinstance(value, list):
        return [replace_source(item, old, new) for item in value]
    return new if value == old else value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ic-binary", type=Path, default=ROOT / "target/release/examples/koblitz_orbit_dlp_fast_online")
    parser.add_argument("--generator-binary", type=Path, default=ROOT / "target/release/examples/koblitz_rho_fixture")
    parser.add_argument("--rho-binary", type=Path, default=ROOT / "target/release/examples/koblitz_rho_batch_ks_strong_online")
    args = parser.parse_args()
    assert subprocess.check_output(["git", "rev-parse", "origin/main"], cwd=ROOT).decode().strip() == MAIN
    parent = json.loads((PARENT / "candidate_control.json").read_text())
    parent_record = parent["record"]
    old_source = parent_record["implementation"]["example_source_sha256"]
    source = sha(SOURCE)
    base = json.loads(BASE.read_text())
    assert (base["n"], base["a"], base["factor_base_points"], base["orbit_columns"]) == (53, 0, 23320, 220)
    assert base["base_hash"] == BASE_HASH
    assert sha(BASE) == parent_record["factor_base"]["enumerated_set_sha256"]
    candidate_ids = {}
    for label, enabled in (("control", False), ("symmetry", True)):
        record = replace_source(copy.deepcopy(parent_record), old_source, source)
        record["implementation"]["library_parent_commit"] = MAIN
        record["implementation"]["flags"]["KIC_PAIR_SYMMETRY"] = int(enabled)
        record["point_decomposition"]["pair_symmetry"] = enabled
        record["point_decomposition"]["equation_order"] = (
            "full ordered left/right/relative/shift/root scan with later swapped pair skipped"
            if enabled else "full ordered left/right/relative/shift/root scan"
        )
        record["point_decomposition"]["index_construction"] = (
            "one S3 solve per unordered pair orbit; restore swapped aliases in original scan slots"
            if enabled else "one S3 solve per ordered pair/relative shift"
        )
        record["implementation"]["cargo_lock_sha256"] = sha(ROOT / "Cargo.lock")
        candidate_id = PREFIX + hashlib.sha256(canonical(record)).hexdigest()[:12]
        candidate_ids[label] = candidate_id
        save_new_or_equal(HERE / f"candidate_{label}.json", {"candidate_id": candidate_id, "record": record})
    frozen = {
        "schema": "n53-pair-symmetry-prep-v1",
        "base_hash_blake3": BASE_HASH,
        "base_file_sha256": sha(BASE),
        "curve_id": "EC1N53Ce0hb097de99be9a",
        "source_parent_commit": MAIN,
        "source_sha256": source,
        "ic_binary_sha256": sha(args.ic_binary),
        "generator_source_sha256": sha(GENERATOR),
        "generator_binary_sha256": sha(args.generator_binary),
        "rho_source_sha256": sha(RHO_SOURCE),
        "rho_binary_sha256": sha(args.rho_binary),
        "cargo_lock_sha256": sha(ROOT / "Cargo.lock"),
        "protocol_sha256": sha(HERE / "PROTOCOL.md"),
        "freeze_source_sha256": sha(Path(__file__)),
        "candidate_ids": candidate_ids,
        "public_hash_seed": 53271307,
        "generator_seed": 700001,
        "generator_args": ["53", "0", "signed_frobenius", "1", "strong", "700001", "hash:53271307"],
        "rank_seeds": list(range(532053, 532059)),
        "rho_seeds": list(range(531153, 531159)),
        "rho_config": {"rung": 3, "lanes": 32, "distinguished_bits": 4, "jump_count": 32, "target_count": 1},
        "resource_envelope": {"threads": 1, "wall_cap_seconds": 60, "observed_rss_cap_bytes": 17179869184},
        "timing_class": "exploratory_shared_host",
    }
    save_new_or_equal(HERE / "FROZEN_PREP.json", frozen)
    print(json.dumps({"status": "FROZEN", "candidate_ids": candidate_ids, "source_sha256": source}, sort_keys=True))


if __name__ == "__main__":
    main()
