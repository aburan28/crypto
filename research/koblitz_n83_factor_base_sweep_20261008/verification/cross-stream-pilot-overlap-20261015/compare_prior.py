#!/usr/bin/env python3
"""Compare every archived-pair overlap with the earlier full N83 census."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / "pilot-overlap-census-20261013/overlap.json"
INPUTS = {
    "prior": (PRIOR, "bcc6b479909fe75992252f50112d284159cd61a3392def648d3100c7e4869a77"),
    "replication": (HERE / "overlap.json", "f395be7b22c6d5a05bf1bc6c725ebc927420bacbc2408ac8d9169fe6ed0e4ea2"),
    "independent_audit": (HERE / "audit.json", "8ee5622cb6da482993d7cd7b2db3e9c5102367dd71805a6f061167758f2c8235"),
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def prior_key(row: dict) -> tuple:
    return (row["base_policy"], row["base_seed"], row["base_orbits"],
            row["pool_policy"], row["pool_seed"], row["pool_orbits"])


def replicate_key(row: dict) -> tuple:
    return (*row["base"], *row["pool"])


def compare() -> dict:
    sources = {}
    for name, (path, digest) in INPUTS.items():
        raw = path.read_bytes()
        require(sha(raw) == digest, f"{name} receipt changed")
        sources[name] = json.loads(raw)
    prior, current, audited = (sources[name] for name in INPUTS)
    require(prior["schema"] == "n83.primary-pilot-exact-overlap-census/v1"
            and prior["status"] == "PASS_archive_backed_exact_orbit_overlap"
            and current["status"] == "PASS_27_ARCHIVES_648_EXACT_PAIRS"
            and audited["status"] == "PASS_27_ARCHIVES_648_PAIRS_MINIMUM_POINT_REPLAY"
            and audited["receipt_sha256"] == INPUTS["replication"][1]
            and prior["ordered_cross_stream_pair_count"] == current["pair_count"] == 648,
            "input status or population changed")
    old = {prior_key(row): row for row in prior["ordered_pairs"]}
    new = {replicate_key(row): row for row in current["pairs"]}
    require(len(old) == len(new) == 648 and set(old) == set(new),
            "ordered pair keys differ")
    old_objects = {row["s3_uri"]: row for row in prior["archive_objects"]}
    new_objects = {row["s3_uri"]: row for row in current["archives"]}
    require(len(old_objects) == len(new_objects) == 27 and set(old_objects) == set(new_objects),
            "archive identity population differs")
    for uri, row in new_objects.items():
        ancestor = old_objects[uri]
        require(row["compressed_blake3"] == ancestor["compressed_blake3"]
                and (row["policy"], row["seed"], row["columns"])
                == (ancestor["policy"], ancestor["seed"], ancestor["orbit_columns"]),
                f"archived object changed: {uri}")
    for key, row in new.items():
        ancestor = old[key]
        require(row["shared_orbits"] == ancestor["shared_signed_orbits"]
                and row["residual_orbits"] == ancestor["residual_signed_orbits"]
                and ancestor["shared_x_values"] == 83 * row["shared_orbits"]
                and row["base_only_orbits"] == key[2] - row["shared_orbits"]
                and row["residual_points"] == 166 * ancestor["residual_signed_orbits"],
                f"prior/new pair differs: {key}")
    disjoint = sum(row["shared_orbits"] == 0 for row in new.values())
    empty = sum(row["residual_orbits"] == 0 for row in new.values())
    singleton = sum(row["residual_orbits"] == 1 for row in new.values())
    require((disjoint, empty, singleton)
            == (prior["disjoint_pair_count"], prior["fully_absorbed_pool_count"],
                prior["singleton_residual_count"]) == (378, 33, 57)
            and 648 - empty - singleton == prior["passes_current_min_two_residual_screen_count"] == 558,
            "summary dispositions changed")
    return {"schema": "n83.cross-stream-overlap-prior-comparison/v1",
            "status": "PASS_EXACT_MATCH_ALL_27_ARCHIVES_AND_648_PAIRS",
            "input_sha256": {name: digest for name, (_, digest) in INPUTS.items()},
            "generator_sha256": sha(Path(__file__).read_bytes()),
            "archive_matches": 27, "pair_matches": 648,
            "disjoint_pairs": disjoint, "empty_residual_pairs": empty,
            "singleton_residual_pairs": singleton,
            "current_worker_min_two_residual_pairs": 558,
            "selected_best_total_runtime": None}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = compare()
    data = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode()
    target = HERE / "comparison.json"
    if args.write:
        with target.open("xb") as stream:
            stream.write(data)
    else:
        require(target.read_bytes() == data, "comparison receipt differs on replay")
    print(json.dumps({"status": result["status"], "pairs": 648,
                      "sha256": sha(data)}))


if __name__ == "__main__":
    main()
