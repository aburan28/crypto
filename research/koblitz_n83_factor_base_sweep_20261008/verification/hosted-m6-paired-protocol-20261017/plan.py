#!/usr/bin/env python3
"""Freeze the next matched N83 4+2 Gray/hash natural-search gate."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
OLD = STUDY / "verification/k600-gray-hash-paired-protocol-20261014/plan.json"
SUPPORT = STUDY / "verification/cross-stream-support-20261017/screen.json"
CAPACITY = STUDY / "verification/paired-capacity-screen-20261013/screen.py"
SAVED = HERE / "plan.json"
HASHES = {
    "old_plan": "fbd3fb943bc0f181797d106c41ee9961b8d70768b5c65a5584002e3c1b5e488f",
    "support": "82ff0dcbd4dbb9fd7ee4269cade509abd1f2d1eefc973dcc64900bea3af36e9e",
}
BASE = "3888e2597fd3519807c9d722f49f83329c93ce5cb0929501cf7cb7cbc3cab1a0"
POOLS = {
    "gray_pool": ("966a0493c96674f0af6d184615ed05309d08fba5c7342957b4186cd50366dc47",
                  ["public_x_gray_prefix", 2026100801, 600], 97),
    "hash_pool": ("b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e",
                  ["public_x_hash", 2026100801, 600], 600),
}


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build() -> dict:
    require(sha256(OLD) == HASHES["old_plan"]
            and sha256(SUPPORT) == HASHES["support"],
            "frozen prior plan or support receipt differs")
    old = json.loads(OLD.read_text())
    support = json.loads(SUPPORT.read_text())
    require(old["schema"] == "n83.k600-gray-hash-paired-source-protocol/v1"
            and old["arm_order"] == ["gray_pool", "hash_pool"]
            and support["status"] == "PASS_EXACT_648_PAIR_SUPPORT_SCREEN"
            and support["pair_count"] == 648,
            "prior source gate identity differs")
    spec = importlib.util.spec_from_file_location("n83_capacity_for_m6", CAPACITY)
    require(spec is not None and spec.loader is not None, "capacity source missing")
    capacity = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(capacity)
    arms = []
    for role, prior in zip(old["arm_order"], old["arms"]):
        pool_hash, pool_key, residual = POOLS[role]
        require(prior["role"] == role and prior["base_object_blake3"] == BASE
                and prior["pool_object_blake3"] == pool_hash
                and prior["residual_signed_orbits"] == residual,
                "frozen archive pair differs")
        case = [row for row in support["cases"]
                if row["base"] == ["public_x_sequential", 2026100802, 600]
                and row["pool"] == pool_key]
        require(len(case) == 1 and case[0]["residual_orbits"] == residual
                and case[0]["current_worker_import_admissible"],
                "whole-subgroup support case differs")
        source_capacity = capacity.requirements(600, 600 + residual, 6, 2)
        require(source_capacity["current_worker_cap_pass"] is True
                and source_capacity["total_variables"] <= 150000
                and source_capacity["domain_clauses"] <= 12000000,
                "4+2 worker source capacity differs")
        arms.append({
            "role": role, "base_object_blake3": BASE,
            "pool_object_blake3": pool_hash,
            "base_compressed_sha256": prior["base_compressed_sha256"],
            "pool_compressed_sha256": prior["pool_compressed_sha256"],
            "base_s3_uri": prior["base_s3_uri"],
            "pool_s3_uri": prior["pool_s3_uri"],
            "shared_signed_orbits": prior["shared_signed_orbits"],
            "residual_signed_orbits": residual,
            "source_capacity": source_capacity,
            "point_multisets_4_plus_2": case[0]["shapes"]["6/2"]["unordered_multisets"],
            "worker_verify_argv": prior["worker_verify_argv"],
            "worker_search_argv": ["/worker", "run-cross-pilot", BASE,
                                   "/pilot-manifest.json", "/base.jsonl.gz",
                                   pool_hash, "/pilot-manifest.json", "/pool.jsonl.gz",
                                   "/corpus.json", "2300", "20000", "64", "4", "2",
                                   "/out/search"],
        })
    return {
        "schema": "n83.hosted-m6-gray-hash-paired-protocol/v1",
        "status": "FROZEN_BOUNDED_NATURAL_SEARCH",
        "source_receipt_sha256": HASHES,
        "curve_slug": old["curve_slug"],
        "public_target": old["public_target"],
        "first_trial_a": old["first_trial_a"],
        "first_trial_b": old["first_trial_b"],
        "arm_order": old["arm_order"],
        "arms": arms,
        "limits": {
            "base_summands": 4, "residual_summands": 2,
            "order_first_pair": True, "conflicts_per_arm": 20000,
            "max_models_per_arm": 64, "worker_wall_seconds_per_arm": 2300,
            "external_wall_seconds_per_arm": 2400,
            "max_variables": 150000, "max_domain_clauses": 12000000,
            "run_cpus": 1, "run_memory_bytes": 4294967296,
            "swap_bytes": 0, "job_cap_seconds": 7500,
            "conservative_full_charge_seconds": 7500,
        },
        "scope_limit": "paired one-query SAT access only; not relation yield, graph rank or total IC runtime",
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    encoded = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if args.write:
        with SAVED.open("x") as stream:
            stream.write(encoded)
    else:
        require(SAVED.read_text() == encoded, "saved M6 plan differs")
    print(json.dumps({"status": "PASS", "sha256": hashlib.sha256(encoded.encode()).hexdigest()},
                     sort_keys=True))


if __name__ == "__main__":
    main()
