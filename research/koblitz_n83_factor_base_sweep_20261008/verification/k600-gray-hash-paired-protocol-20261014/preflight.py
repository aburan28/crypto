#!/usr/bin/env python3
"""Freeze the K600 Gray/hash matched-pair source gate; launch no solver."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
ROOT = STUDY.parents[1]
PILOT = STUDY / "pilot-01"
LIVE = HERE.parent / "pilot-pair-live-s3-20261014"
CAPACITY = HERE.parent / "paired-capacity-screen-20261013" / "screen.py"
SCREEN = HERE.parent / "pilot-pair-support-screen-20261014" / "screen.json"
PLAN = HERE / "plan.json"
PINNED_SHA256 = {
    "manifest": (PILOT / "manifest.json", "5bd862987844f709e03beeab9fd265e3d2024011bcadcb82f5af9c72018df6a6"),
    "corpus": (PILOT / "probe-corpus.json", "760ae5c2e88615152e2569653a0e0c0cf3b2f84691f2b0731247f9ee3759dcbd"),
    "live_audit": (LIVE / "audit.json", "0c281fcc873eb914686c6e4ef3ef2bba51347c9cd88ddaad117e3b0c3f83b78e"),
    "support_screen": (SCREEN, "2c2229b9bfc78935cbc6c3ca7c6612f8fade2a40a0e0592f93a99b855d41dd45"),
    "capacity_source": (CAPACITY, "116617e9de7552e48b548e632e80121508515e7079878305de0b975c8a924c5d"),
    "worker_source": (ROOT / "examples/koblitz_n83_mixed_orbit_worker.rs", "5f310f7495989a94004a6545bed5682d54d3bc4b4f41db543ed24408a1fa1b7b"),
    "cargo_lock": (ROOT / "Cargo.lock", "ac5cd7a3459be776180afe1eb2975c76b1a88f84bd1e99aa20475927ec260a90"),
}
BASE = "3888e2597fd3519807c9d722f49f83329c93ce5cb0929501cf7cb7cbc3cab1a0"
POOLS = {
    "gray_pool": "966a0493c96674f0af6d184615ed05309d08fba5c7342957b4186cd50366dc47",
    "hash_pool": "b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e",
}
SUBGROUP_ORDER = 2417851639230796216685689
PRIOR_CHARGE = 23100
CHARGE_PER_ARM = {"build": 600, "import": 360, "search": 1200, "cleanup": 240}


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_json(path: Path, expected_sha256: str) -> dict:
    require(sha256(path) == expected_sha256, f"pinned bytes changed: {path.name}")
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"JSON object required: {path.name}")
    return value


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"module unavailable: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_plan(archive_dir: Path) -> dict:
    for label, (path, expected) in PINNED_SHA256.items():
        require(sha256(path) == expected, f"{label} source bytes changed")
    audit_receipt = checked_json(*PINNED_SHA256["live_audit"])
    live_audit = load_module("n83_live_pair_audit", LIVE / "audit.py")
    require(live_audit.audit(archive_dir) == audit_receipt,
            "local archive replay differs from saved live-S3 receipt")
    capacity = load_module("n83_paired_capacity", CAPACITY)
    screen = checked_json(*PINNED_SHA256["support_screen"])
    corpus = checked_json(*PINNED_SHA256["corpus"])
    manifest = checked_json(*PINNED_SHA256["manifest"])
    rows = {row["compressed_blake3"]: row for row in manifest["bases"]}
    require(set((BASE, *POOLS.values())) <= rows.keys(), "selected objects missing from manifest")
    targets = [row for row in corpus["targets"] if row.get("a") == 0 and row.get("fixture") == 0]
    require(len(targets) == 1 and targets[0]["point"] == [
        "1834913121178197242856238", "8797964190732846213981377"
    ], "first public target changed")
    require(audit_receipt["status"] == "PASS_ARCHIVE_BYTES_AND_ORBIT_OVERLAP", "live audit status changed")
    require(audit_receipt["selected_best_total_runtime"] is None, "unexpected runtime selection")

    arms = []
    for role, pool_hash in POOLS.items():
        pool = audit_receipt["objects"][role]
        base = audit_receipt["objects"]["base"]
        overlap = audit_receipt["pairs"][role]
        require(base["compressed_blake3"] == BASE and pool["compressed_blake3"] == pool_hash,
                "audit object identity changed")
        require(rows[pool_hash]["s3_uri"] == pool["s3_uri"], "pool S3 URI mismatch")
        require(rows[BASE]["s3_uri"] == base["s3_uri"], "base S3 URI mismatch")
        require(overlap["shared_signed_orbits"] + overlap["residual_signed_orbits"] == 600,
                "pool orbit partition changed")
        expected_residual = 97 if role == "gray_pool" else 600
        require(overlap["residual_signed_orbits"] == expected_residual, "residual orbit count changed")
        source_capacity = capacity.requirements(600, 600 + expected_residual, 5, 2)
        require(source_capacity["current_worker_cap_pass"] is True, "worker capacity gate failed")
        matching = [row for row in screen["ordered_pairs"]
                    if row["base_object_blake3"] == BASE and row["pool_object_blake3"] == pool_hash]
        require(len(matching) == 1, "support row missing or duplicated")
        shape = [item for item in matching[0]["shape_bounds"]
                 if item["summands"] == 5 and item["residual_summands"] == 2]
        require(len(shape) == 1, "3+2 support shape missing")
        multiset_count = int(shape[0]["point_multisets_including_repetition"])
        rank_floor = (600 * (SUBGROUP_ORDER - 1) + 2 * multiset_count - 1) // (2 * multiset_count)
        require(rank_floor == (33975 if role == "gray_pool" else 889), "rank floor changed")
        arms.append({
            "role": role,
            "base_object_blake3": BASE,
            "pool_object_blake3": pool_hash,
            "base_s3_uri": base["s3_uri"],
            "pool_s3_uri": pool["s3_uri"],
            "base_compressed_sha256": base["compressed_sha256"],
            "pool_compressed_sha256": pool["compressed_sha256"],
            "shared_signed_orbits": overlap["shared_signed_orbits"],
            "residual_signed_orbits": expected_residual,
            "source_capacity": source_capacity,
            "point_multisets_3_plus_2": str(multiset_count),
            "conditional_half_K_row_query_floor": rank_floor,
            "worker_verify_argv": ["/worker", "verify-cross-pilot", BASE,
                                   "/pilot-manifest.json", "/base.jsonl.gz", pool_hash,
                                   "/pilot-manifest.json", "/pool.jsonl.gz"],
            "worker_search_argv": ["/worker", "run-cross-pilot", BASE,
                                   "/pilot-manifest.json", "/base.jsonl.gz", pool_hash,
                                   "/pilot-manifest.json", "/pool.jsonl.gz", "/corpus.json",
                                   "1100", "5000", "32", "3", "2", "/out/search"],
        })

    requested = 2 * sum(CHARGE_PER_ARM.values())
    require(requested == 4800, "paired charge arithmetic changed")
    return {
        "schema": "n83.k600-gray-hash-paired-source-protocol/v1",
        "status": "SOURCE_READY_AUTHORIZATION_PENDING",
        "study": "koblitz_n83_factor_base_sweep_20261008",
        "curve_slug": "icv1-f2m83-tm6151469093347-debefd74",
        "subgroup_order": str(SUBGROUP_ORDER),
        "source_sha256": {label: expected for label, (_, expected) in PINNED_SHA256.items()},
        "live_audit_source_sha256": audit_receipt["source_sha256"],
        "container_image": "sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31",
        "public_target": targets[0]["point"],
        "first_trial_a": "1833078656608906180176350",
        "first_trial_b": "163089659398503559150141",
        "worker_limits": {
            "base_summands": 3, "residual_summands": 2, "order_first_pair": True,
            "conflicts": 5000, "max_models": 32, "worker_wall_seconds": 1100,
            "max_variables": 150000, "max_domain_clauses": 12000000,
            "run_cpus": 1, "run_memory_bytes": 4294967296, "swap_bytes": 0,
            "offline_build_cpus": 4, "offline_build_memory_bytes": 8589934592,
        },
        "arm_order": ["gray_pool", "hash_pool"],
        "charge_per_arm_seconds": CHARGE_PER_ARM,
        "prior_conservative_charged_seconds": PRIOR_CHARGE,
        "additional_seconds_requested_not_authorized": requested,
        "hypothetical_cumulative_if_fully_authorized_and_charged": PRIOR_CHARGE + requested,
        "arms": arms,
        "new_n83_solver_search_executed": False,
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-dir", type=Path, required=True)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    plan = make_plan(args.archive_dir)
    encoded = (json.dumps(plan, indent=2, sort_keys=True) + "\n").encode()
    if args.write:
        PLAN.write_bytes(encoded)
    else:
        require(PLAN.read_bytes() == encoded, "saved paired source protocol differs from replay")
    print(json.dumps({"status": plan["status"], "arms": [
        {"role": row["role"], "residual_signed_orbits": row["residual_signed_orbits"]}
        for row in plan["arms"]
    ]}, sort_keys=True))


if __name__ == "__main__":
    main()
