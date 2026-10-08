#!/usr/bin/env python3
"""Freeze the two source-bound N53 S3-query implementations and workload."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CRYPTO = ROOT
CONTROL = HERE
TARGET = HERE / "public_target.json"
SOURCE_WORKLOAD = HERE / "source_workload.json"
BINARY_DIR = Path("/Volumes/SSD990/crypto/target/release")
PREFIX = "IC1N53Ckb1fb25864PDP4rootRCguidedLAgaussTDdirectISO0h"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(record: dict) -> bytes:
    return json.dumps(record, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode("utf-8")


def write_immutable(path: Path, value: dict) -> None:
    data = json.dumps(value, indent=2, sort_keys=True) + "\n"
    if path.exists():
        if path.read_text() != data:
            raise RuntimeError(f"frozen artifact would change: {path}")
    else:
        path.write_text(data)


def main() -> None:
    original = json.loads((CONTROL / "control_manifest.json").read_text())
    workload = copy.deepcopy(json.loads(SOURCE_WORKLOAD.read_text())["record"])
    assert sha(HERE / "baseline.rs") == (
        "00c0801b2f3e5150237d61e72176465ae5cdeec5e2878b8053d279dbb2ad03fb"
    )
    assert sha(TARGET) == (
        "9c75214e1a276a921ff8223884adb101b358ae1010ccddb289a3e3b13113216d"
    )
    workload["resource_envelope"] = {
        "candidate_processes": 1,
        "relation_and_index_workers": 14,
        "maximum_process_threads_including_main": 15,
        "memory_limit_bytes": None,
        "rho_workers_if_paired": 14,
    }
    workload["cache_policy"] = (
        "one fresh process per run; factor base and S3 index are prepared "
        "before the online clock, but relation collection, matrix rank, "
        "target span, and scalar recovery occur inside the online interval; "
        "no cross-run shared cache"
    )
    workload["input_law"] = (
        "one fixed previously unseen public subgroup point supplied as coordinates; "
        "fixture scalar is absent from both IC inputs"
    )
    workload["comparison_question"] = (
        "same-point before/after S3 target-query implementation cost with 14 workers"
    )
    workload_hash = hashlib.sha256(canonical(workload)).hexdigest()
    workload_id = workload_hash[:12]
    write_immutable(HERE / "workload.json", {
        "identity_sha256": workload_hash,
        "record": workload,
        "workload_id": workload_id,
    })

    components = {
        "binary_curve_point_arithmetic": CRYPTO / "src/binary_ecc/curve.rs",
        "binary_field_arithmetic": CRYPTO / "src/cryptanalysis/semaev_decomp.rs",
        "fast_binary_curve_arithmetic": CRYPTO / "src/cryptanalysis/koblitz_fast_arith.rs",
        "koblitz_curve_and_subgroup_definition": CRYPTO / "src/cryptanalysis/koblitz_index_calculus.rs",
        "cargo_manifest": HERE / "Cargo.toml",
        "cargo_lock": HERE / "Cargo.lock",
        "dependency_cargo_manifest": CRYPTO / "Cargo.toml",
        "dependency_cargo_lock": CRYPTO / "Cargo.lock",
    }
    source_hashes = {name: sha(path) for name, path in components.items()}
    frozen = {}
    for variant in ("baseline", "candidate"):
        source = HERE / f"{variant}.rs"
        binary = BINARY_DIR / f"s3-{variant}-v4"
        if not binary.is_file():
            raise RuntimeError(f"build the release binary first: {binary}")
        source_hash = sha(source)
        binary_hash = sha(binary)
        record = copy.deepcopy(original["identity_record"])
        implementation = record["implementation"]
        implementation["algorithm_affecting_flags"]["s3_pair_shared_inversion"] = (
            variant == "candidate"
        )
        implementation["algorithm_affecting_flags"]["target_s3_call_count_policy"] = (
            "charge both S3 queries when a state is batched"
            if variant == "candidate" else "charge each S3 query as visited"
        )
        implementation["build_command"] = (
            "cargo build --release --offline --bins"
        )
        implementation["executable_sha256"] = binary_hash
        implementation["source_component_sha256"].update(source_hashes)
        implementation["source_component_sha256"]["slice_ic_example"] = source_hash
        record["point_decomposition"]["query_inversion_policy"] = (
            "share one Montgomery inversion between the two generic target S3 "
            "queries per indexed state; use the original scalar solver if either "
            "query is exceptional; preserve root and witness order"
            if variant == "candidate" else
            "solve each target S3 query separately in indexed root order"
        )
        record["point_decomposition"]["implementation_version"] = (
            original["identity_record"]["point_decomposition"]["implementation_version"]
            + ("; two-query shared inversion" if variant == "candidate"
               else "; isolated baseline rebuild")
        )
        record["point_decomposition"]["source_digest"] = source_hash
        for section in ("relation_collection", "relation_linear_algebra", "target_descent"):
            record[section]["source_digest"] = source_hash
        identity_hash = hashlib.sha256(canonical(record)).hexdigest()
        candidate_id = PREFIX + identity_hash[:12]
        manifest = {
            "candidate_id": candidate_id,
            "identity_sha256": identity_hash,
            "identity_record": record,
            "candidate_freeze": {
                "kind": "pre-run-frozen-candidate",
                "source_path": str(source),
                "source_sha256": source_hash,
                "binary_path": str(binary),
                "binary_sha256": binary_hash,
                "workload_path": str(HERE / "workload.json"),
                "workload_id": workload_id,
                "target_path": str(TARGET),
                "target_sha256": sha(TARGET),
                "resource_envelope": workload["resource_envelope"],
                "source_control_candidate_id": original["candidate_id"],
                "source_control_manifest_sha256": sha(CONTROL / "control_manifest.json"),
            },
        }
        path = HERE / f"{variant}_manifest.json"
        write_immutable(path, manifest)
        frozen[variant] = {"candidate_id": candidate_id, "manifest_sha256": sha(path)}
    write_immutable(HERE / "freeze_receipt.json", {
        "kind": "s3_pair_query_pre_run_freeze",
        "workload_id": workload_id,
        "workload_sha256": sha(HERE / "workload.json"),
        "target_sha256": sha(TARGET),
        "variants": frozen,
    })
    print(json.dumps({"workload_id": workload_id, "variants": frozen}, sort_keys=True))


if __name__ == "__main__":
    main()
