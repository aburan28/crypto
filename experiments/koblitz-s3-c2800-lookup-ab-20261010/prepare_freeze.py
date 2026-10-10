"""Freeze three exact S3 lookup modes and a balanced T025 stage diagnostic."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_c2800lookup_ab.rs"
RUNNER = HERE / "run_diagnostic.py"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
LATE_MANIFEST = (ROOT / "experiments/koblitz-s3-c2800-prefetch-20261010/"
                 "identities/candidates/"
                 "IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0h9d7b111b5b0f.json")
REFERENCE = (ROOT / "experiments/koblitz-s3-c2800-prefetch-20261010/local/"
             "verified_t025/report.json.gz")
REFERENCE_RECEIPT = (ROOT / "experiments/koblitz-s3-c2800-prefetch-20261010/local/"
                     "verified_t025/receipt.json")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def replace(value: object, old: str, new: str) -> object:
    if isinstance(value, dict):
        return {key: replace(item, old, new) for key, item in value.items()}
    if isinstance(value, list):
        return [replace(item, old, new) for item in value]
    return new if value == old else value


def prepare(binary: Path) -> None:
    freeze_path = HERE / "DIAGNOSTIC_FREEZE.json"
    schedule_path = HERE / "SCHEDULE.json"
    if freeze_path.exists() or schedule_path.exists():
        raise FileExistsError("refusing to overwrite frozen diagnostic")
    source_sha = sha(SOURCE.read_bytes())
    commit = subprocess.check_output(
        ["git", "log", "-1", "--format=%H", "--", str(SOURCE)],
        cwd=ROOT, text=True).strip()
    assert sha(subprocess.check_output(
        ["git", "show", f"{commit}:examples/{SOURCE.name}"], cwd=ROOT)) == source_sha
    old = json.loads(LATE_MANIFEST.read_bytes())
    old_id = old.pop("candidate_id")
    assert old_id.endswith("h" + sha(canonical(old))[:12])
    old_sha = old["implementation"]["example_source_sha256"]
    base = replace(old, old_sha, source_sha)
    base["implementation"]["frozen_source_commit"] = commit
    policies = {
        "none": "none",
        "late": dict(old["point_decomposition"]["target_lookup_prefetch"]),
        "early": dict(old["point_decomposition"]["target_lookup_prefetch"]),
    }
    policies["early"]["prefetch_policy"] = (
        "prefetch_each_partner_as_canonicalized_then_lookup_in_original_order")
    candidate_ids = {}
    candidate_paths = {}
    identities = HERE / "identities/candidates"
    identities.mkdir(parents=True, exist_ok=True)
    for mode in ("none", "late", "early"):
        candidate = json.loads(canonical(base))
        candidate["point_decomposition"]["target_lookup_prefetch"] = policies[mode]
        candidate["implementation"]["algorithm_flags"]["target_lookup_prefetch"] = policies[mode]
        candidate["point_decomposition"]["implementation"] = (
            "parallel_sharded_s3_index_prelifted_group_sum_roots_batch8_"
            "inverse_target_conjugates_exact_u64_mixed_state_index_slots_"
            f"compact_poly_states_lookup_ab_{mode}")
        candidate["implementation"]["algorithm_flags"]["secondary_diagnostic_modes"] = [
            "none", "late", "early"]
        candidate_id = ("IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0h"
                        + sha(canonical(candidate))[:12])
        candidate["candidate_id"] = candidate_id
        manifest = identities / (candidate_id + ".json")
        manifest.write_bytes(canonical(candidate) + b"\n")
        candidate_ids[mode] = candidate_id
        candidate_paths[mode] = str(manifest.relative_to(ROOT))
    parent = json.loads(PARENT.read_bytes())
    target = next(item for item in parent["targets"]
                  if item["n"] == 53 and item["label"] == "T025")
    point = PARENT.parent / target["point_file_relative"]
    assert sha(point.read_bytes()) == target["point_sha256"]
    assert json.loads(point.read_bytes()) == target["public_point"]
    workload = json.loads((PARENT.parent / "workloads" /
                           ("W" + target["workload_id"] + ".json")).read_bytes())
    assert workload["workload_id"] == sha(canonical(workload["record"]))[:12]
    assert workload["workload_id"] == target["workload_id"]
    reference = json.loads(gzip.decompress(REFERENCE.read_bytes()))
    ref_receipt = json.loads(REFERENCE_RECEIPT.read_bytes())
    assert ref_receipt["verified"] and ref_receipt["memory_within_envelope"]
    assert reference["recovered_scalar"] == target["fixture_scalar"]
    assert reference["target"] == target["public_point"]
    assert reference["prepared_rank"] == 2800
    # Six permutations balance every mode over all three positions. Repeat
    # the complete six-permutation block three times: 18 queries per mode.
    schedule = [mode for _ in range(3)
                for permutation in itertools.permutations((0, 1, 2))
                for mode in permutation]
    assert len(schedule) == 54
    schedule_path.write_bytes(canonical(schedule) + b"\n")
    plan = {
        "kind": "s3_c2800_lookup_ab_t025_diagnostic_freeze_v1",
        "source_sha256": source_sha,
        "source_commit": commit,
        "binary_sha256": sha(binary.read_bytes()),
        "runner_sha256": sha(RUNNER.read_bytes()),
        "candidate_ids": candidate_ids,
        "candidate_manifests": candidate_paths,
        "workload_id": target["workload_id"],
        "parent_freeze_sha256": sha(PARENT.read_bytes()),
        "point_file_relative": target["point_file_relative"],
        "point_sha256": target["point_sha256"],
        "public_point": target["public_point"],
        "fixture_scalar": target["fixture_scalar"],
        "reference_receipt_sha256": sha(REFERENCE_RECEIPT.read_bytes()),
        "reference_gzip_sha256": sha(REFERENCE.read_bytes()),
        "reference_raw_sha256": sha(gzip.decompress(REFERENCE.read_bytes())),
        "schedule_sha256": sha(schedule_path.read_bytes()),
        "schedule_length": 54,
        "queries_per_mode": 18,
        "cache_flush_bytes_per_query": 64 * 1024 * 1024,
        "secondary_metric": "same_preparation_repeated_target_extract_wall_and_thread_cpu",
        "primary_target_count": 1,
        "relation_seed": 20260928,
        "relation_workers": 14,
        "target_workers": 1,
        "timeout_seconds": 1200,
        "memory_envelope_bytes": 34359738368,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    freeze_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"freeze_sha256": sha(freeze_path.read_bytes()),
                      "candidate_ids": candidate_ids,
                      "source_sha256": source_sha,
                      "binary_sha256": plan["binary_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.binary.resolve())
