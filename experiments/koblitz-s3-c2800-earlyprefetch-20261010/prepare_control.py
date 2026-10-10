"""Freeze the exact early-prefetch candidate and one-point T025 control."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_c2800earlyprefetch_targetspan.rs"
RUNNER = HERE / "run_control.py"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
OLD = (ROOT / "experiments/koblitz-s3-c2800-prefetch-20261010/"
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


def replace(value: object, substitutions: dict) -> object:
    if isinstance(value, dict):
        return {key: replace(item, substitutions) for key, item in value.items()}
    if isinstance(value, list):
        return [replace(item, substitutions) for item in value]
    return substitutions.get(value, value)


def prepare(binary: Path) -> None:
    freeze_path = HERE / "CONTROL_FREEZE.json"
    if freeze_path.exists():
        raise FileExistsError(freeze_path)
    old = json.loads(OLD.read_bytes())
    old_id = old.pop("candidate_id")
    assert old_id.endswith("h" + sha(canonical(old))[:12])
    old_sha = old["implementation"]["example_source_sha256"]
    new_sha = sha(SOURCE.read_bytes())
    commit = subprocess.check_output(
        ["git", "log", "-1", "--format=%H", "--", str(SOURCE)],
        cwd=ROOT, text=True).strip()
    assert sha(subprocess.check_output(
        ["git", "show", f"{commit}:examples/{SOURCE.name}"], cwd=ROOT)) == new_sha
    old_policy = "compute_all_batch_roots_then_prefetch_then_lookup"
    new_policy = "prefetch_each_partner_as_canonicalized_then_lookup_in_original_order"
    candidate = replace(old, {old_sha: new_sha,
                              old["implementation"]["frozen_source_commit"]: commit,
                              old_policy: new_policy})
    candidate["point_decomposition"]["implementation"] = (
        candidate["point_decomposition"]["implementation"]
        .replace("_batch8_root_prefetch", "_batch8_early_root_prefetch"))
    assert candidate["point_decomposition"]["target_lookup_prefetch"]["prefetch_policy"] == new_policy
    assert candidate["implementation"]["algorithm_flags"]["target_lookup_prefetch"]["prefetch_policy"] == new_policy
    candidate_id = ("IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0h"
                    + sha(canonical(candidate))[:12])
    candidate["candidate_id"] = candidate_id
    identities = HERE / "identities/candidates"
    identities.mkdir(parents=True, exist_ok=True)
    manifest = identities / (candidate_id + ".json")
    manifest.write_bytes(canonical(candidate) + b"\n")
    parent = json.loads(PARENT.read_bytes())
    target = next(item for item in parent["targets"]
                  if item["n"] == 53 and item["label"] == "T025")
    point = PARENT.parent / target["point_file_relative"]
    assert sha(point.read_bytes()) == target["point_sha256"]
    assert json.loads(point.read_bytes()) == target["public_point"]
    ref = json.loads(gzip.decompress(REFERENCE.read_bytes()))
    reference_receipt = json.loads(REFERENCE_RECEIPT.read_bytes())
    assert reference_receipt["verified"] and reference_receipt["memory_within_envelope"]
    assert ref["recovered_scalar"] == target["fixture_scalar"]
    assert ref["group_verified"] is True and ref["target"] == target["public_point"]
    assert ref["factor_base_points"] == 296800 and ref["prepared_rank"] == 2800
    plan = {
        "kind": "s3_c2800_earlyprefetch_t025_control_freeze_v1",
        "candidate_id": candidate_id,
        "candidate_manifest": str(manifest.relative_to(ROOT)),
        "source_sha256": new_sha,
        "source_commit": commit,
        "binary_sha256": sha(binary.read_bytes()),
        "runner_sha256": sha(RUNNER.read_bytes()),
        "workload_id": target["workload_id"],
        "parent_freeze_sha256": sha(PARENT.read_bytes()),
        "point_file_relative": target["point_file_relative"],
        "point_sha256": target["point_sha256"],
        "public_point": target["public_point"],
        "fixture_scalar": target["fixture_scalar"],
        "reference_receipt_sha256": sha(REFERENCE_RECEIPT.read_bytes()),
        "reference_gzip_sha256": sha(REFERENCE.read_bytes()),
        "reference_raw_sha256": sha(gzip.decompress(REFERENCE.read_bytes())),
        "relation_seed": 20260928,
        "relation_workers": 14,
        "timeout_seconds": 1200,
        "memory_envelope_bytes": 34359738368,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
        "stage_gate": "verified_scalar_same_first_witness_full_rank_and_rss_at_most_32gib",
    }
    freeze_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"candidate_id": candidate_id,
                      "freeze_sha256": sha(freeze_path.read_bytes()),
                      "source_sha256": new_sha,
                      "binary_sha256": plan["binary_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.binary.resolve())
