"""Freeze the T025 control and exact C2800 prefetch candidate identity."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_c2800prefetch_targetspan.rs"
RUNNER = HERE / "run_control.py"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
OLD = (ROOT / "experiments/koblitz-s3-c2800-mixed-20261009/identities/"
       "candidates/IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0hbbd5b5a111fb.json")
REFERENCE = (ROOT / "experiments/koblitz-s3-c2800-mixed-20261009/local/"
             "feasibility_report.json.gz")


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
    freeze_path = HERE / "CONTROL_FREEZE.json"
    if freeze_path.exists():
        raise FileExistsError(freeze_path)
    source_sha = sha(SOURCE.read_bytes())
    binary_sha = sha(binary.read_bytes())
    old = json.loads(OLD.read_bytes())
    old_id = old.pop("candidate_id")
    assert old_id.endswith("h" + sha(canonical(old))[:12])
    old_sha = old["implementation"]["example_source_sha256"]
    candidate = replace(old, old_sha, source_sha)
    commit = subprocess.check_output(
        ["git", "log", "-1", "--format=%H", "--", str(SOURCE)],
        cwd=ROOT, text=True).strip()
    assert commit and sha(subprocess.check_output(
        ["git", "show", f"{commit}:examples/{SOURCE.name}"],
        cwd=ROOT)) == source_sha
    candidate["implementation"]["frozen_source_commit"] = commit
    prefetch = {
        "batch_width": 8,
        "max_candidates": 32,
        "prefetch_hint": "read_first_slot_L1",
        "prefetch_policy": "compute_all_batch_roots_then_prefetch_then_lookup",
        "lookup_order": "original_state_and_root_order",
        "architecture": "x86_64_mm_prefetch_or_aarch64_prfm_with_portable_noop",
    }
    candidate["point_decomposition"]["target_lookup_prefetch"] = prefetch
    candidate["point_decomposition"]["implementation"] += "_batch8_root_prefetch"
    candidate["implementation"]["algorithm_flags"]["target_lookup_prefetch"] = prefetch
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
    workload = PARENT.parent / "workloads" / ("W" + target["workload_id"] + ".json")
    work = json.loads(workload.read_bytes())
    ident = work["workload_id"]
    assert ident == target["workload_id"] == sha(canonical(work["record"]))[:12]
    assert work["identity_sha256"] == sha(canonical(work["record"]))
    reference_raw = __import__("gzip").decompress(REFERENCE.read_bytes())
    reference = json.loads(reference_raw)
    assert reference["target"] == target["public_point"]
    assert reference["recovered_scalar"] == target["fixture_scalar"]
    assert reference["factor_base_points"] == 296800
    plan = {
        "kind": "s3_c2800_prefetch_t025_control_freeze_v1",
        "source_sha256": source_sha,
        "source_commit": commit,
        "binary_sha256": binary_sha,
        "runner_sha256": sha(RUNNER.read_bytes()),
        "candidate_id": candidate_id,
        "candidate_manifest": str(manifest.relative_to(ROOT)),
        "workload_id": ident,
        "parent_freeze_sha256": sha(PARENT.read_bytes()),
        "point_file_relative": target["point_file_relative"],
        "point_sha256": target["point_sha256"],
        "public_point": target["public_point"],
        "fixture_scalar": target["fixture_scalar"],
        "reference_gzip_sha256": sha(REFERENCE.read_bytes()),
        "reference_raw_sha256": sha(reference_raw),
        "relation_seed": 20260928,
        "relation_workers": 14,
        "memory_envelope_bytes": 34359738368,
        "timeout_seconds": 1200,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
        "stage_gate": "verified_t025_same_first_witness_full_rank_and_rss_at_most_32gib",
    }
    freeze_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"candidate_id": candidate_id, "source_sha256": source_sha,
                      "binary_sha256": binary_sha,
                      "freeze_sha256": sha(freeze_path.read_bytes())},
                     sort_keys=True))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.binary.resolve())
