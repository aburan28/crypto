"""Assign an IC1 identity to four-worker single-target C1800 extraction."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OLD = ROOT / "experiments/koblitz-s3-c1800-norow-20261010/identities/candidates/IC1N53Ckb1fb190800PDP4rootRCsampleLAgaussTDdirectISO0h671db6e1f6fd.json"
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan_norow.rs"
BINARY = ROOT / "target/release/examples/s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan_norow"
CONTROL = HERE / "local/control/workers4_receipt.json"
CONTROL_RAW = HERE / "local/control/workers4_report.json.gz"
OUTPUT = HERE / "identities/candidates"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def main() -> None:
    manifest = json.loads(OLD.read_bytes())
    old_id = manifest.pop("candidate_id")
    assert old_id.endswith("h" + sha(canonical(manifest))[:12])
    control = json.loads(CONTROL.read_bytes())
    raw = gzip.decompress(CONTROL_RAW.read_bytes())
    report = json.loads(raw)
    assert control["scalar_verified"] and control["full_rank_before_target"]
    assert control["memory_within_envelope"] and control["target_workers"] == 4
    assert control["raw_report_sha256"] == sha(raw)
    assert control["source_sha256"] == sha(SOURCE.read_bytes())
    assert control["binary_sha256"] == sha(BINARY.read_bytes())
    assert report["target_query_workers"] == 4
    assert report["factor_base_digest"] == manifest["factor_base"]["enumerated_set_digest"]
    assert report["factor_base_points"] == manifest["factor_base"]["actual_usable_points_B"]
    manifest["implementation"]["algorithm_flags"]["target_query_workers"] = 4
    manifest["implementation"]["algorithm_flags"]["target_query_worker_spacing"] = (
        "index_states_len_divided_by_workers")
    manifest["implementation"]["algorithm_flags"]["target_query_stop_policy"] = (
        "first_verified_relation_cooperative_cancellation_and_worker_drain")
    manifest["point_decomposition"]["target_query_workers"] = 4
    manifest["point_decomposition"]["target_query_worker_spacing"] = (
        "index_states_len_divided_by_workers")
    manifest["point_decomposition"]["target_query_stop_policy"] = (
        "first_verified_relation_cooperative_cancellation_and_worker_drain")
    manifest["target_descent"]["target_query_workers"] = 4
    candidate_id = "IC1N53Ckb1fb190800PDP4rootRCsampleLAgaussTDdirectISO0h" + sha(
        canonical(manifest))[:12]
    OUTPUT.mkdir(parents=True, exist_ok=True)
    path = OUTPUT / (candidate_id + ".json")
    if path.exists():
        raise FileExistsError(path)
    path.write_bytes(canonical({**manifest, "candidate_id": candidate_id}) + b"\n")
    print(json.dumps({"candidate_id": candidate_id,
                      "manifest_sha256": sha(path.read_bytes())}, sort_keys=True))


if __name__ == "__main__":
    main()
