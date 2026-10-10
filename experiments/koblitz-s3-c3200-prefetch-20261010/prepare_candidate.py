"""Assign an IC1 identity to the verified C3200 quotient-root prefetch path."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OLD = ROOT / "experiments/koblitz-s3-c3200-projective-20261010/identities/candidates/IC1N53Ckb1fb339200PDP4rootRCsampleLAgaussTDdirectISO0h0e595bd6143a.json"
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_c3200projective_prefetch_targetspan.rs"
BINARY = ROOT / "target/release/examples/s3_partial_prepared_fullrank_c3200projective_prefetch_targetspan"
CONTROL = HERE / "local/control/receipt.json"
REPORT = HERE / "local/control/report.json.gz"
SAGE_REPLAY = HERE / "independent_sage_replay.json"
OUTPUT = HERE / "identities/candidates"
PREFETCH_POLICY = (
    "prefetch_all_mixed_quotient_home_slots_after_eight_state_root_batch_then_original_lookup_order")
SOURCE_COMMIT = "d860557cea446fe541146aa8363eb68ac2f1381a"


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
    raw = gzip.decompress(REPORT.read_bytes())
    report = json.loads(raw)
    replay = json.loads(SAGE_REPLAY.read_bytes())
    source_sha = sha(SOURCE.read_bytes())
    assert control["scalar_verified"] and control["full_rank_before_target"]
    assert control["index_layout_ok"] and control["memory_within_envelope"]
    assert control["same_first_relation_witness_as_reference"]
    assert control["source_sha256"] == source_sha
    assert control["binary_sha256"] == sha(BINARY.read_bytes())
    assert control["raw_report_sha256"] == sha(raw)
    assert report["factor_base_points"] == 339200
    assert report["target_lookup_prefetch_policy"] == PREFETCH_POLICY
    assert replay["generator_order_verified"]
    assert replay["scalar_times_generator_equals_public_target"]
    assert replay["report_raw_sha256"] == sha(raw)
    old_source_sha = manifest["implementation"]["example_source_sha256"]
    assert old_source_sha != source_sha
    for section, field in (("factor_base", "source_sha256"),
                           ("implementation", "example_source_sha256"),
                           ("point_decomposition", "implementation_source_sha256"),
                           ("relation_collection", "source_sha256"),
                           ("relation_linear_algebra", "source_sha256"),
                           ("target_descent", "source_sha256")):
        assert manifest[section][field] == old_source_sha
        manifest[section][field] = source_sha
    manifest["implementation"]["frozen_source_commit"] = SOURCE_COMMIT
    manifest["implementation"]["algorithm_flags"]["target_lookup_prefetch_policy"] = (
        PREFETCH_POLICY)
    pdp = manifest["point_decomposition"]
    pdp["implementation"] += "_target_mixed_quotient_home_slot_prefetch"
    pdp["target_lookup_prefetch_policy"] = PREFETCH_POLICY
    candidate_id = "IC1N53Ckb1fb339200PDP4rootRCsampleLAgaussTDdirectISO0h" + sha(
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
