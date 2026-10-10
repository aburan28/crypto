"""Freeze the exact same-base N53 C2800 mixed/prefetch paired screen."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
CONTROL = HERE / "local/verified_t025/receipt.json"
RUNNER = HERE / "run_paired.py"
BASELINE_MANIFEST = (ROOT / "experiments/koblitz-s3-c2800-mixed-20261009/"
                     "identities/candidates/"
                     "IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0hbbd5b5a111fb.json")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def arm(binary: Path, manifest_path: Path) -> dict:
    manifest = json.loads(manifest_path.read_bytes())
    ident = manifest.pop("candidate_id")
    assert ident.endswith("h" + sha(canonical(manifest))[:12])
    source_hash = manifest["implementation"]["example_source_sha256"]
    source = (ROOT / "examples/s3_partial_prepared_fullrank_c2800mixed_targetspan.rs"
              if manifest_path == BASELINE_MANIFEST else
              ROOT / "examples/s3_partial_prepared_fullrank_c2800prefetch_targetspan.rs")
    assert sha(source.read_bytes()) == source_hash
    assert manifest["factor_base"]["actual_usable_points_B"] == 296800
    assert manifest["factor_base"]["effective_relation_columns"] == 2800
    assert manifest["implementation"]["algorithm_flags"]["root_table_slots"] == 1073741824
    return {"candidate_id": ident,
            "candidate_manifest": str(manifest_path.relative_to(ROOT)),
            "source_relative": str(source.relative_to(ROOT)),
            "source_sha256": source_hash,
            "binary_sha256": sha(binary.read_bytes()),
            "factor_base_digest": manifest["factor_base"]["enumerated_set_digest"]}


def prepare(baseline: Path, prefetch: Path) -> None:
    freeze_path = HERE / "SCREEN_FREEZE.json"
    if freeze_path.exists():
        raise FileExistsError(freeze_path)
    control = json.loads(CONTROL.read_bytes())
    assert control["verified"] and control["memory_within_envelope"]
    parent = json.loads(PARENT.read_bytes())
    labels = [f"T{i:03}" for i in range(173, 181)]
    n53 = {t["label"]: t for t in parent["targets"] if t["n"] == 53}
    assert all(label in n53 for label in labels)
    new_manifest = next((HERE / "identities/candidates").glob(
        "IC1N53Ckb1fb296800PDP4rootRCsampleLAgaussTDdirectISO0h*.json"))
    arms = {"baseline": arm(baseline, BASELINE_MANIFEST),
            "prefetch": arm(prefetch, new_manifest)}
    assert arms["prefetch"]["candidate_id"] == control["candidate_id"]
    assert arms["prefetch"]["binary_sha256"] == control["binary_sha256"]
    assert arms["baseline"]["factor_base_digest"] == arms["prefetch"]["factor_base_digest"]
    plan = {
        "kind": "s3_n53_c2800_prefetch_paired_screen_freeze_v1",
        "parent_freeze_sha256": sha(PARENT.read_bytes()),
        "control_receipt_sha256": sha(CONTROL.read_bytes()),
        "runner_sha256": sha(RUNNER.read_bytes()),
        "arms": arms,
        "target_labels": labels,
        "repetitions_per_arm_per_target": 1,
        "pair_order": "baseline_prefetch_for_even_target_index_else_prefetch_baseline",
        "relation_seed": 20260928,
        "relation_workers": 14,
        "target_workers": 1,
        "timeout_seconds_per_run": 1200,
        "memory_envelope_bytes": 34359738368,
        "stop_after_first_invalid_or_over_budget_run": True,
        "exploratory_geomean_speedup_threshold": "2",
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    freeze_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"freeze_sha256": sha(freeze_path.read_bytes()),
                      "candidates": {key: value["candidate_id"]
                                     for key, value in arms.items()},
                      "targets": len(labels)}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--prefetch", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.baseline.resolve(), args.prefetch.resolve())
