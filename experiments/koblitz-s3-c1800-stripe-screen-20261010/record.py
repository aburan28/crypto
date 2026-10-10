"""Bind verified striped-schedule call diagnostics to candidate and input IDs."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE_NAME = "s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan"
STRIPE_NAME = "s3_partial_prepared_fullrank_prefetch_projective_c1800stripe_targetspan"
BASE_SCREEN = ROOT / "experiments/koblitz-s3-batch32-screen-20261010"
PARENT = ROOT / "experiments/koblitz-s3-c1800-matchedsource-n53-20261010"
LABELS = ("T205", "T209", "T211", "T212")
SCHEDULE = "c1200_prefix256_then_striped8_c1200_then_complete_c1800"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(record: dict) -> bytes:
    return json.dumps(record, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def replace_source(value, old: str, new: str):
    if isinstance(value, dict):
        return {key: replace_source(item, old, new) for key, item in value.items()}
    if isinstance(value, list):
        return [replace_source(item, old, new) for item in value]
    return new if value == old else value


def main(stripe_binary: Path, baseline_dir: Path, stripe_dir: Path) -> None:
    if (HERE / "CALL_DIAGNOSTIC.json").exists():
        raise FileExistsError(HERE / "CALL_DIAGNOSTIC.json")
    parent_freeze_path = PARENT / "FREEZE.json"
    parent_freeze = json.loads(parent_freeze_path.read_bytes())
    baseline_freeze_path = BASE_SCREEN / "FREEZE.json"
    baseline_freeze = json.loads(baseline_freeze_path.read_bytes())
    baseline_arm = baseline_freeze["curves"]["53"]["arms"]["batch8"]
    reference_arm = parent_freeze["arms"]["expanded"]
    assert baseline_arm["candidate_id"] == reference_arm["candidate_id"]
    assert baseline_arm["binary_name"] == BASE_NAME
    assert baseline_freeze["binary_sha256"][BASE_NAME] == reference_arm["binary_sha256"]
    source_path = ROOT / "examples" / f"{STRIPE_NAME}.rs"
    source_raw = source_path.read_bytes()
    source_sha = sha(source_raw)
    binary_sha = sha(stripe_binary.read_bytes())
    commit = subprocess.check_output(
        ["git", "-C", str(ROOT), "log", "-1", "--format=%H", "--",
         f"examples/{STRIPE_NAME}.rs"], text=True).strip()
    assert subprocess.check_output(
        ["git", "-C", str(ROOT), "show", f"{commit}:examples/{STRIPE_NAME}.rs"]) == source_raw
    tree_oid = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip()
    assert tree_oid == "0a574800502e67882d8f8f88a846c1e80ee795fd"
    reference_path = PARENT / reference_arm["candidate_manifest_relative"]
    reference = json.loads(reference_path.read_bytes())
    assert sha(reference_path.read_bytes()) == reference_arm["candidate_manifest_sha256"]
    candidate = replace_source(reference, reference_arm["source_sha256"], source_sha)
    del candidate["candidate_id"]
    candidate["implementation"]["frozen_source_commit"] = commit
    flags = candidate["implementation"]["algorithm_flags"]
    flags["target_search_schedule"] = SCHEDULE
    flags["target_search_stripes"] = 8
    flags["target_search_lead_states"] = 256
    pdp = candidate["point_decomposition"]
    pdp["target_search_schedule"] = SCHEDULE
    pdp["target_search_stripes"] = 8
    pdp["target_search_lead_states"] = 256
    pdp["target_query_policy"] = (
        "point_seeded_c1200_prefix256_then_striped8_c1200_then_complete_c1800_"
        "with_first_verified_relation_and_worker_drain")
    candidate_id = reference_arm["candidate_id"].rsplit("h", 1)[0] + "h" + sha(canonical(candidate))[:12]
    candidate["candidate_id"] = candidate_id
    identities = HERE / "identities/candidates"
    identities.mkdir(parents=True, exist_ok=True)
    manifest_path = identities / f"{candidate_id}.json"
    manifest_path.write_bytes(canonical(candidate) + b"\n")
    ledger_path = baseline_dir / "raw.jsonl"
    ledger_raw = ledger_path.read_bytes()
    ledger = [json.loads(line) for line in ledger_raw.splitlines()]
    baseline_rows = {(row["n"], row["label"], row["arm"]): row for row in ledger}
    targets = {target["label"]: target for target in parent_freeze["selected_targets"]}
    rows = []
    for label in LABELS:
        target = targets[label]
        base_row = baseline_rows[(53, label, "batch8")]
        baseline_raw = (baseline_dir / f"N53_{label}_batch8_R1.json").read_bytes()
        stripe_raw = (stripe_dir / f"s3-c1800-stripe-smoke-{label}.json").read_bytes()
        stdout_raw = (stripe_dir / f"s3-c1800-stripe-smoke-{label}.stdout").read_bytes()
        stderr_raw = (stripe_dir / f"s3-c1800-stripe-smoke-{label}.stderr").read_bytes()
        baseline = json.loads(baseline_raw)
        stripe = json.loads(stripe_raw)
        assert base_row["verified_native_and_fixture"] and base_row["report_sha256"] == sha(baseline_raw)
        assert base_row["candidate_id"] == reference_arm["candidate_id"]
        assert base_row["workload_id"] == target["workload_id"]
        assert baseline["group_verified"] and stripe["group_verified"]
        assert baseline["target"] == stripe["target"] == target["public_point"]
        assert baseline["recovered_scalar"] == stripe["recovered_scalar"] == target["fixture_scalar"]
        for key in ("n", "a", "factor_base_points", "factor_base_digest", "prepared_rank",
                    "root_table_slots", "target_query_workers", "relation_collection_workers"):
            assert baseline[key] == stripe[key], (label, key)
        assert stripe["target_search_schedule"] == SCHEDULE
        assert stripe["target_search_stripes"] == 8
        assert stripe["target_query_batch_width"] == 8
        assert stripe_raw == stdout_raw
        assert stripe["target_count"] == 1 and len(stripe["target_relation_indices"]) == 4
        assert baseline["target_s3_calls"] > 0 and stripe["target_s3_calls"] > 0
        rows.append({
            "label": label, "workload_id": target["workload_id"],
            "point_sha256": target["point_sha256"],
            "reference_candidate_id": reference_arm["candidate_id"],
            "reference_run_id": base_row["run_id"],
            "reference_report_sha256": sha(baseline_raw),
            "candidate_id": candidate_id,
            "run_id": f"{candidate_id}W{target['workload_id']}R1",
            "report_sha256": sha(stripe_raw),
            "stdout_sha256": sha(stdout_raw), "stderr_sha256": sha(stderr_raw),
            "reference_s3_calls": baseline["target_s3_calls"],
            "stripe_s3_calls": stripe["target_s3_calls"],
            "reference_over_stripe_calls": baseline["target_s3_calls"] / stripe["target_s3_calls"],
            "reference_witness": baseline["target_relation_indices"],
            "stripe_witness": stripe["target_relation_indices"],
            "recovered_scalar": stripe["recovered_scalar"],
            "native_group_verified": True,
        })
    ratio = math.exp(statistics.mean(math.log(row["reference_over_stripe_calls"]) for row in rows))
    receipt = {
        "kind": "s3_c1800_stripe_call_diagnostic_v1",
        "status": "complete_native_verified",
        "candidate_id": candidate_id,
        "candidate_manifest_relative": str(manifest_path.relative_to(HERE)),
        "candidate_manifest_sha256": sha(manifest_path.read_bytes()),
        "reference_candidate_id": reference_arm["candidate_id"],
        "parent_freeze_sha256": sha(parent_freeze_path.read_bytes()),
        "baseline_screen_freeze_sha256": sha(baseline_freeze_path.read_bytes()),
        "baseline_screen_ledger_sha256": sha(ledger_raw),
        "source_sha256": source_sha, "binary_sha256": binary_sha,
        "source_commit": commit, "source_tree_git_oid": tree_oid,
        "cargo_lock_sha256": sha((ROOT / "Cargo.lock").read_bytes()),
        "recorder_sha256": sha(Path(__file__).read_bytes()),
        "target_labels": list(LABELS), "rows": rows,
        "geometric_reference_over_stripe_s3_calls": ratio,
        "online_wall_ratio": None, "controlled_speedup": None,
        "timing_use": "correctness and deterministic S3 call counts only; smoke timings excluded",
    }
    (HERE / "CALL_DIAGNOSTIC.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"candidate_id": candidate_id, "verified": len(rows),
                      "call_ratio": ratio}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stripe-binary", type=Path, required=True)
    parser.add_argument("--baseline-dir", type=Path, required=True)
    parser.add_argument("--stripe-dir", type=Path, required=True)
    args = parser.parse_args()
    main(args.stripe_binary.resolve(), args.baseline_dir.resolve(), args.stripe_dir.resolve())
