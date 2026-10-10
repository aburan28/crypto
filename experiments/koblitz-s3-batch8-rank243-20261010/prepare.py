"""Freeze a one-worker rank-243 batch-inversion comparison after smoke checks."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009"
PRIOR = ROOT / "experiments/koblitz-s3-prepared-20261009"
TEMPLATE_COMMIT = "e7ae9012d02a13ab48851c00437e2f432c7775f1"
OLD_TEMPLATE = "experiments/koblitz-s3-batch8-20261009/identities/candidates"
NAMES = {"baseline": "s3_partial_prepared_parallel_targetspan",
         "batch8": "s3_partial_prepared_batch8_targetspan"}
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(record: dict) -> bytes:
    return json.dumps(record, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def read_template(arm: str) -> dict:
    suffix = ("h1e2914fbd8c2" if arm == "baseline" else "hff55c3b35e32")
    name = f"IC1N53Ckb1fb41976PDP4rootRCsampleLAgaussTDdirectISO0{suffix}.json"
    raw = subprocess.check_output(["git", "-C", str(ROOT), "show",
                                   f"{TEMPLATE_COMMIT}:{OLD_TEMPLATE}/{name}"])
    return json.loads(raw)


def make_candidate(n: int, arm: str, prior: dict, smoke: dict,
                   source_hash: str, source_tree: str, head: str) -> dict:
    record = read_template(arm)
    record.pop("candidate_id")
    for key in ("curve", "field", "endomorphism", "factor_base"):
        record[key] = prior[key]
    point_decomposition = record["point_decomposition"]
    point_decomposition["root_table_slots"] = smoke["root_table_slots"]
    point_decomposition["target_query_workers"] = 1
    record["implementation"]["algorithm_flags"].update(
        columns=244, prepare_rank=243, target_query_workers=1)
    record["implementation"]["frozen_source_commit"] = head
    record["implementation"]["rust_library_source_tree_git_oid"] = source_tree
    record["implementation"]["example_source_sha256"] = source_hash
    record["factor_base"]["source_sha256"] = source_hash
    point_decomposition["implementation_source_sha256"] = source_hash
    record["relation_collection"]["pre_target_rank_threshold"] = 243
    record["relation_collection"]["source_sha256"] = source_hash
    record["relation_linear_algebra"]["rank_criterion"] = (
        "fixed_pre_target_rank_then_target_row_span_or_full_244_fallback")
    record["relation_linear_algebra"]["source_sha256"] = source_hash
    record["target_descent"]["target_query_workers"] = 1
    record["target_descent"]["source_sha256"] = source_hash
    if arm == "batch8":
        assert point_decomposition["target_query_batch_width"] == 8
        assert record["implementation"]["algorithm_flags"]["target_query_batch_width"] == 8
    prefix = f"IC1N{n}Ckb1fb{smoke['factor_base_points']}PDP4rootRCsampleLAgaussTDdirectISO0h"
    record["candidate_id"] = prefix + sha(canonical(record))[:12]
    return record


def main(binaries: Path, smoke_dir: Path) -> None:
    if (HERE / "FREEZE.json").exists():
        raise FileExistsError(HERE / "FREEZE.json")
    prior_freeze_path = PARENT / "FREEZE.json"
    prior_freeze = json.loads(prior_freeze_path.read_bytes())
    head = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
    source_tree = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip()
    source_hashes = {arm: sha((ROOT / "examples" / f"{name}.rs").read_bytes())
                     for arm, name in NAMES.items()}
    binary_hashes = {arm: sha((binaries / name).read_bytes()) for arm, name in NAMES.items()}
    manifests = {}
    smoke_hashes = {}
    for n in (41, 53):
        prior_name = prior_freeze["curves"][str(n)]["reference_candidate_id"]
        prior_path = PRIOR / "identities/candidates" / f"{prior_name}.json"
        prior = json.loads(prior_path.read_bytes())
        assert prior["candidate_id"] == prior_name
        assert prior["factor_base"]["actual_usable_points_B"] == prior_freeze["curves"][str(n)]["factor_base_points"]
        reports = {}
        for arm in NAMES:
            smoke_path = smoke_dir / f"n{n}_{arm}.json"
            smoke_hashes[f"n{n}_{arm}"] = sha(smoke_path.read_bytes())
            report = json.loads(smoke_path.read_bytes())
            assert report["n"] == n and report["group_verified"] is True
            assert report["factor_base_points"] == prior["factor_base"]["actual_usable_points_B"]
            assert report["factor_base_digest"] == prior["factor_base"]["enumerated_set_digest"]
            assert report["prepared_rank"] == 243 and report["target_query_workers"] == 1
            assert report["target_count"] == 1
            assert report["target_query_batch_width"] == 8 if arm == "batch8" else True
            assert abs(sum(report["timing_ms"][k] for k in PHASES) -
                       report["timing_ms"]["target_online_after_reusable_setup"]) < .02
            reports[arm] = report
        assert reports["baseline"]["target"] == reports["batch8"]["target"]
        assert reports["baseline"]["recovered_scalar"] == reports["batch8"]["recovered_scalar"]
        assert reports["baseline"]["root_table_slots"] == reports["batch8"]["root_table_slots"]
        for arm in NAMES:
            manifest = make_candidate(n, arm, prior, reports[arm],
                                      source_hashes[arm], source_tree, head)
            path = HERE / "identities/candidates" / f"{manifest['candidate_id']}.json"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(canonical(manifest) + b"\n")
            manifests[f"n{n}_{arm}"] = {"candidate_id": manifest["candidate_id"],
                                        "path": str(path.relative_to(HERE)),
                                        "sha256": sha(path.read_bytes())}
    selected = [t for t in prior_freeze["targets"]
                if t["label"] in [f"T{i:03}" for i in range(25, 29)]]
    assert len(selected) == 8
    freeze = {
        "kind": "s3_rank243_one_worker_batch8_screen_v1",
        "source_commit": head, "source_tree_git_oid": source_tree,
        "cargo_lock_sha256": sha((ROOT / "Cargo.lock").read_bytes()),
        "parent_freeze_sha256": sha(prior_freeze_path.read_bytes()),
        "prior_candidate_manifest_sha256": {
            str(n): sha((PRIOR / "identities/candidates" /
                         f"{prior_freeze['curves'][str(n)]['reference_candidate_id']}.json").read_bytes())
            for n in (41, 53)},
        "source_sha256": source_hashes, "binary_sha256": binary_hashes,
        "binary_names": NAMES, "candidate_manifests": manifests,
        "smoke_report_sha256": smoke_hashes,
        "selected_targets": selected, "selected_indices": list(range(25, 29)),
        "pre_target_rank": 243, "orbit_columns": 244,
        "relation_seed": 20260928, "relation_workers": 14,
        "target_workers": 1, "repetitions_per_arm_per_target": 2,
        "planned_run_count": 32, "per_run_timeout_seconds": 180,
        "target_point_input": "one frozen public point per fresh process",
        "paired_order": "ABBA_or_BAAB_alternating_by_curve_and_target",
        "online_phases": list(PHASES),
        "development_advance_gate": "verified_online_geometric_baseline_over_batch8_at_least_2_on_both_curves",
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
        "runner_sha256": sha((HERE / "run.py").read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze.py").read_bytes()),
    }
    (HERE / "FREEZE.json").write_text(json.dumps(freeze, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"freeze_sha256": sha((HERE / "FREEZE.json").read_bytes()),
                      "candidate_ids": {key: value["candidate_id"] for key, value in manifests.items()},
                      "planned_run_count": 32}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binaries", type=Path, required=True)
    parser.add_argument("--smoke-dir", type=Path, required=True)
    args = parser.parse_args()
    main(args.binaries.resolve(), args.smoke_dir.resolve())
