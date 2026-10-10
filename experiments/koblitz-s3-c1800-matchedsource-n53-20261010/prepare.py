"""Freeze the eight N53 targets for a same-source C1200/C1800 screen."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-c1800-combined-20261010"
SMOKE = Path("/private/tmp/s3-n53-source-matched-c1200-smoke-20261010")
BINARIES = Path("/private/tmp/crypto-s3-release-examples-20261010")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(record: dict) -> bytes:
    return json.dumps(record, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def write_json(path: Path, record: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")


def make_baseline(source_candidate: dict, smoke: dict) -> dict:
    candidate = copy.deepcopy(source_candidate)
    candidate.pop("candidate_id")
    columns = smoke["orbit_columns"]
    assert columns == smoke["prepared_rank"] == 1200
    base = candidate["factor_base"]
    base.update({
        "actual_usable_points_B": smoke["factor_base_points"],
        "effective_relation_columns": columns,
        "enumerated_set_digest": smoke["factor_base_digest"],
        "last_scanned_x": smoke["factor_base_scanned_x"],
        "requested_orbit_columns": columns,
        "start_x": smoke["factor_base_scan_start_x"],
    })
    point = candidate["point_decomposition"]
    point.update({
        "indexed_state_payload_bytes": smoke["index_state_payload_bytes"],
        "root_table_payload_bytes": smoke["root_table_slots"] * 13,
        "root_table_slots": smoke["root_table_slots"],
        "target_query_policy": "point_seeded_native_rotating_state_order_with_first_verified_relation_and_worker_drain",
        "target_search_prefix_columns": columns,
        "target_search_schedule": "native_rotating_state_order",
    })
    relation = candidate["relation_collection"]
    relation["pre_target_rank_threshold"] = columns
    relation["stop_criterion"] = (
        f"complete_{columns}_rank_before_target_input_if_attainable_then_online_span_fallback_if_incomplete")
    candidate["relation_linear_algebra"]["rank_criterion"] = (
        f"full_{columns}_pre_target_if_attained_then_existing_fallback")
    flags = candidate["implementation"]["algorithm_flags"]
    flags.update({
        "columns": columns,
        "index_state_payload_bytes": smoke["index_state_payload_bytes"],
        "prepare_rank": columns,
        "root_table_payload_bytes": smoke["root_table_slots"] * 13,
        "root_table_slots": smoke["root_table_slots"],
        "target_search_prefix_columns": columns,
        "target_search_schedule": "native_rotating_state_order",
    })
    candidate["candidate_id"] = (
        f"IC1N53Ckb1fb{smoke['factor_base_points']}PDP4rootRCsampleLAgaussTDdirectISO0h"
        + sha(canonical(candidate))[:12])
    return candidate


def prepare() -> None:
    freeze_path = HERE / "FREEZE.json"
    assert not freeze_path.exists(), "Freeze is immutable"
    prior_path = PARENT / "CONFIRM_FREEZE.json"
    prior = json.loads(prior_path.read_bytes())
    targets = prior["targets"]
    labels = [f"T{number:03}" for number in range(205, 213)]
    assert [target["label"] for target in targets] == labels
    assert len({target["workload_id"] for target in targets}) == 8
    for target in targets:
        point = ROOT / target["point_file_relative"]
        assert sha(point.read_bytes()) == target["point_sha256"]
        assert json.loads(point.read_bytes()) == target["public_point"]
        workload = ROOT / target["workload_manifest_relative"]
        assert sha(workload.read_bytes()) == target["workload_manifest_sha256"]
        identity = json.loads(workload.read_bytes())
        assert identity["workload_id"] == target["workload_id"]
        assert identity["identity_sha256"] == target["workload_identity_sha256"]
        assert sha(canonical(identity["record"])) == identity["identity_sha256"]
        assert identity["record"]["target_point"] == target["public_point"]
    smoke_raw = SMOKE.with_suffix(".json").read_bytes()
    assert smoke_raw == SMOKE.with_suffix(".stdout").read_bytes()
    smoke = json.loads(smoke_raw)
    assert smoke["n"] == 53 and smoke["a"] == 0
    assert smoke["group_verified"] is True
    assert smoke["target"] == targets[0]["public_point"]
    assert smoke["recovered_scalar"] == targets[0]["fixture_scalar"]
    assert smoke["final_rank_complete"] is True and smoke["fallback_attempts"] == 0
    assert smoke["target_query_workers"] == 1 and smoke["relation_collection_workers"] == 14
    assert smoke["target_search_schedule"] == "native_rotating_state_order"
    assert smoke["factor_base_points"] == prior["arms"]["baseline"]["factor_base_points"]
    assert smoke["factor_base_digest"] == prior["arms"]["baseline"]["factor_base_digest"]
    source_arm = prior["arms"]["combined"]
    source = ROOT / source_arm["source_relative"]
    binary_name = source.stem
    binary = BINARIES / binary_name
    assert sha(source.read_bytes()) == source_arm["source_sha256"]
    binary_sha = sha(binary.read_bytes())
    source_manifest = ROOT / source_arm["candidate_manifest_relative"]
    assert sha(source_manifest.read_bytes()) == source_arm["candidate_manifest_sha256"]
    baseline_candidate = make_baseline(json.loads(source_manifest.read_bytes()), smoke)
    baseline_manifest = HERE / "identities/candidates" / f"{baseline_candidate['candidate_id']}.json"
    write_json(baseline_manifest, baseline_candidate)
    baseline = {
        "binary_name": binary_name,
        "binary_sha256": binary_sha,
        "source_sha256": source_arm["source_sha256"],
        "candidate_id": baseline_candidate["candidate_id"],
        "candidate_manifest_relative": str(baseline_manifest.relative_to(HERE)),
        "candidate_manifest_sha256": sha(baseline_manifest.read_bytes()),
        "columns": 1200, "prepare_rank": 1200,
        "factor_base_points": smoke["factor_base_points"],
        "factor_base_digest": smoke["factor_base_digest"],
        "root_table_slots": smoke["root_table_slots"],
        "target_search_schedule": "native_rotating_state_order",
        "target_search_prefix_columns": 1200,
    }
    expanded = {
        "binary_name": binary_name,
        "binary_sha256": binary_sha,
        "source_sha256": source_arm["source_sha256"],
        "candidate_id": source_arm["candidate_id"],
        "candidate_manifest_relative": "../" + str(source_manifest.relative_to(ROOT / "experiments")),
        "candidate_manifest_sha256": source_arm["candidate_manifest_sha256"],
        "columns": source_arm["columns"],
        "prepare_rank": source_arm["prepare_rank"],
        "factor_base_points": source_arm["factor_base_points"],
        "factor_base_digest": source_arm["factor_base_digest"],
        "root_table_slots": source_arm["root_table_slots"],
        "target_search_schedule": source_arm["target_search_schedule"],
        "target_search_prefix_columns": source_arm["target_search_prefix_columns"],
    }
    run_order = []
    for index, label in enumerate(labels):
        order = (("baseline", 1), ("expanded", 1), ("expanded", 2), ("baseline", 2))
        if index % 2:
            order = (("expanded", 1), ("baseline", 1), ("baseline", 2), ("expanded", 2))
        run_order.extend({"label": label, "arm": arm, "repetition": repetition}
                         for arm, repetition in order)
    freeze = {
        "kind": "s3_n53_c1200_c1800_matchedsource8_v1",
        "source_tree_git_oid": subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip(),
        "cargo_lock_sha256": sha((ROOT / "Cargo.lock").read_bytes()),
        "parent_freeze_sha256": sha(prior_path.read_bytes()),
        "parent_analysis_sha256": sha((PARENT / "analysis.json").read_bytes()),
        "parent_sage_replay_sha256": sha(
            (PARENT / "local/paired_confirm/independent_sage_replay.json").read_bytes()),
        "selected_labels": labels, "selected_targets": targets,
        "run_order": run_order, "planned_run_count": len(run_order),
        "relation_seed": prior["relation_seed"],
        "relation_workers": 14, "target_workers": 1,
        "per_run_timeout_seconds": 180,
        "observed_peak_rss_envelope_bytes": 32 * 1024**3,
        "repetitions_per_arm_per_target": 2,
        "one_public_target_per_process": True,
        "measured_interval": "target_online_after_reusable_setup; all target-dependent phases through scalar replay",
        "comparison_variable": "C1200 versus C1800 factor-base and S3 index with identical source/executable; C1800 first searches its embedded C1200 prefix",
        "development_gate_online_ratio_at_least": "2.0",
        "matched_source_gate_repeat_aware_lower_at_least": "2.0",
        "bootstrap_resamples": 20000,
        "bootstrap_seed": 20261010,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
        "smoke_report_sha256": sha(smoke_raw),
        "smoke_archive_relative": "local/smoke_raw.tar.zst",
        "smoke_archive_sha256": sha((HERE / "local/smoke_raw.tar.zst").read_bytes()),
        "runner_sha256": sha((HERE / "run.py").read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze.py").read_bytes()),
        "replay_sha256": sha((HERE / "replay_sage.py").read_bytes()),
        "preparer_sha256": sha(Path(__file__).read_bytes()),
        "arms": {"baseline": baseline, "expanded": expanded},
    }
    write_json(freeze_path, freeze)
    print(json.dumps({"freeze_sha256": sha(freeze_path.read_bytes()),
                      "candidate_ids": {arm: data["candidate_id"] for arm, data in freeze["arms"].items()},
                      "targets": labels, "planned_runs": len(run_order)}))


if __name__ == "__main__":
    prepare()
