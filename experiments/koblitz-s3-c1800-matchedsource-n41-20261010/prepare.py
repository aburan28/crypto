"""Freeze fresh N41 targets for a same-source C1200/C1800 ablation."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009"
PILOT = ROOT / "experiments/koblitz-s3-c1800-n41-20261010"
CONFIRMATION = ROOT / "experiments/koblitz-s3-c1800-n41-confirmation24-20261010"
SMOKE = Path("/private/tmp/s3-n41-source-matched-c1200-smoke-20261010")
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
    factor_base = candidate["factor_base"]
    factor_base.update({
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
        "target_search_prefix_columns": columns,
        "target_search_schedule": "native_rotating_state_order",
        "target_query_policy": "point_seeded_native_rotating_state_order_with_first_verified_relation_and_worker_drain",
    })
    relation = candidate["relation_collection"]
    relation["pre_target_rank_threshold"] = columns
    relation["stop_criterion"] = (
        f"complete_{columns}_rank_before_target_input_if_attainable_then_online_span_fallback_if_incomplete")
    algebra = candidate["relation_linear_algebra"]
    algebra["rank_criterion"] = f"full_{columns}_pre_target_if_attained_then_existing_fallback"
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
        f"IC1N41Ckb1fb{smoke['factor_base_points']}PDP4rootRCsampleLAgaussTDdirectISO0h"
        + sha(canonical(candidate))[:12])
    return candidate


def prepare() -> None:
    freeze_path = HERE / "FREEZE.json"
    assert not freeze_path.exists(), "Freeze is immutable"
    parent_path = PARENT / "FREEZE.json"
    pilot_path = PILOT / "FREEZE.json"
    parent = json.loads(parent_path.read_bytes())
    pilot = json.loads(pilot_path.read_bytes())
    prior_confirmation = json.loads((CONFIRMATION / "local/analysis.json").read_bytes())
    assert prior_confirmation["complete"] and prior_confirmation["verified_runs"] == 48
    selected_labels = [f"T{number:03}" for number in range(66, 90)]
    targets = [target for target in parent["targets"]
               if target["n"] == 41 and target["label"] in selected_labels]
    assert [target["label"] for target in targets] == selected_labels
    assert not set(selected_labels) & set(pilot["selected_labels"])
    assert not set(selected_labels) & set(json.loads(
        (CONFIRMATION / "FREEZE.json").read_bytes())["selected_labels"])
    assert len({target["workload_id"] for target in targets}) == 24
    for target in targets:
        point = PARENT / target["point_file_relative"]
        assert sha(point.read_bytes()) == target["point_sha256"]
        assert json.loads(point.read_bytes()) == target["public_point"]
        identity = json.loads((PARENT / "workloads" /
                               f"W{target['workload_id']}.json").read_bytes())
        assert identity["workload_id"] == target["workload_id"]
        assert identity["identity_sha256"] == target["workload_identity_sha256"]
        assert sha(canonical(identity["record"])) == identity["identity_sha256"]
        assert identity["record"]["target_point"] == target["public_point"]

    smoke_raw = SMOKE.with_suffix(".json").read_bytes()
    smoke_stdout = SMOKE.with_suffix(".stdout").read_bytes()
    assert smoke_raw == smoke_stdout
    smoke = json.loads(smoke_raw)
    control = next(target for target in parent["targets"]
                   if target["n"] == 41 and target["label"] == "T025")
    assert smoke["group_verified"] is True
    assert smoke["target"] == control["public_point"]
    assert smoke["recovered_scalar"] == control["fixture_scalar"]
    assert smoke["final_rank_complete"] is True and smoke["fallback_attempts"] == 0
    assert smoke["target_query_workers"] == 1 and smoke["relation_collection_workers"] == 14
    assert smoke["target_search_schedule"] == "native_rotating_state_order"
    assert smoke["factor_base_points"] == pilot["arms"]["baseline"]["factor_base_points"]
    assert smoke["factor_base_digest"] == pilot["arms"]["baseline"]["factor_base_digest"]
    source_arm = pilot["arms"]["expanded"]
    source_path = ROOT / "examples" / f"{source_arm['binary_name']}.rs"
    binary = BINARIES / source_arm["binary_name"]
    assert sha(source_path.read_bytes()) == source_arm["source_sha256"]
    assert sha(binary.read_bytes()) == source_arm["binary_sha256"]
    source_manifest = PILOT / source_arm["candidate_manifest_relative"]
    assert sha(source_manifest.read_bytes()) == source_arm["candidate_manifest_sha256"]
    baseline_candidate = make_baseline(json.loads(source_manifest.read_bytes()), smoke)
    baseline_manifest = HERE / "identities/candidates" / f"{baseline_candidate['candidate_id']}.json"
    write_json(baseline_manifest, baseline_candidate)
    baseline = {
        "binary_name": source_arm["binary_name"],
        "binary_sha256": source_arm["binary_sha256"],
        "source_sha256": source_arm["source_sha256"],
        "candidate_id": baseline_candidate["candidate_id"],
        "candidate_manifest_relative": str(baseline_manifest.relative_to(HERE)),
        "candidate_manifest_sha256": sha(baseline_manifest.read_bytes()),
        "columns": 1200, "prepare_rank": 1200,
        "factor_base_points": smoke["factor_base_points"],
        "factor_base_digest": smoke["factor_base_digest"],
        "root_table_slots": smoke["root_table_slots"],
    }
    expanded = copy.deepcopy(source_arm)
    expanded_manifest = PILOT / expanded["candidate_manifest_relative"]
    expanded["candidate_manifest_relative"] = "../" + str(
        expanded_manifest.relative_to(ROOT / "experiments"))
    run_order = []
    for index, label in enumerate(selected_labels):
        order = ("baseline", "expanded") if index % 2 == 0 else ("expanded", "baseline")
        run_order.extend({"label": label, "arm": arm, "repetition": 1}
                         for arm in order)
    freeze = {
        "kind": "s3_n41_c1200_c1800_matchedsource24_v1",
        "source_tree_git_oid": subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip(),
        "cargo_lock_sha256": sha((ROOT / "Cargo.lock").read_bytes()),
        "parent_freeze_sha256": sha(parent_path.read_bytes()),
        "pilot_freeze_sha256": sha(pilot_path.read_bytes()),
        "pilot_analysis_sha256": sha((PILOT / "local/analysis.json").read_bytes()),
        "pilot_sage_replay_sha256": sha((PILOT / "local/sage_replay.json").read_bytes()),
        "prior_confirmation_analysis_sha256": sha(
            (CONFIRMATION / "local/analysis.json").read_bytes()),
        "input_law": parent["input_law"],
        "selected_labels": selected_labels,
        "selected_targets": targets,
        "run_order": run_order,
        "planned_run_count": len(run_order),
        "relation_seed": pilot["relation_seed"],
        "relation_workers": 14, "target_workers": 1,
        "per_run_timeout_seconds": 180,
        "observed_peak_rss_envelope_bytes": 32 * 1024**3,
        "repetitions_per_arm_per_target": 1,
        "one_public_target_per_process": True,
        "measured_interval": pilot["measured_interval"],
        "comparison_variable": "actual base size, rank, and S3 index size; same source, executable, target kernel, and N41 native search schedule",
        "development_gate_online_ratio_at_least": "2.0",
        "matched_source_gate_bootstrap_lower_at_least": "2.0",
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
        "smoke_report_sha256": sha(smoke_raw),
        "smoke_archive_relative": "local/smoke_raw.tar.zst",
        "smoke_archive_sha256": sha((HERE / "local/smoke_raw.tar.zst").read_bytes()),
        "runner_sha256": sha((HERE / "run.py").read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze.py").read_bytes()),
        "preparer_sha256": sha(Path(__file__).read_bytes()),
        "arms": {"baseline": baseline, "expanded": expanded},
    }
    write_json(freeze_path, freeze)
    print(json.dumps({"freeze_sha256": sha(freeze_path.read_bytes()),
                      "baseline_candidate_id": baseline["candidate_id"],
                      "expanded_candidate_id": expanded["candidate_id"],
                      "target_count": len(targets), "planned_runs": len(run_order)}))


if __name__ == "__main__":
    prepare()
