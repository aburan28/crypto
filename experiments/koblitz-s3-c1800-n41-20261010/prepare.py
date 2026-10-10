"""Freeze a fresh N41 C1200/C1800 S3 screen from checked smoke reports."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009"
N41_TEMPLATE = ROOT / "experiments/koblitz-s3-prepared-20261009/identities/candidates/IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h3f64f5604376.json"
TEMPLATES = {
    "baseline": ROOT / "experiments/koblitz-s3-c2800-mixed-20261009/identities/candidates/IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0hf8b850f65ece.json",
    "expanded": ROOT / "experiments/koblitz-s3-c1800-combined-20261010/identities/candidates/IC1N53Ckb1fb190800PDP4rootRCsampleLAgaussTDdirectISO0h1ab23390e700.json",
}
BINARIES = Path("/private/tmp/crypto-s3-release-examples-20261010")
SMOKE = Path("/private/tmp/s3-n41-c1800-smoke-20261010")
BINARY_NAMES = {
    "baseline": "s3_partial_prepared_fullrank_c2800mixed_targetspan",
    "expanded": "s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan",
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(record: dict) -> bytes:
    return json.dumps(record, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def write_json(path: Path, record: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")


def make_candidate(arm: str, report: dict, n41: dict) -> dict:
    candidate = copy.deepcopy(json.loads(TEMPLATES[arm].read_bytes()))
    candidate.pop("candidate_id")
    for key in ("field", "curve", "endomorphism"):
        candidate[key] = copy.deepcopy(n41[key])
    base = candidate["factor_base"]
    base.update({
        "actual_usable_points_B": report["factor_base_points"],
        "effective_relation_columns": report["orbit_columns"],
        "enumerated_set_digest": report["factor_base_digest"],
        "last_scanned_x": report["factor_base_scanned_x"],
        "requested_orbit_columns": report["orbit_columns"],
        "start_x": report["factor_base_scan_start_x"],
    })
    point = candidate["point_decomposition"]
    point.update({
        "indexed_state_payload_bytes": report["index_state_payload_bytes"],
        "root_table_payload_bytes": report["root_table_slots"] * 13,
        "root_table_slots": report["root_table_slots"],
        "target_conjugate_count": 41,
        "target_query_policy": "point_seeded_native_rotating_state_order_with_first_verified_relation_and_worker_drain",
        "target_search_prefix_columns": report["orbit_columns"],
        "target_search_schedule": "native_rotating_state_order",
    })
    relation = candidate["relation_collection"]
    relation["pre_target_rank_threshold"] = report["orbit_columns"]
    relation["stop_criterion"] = (
        f"complete_{report['orbit_columns']}_rank_before_target_input_if_attainable_then_online_span_fallback_if_incomplete")
    algebra = candidate["relation_linear_algebra"]
    algebra["modulus"] = report["subgroup_order"]
    algebra["rank_criterion"] = (
        f"full_{report['orbit_columns']}_pre_target_if_attained_then_existing_fallback")
    flags = candidate["implementation"]["algorithm_flags"]
    flags.update({
        "columns": report["orbit_columns"],
        "index_state_payload_bytes": report["index_state_payload_bytes"],
        "inverse_frobenius_target_conjugates": 41,
        "prepare_rank": report["orbit_columns"],
        "root_table_payload_bytes": report["root_table_slots"] * 13,
        "root_table_slots": report["root_table_slots"],
        "target_search_prefix_columns": report["orbit_columns"],
        "target_search_schedule": "native_rotating_state_order",
    })
    digest = sha(canonical(candidate))[:12]
    candidate["candidate_id"] = (
        f"IC1N41Ckb1fb{report['factor_base_points']}PDP4rootRCsampleLAgaussTDdirectISO0h{digest}")
    return candidate


def prepare() -> None:
    assert not (HERE / "FREEZE.json").exists(), "Freeze is immutable"
    parent = json.loads((PARENT / "FREEZE.json").read_bytes())
    n41 = json.loads(N41_TEMPLATE.read_bytes())
    selected_labels = [f"T{number:03}" for number in range(34, 42)]
    selected_targets = [target for target in parent["targets"]
                        if target["n"] == 41 and target["label"] in selected_labels]
    assert [target["label"] for target in selected_targets] == selected_labels
    assert len({target["workload_id"] for target in selected_targets}) == 8
    for target in selected_targets:
        point = PARENT / target["point_file_relative"]
        assert sha(point.read_bytes()) == target["point_sha256"]
        assert json.loads(point.read_bytes()) == target["public_point"]
        workload = PARENT / "workloads" / f"W{target['workload_id']}.json"
        identity = json.loads(workload.read_bytes())
        assert identity["workload_id"] == target["workload_id"]
        assert identity["identity_sha256"] == target["workload_identity_sha256"]
        assert sha(canonical(identity["record"])) == identity["identity_sha256"]
        assert identity["record"]["target_point"] == target["public_point"]
    smoke_target = next(target for target in parent["targets"]
                        if target["n"] == 41 and target["label"] == "T025")
    arms = {}
    smoke = {}
    for arm in ("baseline", "expanded"):
        binary_name = BINARY_NAMES[arm]
        source = ROOT / "examples" / f"{binary_name}.rs"
        report_path = SMOKE / f"{arm}.json"
        report = json.loads(report_path.read_bytes())
        template = json.loads(TEMPLATES[arm].read_bytes())
        assert sha(source.read_bytes()) == template["implementation"]["example_source_sha256"]
        assert report["group_verified"] is True
        assert report["n"] == 41 and report["a"] == 0
        assert report["target"] == smoke_target["public_point"]
        assert report["recovered_scalar"] == smoke_target["fixture_scalar"]
        assert report["prepared_rank"] == report["orbit_columns"]
        assert report["final_rank_complete"] is True
        assert report["fallback_attempts"] == 0
        assert report["relation_collection_workers"] == 14
        assert report["target_query_workers"] == 1
        assert report["target_search_schedule"] == "native_rotating_state_order"
        candidate = make_candidate(arm, report, n41)
        candidate_path = HERE / "identities/candidates" / f"{candidate['candidate_id']}.json"
        write_json(candidate_path, candidate)
        arms[arm] = {
            "binary_name": binary_name,
            "binary_sha256": sha((BINARIES / binary_name).read_bytes()),
            "source_sha256": sha(source.read_bytes()),
            "candidate_id": candidate["candidate_id"],
            "candidate_manifest_relative": str(candidate_path.relative_to(HERE)),
            "candidate_manifest_sha256": sha(candidate_path.read_bytes()),
            "columns": report["orbit_columns"],
            "prepare_rank": report["prepared_rank"],
            "factor_base_points": report["factor_base_points"],
            "factor_base_digest": report["factor_base_digest"],
            "root_table_slots": report["root_table_slots"],
        }
        smoke[arm] = {
            "report_sha256": sha(report_path.read_bytes()),
            "target_label": "T025",
            "group_verified": report["group_verified"],
            "prepared_rank": report["prepared_rank"],
            "target_online_ms": report["timing_ms"]["target_online_after_reusable_setup"],
        }
    run_order = []
    for index, label in enumerate(selected_labels):
        order = ("baseline", "expanded") if index % 2 == 0 else ("expanded", "baseline")
        run_order.extend({"label": label, "arm": arm, "repetition": 1}
                         for arm in order)
    freeze = {
        "kind": "s3_n41_c1200_c1800_fresh_screen_v1",
        "source_tree_git_oid": subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip(),
        "cargo_lock_sha256": sha((ROOT / "Cargo.lock").read_bytes()),
        "parent_freeze_sha256": sha((PARENT / "FREEZE.json").read_bytes()),
        "input_law": parent["input_law"],
        "selected_labels": selected_labels,
        "selected_targets": selected_targets,
        "run_order": run_order,
        "planned_run_count": len(run_order),
        "relation_seed": 20260928,
        "relation_workers": 14,
        "target_workers": 1,
        "per_run_timeout_seconds": 180,
        "observed_peak_rss_envelope_bytes": 32 * 1024**3,
        "repetitions_per_arm_per_target": 1,
        "one_public_target_per_process": True,
        "measured_interval": "target_online_after_reusable_setup; all target-dependent phases through independent scalar replay",
        "comparison_variable": "complete C1200 versus C1800 candidate, including factor-base, index, and source-specific target checks",
        "development_gate_online_ratio_at_least": "2.0",
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
        "runner_sha256": sha((HERE / "run.py").read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze.py").read_bytes()),
        "preparer_sha256": sha(Path(__file__).read_bytes()),
        "arms": arms,
        "smoke": smoke,
        "smoke_archive_relative": "local/smoke_raw.tar.gz",
        "smoke_archive_sha256": sha((HERE / "local/smoke_raw.tar.gz").read_bytes()),
    }
    write_json(HERE / "FREEZE.json", freeze)
    print(json.dumps({"freeze_sha256": sha((HERE / "FREEZE.json").read_bytes()),
                      "selected_labels": selected_labels,
                      "candidate_ids": {arm: data["candidate_id"] for arm, data in arms.items()}},
                     indent=2))


if __name__ == "__main__":
    prepare()
