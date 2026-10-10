"""Freeze public-only inputs for a future Modal N53 confirmation run."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PANEL = HERE.parent / "PANEL_FREEZE.json"
EXECUTION = HERE.parent / "EXECUTION_FREEZE.json"
GENERATED = HERE.parent / "GENERATED_TARGETS.json"
BUNDLE = HERE / "PUBLIC_BUNDLE.json"
RECEIPT = HERE / "PUBLIC_BUNDLE_RECEIPT.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def prepare(verify: bool = False) -> dict:
    if not verify and (BUNDLE.exists() or RECEIPT.exists()):
        raise FileExistsError("public bundle or receipt already exists")
    panel_raw = PANEL.read_bytes()
    panel = json.loads(panel_raw)
    execution_raw = EXECUTION.read_bytes()
    execution = json.loads(execution_raw)
    generated_raw = GENERATED.read_bytes()
    assert panel["kind"] == "s3_n53_c1800_confirmation34_panel_freeze_v1"
    assert execution["kind"] == "s3_n53_c1800_confirmation34_execution_freeze_v1"
    assert execution["panel_freeze_sha256"] == sha(panel_raw)
    assert execution["generated_targets_sha256"] == sha(generated_raw)
    assert panel["generated_targets_sha256"] == sha(generated_raw)
    assert len(panel["targets"]) == 34
    assert len(panel["run_order"]) == 136
    assert len(set(panel["target_labels"])) == 34
    assert panel["relation_workers"] == 14 and panel["target_workers"] == 1
    assert panel["repetitions_per_arm_per_target"] == 2
    assert panel["observed_peak_rss_envelope_bytes"] == 34359738368
    assert {row["label"] for row in panel["targets"]} == set(panel["target_labels"])
    assert {row["n"] for row in panel["targets"]} == {53}

    arms = {}
    for name in ("baseline", "combined"):
        frozen = panel["arms"][name]
        binary = execution["binaries"][name]
        assert binary["sha256"] and binary["architecture"] == "x86_64-linux-elf"
        assert sha((ROOT / frozen["source_relative"]).read_bytes()) == frozen[
            "source_sha256"]
        assert sha((ROOT / frozen["candidate_manifest_relative"]).read_bytes()) == frozen[
            "candidate_manifest_sha256"]
        assert frozen["target_workers"] == 1
        arms[name] = {
            "binary_sha256": binary["sha256"],
            "candidate_id": frozen["candidate_id"],
            "columns": frozen["columns"],
            "factor_base_digest": frozen["factor_base_digest"],
            "factor_base_points": frozen["factor_base_points"],
            "prepare_rank": frozen["prepare_rank"],
            "root_table_capacity_policy": frozen["root_table_capacity_policy"],
            "root_table_layout": frozen["root_table_layout"],
            "root_table_slots": frozen["root_table_slots"],
            "target_search_prefix_columns": frozen["target_search_prefix_columns"],
            "target_search_schedule": frozen["target_search_schedule"],
        }

    targets = []
    seen_points = set()
    for frozen in panel["targets"]:
        point_raw = (ROOT / frozen["point_file_relative"]).read_bytes()
        workload_raw = (ROOT / frozen["workload_manifest_relative"]).read_bytes()
        assert sha(point_raw) == frozen["point_sha256"]
        assert json.loads(point_raw) == frozen["public_point"]
        assert sha(workload_raw) == frozen["workload_manifest_sha256"]
        workload = json.loads(workload_raw)
        assert workload["workload_id"] == frozen["workload_id"]
        assert workload["record"]["target_point"] == frozen["public_point"]
        point = tuple(frozen["public_point"])
        assert point not in seen_points
        seen_points.add(point)
        targets.append({
            "label": frozen["label"],
            "n": 53,
            "public_point": frozen["public_point"],
            "point_sha256": frozen["point_sha256"],
            "workload_id": frozen["workload_id"],
        })
    schedule = []
    counts = {}
    for row in panel["run_order"]:
        assert row["label"] in panel["target_labels"]
        assert row["arm"] in arms and row["repetition"] in (1, 2)
        key = (row["label"], row["arm"], row["repetition"])
        counts[key] = counts.get(key, 0) + 1
        schedule.append({"label": row["label"], "arm": row["arm"],
                         "repetition": row["repetition"]})
    assert len(counts) == 136 and all(count == 1 for count in counts.values())

    bundle = {
        "kind": "s3_modal_c1800_confirmation34_public_bundle_v1",
        "panel_freeze_sha256": sha(panel_raw),
        "execution_freeze_sha256": sha(execution_raw),
        "generated_targets_sha256": sha(generated_raw),
        "source_tree_git_oid": execution["linux_build_source_tree_git_oid"],
        "arms": arms,
        "targets": targets,
        "target_count": 34,
        "run_order": schedule,
        "planned_run_count": 136,
        "relation_seed": panel["relation_seed"],
        "relation_workers": panel["relation_workers"],
        "target_workers": panel["target_workers"],
        "memory_envelope_bytes": panel["observed_peak_rss_envelope_bytes"],
        "per_process_timeout_seconds": panel["timeout_seconds_per_process"],
    }
    raw = canonical(bundle) + b"\n"
    assert b"fixture_scalar" not in raw
    assert b"scalar_derivation" not in raw
    assert all(str(row["fixture_scalar"]).encode() not in raw
               for row in panel["targets"])
    receipt = {
        "kind": "s3_modal_c1800_confirmation34_public_bundle_receipt_v1",
        "prepare_source_sha256": sha(Path(__file__).read_bytes()),
        "panel_freeze_sha256": sha(panel_raw),
        "execution_freeze_sha256": sha(execution_raw),
        "public_bundle_sha256": sha(raw),
        "public_bundle_bytes": len(raw),
        "target_count": 34,
        "planned_run_count": 136,
        "fixture_scalars_in_bundle": False,
        "source_repository_in_bundle": False,
        "upload_performed_by_preparation_script": False,
    }
    if verify:
        assert BUNDLE.read_bytes() == raw
        assert json.loads(RECEIPT.read_bytes()) == receipt
    else:
        BUNDLE.write_bytes(raw)
        RECEIPT.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    print(json.dumps(prepare(args.verify), sort_keys=True))
