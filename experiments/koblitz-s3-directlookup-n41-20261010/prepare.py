"""Freeze the N41 C1200 PMULL prefetch versus direct-lookup screen."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-c1200-fold-n41-confirmation48-20261010"
DIRECT_NAME = "s3_partial_prepared_fullrank_projective_directlookup_targetspan"
LABELS = [f"T{index:03}" for index in range(210, 222)]


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: dict) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def replace_source(value, before: str, after: str):
    if isinstance(value, dict):
        return {key: replace_source(item, before, after) for key, item in value.items()}
    if isinstance(value, list):
        return [replace_source(item, before, after) for item in value]
    return after if value == before else value


def main(direct_binary: Path, prefetch_binary: Path, probe: Path,
         timeout_seconds: int, rss_envelope_bytes: int,
         preview: bool) -> None:
    freeze_path = HERE / "FREEZE.json"
    assert not freeze_path.exists(), "The freeze is immutable"
    assert platform.system() == "Darwin" and platform.machine() == "arm64"
    assert timeout_seconds > 0 and rss_envelope_bytes > 0
    parent_path = PARENT / "FREEZE.json"
    parent = json.loads(parent_path.read_bytes())
    assert json.loads((PARENT / "local/analysis.json").read_bytes())["complete"]
    assert json.loads((PARENT / "local/sage_replay.json").read_bytes())[
        "all_scalars_independently_replayed"] is True
    selected = [target for target in parent["selected_targets"]
                if target["label"] in LABELS]
    assert [target["label"] for target in selected] == LABELS
    assert len({target["workload_id"] for target in selected}) == len(LABELS)
    for target in selected:
        point = ROOT / target["point_file_relative"]
        assert sha(point.read_bytes()) == target["point_sha256"]
        assert json.loads(point.read_bytes()) == target["public_point"]
        workload = ROOT / target["workload_manifest_relative"]
        assert sha(workload.read_bytes()) == target["workload_manifest_sha256"]
        identity = json.loads(workload.read_bytes())
        assert identity["workload_id"] == target["workload_id"]
        assert sha(canonical(identity["record"])) == identity["identity_sha256"]

    baseline = copy.deepcopy(parent["arms"]["expanded"])
    assert sha(prefetch_binary.read_bytes()) == baseline["binary_sha256"]
    old_source = ROOT / "examples" / f"{baseline['binary_name']}.rs"
    assert sha(old_source.read_bytes()) == baseline["source_sha256"]
    old_manifest_path = (PARENT / baseline["candidate_manifest_relative"]).resolve()
    assert sha(old_manifest_path.read_bytes()) == baseline["candidate_manifest_sha256"]
    baseline["candidate_manifest_relative"] = os.path.relpath(old_manifest_path, HERE)

    new_source = ROOT / "examples" / f"{DIRECT_NAME}.rs"
    new_source_sha = sha(new_source.read_bytes())
    new_binary_sha = sha(direct_binary.read_bytes())
    assert new_source_sha != baseline["source_sha256"]
    assert new_binary_sha != baseline["binary_sha256"]
    source_commit = subprocess.check_output(
        ["git", "-C", str(ROOT), "log", "-1", "--format=%H", "--",
         f"examples/{DIRECT_NAME}.rs"], text=True).strip()
    assert len(source_commit) == 40
    committed_source = subprocess.check_output(
        ["git", "-C", str(ROOT), "show",
         f"{source_commit}:examples/{DIRECT_NAME}.rs"])
    assert sha(committed_source) == new_source_sha
    candidate = json.loads(old_manifest_path.read_bytes())
    old_candidate_id = candidate.pop("candidate_id")
    assert old_candidate_id == baseline["candidate_id"]
    candidate = replace_source(candidate, baseline["source_sha256"], new_source_sha)
    flags = candidate["implementation"]["algorithm_flags"]
    flags["target_lookup_prefetch_policy"] = "none_immediate_root_order_lookup"
    candidate["point_decomposition"]["target_lookup_prefetch_policy"] = (
        "none_immediate_root_order_lookup")
    candidate["point_decomposition"]["implementation"] = (
        "parallel_sharded_s3_index_prelifted_group_sum_roots_batch8_inverse_"
        "target_conjugates_dense_packed_u40_root_slots_compact_poly_states_"
        "direct_lookup_projective_target_checks")
    candidate["implementation"]["frozen_source_commit"] = source_commit
    assert candidate["implementation"]["example_source_sha256"] == new_source_sha
    candidate_id = old_candidate_id.rsplit("h", 1)[0] + "h" + sha(canonical(candidate))[:12]
    candidate["candidate_id"] = candidate_id
    manifest_path = HERE / "identities/candidates" / f"{candidate_id}.json"
    manifest_bytes = canonical(candidate) + b"\n"
    if not preview:
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        assert not manifest_path.exists()
        manifest_path.write_bytes(manifest_bytes)

    direct = copy.deepcopy(baseline)
    direct.update({
        "binary_name": DIRECT_NAME,
        "binary_sha256": new_binary_sha,
        "source_sha256": new_source_sha,
        "candidate_id": candidate_id,
        "candidate_manifest_relative": os.path.relpath(manifest_path, HERE),
        "candidate_manifest_sha256": sha(manifest_bytes),
    })
    for key in ("columns", "factor_base_points", "factor_base_digest",
                "prepare_rank", "root_table_slots"):
        assert direct[key] == baseline[key]
    expected_probe = subprocess.check_output(
        [str(probe)], env={**os.environ, "GF2_ARM_PMULL_FOLD": "1"}, text=True)
    assert expected_probe == (
        "n=41 kernel=pmull folds=true\nn=53 kernel=pmull folds=true\n")
    order = []
    for index, label in enumerate(LABELS):
        arms = ("prefetch", "direct") if index % 2 == 0 else ("direct", "prefetch")
        order.extend({"label": label, "arm": arm, "repetition": 1} for arm in arms)
    freeze = {
        "kind": "s3_n41_c1200_directlookup_screen_v1",
        "comparison_variable": "immediate root-order lookup versus eight-state speculative prefetch at fixed N41 C1200 base and PMULL fold",
        "screen_status": "exploratory_reused_public_points",
        "source_tree_git_oid": parent["source_tree_git_oid"],
        "source_commit": source_commit,
        "cargo_lock_sha256": sha((ROOT / "Cargo.lock").read_bytes()),
        "parent_freeze_sha256": sha(parent_path.read_bytes()),
        "parent_analysis_sha256": sha((PARENT / "local/analysis.json").read_bytes()),
        "parent_sage_replay_sha256": sha((PARENT / "local/sage_replay.json").read_bytes()),
        "sage_runtime_info_sha256": sha((HERE / "local/sage_runtime_info.json").read_bytes()),
        "probe_binary_sha256": sha(probe.read_bytes()),
        "probe_expected_stdout": expected_probe,
        "host_os": platform.system(), "host_architecture": platform.machine(),
        "selected_labels": LABELS, "selected_targets": selected,
        "run_order": order, "planned_run_count": len(order),
        "arms": {"prefetch": baseline, "direct": direct},
        "arms_env": {"prefetch": "1", "direct": "1"},
        "relation_seed": parent["relation_seed"],
        "relation_workers": parent["relation_workers"],
        "target_workers": parent["target_workers"],
        "per_run_timeout_seconds": timeout_seconds,
        "observed_peak_rss_envelope_bytes": rss_envelope_bytes,
        "one_public_target_per_process": True,
        "repetitions_per_arm_per_target": 1,
        "measured_interval": parent["measured_interval"],
        "bootstrap_resamples": 20000,
        "bootstrap_seed": 20261010,
        "development_gate_online_ratio_at_least": "2.0",
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
        "runner_sha256": sha((HERE / "run.py").read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze.py").read_bytes()),
        "replay_sha256": sha((HERE / "replay_sage.py").read_bytes()),
        "preparer_sha256": sha(Path(__file__).read_bytes()),
    }
    assert freeze["planned_run_count"] == 24
    freeze_bytes = (json.dumps(freeze, indent=2, sort_keys=True) + "\n").encode()
    if not preview:
        freeze_path.write_bytes(freeze_bytes)
    print(json.dumps({"freeze_sha256": sha(freeze_bytes), "preview": preview,
                      "candidate_id": candidate_id,
                      "selected_targets": [LABELS[0], LABELS[-1]],
                      "planned_runs": len(order)}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--direct-binary", type=Path, required=True)
    parser.add_argument("--prefetch-binary", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--per-run-timeout-seconds", type=int, required=True)
    parser.add_argument("--rss-envelope-bytes", type=int, required=True)
    parser.add_argument("--preview", action="store_true")
    args = parser.parse_args()
    main(args.direct_binary.resolve(), args.prefetch_binary.resolve(),
         args.probe.resolve(), args.per_run_timeout_seconds,
         args.rss_envelope_bytes, args.preview)
