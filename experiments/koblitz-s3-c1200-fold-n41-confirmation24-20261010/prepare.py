"""Freeze N41 old C1200/table versus improved C1200/PMULL-fold pairs."""

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
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009"
TARGETKERNEL = ROOT / "experiments/koblitz-s3-targetkernel-n41-20261010"
FOLD_SCREEN = ROOT / "experiments/koblitz-s3-pmull-fold-screen-20261010"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(record: dict) -> bytes:
    return json.dumps(record, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def prepare(binaries: Path, probe: Path) -> None:
    freeze_path = HERE / "FREEZE.json"
    assert not freeze_path.exists(), "Freeze is immutable"
    assert platform.system() == "Darwin" and platform.machine() == "arm64"
    parent_path = PARENT / "FREEZE.json"
    parent = json.loads(parent_path.read_bytes())
    generated_path = HERE / "GENERATED_TARGETS.json"
    generated = json.loads(generated_path.read_bytes())
    assert generated["kind"] == "s3_n41_c1200_fold_confirmation24_public_targets_v1"
    assert generated["parent_freeze_sha256"] == sha(parent_path.read_bytes())
    assert generated["generator_order_verified"] is True
    assert generated["existing_t200_fixture_replayed"] is True
    labels = [f"T{index:03}" for index in range(186, 210)]
    targets = generated["targets"]
    assert [target["label"] for target in targets] == labels
    assert len({target["workload_id"] for target in targets}) == 24
    parent_by_label = {target["label"]: target for target in parent["targets"]
                       if target["n"] == 41}
    for target in targets:
        assert target["n"] == 41
        if target["label"] in parent_by_label:
            prior = parent_by_label[target["label"]]
            for key in ("fixture_scalar", "scalar_derivation_sha256", "public_point",
                        "workload_id", "workload_identity_sha256"):
                assert target[key] == prior[key]
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

    kernel = json.loads((TARGETKERNEL / "FREEZE.json").read_bytes())
    fold_screen = json.loads((FOLD_SCREEN / "FREEZE.json").read_bytes())
    assert json.loads((TARGETKERNEL / "local/analysis.json").read_bytes())["complete"]
    assert json.loads((TARGETKERNEL / "local/sage_replay.json").read_bytes())[
        "all_scalars_independently_replayed"] is True
    assert json.loads((FOLD_SCREEN / "local/analysis.json").read_bytes())["complete"]
    assert not set(labels) & set(kernel["selected_labels"])
    arms = {"baseline": copy.deepcopy(kernel["arms"]["baseline"]),
            "expanded": copy.deepcopy(kernel["arms"]["improved"])}
    fold = fold_screen["curves"]["41"]["arms"]["fold"]
    arms["expanded"]["candidate_id"] = fold["candidate_id"]
    arms["expanded"]["candidate_manifest_sha256"] = fold["candidate_manifest_sha256"]
    fold_manifest = (FOLD_SCREEN / fold["candidate_manifest_relative"]).resolve()
    arms["expanded"]["candidate_manifest_relative"] = os.path.relpath(fold_manifest, HERE)
    for arm, record in arms.items():
        if arm == "baseline":
            path = (TARGETKERNEL / record["candidate_manifest_relative"]).resolve()
            record["candidate_manifest_relative"] = os.path.relpath(path, HERE)
        else:
            path = fold_manifest
        assert sha(path.read_bytes()) == record["candidate_manifest_sha256"]
        candidate = json.loads(path.read_bytes())
        candidate_id = candidate.pop("candidate_id")
        assert candidate_id == record["candidate_id"]
        assert candidate_id.endswith("h" + sha(canonical(candidate))[:12])
        binary = binaries / record["binary_name"]
        source = ROOT / "examples" / f"{record['binary_name']}.rs"
        assert sha(binary.read_bytes()) == record["binary_sha256"]
        assert sha(source.read_bytes()) == record["source_sha256"]
    for key in ("columns", "factor_base_points", "factor_base_digest", "prepare_rank",
                "root_table_slots"):
        assert arms["baseline"][key] == arms["expanded"][key]
    assert sha(probe.read_bytes()) == fold_screen["probe_binary_sha256"]
    probe_expected = {}
    for arm, flag in (("baseline", "0"), ("expanded", "1")):
        env = dict(os.environ)
        env["GF2_ARM_PMULL_FOLD"] = flag
        probe_expected[arm] = subprocess.check_output([str(probe)], env=env, text=True)
    assert probe_expected == {
        "baseline": "n=41 kernel=pmull folds=false\nn=53 kernel=pmull folds=false\n",
        "expanded": "n=41 kernel=pmull folds=true\nn=53 kernel=pmull folds=true\n",
    }
    run_order = []
    for index, label in enumerate(labels):
        order = ("baseline", "expanded") if index % 2 == 0 else ("expanded", "baseline")
        run_order.extend({"label": label, "arm": arm, "repetition": 1}
                         for arm in order)
    freeze = {
        "kind": "s3_n41_c1200_targetkernel_pmullfold_confirmation24_v1",
        "source_tree_git_oid": subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip(),
        "cargo_lock_sha256": sha((ROOT / "Cargo.lock").read_bytes()),
        "parent_freeze_sha256": sha(parent_path.read_bytes()),
        "generated_targets_sha256": sha(generated_path.read_bytes()),
        "generator_sha256": sha((HERE / "make_targets_sage.py").read_bytes()),
        "sage_runtime_info_sha256": sha((HERE / "local/sage_runtime_info.json").read_bytes()),
        "targetkernel_freeze_sha256": sha((TARGETKERNEL / "FREEZE.json").read_bytes()),
        "targetkernel_analysis_sha256": sha((TARGETKERNEL / "local/analysis.json").read_bytes()),
        "targetkernel_sage_replay_sha256": sha((TARGETKERNEL / "local/sage_replay.json").read_bytes()),
        "fold_screen_freeze_sha256": sha((FOLD_SCREEN / "FREEZE.json").read_bytes()),
        "fold_screen_analysis_sha256": sha((FOLD_SCREEN / "local/analysis.json").read_bytes()),
        "probe_source_sha256": sha((FOLD_SCREEN / "probe.rs").read_bytes()),
        "probe_binary_sha256": sha(probe.read_bytes()),
        "probe_expected_stdout": probe_expected,
        "arms_env": {"baseline": "0", "expanded": "1"},
        "host_os": platform.system(), "host_architecture": platform.machine(),
        "input_law": parent["input_law"],
        "selected_labels": labels, "selected_targets": targets,
        "run_order": run_order, "planned_run_count": len(run_order),
        "relation_seed": kernel["relation_seed"],
        "relation_workers": kernel["relation_workers"],
        "target_workers": kernel["target_workers"],
        "per_run_timeout_seconds": kernel["per_run_timeout_seconds"],
        "observed_peak_rss_envelope_bytes": kernel["observed_peak_rss_envelope_bytes"],
        "repetitions_per_arm_per_target": 1,
        "one_public_target_per_process": True,
        "measured_interval": kernel["measured_interval"],
        "comparison_variable": "old C1200/table versus improved C1200/ARM PMULL fold with identical N41 curve, factor base, rank, root capacity, target and resources",
        "development_gate_online_ratio_at_least": "2.0",
        "confirmation_gate_bootstrap_lower_at_least": "2.0",
        "bootstrap_resamples": 20000, "bootstrap_seed": 20261010,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
        "runner_sha256": sha((HERE / "run.py").read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze.py").read_bytes()),
        "replay_sha256": sha((HERE / "replay_sage.py").read_bytes()),
        "preparer_sha256": sha(Path(__file__).read_bytes()),
        "arms": arms,
    }
    assert freeze["planned_run_count"] == 48
    freeze_path.write_text(json.dumps(freeze, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"freeze_sha256": sha(freeze_path.read_bytes()),
                      "selected_targets": [labels[0], labels[-1]],
                      "candidate_ids": {arm: record["candidate_id"] for arm, record in arms.items()},
                      "planned_runs": len(run_order)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binaries", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.binaries.resolve(), args.probe.resolve())
