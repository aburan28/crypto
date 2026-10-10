"""Freeze 34 fresh N53 C1200/table versus C1800/PMULL-fold pairs."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-c1800-confirmation34-20261010"
MATCHED = ROOT / "experiments/koblitz-s3-c1800-matchedsource-n53-20261010"
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
    panel_path = PARENT / "PANEL_FREEZE.json"
    panel = json.loads(panel_path.read_bytes())
    assert panel["kind"] == "s3_n53_c1800_confirmation34_panel_freeze_v1"
    labels = [f"T{index:03}" for index in range(213, 247)]
    assert panel["target_labels"] == labels
    targets = panel["targets"]
    assert [target["label"] for target in targets] == labels
    assert len({target["workload_id"] for target in targets}) == 34
    for target in targets:
        assert target["n"] == 53
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

    matched = json.loads((MATCHED / "FREEZE.json").read_bytes())
    fold_screen = json.loads((FOLD_SCREEN / "FREEZE.json").read_bytes())
    assert json.loads((MATCHED / "local/analysis.json").read_bytes())["complete"]
    assert json.loads((MATCHED / "local/sage_replay.json").read_bytes())[
        "all_scalars_independently_replayed"] is True
    assert json.loads((FOLD_SCREEN / "local/analysis.json").read_bytes())["complete"]
    assert not set(labels) & set(matched["selected_labels"])
    assert panel["relation_seed"] == matched["relation_seed"]
    arms = copy.deepcopy(matched["arms"])
    fold = fold_screen["curves"]["53"]["arms"]["fold"]
    expanded = arms["expanded"]
    expanded["candidate_id"] = fold["candidate_id"]
    expanded["candidate_manifest_relative"] = str(
        (FOLD_SCREEN / fold["candidate_manifest_relative"]).resolve().relative_to(
            ROOT / "experiments"))
    expanded["candidate_manifest_relative"] = "../" + expanded["candidate_manifest_relative"]
    expanded["candidate_manifest_sha256"] = fold["candidate_manifest_sha256"]
    for arm, record in arms.items():
        if arm == "baseline":
            manifest_path = (MATCHED / record["candidate_manifest_relative"]).resolve()
            record["candidate_manifest_relative"] = "../" + str(
                manifest_path.relative_to(ROOT / "experiments"))
        else:
            manifest_path = (HERE / record["candidate_manifest_relative"]).resolve()
        assert sha(manifest_path.read_bytes()) == record["candidate_manifest_sha256"]
        candidate = json.loads(manifest_path.read_bytes())
        candidate_id = candidate.pop("candidate_id")
        assert candidate_id == record["candidate_id"]
        assert candidate_id.endswith("h" + sha(canonical(candidate))[:12])
        binary = binaries / record["binary_name"]
        source = ROOT / "examples" / f"{record['binary_name']}.rs"
        assert sha(binary.read_bytes()) == record["binary_sha256"]
        assert sha(source.read_bytes()) == record["source_sha256"]
    assert arms["baseline"]["binary_sha256"] == arms["expanded"]["binary_sha256"]
    assert arms["baseline"]["source_sha256"] == arms["expanded"]["source_sha256"]
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
    shuffled = labels[:]
    random.Random(20261010).shuffle(shuffled)
    run_order = []
    for index, label in enumerate(shuffled):
        order = ("baseline", "expanded") if index % 2 == 0 else ("expanded", "baseline")
        run_order.extend({"label": label, "arm": arm, "repetition": 1}
                         for arm in order)
    smoke = (MATCHED / matched["smoke_archive_relative"]).resolve()
    freeze = {
        "kind": "s3_n53_c1200_c1800_pmullfold_confirmation34_v1",
        "source_tree_git_oid": subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip(),
        "cargo_lock_sha256": sha((ROOT / "Cargo.lock").read_bytes()),
        "parent_panel_freeze_sha256": sha(panel_path.read_bytes()),
        "matched_freeze_sha256": sha((MATCHED / "FREEZE.json").read_bytes()),
        "matched_analysis_sha256": sha((MATCHED / "local/analysis.json").read_bytes()),
        "matched_sage_replay_sha256": sha((MATCHED / "local/sage_replay.json").read_bytes()),
        "fold_screen_freeze_sha256": sha((FOLD_SCREEN / "FREEZE.json").read_bytes()),
        "fold_screen_analysis_sha256": sha((FOLD_SCREEN / "local/analysis.json").read_bytes()),
        "probe_source_sha256": sha((FOLD_SCREEN / "probe.rs").read_bytes()),
        "probe_binary_sha256": sha(probe.read_bytes()),
        "probe_expected_stdout": probe_expected,
        "arms_env": {"baseline": "0", "expanded": "1"},
        "host_os": platform.system(), "host_architecture": platform.machine(),
        "selected_labels": labels, "selected_targets": targets,
        "run_order": run_order, "run_order_policy": "seeded target shuffle then alternating AB and BA",
        "run_order_seed": 20261010, "planned_run_count": len(run_order),
        "relation_seed": matched["relation_seed"],
        "relation_workers": matched["relation_workers"],
        "target_workers": matched["target_workers"],
        "per_run_timeout_seconds": matched["per_run_timeout_seconds"],
        "observed_peak_rss_envelope_bytes": matched["observed_peak_rss_envelope_bytes"],
        "repetitions_per_arm_per_target": 1,
        "one_public_target_per_process": True,
        "measured_interval": matched["measured_interval"],
        "comparison_variable": "C1200/table versus C1800/ARM64 PMULL sparse-tail fold in one source and executable; exact N53 curve and target paired",
        "development_gate_online_ratio_at_least": "2.0",
        "confirmation_gate_bootstrap_lower_at_least": "2.0",
        "bootstrap_resamples": 20000, "bootstrap_seed": 20261010,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
        "smoke_archive_relative": "../" + str(smoke.relative_to(ROOT / "experiments")),
        "smoke_archive_sha256": sha(smoke.read_bytes()),
        "runner_sha256": sha((HERE / "run.py").read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze.py").read_bytes()),
        "replay_sha256": sha((HERE / "replay_sage.py").read_bytes()),
        "preparer_sha256": sha(Path(__file__).read_bytes()),
        "arms": arms,
    }
    assert freeze["planned_run_count"] == 68
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
