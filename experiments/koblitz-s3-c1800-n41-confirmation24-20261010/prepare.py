"""Freeze 24 untouched N41 targets for C1200/C1800 confirmation."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PILOT = ROOT / "experiments/koblitz-s3-c1800-n41-20261010"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(record: dict) -> bytes:
    return json.dumps(record, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def prepare() -> None:
    freeze_path = HERE / "FREEZE.json"
    assert not freeze_path.exists(), "Freeze is immutable"
    pilot = json.loads((PILOT / "FREEZE.json").read_bytes())
    pilot_analysis = json.loads((PILOT / "local/analysis.json").read_bytes())
    pilot_replay = json.loads((PILOT / "local/sage_replay.json").read_bytes())
    assert pilot_analysis["complete"] and pilot_analysis["verified_runs"] == 16
    assert pilot_replay["all_scalars_independently_replayed"] is True
    parent = json.loads((PARENT / "FREEZE.json").read_bytes())
    labels = [f"T{number:03}" for number in range(42, 66)]
    targets = [target for target in parent["targets"]
               if target["n"] == 41 and target["label"] in labels]
    assert [target["label"] for target in targets] == labels
    assert len({target["workload_id"] for target in targets}) == 24
    assert not ({target["label"] for target in targets} &
                set(pilot["selected_labels"]))
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
    arms = copy.deepcopy(pilot["arms"])
    for arm in arms.values():
        candidate_path = PILOT / arm["candidate_manifest_relative"]
        assert sha(candidate_path.read_bytes()) == arm["candidate_manifest_sha256"]
        arm["candidate_manifest_relative"] = str(candidate_path.relative_to(ROOT / "experiments"))
        arm["candidate_manifest_relative"] = "../" + arm["candidate_manifest_relative"]
    run_order = []
    for index, label in enumerate(labels):
        order = ("baseline", "expanded") if index % 2 == 0 else ("expanded", "baseline")
        run_order.extend({"label": label, "arm": arm, "repetition": 1}
                         for arm in order)
    freeze = {
        "kind": "s3_n41_c1200_c1800_confirmation24_v1",
        "source_tree_git_oid": subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip(),
        "cargo_lock_sha256": sha((ROOT / "Cargo.lock").read_bytes()),
        "parent_freeze_sha256": sha((PARENT / "FREEZE.json").read_bytes()),
        "pilot_freeze_sha256": sha((PILOT / "FREEZE.json").read_bytes()),
        "pilot_analysis_sha256": sha((PILOT / "local/analysis.json").read_bytes()),
        "pilot_sage_replay_sha256": sha((PILOT / "local/sage_replay.json").read_bytes()),
        "input_law": parent["input_law"],
        "selected_labels": labels,
        "selected_targets": targets,
        "run_order": run_order,
        "planned_run_count": len(run_order),
        "relation_seed": pilot["relation_seed"],
        "relation_workers": pilot["relation_workers"],
        "target_workers": pilot["target_workers"],
        "per_run_timeout_seconds": pilot["per_run_timeout_seconds"],
        "observed_peak_rss_envelope_bytes": pilot["observed_peak_rss_envelope_bytes"],
        "repetitions_per_arm_per_target": 1,
        "one_public_target_per_process": True,
        "measured_interval": pilot["measured_interval"],
        "comparison_variable": pilot["comparison_variable"],
        "development_gate_online_ratio_at_least": "2.0",
        "confirmation_gate_bootstrap_lower_at_least": "2.0",
        "bootstrap_resamples": 20000,
        "bootstrap_seed": 20261010,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
        "smoke_archive_relative": "../koblitz-s3-c1800-n41-20261010/local/smoke_raw.tar.gz",
        "smoke_archive_sha256": pilot["smoke_archive_sha256"],
        "runner_sha256": sha((HERE / "run.py").read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze.py").read_bytes()),
        "preparer_sha256": sha(Path(__file__).read_bytes()),
        "arms": arms,
    }
    freeze_path.write_text(json.dumps(freeze, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"freeze_sha256": sha(freeze_path.read_bytes()),
                      "target_count": len(targets), "planned_runs": len(run_order),
                      "first_target": labels[0], "last_target": labels[-1]}))


if __name__ == "__main__":
    prepare()
