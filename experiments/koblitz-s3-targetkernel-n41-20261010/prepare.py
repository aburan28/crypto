"""Freeze fresh N41 targets for a fixed-base old/new target-kernel test."""

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
MATCHED = ROOT / "experiments/koblitz-s3-c1800-matchedsource-n41-20261010"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(record: dict) -> bytes:
    return json.dumps(record, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def prepare() -> None:
    freeze_path = HERE / "FREEZE.json"
    assert not freeze_path.exists(), "Freeze is immutable"
    parent_path = PARENT / "FREEZE.json"
    parent = json.loads(parent_path.read_bytes())
    pilot = json.loads((PILOT / "FREEZE.json").read_bytes())
    confirmation = json.loads((CONFIRMATION / "FREEZE.json").read_bytes())
    matched = json.loads((MATCHED / "FREEZE.json").read_bytes())
    assert json.loads((MATCHED / "local/analysis.json").read_bytes())["complete"]
    labels = [f"T{number:03}" for number in range(90, 114)]
    targets = [target for target in parent["targets"]
               if target["n"] == 41 and target["label"] in labels]
    assert [target["label"] for target in targets] == labels
    assert len({target["workload_id"] for target in targets}) == 24
    assert not set(labels) & (set(pilot["selected_labels"]) |
                              set(confirmation["selected_labels"]) |
                              set(matched["selected_labels"]))
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
    arms = {
        "baseline": copy.deepcopy(pilot["arms"]["baseline"]),
        "improved": copy.deepcopy(matched["arms"]["baseline"]),
    }
    for arm, original in (("baseline", PILOT), ("improved", MATCHED)):
        frozen_arm = arms[arm]
        candidate_path = original / frozen_arm["candidate_manifest_relative"]
        assert sha(candidate_path.read_bytes()) == frozen_arm["candidate_manifest_sha256"]
        frozen_arm["candidate_manifest_relative"] = "../" + str(
            candidate_path.relative_to(ROOT / "experiments"))
        source = ROOT / "examples" / f"{frozen_arm['binary_name']}.rs"
        assert sha(source.read_bytes()) == frozen_arm["source_sha256"]
        binary = Path("/private/tmp/crypto-s3-release-examples-20261010") / frozen_arm["binary_name"]
        assert sha(binary.read_bytes()) == frozen_arm["binary_sha256"]
    for key in ("columns", "prepare_rank", "factor_base_points", "factor_base_digest", "root_table_slots"):
        assert arms["baseline"][key] == arms["improved"][key]
    assert arms["baseline"]["source_sha256"] != arms["improved"]["source_sha256"]
    run_order = []
    for index, label in enumerate(labels):
        order = ("baseline", "improved") if index % 2 == 0 else ("improved", "baseline")
        run_order.extend({"label": label, "arm": arm, "repetition": 1}
                         for arm in order)
    freeze = {
        "kind": "s3_n41_c1200_targetkernel24_v1",
        "source_tree_git_oid": subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip(),
        "cargo_lock_sha256": sha((ROOT / "Cargo.lock").read_bytes()),
        "parent_freeze_sha256": sha(parent_path.read_bytes()),
        "pilot_freeze_sha256": sha((PILOT / "FREEZE.json").read_bytes()),
        "pilot_analysis_sha256": sha((PILOT / "local/analysis.json").read_bytes()),
        "pilot_sage_replay_sha256": sha((PILOT / "local/sage_replay.json").read_bytes()),
        "prior_confirmation_analysis_sha256": sha((CONFIRMATION / "local/analysis.json").read_bytes()),
        "prior_matchedsource_analysis_sha256": sha((MATCHED / "local/analysis.json").read_bytes()),
        "input_law": parent["input_law"],
        "selected_labels": labels, "selected_targets": targets,
        "run_order": run_order, "planned_run_count": len(run_order),
        "relation_seed": 20260928,
        "relation_workers": 14, "target_workers": 1,
        "per_run_timeout_seconds": 180,
        "observed_peak_rss_envelope_bytes": 32 * 1024**3,
        "repetitions_per_arm_per_target": 1,
        "one_public_target_per_process": True,
        "measured_interval": pilot["measured_interval"],
        "comparison_variable": "target-kernel implementation source; identical curve, base digest, columns, prepared rank, and root-table capacity",
        "development_gate_online_ratio_at_least": "2.0",
        "target_kernel_gate_bootstrap_lower_at_least": "2.0",
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
        "smoke_archive_relative": "../koblitz-s3-c1800-n41-20261010/local/smoke_raw.tar.gz",
        "smoke_archive_sha256": pilot["smoke_archive_sha256"],
        "matched_smoke_archive_relative": "../koblitz-s3-c1800-matchedsource-n41-20261010/local/smoke_raw.tar.zst",
        "matched_smoke_archive_sha256": matched["smoke_archive_sha256"],
        "runner_sha256": sha((HERE / "run.py").read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze.py").read_bytes()),
        "preparer_sha256": sha(Path(__file__).read_bytes()),
        "arms": arms,
    }
    freeze_path.write_text(json.dumps(freeze, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"freeze_sha256": sha(freeze_path.read_bytes()),
                      "selected_targets": [labels[0], labels[-1]],
                      "candidate_ids": {arm: item["candidate_id"] for arm, item in arms.items()},
                      "planned_runs": len(run_order)}))


if __name__ == "__main__":
    prepare()
