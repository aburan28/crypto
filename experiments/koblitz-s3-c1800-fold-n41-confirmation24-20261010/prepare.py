"""Freeze fresh N41 C1200/table versus C1800/PMULL-fold targets."""

from __future__ import annotations

import copy
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PILOT = ROOT / "experiments/koblitz-s3-c1800-n41-20261010"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009"
CONFIRMATION = ROOT / "experiments/koblitz-s3-c1800-n41-confirmation24-20261010"
FOLD_SCREEN = ROOT / "experiments/koblitz-s3-pmull-fold-screen-20261010"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(record: dict) -> bytes:
    return json.dumps(record, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def prepare(binaries: Path, probe: Path) -> None:
    freeze_path = HERE / "FREEZE.json"
    assert not freeze_path.exists(), "Freeze is immutable"
    pilot = json.loads((PILOT / "FREEZE.json").read_bytes())
    pilot_analysis = json.loads((PILOT / "local/analysis.json").read_bytes())
    pilot_replay = json.loads((PILOT / "local/sage_replay.json").read_bytes())
    assert pilot_analysis["complete"] and pilot_analysis["verified_runs"] == 16
    assert pilot_replay["all_scalars_independently_replayed"] is True
    parent = json.loads((PARENT / "FREEZE.json").read_bytes())
    assert platform.system() == "Darwin" and platform.machine() == "arm64"
    labels = [f"T{number:03}" for number in range(114, 138)]
    targets = [target for target in parent["targets"]
               if target["n"] == 41 and target["label"] in labels]
    assert [target["label"] for target in targets] == labels
    assert len({target["workload_id"] for target in targets}) == 24
    assert not (set(labels) & (set(pilot["selected_labels"]) |
                              set(json.loads((CONFIRMATION / "FREEZE.json").read_bytes())["selected_labels"])))
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
        source = ROOT / "examples" / f"{arm['binary_name']}.rs"
        binary = binaries / arm["binary_name"]
        assert sha(source.read_bytes()) == arm["source_sha256"]
        assert sha(binary.read_bytes()) == arm["binary_sha256"]
    assert sha(probe.read_bytes()) == "f0df21eab46ea0c0be6ab3cce70fb6972a92963b72f8dd724fb8e30e2e273a3b"
    probe_expected = {}
    for arm, flag in (("baseline", "0"), ("expanded", "1")):
        env = dict(os.environ)
        env["GF2_ARM_PMULL_FOLD"] = flag
        probe_expected[arm] = subprocess.check_output([str(probe)], env=env, text=True)
    assert probe_expected == {
        "baseline": "n=41 kernel=pmull folds=false\nn=53 kernel=pmull folds=false\n",
        "expanded": "n=41 kernel=pmull folds=true\nn=53 kernel=pmull folds=true\n",
    }
    expanded_manifest = json.loads((HERE / arms["expanded"]["candidate_manifest_relative"]).read_bytes())
    assert expanded_manifest.pop("candidate_id") == arms["expanded"]["candidate_id"]
    flags = expanded_manifest["implementation"]["algorithm_flags"]
    flags["gf2_arm_pmull_fold_env"] = "1"
    flags["gf2_arithmetic_reduction"] = "aarch64_pmull_sparse_tail_fold"
    flags["required_cpu_architecture"] = "aarch64_pmull"
    for stage in ("factor_base", "point_decomposition", "relation_collection",
                  "relation_linear_algebra", "target_descent"):
        expanded_manifest[stage]["gf2_arithmetic_reduction"] = "aarch64_pmull_sparse_tail_fold"
    prefix = arms["expanded"]["candidate_id"].rsplit("h", 1)[0]
    candidate_id = prefix + "h" + sha(canonical(expanded_manifest))[:12]
    expanded_manifest["candidate_id"] = candidate_id
    manifest_path = HERE / "identities/candidates" / f"{candidate_id}.json"
    manifest_path.write_bytes(canonical(expanded_manifest) + b"\n")
    arms["expanded"]["candidate_id"] = candidate_id
    arms["expanded"]["candidate_manifest_relative"] = str(manifest_path.relative_to(HERE))
    arms["expanded"]["candidate_manifest_sha256"] = sha(manifest_path.read_bytes())
    run_order = []
    for index, label in enumerate(labels):
        order = ("baseline", "expanded") if index % 2 == 0 else ("expanded", "baseline")
        run_order.extend({"label": label, "arm": arm, "repetition": 1}
                         for arm in order)
    freeze = {
        "kind": "s3_n41_c1200_c1800_pmullfold_confirmation24_v1",
        "source_tree_git_oid": subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip(),
        "cargo_lock_sha256": sha((ROOT / "Cargo.lock").read_bytes()),
        "parent_freeze_sha256": sha((PARENT / "FREEZE.json").read_bytes()),
        "pilot_freeze_sha256": sha((PILOT / "FREEZE.json").read_bytes()),
        "pilot_analysis_sha256": sha((PILOT / "local/analysis.json").read_bytes()),
        "pilot_sage_replay_sha256": sha((PILOT / "local/sage_replay.json").read_bytes()),
        "prior_confirmation_analysis_sha256": sha((CONFIRMATION / "local/analysis.json").read_bytes()),
        "fold_screen_freeze_sha256": sha((FOLD_SCREEN / "FREEZE.json").read_bytes()),
        "fold_screen_analysis_sha256": sha((FOLD_SCREEN / "local/analysis.json").read_bytes()),
        "probe_source_sha256": sha((FOLD_SCREEN / "probe.rs").read_bytes()),
        "probe_binary_sha256": sha(probe.read_bytes()),
        "probe_expected_stdout": probe_expected,
        "arms_env": {"baseline": "0", "expanded": "1"},
        "host_os": platform.system(), "host_architecture": platform.machine(),
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
        "comparison_variable": "C1200/table versus C1800/ARM64 PMULL sparse-tail fold, with identical curve, public target, input law, resource envelope, and run order pairing",
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
        "replay_sha256": sha((HERE / "replay_sage.py").read_bytes()),
        "preparer_sha256": sha(Path(__file__).read_bytes()),
        "arms": arms,
    }
    freeze_path.write_text(json.dumps(freeze, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"freeze_sha256": sha(freeze_path.read_bytes()),
                      "target_count": len(targets), "planned_runs": len(run_order),
                      "first_target": labels[0], "last_target": labels[-1]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binaries", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.binaries.resolve(), args.probe.resolve())
