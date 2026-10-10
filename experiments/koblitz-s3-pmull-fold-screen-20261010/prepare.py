"""Freeze paired ARM PMULL fold versus table S3 target solves."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE_NAME = "s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan"
PARENTS = {
    41: (ROOT / "experiments/koblitz-s3-targetkernel-n41-20261010", "improved"),
    53: (ROOT / "experiments/koblitz-s3-c1800-matchedsource-n53-20261010", "expanded"),
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: dict) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def main(binaries: Path, probe: Path, equivalence: Path, rlib: Path) -> None:
    if (HERE / "FREEZE.json").exists():
        raise FileExistsError(HERE / "FREEZE.json")
    source_hashes = {BASE_NAME: sha((ROOT / "examples" / f"{BASE_NAME}.rs").read_bytes())}
    binary_hashes = {BASE_NAME: sha((binaries / BASE_NAME).read_bytes())}
    assert source_hashes[BASE_NAME] == "79f77fef1857c9aed10b8d8913ea796e501f9b0afd615345746c453f8d4907b2"
    assert binary_hashes[BASE_NAME] == "d10625fe3aba0e7b61c44f001f119c1bb3b7f4fdde4ee558df1501f1f1ab1164"
    assert platform.system() == "Darwin" and platform.machine() == "arm64"
    assert sha(probe.read_bytes()) == "f0df21eab46ea0c0be6ab3cce70fb6972a92963b72f8dd724fb8e30e2e273a3b"
    assert sha(equivalence.read_bytes()) == "9008d764bf8dc28a2d222df64a32b858d6af6dfdf892a095f360c29dda6f99cf"
    assert sha(rlib.read_bytes()) == "8cb392a206f32e71670c36d5ac31d59e8c8accd1fcb17a8584870b2a235add28"
    equivalence_output = (HERE / "local/equivalence.txt").read_bytes()
    assert equivalence_output == subprocess.check_output([str(equivalence)])
    expected_probe = {}
    for arm, flag in (("table", "0"), ("fold", "1")):
        env = dict(os.environ)
        env["GF2_ARM_PMULL_FOLD"] = flag
        expected_probe[arm] = subprocess.check_output([str(probe)], env=env, text=True)
    assert expected_probe == {
        "table": "n=41 kernel=pmull folds=false\nn=53 kernel=pmull folds=false\n",
        "fold": "n=41 kernel=pmull folds=true\nn=53 kernel=pmull folds=true\n",
    }
    tree_oid = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip()
    assert tree_oid == "0a574800502e67882d8f8f88a846c1e80ee795fd"
    freeze = {
        "kind": "s3_pmull_fold_screen_v1",
        "comparison_variable": "GF2_ARM_PMULL_FOLD=0 versus 1 in one executable on physical ARM64",
        "measured_interval": "one public target, target-dependent query through independently checked scalar",
        "controlled_speedup": None,
        "host_exclusive_partition_receipt": None,
        "big_speedup_gate_online_ratio_at_least": "2.0",
        "per_run_timeout_seconds": 180,
        "observed_peak_rss_envelope_bytes": 34359738368,
        "relation_seed": 20260928,
        "relation_workers": 14,
        "target_workers": 1,
        "host_os": platform.system(), "host_architecture": platform.machine(),
        "source_tree_git_oid": tree_oid,
        "cargo_lock_sha256": sha((ROOT / "Cargo.lock").read_bytes()),
        "probe_source_sha256": sha((HERE / "probe.rs").read_bytes()),
        "probe_binary_sha256": sha(probe.read_bytes()),
        "probe_expected_stdout": expected_probe,
        "equivalence_source_sha256": sha((HERE / "equivalence.rs").read_bytes()),
        "equivalence_binary_sha256": sha(equivalence.read_bytes()),
        "equivalence_output_sha256": sha(equivalence_output),
        "library_rlib_sha256": sha(rlib.read_bytes()),
        "arms_env": {"table": "0", "fold": "1"},
        "preparer_sha256": sha(Path(__file__).read_bytes()),
        "runner_sha256": sha((HERE / "run.py").read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze.py").read_bytes()),
        "replay_sha256": sha((HERE / "replay_sage.py").read_bytes()),
        "source_sha256": source_hashes,
        "binary_sha256": binary_hashes,
        "parent_freeze_sha256": {},
        "curves": {},
        "selected_targets": [],
        "run_order": [],
    }
    identities = HERE / "identities/candidates"
    identities.mkdir(parents=True, exist_ok=True)
    for n, (parent, parent_arm) in PARENTS.items():
        parent_freeze_path = parent / "FREEZE.json"
        parent_freeze = json.loads(parent_freeze_path.read_bytes())
        parent_arm_record = parent_freeze["arms"][parent_arm]
        assert parent_arm_record["binary_name"] == BASE_NAME
        assert parent_arm_record["binary_sha256"] == binary_hashes[BASE_NAME]
        assert parent_arm_record["source_sha256"] == source_hashes[BASE_NAME]
        freeze["parent_freeze_sha256"][str(n)] = sha(parent_freeze_path.read_bytes())
        baseline_path = parent / parent_arm_record["candidate_manifest_relative"]
        baseline = json.loads(baseline_path.read_bytes())
        assert baseline["candidate_id"] == parent_arm_record["candidate_id"]
        candidate = json.loads(baseline_path.read_bytes())
        del candidate["candidate_id"]
        flags = candidate["implementation"]["algorithm_flags"]
        flags["gf2_arm_pmull_fold_env"] = "1"
        flags["gf2_arithmetic_reduction"] = "aarch64_pmull_sparse_tail_fold"
        flags["required_cpu_architecture"] = "aarch64_pmull"
        for stage in ("factor_base", "point_decomposition", "relation_collection",
                      "relation_linear_algebra", "target_descent"):
            candidate[stage]["gf2_arithmetic_reduction"] = "aarch64_pmull_sparse_tail_fold"
        prefix = parent_arm_record["candidate_id"].rsplit("h", 1)[0]
        candidate_id = prefix + "h" + sha(canonical(candidate))[:12]
        candidate["candidate_id"] = candidate_id
        manifest_path = identities / f"{candidate_id}.json"
        manifest_path.write_bytes(canonical(candidate) + b"\n")
        freeze["curves"][str(n)] = {
            "columns": parent_arm_record["columns"],
            "factor_base_points": parent_arm_record["factor_base_points"],
            "factor_base_digest": parent_arm_record["factor_base_digest"],
            "prepare_rank": parent_arm_record["prepare_rank"],
            "root_table_slots": parent_arm_record["root_table_slots"],
            "target_search_schedule": parent_arm_record.get(
                "target_search_schedule", "native_rotating_state_order"),
            "target_search_prefix_columns": parent_arm_record.get(
                "target_search_prefix_columns", parent_arm_record["columns"]),
            "arms": {
                "table": {"binary_name": BASE_NAME,
                           "candidate_id": parent_arm_record["candidate_id"],
                           "candidate_manifest_sha256": sha(baseline_path.read_bytes())},
                "fold": {"binary_name": BASE_NAME, "candidate_id": candidate_id,
                            "candidate_manifest_relative": str(manifest_path.relative_to(HERE)),
                            "candidate_manifest_sha256": sha(manifest_path.read_bytes())},
            },
        }
        for target in parent_freeze["selected_targets"][:8]:
            point_path = (parent.parent / "koblitz-s3-crosscurve-20261009" /
                          target["point_file_relative"] if n == 41 else
                          ROOT / target["point_file_relative"])
            assert sha(point_path.read_bytes()) == target["point_sha256"]
            freeze["selected_targets"].append({
                "n": n, "label": target["label"],
                "point_file_relative": str(point_path.relative_to(ROOT)),
                "point_sha256": target["point_sha256"],
                "public_point": target["public_point"],
                "fixture_scalar": target["fixture_scalar"],
                "workload_id": target["workload_id"],
                "workload_identity_sha256": target["workload_identity_sha256"],
            })
    by_curve = {n: [t for t in freeze["selected_targets"] if t["n"] == n]
                for n in (41, 53)}
    for index in range(8):
        for n in (41, 53):
            target = by_curve[n][index]
            order = ("table", "fold") if (index + n) % 2 else ("fold", "table")
            for arm in order:
                freeze["run_order"].append({"n": n, "label": target["label"],
                                            "arm": arm, "repetition": 1})
    freeze["planned_run_count"] = len(freeze["run_order"])
    assert freeze["planned_run_count"] == 32
    (HERE / "FREEZE.json").write_text(json.dumps(freeze, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"freeze_sha256": sha((HERE / "FREEZE.json").read_bytes()),
                      "planned_runs": 32,
                      "candidate_ids": {n: freeze["curves"][str(n)]["arms"]["fold"]["candidate_id"]
                                        for n in (41, 53)}}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binaries", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--equivalence", type=Path, required=True)
    parser.add_argument("--rlib", type=Path, required=True)
    args = parser.parse_args()
    main(args.binaries.resolve(), args.probe.resolve(), args.equivalence.resolve(), args.rlib.resolve())
