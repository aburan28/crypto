"""Freeze the source-identical T025-T033 N41/N53 extension before timing."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OLD = ROOT / "experiments/koblitz-s3-crosscurve-20261009"
FREEZE = HERE / "FREEZE.json"
SOURCE_COMMIT = "b562cdaed14254d4f6c5f2a4640b302592e8882b"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def make(binaries: Path) -> dict:
    old_raw = (OLD / "FREEZE.json").read_bytes()
    plan = json.loads(old_raw)
    assert plan["kind"] == "s3_rank243_matched_curve_fresh_panel_v1"
    assert plan["initial_pilot_count_per_curve"] == 24
    assert plan["workers"] == 14 and plan["timeout_seconds_per_run"] == 180
    old_analysis_raw = (OLD / "local/analysis_pilot24.json").read_bytes()
    old_analysis = json.loads(old_analysis_raw)
    assert old_analysis["matched_index_count"] == 24
    assert old_analysis["run_count"] == 192
    assert all(old_analysis["curves"][str(n)]["verified_target_count"] == 24
               for n in (41, 53))
    expected_binaries = {
        "reference": ("s3_pair_candidate_frozen", plan["reference_binary_sha256_local"]),
        "prepared": ("s3_partial_prepared_targetspan", plan["partial_binary_sha256_local"]),
    }
    for arm, (name, digest) in expected_binaries.items():
        assert (binaries / name).is_file()
        assert sha((binaries / name).read_bytes()) == digest
        assert old_analysis["binary_sha256"][arm] == digest
    assert sha((ROOT / "experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs").read_bytes()) == plan[
        "reference_source_sha256"]
    assert sha((ROOT / "examples/s3_partial_prepared_targetspan.rs").read_bytes()) == plan[
        "partial_source_sha256"]
    assert subprocess.check_output(["git", "rev-parse", SOURCE_COMMIT], cwd=ROOT,
                                   text=True).strip() == SOURCE_COMMIT
    tree = subprocess.check_output(["git", "rev-parse", SOURCE_COMMIT + "^{tree}"],
                                   cwd=ROOT, text=True).strip()
    spec = importlib.util.spec_from_file_location("crosscurve_runner", OLD / "run_local_pilot.py")
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    selected = runner.selected_targets(plan, 25, 9)
    assert len(selected) == 18
    assert {int(row["label"][1:]) for row in selected} == set(range(25, 34))
    order = []
    for target in selected:
        repeats = Counter()
        for arm in runner.block_order(target["n"], target["label"]):
            repeats[arm] += 1
            candidate = plan["curves"][str(target["n"])][
                "reference_candidate_id" if arm == "reference" else "rank243_candidate_id"]
            order.append({
                "n": target["n"], "label": target["label"], "arm": arm,
                "repetition": repeats[arm],
                "run_id": f"{candidate}W{target['workload_id']}R{repeats[arm]}",
            })
    assert len(order) == 72 and len({row["run_id"] for row in order}) == 72
    return {
        "kind": "s3_crosscurve_extension33_freeze_v1",
        "old_freeze_sha256": sha(old_raw),
        "old_raw_archive_sha256": sha((OLD / "local/raw_pilot24.tar.gz").read_bytes()),
        "old_analysis_sha256": sha(old_analysis_raw),
        "old_runner_sha256": sha((OLD / "run_local_pilot.py").read_bytes()),
        "old_analyzer_sha256": sha((OLD / "analyze_local_pilot.py").read_bytes()),
        "old_sage_replay_sha256": sha((OLD / "replay_sage.py").read_bytes()),
        "source_commit": SOURCE_COMMIT,
        "source_tree_oid": tree,
        "rebuild_command": "CARGO_TARGET_DIR=/private/tmp/s3-crosscurve-b562-build cargo build --release --offline --locked --example s3_pair_candidate_frozen --example s3_partial_prepared_targetspan",
        "cargo_version": subprocess.check_output(["cargo", "--version"], text=True).strip(),
        "rustc_version": subprocess.check_output(["rustc", "--version"], text=True).strip(),
        "binary_sha256": {arm: digest for arm, (_, digest) in expected_binaries.items()},
        "binary_bytes": {arm: (binaries / name).stat().st_size
                         for arm, (name, _) in expected_binaries.items()},
        "selected_indices": list(range(25, 34)),
        "selected_targets": [{key: row[key] for key in
                              ("n", "label", "public_point", "point_sha256", "workload_id")}
                             for row in selected],
        "run_order": order,
        "planned_run_count": 72,
        "workers": 14,
        "per_run_timeout_seconds": 180,
        "comparison": "rank243_prepared_versus_target_dependent_reference_on_same_point",
        "operation_count_boundary": "reference_rank_target_s3_calls_plus_target_s3_calls_over_prepared_fallback_s3_calls_plus_target_s3_calls",
        "bootstrap": {"draws": 20000, "seed": 20261010,
                      "unit": "matched_target_index", "interval": "two_sided_95_percent"},
        "prepare_source_sha256": sha(Path(__file__).read_bytes()),
        "runner_sha256": sha((HERE / "run.py").read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze.py").read_bytes()),
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binaries", required=True, type=Path)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    plan = make(args.binaries.resolve())
    raw = json.dumps(plan, indent=2, sort_keys=True).encode() + b"\n"
    if args.verify:
        assert FREEZE.read_bytes() == raw
    else:
        if FREEZE.exists():
            raise FileExistsError(FREEZE)
        FREEZE.write_bytes(raw)
    print(json.dumps({"freeze_sha256": sha(raw),
                      "run_count": plan["planned_run_count"],
                      "source_commit": plan["source_commit"],
                      "binary_sha256": plan["binary_sha256"]}, sort_keys=True))
