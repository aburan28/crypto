#!/usr/bin/env python3
"""Rebuild the independent C1200/C1800-striped confirmation."""

import gzip
import hashlib
import json
import math
import random
import statistics
import tarfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FREEZE = HERE / "CONFIRM_FREEZE.json"
LEDGER = HERE / "local/raw_paired32.jsonl.gz"
REPORTS = HERE / "local/raw_reports.tar.gz"
HOST = HERE / "local/host.json"


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def geomean(values):
    return math.exp(statistics.mean(math.log(v) for v in values))


def percentile(values, fraction):
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = int(position)
    return ordered[lower] + (ordered[min(lower + 1, len(ordered) - 1)] - ordered[lower]) * (position - lower)


def main():
    freeze = json.loads(FREEZE.read_text())
    host = json.loads(HOST.read_text())
    with gzip.open(LEDGER, "rt") as source:
        rows = [json.loads(line) for line in source]
    assert len(rows) == 32
    assert host["host_exclusive_partition_receipt"] is None
    by_key = {}
    with tarfile.open(REPORTS, "r:gz") as archive:
        archive_names = set(archive.getnames())
        for row in rows:
            key = (row["label"], row["arm"], row["repetition"])
            assert key not in by_key
            by_key[key] = row
            assert row["label"] in freeze["target_labels"]
            assert row["arm"] in freeze["arms"]
            assert row["repetition"] in (1, 2)
            assert row["candidate_id"] == freeze["arms"][row["arm"]]["candidate_id"]
            assert row["binary_sha256"] == freeze["binary_sha256"]
            assert row["source_sha256"] == freeze["source_sha256"]
            assert row["exit_code"] == 0 and row["error"] is None
            assert row["verified_native_and_fixture"] and row["phase_ok"]
            report = row["report"]
            assert report["group_verified"] and report["target_relation_found"]
            assert report["final_rank_complete"]
            assert report["prepared_rank"] == freeze["arms"][row["arm"]]["prepare_rank"]
            assert report["factor_base_digest"] == freeze["arms"][row["arm"]]["factor_base_digest"]
            assert report["factor_base_points"] == freeze["arms"][row["arm"]]["factor_base_points"]
            assert report["target_search_schedule"] == freeze["arms"][row["arm"]]["target_search_schedule"]
            assert report["target_search_prefix_columns"] == freeze["arms"][row["arm"]]["target_search_prefix_columns"]
            assert report["target_search_lead_states"] == freeze["arms"][row["arm"]]["target_search_lead_states"]
            assert report["target_search_stripes"] == freeze["arms"][row["arm"]]["target_search_stripes"]
            assert report["target_count"] == 1 and report["target_query_workers"] == 1
            assert report["timing_ms"]["target_online_after_reusable_setup"] > 0
            assert report["target_s3_calls"] > 0
            name = f"{row['label']}_{row['arm']}_R{row['repetition']}.json"
            assert name in archive_names
            with archive.extractfile(name) as member:
                raw = member.read()
            assert sha256(raw) == row["raw_sha256"]
            assert json.loads(raw) == report
    assert len(by_key) == 32 and len(archive_names) == 32

    target_rows = []
    call_ratios = []
    wall_ratios = []
    aa = {"baseline": [], "expanded": []}
    for label in freeze["target_labels"]:
        baseline = [by_key[label, "baseline", repetition] for repetition in (1, 2)]
        expanded = [by_key[label, "expanded", repetition] for repetition in (1, 2)]
        point_hashes = {row["point_sha256"] for row in baseline + expanded}
        workload_ids = {row["workload_id"] for row in baseline + expanded}
        assert len(point_hashes) == len(workload_ids) == 1
        calls = {}
        walls = {}
        for arm, pair in (("baseline", baseline), ("expanded", expanded)):
            arm_calls = {row["report"]["target_s3_calls"] for row in pair}
            assert len(arm_calls) == 1
            calls[arm] = arm_calls.pop()
            wall_samples = [row["report"]["timing_ms"]["target_online_after_reusable_setup"] for row in pair]
            walls[arm] = geomean(wall_samples)
            aa[arm].append(abs(math.log(wall_samples[0] / wall_samples[1])))
        call_ratio = calls["baseline"] / calls["expanded"]
        wall_ratio = walls["baseline"] / walls["expanded"]
        call_ratios.append(call_ratio)
        wall_ratios.append(wall_ratio)
        target_rows.append({
            "label": label,
            "point_sha256": point_hashes.pop(),
            "workload_id": workload_ids.pop(),
            "baseline_target_s3_calls": calls["baseline"],
            "expanded_target_s3_calls": calls["expanded"],
            "baseline_online_ms_geomean": walls["baseline"],
            "expanded_online_ms_geomean": walls["expanded"],
            "baseline_over_expanded_call_ratio": call_ratio,
            "baseline_over_expanded_online_ratio": wall_ratio,
        })

    rng = random.Random(20261009)
    bootstrap = {"call": [], "online": []}
    for _ in range(20000):
        indices = [rng.randrange(len(target_rows)) for _ in target_rows]
        bootstrap["call"].append(geomean([call_ratios[i] for i in indices]))
        bootstrap["online"].append(geomean([wall_ratios[i] for i in indices]))
    call_geomean = geomean(call_ratios)
    online_geomean = geomean(wall_ratios)
    result = {
        "kind": "s3_n53_c1800stripe_paired_confirmation_analysis_v1",
        "freeze_sha256": sha256(FREEZE.read_bytes()),
        "ledger_decompressed_sha256": sha256(gzip.decompress(LEDGER.read_bytes())),
        "raw_reports_archive_sha256": sha256(REPORTS.read_bytes()),
        "target_count": len(target_rows),
        "run_count": len(rows),
        "native_and_fixture_verified_count": len(rows),
        "root_call_geomean_baseline_over_expanded": call_geomean,
        "root_call_bootstrap95_target_resample": [percentile(bootstrap["call"], 0.025), percentile(bootstrap["call"], 0.975)],
        "root_call_targets_at_least_2x": sum(r >= 2 for r in call_ratios),
        "online_geomean_baseline_over_expanded": online_geomean,
        "online_bootstrap95_target_resample": [percentile(bootstrap["online"], 0.025), percentile(bootstrap["online"], 0.975)],
        "online_targets_at_least_2x": sum(r >= 2 for r in wall_ratios),
        "aa_p95_absolute_log_ratio": {arm: percentile(values, 0.95) for arm, values in aa.items()},
        "confirmation_root_call_geomean_ratio_min": freeze["confirmation_root_call_geomean_ratio_min"],
        "confirmation_root_call_bootstrap95_lower_strictly_gt": freeze["confirmation_root_call_bootstrap95_lower_strictly_gt"],
        "confirmation_passed": (
            call_geomean >= freeze["confirmation_root_call_geomean_ratio_min"]
            and percentile(bootstrap["call"], 0.025)
            > freeze["confirmation_root_call_bootstrap95_lower_strictly_gt"]
            and len(rows) == 32
        ),
        "exploratory_online_geomean_ratio_min": freeze["exploratory_online_geomean_ratio_min"],
        "exploratory_online_bootstrap95_lower_strictly_gt": freeze["exploratory_online_bootstrap95_lower_strictly_gt"],
        "exploratory_online_gate_passed": (
            online_geomean >= freeze["exploratory_online_geomean_ratio_min"]
            and percentile(bootstrap["online"], 0.025)
            > freeze["exploratory_online_bootstrap95_lower_strictly_gt"]
            and len(rows) == 32
        ),
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
        "targets": target_rows,
    }
    (HERE / "analysis.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "targets"}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
