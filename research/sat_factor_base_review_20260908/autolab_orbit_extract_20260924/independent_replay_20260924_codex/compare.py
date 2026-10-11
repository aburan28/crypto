#!/usr/bin/env python3
"""Summarize frozen compact-orbit stage receipts without a DLP speed claim."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
import tarfile

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "all_raw_runs.tar.gz"
GROUPS = {
    "frozen_unsorted": ("frozen_unsorted_run_001", "frozen_unsorted_run_002"),
    "sorted_key_prototype": ("deterministic_run_004", "deterministic_run_005"),
    "scan_order_prototype": ("scan_order_run_001",),
    "deterministic_direct_map": ("direct_map_run_001", "direct_map_run_002"),
}
HOLDOUT = {
    "frozen_unsorted": ("holdout_unsorted_001", "holdout_unsorted_002"),
    "deterministic_direct_map": ("holdout_deterministic_001", "holdout_deterministic_002"),
}
BASE_HASH = "d859319015ea405fd18aee41b51396ce4edcab64ef66265d8edcdeb5e040eb71"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "preserve earlier evidence"
    with tarfile.open(ARCHIVE, "r:gz") as tar:
        def read_json(name):
            return json.load(tar.extractfile(name))
        def replay(name):
            receipt = read_json(f"{name}/validation.json")
            assert receipt["accepted"] and receipt["replayed_targets"] == 12
            assert receipt["independently_checked_archived_relations"] == 12
            assert receipt["natural_targets"] == 8
            assert receipt["base_hash"] == BASE_HASH
            assert receipt["pair_table_entries"] == receipt["edge_selectors"] == 0
            rows = {}
            for kind, count in (("natural", 8), ("planted", 4)):
                for seed in range(1, count + 1):
                    producer = read_json(f"{name}/{kind}_{seed:02d}.stdout.jsonl")
                    extraction = producer["compact_orbit_extraction"]
                    assert extraction["group_valid"] and producer["invalid_group_lifts"] == 0
                    assert producer["decomposition_verdict"] == "SAT"
                    assert producer["pair_table_entries"] == extraction["pair_table_entries"] == 0
                    assert extraction["index_entries"] == 5081560
                    rows[(kind, seed)] = extraction
            return receipt, rows
        def holdout(name):
            report = read_json(f"{name}/results.json")
            assert report["status"] == "PASS" and report["verified_hits"] == 8
            assert report["base_hash"] == BASE_HASH
            assert report["seeds"] == list(range(9, 17))
            rows = {row["seed"]: row for row in report["rows"]}
            assert len(rows) == 8 and all(row["group_valid"] and row["independent_lift"] for row in rows.values())
            return report, rows
        stage = {}
        replay_data = {}
        for label, names in GROUPS.items():
            records = [replay(name) for name in names]
            receipts = [item[0] for item in records]
            assert len({r["producer_source_sha256"] for r in receipts}) == 1
            assert len({r["producer_executable_sha256"] for r in receipts}) == 1
            assert len({r["independent_verifier_sha256"] for r in receipts}) == 1
            per_run = []
            for receipt, rows in records:
                per_run.append({
                    "name": names[len(per_run)],
                    "natural_median_extract_ms": statistics.median(rows["natural", seed]["extract_ms"] for seed in range(1, 9)),
                    "planted_median_extract_ms": statistics.median(rows["planted", seed]["extract_ms"] for seed in range(1, 5)),
                    "natural_median_trials": statistics.median(rows["natural", seed]["trials"] for seed in range(1, 9)),
                    "planted_median_trials": statistics.median(rows["planted", seed]["trials"] for seed in range(1, 5)),
                })
            stage[label] = {
                "source_sha256": receipts[0]["producer_source_sha256"],
                "executable_sha256": receipts[0]["producer_executable_sha256"],
                "runs": per_run,
                "independently_verified_relations": 12 * len(names),
            }
            replay_data[label] = records
        direct = replay_data["deterministic_direct_map"]
        assert all(
            direct[0][1][key]["x_codes"] == direct[1][1][key]["x_codes"]
            and direct[0][1][key]["pinned_intermediates"] == direct[1][1][key]["pinned_intermediates"]
            and direct[0][1][key]["trials"] == direct[1][1][key]["trials"]
            for key in direct[0][1]
        )
        unsorted = replay_data["frozen_unsorted"]
        changed_unsorted_witnesses = sum(
            unsorted[0][1][key]["x_codes"] != unsorted[1][1][key]["x_codes"]
            for key in unsorted[0][1]
        )
        holdout_stage = {}
        holdout_data = {}
        for label, names in HOLDOUT.items():
            records = [holdout(name) for name in names]
            reports = [item[0] for item in records]
            assert len({r["source_sha256"] for r in reports}) == 1
            assert len({r["executable_sha256"] for r in reports}) == 1
            per_run = []
            for name, (_, rows) in zip(names, records):
                per_run.append({
                    "name": name,
                    "median_extract_ms": statistics.median(row["extract_ms"] for row in rows.values()),
                    "median_process_wall_seconds": statistics.median(row["process_wall_seconds"] for row in rows.values()),
                    "median_trials": statistics.median(row["trials"] for row in rows.values()),
                })
            holdout_stage[label] = {
                "source_sha256": reports[0]["source_sha256"],
                "executable_sha256": reports[0]["executable_sha256"],
                "runs": per_run,
                "independently_verified_relations": 8 * len(names),
            }
            holdout_data[label] = records
        candidate_holdout = holdout_data["deterministic_direct_map"]
        assert all(
            candidate_holdout[0][1][seed]["raw_producer"]["compact_orbit_extraction"]["x_codes"]
            == candidate_holdout[1][1][seed]["raw_producer"]["compact_orbit_extraction"]["x_codes"]
            and candidate_holdout[0][1][seed]["trials"] == candidate_holdout[1][1][seed]["trials"]
            for seed in range(9, 17)
        )
        assert stage["deterministic_direct_map"]["source_sha256"] == holdout_stage["deterministic_direct_map"]["source_sha256"]
        assert stage["frozen_unsorted"]["source_sha256"] == holdout_stage["frozen_unsorted"]["source_sha256"]
    result = {
        "schema": "compact-orbit-stage-comparison-v1",
        "status": "PASS",
        "raw_archive_sha256": hashlib.sha256(ARCHIVE.read_bytes()).hexdigest(),
        "base_hash": BASE_HASH,
        "fixed_fixtures": stage,
        "natural_holdout_seeds_9_to_16": holdout_stage,
        "deterministic_candidate_repeated_witnesses_equal": True,
        "unsorted_repeated_witnesses_changed_count": changed_unsorted_witnesses,
        "full_dlp_total_operations": None,
        "rho_ratio": None,
        "complete_negative_coverage": None,
        "claim_boundary": "A reproducible positive compact-orbit extractor and stage-only comparison. No useful-rank collection, final scalar, or matched-rho result; wall timings are exploratory.",
    }
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": result["status"],
        "candidate_repeated_witnesses_equal": result["deterministic_candidate_repeated_witnesses_equal"],
        "unsorted_changed_witnesses": changed_unsorted_witnesses,
        "fixed_natural_extract_medians_ms": [r["natural_median_extract_ms"] for r in stage["deterministic_direct_map"]["runs"]],
        "holdout_hits": holdout_stage["deterministic_direct_map"]["independently_verified_relations"],
    }))


if __name__ == "__main__":
    main()
