"""Write EC1/IC1/W/R manifests for the frozen partial-rank holdout panel."""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR = ROOT / "experiments/koblitz-s3-prepared-20261009/make_identities.py"
spec = importlib.util.spec_from_file_location("prior_s3_identities", PRIOR)
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)
PARTIAL_SOURCE = ROOT / "examples/s3_partial_prepared_targetspan.rs"
SNAPSHOT = "833ef86191947ba8e4def56a78ee6c51a4ea0ff2"


def partial_candidate(n: int, report: dict, curve: dict,
                      source_tree: str, cargo_lock_sha: str, rank: int) -> dict:
    record = copy.deepcopy(prior.candidate_record(
        n, "candidate", report, curve, SNAPSHOT, source_tree, cargo_lock_sha))
    source_sha = prior.sha(PARTIAL_SOURCE)
    for section in ("factor_base", "relation_collection", "relation_linear_algebra",
                    "target_descent"):
        record[section]["source_sha256"] = source_sha
    record["point_decomposition"]["implementation_source_sha256"] = source_sha
    record["point_decomposition"]["limits"] = {
        "max_total_relation_attempts": "max(8*columns,columns+32)"}
    record["relation_collection"]["stop_criterion"] = (
        "fixed_target_independent_rank_before_target_input_then_"
        "online_fallback_until_target_coefficient_row_in_span")
    record["relation_collection"]["pre_target_rank_threshold"] = rank
    record["relation_collection"]["fallback_charging"] = "all_fallback_queries_and_attempts_online"
    record["relation_linear_algebra"]["rank_criterion"] = (
        "fixed_pre_target_rank_then_target_row_span_or_full_244_fallback")
    record["implementation"]["example_source_sha256"] = source_sha
    record["implementation"]["algorithm_flags"]["prepare_rank"] = rank
    record["implementation"]["algorithm_flags"]["online_fallback"] = True
    return record


def workload_record(n: int, fixture: dict, curve_id: str, seed_label: str,
                    relation_seed: int) -> dict:
    return {"curve_id": curve_id, "subgroup_order_r": prior.CURVES[n]["r"],
            "target_count": 1, "target_point": fixture["public_point"],
            "target_label": fixture["label"],
            "target_input_law": "known_scalar_sha256_public_point_fixture_v1",
            "target_seed_label": seed_label,
            "target_index": int(fixture["label"][1:]),
            "scalar_derivation_sha256": fixture["scalar_derivation_sha256"],
            "relation_seed": relation_seed,
            "preparation_state_at_online_start": "candidate_policy_reusable_setup_ready"}


def main(ledger: Path, output: Path) -> None:
    entries = [json.loads(line) for line in ledger.read_text().splitlines()]
    assert len(entries) in (96, 304)
    freeze = json.loads((HERE / "FREEZE.json").read_text())
    if len(entries) == 304:
        extension = json.loads((HERE / "EXTENSION_FREEZE.json").read_text())
        assert prior.sha(HERE / "FREEZE.json") == extension["parent_freeze_sha256"]
        selected_labels = {f"T{index:03}" for index in range(25, 101)}
    else:
        selected_labels = {f"T{index:03}" for index in range(1, 25)}
    assert {(row["n"], row["label"], row["arm"]) for row in entries} == {
        (n, label, arm) for n in (41, 53)
        for label in selected_labels for arm in ("reference", "partial")}
    assert prior.sha(PARTIAL_SOURCE) == freeze["partial_source_sha256"]
    source_tree = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", f"{SNAPSHOT}:src"], text=True).strip()
    current_tree = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip()
    assert source_tree == current_tree
    cargo_lock_sha = prior.sha(ROOT / "Cargo.lock")
    curves = {n: prior.curve_record(n) for n in (41, 53)}
    candidates = {}
    workloads = {}
    for n in (41, 53):
        curve_id = curves[n]["curve"]["curve_id"]
        prior.write(output / "curves" / f"{curve_id}.json", curves[n])
        reference_report = next(row["report"] for row in entries
                                if row["n"] == n and row["arm"] == "reference")
        partial_report = next(row["report"] for row in entries
                              if row["n"] == n and row["arm"] == "partial")
        for arm, report in (("reference", reference_report), ("partial", partial_report)):
            if arm == "reference":
                record = prior.candidate_record(
                    n, "candidate", report, curves[n], SNAPSHOT, source_tree, cargo_lock_sha)
            else:
                record = partial_candidate(
                    n, report, curves[n], source_tree, cargo_lock_sha,
                    freeze["rank_thresholds"][str(n)])
            candidate_id = (f"IC1N{n}Ckb1fb{report['factor_base_points']}PDP4rootRCsample"
                            f"LAgaussTDdirectISO0h{prior.digest(record)}")
            candidates[n, arm] = candidate_id
            prior.write(output / "candidates" / f"{candidate_id}.json",
                        {**record, "candidate_id": candidate_id})
        fixtures_path = HERE / f"holdout/n{n}/fixtures.json"
        fixtures = json.loads(fixtures_path.read_text())
        assert prior.sha(fixtures_path) == freeze["curves"][str(n)]["fixtures_sha256"]
        for fixture in fixtures["fixtures"]:
            if fixture["label"] not in selected_labels:
                continue
            record = workload_record(n, fixture, curve_id, fixtures["seed_label"],
                                     freeze["relation_seed"])
            workload_id = prior.digest(record)
            workloads[n, fixture["label"]] = workload_id
            prior.write(output / "workloads" / f"W{workload_id}.json",
                        {**record, "workload_id": workload_id})
    run_rows = []
    for order, row in enumerate(entries, 1):
        report = row["report"]
        candidate_id = candidates[row["n"], row["arm"]]
        workload_id = workloads[row["n"], row["label"]]
        run_rows.append({
            "candidate_id": candidate_id, "workload_id": workload_id,
            "run_id": f"{candidate_id}W{workload_id}R1",
            "run_order": order, "arm": row["arm"], "n": row["n"],
            "target_label": row["label"],
            "target": report["target"] if report else None,
            "status": ("verified_exploratory_mac_contended"
                       if row["verified_native_and_fixture"] else row["error"] or "failure"),
            "strict_host_isolation_receipt": None,
            "verified_native_and_fixture": row["verified_native_and_fixture"],
            "exit_code": row["exit_code"], "error": row["error"],
            "binary_sha256": row["binary_sha256"],
            "target_sha256": row["target_sha256"],
            "raw_sha256": row["raw_sha256"],
            "resource_envelope": {"relation_workers": freeze["workers"],
                                  "timeout_seconds": freeze["timeout_seconds_per_run"]},
            "prepared_rank": (report["prepared_rank"] if row["arm"] == "partial" else None)
                             if report else None,
            "rank_attempts": report["rank_attempts"] if report else None,
            "fallback_attempts": report.get("fallback_attempts") if report else None,
            "rank_failures": report["rank_failures"] if report else None,
            "fallback_failures": report.get("fallback_failures") if report else None,
            "verified_relations": (report.get("total_verified_relations",
                                               report["rank_verified_relations"]) if report else None),
            "novel_rows": (report.get("total_new_rows", report["rank_new_rows"])
                           if report else None),
            "final_rank": (report.get("final_rank", report["rank"]) if report else None),
            "actual_factor_base_points": report["factor_base_points"] if report else None,
            "folded_columns": report["orbit_columns"] if report else None,
            "target_relation_found": report["target_relation_found"] if report else None,
            "recovered_scalar": report["recovered_scalar"] if report else None,
            "peak_rss_bytes": report["peak_rss_bytes"] if report else None,
            "timing_ms": report["timing_ms"] if report else None,
        })
    measurement_file = output / "measurements.jsonl"
    measurement_data = b"".join(prior.canonical(row) + b"\n" for row in run_rows)
    if measurement_file.exists():
        assert measurement_file.read_bytes() == measurement_data
    else:
        measurement_file.write_bytes(measurement_data)
    print(json.dumps({"curves": [curves[n]["curve"]["curve_id"] for n in (41, 53)],
                      "candidates": list(candidates.values()),
                      "workload_count": len(workloads), "run_count": len(run_rows),
                      "measurements_sha256": prior.sha(measurement_file)}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    main(args.ledger, args.output)
