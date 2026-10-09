"""Write IC1 candidates and one-target workloads for the N53 rank-242 panel."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess

import make_identities as identities


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
prior = identities.prior


def main(ledger: Path, output: Path, plan_path: Path) -> None:
    plan = json.loads(plan_path.read_text())
    rows = [json.loads(line) for line in ledger.read_text().splitlines()]
    assert len(rows) == 81
    assert {(row["label"], row["arm"]) for row in rows} == {
        (target["label"], arm) for target in plan["targets"]
        for arm in ("reference", "shallow", "deep")}
    source_tree = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", f"{identities.SNAPSHOT}:src"],
        text=True).strip()
    assert source_tree == subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD:src"], text=True).strip()
    cargo_lock_sha = prior.sha(ROOT / "Cargo.lock")
    curve = prior.curve_record(53)
    curve_id = curve["curve"]["curve_id"]
    prior.write(output / "curves" / f"{curve_id}.json", curve)
    candidates = {}
    for arm in ("reference", "shallow", "deep"):
        report = next(row["report"] for row in rows if row["arm"] == arm)
        if arm == "reference":
            record = prior.candidate_record(
                53, "candidate", report, curve, identities.SNAPSHOT,
                source_tree, cargo_lock_sha)
        else:
            record = identities.partial_candidate(
                53, report, curve, source_tree, cargo_lock_sha,
                plan["rank_policies"][arm])
        candidate_id = (f"IC1N53Ckb1fb{report['factor_base_points']}PDP4rootRCsample"
                        f"LAgaussTDdirectISO0h{prior.digest(record)}")
        candidates[arm] = candidate_id
        prior.write(output / "candidates" / f"{candidate_id}.json",
                    {**record, "candidate_id": candidate_id})
    workloads = {}
    for target in plan["targets"]:
        seed = (plan["fresh_seed_label"] if target["family"] == "fresh"
                else "s3-partial-prepared-holdout-20261009-v1")
        fixture = {"label": target["label"],
                   "public_point": target["public_point"],
                   "scalar_derivation_sha256": target["scalar_derivation_sha256"]}
        record = identities.workload_record(
            53, fixture, curve_id, seed, plan["relation_seed"])
        workload_id = prior.digest(record)
        workloads[target["label"]] = workload_id
        prior.write(output / "workloads" / f"W{workload_id}.json",
                    {**record, "workload_id": workload_id})
    measurement_rows = []
    for order, row in enumerate(rows, 1):
        report = row["report"]
        candidate_id = candidates[row["arm"]]
        workload_id = workloads[row["label"]]
        measurement_rows.append({
            "candidate_id": candidate_id, "workload_id": workload_id,
            "run_id": f"{candidate_id}W{workload_id}R1",
            "run_order": order, "arm": row["arm"], "n": 53,
            "target_label": row["label"], "target_family": row["family"],
            "target": report["target"] if report else None,
            "status": ("verified_exploratory_mac_contended"
                       if row["verified_native_and_fixture"] else row["error"] or "failure"),
            "strict_host_isolation_receipt": None,
            "verified_native_and_fixture": row["verified_native_and_fixture"],
            "exit_code": row["exit_code"], "error": row["error"],
            "binary_sha256": row["binary_sha256"],
            "target_sha256": row["target_sha256"],
            "raw_sha256": row["raw_sha256"],
            "resource_envelope": {"relation_workers": plan["workers"],
                                  "timeout_seconds": plan["timeout_seconds_per_run"]},
            "prepared_rank": report.get("prepared_rank") if report else None,
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
    measurement_data = b"".join(prior.canonical(row) + b"\n" for row in measurement_rows)
    if measurement_file.exists():
        assert measurement_file.read_bytes() == measurement_data
    else:
        measurement_file.write_bytes(measurement_data)
    print(json.dumps({"curve": curve_id, "candidates": candidates,
                      "workload_count": len(workloads), "run_count": len(rows),
                      "measurements_sha256": prior.sha(measurement_file)}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--plan", type=Path,
                        default=HERE / "rank242/FREEZE.json")
    args = parser.parse_args()
    main(args.ledger, args.output, args.plan)
