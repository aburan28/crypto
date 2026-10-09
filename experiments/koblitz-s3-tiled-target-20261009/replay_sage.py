"""Independently replay N53 tiled target relations and logs in Sage."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

from sage.all import GF, matrix, vector
from sage.env import SAGE_VERSION


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR_DIR = ROOT / "experiments/koblitz-s3-prepared-20261009"
spec = importlib.util.spec_from_file_location("prior_s3_replay", PRIOR_DIR / "replay_local_sage.py")
prior_replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior_replay)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main(ledger_path: Path, runtime_path: Path, receipt_path: Path) -> None:
    if receipt_path.exists():
        raise FileExistsError(receipt_path)
    runtime = json.loads(runtime_path.read_text())
    assert runtime["status"] == "verified" and runtime["sage_version"] == SAGE_VERSION
    plan_path = HERE / "FREEZE.json"
    plan = json.loads(plan_path.read_text())
    parent_path = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
    parent_bytes = parent_path.read_bytes()
    assert plan["parent_plan_sha256"] == sha(parent_bytes)
    parent = json.loads(parent_bytes)
    targets = {(item["n"], item["label"]): item for item in parent["targets"]}
    prior_path = PRIOR_DIR / "prior_base_replay.json"
    prior = json.loads(prior_path.read_text())
    assert prior["all_factor_base_points_and_labels_checked"] is True
    ledger_bytes = ledger_path.read_bytes()
    rows = [json.loads(line) for line in ledger_bytes.splitlines()]
    assert len(rows) == 8 * plan["target_count"]
    checked = []
    distinct_witnesses = set()
    per_target = {}
    for n in (53,):
        cfg, base, curve, generator, representatives, points, labels, point, _ = (
            prior_replay.setup_curve(n, prior))
        r = cfg["r"]
        field = GF(r)
        for entry in (row for row in rows if row["n"] == n):
            target = targets[(n, entry["label"])]
            report = entry["report"]
            assert entry["verified_native_and_fixture"] and entry["phase_ok"]
            assert report["group_verified"] and report["target_count"] == 1
            assert report["factor_base_digest"] == base["base_hash"]
            assert report["factor_base_points"] == len(points)
            assert report["target"] == target["public_point"]
            assert report["recovered_scalar"] == target["fixture_scalar"]
            public_point = point(target["public_point"])
            assert target["fixture_scalar"] * generator == public_point
            indices = report["target_relation_indices"]
            assert sum((points[index] for index in indices), curve(0)) == public_point
            target_row = prior_replay.row_coefficients(indices, labels, r)
            witnesses = report["rank_relation_witnesses"]
            assert entry["arm"] in ("parallel", "tile_2048", "tile_4096",
                                    "tile_8192")
            candidate = (plan["reference_candidate_id"]
                         if entry["tile_states"] is None else
                         plan["tiled_candidate_ids_by_tile_states"][str(entry["tile_states"])])
            assert entry["candidate_id"] == candidate
            assert entry["run_id"] == (
                f"{candidate}W{target['workload_id']}R{entry['repetition']}")
            assert report["prepared_rank_reached"] and report["prepared_rank"] == 243
            assert report["group_lift_valid_states"] == report["regular_states"]
            assert report["group_lift_fallback_states"] == 0
            witnesses = witnesses + report["fallback_relation_witnesses"]
            declared_rank = report["final_rank"]
            rows_matrix = []
            rhs = []
            for witness in witnesses:
                scalar = int(witness["relation_scalar"])
                witness_indices = tuple(witness["point_indices"])
                key = (n, scalar, witness_indices)
                if key not in distinct_witnesses:
                    assert sum((points[index] for index in witness_indices), curve(0)) == scalar * generator
                    distinct_witnesses.add(key)
                rows_matrix.append(prior_replay.row_coefficients(witness_indices, labels, r))
                rhs.append(scalar)
            relation_matrix = matrix(field, rows_matrix)
            relation_rhs = vector(field, rhs)
            assert relation_matrix.rank() == declared_rank
            weights = relation_matrix.transpose().solve_right(vector(field, target_row))
            recovered = int(weights.dot_product(relation_rhs)) % r
            assert recovered == target["fixture_scalar"] and recovered * generator == public_point
            per_target.setdefault((n, entry["label"]), {})[
                (entry["arm"], entry["run_id"])] = [
                    (int(w["relation_scalar"]), tuple(w["point_indices"]))
                    for w in witnesses]
            checked.append({"n": n, "label": entry["label"], "arm": entry["arm"],
                            "run_id": entry["run_id"], "rank": declared_rank,
                            "independent_recovery": True})
    assert len(checked) == len(rows)
    for (n, label), arms in per_target.items():
        assert len(arms) == 8
        witnesses = list(arms.values())
        assert all(a[:min(len(a), len(b))] == b[:min(len(a), len(b))]
                   for a in witnesses for b in witnesses)
    receipt = {"kind": "s3_tiled_size_independent_sage_replay_v1",
               "sage_version": SAGE_VERSION,
               "sage_runtime_info_sha256": sha(runtime_path.read_bytes()),
               "plan_sha256": sha(plan_path.read_bytes()),
               "parent_plan_sha256": sha(parent_bytes),
               "ledger_sha256": sha(ledger_bytes),
               "run_count": len(rows), "target_count": len(per_target),
               "distinct_relation_witnesses_group_checked": len(distinct_witnesses),
               "all_scalars_independently_recovered": True,
               "runs": checked}
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"receipt": str(receipt_path), "runs": len(rows),
                      "witnesses": len(distinct_witnesses)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--runtime-info", type=Path, default=HERE / "sage_runtime_info.json")
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    main(args.ledger, args.runtime_info, args.receipt)
