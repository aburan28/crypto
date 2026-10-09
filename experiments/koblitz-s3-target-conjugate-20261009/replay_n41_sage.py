"""Independently replay the N41 target-conjugate relations in checked Sage."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

from sage.all import GF, matrix, vector
from sage.env import SAGE_VERSION


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR = ROOT / "experiments/koblitz-s3-prepared-20261009"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
spec = importlib.util.spec_from_file_location("prepared_s3_replay", PRIOR / "replay_local_sage.py")
prior_replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior_replay)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def replay(ledger_path: Path, runtime_path: Path, receipt_path: Path) -> None:
    if receipt_path.exists():
        raise FileExistsError(receipt_path)
    runtime = json.loads(runtime_path.read_text())
    assert runtime["status"] == "verified" and runtime["sage_version"] == SAGE_VERSION
    plan_bytes = (HERE / "N41_FREEZE.json").read_bytes()
    plan = json.loads(plan_bytes)
    parent_bytes = PARENT.read_bytes()
    assert sha(parent_bytes) == plan["parent_plan_sha256"]
    targets = {(item["n"], item["label"]): item
               for item in json.loads(parent_bytes)["targets"]}
    prior = json.loads((PRIOR / "prior_base_replay.json").read_text())
    assert prior["all_factor_base_points_and_labels_checked"] is True
    cfg, base, curve, generator, _, points, labels, point, _ = (
        prior_replay.setup_curve(41, prior))
    assert base["base_hash"] == plan["factor_base_digest"]
    assert len(points) == plan["factor_base_points"]
    r = cfg["r"]
    payload = ledger_path.read_bytes()
    raw = gzip.decompress(payload) if payload[:2] == b"\x1f\x8b" else payload
    rows = [json.loads(line) for line in raw.splitlines()]
    assert len(rows) == 32
    checked = []
    seen_witnesses = set()
    for entry in rows:
        target = targets[(41, entry["label"])]
        frozen = plan["arms"][entry["arm"]]
        assert entry["run_id"] == (
            f"{frozen['candidate_id']}W{target['workload_id']}R{entry['repetition']}")
        assert entry["verified_native_and_fixture"] and entry["phase_ok"]
        report = entry["report"]
        assert report["group_verified"] and report["target_count"] == 1
        assert report["factor_base_digest"] == base["base_hash"]
        assert report["target"] == target["public_point"]
        assert report["recovered_scalar"] == target["fixture_scalar"]
        public_point = point(target["public_point"])
        assert target["fixture_scalar"] * generator == public_point
        indices = tuple(report["target_relation_indices"])
        assert len(indices) == 4
        assert sum((points[index] for index in indices), curve(0)) == public_point
        target_row = prior_replay.row_coefficients(indices, labels, r)
        witnesses = report["rank_relation_witnesses"] + report["fallback_relation_witnesses"]
        coefficients, right_hand_side = [], []
        for witness in witnesses:
            scalar = int(witness["relation_scalar"])
            witness_indices = tuple(witness["point_indices"])
            assert len(witness_indices) == 4
            key = (scalar, witness_indices)
            if key not in seen_witnesses:
                assert sum((points[index] for index in witness_indices),
                           curve(0)) == scalar * generator
                seen_witnesses.add(key)
            coefficients.append(prior_replay.row_coefficients(witness_indices, labels, r))
            right_hand_side.append(scalar)
        field = GF(r)
        relation_matrix = matrix(field, coefficients)
        assert relation_matrix.rank() == report["final_rank"]
        weights = relation_matrix.transpose().solve_right(vector(field, target_row))
        recovered = int(weights.dot_product(vector(field, right_hand_side))) % r
        assert recovered == target["fixture_scalar"]
        assert recovered * generator == public_point
        checked.append({"label": entry["label"], "arm": entry["arm"],
                        "run_id": entry["run_id"], "rank": report["final_rank"],
                        "independent_recovery": True})
    receipt = {"kind": "s3_target_conjugate_n41_independent_sage_replay_v1",
               "sage_version": SAGE_VERSION,
               "runtime_info_sha256": sha(runtime_path.read_bytes()),
               "plan_sha256": sha(plan_bytes),
               "parent_plan_sha256": sha(parent_bytes),
               "ledger_sha256": sha(raw),
               "run_count": len(checked), "target_count": 8,
               "factor_base_points_and_labels_checked": len(points),
               "distinct_relation_witnesses_group_checked": len(seen_witnesses),
               "all_scalars_independently_recovered": True,
               "runs": checked}
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"receipt": str(receipt_path), "runs": len(checked),
                      "witnesses": len(seen_witnesses)}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--runtime-info", type=Path,
                        default=HERE / "local/sage_runtime_info.json")
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    replay(args.ledger, args.runtime_info, args.receipt)
