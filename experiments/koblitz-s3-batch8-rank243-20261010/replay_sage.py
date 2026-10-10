"""Independently replay rank-243 batch-inversion scalars and relation rows."""

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


def main(ledger: Path, runtime_path: Path, receipt_path: Path) -> None:
    if receipt_path.exists():
        raise FileExistsError(receipt_path)
    runtime = json.loads(runtime_path.read_bytes())
    assert runtime["status"] == "verified" and runtime["sage_version"] == SAGE_VERSION
    freeze_path = HERE / "FREEZE.json"
    freeze = json.loads(freeze_path.read_bytes())
    targets = {(target["n"], target["label"]): target
               for target in freeze["selected_targets"]}
    prior = json.loads((PRIOR_DIR / "prior_base_replay.json").read_bytes())
    assert prior["all_factor_base_points_and_labels_checked"] is True
    raw = ledger.read_bytes()
    rows = [json.loads(line) for line in raw.splitlines()]
    assert len(rows) == freeze["planned_run_count"] and all(
        row["verified_native_and_fixture"] for row in rows)
    distinct_witnesses = set()
    checked = []
    per_target = {}
    for n in (41, 53):
        config, base, curve, generator, _, points, labels, point, _ = (
            prior_replay.setup_curve(n, prior))
        r = config["r"]
        field = GF(r)
        for entry in (row for row in rows if row["n"] == n):
            target = targets[(n, entry["label"])]
            report = entry["report"]
            assert entry["phase_ok"] and report["group_verified"]
            assert report["target_count"] == 1
            assert report["factor_base_digest"] == base["base_hash"]
            assert report["factor_base_points"] == len(points)
            assert report["target"] == target["public_point"]
            assert report["recovered_scalar"] == target["fixture_scalar"]
            assert report["prepared_rank_reached"] and report["prepared_rank"] == 243
            assert entry["candidate_id"] == freeze["candidate_manifests"][
                f"n{n}_{entry['arm']}"]["candidate_id"]
            public_point = point(target["public_point"])
            assert target["fixture_scalar"] * generator == public_point
            indices = report["target_relation_indices"]
            assert sum((points[index] for index in indices), curve(0)) == public_point
            target_row = prior_replay.row_coefficients(indices, labels, r)
            witnesses = (report["rank_relation_witnesses"] +
                         report["fallback_relation_witnesses"])
            matrix_rows = []
            rhs = []
            for witness in witnesses:
                scalar = int(witness["relation_scalar"])
                witness_indices = tuple(witness["point_indices"])
                key = (n, scalar, witness_indices)
                if key not in distinct_witnesses:
                    assert sum((points[index] for index in witness_indices), curve(0)) == (
                        scalar * generator)
                    distinct_witnesses.add(key)
                matrix_rows.append(prior_replay.row_coefficients(witness_indices, labels, r))
                rhs.append(scalar)
            relation_matrix = matrix(field, matrix_rows)
            relation_rhs = vector(field, rhs)
            assert relation_matrix.rank() == report["final_rank"]
            weights = relation_matrix.transpose().solve_right(vector(field, target_row))
            recovered = int(weights.dot_product(relation_rhs)) % r
            assert recovered == target["fixture_scalar"]
            assert recovered * generator == public_point
            per_target.setdefault((n, entry["label"]), {})[
                (entry["arm"], entry["repetition"])] = witnesses
            checked.append({"n": n, "label": entry["label"],
                            "arm": entry["arm"], "repetition": entry["repetition"],
                            "rank": report["final_rank"], "independent_recovery": True})
    assert len(checked) == len(rows)
    for arms in per_target.values():
        assert len(arms) == 4
        witness_lists = list(arms.values())
        assert all(a[:min(len(a), len(b))] == b[:min(len(a), len(b))]
                   for a in witness_lists for b in witness_lists)
    receipt = {
        "kind": "s3_rank243_batch8_independent_sage_replay_v1",
        "sage_version": SAGE_VERSION,
        "runtime_info_sha256": sha(runtime_path.read_bytes()),
        "freeze_sha256": sha(freeze_path.read_bytes()),
        "ledger_sha256": sha(raw),
        "replay_source_sha256": sha(Path(__file__).read_bytes()),
        "run_count": len(rows), "target_count": len(per_target),
        "distinct_relation_witnesses_group_checked": len(distinct_witnesses),
        "all_scalars_independently_recovered": True,
        "runs": checked,
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"receipt": str(receipt_path), "runs": len(rows),
                      "witnesses": len(distinct_witnesses)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--runtime-info", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    main(args.ledger.resolve(), args.runtime_info.resolve(), args.receipt.resolve())
