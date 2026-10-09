"""Independently replay the frozen partial-rank S3 holdout with checked Sage.

Run this only with /Volumes/SSD990/cryptanalysis/sage -python after saving
the launcher's --runtime-info receipt. The prior base replay checks all base
points; this run checks every fresh target and each distinct relation witness.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import time

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


def replay_curve(n: int, rows: list[dict], prior: dict, freeze: dict,
                 first_target: int, target_count_per_curve: int) -> dict:
    cfg, base, curve, generator, representatives, points, labels, point, base_file_sha = (
        prior_replay.setup_curve(n, prior))
    modulus = cfg["r"]
    field = GF(modulus)
    fixtures_path = HERE / f"holdout/n{n}/fixtures.json"
    assert sha(fixtures_path.read_bytes()) == freeze["curves"][str(n)]["fixtures_sha256"]
    fixtures = {item["label"]: item for item in json.loads(fixtures_path.read_text())["fixtures"]}
    assert len(fixtures) == 100 and len(rows) == 2 * target_count_per_curve
    traces = {}
    checked_witnesses = set()
    checked_labels = set()
    checked_targets = []
    by_target = {}

    for entry in rows:
        label, arm = entry["label"], entry["arm"]
        fixture = fixtures[label]
        report = entry["report"]
        assert entry["error"] is None and entry["exit_code"] == 0
        assert entry["verified_native_and_fixture"] is True and entry["phase_ok"] is True
        assert report["group_verified"] is True and report["target_count"] == 1
        assert report["factor_base_digest"] == base["base_hash"]
        assert report["factor_base_points"] == len(points) and report["orbit_columns"] == 244
        assert report["target"] == fixture["public_point"]
        assert report["recovered_scalar"] == int(fixture["fixture_scalar"])
        target = point(fixture["public_point"])
        assert int(fixture["fixture_scalar"]) * generator == target

        target_indices = list(map(int, report["target_relation_indices"]))
        assert len(target_indices) == 4
        assert sum((points[index] for index in target_indices), curve(0)) == target
        checked_labels.update(target_indices)
        target_row = prior_replay.row_coefficients(target_indices, labels, modulus)

        witnesses = report["rank_relation_witnesses"]
        if arm == "partial":
            assert report["prepared_rank_reached"] is True
            assert report["rank"] == freeze["rank_thresholds"][str(n)]
            assert len(witnesses) == report["rank_verified_relations"]
            assert len(report["fallback_relation_witnesses"]) == report["fallback_verified_relations"]
            witnesses = witnesses + report["fallback_relation_witnesses"]
            declared_rank = report["final_rank"]
            assert report["total_new_rows"] == report["rank_new_rows"] + report["fallback_new_rows"]
            assert report["total_verified_relations"] == len(witnesses)
        else:
            assert len(witnesses) == report["rank_verified_relations"]
            declared_rank = report["rank"]
        trace_hash = sha(json.dumps(witnesses, sort_keys=True, separators=(",", ":")).encode())
        if trace_hash not in traces:
            matrix_rows = []
            right_hand_side = []
            for witness in witnesses:
                indices = tuple(map(int, witness["point_indices"]))
                scalar = int(witness["relation_scalar"])
                assert len(indices) == 4 and 0 <= scalar < modulus
                witness_key = (scalar, indices)
                if witness_key not in checked_witnesses:
                    assert sum((points[index] for index in indices), curve(0)) == scalar * generator
                    checked_witnesses.add(witness_key)
                checked_labels.update(indices)
                matrix_rows.append(prior_replay.row_coefficients(indices, labels, modulus))
                right_hand_side.append(scalar)
            relation_matrix = matrix(field, matrix_rows)
            relation_rhs = vector(field, right_hand_side)
            rank = int(relation_matrix.rank())
            assert rank == declared_rank
            traces[trace_hash] = (relation_matrix, relation_rhs, rank)
        relation_matrix, relation_rhs, rank = traces[trace_hash]
        weights = relation_matrix.transpose().solve_right(vector(field, target_row))
        recovered = int(weights.dot_product(relation_rhs)) % modulus
        assert recovered == int(fixture["fixture_scalar"])
        assert recovered * generator == target
        if report["solved_base_logs"] is not None:
            assert rank == 244
            logs = relation_matrix.solve_right(relation_rhs)
            assert [int(value) for value in logs] == report["solved_base_logs"]
            assert int(vector(field, target_row).dot_product(logs)) % modulus == recovered

        by_target.setdefault(label, {})[arm] = witnesses
        checked_targets.append({"label": label, "arm": arm,
                                "raw_sha256": entry["raw_sha256"],
                                "trace_sha256": trace_hash,
                                "rank": rank, "recovered_scalar": recovered})
        print(json.dumps({"n": n, "target": label, "arm": arm,
                          "rank": rank, "independent_recovery": True}), flush=True)

    assert set(by_target) == {f"T{index:03}" for index in range(
        first_target, first_target + target_count_per_curve)}
    for label, arms in by_target.items():
        assert set(arms) == {"reference", "partial"}
        ref, part = arms["reference"], arms["partial"]
        assert ref[:min(len(ref), len(part))] == part[:min(len(ref), len(part))], label
    for index in checked_labels:
        column, coefficient = labels[index]
        assert points[index] == coefficient * representatives[column]
    assert len({(item["label"], item["arm"]) for item in checked_targets}) == 2 * target_count_per_curve
    return {"n": n, "field_modulus_powers": cfg["modulus_powers"],
            "subgroup_order": modulus, "generator": cfg["generator"],
            "base_fixture_sha256": base_file_sha, "base_digest": base["base_hash"],
            "factor_base_points_prior_checked": len(points),
            "point_labels_rechecked": len(checked_labels),
            "distinct_relation_traces": len(traces),
            "distinct_relation_witnesses_group_checked": len(checked_witnesses),
            "all_common_relation_prefixes_equal": True,
            "targets": checked_targets}


def main(ledger_path: Path, runtime_path: Path, output_path: Path,
         first_target: int, target_count_per_curve: int) -> None:
    if output_path.exists():
        raise FileExistsError(output_path)
    runtime = json.loads(runtime_path.read_text())
    assert runtime["status"] == "verified" and runtime["sage_version"] == SAGE_VERSION
    prior_path = PRIOR_DIR / "prior_base_replay.json"
    prior = json.loads(prior_path.read_text())
    assert prior["all_factor_base_points_and_labels_checked"] is True
    freeze_path = HERE / "FREEZE.json"
    freeze = json.loads(freeze_path.read_text())
    started = time.perf_counter_ns()
    ledger_bytes = ledger_path.read_bytes()
    rows = [json.loads(line) for line in ledger_bytes.splitlines()]
    assert len(rows) == 4 * target_count_per_curve
    for entry in rows:
        raw_path = ledger_path.parent / f"n{entry['n']}_{entry['label']}_{entry['arm']}.jsonl"
        assert sha(raw_path.read_bytes()) == entry["raw_sha256"]
        assert json.loads(raw_path.read_text()) == entry["report"]
        point_path = HERE / f"holdout/n{entry['n']}/{entry['label']}/public_target.json"
        assert sha(point_path.read_bytes()) == entry["target_sha256"]
    curves = [replay_curve(n, [row for row in rows if row["n"] == n], prior, freeze,
                           first_target, target_count_per_curve)
              for n in (41, 53)]
    result = {"kind": "s3_partial_prepared_fresh_holdout_independent_sage_replay_v1",
              "sage_version": SAGE_VERSION,
              "sage_runtime_info_sha256": sha(runtime_path.read_bytes()),
              "replay_script_sha256": sha(Path(__file__).read_bytes()),
              "prior_base_replay_sha256": sha(prior_path.read_bytes()),
              "freeze_sha256": sha(freeze_path.read_bytes()),
              "ledger_sha256": sha(ledger_bytes), "run_count": len(rows),
              "target_count": 2 * target_count_per_curve,
              "first_target_index": first_target,
              "target_count_per_curve": target_count_per_curve,
              "all_scalars_independently_recovered": True,
              "distinct_relation_witnesses_group_checked": sum(
                  curve["distinct_relation_witnesses_group_checked"] for curve in curves),
              "curves": curves, "replay_wall_ns": time.perf_counter_ns() - started}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"receipt": str(output_path), "runs": len(rows),
                      "witnesses": result["distinct_relation_witnesses_group_checked"]}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", required=True, type=Path)
    parser.add_argument("--runtime-info", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    parser.add_argument("--first-target", type=int, default=1)
    parser.add_argument("--target-count-per-curve", type=int, default=24)
    args = parser.parse_args()
    main(args.ledger, args.runtime_info, args.receipt,
         args.first_target, args.target_count_per_curve)
