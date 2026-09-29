#!/usr/bin/env python3
"""Timed base-log solving and scalar-blind public-point recovery for native L384.

This is the work needed to turn the producer's x-coordinate witnesses into
actual logarithms. The independent audit uses a different field/group law.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SHARED = HERE.parent / "autolab_shared_log_n53_20260925"
BASE_HASH = "2de8ec46916999dfe5e1ac95e68ea15a3f001eec27e013fa96d9525719b4f684"
POINTS_HASH = "0de93f652949981367385818685a255c98c4d230a0319fc658b7e91836d36cce"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


math = load(SHARED / "generate_targets.py", "operational_shift_reduce_group")


class SparseElim:
    """Operational modular row reduction; independent audit uses cold_batch_rank.Echelon."""

    def __init__(self, columns: int, modulus: int):
        self.columns = columns
        self.modulus = modulus
        self.pivots = {}

    def insert(self, coefficients: dict[int, int], rhs: int) -> bool:
        modulus = self.modulus
        row = {column: value % modulus for column, value in coefficients.items() if value % modulus}
        rhs %= modulus
        while row:
            pivot = min(row)
            if pivot not in self.pivots:
                scale = pow(row[pivot], -1, modulus)
                self.pivots[pivot] = ({column: value * scale % modulus for column, value in row.items()}, rhs * scale % modulus)
                return True
            old, old_rhs = self.pivots[pivot]
            factor = row[pivot]
            for column, value in old.items():
                updated = (row.get(column, 0) - factor * value) % modulus
                if updated:
                    row[column] = updated
                else:
                    row.pop(column, None)
            rhs = (rhs - factor * old_rhs) % modulus
        assert rhs == 0, "inconsistent relation"
        return False

    def solution(self):
        if len(self.pivots) != self.columns:
            return None
        result = [0] * self.columns
        for pivot in sorted(self.pivots, reverse=True):
            row, rhs = self.pivots[pivot]
            result[pivot] = (rhs - sum(value * result[column] for column, value in row.items() if column > pivot)) % self.modulus
        return result


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def base_maps(header: dict):
    assert header["base_hash"] == BASE_HASH
    assert header["n"] == 53 and header["a"] == 0
    assert header["subgroup_order"] == math.ORDER
    assert tuple(header["generator"]) == math.GENERATOR
    points = [tuple(p) for p in header["factor_base_point_coordinates"]]
    labels = [tuple(label) for label in header["factor_base_point_labels"]]
    reps = [tuple(p) for p in header["factor_base_representatives"]]
    assert len(points) == len(labels) == 23320
    assert len(reps) == 220 and len(set(points)) == len(points)
    by_point = dict(zip(points, labels))
    by_x = defaultdict(list)
    for point in points:
        by_x[point[0]].append(point)
    assert len(by_x) == 11660 and all(len(pair) == 2 for pair in by_x.values())
    assert all(0 <= col < 220 and 0 < coefficient < math.ORDER for col, coefficient in labels)
    return by_point, by_x, reps


def recover_witness(by_x, target, codes, intermediates):
    """Resolve four factor-base signs using the pinned pair x values and Q."""
    assert len(codes) == 4 and len(intermediates) == 2
    assert all(code in by_x for code in codes)
    pair_choices = []
    for offset, pinned in ((0, intermediates[0]), (2, intermediates[1])):
        choices = []
        for left, right in itertools.product(by_x[codes[offset]], by_x[codes[offset + 1]]):
            total = math.add(left, right)
            if total is not None and total[0] == pinned:
                choices.append((left, right, total))
        assert choices, "no pinned pair lift"
        pair_choices.append(choices)
    for left, right in itertools.product(*pair_choices):
        if math.add(left[2], right[2]) == target:
            return left[:2] + right[:2]
    raise AssertionError("no four-point lift matches public target")


def training(args):
    header_file = args.training / "base_header.jsonl"
    raw_file = args.training / "producer.stdout.jsonl"
    schedule_file = args.training / "target_scalars.txt"
    header = json.loads(header_file.read_text())
    by_point, by_x, reps = base_maps(header)
    schedule = [int(line) for line in schedule_file.read_text().splitlines()]
    assert len(schedule) == len(set(schedule)) == 512
    observed = json.loads(raw_file.read_text())
    assert observed["factor_base_input_hash"] is None
    assert observed["factor_base_input_path"] is None
    batch = observed["compact_orbit_batch"]
    assert batch["targets_requested"] == 512
    assert batch["regular_scan_policy"] == "target_cyclic_v1"
    relations, failures, queries = (
        batch["relations"], batch["failed_target_scalars"], batch["query_observations"]
    )
    assert len(queries) == 512 and len(relations) + len(failures) == 512
    matrix = SparseElim(220, math.ORDER)
    first_full_rank = None
    gains = []
    post_rank_predictions = 0
    solution = None
    relation_index = failure_index = 0
    for index, (scalar, query) in enumerate(zip(schedule, queries), start=1):
        assert query["scalar"] == scalar
        if not query["hit"]:
            assert failures[failure_index] == scalar
            failure_index += 1
            continue
        relation = relations[relation_index]
        relation_index += 1
        assert relation["scalar"] == scalar
        target = math.scalar(math.GENERATOR, scalar)
        assert target is not None
        lift = recover_witness(by_x, target, relation["x_codes"], relation["pinned_intermediates"])
        row = {}
        for point in lift:
            column, coefficient = by_point[point]
            row[column] = (row.get(column, 0) + coefficient) % math.ORDER
        if solution is not None:
            assert sum(value * solution[col] for col, value in row.items()) % math.ORDER == scalar
            post_rank_predictions += 1
        if matrix.insert(row, scalar):
            gains.append(index)
            if len(matrix.pivots) == 220 and first_full_rank is None:
                first_full_rank = index
                solution = matrix.solution()
    assert relation_index == len(relations) and failure_index == len(failures)
    assert solution == matrix.solution()
    if solution is not None:
        for rep, log in zip(reps, solution):
            assert math.scalar(math.GENERATOR, log) == rep
    report = {
        "classification": "OPERATIONAL_TRAINING_LOGS_RECOVERED" if solution is not None else "OPERATIONAL_RANK_DEFICIENT",
        "training_raw_sha256": sha(raw_file),
        "base_header_sha256": sha(header_file),
        "schedule_sha256": sha(schedule_file),
        "rank": len(matrix.pivots), "columns": 220,
        "relations_used": len(relations), "failed_target_scalars": failures,
        "first_full_rank_at": first_full_rank,
        "rank_gain_indices": gains, "post_rank_predictions": post_rank_predictions,
        "base_log_solution": solution,
    }
    save(args.out, report)


def point(args):
    assert sha(args.points) == POINTS_HASH
    points = [tuple(json.loads(line)) for line in args.points.read_text().splitlines()]
    assert len(points) == len(set(points)) == 384
    header_file = args.training / "base_header.jsonl"
    by_point, by_x, _ = base_maps(json.loads(header_file.read_text()))
    logs_file = args.training / "operational_solution.json"
    logs_report = json.loads(logs_file.read_text())
    assert logs_report["classification"] == "OPERATIONAL_TRAINING_LOGS_RECOVERED"
    logs = logs_report["base_log_solution"]
    assert len(logs) == 220
    raw_file = args.raw / "ic.stdout.jsonl"
    observed = json.loads(raw_file.read_text())
    assert observed["factor_base_input_hash"] == BASE_HASH
    batch = observed["compact_orbit_point_batch"]
    assert batch["targets_requested"] == 384
    assert batch["regular_scan_policy"] == "target_cyclic_v1"
    relations, failures, queries = (
        batch["relations"], batch["failed_target_points"], batch["query_observations"]
    )
    assert len(queries) == 384 and len(relations) + len(failures) == 384
    recovered = []
    relation_index = failure_index = 0
    for index, (target, query) in enumerate(zip(points, queries)):
        assert tuple(query["target_point"]) == target, index
        if not query["hit"]:
            assert tuple(failures[failure_index]) == target
            failure_index += 1
            recovered.append(None)
            continue
        relation = relations[relation_index]
        relation_index += 1
        assert tuple(relation["target_point"]) == target, index
        assert query["partner_roots"] <= 2 * query["s3_calls"], index
        assert 0 <= query["group_lift_attempts"] <= query["indexed_partner_hits"] <= query["partner_roots"], index
        lift = recover_witness(by_x, target, relation["x_codes"], relation["pinned_intermediates"])
        value = sum(by_point[p][1] * logs[by_point[p][0]] for p in lift) % math.ORDER
        assert math.scalar(math.GENERATOR, value) == target, index
        recovered.append(value)
    assert relation_index == len(relations) and failure_index == len(failures)
    report = {
        "classification": "OPERATIONAL_384_POINT_LOGS_RECOVERED" if not failures else "OPERATIONAL_PARTIAL_POINT_LOGS",
        "points_sha256": sha(args.points), "raw_sha256": sha(raw_file),
        "base_header_sha256": sha(header_file), "training_solution_sha256": sha(logs_file),
        "count": len(relations), "failed_target_points": failures,
        "recovered_scalars": recovered,
        "regular_states": batch["regular_states"], "index_entries": batch["index_entries"],
        "s3_calls_sum": sum(item["s3_calls"] for item in queries),
        "partner_roots_sum": sum(item["partner_roots"] for item in queries),
        "indexed_partner_hits_sum": sum(item["indexed_partner_hits"] for item in queries),
        "group_lift_attempts_sum": sum(item["group_lift_attempts"] for item in queries),
        "query_observations": queries,
    }
    save(args.out, report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("training", "point"), required=True)
    parser.add_argument("--training", type=Path, required=True)
    parser.add_argument("--raw", type=Path)
    parser.add_argument("--points", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.stage == "training":
        training(args)
    else:
        assert args.raw is not None and args.points is not None
        point(args)
    print(json.dumps({"stage": args.stage, "verdict": "PASS", "output_sha256": sha(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
