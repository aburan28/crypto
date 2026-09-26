#!/usr/bin/env python3
"""Independent field/group and modular-rank replay of the four rank arms."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ORBIT = HERE.parent / "autolab_orbit_extract_20260924"
REPLAY_PATH = ORBIT / "independent_replay_20260924_codex/replay.py"
RANK_PATH = ORBIT / "cold_batch_rank.py"
SCALARS = HERE / "target_scalars.txt"
POINTS = HERE / "target_points.jsonl"
ARMS = (
    ("native_lex", "native", "lex"),
    ("certified_cyclic", "certified", "target_cyclic_v1"),
    ("certified_lex", "certified", "lex"),
    ("native_cyclic", "native", "target_cyclic_v1"),
)
MASK = (1 << 64) - 1


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def splitmix64(value: int) -> int:
    value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9 & MASK
    value = (value ^ (value >> 27)) * 0x94D049BB133111EB & MASK
    return value ^ (value >> 31)


def expected_start(point: tuple[int, int], policy: str, keys: int) -> int:
    if policy == "lex" or keys == 0:
        return 0
    assert policy == "target_cyclic_v1"
    x, y = point
    mixed = (splitmix64(x ^ 0x6B69632D78353321)
             ^ splitmix64(y ^ 0x6B69632D79353321)
             ^ splitmix64(53 ^ 0x6B69632D6E353321))
    return splitmix64(mixed ^ 0x7461726765747631) % keys


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replay_arm(root: Path, name: str, base: str, policy: str, scalars: list[int],
               points: list[tuple[int, int]], verifier, rank_module) -> dict:
    raw_path = root / name / "producer.stdout.jsonl"
    raw = json.loads(raw_path.read_text())
    assert raw["n"] == 53 and raw["a"] == 0
    assert raw["pair_table_entries"] == raw["pair_selector_variables"] == 0
    assert raw["models_examined"] == 0
    assert raw["decomposition_verdict"] == "UNKNOWN"
    if base == "native":
        assert raw["factor_base_input_hash"] is None
        assert raw["factor_base_input_path"] is None
    else:
        assert raw["factor_base_input_hash"] == "d859319015ea405fd18aee41b51396ce4edcab64ef66265d8edcdeb5e040eb71"
    header = raw["compact_orbit_base_header"]
    assert header["kind"] == "point_defined_factor_base"
    assert header["n"] == 53 and header["a"] == 0
    assert header["orbit_columns"] == 220 and header["factor_base_points"] == 23320
    curve = verifier.Curve(header)
    assert curve.on_curve(curve.generator)
    assert curve.scalar(curve.generator, curve.order) is None
    by_point, reps, lam = rank_module.verify_orbit_labels(curve, header)
    assert len(by_point) == 23320 and len(reps) == 220
    by_x = defaultdict(list)
    for point in by_point:
        by_x[point[0]].append(point)
    assert len(by_x) == 11660 and all(len(pair) == 2 for pair in by_x.values())
    batch = raw["compact_orbit_batch"]
    assert batch["regular_scan_policy"] == policy
    assert batch["targets_requested"] == len(scalars) == len(points) == 512
    queries = batch["query_observations"]
    relations = batch["relations"]
    failures = batch["failed_target_scalars"]
    assert len(queries) == 512
    assert len(relations) + len(failures) == 512
    matrix = rank_module.Echelon(220, curve.order)
    occurrence = Counter()
    first_pair = Counter()
    rank_gains = []
    full_rank_at = None
    solution = None
    relation_index = 0
    failed_index = 0
    nonzero_starts = 0
    all_starts = []
    for target_index, (scalar, point, query) in enumerate(zip(scalars, points, queries), start=1):
        assert query["scalar"] == scalar
        assert query["hit"] in (True, False)
        assert query["regular_scan_start"] == expected_start(point, policy, batch["regular_states"])
        all_starts.append(query["regular_scan_start"])
        nonzero_starts += query["regular_scan_start"] != 0
        assert 0 <= query["regular_keys_visited"] <= batch["regular_states"]
        assert 0 <= query["group_lift_attempts"] <= query["indexed_partner_hits"] <= query["partner_roots"] <= 2 * query["s3_calls"]
        target = curve.scalar(curve.generator, scalar)
        assert target == point and curve.on_curve(target)
        if not query["hit"]:
            assert failures[failed_index] == scalar
            failed_index += 1
            assert query["regular_keys_visited"] == batch["regular_states"]
            continue
        relation = relations[relation_index]
        relation_index += 1
        assert relation["scalar"] == scalar
        lifted = verifier.check_witness(curve, by_x, target,
                                        relation["x_codes"], relation["pinned_intermediates"])
        row = {}
        for item in map(tuple, lifted):
            column, coefficient = by_point[item]
            row[column] = (row.get(column, 0) + coefficient) % curve.order
            occurrence[column] += 1
        first_pair[by_point[tuple(lifted[0])][0]] += 1
        if solution is not None:
            assert sum(coefficient * solution[column] for column, coefficient in row.items()) % curve.order == scalar
        if matrix.insert(row, scalar):
            rank_gains.append(target_index)
            if len(matrix.pivots) == 220 and full_rank_at is None:
                full_rank_at = target_index
                solution = matrix.solution()
    assert relation_index == len(relations) and failed_index == len(failures)
    assert matrix.solution() == solution
    if solution is not None:
        for rep, log in zip(reps, solution):
            assert curve.scalar(rep, curve.order) is None
            assert curve.scalar(curve.generator, log) == rep
    return {
        "name": name, "base": base, "policy": policy,
        "raw_sha256": sha(raw_path),
        "base_arrays_sha256": hashlib.sha256(json.dumps([
            header["factor_base_point_coordinates"],
            header["factor_base_point_labels"],
            header["factor_base_representatives"]
        ], separators=(",", ":")).encode()).hexdigest(),
        "base_header": header,
        "targets_requested": 512,
        "relations_independently_group_replayed": len(relations),
        "failed_target_scalars": failures,
        "rank": len(matrix.pivots), "columns": 220,
        "first_full_rank_target_index": full_rank_at,
        "rank_gain_indices": rank_gains,
        "zero_occurrence_columns": sorted(set(range(220)) - set(occurrence)),
        "least_positive_column_occurrences": min(occurrence.values()) if occurrence else None,
        "first_pair_column0_count": first_pair[0],
        "first_pair_distinct_columns": len(first_pair),
        "nonzero_scan_starts": nonzero_starts,
        "distinct_scan_starts": len(set(all_starts)),
        "s3_calls": sum(q["s3_calls"] for q in queries),
        "regular_keys_visited": sum(q["regular_keys_visited"] for q in queries),
        "query_ms_sum": batch["query_ms_sum"],
        "regular_scan_ms": batch["regular_state_scan_ms"],
        "root_index_build_ms": batch["root_index_build_ms"],
        "base_ms": raw["base_ms"],
        "producer_stage_ms": batch["charged_total_ms"],
        "orbit_labels_replayed": len(by_point),
        "frobenius_eigenvalue": lam,
    }


def run(out: Path):
    verifier = load(REPLAY_PATH, "n53_rank_rotation_group_replay")
    rank_module = load(RANK_PATH, "n53_rank_rotation_rank_replay")
    scalars = [int(line) for line in SCALARS.read_text().splitlines()]
    points = [tuple(json.loads(line)) for line in POINTS.read_text().splitlines()]
    rows = [replay_arm(out, name, base, policy, scalars, points, verifier, rank_module)
            for name, base, policy in ARMS]
    by_name = {row["name"]: row for row in rows}
    assert by_name["native_lex"]["base_arrays_sha256"] == by_name["native_cyclic"]["base_arrays_sha256"]
    assert by_name["certified_lex"]["base_arrays_sha256"] == by_name["certified_cyclic"]["base_arrays_sha256"]
    assert by_name["native_lex"]["base_arrays_sha256"] != by_name["certified_lex"]["base_arrays_sha256"]
    for row in rows:
        del row["base_header"]
    report = {
        "classification": "INDEPENDENT_GROUP_AND_RANK_REPLAY_STAGE_ONLY",
        "scope": "Known-scalar n53 rank/selection control, not logarithm recovery or an IC/rho comparison",
        "target_scalars_sha256": sha(SCALARS),
        "target_points_sha256": sha(POINTS),
        "independent_group_source_sha256": sha(REPLAY_PATH),
        "independent_rank_source_sha256": sha(RANK_PATH),
        "arms": rows,
        "native_rank_delta_rotated_minus_lex": by_name["native_cyclic"]["rank"] - by_name["native_lex"]["rank"],
        "certified_rank_delta_rotated_minus_lex": by_name["certified_cyclic"]["rank"] - by_name["certified_lex"]["rank"],
        "all_four_bases_pairwise_identical_within_base": True,
        "attack_speed_crossover": None,
        "common_operation_unit": None,
    }
    path = out / "audit.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"classification": report["classification"],
                      "native_rank_delta": report["native_rank_delta_rotated_minus_lex"],
                      "certified_rank_delta": report["certified_rank_delta_rotated_minus_lex"]}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.out)
