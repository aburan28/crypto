#!/usr/bin/env python3
"""Independent full-pair/group-law replay of compact four-sum oracle receipts."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import random
import types

HERE = Path(__file__).resolve().parent
SWEEP = HERE.parent / "compact_base_sweep_20260925"
MATH = SWEEP / "evidence/source/pr737_independent_math.py"
BASE_HASHES = {
    (37, 3): "2722f3c7271ea410e47ff9d20b67496a9a28120804f8b790e0729d8c86ea6ddb",
    (41, 8): "cbdd871504ebac2f36e16ce652bc6dfe9623d3fd62723b05f7112d3d94614520",
    (41, 12): "a3856771a5ea1ad3a1e2e262eebcc1b2ccde8293e53af4360268cdcf72ab7fa8",
}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def compact(value) -> bytes:
    return json.dumps(value, separators=(",", ":")).encode()


def read_lines(path: Path) -> list[dict]:
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as stream:
            data = stream.read()
    else:
        data = path.read_bytes()
    return [json.loads(line) for line in data.splitlines() if line.strip()]


def math_module():
    module = types.ModuleType("pr737_independent_math")
    module.__file__ = str(MATH)
    exec(compile(MATH.read_bytes(), str(MATH), "exec"), module.__dict__)
    return module


def euclid_inverse(curve, value: int) -> int:
    """Independent binary-polynomial extended Euclid, reduced mod f."""
    assert value
    u, v = value, curve.modulus
    g, h = 1, 0
    while u != 1:
        shift = u.bit_length() - v.bit_length()
        if shift < 0:
            u, v, g, h, shift = v, u, h, g, -shift
        u ^= v << shift
        g ^= h << shift
    while g.bit_length() > curve.n:
        g ^= curve.modulus << (g.bit_length() - curve.n - 1)
    assert curve.mul(value, g) == 1
    return g


def neg(point):
    return None if point is None else (point[0], point[0] ^ point[1])


def s3(curve, x: int, y: int, z: int) -> int:
    cross = curve.mul(x, y) ^ curve.mul(x, z) ^ curve.mul(y, z)
    return curve.square(cross) ^ curve.mul(curve.mul(x, y), z) ^ 1


def check_group_law(curve_class, n: int, header: dict, points: list[tuple[int, int]]) -> object:
    original = curve_class(n, 0, header["field_modulus_low_terms"])
    accelerated = curve_class(n, 0, header["field_modulus_low_terms"])
    accelerated.inverse = types.MethodType(euclid_inverse, accelerated)
    rng = random.Random(0x4F5241434C45 + n)
    for _ in range(128):
        value = rng.randrange(1, 1 << n)
        assert original.inverse(value) == accelerated.inverse(value)
    order2 = (0, 1)
    assert original.is_on_curve(order2)
    inputs = [(None, None), (None, points[0]), (points[0], None),
              (points[0], neg(points[0])), (points[0], points[0]),
              (order2, order2), (order2, points[0]), (points[0], order2)]
    inputs += [(points[rng.randrange(len(points))], points[rng.randrange(len(points))])
               for _ in range(128)]
    for left, right in inputs:
        assert original.add(left, right) == accelerated.add(left, right)
    return accelerated


def classify_witness(curve, points, tuple4, target, root_owner, root_collisions):
    chosen = [points[i] for i in tuple4]
    assert curve.add(curve.add(chosen[0], chosen[1]),
                     curve.add(chosen[2], chosen[3])) == target
    finite = 0
    collision_proxy = False
    for left_pair, right_pair in (((0, 1), (2, 3)),
                                  ((0, 2), (1, 3)),
                                  ((0, 3), (1, 2))):
        a, b = (chosen[i] for i in left_pair)
        c, d = (chosen[i] for i in right_pair)
        u, v = curve.add(a, b), curve.add(c, d)
        if u is None or v is None:
            continue
        finite += 1
        assert s3(curve, a[0], b[0], u[0]) == 0
        assert s3(curve, c[0], d[0], v[0]) == 0
        assert s3(curve, u[0], v[0], target[0]) == 0
        for p, q, sum_point in ((a, b, u), (c, d, v)):
            endpoint_key = tuple(sorted((p[0], q[0])))
            assert root_owner[sum_point[0]] == endpoint_key or sum_point[0] in root_collisions
            collision_proxy |= sum_point[0] in root_collisions
    opposite_pairs = sum(curve.add(chosen[i], chosen[j]) is None
                         for i in range(4) for j in range(i + 1, 4))
    return {"finite_balanced_partitions": finite,
            "repeated_indices": len(set(tuple4)) < 4,
            "opposite_point_pairs": opposite_pairs,
            "root_x_endpoint_collision_proxy": collision_proxy}


def verify_arm(n: int, r: int, oracle_path: Path, extractor_path: Path) -> dict:
    header = json.loads((HERE / f"base_n{n}_R{r}.json").read_bytes())
    points = [tuple(row) for row in header["factor_base_point_coordinates"]]
    assert sha(compact(header["factor_base_point_coordinates"])) == BASE_HASHES[(n, r)]
    assert len(points) == 2 * n * r and len(set(points)) == len(points)
    math = math_module()
    curve = check_group_law(math.Curve, n, header, points)
    order = int(header["subgroup_order"])
    generator = tuple(header["generator"])
    assert curve.is_on_curve(generator) and curve.scalar(generator, order) is None
    for point in points:
        assert curve.is_on_curve(point) and curve.scalar(point, order) is None
    target_count = 128 if (n, r) == (41, 12) else 512
    full_target_data = (HERE / f"target_points_n{n}.jsonl").read_bytes()
    if (n, r) == (41, 12):
        prefix_data = (HERE / "target_points_n41_R12.jsonl").read_bytes()
        assert prefix_data == b"".join(full_target_data.splitlines(keepends=True)[:128])
        assert sha(prefix_data) == "1b154bdd8aa9dabdb37d2dd5a7bfee69a1fa8284ac5db678ef73eaf2c2f297a5"
        effective_target_file = HERE / "target_points_n41_R12.jsonl"
    else:
        effective_target_file = HERE / f"target_points_n{n}.jsonl"
    target_rows = read_lines(effective_target_file)
    scalar_rows = list(map(int, (HERE / f"target_scalars_n{n}.txt").read_text().split()))[:target_count]
    targets = [tuple(row) for row in target_rows]
    assert len(targets) == len(scalar_rows) == target_count
    for point, scalar in zip(targets, scalar_rows):
        assert curve.is_on_curve(point) and curve.scalar(generator, scalar) == point
    oracle = read_lines(oracle_path)
    extractor = read_lines(extractor_path)
    assert len(oracle) == target_count + 1 and len(extractor) == 1
    oracle_header, outcome_rows = oracle[0], oracle[1:]
    assert oracle_header["kind"] == "complete_four_sum_header"
    assert (oracle_header["n"], oracle_header["R"], oracle_header["target_count"]) == (n, r, target_count)
    produced = extractor[0]
    assert (produced["n"], produced["a"], produced["orbit_columns"]) == (n, 0, r)
    assert sha(compact(produced["compact_orbit_base_header"]["factor_base_point_coordinates"])) == BASE_HASHES[(n, r)]
    batch = produced["compact_orbit_point_batch"]
    observations = batch["query_observations"]
    assert batch["targets_requested"] == len(observations) == target_count
    assert batch["sat_verification_included"] is False
    assert [tuple(row["target_point"]) for row in observations] == targets
    assert [tuple(row["target"]) for row in outcome_rows] == targets
    assert [row["index"] for row in outcome_rows] == list(range(target_count))

    # Independent full construction makes the oracle's pair-count and all
    # occupancy/collision statistics replayable, including infinity.
    pair_counts = Counter()
    root_owner = {}
    root_collisions = set()
    for i, p in enumerate(points):
        for q in points[i:]:
            result = curve.add(p, q)
            pair_counts[result] += 1
            if result is not None:
                x_pair = tuple(sorted((p[0], q[0])))
                prior = root_owner.setdefault(result[0], x_pair)
                if prior != x_pair:
                    root_collisions.add(result[0])
    pair_entries = len(points) * (len(points) + 1) // 2
    assert sum(pair_counts.values()) == oracle_header["pair_entries"] == pair_entries
    assert len(pair_counts) == oracle_header["unique_pair_sums"]
    assert pair_entries - len(pair_counts) == oracle_header["pair_collisions"]
    assert pair_counts[None] == oracle_header["infinity_pair_entries"]
    assert {str(k): v for k, v in sorted(Counter(pair_counts.values()).items())} == oracle_header["bucket_histogram"]
    assert all(row["unique_sum_probes"] == len(pair_counts) == row["lookup_count"] for row in outcome_rows)

    by_x = defaultdict(list)
    for point in points:
        by_x[point[0]].append(point)
    from importlib.util import module_from_spec, spec_from_file_location
    vpath = SWEEP / "verify.py"
    spec = spec_from_file_location("compact_sweep_verify", vpath)
    replay = module_from_spec(spec)
    spec.loader.exec_module(replay)
    rels = batch["relations"]
    assert len(rels) == sum(bool(o["hit"]) for o in observations)
    assert len(batch["failed_target_points"]) == target_count - len(rels)
    for rel in rels:
        point = tuple(rel["target_point"])
        assert point in targets
        replay.check_witness(curve, by_x, point, rel["x_codes"], rel["pinned_intermediates"])

    detail = []
    table = Counter()
    witness_count = 0
    finite_misses = 0
    zero_finite_misses = 0
    collision_proxy_misses = 0
    for i, (row, observation, target) in enumerate(zip(outcome_rows, observations, targets)):
        assert row["kind"] == "complete_four_sum_target" and row["index"] == i
        witnesses = row["witnesses"]
        assert row["member"] == bool(witnesses)
        assert len(witnesses) == row["distinct_four_multisets"]
        assert row["matched_pair_partition_products"] >= len(witnesses)
        assert row["duplicate_partition_products"] == row["matched_pair_partition_products"] - len(witnesses)
        assert all(len(t) == 4 and t == sorted(t) and all(0 <= x < len(points) for x in t)
                   for t in witnesses)
        assert len({tuple(t) for t in witnesses}) == len(witnesses)
        witness_classes = [classify_witness(curve, points, t, target, root_owner, root_collisions)
                           for t in witnesses]
        witness_count += len(witnesses)
        oracle_hit = bool(witnesses)
        extractor_hit = bool(observation["hit"])
        assert not extractor_hit or oracle_hit
        table[(oracle_hit, extractor_hit)] += 1
        if oracle_hit and not extractor_hit:
            if any(c["finite_balanced_partitions"] for c in witness_classes):
                finite_misses += 1
            else:
                zero_finite_misses += 1
            if any(c["root_x_endpoint_collision_proxy"] for c in witness_classes):
                collision_proxy_misses += 1
        detail.append({"index": i, "oracle_hit": oracle_hit,
                       "extractor_hit": extractor_hit,
                       "witness_count": len(witnesses),
                       "witness_classes": witness_classes})

    miss_indices = [i for i, row in enumerate(outcome_rows) if not row["member"]]
    def sample_score(i):
        return hashlib.sha256(f"ECC2K-FOUR-SUM-MISS-SAMPLE-v1/{n}/{r}/{i}".encode()).digest()
    sampled = sorted(miss_indices, key=sample_score)[:16]
    assert len(sampled) == 16
    for index in sampled:
        target = targets[index]
        for pair_sum in pair_counts:
            assert curve.add(target, neg(pair_sum)) not in pair_counts, (n, r, index)
    return {"schema_version": "1.0", "n": n, "R": r,
            "target_count": target_count, "pair_entries": pair_entries,
            "unique_pair_sums": len(pair_counts),
            "pair_collisions": pair_entries - len(pair_counts),
            "infinity_pair_entries": pair_counts[None],
            "distinct_root_x_collision_proxy_count": len(root_collisions),
            "oracle_members": table[(True, True)] + table[(True, False)],
            "extractor_hits": table[(True, True)],
            "oracle_positive_extractor_negative": table[(True, False)],
            "both_negative": table[(False, False)],
            "forbidden_extractor_positive_oracle_negative": table[(False, True)],
            "finite_balanced_missed_targets": finite_misses,
            "zero_finite_balanced_missed_targets": zero_finite_misses,
            "root_x_collision_proxy_missed_targets": collision_proxy_misses,
            "total_distinct_four_multiset_witnesses": witness_count,
            "full_independent_miss_indices": sampled,
            "all_positive_witnesses_replayed": True,
            "all_pair_sums_independently_rebuilt": True,
            "crosschecked_128_euclid_inverses_and_136_group_adds": True,
            "target_detail": detail}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--r", type=int, required=True)
    parser.add_argument("--oracle", type=Path, required=True)
    parser.add_argument("--extractor", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = verify_arm(args.n, args.r, args.oracle, args.extractor)
    data = json.dumps(report, indent=2, sort_keys=True).encode() + b"\n"
    if args.out.exists():
        assert args.out.read_bytes() == data
    else:
        args.out.write_bytes(data)
    print(json.dumps({k: v for k, v in report.items() if k != "target_detail"}, sort_keys=True))


if __name__ == "__main__":
    main()
