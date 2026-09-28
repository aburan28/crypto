#!/usr/bin/env python3
"""Independent bit-serial/Fermat replay of a near-threshold slot census."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import resource
import signal
import sys
import time
from array import array
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "rotated_pdp_corpus_20260925" / "verify.py"
spec = importlib.util.spec_from_file_location("pdp_independent_verify", SOURCE)
assert spec and spec.loader
pdp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pdp)
Q = 130873


def highwater() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def target_scalars() -> list[int]:
    result, used, counter = [], set(), 0
    while len(result) != 64:
        phrase = "ECC2K130-ROTATED-SLOT-PLACEMENT-20260925-v1/target/" + str(counter)
        candidate = 1 + int.from_bytes(hashlib.sha256(phrase.encode("ascii")).digest(), "big") % (Q - 1)
        if candidate not in used:
            used.add(candidate)
            result.append(candidate)
        counter += 1
    return result


def independent_rank(rows: list[list[int]]) -> int:
    if not rows:
        return 0
    matrix = [[x % Q for x in row] for row in rows]
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next((i for i in range(rank, len(matrix)) if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], Q - 2, Q)
        matrix[rank] = [(x * inverse) % Q for x in matrix[rank]]
        for i in range(rank + 1, len(matrix)):
            coefficient = matrix[i][column]
            if coefficient:
                matrix[i] = [(x - coefficient * y) % Q
                             for x, y in zip(matrix[i], matrix[rank])]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def replay(position: int, data_dir: Path) -> dict:
    started, cpu_started = time.perf_counter(), time.process_time()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError("verifier wall cap")))
    signal.setitimer(signal.ITIMER_REAL, 600)
    try:
        result = json.loads((data_dir / "result.json").read_text())
        assert result["position"] == position and result["q"] == Q and result["joint_rank"] == 13
        f = pdp.parent_verify.GF(19, 0x80027)
        pdp.rabin_prime_degree(f)
        curve = pdp.parent_verify.E(f)
        assert pdp.source_order(19) == 4 * Q
        h = pdp.generator(curve, Q)
        assert h == (385982, 301867)
        normal = []
        conjugate = 3
        for _ in range(19):
            normal.append(conjugate)
            conjugate = f.square(conjugate)
        assert conjugate == 3 and pdp.rank(normal) == 19 and f.trace(3) == 1
        bases = []
        for i in range(6):
            basis = [normal[i], normal[i + 6]]
            if i == position:
                basis.append(normal[18])
            bases.append(basis)
        assert [pdp.rank(basis) for basis in bases] == result["slot_ranks"]
        assert pdp.rank([element for basis in bases for element in basis]) == 13
        factors = []
        for basis in bases:
            xs = {0}
            for element in basis:
                xs |= {x ^ element for x in list(xs)}
            assert len(xs) == 1 << len(basis)
            factors.append(sorted(p for x in xs for p in pdp.lifts(curve, x)))
        assert [len(factor) for factor in factors] == result["factor_sizes"]

        point_to_scalar = {}
        points = []
        current = None
        for scalar in range(Q):
            assert current not in point_to_scalar
            points.append(current)
            point_to_scalar[current] = scalar
            current = curve.add(current, h)
        assert current is None and len(point_to_scalar) == Q
        factors_as_scalars = [[point_to_scalar[curve.scalar(point, 4)]
                               for point in factor] for factor in factors]
        assert [len(set(values)) for values in factors_as_scalars] == result["slot_signed_projected_counts"]
        slot_budgets = [len({min(value, Q - value) for value in values if value})
                        for values in factors_as_scalars]
        assert slot_budgets == result["slot_pair_column_counts"]
        assert sum(slot_budgets) == result["column_budget"]
        columns = sorted({min(value, Q - value) for values in factors_as_scalars
                          for value in values if value})
        assert columns == result["global_columns"]
        assert len(columns) == result["global_column_count"]

        # Deliberately traverse dense scalar indices and factor choices in the
        # opposite nesting order to the producer's sparse-support recurrence.
        counts = [0] * Q
        counts[0] = 1
        for values in factors_as_scalars:
            new_counts = [0] * Q
            for shift in values:
                for old, multiplicity in enumerate(counts):
                    if multiplicity:
                        index = old + shift
                        if index >= Q:
                            index -= Q
                        new_counts[index] += multiplicity
            counts = new_counts
        assert sum(counts) == math.prod(map(len, factors)) == result["tuple_count"]
        assert sys.byteorder == "little" and array("I").itemsize == 4
        raw = array("I", counts).tobytes()
        assert raw == (data_dir / "counts.u32le").read_bytes()
        assert hashlib.sha256(raw).hexdigest() == result["counts_sha256"]
        support = sum(bool(count) for count in counts)
        m2 = sum(count * count for count in counts)
        assert result["support"] == support and result["misses"] == Q - support
        assert result["second_moment"] == m2
        assert result["collision_pairs"] == (m2 - sum(counts)) // 2
        assert result["cauchy_numerator"] == sum(counts) ** 2
        assert result["cauchy_denominator"] == m2
        assert sum(counts) ** 2 <= m2 * support

        rows = []
        assert len(result["sample_targets"]) == 64
        assert [item["k"] for item in result["sample_targets"]] == target_scalars()
        for item in result["sample_targets"]:
            scalar = item["k"]
            assert item["count"] == counts[scalar]
            if counts[scalar] == 0:
                assert item["indices"] is None and item["row"] is None
                continue
            indices = item["indices"]
            assert len(indices) == 6
            total = None
            coefficient = {column: 0 for column in columns}
            for i, j in enumerate(indices):
                assert 0 <= j < len(factors[i])
                total = curve.add(total, factors[i][j])
                value = factors_as_scalars[i][j]
                if value:
                    col = min(value, Q - value)
                    coefficient[col] += 1 if value == col else -1
            assert curve.scalar(total, 4) == points[scalar]
            sparse = [[column, coefficient[column]] for column in columns if coefficient[column]]
            assert item["row"] == sparse
            dense = [coefficient[column] for column in columns]
            assert sum(a * b for a, b in zip(dense, columns)) % Q == scalar
            rows.append(dense)
        assert len(rows) == result["sample_support"]
        assert independent_rank(rows) == result["sample_rank"]
        assert result["wall_seconds"] <= 300 and result["peak_rss_bytes"] <= 536870912
        verification = {"position": position, "status": "PASS",
                        "all_q_multiplicities_equal": True,
                        "sample_targets_checked": 64,
                        "sample_support": len(rows), "sample_rank": independent_rank(rows),
                        "second_moment": m2,
                        "field_operations": dict(f.ops), "curve_operations": dict(curve.ops),
                        "wall_seconds": time.perf_counter() - started,
                        "cpu_seconds": time.process_time() - cpu_started,
                        "peak_rss_bytes": highwater()}
        assert verification["wall_seconds"] <= 600 and verification["peak_rss_bytes"] <= 536870912
        return verification
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--position", type=int, required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    output = replay(args.position, args.data)
    args.out.write_text(json.dumps(output, sort_keys=True, separators=(",", ":")) + "\n")
