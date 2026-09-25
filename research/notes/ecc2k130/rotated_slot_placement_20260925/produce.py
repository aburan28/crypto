#!/usr/bin/env python3
"""Exact cyclic group convolution for a frozen n19 slot placement."""
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
PARENT = HERE.parent / "rotated_pdp_corpus_20260925" / "producer.py"
spec = importlib.util.spec_from_file_location("pdp_producer", PARENT)
assert spec and spec.loader
pdp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pdp)
Q = 130873
Point = tuple[int, int] | None


def rss() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def selected_targets() -> list[int]:
    domain = "ECC2K130-ROTATED-SLOT-PLACEMENT-20260925-v1"
    seen, targets, counter = set(), [], 0
    while len(targets) < 64:
        digest = hashlib.sha256(f"{domain}/target/{counter}".encode()).digest()
        k = 1 + int.from_bytes(digest, "big") % (Q - 1)
        if k not in seen:
            seen.add(k)
            targets.append(k)
        counter += 1
    return targets


def rank_mod(rows: list[list[int]], q: int) -> int:
    pivots: dict[int, list[int]] = {}
    for original in rows:
        row = [v % q for v in original]
        for pivot in sorted(pivots):
            if row[pivot]:
                factor = row[pivot]
                row = [(x - factor * y) % q for x, y in zip(row, pivots[pivot])]
        pivot = next((j for j, v in enumerate(row) if v), None)
        if pivot is not None:
            inverse = pow(row[pivot], -1, q)
            pivots[pivot] = [(v * inverse) % q for v in row]
    return len(pivots)


def run(position: int, out: Path) -> None:
    config = json.loads((HERE / "INPUT.json").read_text())
    assert config["positions"] == list(range(6)) and position in config["positions"]
    started, cpu_started = time.perf_counter(), time.process_time()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError("producer wall cap")))
    signal.setitimer(signal.ITIMER_REAL, 300)
    try:
        f = pdp.parent.Field(19, [0, 1, 2, 5])
        f.rabin_prime_degree()
        curve = pdp.parent.Curve(f)
        assert pdp.parent.source_group_order(19) == 4 * Q
        h = pdp.generator(curve, Q)
        assert h == (385982, 301867)
        conjugates = pdp.parent.normal_conjugates(f, 3)
        assert len(conjugates) == 19 and pdp.parent.rank(conjugates) == 19
        bases = [[conjugates[i + 6 * j] for j in range(3)] +
                 ([conjugates[18]] if i == position else []) for i in range(6)]
        assert [pdp.parent.rank(b) for b in bases] == [3 + (i == position) for i in range(6)]
        assert pdp.parent.rank([b for basis in bases for b in basis]) == 19
        factors = [sorted(p for x in pdp.parent.all_x(b) for p in pdp.lifts(curve, x))
                   for b in bases]
        assert all(factors)
        point_to_scalar: dict[Point, int] = {}
        current: Point = None
        for k in range(Q):
            assert current not in point_to_scalar
            point_to_scalar[current] = k
            current = curve.add(current, h)
        assert current is None and len(point_to_scalar) == Q
        deltas = [[point_to_scalar[curve.scalar(p, 4)] for p in factor]
                  for factor in factors]
        slot_signed = [len(set(items)) for items in deltas]
        slot_columns = [len({min(k, Q - k) for k in items if k}) for items in deltas]
        budget = sum(slot_columns)
        columns = sorted({min(k, Q - k) for items in deltas for k in items if k})
        col_index = {k: i for i, k in enumerate(columns)}
        total = math.prod(map(len, factors))
        assert total < 2**32

        counts = [0] * Q
        counts[0] = 1
        active = [0]
        predecessors: list[array] = []
        choices: list[array] = []
        for factor_deltas in deltas:
            next_counts = [0] * Q
            predecessor = array("i", [-1]) * Q
            choice = array("i", [-1]) * Q
            next_active = []
            for old in active:
                multiplicity = counts[old]
                for j, delta in enumerate(factor_deltas):
                    new = old + delta
                    if new >= Q:
                        new -= Q
                    if next_counts[new] == 0:
                        predecessor[new], choice[new] = old, j
                        next_active.append(new)
                    next_counts[new] += multiplicity
            counts, active = next_counts, next_active
            predecessors.append(predecessor)
            choices.append(choice)
        assert sum(counts) == total
        assert sys.byteorder == "little" and array("I").itemsize == 4
        count_bytes = array("I", counts).tobytes()
        count_path = out / "counts.u32le"
        count_path.write_bytes(count_bytes)

        target_rows, dense_rows = [], []
        for target in selected_targets():
            multiplicity = counts[target]
            if multiplicity == 0:
                target_rows.append({"k": target, "count": 0, "indices": None, "row": None})
                continue
            cursor = target
            indices = [0] * 6
            for stage in range(5, -1, -1):
                indices[stage] = choices[stage][cursor]
                cursor = predecessors[stage][cursor]
                assert cursor >= 0 and indices[stage] >= 0
            assert cursor == 0
            acc: Point = None
            row = [0] * len(columns)
            for stage, j in enumerate(indices):
                p = factors[stage][j]
                acc = curve.add(acc, p)
                delta = deltas[stage][j]
                if delta:
                    canonical = min(delta, Q - delta)
                    row[col_index[canonical]] += 1 if delta == canonical else -1
            assert curve.scalar(acc, 4) == curve.scalar(h, target)
            assert sum(v * k for v, k in zip(row, columns)) % Q == target
            dense_rows.append(row)
            target_rows.append({"k": target, "count": multiplicity,
                                "indices": indices,
                                "row": [[columns[j], v] for j, v in enumerate(row) if v]})
        m2 = sum(v * v for v in counts)
        assert (m2 - total) % 2 == 0
        result = {
            "position": position, "n": 19, "q": Q, "beta": 3, "joint_rank": 19,
            "slot_ranks": [len(b) for b in bases],
            "factor_sizes": list(map(len, factors)),
            "slot_signed_projected_counts": slot_signed,
            "slot_pair_column_counts": slot_columns,
            "column_budget": budget, "global_columns": columns,
            "global_column_count": len(columns), "tuple_count": total,
            "support": len(active), "misses": Q - len(active),
            "second_moment": m2, "collision_pairs": (m2 - total) // 2,
            "cauchy_numerator": total * total, "cauchy_denominator": m2,
            "sample_support": len(dense_rows), "sample_rank": rank_mod(dense_rows, Q),
            "sample_targets": target_rows,
            "counts_sha256": hashlib.sha256(count_bytes).hexdigest(),
            "field_operations": dict(f.operations), "curve_operations": dict(curve.operations),
            "wall_seconds": time.perf_counter() - started,
            "cpu_seconds": time.process_time() - cpu_started, "peak_rss_bytes": rss(),
        }
        assert result["wall_seconds"] <= 300 and result["peak_rss_bytes"] <= 536870912
        (out / "result.json").write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--position", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    run(args.position, args.out)
