#!/usr/bin/env python3
"""Independent GF(2^n) and modular-LA replay for the paired full-rank fixture."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def rows(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def only(items: list[dict], kind: str) -> dict:
    found = [item for item in items if item["kind"] == kind]
    assert len(found) == 1, (kind, len(found))
    return found[0]


class Curve:
    def __init__(self, n: int, a: int, low_terms: list[int]):
        self.n = n
        self.a = a
        self.mask = (1 << n) - 1
        self.reduction = sum(1 << index for index in low_terms)
        self.modulus = (1 << n) | self.reduction

    def mul(self, x: int, y: int) -> int:
        assert 0 <= x <= self.mask and 0 <= y <= self.mask
        z = 0
        for _ in range(self.n):
            if y & 1:
                z ^= x
            y >>= 1
            carry = x >> (self.n - 1)
            x = (x << 1) & self.mask
            if carry:
                x ^= self.reduction
        return z

    def square(self, x: int) -> int:
        return self.mul(x, x)

    def power(self, x: int, exponent: int) -> int:
        value = 1
        while exponent:
            if exponent & 1:
                value = self.mul(value, x)
            x = self.square(x)
            exponent >>= 1
        return value

    def inverse(self, x: int) -> int:
        assert x != 0
        result = self.power(x, (1 << self.n) - 2)
        assert self.mul(result, x) == 1
        return result

    def is_on_curve(self, point: tuple[int, int] | None) -> bool:
        if point is None:
            return True
        x, y = point
        if not (0 <= x <= self.mask and 0 <= y <= self.mask):
            return False
        left = self.square(y) ^ self.mul(x, y)
        right = self.mul(self.square(x), x) ^ self.mul(self.a, self.square(x)) ^ 1
        return left == right

    def add(self, left: tuple[int, int] | None, right: tuple[int, int] | None):
        if left is None:
            return right
        if right is None:
            return left
        x1, y1 = left
        x2, y2 = right
        if x1 == x2:
            if y1 ^ y2 == x1:
                return None
            assert y1 == y2
            if x1 == 0:
                return None
            lam = x1 ^ self.mul(y1, self.inverse(x1))
            x3 = self.square(lam) ^ lam ^ self.a
            y3 = self.square(x1) ^ self.mul(lam ^ 1, x3)
            return (x3, y3)
        lam = self.mul(y1 ^ y2, self.inverse(x1 ^ x2))
        x3 = self.square(lam) ^ lam ^ x1 ^ x2 ^ self.a
        y3 = self.mul(lam, x1 ^ x3) ^ x3 ^ y1
        return (x3, y3)

    def scalar(self, point: tuple[int, int] | None, k: int):
        assert k >= 0
        acc = None
        while k:
            if k & 1:
                acc = self.add(acc, point)
            point = self.add(point, point)
            k >>= 1
        return acc


def pt(raw):
    return None if raw is None else (int(raw[0]), int(raw[1]))


def compact_pt(raw):
    """Producer reserves (0, 0) for infinity and stores affine x+1."""
    x, y = int(raw[0]), int(raw[1])
    return None if x == 0 else (x - 1, y)


def rank_and_solve(matrix: list[list[int]], rhs: list[int], modulus: int, columns: int):
    a = [[v % modulus for v in row] + [b % modulus] for row, b in zip(matrix, rhs)]
    pivot_row = 0
    pivots = []
    for col in range(columns):
        found = next((i for i in range(pivot_row, len(a)) if a[i][col]), None)
        if found is None:
            continue
        a[pivot_row], a[found] = a[found], a[pivot_row]
        inv = pow(a[pivot_row][col], -1, modulus)
        a[pivot_row] = [(v * inv) % modulus for v in a[pivot_row]]
        for i in range(len(a)):
            if i == pivot_row or not a[i][col]:
                continue
            factor = a[i][col]
            a[i] = [(x - factor * y) % modulus for x, y in zip(a[i], a[pivot_row])]
        pivots.append(col)
        pivot_row += 1
    assert all(any(row[:-1]) or row[-1] == 0 for row in a), "inconsistent system"
    if pivot_row < columns:
        return pivot_row, None
    result = [0] * columns
    for i, col in enumerate(pivots):
        result[col] = a[i][-1]
    return pivot_row, result


class IncrementalRank:
    def __init__(self, columns: int, modulus: int):
        self.columns = columns
        self.modulus = modulus
        self.pivots: dict[int, list[int]] = {}

    def add(self, row: list[int]) -> int:
        value = [x % self.modulus for x in row]
        for col in sorted(self.pivots):
            factor = value[col]
            if factor:
                pivot = self.pivots[col]
                value = [(x - factor * y) % self.modulus for x, y in zip(value, pivot)]
        col = next((i for i, x in enumerate(value) if x), None)
        if col is not None:
            inv = pow(value[col], -1, self.modulus)
            self.pivots[col] = [(x * inv) % self.modulus for x in value]
        return len(self.pivots)


def verify(ic_path: Path, rho_path: Path) -> dict:
    ic = rows(ic_path)
    rho_rows = rows(rho_path)
    base = only(ic, "point_defined_factor_base")
    summary = only(ic, "relation_rank_summary")
    rho = only(rho_rows, "rho_public_fixture")
    receipts = [row for row in ic if row["kind"] == "relation_rank_receipt"]
    n, a, r = int(base["n"]), int(base["a"]), int(base["subgroup_order"])
    curve = Curve(n, a, base["field_modulus_low_terms"])
    assert rho["n"] == n and rho["a"] == a and rho["subgroup_order"] == r
    assert rho["field_modulus_low_terms"] == base["field_modulus_low_terms"]
    assert rho["automorphism_size"] == 2 * n
    assert base["base_hash"] == summary["base_hash"]
    assert summary["status"] == "FULL_RANK"
    assert summary["required_surplus_relations"] == 0
    assert len(receipts) == summary["admitted_relations"]
    assert summary["terminal_rank"] == summary["matrix_columns"]
    assert summary["matrix_columns"] == base["orbit_columns"] + 1
    assert summary["target_scalar_constructed"] is False
    assert summary["published_fixture_scalar"] is None
    assert rho["target_scalar_constructed"] is False
    assert rho["published_fixture_scalar"] is None
    assert summary["public_hash_seed"] == rho["public_hash_seed"]
    assert summary["public_hash_counter"] == rho["public_hash_counter"]
    generator = pt(base["generator"])
    q = pt(summary["published_q"])
    assert q == pt(rho["published_q"])
    assert generator == pt(rho["generator"])
    assert curve.is_on_curve(generator) and curve.is_on_curve(q)
    assert curve.scalar(generator, r) is None and curve.scalar(q, r) is None
    recovered = int(summary["recovered_fixture_scalar"])
    assert recovered == int(rho["recovered_fixture_scalar"])
    assert curve.scalar(generator, recovered) == q
    assert rho["verified"] and rho["reference_group_validation"]
    assert rho["walk_steps"] > 0 and rho["ideal_steps"] > 0

    points = [pt(item) for item in base["factor_base_point_coordinates"]]
    labels = [(int(col), int(coeff)) for col, coeff in base["factor_base_point_labels"]]
    reps = [pt(item) for item in base["factor_base_representatives"]]
    assert len(points) == len(labels) == base["factor_base_points"]
    assert len(reps) == base["orbit_columns"]
    assert all(curve.is_on_curve(p) for p in points)
    assert all(curve.is_on_curve(p) for p in reps)
    lam = int(rho["lambda"])
    assert curve.scalar(generator, lam) == (curve.square(generator[0]), curve.square(generator[1]))
    orbit_image = {}
    for col, rep in enumerate(reps):
        current, coefficient = rep, 1
        for _ in range(n):
            for value, image in ((coefficient, current), ((-coefficient) % r, (current[0], current[1] ^ current[0]))):
                key = (col, value)
                assert key not in orbit_image or orbit_image[key] == image
                orbit_image[key] = image
            current = (curve.square(current[0]), curve.square(current[1]))
            coefficient = (coefficient * lam) % r
        assert current == rep and coefficient == 1
    for index, point in enumerate(points):
        col, coeff = labels[index]
        assert 0 <= col < len(reps) and 0 <= coeff < r
        assert orbit_image[(col, coeff)] == point, ("orbit label", index)

    matrix = []
    rhs = []
    incremental = IncrementalRank(len(reps) + 1, r)
    for index, row in enumerate(receipts):
        assert row["accepted_relation"] == index + 1
        assert row["base_hash"] == base["base_hash"]
        assert row["public_hash_seed"] == summary["public_hash_seed"]
        assert row["public_hash_counter"] == summary["public_hash_counter"]
        indices = [int(i) for i in row["factor_point_indices"]]
        assert all(0 <= i < len(points) for i in indices)
        assert [(int(x["column"]), int(x["coefficient"])) for x in row["orbit_labels"]] == [labels[i] for i in indices]
        a_coeff = int(row["coefficient_a"])
        b_coeff = int(row["coefficient_b"])
        target = curve.add(curve.scalar(generator, a_coeff), curve.scalar(q, b_coeff))
        assert target == compact_pt(row["target_point_key"]), ("target", index)
        summed = None
        for i in indices:
            summed = curve.add(summed, points[i])
        assert summed == target, ("relation", index)
        reconstructed = [0] * (len(reps) + 1)
        for i in indices:
            col, coeff = labels[i]
            reconstructed[col] = (reconstructed[col] + coeff) % r
        reconstructed[-1] = (-b_coeff) % r
        assert reconstructed == row["sparse_row"], ("row", index)
        matrix.append(reconstructed)
        rhs.append(a_coeff)
        rank = incremental.add(reconstructed)
        assert row["rank_before"] == (0 if index == 0 else receipts[index-1]["rank_after"])
        assert row["rank_after"] == rank, ("rank", index)
        assert row["rank_incremented"] == (rank > row["rank_before"])
    rank, solution = rank_and_solve(matrix, rhs, r, len(reps) + 1)
    assert rank == len(reps) + 1 and solution is not None
    assert solution[:-1] == summary["factor_base_log_solution"]
    assert solution[-1] == recovered
    for rep, log in zip(reps, solution[:-1]):
        assert curve.scalar(generator, int(log)) == rep
    return {
        "verdict": "PASS",
        "n": n,
        "a": a,
        "public_hash_seed": summary["public_hash_seed"],
        "public_hash_counter": summary["public_hash_counter"],
        "q": list(q),
        "recovered_scalar": recovered,
        "relations_replayed": len(receipts),
        "base_points_orbit_labels_replayed": len(points),
        "terminal_rank": rank,
        "matrix_columns": len(reps) + 1,
        "rho_automorphism_size": rho["automorphism_size"],
        "rho_walk_steps": rho["walk_steps"],
        "rho_ideal_steps": rho["ideal_steps"],
        "rho_steps_over_ideal": rho["walk_steps"] / rho["ideal_steps"],
        "ic_stdout_sha256": hashlib.sha256(ic_path.read_bytes()).hexdigest(),
        "rho_stdout_sha256": hashlib.sha256(rho_path.read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("ic_jsonl", type=Path)
    parser.add_argument("rho_jsonl", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    report = verify(args.ic_jsonl, args.rho_jsonl)
    result = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(result)
    print(result, end="")


if __name__ == "__main__":
    main()
