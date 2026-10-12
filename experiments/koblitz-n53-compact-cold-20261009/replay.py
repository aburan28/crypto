#!/usr/bin/env python3
"""Independent group, orbit, modular-rank, and scalar replay for compact IC."""

import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open() as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def one(path, kind):
    found = [row for row in rows(path) if row.get("kind") == kind]
    assert len(found) == 1, (path, kind, len(found))
    return found[0]


class BinaryCurve:
    def __init__(self, n, a, low_terms):
        self.n = n
        self.a = a
        self.modulus = (1 << n) | sum(1 << int(i) for i in low_terms)
        self.mask = (1 << n) - 1

    def mul(self, a, b):
        out = 0
        while b:
            if b & 1:
                out ^= a
            b >>= 1
            a <<= 1
            if a & (1 << self.n):
                a ^= self.modulus
        return out & self.mask

    def inv(self, a):
        assert a
        u, v, left, right = a, self.modulus, 1, 0
        while u != 1:
            assert u, "field modulus is not irreducible"
            shift = u.bit_length() - v.bit_length()
            if shift < 0:
                u, v, left, right = v, u, right, left
                shift = -shift
            u ^= v << shift
            left ^= right << shift
        while left.bit_length() > self.n:
            left ^= self.modulus << (left.bit_length() - self.n - 1)
        return left

    def on_curve(self, point):
        if point is None:
            return True
        x, y = point
        return (self.mul(y, y) ^ self.mul(x, y)) == (
            self.mul(self.mul(x, x), x) ^ (self.mul(x, x) if self.a else 0) ^ 1)

    def neg(self, point):
        if point is None:
            return None
        x, y = point
        return x, x ^ y

    def frob(self, point):
        if point is None:
            return None
        x, y = point
        return self.mul(x, x), self.mul(y, y)

    def add(self, p, q):
        if p is None:
            return q
        if q is None:
            return p
        x1, y1 = p
        x2, y2 = q
        if x1 == x2:
            if y1 ^ y2 == x1:
                return None
            assert y1 == y2 and x1
            slope = x1 ^ self.mul(y1, self.inv(x1))
            x3 = self.mul(slope, slope) ^ slope ^ self.a
            y3 = self.mul(x1, x1) ^ self.mul(slope ^ 1, x3)
        else:
            slope = self.mul(y1 ^ y2, self.inv(x1 ^ x2))
            x3 = self.mul(slope, slope) ^ slope ^ x1 ^ x2 ^ self.a
            y3 = self.mul(slope, x1 ^ x3) ^ x3 ^ y1
        return x3, y3

    def scale(self, value, point):
        out = None
        while value:
            if value & 1:
                out = self.add(out, point)
            point = self.add(point, point)
            value >>= 1
        return out


def point(value):
    return None if value is None else tuple(map(int, value))


def check_pair_intermediates(curve, points, codes, intermediate):
    assert [p[0] for p in points] == list(map(int, codes))
    # The x-only S3 index stores both sign branches. `lift` may select the
    # opposite branch when it chooses affine signs, so check root membership
    # for each pair rather than forcing the selected pair sum to one root.
    for offset, root in ((0, intermediate[0]), (2, intermediate[1])):
        roots = {
            total[0]
            for total in (
                curve.add(points[offset], points[offset + 1]),
                curve.add(points[offset], curve.neg(points[offset + 1])),
            )
            if total is not None
        }
        assert int(root) in roots


def check_base(base, rank_header):
    n, a = int(base["n"]), int(base["a"])
    curve = BinaryCurve(n, a, base["field_modulus_low_terms"])
    r = int(base["subgroup_order"])
    pts = [point(p) for p in base["factor_base_point_coordinates"]]
    labels = [tuple(map(int, label)) for label in base["factor_base_point_labels"]]
    reps = [point(p) for p in base["factor_base_representatives"]]
    columns = int(base["orbit_columns"])
    assert len(reps) == columns and len(pts) == len(labels) == 2 * n * columns
    assert len(set(pts)) == len(pts) and None not in pts
    assert rank_header["base_hash"] == base["base_hash"]
    assert rank_header["factor_base_points"] == len(pts)
    assert rank_header["orbit_columns"] == columns
    assert int(rank_header["subgroup_order"]) == r
    assert int(rank_header["n"]) == n and int(rank_header["a"]) == a
    generator = point(rank_header["generator"])
    assert curve.on_curve(generator) and curve.scale(r, generator) is None
    lam = labels[2][1]
    assert 0 < lam < r
    for col, rep in enumerate(reps):
        assert curve.on_curve(rep) and curve.scale(r, rep) is None
        assert curve.scale(lam, rep) == curve.frob(rep)
        current, coefficient = rep, 1
        for k in range(n):
            i = 2 * (col * n + k)
            assert pts[i] == current and pts[i + 1] == curve.neg(current)
            assert labels[i] == (col, coefficient)
            assert labels[i + 1] == (col, (-coefficient) % r)
            assert curve.on_curve(current)
            current = curve.frob(current)
            coefficient = coefficient * lam % r
        assert current == rep and coefficient == 1
    return curve, r, generator, pts, labels, reps


def insert_row(row, pivots, r, columns):
    row = [int(value) % r for value in row]
    assert len(row) == columns + 1
    for col in range(columns):
        value = row[col]
        if value == 0:
            continue
        if col in pivots:
            pivot = pivots[col]
            for j in range(col, columns + 1):
                row[j] = (row[j] - value * pivot[j]) % r
        else:
            inverse = pow(value, -1, r)
            for j in range(col, columns + 1):
                row[j] = row[j] * inverse % r
            pivots[col] = row
            return True
    assert row[columns] == 0, "inconsistent relation system"
    return False


def solve(pivots, r, columns):
    assert len(pivots) == columns
    logs = [0] * columns
    for col in range(columns - 1, -1, -1):
        row = pivots[col]
        logs[col] = (row[columns] - sum(row[j] * logs[j] for j in range(col + 1, columns))) % r
    return logs


def group_sum(curve, chosen):
    out = None
    for item in chosen:
        out = curve.add(out, item)
    return out


def replay_rank(path, curve, r, generator, pts, labels, reps, rank_seed):
    iterator = rows(path)
    header = next(iterator)
    assert header["kind"] == "compact_orbit_rank_header"
    columns = len(reps)
    pivots = {}
    equations = []
    attempts = relations = failures = without_gain = 0
    solution = None
    for row in iterator:
        if row["kind"] == "compact_orbit_rank_solution":
            assert solution is None
            solution = row
            continue
        assert solution is None and row["kind"] == "compact_orbit_rank_attempt"
        assert int(row["attempt_index"]) == attempts
        rank_seed = (rank_seed * 6364136223846793005 + 1442695040888963407) & ((1 << 64) - 1)
        scalar = (rank_seed >> 11) % (r - 1) + 1
        assert int(row["scalar"]) == scalar
        column = next(col for col in range(columns) if col not in pivots)
        assert int(row["pivotless_column"]) == column
        assert int(row["rank_before"]) == len(pivots)
        if row["found"]:
            indices = list(map(int, row["point_indices"]))
            assert len(indices) == 4 and all(0 <= i < len(pts) for i in indices)
            chosen = [pts[i] for i in indices]
            target = point(row["target"])
            assert curve.add(target, reps[column]) == curve.scale(scalar, generator)
            assert group_sum(curve, chosen) == target
            check_pair_intermediates(curve, chosen, row["x_codes"], row["pinned_intermediates"])
            expected = [0] * (columns + 1)
            for i in indices:
                col, coeff = labels[i]
                expected[col] = (expected[col] + coeff) % r
            expected[column] = (expected[column] + 1) % r
            expected[columns] = scalar
            assert list(map(int, row["row"])) == expected
            gained = insert_row(expected, pivots, r, columns)
            assert row["gained"] is gained
            relations += 1
            without_gain += not gained
            equations.append(expected)
        else:
            assert int(row["rank_after"]) == len(pivots)
            failures += 1
        attempts += 1
        assert int(row["rank_after"]) == len(pivots)
    assert solution is not None and len(pivots) == columns
    assert int(solution["rank"]) == columns
    assert int(solution["attempts"]) == attempts
    assert int(solution["relations"]) == relations
    assert int(solution["failures"]) == failures
    assert int(solution["rows_without_gain"]) == without_gain
    logs = solve(pivots, r, columns)
    assert list(map(int, solution["logs"])) == logs
    for equation in equations:
        assert sum(a * b for a, b in zip(equation[:columns], logs)) % r == equation[columns]
    for col, rep in enumerate(reps):
        assert curve.scale(logs[col], generator) == rep
    return {"attempts": attempts, "relations": relations, "failures": failures,
            "rows_without_gain": without_gain, "rank": len(pivots)}, logs


def replay_target(path, curve, r, generator, pts, labels, logs, workload=None):
    target_row = one(path, "compact_orbit_dlp_target")
    target = point(target_row["published_q"])
    assert target == point(target_row["target"])
    if workload is not None:
        assert target == point(workload["primary_target"])
    indices = list(map(int, target_row["point_indices"]))
    assert len(indices) == 4 and all(0 <= i < len(pts) for i in indices)
    chosen = [pts[i] for i in indices]
    assert group_sum(curve, chosen) == target
    check_pair_intermediates(curve, chosen, target_row["x_codes"], target_row["pinned_intermediates"])
    scalar = sum(labels[i][1] * logs[labels[i][0]] for i in indices) % r
    assert scalar == int(target_row["recovered_scalar"])
    assert curve.scale(scalar, generator) == target
    assert target_row["group_verified"] is True
    if workload is not None:
        assert scalar == int(workload["verification_scalar"])
    return target, scalar


def replay(run_dir, workload_path=None, require_receipt=False):
    run_dir = Path(run_dir)
    base = one(run_dir / "base.jsonl", "point_defined_factor_base")
    rank_header = one(run_dir / "rank.jsonl", "compact_orbit_rank_header")
    curve, r, generator, pts, labels, reps = check_base(base, rank_header)
    workload = json.loads(Path(workload_path).read_text()) if workload_path else None
    if workload is not None:
        assert int(workload["n"]) == curve.n and int(workload["a"]) == curve.a
        assert int(workload["subgroup_order"]) == r
        assert point(workload["generator"]) == generator
        assert list(map(int, workload["field_modulus_low_terms"])) == base["field_modulus_low_terms"]
    seed = int(workload["rank_seed"]) if workload is not None else 20261009
    rank, logs = replay_rank(run_dir / "rank.jsonl", curve, r, generator, pts, labels, reps, seed)
    target, scalar = replay_target(run_dir / "ic-target.jsonl", curve, r, generator, pts, labels, logs, workload)
    summary = one(run_dir / "ic.stdout.jsonl", "compact_orbit_dlp_summary")
    assert int(summary["rank"]) == rank["rank"]
    assert int(summary["rank_attempts"]) == rank["attempts"]
    assert int(summary["rank_relations"]) == rank["relations"]
    assert int(summary["rank_failures"]) == rank["failures"]
    assert int(summary["factor_base_points"]) == len(pts)
    assert summary["base_hash"] == base["base_hash"]
    if (run_dir / "rho.stdout.jsonl").exists():
        rho = one(run_dir / "rho.stdout.jsonl", "rho_public_fixture")
        assert point(rho["published_q"]) == target
        assert int(rho["recovered_fixture_scalar"]) == scalar and rho["verified"] is True
        assert curve.scale(int(rho["recovered_fixture_scalar"]), generator) == target
    if require_receipt:
        receipt = json.loads((run_dir / "receipt.json").read_text())
        assert receipt["ic"]["status"] == receipt["rho"]["status"] == "success"
        assert receipt["analysis"]["status"] == "pending_independent_relation_replay"
        for name, digest in receipt["files_sha256"].items():
            assert sha(run_dir / name) == digest, name
    return {"status": "PASS", "curve_n": curve.n, "subgroup_order": r,
            "factor_base_points": len(pts), "orbit_columns": len(reps),
            "rank": rank, "target": target, "recovered_scalar": scalar,
            "base_sha256": sha(run_dir / "base.jsonl"),
            "rank_trace_sha256": sha(run_dir / "rank.jsonl"),
            "target_sha256": sha(run_dir / "ic-target.jsonl")}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--workload", type=Path)
    parser.add_argument("--require-receipt", action="store_true")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = replay(args.run_dir, args.workload, args.require_receipt)
    output = json.dumps(result, sort_keys=True, indent=2) + "\n"
    if args.out:
        with args.out.open("x") as handle:
            handle.write(output)
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
