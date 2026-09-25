#!/usr/bin/env python3
"""Independent exact replay of the frozen beta=3 F0 column census.

This verifier intentionally imports no producer arithmetic.  In particular,
it reconstructs every Gray-ordered coefficient mask from its bits instead of
using the producer's incremental x update, and obtains inverses by polynomial
long-division Euclid rather than the producer's inverse algorithm.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import resource
import signal
import sys
import time
from collections import Counter
from pathlib import Path


Point = tuple[int, int] | None
CAP_WALL = {"pilot": 150, "full": 3600}
CAP_RSS = 768 * 1024 * 1024


def peak_rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def assert_rss_cap() -> None:
    if peak_rss_bytes() > CAP_RSS:
        raise MemoryError(f"independent replay exceeded {CAP_RSS} bytes RSS")


def polynomial_product(a: int, b: int) -> int:
    result = 0
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        b >>= 1
    return result


def polynomial_divmod(a: int, b: int) -> tuple[int, int]:
    if b == 0:
        raise ZeroDivisionError("polynomial division by zero")
    quotient = 0
    degree = b.bit_length() - 1
    while a and a.bit_length() - 1 >= degree:
        shift = a.bit_length() - 1 - degree
        quotient ^= 1 << shift
        a ^= b << shift
    return quotient, a


class Field:
    def __init__(self, n: int, polynomial: int):
        assert n in (13, 19, 131)
        assert polynomial.bit_length() == n + 1 and polynomial & 1
        self.n = n
        self.polynomial = polynomial
        self.mask = (1 << n) - 1
        self.operations: Counter[str] = Counter()
        self.square_byte = tuple(
            sum(((byte >> bit) & 1) << (2 * bit) for bit in range(8))
            for byte in range(256)
        )

    def reduce(self, value: int) -> int:
        self.operations["field_reduction"] += 1
        return polynomial_divmod(value, self.polynomial)[1]

    def mul(self, a: int, b: int) -> int:
        self.operations["field_mul"] += 1
        return self.reduce(polynomial_product(a, b))

    def square(self, a: int) -> int:
        self.operations["field_square"] += 1
        expanded = 0
        byte_index = 0
        while a:
            expanded ^= self.square_byte[a & 255] << (16 * byte_index)
            a >>= 8
            byte_index += 1
        return self.reduce(expanded)

    def inverse(self, a: int) -> int:
        if a == 0:
            raise ZeroDivisionError("zero has no field inverse")
        self.operations["field_inverse"] += 1
        r0, r1 = self.polynomial, a
        t0, t1 = 0, 1
        while r1:
            self.operations["inverse_division"] += 1
            quotient, remainder = polynomial_divmod(r0, r1)
            r0, r1 = r1, remainder
            t0, t1 = t1, t0 ^ polynomial_product(quotient, t1)
        assert r0 == 1, "field polynomial or element is not invertible"
        result = self.reduce(t0)
        return result

    def direct_trace(self, a: int) -> int:
        self.operations["direct_trace"] += 1
        total = 0
        power = a
        for _ in range(self.n):
            total ^= power
            power = self.square(power)
        assert power == a and total in (0, 1)
        return total

    def trace_mask(self) -> int:
        return sum(self.direct_trace(1 << bit) << bit for bit in range(self.n))


def trace_with_mask(trace_mask: int, x: int) -> int:
    return (trace_mask & x).bit_count() & 1


def rank(vectors: list[int]) -> int:
    pivots: dict[int, int] = {}
    for original in vectors:
        value = original
        while value:
            bit = value.bit_length() - 1
            if bit in pivots:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                break
    return len(pivots)


def normal_slice(field: Field, beta: int, m: int, d: int) -> list[int]:
    assert m * d <= field.n
    conjugates: list[int] = []
    value = beta
    for _ in range(field.n):
        conjugates.append(value)
        value = field.square(value)
    assert value == beta and rank(conjugates) == field.n
    assert field.direct_trace(beta) == 1
    assert 1 == xor_values(conjugates)
    basis = [conjugates[m * j] for j in range(d)]
    assert rank(basis) == d
    assert rank(basis + [1]) == d + 1  # x=1 is outside F0.
    return basis


def xor_values(values: list[int]) -> int:
    result = 0
    for value in values:
        result ^= value
    return result


def pivot_rows(basis: list[int]) -> dict[int, int]:
    rows: dict[int, int] = {}
    for original in basis:
        value = original
        while value:
            bit = value.bit_length() - 1
            if bit in rows:
                value ^= rows[bit]
            else:
                rows[bit] = value
                break
        else:
            raise AssertionError("dependent factor basis")
    return rows


def in_span(value: int, rows: dict[int, int]) -> bool:
    while value:
        bit = value.bit_length() - 1
        if bit not in rows:
            return False
        value ^= rows[bit]
    return True


def x_from_natural_mask(basis: list[int], mask: int) -> int:
    value = 0
    bit = 0
    while mask:
        if mask & 1:
            value ^= basis[bit]
        bit += 1
        mask >>= 1
    return value


def projected_x(field: Field, x: int, inverse_x: int) -> int:
    """x([4]P) from x(P); valid for nonzero x other than 1."""
    assert x not in (0, 1)
    inverse_x_squared = field.square(inverse_x)
    x_two = field.square(x) ^ inverse_x_squared
    assert x_two != 0
    inverse_x_two = field.inverse(x_two)
    return field.square(x_two) ^ field.square(inverse_x_two)


def classify_x(field: Field, trace_mask: int, x: int) -> tuple[int, int, int | None]:
    """Return (lift count, projected x sentinel or value, inverse x)."""
    if x == 0:
        return 1, -1, None  # (0,1) is the only affine lift; [4](0,1)=O.
    assert x != 1, "the normal slice must exclude the order-four x=1 fibre"
    inverse_x = field.inverse(x)
    rhs = x ^ field.square(inverse_x)
    if trace_with_mask(trace_mask, rhs):
        return 0, -2, inverse_x
    column = projected_x(field, x, inverse_x)
    assert column != 0  # A nonzero prime-subgroup [4] image has x != 0.
    return 2, column, inverse_x


def half_trace(field: Field, rhs: int) -> int:
    assert field.n & 1 and field.direct_trace(rhs) == 0
    total = 0
    power = rhs
    for _ in range((field.n + 1) // 2):
        total ^= power
        power = field.square(field.square(power))
    assert field.square(total) ^ total == rhs
    return total


class Curve:
    total_operations: Counter[str] = Counter()

    def __init__(self, field: Field):
        self.field = field
        self.operations: Counter[str] = Counter()

    def on_curve(self, point: Point) -> bool:
        if point is None:
            return True
        x, y = point
        f = self.field
        return f.square(y) ^ f.mul(x, y) == f.mul(f.square(x), x) ^ 1

    @staticmethod
    def negate(point: Point) -> Point:
        return None if point is None else (point[0], point[0] ^ point[1])

    def add(self, a: Point, b: Point) -> Point:
        self.operations["point_add"] += 1
        Curve.total_operations["point_add"] += 1
        if a is None:
            return b
        if b is None:
            return a
        x, y = a
        u, v = b
        f = self.field
        if x == u:
            if y ^ v == x:
                return None
            assert y == v and x != 0
            slope = x ^ f.mul(y, f.inverse(x))
            new_x = f.square(slope) ^ slope
            new_y = f.square(x) ^ f.mul(slope ^ 1, new_x)
        else:
            slope = f.mul(y ^ v, f.inverse(x ^ u))
            new_x = f.square(slope) ^ slope ^ x ^ u
            new_y = f.mul(slope, x ^ new_x) ^ new_x ^ y
        return new_x, new_y

    def four(self, point: Point) -> Point:
        self.operations["point_four"] += 1
        Curve.total_operations["point_four"] += 1
        twice = self.add(point, point)
        return self.add(twice, twice)


def torsion_kernel_check(curve: Curve) -> list[Point]:
    torsion: list[Point] = [None, (0, 1), (1, 0), (1, 1)]
    assert len(set(torsion)) == 4
    assert all(curve.on_curve(point) and curve.four(point) is None
               for point in torsion)
    assert curve.add((1, 0), (1, 0)) == (0, 1)
    assert curve.add((1, 1), (1, 1)) == (0, 1)
    assert curve.add((1, 0), (0, 1)) == (1, 1)
    return torsion


def group_law_check(field: Field, trace_mask: int, x: int, lift: int, column: int) -> bool:
    curve = Curve(field)
    torsion = torsion_kernel_check(curve)
    torsion_two = (0, 1)
    if x == 0:
        assert lift == 1 and column == -1
        assert curve.four(torsion_two) is None
        return False
    assert x != 1
    inverse_x = field.inverse(x)
    assert field.mul(x, inverse_x) == 1
    rhs = x ^ field.square(inverse_x)
    assert trace_with_mask(trace_mask, rhs) == (lift == 0)
    assert field.direct_trace(rhs) == (lift == 0)
    if lift == 0:
        assert column == -2
        return False
    assert lift == 2 and column >= 0
    z = half_trace(field, rhs)
    a = (x, field.mul(x, z))
    b = curve.negate(a)
    assert a != b and curve.on_curve(a) and curve.on_curve(b)
    four_a, four_b = curve.four(a), curve.four(b)
    assert four_a is not None and four_b == curve.negate(four_a)
    assert four_a[0] == four_b[0] == column
    assert curve.add(a, b) is None
    translated = curve.add(a, torsion_two)
    assert translated is not None and translated[0] == inverse_x
    assert curve.four(translated) == four_a
    translations = [curve.add(a, shift) for shift in torsion]
    assert len(set(translations)) == 4
    assert all(curve.four(point) == four_a for point in translations)
    return True


def inspect_model(model: dict) -> tuple[Field, list[int], int]:
    field = Field(model["n"], model["poly"])
    assert model["beta"] == 3
    basis = normal_slice(field, model["beta"], model["m"], model["d"])
    return field, basis, field.trace_mask()


def small_control(model: dict) -> dict:
    field, basis, trace_mask = inspect_model(model)
    curve = Curve(field)
    torsion = torsion_kernel_check(curve)
    mask_count = 1 << model["d"]
    projected_x_counts: Counter[int] = Counter()
    signed_partition: dict[int, set[Point]] = {}
    rational_nonzero = 0
    for mask in range(mask_count):
        x = x_from_natural_mask(basis, mask)
        lift, column, inverse_x = classify_x(field, trace_mask, x)
        assert group_law_check(field, trace_mask, x, lift, column) == (lift == 2)
        if lift == 1:
            signed_partition.setdefault(-1, set()).add((0, 1))
        if lift == 2:
            rational_nonzero += 1
            projected_x_counts[column] += 1
            assert inverse_x is not None
            z = half_trace(field, x ^ field.square(inverse_x))
            point = (x, field.mul(x, z))
            signed_partition.setdefault(column, set()).update((point, curve.negate(point)))
    assert signed_partition[-1] == {(0, 1)}
    assert sum(projected_x_counts.values()) == rational_nonzero
    assert max(projected_x_counts.values(), default=0) <= 4
    assert sum(len(points) for points in signed_partition.values()) == 1 + 2 * rational_nonzero
    factor_points = set().union(*signed_partition.values())
    # Two factor points have the same canonical signed [4] column exactly when
    # they lie in one signed E[4] torsion orbit.  Check all four translates.
    for column, points in signed_partition.items():
        representative = next(iter(points))
        signed_orbit = {curve.add(sign, shift)
                        for sign in (representative, curve.negate(representative))
                        for shift in torsion}
        assert signed_orbit & factor_points == points
        for point in points:
            projected = curve.four(point)
            assert (projected is None and column == -1) or (
                projected is not None and projected[0] == column)
    # These two x=1 order-four points are outside F0 but must remain explicit
    # in the full rational [4] kernel and zero-column torsion accounting.
    assert all(point not in factor_points and curve.four(point) is None
               for point in ((1, 0), (1, 1)))
    # The fully archived n13/n19 beta=3 F0 factors have 5 and 7 points.
    assert 1 + 2 * rational_nonzero == {13: 5, 19: 7}[field.n]
    return {"n": field.n, "mask_count": mask_count,
            "rational_nonzero_x": rational_nonzero,
            "physical_points": 1 + 2 * rational_nonzero,
            "nonzero_column_pairs": len(projected_x_counts),
            "max_column_preimages": max(projected_x_counts.values(), default=0),
            "torsion_orbit_partition_checks": len(signed_partition),
            "field_operations": dict(field.operations)}


def verify_samples(input_data: dict, field: Field, basis: list[int], trace_mask: int) -> dict:
    sample_masks = input_data["sample_masks"]
    assert len(sample_masks) == len(set(sample_masks)) == 64
    expected = [0, 1, 1 << 20, (1 << 21) - 1]
    used = set(expected)
    counter = 0
    while len(expected) < 64:
        payload = f"{input_data['domain']}/sample/{counter}".encode("ascii")
        candidate = int.from_bytes(hashlib.sha256(payload).digest(), "big") % (1 << 21)
        counter += 1
        if candidate not in used:
            expected.append(candidate)
            used.add(candidate)
    assert sample_masks == expected
    checked_lifts = 0
    torsion_translation_collisions = 0
    for mask in sample_masks:
        x = x_from_natural_mask(basis, mask)
        lift, column, inverse_x = classify_x(field, trace_mask, x)
        if group_law_check(field, trace_mask, x, lift, column):
            checked_lifts += 1
            assert inverse_x is not None
            if rank(basis + [inverse_x]) == len(basis):
                inverse_lift, inverse_column, _ = classify_x(field, trace_mask, inverse_x)
                assert inverse_lift == 2 and inverse_column == column
                torsion_translation_collisions += 1
    return {"sample_masks": len(sample_masks),
            "sample_rational_nonzero_x": checked_lifts,
            "sample_inverse_x_in_F0": torsion_translation_collisions}


def self_test(input_data: dict) -> dict:
    assert input_data["domain"] == "ECC2K130-N131-BETA3-F0-CENSUS-20260925-v1"
    by_n = {model["n"]: model for model in input_data["models"]}
    assert set(by_n) == {13, 19, 131}
    assert (by_n[13]["poly"], by_n[19]["poly"], by_n[131]["poly"]) == (
        0x201B, 0x80027, 0x800000000000000000000000000002007)
    Curve.total_operations.clear()
    controls = [small_control(by_n[n]) for n in (13, 19)]
    field, basis, trace_mask = inspect_model(by_n[131])
    assert (by_n[131]["m"], by_n[131]["d"]) == (6, 21)
    samples = verify_samples(input_data, field, basis, trace_mask)
    assert_rss_cap()
    return {"status": "pass", "small_controls": controls,
            "n131_group_law_samples": samples,
            "n131_field_operations": dict(field.operations),
            "point_operations": dict(Curve.total_operations),
            "peak_rss_bytes": peak_rss_bytes()}


def replay(input_data: dict, summary: dict) -> dict:
    assert input_data["domain"] == "ECC2K130-N131-BETA3-F0-CENSUS-20260925-v1"
    by_n = {model["n"]: model for model in input_data["models"]}
    assert set(by_n) == {13, 19, 131}
    assert (by_n[13]["poly"], by_n[19]["poly"], by_n[131]["poly"]) == (
        0x201B, 0x80027, 0x800000000000000000000000000002007)
    preflight = self_test(input_data)
    field, basis, trace_mask = inspect_model(by_n[131])
    assert (by_n[131]["m"], by_n[131]["d"]) == (6, 21)
    total = summary["total_masks"]
    assert total in (input_data["pilot_masks"], input_data["full_masks"])
    assert input_data["full_masks"] == 1 << by_n[131]["d"]
    assert input_data["pilot_masks"] == 1 << 15
    digest = hashlib.sha256()
    projected_x_counts: Counter[int] = Counter()
    rows = pivot_rows(basis)
    zero_x_count = one_x_count = liftable_nonzero_x = inverse_pair_collisions = 0
    for ordinal in range(total):
        mask = ordinal ^ (ordinal >> 1)
        x = x_from_natural_mask(basis, mask)
        lift, column, inverse_x = classify_x(field, trace_mask, x)
        digest.update(f"{mask},{x},{lift},{column}\n".encode("ascii"))
        if x == 0:
            zero_x_count += 1
        elif x == 1:
            one_x_count += 1
        if lift == 2:
            liftable_nonzero_x += 1
            projected_x_counts[column] += 1
            assert inverse_x is not None
            if x < inverse_x and in_span(inverse_x, rows):
                inverse_pair_collisions += 1
        if (ordinal + 1) % 8192 == 0:
            assert_rss_cap()
    assert zero_x_count == 1 and one_x_count == 0
    assert sum(projected_x_counts.values()) == liftable_nonzero_x
    assert max(projected_x_counts.values(), default=0) <= 4
    multiplicities = Counter(projected_x_counts.values())
    measured = {"total_masks": total,
                "zero_x_count": zero_x_count,
                "one_x_count": one_x_count,
                "liftable_nonzero_x": liftable_nonzero_x,
                "inverse_pair_collisions": inverse_pair_collisions,
                "nonzero_signed_columns": len(projected_x_counts),
                "multiplicity_histogram": {str(value): multiplicities[value]
                                           for value in range(1, 5)},
                "max_column_preimages": max(projected_x_counts.values(), default=0),
                "row_sha256": digest.hexdigest()}
    for key, value in measured.items():
        assert summary[key] == value, (key, summary[key], value)
    assert summary["physical_factor_points"] == 1 + 2 * liftable_nonzero_x
    assert summary["distinct_projected_points"] == 1 + 2 * len(projected_x_counts)
    assert summary["total_signed_columns_including_O"] == 1 + len(projected_x_counts)
    assert summary["zero_column_count"] == 1
    assert_rss_cap()
    return {"status": "pass", "measured": measured,
            "derived": {"physical_factor_points": 1 + 2 * liftable_nonzero_x,
                        "distinct_projected_points": 1 + 2 * len(projected_x_counts),
                        "total_signed_columns_including_O": 1 + len(projected_x_counts),
                        "zero_column_count": 1},
            "preflight": preflight,
            "small_controls": preflight["small_controls"],
            "n131_group_law_samples": preflight["n131_group_law_samples"],
            "field_operations_enumeration_and_setup": dict(field.operations),
            "point_operations": dict(Curve.total_operations),
            "peak_rss_bytes": peak_rss_bytes(),
            "nonzero_column_definition": "one projected x represents a signed +/-[4] point pair"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if not args.self_test and args.summary is None:
        parser.error("--summary is required unless --self-test is set")
    if args.self_test and args.summary is not None:
        parser.error("--summary is not used with --self-test")
    started = time.perf_counter()
    cpu_started = time.process_time()
    data = json.loads(args.input.read_text())
    summary = None if args.self_test else json.loads(args.summary.read_text())
    if summary is None:
        cap = CAP_WALL["pilot"]
    else:
        assert summary["total_masks"] in (data["pilot_masks"], data["full_masks"])
        cap = CAP_WALL["pilot" if summary["total_masks"] == data["pilot_masks"] else "full"]
    def expire(_signum, _frame):
        raise TimeoutError(f"independent replay exceeded {cap} seconds")
    signal.signal(signal.SIGALRM, expire)
    signal.alarm(cap)
    try:
        result = self_test(data) if args.self_test else replay(data, summary)
        result["wall_seconds"] = time.perf_counter() - started
        result["cpu_seconds"] = time.process_time() - cpu_started
        result["peak_rss_bytes"] = peak_rss_bytes()
        result["input_sha256"] = hashlib.sha256(args.input.read_bytes()).hexdigest()
        if args.summary is not None:
            result["summary_sha256"] = hashlib.sha256(args.summary.read_bytes()).hexdigest()
        assert result["wall_seconds"] <= cap
        assert_rss_cap()
    except BaseException as error:
        result = {"status": "failed", "error": repr(error),
                  "wall_seconds": time.perf_counter() - started,
                  "cpu_seconds": time.process_time() - cpu_started,
                  "peak_rss_bytes": peak_rss_bytes()}
        args.output.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")
        raise
    finally:
        signal.alarm(0)
    args.output.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
