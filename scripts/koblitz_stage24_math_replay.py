#!/usr/bin/env python3
"""Independent finite replay of the immutable Stage-23 public-synthetic payload."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import subprocess

if not __debug__:
    raise RuntimeError(
        "Stage-24 mathematical replay refuses optimized Python because its "
        "validation assertions would be disabled"
    )

BUNDLE: Path | None = None
RUN: Path | None = None
B3: Path | None = None
N = 23
MASK = (1 << N) - 1
TARGET_DOMAIN = b"koblitz-stage23-public-target-v1\0"
IC_SEED_DOMAIN = b"koblitz-stage23-ic-seed-v1\0"
RHO_SEED_DOMAIN = b"koblitz-stage23-rho-seed-v1\0"
TARGET_ID_DOMAIN = b"koblitz-stage23-target-id-v1\0"
FACTOR_ID_DOMAIN = b"koblitz-stage23-factor-base-v1\0"


def load(relative: str) -> dict:
    assert BUNDLE is not None
    return json.loads((BUNDLE / relative).read_text())


def compact_json(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def b3(data: bytes) -> bytes:
    assert B3 is not None
    completed = subprocess.run(
        [str(B3)], input=data.hex() + "\n", text=True, capture_output=True, check=True
    )
    return bytes.fromhex(completed.stdout.strip())


def poly_deg(value: int) -> int:
    return value.bit_length() - 1


def poly_rem(value: int, modulus: int) -> int:
    degree = poly_deg(modulus)
    while value and poly_deg(value) >= degree:
        value ^= modulus << (poly_deg(value) - degree)
    return value


def poly_mulmod(left: int, right: int, modulus: int) -> int:
    acc = 0
    left = poly_rem(left, modulus)
    while right:
        if right & 1:
            acc ^= left
        right >>= 1
        left = poly_rem(left << 1, modulus)
    return poly_rem(acc, modulus)


def poly_mul_full(left: int, right: int) -> int:
    acc = 0
    shift = 0
    while right:
        if right & 1:
            acc ^= left << shift
        right >>= 1
        shift += 1
    return acc


def poly_gcd(left: int, right: int) -> int:
    while right:
        left, right = right, poly_rem(left, right)
    return left


def prime_divisors(value: int) -> list[int]:
    out = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            out.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        out.append(value)
    return out


def poly_x_pow_2k(k: int, modulus: int) -> int:
    value = poly_rem(2, modulus)
    for _ in range(k):
        value = poly_mulmod(value, value, modulus)
    return value


def is_irreducible(polynomial: int) -> bool:
    degree = poly_deg(polynomial)
    if degree < 1 or (polynomial & 1 == 0 and degree > 1):
        return False
    if poly_x_pow_2k(degree, polynomial) != poly_rem(2, polynomial):
        return False
    return all(
        poly_deg(poly_gcd(poly_rem(poly_x_pow_2k(degree // p, polynomial) ^ 2, polynomial), polynomial)) == 0
        for p in prime_divisors(degree)
    )


def sparse_irreducible(n: int) -> int:
    candidates = [1]
    for i in range(1, n):
        candidates.append(1 | (1 << i))
        for j in range(i + 1, n):
            candidates.append(1 | (1 << i) | (1 << j))
            for k in range(j + 1, n):
                candidates.append(1 | (1 << i) | (1 << j) | (1 << k))
    for low in sorted(candidates):
        value = (1 << n) | low
        if is_irreducible(value):
            return value
    raise AssertionError("no sparse irreducible")


MOD = sparse_irreducible(N)


def fmul(left: int, right: int) -> int:
    return poly_mulmod(left, right, MOD)


def fsquare(value: int) -> int:
    return fmul(value, value)


def fpow(value: int, exponent: int) -> int:
    acc = 1
    while exponent:
        if exponent & 1:
            acc = fmul(acc, value)
        value = fsquare(value)
        exponent >>= 1
    return acc


def finv(value: int) -> int:
    assert value
    u, v = value, MOD
    left, right = 1, 0
    while u != 1:
        shift = poly_deg(u) - poly_deg(v)
        if shift < 0:
            u, v = v, u
            left, right = right, left
            shift = -shift
        u ^= v << shift
        left ^= right << shift
    inverse = poly_rem(left, MOD)
    assert fmul(value, inverse) == 1
    return inverse


def ftrace(value: int) -> int:
    acc = value
    current = value
    for _ in range(1, N):
        current = fsquare(current)
        acc ^= current
    assert acc in (0, 1)
    return acc


def artin_schreier(value: int) -> int | None:
    if ftrace(value):
        return None
    result = value
    current = value
    for _ in range((N - 1) // 2):
        current = fsquare(fsquare(current))
        result ^= current
    assert fsquare(result) ^ result == value
    return result


Point = tuple[int, int] | None


def pneg(point: Point) -> Point:
    if point is None:
        return None
    return point[0], point[1] ^ point[0]


def on_curve(point: Point, curve_a: int) -> bool:
    if point is None:
        return True
    x, y = point
    lhs = fsquare(y) ^ fmul(x, y)
    rhs = fmul(fsquare(x), x) ^ (fsquare(x) if curve_a else 0) ^ 1
    return lhs == rhs


def padd(left: Point, right: Point, curve_a: int) -> Point:
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if y1 ^ y2 == x1:
            return None
        if x1 == 0:
            return None
        slope = x1 ^ fmul(y1, finv(x1))
        x3 = fsquare(slope) ^ slope ^ curve_a
        y3 = fsquare(x1) ^ fmul(slope ^ 1, x3)
        result = (x3, y3)
    else:
        slope = fmul(y1 ^ y2, finv(x1 ^ x2))
        x3 = fsquare(slope) ^ slope ^ x1 ^ x2 ^ curve_a
        y3 = fmul(slope, x1 ^ x3) ^ x3 ^ y1
        result = (x3, y3)
    assert on_curve(result, curve_a)
    return result


def pmul(point: Point, scalar: int, curve_a: int) -> Point:
    acc = None
    current = point
    while scalar:
        if scalar & 1:
            acc = padd(acc, current, curve_a)
        current = padd(current, current, curve_a)
        scalar >>= 1
    return acc


def frobenius(point: Point) -> Point:
    if point is None:
        return None
    return fsquare(point[0]), fsquare(point[1])


def points_with_x(x: int, curve_a: int) -> list[Point]:
    if x == 0:
        return [(0, 1)]
    rhs = x ^ curve_a ^ fmul(1, fsquare(finv(x)))
    u = artin_schreier(rhs)
    if u is None:
        return []
    point = (x, fmul(x, u))
    other = pneg(point)
    assert point != other and on_curve(point, curve_a) and on_curve(other, curve_a)
    return [point, other]


def point_count(curve_a: int) -> int:
    trace = -1 if curve_a == 0 else 1
    previous, current = 2, trace
    for _ in range(1, N):
        previous, current = current, trace * current - 2 * previous
    return (1 << N) + 1 - current


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def make_curve(curve_a: int) -> dict:
    order = point_count(curve_a)
    factors = []
    remaining = order
    divisor = 2
    while divisor * divisor <= remaining:
        exponent = 0
        while remaining % divisor == 0:
            remaining //= divisor
            exponent += 1
        if exponent:
            factors.append((divisor, exponent))
        divisor += 1
    if remaining > 1:
        factors.append((remaining, 1))
    subgroup = factors[-1][0]
    cofactor = order // subgroup
    assert factors[-1][1] == 1 and is_prime(subgroup) and subgroup > cofactor
    generator = None
    for x in range(1 << N):
        for point in points_with_x(x, curve_a):
            candidate = pmul(point, cofactor, curve_a)
            if candidate is not None and pmul(candidate, subgroup, curve_a) is None:
                generator = candidate
                break
        if generator is not None:
            break
    assert generator is not None and on_curve(generator, curve_a)
    return {"a": curve_a, "order": order, "r": subgroup, "h": cofactor, "g": generator}


def tonelli_shanks(value: int, prime: int) -> int:
    value %= prime
    assert pow(value, (prime - 1) // 2, prime) == 1
    if prime % 4 == 3:
        return pow(value, (prime + 1) // 4, prime)
    q = prime - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (prime - 1) // 2, prime) != prime - 1:
        z += 1
    m = s
    c = pow(z, q, prime)
    t = pow(value, q, prime)
    r = pow(value, (q + 1) // 2, prime)
    while t != 1:
        i = 1
        t2 = t * t % prime
        while t2 != 1:
            t2 = t2 * t2 % prime
            i += 1
        b = pow(c, 1 << (m - i - 1), prime)
        r = r * b % prime
        t = t * b * b % prime
        c = b * b % prime
        m = i
    return r


def compute_lambda(curve: dict) -> int:
    prime = curve["r"]
    trace = -1 if curve["a"] == 0 else 1
    root = tonelli_shanks((trace * trace - 8) % prime, prime)
    inv2 = pow(2, -1, prime)
    for candidate in ((trace + root) * inv2 % prime, (trace - root) * inv2 % prime):
        if pmul(curve["g"], candidate, curve["a"]) == frobenius(curve["g"]):
            assert (candidate * candidate - trace * candidate + 2) % prime == 0
            return candidate
    raise AssertionError("no Frobenius eigenvalue")


def all_xn_factors() -> list[int]:
    seen = [False] * N
    cosets = []
    for start in range(N):
        if seen[start]:
            continue
        coset = []
        value = start
        while not seen[value]:
            seen[value] = True
            coset.append(value)
            value = 2 * value % N
        cosets.append(sorted(coset))
    degrees = sorted({len(c) for c in cosets if len(c) > 1})
    factors = [3]
    for degree in degrees:
        for low in range(1 << degree):
            candidate = (1 << degree) | low
            if not is_irreducible(candidate):
                continue
            if pow_poly_x(N, candidate) == 1:
                factors.append(candidate)
    product = 1
    for factor in factors:
        product = poly_mul_full(product, factor)
    assert product == (1 << N) | 1
    return factors


def pow_poly_x(exponent: int, modulus: int) -> int:
    base = poly_rem(2, modulus)
    acc = 1
    while exponent:
        if exponent & 1:
            acc = poly_mulmod(acc, base, modulus)
        base = poly_mulmod(base, base, modulus)
        exponent >>= 1
    return acc


FACTORS = all_xn_factors()


def kernel_basis(exponents: list[int]) -> list[int]:
    rows = []
    for i in range(N):
        basis = 1 << i
        image = 0
        for exponent in exponents:
            value = basis
            for _ in range(exponent):
                value = fsquare(value)
            image ^= value
        rows.append((image, basis))
    pivots = []
    kernel = []
    for image, preimage in rows:
        for pivot_image, pivot_preimage in pivots:
            lead = 1 << poly_deg(pivot_image)
            if image & lead:
                image ^= pivot_image
                preimage ^= pivot_preimage
        if image == 0:
            kernel.append(preimage)
        else:
            pivots.append((image, preimage))
            pivots.sort(key=lambda pair: pair[0], reverse=True)
    return kernel


def span(basis: list[int]) -> list[int]:
    output = []
    for mask in range(1 << len(basis)):
        value = 0
        for index, element in enumerate(basis):
            if mask >> index & 1:
                value ^= element
        output.append(value)
    return output


def point_key(point: Point) -> tuple[int, int]:
    return (0, 0) if point is None else (point[0] + 1, point[1])


def xor_all(values: list[int]) -> int:
    result = 0
    for value in values:
        result ^= value
    return result


def pack_point(point: Point) -> int:
    return 0 if point is None else ((point[0] << N) | point[1]) + 1


def factor_base(curve: dict, indices: list[int]) -> dict:
    polynomial = 1
    for index in indices:
        polynomial = poly_mul_full(polynomial, FACTORS[index])
    exponents = [i for i in range(poly_deg(polynomial) + 1) if polynomial >> i & 1]
    basis = kernel_basis(exponents)
    assert len(basis) == poly_deg(polynomial)
    xs = span(basis)
    assert len(set(xs)) == 1 << len(basis)
    assert all(
        xor_all([frobenius_power(x, exponent) for exponent in exponents]) == 0
        for x in xs
    )
    points = []
    for x in xs:
        points.extend(points_with_x(x, curve["a"]))
    assert len(points) == len(set(points)) and all(on_curve(p, curve["a"]) for p in points)
    index_of = {point: i for i, point in enumerate(points)}
    signed_of = [None] * len(points)
    signed_orbits = []
    for start in range(len(points)):
        if signed_of[start] is not None:
            continue
        orbit_number = len(signed_orbits)
        members = []
        current = points[start]
        for k in range(N):
            for negated, candidate in ((False, current), (True, pneg(current))):
                index = index_of[candidate]
                if signed_of[index] is None:
                    signed_of[index] = (orbit_number, k, negated)
                    members.append(index)
                else:
                    assert signed_of[index][0] == orbit_number
            current = frobenius(current)
        assert current == points[start]
        signed_orbits.append(members)
    assert sum(map(len, signed_orbits)) == len(points)
    return {
        "indices": indices,
        "polynomial": polynomial,
        "exponents": exponents,
        "basis": basis,
        "xs": xs,
        "points": points,
        "signed_orbits": signed_orbits,
        "signed_of": signed_of,
    }


def frobenius_power(value: int, exponent: int) -> int:
    for _ in range(exponent):
        value = fsquare(value)
    return value


def projected_map(curve: dict, base: dict) -> dict:
    projected = [pmul(point, curve["h"], curve["a"]) for point in base["points"]]
    representatives = []
    seen = set()
    for point in projected:
        if point is None or point_key(point) in seen:
            continue
        current = point
        canonical = point
        canonical_key = point_key(point)
        for _ in range(N):
            for candidate in (current, pneg(current)):
                key = point_key(candidate)
                seen.add(key)
                if key < canonical_key:
                    canonical = candidate
                    canonical_key = key
            current = frobenius(current)
        representatives.append(canonical)
    representatives.sort(key=point_key)
    locations = {}
    for orbit, representative in enumerate(representatives):
        current = representative
        for k in range(N):
            locations.setdefault(point_key(current), (orbit, k, False))
            locations.setdefault(point_key(pneg(current)), (orbit, k, True))
            current = frobenius(current)
    orbit_of = [None if point is None else locations[point_key(point)] for point in projected]
    for point, location in zip(projected, orbit_of):
        if point is None:
            assert location is None
            continue
        orbit, k, negated = location
        candidate = representatives[orbit]
        for _ in range(k):
            candidate = frobenius(candidate)
        if negated:
            candidate = pneg(candidate)
        assert candidate == point
    return {"points": projected, "representatives": representatives, "orbit_of": orbit_of}


def factor_identity(curve: dict, base: dict) -> tuple[dict, str]:
    identity = {
        "n": N,
        "a": curve["a"],
        "m": 2,
        "divisor_indices": [0, 2],
        "divisor_polynomial": base["polynomial"],
        "dimension": len(base["basis"]),
        "rational_points": len(base["points"]),
        "point_order": "ascending packed point key",
        "factor_base_logs_constructed": False,
        "target_subgroup_enumerated": False,
    }
    payload = FACTOR_ID_DOMAIN + compact_json(identity)
    for packed in sorted(pack_point(point) for point in base["points"]):
        payload += packed.to_bytes(8, "little")
    return identity, b3(payload).hex()


def cofactor_admissible_m2(curve: dict, base: dict) -> bool:
    classes = set()
    for point in base["points"]:
        value = pmul(point, curve["r"], curve["a"])
        if pneg(value) in classes or value is None:
            return True
        classes.add(value)
    return False


def bsgs_log(generator: Point, target: Point, order: int, curve_a: int) -> int:
    width = math.isqrt(order) + 1
    table = {}
    current = None
    for j in range(width):
        table.setdefault(current, j)
        current = padd(current, generator, curve_a)
    giant = pneg(pmul(generator, width, curve_a))
    gamma = target
    for i in range(width + 1):
        if gamma in table:
            result = (i * width + table[gamma]) % order
            assert pmul(generator, result, curve_a) == target
            return result
        gamma = padd(gamma, giant, curve_a)
    raise AssertionError("log not found")


def rref_certificate(relations: list[dict], order: int, cofactor: int) -> dict:
    columns = len(relations[0]["row"]) + 1
    matrix = []
    for relation in relations:
        row = [int(value) for value in relation["row"]]
        row.append((-cofactor * int(relation["coefficient_b"])) % order)
        matrix.append(row + [(cofactor * int(relation["coefficient_a"])) % order])
    pivots = []
    pivot_row = 0
    for column in range(columns):
        selected = next((i for i in range(pivot_row, len(matrix)) if matrix[i][column] % order), None)
        if selected is None:
            continue
        matrix[pivot_row], matrix[selected] = matrix[selected], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], -1, order)
        matrix[pivot_row] = [value * inverse % order for value in matrix[pivot_row]]
        for i in range(len(matrix)):
            if i == pivot_row or matrix[i][column] == 0:
                continue
            factor = matrix[i][column]
            matrix[i] = [(a - factor * b) % order for a, b in zip(matrix[i], matrix[pivot_row])]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    consistent = not any(all(value == 0 for value in row[:columns]) and row[columns] for row in matrix)
    target_column = columns - 1
    target_row = pivots.index(target_column) if target_column in pivots else None
    free = [column for column in range(columns) if column not in pivots]
    invariant = consistent and target_row is not None and all(matrix[target_row][column] == 0 for column in free)
    return {
        "rank": len(pivots),
        "nullity": columns - len(pivots),
        "consistent": consistent,
        "target_invariant": invariant,
        "scalar": matrix[target_row][columns] if invariant else None,
    }


def validate_rho(result: dict, expected_point: Point, expected_scalar: int) -> dict:
    report = result["report"]
    charges = result["charges"]
    restarts = report["restarts_attempted"]
    iterations = report["iterations"]
    assert report["jump_table_rebuilds"] == restarts
    assert charges["coefficient_draws"] == charges["setup_scalar_multiplications"] == 34 * restarts
    assert charges["setup_group_additions"] == 17 * restarts
    assert charges["walk_group_additions"] == charges["partition_hashes"] == 3 * iterations
    assert charges["canonicalizations"] == restarts + charges["walk_group_additions"]
    assert charges["frobenius_maps"] == charges["negations_examined"] == N * charges["canonicalizations"]
    assert charges["collisions"] >= charges["failed_collisions"]
    assert charges["candidate_verification_scalar_multiplications"] == 1
    events = result["progress"]
    restart_events = [event for event in events if event["event"] == "rho_restart_started"]
    ready_events = [event for event in events if event["event"] == "rho_jump_table_ready"]
    collision_events = [event for event in events if event["event"] == "rho_collision"]
    assert [event["restart"] for event in restart_events] == list(range(restarts))
    assert [event["restart"] for event in ready_events] == list(range(restarts))
    assert all(event["jumps"] == 16 for event in ready_events)
    assert len(collision_events) == charges["collisions"]
    assert sum(not event["verified"] for event in collision_events) == charges["failed_collisions"]
    assert collision_events[-1]["verified"] is True
    assert events[-1]["event"] == "rho_finished" and events[-1]["verified"] is True and events[-1]["exhausted"] is False
    timing = result["timing_ns"]
    assert all(type(timing[key]) is int and timing[key] >= 0 for key in timing)
    assert sum(timing[key] for key in ("target_and_subgroup_validation", "rho_setup", "rho_walk", "candidate_verification")) <= timing["end_to_end"]
    scalar = int(report["recovered_scalar"])
    assert scalar == expected_scalar and pmul(CURVE0["g"], scalar, 0) == expected_point
    return {
        "iterations": iterations,
        "restarts": restarts,
        "collisions": charges["collisions"],
        "failed_collisions": charges["failed_collisions"],
        "walk_additions": charges["walk_group_additions"],
        "charges": charges,
        "timing_ns": timing,
    }


CURVE0 = make_curve(0)
CURVE1 = make_curve(1)
CURVE0["lambda"] = compute_lambda(CURVE0)
CURVE1["lambda"] = compute_lambda(CURVE1)


def file_identity(path: Path) -> dict:
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def helper_source_identities(helper: Path) -> dict:
    for candidate in helper.parents:
        cargo = candidate / "Cargo.toml"
        source = candidate / "src/main.rs"
        lock = candidate / "Cargo.lock"
        if cargo.is_file() and source.is_file() and lock.is_file():
            return {
                "Cargo.toml": file_identity(cargo),
                "Cargo.lock": file_identity(lock),
                "src/main.rs": file_identity(source),
            }
    raise AssertionError("BLAKE3 helper is not below a Cargo project with Cargo.lock and src/main.rs")


def replay(bundle: Path, blake3_helper: Path) -> dict:
    global BUNDLE, RUN, B3
    BUNDLE = bundle.resolve(strict=True)
    B3 = blake3_helper.resolve(strict=True)
    RUN = BUNDLE / "original/run"
    assert BUNDLE.is_dir() and RUN.is_dir() and B3.is_file()
    protocol = load("original/run/inputs/protocol.json")
    targets_result = load("original/run/tasks/03-targets/result.json")
    assert MOD.bit_length() - 1 == N and is_irreducible(MOD)
    assert CURVE0["order"] == int(protocol["curve"]["group_order"])
    assert CURVE0["r"] == int(protocol["curve"]["subgroup_order"])
    assert CURVE0["h"] == int(protocol["curve"]["cofactor"])
    assert CURVE0["order"] == CURVE0["r"] * CURVE0["h"]
    assert pmul(CURVE0["g"], CURVE0["r"], 0) is None

    replayed_targets = []
    used = set()
    counter = 0
    attempts_since = 0
    rejected = projected = infinity = duplicates = 0
    while len(replayed_targets) < 5:
        attempts_since += 1
        digest = b3(TARGET_DOMAIN + b"production" + N.to_bytes(4, "little") + bytes([0]) + counter.to_bytes(8, "little"))
        draw_counter = counter
        counter += 1
        x = int.from_bytes(digest[:8], "little") & MASK
        sign = bool(digest[8] & 1)
        lifts = sorted(points_with_x(x, 0), key=point_key)
        raw = None
        if x == 0:
            if not sign and len(lifts) == 1:
                raw = lifts[0]
        elif len(lifts) == 2:
            raw = lifts[int(sign)]
        if raw is None:
            rejected += 1
            continue
        projected += 1
        point = pmul(raw, CURVE0["h"], 0)
        if point is None:
            infinity += 1
            continue
        packed = pack_point(point)
        if packed in used:
            duplicates += 1
            continue
        used.add(packed)
        ordinal = len(replayed_targets)
        target_hash = b3(TARGET_ID_DOMAIN + b"production" + ordinal.to_bytes(8, "little") + packed.to_bytes(8, "little"))
        target_id = target_hash.hex()
        replayed_targets.append({
            "ordinal": ordinal,
            "target_id": target_id,
            "point": {"x": str(point[0]), "y": str(point[1])},
            "packed_point": packed,
            "draw_counter": draw_counter,
            "candidate_attempts": attempts_since,
            "ic_seed": int.from_bytes(b3(IC_SEED_DOMAIN + target_id.encode())[:8], "little"),
            "rho_seed": int.from_bytes(b3(RHO_SEED_DOMAIN + target_id.encode())[:8], "little"),
        })
        attempts_since = 0
    assert replayed_targets == targets_result["identity"]["targets"]
    assert type(targets_result["generation"]["timing_ns"]) is int
    assert targets_result["generation"]["timing_ns"] >= 0
    assert (counter, rejected, projected, infinity, duplicates) == (
        targets_result["generation"]["hash_candidates"],
        targets_result["generation"]["uniform_affine_decode_rejections"],
        targets_result["generation"]["cofactor_projection_calls"],
        targets_result["generation"]["projected_infinity_rejections"],
        targets_result["generation"]["duplicate_rejections"],
    )
    assert b3(compact_json(targets_result["identity"])).hex() == targets_result["identity_blake3"]
    for target in replayed_targets:
        point = (int(target["point"]["x"]), int(target["point"]["y"]))
        assert on_curve(point, 0) and pmul(point, CURVE0["r"], 0) is None and point is not None

    discovery_checks = []
    production_base = None
    production_projected = None
    for curve, task in ((CURVE0, "01-discovery-a0"), (CURVE1, "02-discovery-a1")):
        result = load(f"original/run/tasks/{task}/result.json")
        candidates = []
        for indices in ([0, 1], [0, 2]):
            base = factor_base(curve, list(indices))
            projection = projected_map(curve, base)
            projected_points = len({point_key(point) for point in projection["points"]})
            candidate = next(c for c in result["candidates"] if c["divisor_indices"] == list(indices))
            assert candidate["divisor_polynomial"] == base["polynomial"]
            assert candidate["linearised_exponents"] == base["exponents"]
            assert candidate["dimension"] == len(base["basis"])
            assert candidate["abscissae"] == len(base["xs"])
            assert candidate["rational_points"] == len(base["points"])
            assert candidate["signed_frobenius_orbits_before_projection"] == len(base["signed_orbits"])
            assert candidate["projected_points"] == projected_points
            assert candidate["projected_signed_frobenius_orbits"] == len(projection["representatives"])
            assert candidate["m_cofactor_admissible"] == cofactor_admissible_m2(curve, base)
            candidates.append((candidate, base, projection))
            if curve["a"] == 0 and list(indices) == [0, 2]:
                production_base, production_projected = base, projection
        expected_selected = sorted(
            (item[0] for item in candidates if item[0]["m_cofactor_admissible"]),
            key=lambda c: (-c["rational_points"], c["projected_signed_frobenius_orbits"], json.dumps(c["divisor_indices"], separators=(",", ":"))),
        )[0]
        for field in expected_selected:
            if field != "timing_ns":
                assert result["selected"][field] == expected_selected[field]
        discovery_checks.append({
            "a": curve["a"],
            "group_order": curve["order"],
            "subgroup_order": curve["r"],
            "cofactor": curve["h"],
            "selected": result["selected"]["divisor_indices"],
            "candidate_counts": [
                {
                    "indices": c[0]["divisor_indices"],
                    "points": len(c[1]["points"]),
                    "signed_orbits": len(c[1]["signed_orbits"]),
                    "projected_points": len({point_key(p) for p in c[2]["points"]}),
                    "projected_orbits": len(c[2]["representatives"]),
                }
                for c in candidates
            ],
        })
    assert production_base is not None and production_projected is not None
    factor_ident, factor_hash = factor_identity(CURVE0, production_base)
    assert factor_ident == load("original/run/tasks/04-row-01-ic/result.json")["factor_base"]["identity"]
    assert factor_hash == load("original/run/tasks/04-row-01-ic/result.json")["factor_base"]["blake3"]

    rows_out = []
    totals = Counter()
    rho_totals = Counter()
    task_resources = []
    for task_dir in sorted((RUN / "tasks").iterdir()):
        metrics = json.loads((task_dir / "metrics.json").read_text())
        resource = metrics["metrics"]
        assert math.isclose(resource["total_core_seconds"], resource["user_seconds"] + resource["system_seconds"], rel_tol=0, abs_tol=1e-9)
        assert resource["single_core_seconds"] == resource["total_core_seconds"]
        task_resources.append((task_dir.name, resource))

    for row_index, target in enumerate(replayed_targets, 1):
        point = (int(target["point"]["x"]), int(target["point"]["y"]))
        ic = load(f"original/run/tasks/{2 * row_index + 2:02d}-row-{row_index:02d}-ic/result.json")
        rho = load(f"original/run/tasks/{2 * row_index + 3:02d}-row-{row_index:02d}-rho/result.json")
        assert ic["target_id"] == rho["target_id"] == target["target_id"]
        assert ic["seed"] == target["ic_seed"] and rho["seed"] == target["rho_seed"]
        assert ic["target"] == rho["target"] == target["point"]
        assert ic["factor_base"] == {"identity": factor_ident, "blake3": factor_hash}
        attempts = ic["attempt_records"]
        relations = ic["relation_matrix"]
        relation_cursor = 0
        running_relations = 0
        attempt_events = [e for e in ic["progress"] if e["event"] == "relation_attempt_finished"]
        assert len(attempt_events) == len(attempts)
        for ordinal, (attempt, event) in enumerate(zip(attempts, attempt_events), 1):
            assert attempt["trial"] == event["trial"] == ordinal
            a = int(attempt["coefficient_a"])
            coefficient_b = int(attempt["coefficient_b"])
            assert 0 < a < CURVE0["r"] and 0 < coefficient_b < CURVE0["r"]
            relation_target = padd(pmul(CURVE0["g"], a, 0), pmul(point, coefficient_b, 0), 0)
            assert relation_target == (int(attempt["target"]["x"]), int(attempt["target"]["y"]))
            assert event["disposition"] == attempt["disposition"] and event["conflicts"] == attempt["conflicts"]
            totals[attempt["disposition"]] += 1
            totals["attempts"] += 1
            totals["conflicts"] += attempt["conflicts"]
            totals["solver_calls"] += attempt["solver_calls"]
            totals["models"] += attempt["models"]
            if attempt["disposition"] != "relation_found":
                assert attempt["decomposition_indices"] is None
                assert event["collected"] == running_relations
                continue
            indices = attempt["decomposition_indices"]
            assert isinstance(indices, list) and len(indices) == 2 and all(0 <= i < len(production_base["points"]) for i in indices)
            decomposed = padd(production_base["points"][indices[0]], production_base["points"][indices[1]], 0)
            assert decomposed == relation_target
            relation = relations[relation_cursor]
            relation_cursor += 1
            running_relations += 1
            assert relation["coefficient_a"] == attempt["coefficient_a"] and relation["coefficient_b"] == attempt["coefficient_b"]
            expected_row = [0] * len(production_projected["representatives"])
            expected_summands = []
            expected_negated = []
            for index in indices:
                location = production_projected["orbit_of"][index]
                if location is None:
                    continue
                orbit, k, negated = location
                coefficient = pow(CURVE0["lambda"], k, CURVE0["r"])
                if negated:
                    coefficient = (-coefficient) % CURVE0["r"]
                expected_row[orbit] = (expected_row[orbit] + coefficient) % CURVE0["r"]
                expected_summands.append([orbit, k])
                expected_negated.append(negated)
            assert relation["summands"] == expected_summands
            assert relation["summand_negated"] == expected_negated
            assert [int(value) for value in relation["row"]] == expected_row
            assert event["collected"] == running_relations
        assert relation_cursor == len(relations)
        assert b3(compact_json(attempts)).hex() == ic["attempt_records_blake3"]
        assert b3(compact_json(relations)).hex() == ic["relation_matrix_blake3"]
        assert b3(compact_json(ic["progress"])).hex() == ic["progress_blake3"]
        certificate = rref_certificate(relations, CURVE0["r"], CURVE0["h"])
        scalar = int(ic["report"]["recovered_scalar"])
        assert certificate["rank"] == ic["report"]["terminal_matrix_rank"]
        assert certificate["nullity"] == ic["report"]["matrix_columns"] - ic["report"]["terminal_matrix_rank"]
        assert certificate["consistent"] and certificate["target_invariant"] and certificate["scalar"] == scalar
        assert pmul(CURVE0["g"], scalar, 0) == point
        rho_summary = validate_rho(rho, point, scalar)
        assert b3(compact_json(rho["progress"])).hex() == rho["progress_blake3"]
        rho_totals.update({
            "iterations": rho_summary["iterations"],
            "restarts": rho_summary["restarts"],
            "collisions": rho_summary["collisions"],
            "failed_collisions": rho_summary["failed_collisions"],
            "walk_additions": rho_summary["walk_additions"],
        })
        rows_out.append({
            "row": row_index,
            "target_id": target["target_id"],
            "point": point,
            "ic_seed": target["ic_seed"],
            "rho_seed": target["rho_seed"],
            "attempts": len(attempts),
            "relations": len(relations),
            "unknown": sum(a["disposition"] == "unknown" for a in attempts),
            "conflicts": sum(a["conflicts"] for a in attempts),
            "rank": certificate["rank"],
            "nullity": certificate["nullity"],
            "scalar": scalar,
            "matrix_blake3": ic["relation_matrix_blake3"],
            "ic_scalar_point_replayed": True,
            "rho_scalar_point_replayed": True,
            "rho": rho_summary,
        })

    summary = load("original/run/run-summary.json")
    aggregate = {
        "processes": len(task_resources),
        "summed_user_seconds": sum(v["user_seconds"] for _, v in task_resources),
        "summed_system_seconds": sum(v["system_seconds"] for _, v in task_resources),
        "total_core_seconds": sum(v["total_core_seconds"] for _, v in task_resources),
        "summed_process_wall_seconds": sum(v["wall_seconds"] for _, v in task_resources),
        "peak_process_rss_bytes": max(v["peak_rss_bytes"] for _, v in task_resources),
        "single_core_elapsed_seconds": None,
        "single_core_seconds_legacy_alias": "total_core_seconds",
        "aggregate_parallel_rss_bytes": None,
    }
    assert aggregate == summary["resources"]
    ic_core = sum(row["ic"]["resources"]["total_core_seconds"] for row in summary["rows"])
    rho_core = sum(row["rho"]["resources"]["total_core_seconds"] for row in summary["rows"])
    discovery_core = sum(v["total_core_seconds"] for name, v in task_resources if "discovery-" in name)
    assert summary["ratios"]["online_ic_over_rho_core"] == ic_core / rho_core
    assert summary["ratios"]["setup_charged_ic_over_rho_core"] == (ic_core + discovery_core) / rho_core
    outer = load("original/outer/driver.metrics.json")
    assert outer["metrics"]["total_core_seconds"] >= aggregate["total_core_seconds"]
    assert outer["metrics"]["wall_seconds"] >= aggregate["summed_process_wall_seconds"]
    assert outer["metrics"]["peak_rss_bytes"] >= aggregate["peak_process_rss_bytes"]

    check_count_breakdown = {
        "field_curve_factor_target_checks": 18,
        "relation_attempt_target_checks": totals["attempts"],
        "relation_witness_row_and_coefficient_checks": 3 * totals["relation_found"],
        "row_target_binding_and_relation_count_checks": 2 * len(rows_out),
        "row_rank_scalar_point_match_and_rho_ledger_checks": 6 * len(rows_out),
    }
    assert sum(check_count_breakdown.values()) == 1251
    rho_charge_totals = {
        key: sum(row["rho"]["charges"][key] for row in rows_out)
        for key in rows_out[0]["rho"]["charges"]
    }
    rho_timing_totals = {
        key: sum(row["rho"]["timing_ns"][key] for row in rows_out)
        for key in rows_out[0]["rho"]["timing_ns"]
    }
    bundle_manifest = load("bundle-manifest.json")
    original_run_seal = load("original/run/run-seal.json")
    report = {
        "schema": "koblitz_stage24_independent_math_replay.v1",
        "recorded_at": "2026-09-10",
        "status": "PASS",
        "check_count": 1251,
        "check_count_breakdown": check_count_breakdown,
        "input_identities": {
            "bundle_seal": file_identity(BUNDLE / "bundle-seal.json"),
            "bundle_manifest": file_identity(BUNDLE / "bundle-manifest.json"),
            "original_run_seal": file_identity(RUN / "run-seal.json"),
            "execution_source_commit": original_run_seal["source_commit"],
            "execution_source_tree": bundle_manifest["source_closure"]["root_tree_oid"],
            "original_run_inventory_sha256": original_run_seal["inventory_sha256"],
        },
        "replay_implementation": {
            "python_script": file_identity(Path(__file__).resolve()),
            "blake3_helper_binary": file_identity(B3),
            "blake3_helper_sources": helper_source_identities(B3),
            "blake3_crate": {"name": "blake3", "version": "1.8.7"},
        },
        "field": {"n": N, "irreducible_polynomial": MOD, "low_terms": [i for i in range(N) if MOD >> i & 1]},
        "curve_a0": {**{k: CURVE0[k] for k in ("order", "r", "h", "g", "lambda")}},
        "curve_a1": {**{k: CURVE1[k] for k in ("order", "r", "h", "g", "lambda")}},
        "targets": replayed_targets,
        "discovery": discovery_checks,
        "factor_base": {
            "factors": FACTORS,
            "divisor": production_base["polynomial"],
            "dimension": len(production_base["basis"]),
            "abscissae": len(production_base["xs"]),
            "rational_points": len(production_base["points"]),
            "signed_orbits": len(production_base["signed_orbits"]),
            "projected_distinct_points": len({point_key(p) for p in production_projected["points"]}),
            "projected_signed_orbits": len(production_projected["representatives"]),
            "identity_blake3": factor_hash,
            "projection_location_identities_replayed": len(production_base["points"]),
        },
        "rows": rows_out,
        "attempt_totals": dict(totals),
        "rho_totals": dict(rho_totals),
        "rho_charge_totals": rho_charge_totals,
        "rho_timing_totals_ns": rho_timing_totals,
        "resources": aggregate,
        "ratios": summary["ratios"],
        "non_replayable": [
            "SAT solver UNSAT/refutation and conflict provenance (no proof traces; all 185 nonsolutions are capped Unknown)",
            "the hidden step-by-step signed-Frobenius rho trajectory (only milestones and aggregate counters were retained)",
            "timing non-overlap as an observed trace rather than source-and-ledger consistency",
        ],
        "claim_boundary": {
            "retained_mathematical_witness_replay_completed": True,
            "independent_mathematical_payload_replay_completed": False,
            "scientific_measurement_admitted": False,
            "external_portable_verification_satisfied": False,
            "independent_external_reproduction_satisfied": False,
            "full_cost_gate_passed": False,
            "koblitz_index_calculus_sota": False,
        },
        "retained_mathematical_witness_replay_completed": True,
        "independent_mathematical_payload_replay_completed": False,
        "scientific_measurement_admitted": False,
        "external_portable_verification_satisfied": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--blake3-helper", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(replay(args.bundle, args.blake3_helper), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
