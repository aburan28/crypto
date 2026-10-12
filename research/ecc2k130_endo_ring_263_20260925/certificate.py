#!/usr/bin/env python3
"""Exact kernel-direction certificate and bounded orbit-action cost diagnostic.

The direction test uses the independent bit-polynomial group law from PR #703;
the PR #717 oriented map is imported only for the separate cost diagnostic.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import statistics
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent
PRE = RESEARCH / "ecc2k130_direction_review_20260924"
VELU = RESEARCH / "ecc2k130_oriented_transport_20260924"
REL = RESEARCH / "ecc2k130_relations"
sys.dont_write_bytecode = True
for directory in (PRE, VELU, REL):
    sys.path.insert(0, str(directory))

from redteam_velu_replay import Field, Replay
from fastfield import FastGF2m, IRR131
from relations import Koblitz, CHALLENGE_ELL, CHALLENGE_PX, CHALLENGE_PY
from oriented_velu import BinaryVeluMap

L = 263
R = 680564733841876926932320129493409985129
LAMBDA = 196511074115861092422032515080945363956
FROZEN = {
    "twist_torsion_results.json": "6f1e7b22f3471214d38ec0d3196c88fe9edaf45791e9c05ea67764831fd75368",
    "redteam_velu_replay.py": "c08fbb8de9d54906d783be967a36a6073b2a1d1a3583292badd3adcdd7126790",
    "oriented_velu.py": "a1a656cc32efd612250b1ecf50dff18af79a8f62f970ee79841bfd390603423b",
    "fastfield.py": "b5990fda51700bbfba363251c41f02febfd0fc92edfea4f5ddd647bab472789e",
    "relations.py": "0180501ef8fe00f5c54f910822a28cd86dc3c2c1af164348976c1c2e25cc763f",
}
FILES = {
    "twist_torsion_results.json": PRE / "twist_torsion_results.json",
    "redteam_velu_replay.py": PRE / "redteam_velu_replay.py",
    "oriented_velu.py": VELU / "oriented_velu.py",
    "fastfield.py": REL / "fastfield.py",
    "relations.py": REL / "relations.py",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ensure(condition, message):
    if not condition:
        raise AssertionError(message)


def add_vec(a, b):
    return ((a[0] + b[0]) % L, (a[1] + b[1]) % L)


def mul_vec(a, n):
    return ((a[0] * n) % L, (a[1] * n) % L)


def normalize(v):
    x, y = (v[0] % L, v[1] % L)
    ensure(x or y, "zero vector is not a line")
    return (1, y * pow(x, -1, L) % L) if x else (0, 1)


def next_line(line, columns):
    u, v = columns
    return normalize(add_vec(mul_vec(u, line[0]), mul_vec(v, line[1])))


def cyclic_subgroup(replay, generator, a=1):
    seen = {None: 0}
    point = None
    for k in range(1, L):
        point = replay.add(point, generator, a)
        ensure(point not in seen, "kernel point repeated before 263")
        seen[point] = k
    ensure(replay.add(point, generator, a) is None, "kernel does not close at 263")
    return seen


def exact_checks(frozen, deadline):
    # Derive both F_2 traces by exhaustive point counts.
    count_source = 1 + sum(
        ((y*y ^ x*y) & 1) == ((x*x*x ^ 1) & 1)
        for x in (0, 1) for y in (0, 1))
    count_twist = 1 + sum(
        ((y*y ^ x*y) & 1) == ((x*x*x ^ x*x ^ 1) & 1)
        for x in (0, 1) for y in (0, 1))
    ensure((count_source, count_twist) == (4, 2), "F_2 point count")
    a, b = 1, 0
    for _ in range(131):
        a, b = -2 * b, a - b
    q = 1 << 131
    trace = 2 * a - b
    source_order, twist_order = q + 1 - trace, q + 1 + trace
    ensure((a, b) == (8123678690673902378, 38531015900842053623), "tau recurrence")
    ensure(source_order == 4 * R, "challenge group order")
    ensure(q % L == 1 and a % L == L-1 and b % L == 0 and b % (L*L), "torsion Frobenius")
    ensure(twist_order % (L*L) == 0 and twist_order % (L**3), "twist torsion valuation")
    ensure(CHALLENGE_ELL == R, "public subgroup order")
    ensure((LAMBDA*LAMBDA + LAMBDA + 2) % R == 0 and pow(LAMBDA, 131, R) == 1, "subgroup tau scalar")
    ensure((-7 * L * L) == -484183, "descending discriminant")
    # N(a+b*tau)=a^2-ab+2b^2. On the leaf b is a multiple of 263.
    # Completing the square proves the minimum occurs at |b|=263,
    # a=131 or 132; larger |b| has a larger lower bound.
    min_non_scalar_norm = min(a0*a0 - a0*L + 2*L*L for a0 in (131, 132))
    ensure(min_non_scalar_norm == 121046 and 7*(2*L)**2//4 > min_non_scalar_norm,
           "least non-scalar leaf endomorphism norm")
    result = {
        "source_points_over_F2": count_source,
        "twist_points_over_F2": count_twist,
        "source_trace_over_F2": 3-count_source,
        "source_tau_polynomial": "X^2+X+2",
        "source_endomorphism_order": "Z[tau]=O_K, K=Q(sqrt(-7))",
        "source_endomorphism_conductor": 1,
        "source_endomorphism_discriminant": -7,
        "tau131_A": a,
        "tau131_B": b,
        "Z_pi_conductor": b,
        "q": q,
        "source_order": source_order,
        "twist_order": twist_order,
        "twist_263_primary_structure": "(Z/263Z)^2",
        "r": R,
        "source_prime_subgroup_tau_scalar": LAMBDA,
        "descending_endomorphism_conductor": L,
        "descending_endomorphism_discriminant": -7*L*L,
        "least_non_scalar_leaf_endomorphism_degree": min_non_scalar_norm,
        "least_non_scalar_leaf_endomorphism_witnesses": ["131+263*tau", "132+263*tau"],
        "runs": [],
    }
    p = (CHALLENGE_PX, CHALLENGE_PY)
    for run in frozen["runs"]:
        ensure(time.monotonic() < deadline, "180-second stop")
        field = Field(131, int(run["field_modulus_hex"], 16))
        group = Replay(field)
        ensure(group.on_curve(p, 0) and not group.on_curve(p, 1), "public point curve")
        ensure(group.mul(p, LAMBDA, 0) == (field.sqr(p[0]), field.sqr(p[1])), "public subgroup tau scalar full point")
        basis = [tuple(int(z, 16) for z in pair) for pair in run["basis_on_twist"]]
        columns = [tuple(c) for c in run["frobenius_matrix_columns_mod263"]]
        ensure(len(basis) == len(columns) == 2, "basis dimension")
        for generator, column in zip(basis, columns):
            ensure(group.on_curve(generator, 1), "basis point on twist")
            span = cyclic_subgroup(group, generator)
            ensure(len(span) == L, "basis point order")
            expected = group.add(group.mul(basis[0], column[0], 1), group.mul(basis[1], column[1], 1), 1)
            actual = (field.sqr(generator[0]), field.sqr(generator[1]))
            ensure(expected == actual, "saved matrix column disagrees with independent Frobenius")
        ensure(basis[1] not in cyclic_subgroup(group, basis[0]), "dependent torsion basis")
        mu, mv = columns
        ensure((mu[0] + mv[1]) % L == 1, "twist tau trace")
        ensure((mu[0]*mv[1] - mv[0]*mu[1]) % L == 2, "twist tau norm")
        lines = [(1, s) for s in range(L)] + [(0, 1)]
        seen, orbit_lengths, fixed = set(), [], []
        for line in lines:
            if line in seen:
                continue
            current, orbit = line, []
            while current not in orbit:
                ensure(current not in seen, "orbit collision")
                orbit.append(current)
                seen.add(current)
                current = next_line(current, columns)
            ensure(current == line, "orbit does not close at representative")
            orbit_lengths.append(len(orbit))
            if len(orbit) == 1:
                fixed.append(line)
        ensure(len(seen) == L+1 and sorted(orbit_lengths) == [1, 1, 131, 131], "all line orbits")
        row = {"seed": run["seed"], "twist_tau_matrix_columns": columns,
               "twist_tau_characteristic_polynomial": "X^2-X+2",
               "fixed_lines": fixed, "all_line_orbit_lengths": sorted(orbit_lengths),
               "representatives": []}
        for saved in run["representative_x_maps"]:
            ensure(time.monotonic() < deadline, "180-second stop")
            G = tuple(int(z, 16) for z in saved["kernel_generator_on_twist"])
            ensure(group.on_curve(G, 1), "saved generator on twist")
            span = cyclic_subgroup(group, G)
            tau_G = (field.sqr(G[0]), field.sqr(G[1]))
            eigenvalue = span.get(tau_G)
            line = tuple(saved["line"])
            expected_G = group.add(group.mul(basis[0], line[0], 1),
                                   group.mul(basis[1], line[1], 1), 1)
            ensure(G == expected_G, "saved generator disagrees with saved line")
            matrix_fixed = next_line(line, columns) == line
            ensure((eigenvalue is not None) == matrix_fixed, "generator and matrix disagree on direction")
            ensure(line in lines, "saved line missing")
            ensure((saved["orbit_length"] == 1) == matrix_fixed, "saved orbit length disagrees")
            dest_b = int(saved["codomain_b"], 16)
            ensure((dest_b == 1) == matrix_fixed, "quotient coefficient disagrees with direction")
            xs = {int(x, 16) for x in saved["kernel_abscissae"]}
            ensure(xs == {point[0] for point in span if point is not None}, "kernel abscissa set")
            half = sorted(xs)
            ensure(len(half) == 131, "half-kernel size")
            t = 0
            for x in half:
                t ^= x
            ensure(dest_b == 1 ^ t ^ field.sqr(t), "quotient b from independent field arithmetic")
            if eigenvalue is not None:
                ensure(eigenvalue in (124, 140), "horizontal twist eigenvalue")
            row["representatives"].append({
                "line": list(line), "generator_x": hex(G[0]),
                "matrix_fixed": matrix_fixed,
                "direct_full_point_twist_tau_eigenvalue": eigenvalue,
                "direction": "horizontal" if matrix_fixed else "descending",
                "codomain_b": hex(dest_b),
                "codomain_endomorphism_conductor": 1 if matrix_fixed else L,
                "codomain_endomorphism_discriminant": -7 if matrix_fixed else -7*L*L,
            })
        ensure(sum(z["direction"] == "horizontal" for z in row["representatives"]) == 2,
               "representative balance")
        result["runs"].append(row)
    return result


class CountedField(FastGF2m):
    def __init__(self):
        self.counts = {}
        self._inside_inv = 0
        super().__init__(131, IRR131)
        self.reset()

    def reset(self):
        self.counts = {"mul": 0, "sqr": 0, "inv": 0,
                       "mul_inside_inv": 0, "sqr_inside_inv": 0}

    def mul(self, x, y):
        self.counts["mul"] += 1
        if self._inside_inv:
            self.counts["mul_inside_inv"] += 1
        return super().mul(x, y)

    def sqr(self, x):
        self.counts["sqr"] += 1
        if self._inside_inv:
            self.counts["sqr_inside_inv"] += 1
        return super().sqr(x)

    def inv(self, x):
        self.counts["inv"] += 1
        self._inside_inv += 1
        try:
            return super().inv(x)
        finally:
            self._inside_inv -= 1


def cost_checks(frozen, deadline):
    p = (CHALLENGE_PX, CHALLENGE_PY)
    begin = time.perf_counter()
    field = CountedField()
    field_init = time.perf_counter() - begin
    source, twist = Koblitz(field, a=0), Koblitz(field, a=1)
    rows = []
    for run in frozen["runs"]:
        for saved in run["representative_x_maps"]:
            ensure(time.monotonic() < deadline, "180-second stop")
            G = tuple(int(z, 16) for z in saved["kernel_generator_on_twist"])
            field.reset()
            begin = time.perf_counter()
            phi = BinaryVeluMap.from_generator(source, twist, G, L)
            setup_time = time.perf_counter() - begin
            setup_counts = dict(field.counts)
            field.reset()
            first = phi(p)
            ensure(first[0] == int(saved["image_P_x"], 16), "cost map x value")
            map_counts = dict(field.counts)
            times = []
            for _ in range(11):
                field.reset()
                begin = time.perf_counter()
                ensure(phi(p) == first, "repeated map changed")
                times.append(time.perf_counter() - begin)
                ensure(field.counts == map_counts, "map counts changed")
            field.reset()
            q = (field.sqr(p[0]), field.sqr(p[1]))
            source_tau_counts = dict(field.counts)
            ensure(source.on_curve(q), "source tau point")
            transported_tau = phi(q)
            field.reset()
            leaf_tau = phi.codomain.mul(first, LAMBDA)
            ensure(leaf_tau == transported_tau, "leaf subgroup scalar action")
            leaf_scalar_counts = dict(field.counts)
            leaf_scalar_times = []
            for _ in range(11):
                ensure(time.monotonic() < deadline, "180-second stop")
                field.reset()
                begin = time.perf_counter()
                ensure(phi.codomain.mul(first, LAMBDA) == leaf_tau, "leaf scalar result changed")
                leaf_scalar_times.append(time.perf_counter() - begin)
                ensure(field.counts == leaf_scalar_counts, "leaf scalar counts changed")
            tau_times = []
            for _ in range(11):
                begin = time.perf_counter()
                current = p
                for _ in range(1024):
                    current = (field.sqr(current[0]), field.sqr(current[1]))
                tau_times.append((time.perf_counter() - begin) / 1024)
            rows.append({"seed": run["seed"], "line": saved["line"],
                         "direction": "horizontal" if saved["orbit_length"] == 1 else "descending",
                         "setup_seconds": setup_time, "setup_field_counts": setup_counts,
                         "map_field_counts": map_counts, "map_seconds_11": times,
                         "map_seconds_median": statistics.median(times),
                         "source_tau_field_counts": source_tau_counts,
                         "source_tau_seconds_per_action_11": tau_times,
                         "source_tau_seconds_median": statistics.median(tau_times),
                         "leaf_generic_lambda_field_counts": leaf_scalar_counts,
                         "leaf_generic_lambda_seconds_11": leaf_scalar_times,
                         "leaf_generic_lambda_seconds_median": statistics.median(leaf_scalar_times)})
    return {"shared_field_initialization_seconds": field_init, "rows": rows,
            "scope": "Python implementation stage diagnostic; 11 timed map calls, 11 generic scalar calls, and 11 batches of 1024 tau actions per row, not matched ECDLP work"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error("preserve previous receipt; choose a fresh --out")
    started = time.monotonic()
    receipt = {"schema": "ecc2k130_endo_ring_263_certificate_v1", "status": "RUNNING",
               "source_sha256": digest(Path(__file__)),
               "input_sha256": {name: digest(path) for name, path in FILES.items()},
               "host": platform.platform(), "python": sys.version,
               "command": [sys.executable, *sys.argv],
               "PDP_cost": None, "full_ECDLP_cost": None, "end_to_end_speedup": None}
    try:
        ensure(receipt["input_sha256"] == FROZEN, "frozen input hash changed")
        frozen = json.loads(FILES["twist_torsion_results.json"].read_text())
        receipt["exact"] = exact_checks(frozen, started + 180)
        receipt["orbit_action_cost"] = cost_checks(frozen, started + 180)
        stable = {"exact": receipt["exact"], "field_counts": [
            {"line": row["line"], "seed": row["seed"],
             "setup": row["setup_field_counts"], "map": row["map_field_counts"],
             "source_tau": row["source_tau_field_counts"],
             "leaf_generic_lambda": row["leaf_generic_lambda_field_counts"]}
            for row in receipt["orbit_action_cost"]["rows"]]}
        receipt["deterministic_sha256"] = hashlib.sha256(
            json.dumps(stable, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        ensure(time.monotonic() < started+180, "180-second stop")
        receipt["status"] = "PASS"
    except Exception:
        receipt["status"] = "FAIL"
        receipt["failure"] = traceback.format_exc()
    receipt["elapsed_seconds"] = time.monotonic() - started
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": receipt["status"], "elapsed_seconds": receipt["elapsed_seconds"],
                      "out": str(args.out)}, indent=2))
    if receipt["status"] != "PASS":
        print(receipt["failure"], file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
