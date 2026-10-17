"""Exact separator-fiber join for a parameterized five-source binary S3 chain.

This solver retains the original affine polynomial-root contract, including
all intermediate branches. It uses no group-point or rational-lift filter.
The join separator is the complete 83-bit (or supplied-degree) field element
u2, not an unchecked fingerprint. Exponential source loops are admitted only
after conservative finite-work and record-count preflights.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

from rootfinder.binaryfield import BinaryField
from research.advance.finite_fiber import BinaryQuadraticFiber, evaluate_s3
from research.advance.projected_sat import _equation_masks, _verify_exported_root


def _integer(value, name):
    if isinstance(value, str):
        try:
            value = int(value, 0)
        except ValueError as error:
            raise ValueError(f"{name} must be a canonical integer or hex string") from error
    if type(value) is not int:
        raise ValueError(f"{name} must be an integer")
    return value


def _count_budget(value, name):
    value = _integer(value, name)
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def _basis_and_field(instance):
    if instance.get("summands", 5) != 5:
        raise ValueError("the separator join requires exactly five summands")
    extension = instance.get("extension_field", {})
    if extension.get("characteristic", 2) != 2:
        raise ValueError("the separator join requires characteristic two")
    degree = _integer(extension.get("degree"), "extension degree")
    modulus = _integer(extension.get("modulus_hex", extension.get("modulus")), "modulus")
    field = BinaryField(degree, modulus)
    dimension = _integer(instance.get("source_dimension"), "source dimension")
    if not 1 <= dimension <= degree:
        raise ValueError("source dimension must be between one and extension degree")
    basis = tuple(field.validate(_integer(value, "source basis")) for value in
                  instance.get("source_basis_hex", [1 << bit for bit in range(dimension)]))
    if len(basis) != dimension:
        raise ValueError("source basis length differs from source dimension")
    pivots = {}
    for original in basis:
        residual = original
        while residual:
            pivot = residual.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = residual
                break
            residual ^= pivots[pivot]
        if residual == 0:
            raise ValueError("source basis must be binary-linearly independent")
    n_vars = 5 * dimension + 3 * degree
    if instance.get("n_vars", n_vars) != n_vars:
        raise ValueError("primary variable count differs from the chain layout")
    if instance.get("source_bits", 5 * dimension) != 5 * dimension:
        raise ValueError("source bit count differs from the chain layout")
    expected_layout = ([{"index": i, "kind": "source_x", "offset": i * dimension,
                         "width": dimension} for i in range(5)] +
                       [{"index": i, "kind": "auxiliary_x", "offset": 5 * dimension + i * degree,
                         "width": degree} for i in range(3)])
    if "variable_layout" in instance and instance["variable_layout"] != expected_layout:
        raise ValueError("unsupported primary variable layout")
    expected_links = [["source_x:0", "source_x:1", "auxiliary_x:0"],
                      ["auxiliary_x:0", "source_x:2", "auxiliary_x:1"],
                      ["auxiliary_x:1", "source_x:3", "auxiliary_x:2"],
                      ["auxiliary_x:2", "source_x:4", "target_x"]]
    if "chain_links" in instance and instance["chain_links"] != expected_links:
        raise ValueError("unsupported S3 chain link order")
    target = field.validate(_integer(instance.get("target", {}).get("x_hex"), "target x"))
    curve_b = field.validate(_integer(instance.get("curve", {}).get("b_hex", "0x1"), "curve B"))
    if curve_b == 0:
        raise ValueError("nonzero curve B is required for the finite-fiber preflight")
    return field, basis, target, curve_b, n_vars


def _decode_domain_value(mask, basis):
    value = 0
    for bit, vector in enumerate(basis):
        if mask & (1 << bit):
            value ^= vector
    return value


def _coordinate_output(mask, field, basis, target, curve_b):
    dimension = len(basis)
    coordinate_mask = (1 << dimension) - 1
    sources = [_decode_domain_value((mask >> (i * dimension)) & coordinate_mask, basis)
               for i in range(5)]
    auxiliary = [(mask >> (5 * dimension + i * field.degree)) & (field.order - 1)
                 for i in range(3)]
    links = ((sources[0], sources[1], auxiliary[0]),
             (auxiliary[0], sources[2], auxiliary[1]),
             (auxiliary[1], sources[3], auxiliary[2]),
             (auxiliary[2], sources[4], target))
    return sum(evaluate_s3(field, *link, curve_b=curve_b) << (i * field.degree)
               for i, link in enumerate(links))


class _TimeExhausted(Exception):
    pass


def solve_chain_join(instance, *, max_source_assignments=1_000_000,
                     max_root_calls=5_000_000, max_records=1_000_000,
                     max_output=10_000, max_seconds=None):
    """Enumerate all original chain roots, or explicitly report unresolved.

    Preflight rejects a bound exceeding a cap before source-domain
    materialization. The right table keeps every separator fiber and every
    duplicate-key record. The output cap limits retention, not counting.
    """
    started = time.perf_counter()
    budgets = {"source_assignments": _count_budget(max_source_assignments, "max_source_assignments"),
               "root_calls": _count_budget(max_root_calls, "max_root_calls"),
               "right_records": _count_budget(max_records, "max_records")}
    max_output = _count_budget(max_output, "max_output")
    if max_seconds is not None and (isinstance(max_seconds, bool) or
                                   not isinstance(max_seconds, (int, float)) or
                                   not math.isfinite(max_seconds) or max_seconds < 0):
        raise ValueError("max_seconds must be finite and nonnegative, or None")
    field, basis, target, curve_b, n_vars = _basis_and_field(instance)
    dimension, source_bits = len(basis), 5 * len(basis)
    k = 1 << dimension
    planned = {"domain_size": k, "source_assignments": k + 2 * k * k + k * k * k,
               "root_calls": k + 3 * k * k + 2 * k * k * k,
               "right_records": 4 * k * k}
    stats = {"planned": planned, "budgets": budgets, "domain_materialized": False,
             "root_calls": 0, "right_records": 0, "right_keys": 0,
             "left_partial_records": 0, "joined_assignments": 0,
             "output_cap": max_output}
    retained, count, complete, reason = [], 0, False, None
    exceeded = [name for name, cap in budgets.items() if planned[name] > cap]
    if exceeded:
        reason = "preflight bound exceeds " + ", ".join(exceeded) + " budget"
    elif max_seconds is not None and time.perf_counter() - started >= max_seconds:
        reason = "charged pipeline time budget exhausted before domain materialization"
    else:
        equations = _equation_masks({**instance, "n_vars": n_vars})
        for receipt in instance.get("native_coordinate_evaluation_checks", []):
            mask = _integer(receipt["assignment_hex"], "receipt assignment")
            if not 0 <= mask < 1 << n_vars:
                raise ValueError("receipt assignment is outside the original coordinate layout")
            expected = _integer(receipt["coordinate_output_hex"], "receipt output")
            if _coordinate_output(mask, field, basis, target, curve_b) != expected:
                raise ValueError("field chain disagrees with native coordinate receipt")
        domain = tuple(_decode_domain_value(mask, basis) for mask in range(k))
        stats["domain_materialized"] = True
        quadratic = BinaryQuadraticFiber(field)

        def check_time():
            if max_seconds is not None and time.perf_counter() - started >= max_seconds:
                raise _TimeExhausted

        def roots(left, right):
            check_time()
            if stats["root_calls"] >= budgets["root_calls"]:
                raise ArithmeticError("preflight quadratic-call bound was violated")
            stats["root_calls"] += 1
            return quadratic.s3_roots(left, right, curve_b)

        right_table = {}
        try:
            for x5_bits, x5 in enumerate(domain):
                for u3 in roots(x5, target):
                    for x4_bits, x4 in enumerate(domain):
                        for u2 in roots(u3, x4):
                            right_table.setdefault(u2, []).append((x4_bits, x5_bits, u3))
                            stats["right_records"] += 1
                            if stats["right_records"] > budgets["right_records"]:
                                raise ArithmeticError("preflight right-record bound was violated")
            stats["right_keys"] = len(right_table)
            for x1_bits, x1 in enumerate(domain):
                for x2_bits, x2 in enumerate(domain):
                    for u1 in roots(x1, x2):
                        for x3_bits, x3 in enumerate(domain):
                            for u2 in roots(u1, x3):
                                stats["left_partial_records"] += 1
                                for x4_bits, x5_bits, u3 in right_table.get(u2, ()):
                                    check_time()
                                    x4, x5 = domain[x4_bits], domain[x5_bits]
                                    links = ((x1, x2, u1), (u1, x3, u2),
                                             (u2, x4, u3), (u3, x5, target))
                                    if any(evaluate_s3(field, *link, curve_b=curve_b) for link in links):
                                        raise ArithmeticError("join assignment failed original S3 verification")
                                    source_mask = (x1_bits | (x2_bits << dimension)
                                                   | (x3_bits << (2 * dimension))
                                                   | (x4_bits << (3 * dimension))
                                                   | (x5_bits << (4 * dimension)))
                                    assignment = (source_mask | (u1 << source_bits)
                                                  | (u2 << (source_bits + field.degree))
                                                  | (u3 << (source_bits + 2 * field.degree)))
                                    count += 1
                                    if len(retained) < max_output:
                                        _verify_exported_root(assignment, equations)
                                        retained.append(assignment)
            complete = True
        except _TimeExhausted:
            reason = "charged pipeline time budget exhausted"
        stats["right_keys"] = len(right_table)
    stats["joined_assignments"] = count
    stats["pipeline_seconds"] = time.perf_counter() - started
    output_complete = complete and len(retained) == count
    return {"schema": "separator-fiber-chain-root-result-v1",
            "algorithm": "exact_s3_separator_fiber_join",
            "scope": "all original Boolean roots of the parameterized five-source S3 chain; no group-decomposition claim",
            "status": "complete" if output_complete else "output_truncated" if complete else "unresolved",
            "search_complete": complete, "output_complete": output_complete,
            "solution_count_exact": count if complete else None,
            "solution_count_lower_bound": count,
            "assignments_hex": [hex(value) for value in sorted(retained)],
            "n_vars": n_vars, "source_bits": source_bits,
            "verified": True, "unresolved_reason": reason, "stats": stats,
            "certificate": "all exact quadratic fibers on both sides of the full-field u2 separator exhausted"
                           if complete else "global exhaustion not established"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("instance")
    parser.add_argument("--max-source-assignments", type=int, default=1_000_000)
    parser.add_argument("--max-root-calls", type=int, default=5_000_000)
    parser.add_argument("--max-records", type=int, default=1_000_000)
    parser.add_argument("--max-output", type=int, default=10_000)
    parser.add_argument("--max-seconds", type=float)
    parser.add_argument("--output")
    args = parser.parse_args()
    raw = Path(args.instance).read_bytes()
    problem = json.loads(raw)
    result = solve_chain_join(problem, max_source_assignments=args.max_source_assignments,
                              max_root_calls=args.max_root_calls, max_records=args.max_records,
                              max_output=args.max_output, max_seconds=args.max_seconds)
    result["input_sha256"] = hashlib.sha256(raw).hexdigest()
    result["instance_id"] = problem.get("instance_id")
    data = json.dumps(result, indent=2) + "\n"
    if args.output:
        Path(args.output).write_text(data)
        print(json.dumps({"status": result["status"], "count": result["solution_count_exact"],
                          "output": args.output, "root_calls": result["stats"]["root_calls"]}))
    else:
        print(data)
    return 0 if result["output_complete"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
