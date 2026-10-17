"""Exact source-projected SAT with deterministic finite-field fiber lifting.

Input is the parameterized neutral S3-chain export, not an arbitrary system
whose equations may differ from its curve metadata. Native export hashes and
evaluation receipts are checked; every retained original Boolean root is
independently checked against the exported equations. No GB is computed.
"""

import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys
import time

from research.advance.field_circuit import compile_s3_chain
from research.advance.finite_fiber import BinaryQuadraticFiber, evaluate_s3, five_source_intermediate_fiber


def load_backend():
    local = Path(__file__).resolve().parent / "deps_5_11_21"
    if local.is_dir() and str(local) not in sys.path:
        sys.path.insert(0, str(local))
    try:
        module = importlib.import_module("pycryptosat")
    except (ImportError, OSError) as error:
        raise RuntimeError("Install the pinned pycryptosat 5.11.21 backend; see ADVANCE.md") from error
    if module.__version__ != "5.11.21":
        raise RuntimeError("This experiment requires the registered pycryptosat 5.11.21 backend")
    return module


def _integer(value):
    if isinstance(value, str):
        return int(value, 0)
    if type(value) is not int:
        raise ValueError("canonical integer or hex string required")
    return value


def _decode_sources(mask, dimension, basis):
    values = []
    for block in range(5):
        value = 0
        for bit, vector in enumerate(basis):
            if mask & (1 << (block * dimension + bit)):
                value ^= vector
        values.append(value)
    return tuple(values)


def _decode_auxiliary(mask, degree, source_bits):
    field_mask = (1 << degree) - 1
    return tuple((mask >> (source_bits + i * degree)) & field_mask for i in range(3))


def _check_links(field, sources, auxiliary, target, curve_b):
    x1, x2, x3, x4, x5 = sources
    u1, u2, u3 = auxiliary
    links = ((x1, x2, u1), (u1, x3, u2), (u2, x4, u3), (u3, x5, target))
    if any(evaluate_s3(field, *link, curve_b=curve_b) for link in links):
        raise ArithmeticError("SAT/fiber result failed original field equation verification")


def _equation_masks(instance):
    n = instance["n_vars"]
    equations = instance.get("equations")
    if equations is None:
        return None
    if len(equations) != 4 * instance["extension_field"]["degree"]:
        raise ValueError("coordinate equation count differs from the four-link chain")
    expected = instance.get("equations_sha256")
    digest = hashlib.sha256(json.dumps(equations, separators=(",", ":")).encode()).hexdigest()
    if expected is not None and digest != expected:
        raise ValueError("exported equation hash mismatch")
    result = []
    for equation in equations:
        monomials = []
        for monomial in equation:
            if (monomial != sorted(set(monomial)) or
                    any(type(i) is not int or not 0 <= i < n for i in monomial)):
                raise ValueError("noncanonical exported monomial")
            monomials.append(sum(1 << i for i in monomial))
        result.append(tuple(monomials))
    return tuple(result)


def _verify_exported_root(mask, equations):
    if equations is None:
        return
    for monomials in equations:
        value = 0
        for monomial in monomials:
            value ^= int(mask & monomial == monomial)
        if value:
            raise ArithmeticError("original exported polynomial rejected a retained root")


def solve_projected_chain(instance, *, projection="sources", max_seconds=60.0,
                          max_sat_calls=1_000_000, max_output=10_000,
                          max_gates=1_000_000, conflict_limit=None):
    """Enumerate all original finite-coordinate roots, or explicitly unresolved.

    Source projection uses one global SAT call per distinct source tuple plus
    the final UNSAT call. Each found source tuple's full (at most eight) fiber
    is exhausted exactly before one source blocking clause is added. The
    comparison mode blocks full original assignments using the same solver.
    Output caps retain exact counting only when global exhaustion completes.
    """
    started = time.perf_counter()
    if instance.get("summands", 5) != 5:
        raise ValueError("this finite-fiber integration currently requires five summands")
    if projection not in ("sources", "all_variables"):
        raise ValueError("unknown projection mode")
    if not isinstance(max_seconds, (int, float)) or max_seconds < 0:
        raise ValueError("max_seconds must be nonnegative")
    if type(max_sat_calls) is not int or max_sat_calls < 0 or type(max_output) is not int or max_output < 0:
        raise ValueError("count budgets must be nonnegative integers")
    compiled = compile_s3_chain(instance, max_gates=max_gates)
    field = compiled.circuit.field
    dimension = instance["source_dimension"]
    source_bits = 5 * dimension
    basis = tuple(field.validate(_integer(v)) for v in instance.get("source_basis_hex", [1 << i for i in range(dimension)]))
    target = field.validate(_integer(instance["target"]["x_hex"]))
    curve_b = field.validate(_integer(instance.get("curve", {}).get("b_hex", "0x1")))
    if curve_b == 0:
        raise ValueError("nonzero curve constant is required for the bounded-fiber certificate")
    equations = _equation_masks(instance)
    for receipt in instance.get("native_coordinate_evaluation_checks", []):
        if compiled.evaluate(_integer(receipt["assignment_hex"])) != _integer(receipt["coordinate_output_hex"]):
            raise ValueError("circuit disagrees with native exported coordinate receipt")
    module = load_backend()
    solver = module.Solver(threads=1)
    # Declare even unconstrained original coordinates to get a full model.
    solver.add_clause([compiled.primary_count, -compiled.primary_count])
    for clause in compiled.circuit.cnf_clauses:
        solver.add_clause(list(clause))
    for variables, parity in compiled.circuit.xor_rows:
        solver.add_xor_clause(list(variables), bool(parity))
    quadratic = BinaryQuadraticFiber(field)
    stats = {**compiled.statistics(), "projection": projection,
             "sat_calls": 0, "source_tuples": 0, "fiber_root_calls": 0,
             "recovered_assignments": 0, "output_cap": max_output,
             "compile_and_setup_seconds": time.perf_counter() - started}
    retained, count, complete, reason = [], 0, False, None
    block_width = source_bits if projection == "sources" else compiled.primary_count
    while True:
        remaining = max_seconds - (time.perf_counter() - started)
        if remaining <= 0:
            reason = "charged pipeline time budget exhausted"
            break
        if stats["sat_calls"] >= max_sat_calls:
            reason = "SAT call budget exhausted"
            break
        options = {"time_limit": remaining}
        if conflict_limit is not None:
            options["confl_limit"] = conflict_limit
        verdict, model = solver.solve(**options)
        stats["sat_calls"] += 1
        if verdict is False:
            complete = True
            break
        if verdict is not True or model is None:
            reason = "native SAT search returned unknown within its budget"
            break
        if len(model) <= compiled.primary_count:
            raise ArithmeticError("SAT backend did not return all original coordinates")
        original = sum(1 << (i - 1) for i in range(1, compiled.primary_count + 1) if model[i])
        sources = _decode_sources(original, dimension, basis)
        auxiliary = _decode_auxiliary(original, field.degree, source_bits)
        _check_links(field, sources, auxiliary, target, curve_b)
        source_mask = original & ((1 << source_bits) - 1)
        if projection == "sources":
            fiber = five_source_intermediate_fiber(field, sources, target, curve_b=curve_b,
                                                   quadratic_solver=quadratic)
            chains = fiber["intermediate_chains"]
            if auxiliary not in chains:
                raise ArithmeticError("exact fiber omitted the native SAT witness")
            stats["source_tuples"] += 1
            stats["fiber_root_calls"] += fiber["stats"]["quadratic_root_calls"]
            roots = [source_mask | sum(u << (source_bits + i * field.degree)
                                     for i, u in enumerate(chain)) for chain in chains]
        else:
            roots = [original]
        for root in roots:
            count += 1
            if len(retained) < max_output:
                _verify_exported_root(root, equations)
                retained.append(root)
        solver.add_clause([-i if original & (1 << (i - 1)) else i
                           for i in range(1, block_width + 1)])
    stats["recovered_assignments"] = count
    stats["pipeline_seconds"] = time.perf_counter() - started
    output_complete = complete and len(retained) == count
    return {"schema": "projected-field-circuit-root-result-v1",
            "algorithm": "shared_field_circuit_projected_sat_fibers",
            "backend": "pycryptosat 5.11.21 / CryptoMiniSat, one thread",
            "scope": "all original Boolean roots of the parameterized finite-coordinate S3 chain; no group-decomposition claim",
            "status": "complete" if output_complete else "output_truncated" if complete else "unresolved",
            "search_complete": complete, "output_complete": output_complete,
            "solution_count_exact": count if complete else None,
            "solution_count_lower_bound": count, "assignments_hex": [hex(v) for v in sorted(retained)],
            "n_vars": compiled.primary_count, "source_bits": source_bits,
            "verified": True, "unresolved_reason": reason, "stats": stats,
            "certificate": "UNSAT after exact source-fiber exhaustion" if complete and projection == "sources"
                           else "UNSAT after full original-coordinate model enumeration" if complete
                           else "global exhaustion not established"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("instance")
    parser.add_argument("--projection", choices=("sources", "all_variables"), default="sources")
    parser.add_argument("--max-seconds", type=float, default=60)
    parser.add_argument("--max-output", type=int, default=10_000)
    parser.add_argument("--max-sat-calls", type=int, default=1_000_000)
    parser.add_argument("--output")
    args = parser.parse_args()
    problem = json.loads(Path(args.instance).read_text())
    result = solve_projected_chain(problem, projection=args.projection,
                                   max_seconds=args.max_seconds, max_output=args.max_output,
                                   max_sat_calls=args.max_sat_calls)
    result["input_sha256"] = hashlib.sha256(Path(args.instance).read_bytes()).hexdigest()
    result["instance_id"] = problem.get("instance_id")
    data = json.dumps(result, indent=2) + "\n"
    if args.output:
        Path(args.output).write_text(data)
        print(json.dumps({"status": result["status"], "count": result["solution_count_exact"],
                          "output": args.output, "sat_calls": result["stats"]["sat_calls"]}))
    else:
        print(data)
    return 0 if result["output_complete"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
