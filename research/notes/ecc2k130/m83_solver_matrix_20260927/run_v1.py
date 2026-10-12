"""Frozen m83 fixed-phase decomposition benchmark; deliberately not a DLP benchmark."""
from __future__ import annotations

import argparse
import collections
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import random
import resource
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE / "source" / "archived"), str(HERE / "source" / "archived" / "vendor"),
                str(HERE / "source" / "f5b")]
import core as c
from scale import Setup

N, S, MODULUS = 83, 2, (1 << 83) | (1 << 45) | 7
SEEDS = (260938, 260939)
ENGINES = ("fes", "sat-z3", "f4-boolean", "f5b-boolean")
LIMIT = 5.0


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def evaleq(equation, assignment):
    return sum(assignment & mask == mask for mask in equation) & 1


def clock():
    return {"wall": time.perf_counter(), "cpu": time.process_time()}


def delta(before):
    after = clock()
    return {"wall_s": after["wall"] - before["wall"], "cpu_s": after["cpu"] - before["cpu"]}


def field_generator(curve):
    for x in itertools.count(1):
        points = curve.lift(x)
        if points:
            g = curve.scale(points[0], curve.cofactor)
            if g:
                assert curve.scale(g, curve.r) is None and curve.valid(g)
                return g


def sqrt_mod(a, prime):
    assert pow(a, (prime - 1) // 2, prime) == 1
    q, s = prime - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (prime - 1) // 2, prime) != prime - 1:
        z += 1
    m, cc, t, x = s, pow(z, q, prime), pow(a, q, prime), pow(a, (q + 1) // 2, prime)
    while t != 1:
        i, v = 1, t * t % prime
        while v != 1:
            v = v * v % prime
            i += 1
        b = pow(cc, 1 << (m - i - 1), prime)
        x, t, cc, m = x * b % prime, t * b * b % prime, b * b % prime, i
    assert x * x % prime == a % prime
    return x


def make_inputs(setup, seed, generator):
    curve = setup.curve
    used = []
    for payload in range(1, 1 << S):
        x = setup.normal.to_field(setup.space.encode(payload))
        for point in curve.lift(x):
            if curve.scale(point, curve.r) is None:
                used.append((payload, point))
    planted = []
    for picks in itertools.product(used, repeat=3):
        points = tuple(curve.frob(p, phase) for phase, (_, p) in enumerate(picks))
        if not c.proper(curve, points):
            continue
        target = c.group_sum(curve, points)
        middle = curve.add(points[0], points[1])
        if target is None or middle is None or not target[0] or not middle[0]:
            continue
        if target not in [q["target"] for q in planted]:
            planted.append({"kind": "planted", "index": len(planted), "target": target,
                            "payloads": [p for p, _ in picks], "witness": points})
        if len(planted) == 2:
            break
    rng = random.Random(seed * 1000 + 832)
    natural = [{"kind": "natural", "index": i, "scalar": (scalar := rng.randrange(1, curve.r)),
                "target": curve.scale(generator, scalar)} for i in range(4)]
    return planted, natural, used


def verify_roots(setup, maps, target, roots):
    verified, false_positive, examples = [], [], {}
    for bits in roots:
        xs = [c.evaluate(poly, bits) for poly in maps]
        witness = c.lift(setup.curve, xs, target)
        if witness is None:
            false_positive.append(bits)
        else:
            assert all(setup.curve.valid(p) and setup.curve.scale(p, setup.curve.r) is None
                       for p in witness)
            assert c.group_sum(setup.curve, witness) == target
            verified.append(bits)
            examples[str(bits)] = witness
    return {"verified": verified, "false_positive": false_positive, "witnesses": examples}


def f4_rows(equations, nvars, deadline):
    """Degree-batched squarefree Macaulay matrices, complete through degree nvars."""
    counters = collections.Counter()
    pivots = {}
    seen = set()
    for degree in range(nvars + 1):
        batch = []
        for eq in equations:
            for mask in range(1 << nvars):
                if time.perf_counter() > deadline:
                    raise TimeoutError("F4 Macaulay construction deadline")
                row = c.times_monomial(eq, mask)
                if row and max(m.bit_count() for m in row) == degree:
                    bits = sum(1 << m for m in row)
                    if bits not in seen:
                        seen.add(bits)
                        batch.append(bits)
        counters["batches"] += 1
        counters["rows"] += len(batch)
        for row in batch:
            while row:
                lead = row.bit_length() - 1
                if lead not in pivots:
                    pivots[lead] = row
                    counters["independent_rows"] += 1
                    break
                row ^= pivots[lead]
                counters["row_xors"] += 1
    return [set(i for i in range(1 << nvars) if row >> i & 1) for row in pivots.values()], dict(counters)


def solve(equations, engine):
    deadline = time.perf_counter() + LIMIT
    nvars = 3 * S
    stats = {}
    if engine == "fes":
        roots = [bits for bits in range(1 << nvars)
                 if all(evaleq(eq, bits) == 0 for eq in equations)]
        stats = {"assignments_tested": 1 << nvars}
    elif engine == "f4-boolean":
        basis, stats = f4_rows(equations, nvars, deadline)
        roots = [bits for bits in range(1 << nvars)
                 if all(evaleq(eq, bits) == 0 for eq in basis)]
        stats["root_extraction_assignments"] = 1 << nvars
    elif engine == "f5b-boolean":
        from boolean_f5b import BooleanF5B
        solver = BooleanF5B(nvars, timeout=LIMIT, max_pairs=250000)
        basis = solver.basis([solver.from_terms(eq) for eq in equations])
        roots = [bits for bits in range(1 << nvars)
                 if all(solver.evaluate(eq, bits) == 0 for eq in basis)]
        stats = dict(solver.stats, root_extraction_assignments=1 << nvars)
    elif engine == "sat-z3":
        import z3
        variables = [z3.Bool(f"v{i}") for i in range(nvars)]
        all_monomials = set().union(*equations)
        monomials = {mask: z3.And(*[variables[i] for i in range(nvars) if mask >> i & 1])
                     if mask else z3.BoolVal(True) for mask in all_monomials}
        solver = z3.Solver()
        solver.set(random_seed=13)
        for eq in equations:
            solver.add(z3.Not(c.xor_tree(z3, [monomials[mask] for mask in sorted(eq)])))
        roots = []
        while True:
            remaining = deadline - time.perf_counter()
            if remaining <= 0:
                raise TimeoutError(f"SAT enumeration deadline after {len(roots)} models")
            solver.set(timeout=max(1, int(remaining * 1000)))
            result = solver.check()
            if result == z3.unsat:
                break
            if result != z3.sat:
                raise TimeoutError(f"SAT {solver.reason_unknown()} after {len(roots)} models")
            model = solver.model()
            bits = sum(int(z3.is_true(model.eval(v, model_completion=True))) << i
                       for i, v in enumerate(variables))
            assert all(evaleq(eq, bits) == 0 for eq in equations)
            roots.append(bits)
            solver.add(z3.Or(*[v != z3.BoolVal(bool(bits >> i & 1))
                               for i, v in enumerate(variables)]))
        stats = {"z3_version": z3.get_version_string(), "models_enumerated": len(roots),
                 "statistics": {k: v for k, v in solver.statistics()}}
    else:
        raise ValueError(engine)
    roots = sorted(roots)
    assert all(all(evaleq(eq, bits) == 0 for eq in equations) for bits in roots)
    return roots, stats


def fixed_child(seed, index, engine):
    resource.setrlimit(resource.RLIMIT_AS, (1 << 30, 1 << 30))
    timings = {}
    start = clock()
    setup = Setup(N, S, seed, "ternary_inline", MODULUS)
    curve = setup.curve
    assert curve.order == 4 * 2417851639230796216685689 and curve.cofactor == 4
    g = field_generator(curve)
    root = sqrt_mod(-7 % curve.r, curve.r)
    lambdas = [((-1 + root) * pow(2, -1, curve.r)) % curve.r,
               ((-1 - root) * pow(2, -1, curve.r)) % curve.r]
    lam, = [v for v in lambdas if curve.scale(g, v) == curve.frob(g)]
    assert pow(lam, N, curve.r) == 1 and lam != 1
    planted, natural, used = make_inputs(setup, seed, g)
    inputs = (planted + [None] * (2 - len(planted))) + natural
    assert len(planted) >= 1
    if inputs[index] is None:
        raise ValueError(f"setup_failed: input {index} missing (only {len(planted)} planted)")
    record = inputs[index]
    target = record["target"]
    assert curve.scale(target, curve.r) is None and curve.valid(target)
    timings["preparation"] = delta(start)
    start = clock()
    maps = [setup.phase_maps[i][i] for i in range(3)]
    field_poly = c.Algebra(curve.f, setup.counts, time.perf_counter() + 25).s4(
        *maps, target[0])
    equations = c.boolean_equations(field_poly, N)
    if record["kind"] == "planted":
        bits = sum(payload << (i * S) for i, payload in enumerate(record["payloads"]))
        assert all(evaleq(eq, bits) == 0 for eq in equations)
        assert c.lift(curve, [c.evaluate(m, bits) for m in maps], target) is not None
    timings["equation_build"] = delta(start)
    start = clock()
    roots, statistics = solve(equations, engine)
    timings["solver_including_setup"] = delta(start)
    start = clock()
    outcome = verify_roots(setup, maps, target, roots)
    timings["lifting_and_verification"] = delta(start)
    return {"kind": "fixed_phase", "seed": seed, "index": index, "engine": engine,
            "status": "complete", "input": record, "input_sha256": digest(record),
            "equation_sha256": digest([sorted(eq) for eq in equations]),
            "curve": {"n": N, "modulus": hex(MODULUS), "order": curve.order,
                      "r": curve.r, "cofactor": curve.cofactor, "generator": g,
                      "frobenius_eigenvalue": lam, "normal_generator": setup.normal.beta,
                      "normal_basis": setup.normal.basis,
                      "factor_base_representatives": sorted(set(min(p, curve.neg(p)) for _, p in used)),
                      "phases": list(setup.phases)},
            "nvars": 3 * S, "field_polynomial_terms": len(field_poly),
            "boolean_equations": len(equations), "boolean_terms": sum(map(len, equations)),
            "algebraic_roots": roots, "group_verification": outcome,
            "statistics": statistics, "field_anf_counts": dict(setup.counts), "timing": timings,
            "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
            "python": sys.version, "platform": platform.platform()}


def parent(out):
    out.mkdir(parents=True, exist_ok=False)
    files = sorted((HERE / "source").rglob("*.py")) + [HERE / "run.py", HERE / "PROTOCOL.md"]
    manifest = {"protocol": "m83_solver_matrix_20260927", "seeds": SEEDS,
                "engines": ENGINES, "limits": {"virtual_bytes": 1 << 30,
                "child_wall_seconds": 70, "equation_seconds": 25, "solve_seconds": LIMIT},
                "source_sha256": {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in files}, "python": sys.version, "platform": platform.platform()}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    for seed in SEEDS:
        for index in range(6):
            for engine in ENGINES:
                key = f"fixed-s{seed}-i{index}-{engine}"
                cmd = [sys.executable, str(HERE / "run.py"), "--child", "--seed", str(seed),
                       "--index", str(index), "--engine", engine]
                run_child(out, key, cmd)
        for arm in ("ternary_inline", "implicit"):
            key = f"phase-s{seed}-{arm}"
            cmd = [sys.executable, str(HERE / "source" / "archived" / "scale.py"),
                   "--child", "--n", "83", "--s", "2", "--seed", str(seed),
                   "--arm", arm, "--budget", "25"]
            run_child(out, key, cmd)


def run_child(out, key, cmd):
    t = clock()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=70,
                                env=dict(os.environ, PYTHONHASHSEED="0"))
        record = {"returncode": result.returncode, "stdout": result.stdout,
                  "stderr": result.stderr, "timing": delta(t)}
    except subprocess.TimeoutExpired as error:
        record = {"returncode": None, "timeout": True,
                  "stdout": error.stdout.decode(errors="replace") if error.stdout else "",
                  "stderr": error.stderr.decode(errors="replace") if error.stderr else "",
                  "timing": delta(t)}
    (out / f"{key}.json").write_text(json.dumps(record, sort_keys=True, indent=2) + "\n")
    print(key, record["returncode"], round(record["timing"]["wall_s"], 2), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path)
    p.add_argument("--child", action="store_true")
    p.add_argument("--seed", type=int)
    p.add_argument("--index", type=int)
    p.add_argument("--engine", choices=ENGINES)
    a = p.parse_args()
    if a.child:
        try:
            print(json.dumps(fixed_child(a.seed, a.index, a.engine), sort_keys=True))
        except Exception as error:
            print(json.dumps({"status": "exception", "type": type(error).__name__,
                              "message": str(error), "seed": a.seed, "index": a.index,
                              "engine": a.engine}), flush=True)
            raise
    else:
        parent(a.out)
