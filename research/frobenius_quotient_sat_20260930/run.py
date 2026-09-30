#!/usr/bin/env python3
"""Run the target-only decomposition grid.  One fresh single-thread
CryptoMiniSat solver per attempt; wall time includes model construction,
loading, solving and verification.  Writes one JSON line per attempt."""
import argparse, json, multiprocessing as mp, random, time

import pycryptosat

from encode import Curve, build, decode, verify

FORMULATIONS = {
    "specialized_xor": ("specialized", False, True),
    "specialized_cnf": ("specialized", False, False),
    "subgroup_xor": ("specialized", True, True),
    "line_xor": ("line", False, True),
    "line_subgroup": ("line", True, True),
    "line_subgroup_cnf": ("line", True, False),
}


def factor_space(E, s, rng):
    """Seeded random s-dim subspace of ker(Tr) holding >= 4 valid x."""
    f = E.f
    while True:
        basis = []
        while len(basis) < s:
            v = rng.getrandbits(f.n)
            if f.tr(v):
                continue
            span = {0}
            for b in basis:
                span |= {w ^ b for w in span}
            if v not in span:
                basis.append(v)
        span = {0}
        for b in basis:
            span |= {w ^ b for w in span}
        valid = sorted(x for x in span if x and (P := E.lift(x)) is not None
                       and E.mul(E.h, P) is None)
        if len(valid) >= 4:
            return basis, valid


def targets(E, valid, kind, k, rng):
    out = []
    while len(out) < k:
        if kind == "planted":
            xs = rng.sample(valid, 4)
            pts = [E.lift(x) for x in xs]
            pts = [E.neg(P) if rng.getrandbits(1) else P for P in pts]
            R = None
            for P in pts:
                R = E.add(R, P)
            witness = xs
        else:
            g = None
            while g is None:
                g = E.mul(4, E.random_point(rng))
            R, witness = E.mul(rng.randrange(1, E.h), g), None
        if R is not None:
            out.append((R, witness))
    return out


def attempt(job):
    n, basis, R, form, budget, witness = job
    final, sub, nx = FORMULATIONS[form]
    E = Curve(n)
    t0 = time.perf_counter()
    M, cvars = build(n, basis, R, final, sub, nx)
    solver = pycryptosat.Solver(threads=1, time_limit=budget)
    M.load(solver)
    t1 = time.perf_counter()
    # Each proposed tuple is checked on the curve; a rejected tuple is blocked
    # and the same solver continues inside the same time budget.
    rejected, found = 0, None
    while True:
        sat, sol = solver.solve()
        if not sat:
            status = "timeout" if sat is None else "unsat"
            break
        xs = decode(sol, cvars, basis)
        if verify(E, xs, R):
            status, found = "verified", sorted(xs)
            break
        rejected += 1
        solver.add_clause([-v if sol[v] else v for ci in cvars for v in ci])
    t2 = t3 = time.perf_counter()
    return {"n": n, "formulation": form, "budget_s": budget,
            "target": [hex(R[0]), hex(R[1])], "status": status,
            "found_x": [hex(x) for x in found] if found else None,
            "build_load_s": round(t1 - t0, 4), "solve_s": round(t2 - t1, 4),
            "total_s": round(t3 - t0, 4), "size": M.size(),
            "planted": witness is not None, "rejected_tuples": rejected}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--s", type=int, required=True)
    ap.add_argument("--seed", type=int, default=20260930)
    ap.add_argument("--targets", type=int, default=4)
    ap.add_argument("--kinds", default="planted,random")
    ap.add_argument("--forms", default=",".join(FORMULATIONS))
    ap.add_argument("--budget", type=float, default=2.0)
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = random.Random(a.seed * 1000 + a.n)
    E = Curve(a.n)
    basis, valid = factor_space(E, a.s, rng)
    jobs, meta = [], {"n": a.n, "s": a.s, "seed": a.seed, "odd_subgroup": E.h,
                      "basis": [hex(b) for b in basis],
                      "valid_x_in_space": len(valid), "targets": {}}
    for kind in a.kinds.split(","):
        ts = targets(E, valid, kind, a.targets, rng)
        meta["targets"][kind] = [[hex(R[0]), hex(R[1])] for R, _ in ts]
        for R, w in ts:
            for form in a.forms.split(","):
                jobs.append((a.n, basis, R, form, a.budget, w))
    with open(a.out, "w") as fh:
        fh.write(json.dumps({"meta": meta}) + "\n")
        with mp.Pool(a.procs) as pool:
            for r in pool.imap(attempt, jobs):
                fh.write(json.dumps(r) + "\n")
                fh.flush()
                print(r["formulation"], r["planted"], r["status"], r["total_s"], flush=True)


if __name__ == "__main__":
    main()
