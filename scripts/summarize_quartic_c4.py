#!/usr/bin/env python3
"""Section 11.17's table: C4, the cost of one S5 solve over F_{p^4}, and the
k = 4 crossover it implies, from the frozen runs.

`research/notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md` §11.16 registered
this measurement before any F_{p^4} code existed: the asymptotic constant
r_inf = c_LA / c_rho is derived and does not depend on the solve; C4 sets
only where the relation phase hands over, at p* = 0.0952 C4 / (1 - r_inf),
and the charpoly of the 4,096-dimensional multiplication matrix bounds C4
below at 7.1e10.  This prints the measured C4 per residual, its phases, the
cross-check against the meet-in-the-middle oracle, and n* = p*^4 at both
ends of §11.16's r_inf range.  n* is an extrapolation on the exponents
n^(1/4) (relations) and n^(1/2) (linear algebra and rho), and is printed as
one.

Usage:
    python3 scripts/summarize_quartic_c4.py C4.json [GENERIC.json]
"""

from __future__ import annotations

import json
import math
import sys

FLOOR = 1.0391 * 4096**3  # section 11.16: the k = 3 charpoly coefficient at D = 4,096
C_ADD = 97.0  # F_p multiplications per F_{p^4} addition (derived; measured below)
S_RHO = 1.3
R_INF = {"phi = 0.73": 0.338, "phi = 1.0": 0.634}


def main() -> None:
    rep = json.load(open(sys.argv[1]))
    rows = rep["residuals"]
    print(f"p = {rep['p']}, seed {rep['seed']}, base {rep['base']} points, "
          f"measured {rep['fp_muls_per_add']:.1f} F_p muls per addition, "
          f"symmetrised S5 once per curve: {rep['precompute_muls']:.3e} muls")
    print()
    print("| residual | dim | rows x cols | C4 | echelon | normal forms | charpoly | eigenvectors "
          "| rational eigenvalues | solver = oracle | planted found | wall |")
    print("|---|---:|---|---:|---:|---:|---:|---:|---:|---|---|---:|")
    c4 = []
    for r in rows:
        s = r["solve"]
        assert s["outcome"] == "solved", s["outcome"]
        assert r["agree"], r
        c4.append(s["total_muls"])
        tag = "constructed" if r["constructed"] else "random"
        print(f"| {r['index']} ({tag}) | {s['dim']} | {s['rows']:,} x {s['cols']:,} "
              f"| {s['total_muls']:.4e} | {s['echelon_muls']:.3e} | {s['nf_muls']:.3e} "
              f"| {s['charpoly_muls']:.3e} | {s['eigenvector_muls']:.3e} "
              f"| {s['rational_eigenvalues']} | {r['solver_quadruples']} = {r['oracle_quadruples']} "
              f"| {r['planted_found']} | {s['total_ms'] / 1e3:.0f} s |")
    mean = sum(c4) / len(c4)
    spread = (max(c4) - min(c4)) / mean
    print()
    print(f"C4 mean {mean:.4e} over {len(c4)} residuals (spread {100 * spread:.2f} %), "
          f"{mean / FLOOR:.1f}x the registered floor {FLOOR:.3e}")
    a = 12.0 / (S_RHO * C_ADD)
    for name, r_inf in R_INF.items():
        pstar = a * mean / (1.0 - r_inf)
        print(f"  {name}: r_inf {r_inf}, p* = {pstar:.3e} = 2^{math.log2(pstar):.2f}, "
              f"n* = 2^{4 * math.log2(pstar):.1f} (extrapolated); "
              f"S/rho at p = {rep['p']}: {a * mean / rep['p'] + r_inf:.3e}")
    if len(sys.argv) > 2:
        gen = json.load(open(sys.argv[2]))
        xs = [math.log(g["solve"]["dim"]) for g in gen if g["solve"]["dim"] >= 256]
        ys = [math.log(g["solve"]["total_muls"]) for g in gen if g["solve"]["dim"] >= 256]
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
        pred = math.exp(my + b * (math.log(4096) - mx))
        print()
        print(f"generic planted systems, d = 4..{max(g['equation_degree'] for g in gen)}: "
              f"C ~ dim^{b:.3f}, predicting {pred:.3e} at dim 4,096 "
              f"({100 * (mean / pred - 1):+.1f} % from the S5 measurement)")


if __name__ == "__main__":
    main()
