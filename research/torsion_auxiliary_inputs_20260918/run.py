"""Driver: one JSON per (size, case) under results/, every row verified.

    python3 run.py --bits 24 28 32 36 40 --seeds 3

For each size it searches a prime-order curve whose p − 1 has a divisor
d ≈ p^{1/2} (the p − 1 case) and another whose p + 1 has a divisor
d ≈ p^{1/3} (the p + 1 case), plants α, builds the auxiliary inputs, and
runs every variant with a fresh counter.  The rho reference runs on the
same curve, subgroup and target.  Nothing is timed; everything is counted.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import time

import cheon
import ec
import rho

HERE = os.path.dirname(os.path.abspath(__file__))


def fresh(inst: ec.CurveInstance):
    ctr = ec.Counter()
    return inst.curve(ctr), ctr


def run_case(bits: int, case: str, seeds: int, rho_seeds: int, rho_cap_bits: int) -> dict:
    exponent = 0.5 if case == "p-1" else 1 / 3
    inst = ec.find_curve(bits, seed=1000 * bits + (0 if case == "p-1" else 1),
                         want=case, exponent=exponent)
    p = inst.p
    fac = inst.fac_p_minus_1 if case == "p-1" else inst.fac_p_plus_1
    d = ec.best_divisor(fac, p**exponent)
    out = {
        "bits": bits, "case": case, "curve": inst.to_json(), "d": d,
        "log_p_d": math.log(d, p), "sqrt_p": math.sqrt(p),
        "floor_dlpwai": math.sqrt(p / d),  # generic lower bound Ω(√(p/d))
        "rows": [],
    }
    rng = random.Random(7 * bits + len(case))
    for s in range(seeds):
        alpha = rng.randrange(2, p - 1)
        E0 = inst.curve()
        G = inst.G
        n_aux = 2 * d + 1 if case == "p+1" else d + 1
        # auxiliary inputs are given by the problem: not charged
        G_pows = [G]
        for _ in range(n_aux - 1):
            G_pows.append(E0.mul(alpha, G_pows[-1]))
        G1, Gd = G_pows[1], G_pows[d]

        def record(name: str, res: dict, ctr: ec.Counter, extra=None):
            snap = ctr.snapshot()
            row = {
                "variant": name, "seed": s, "alpha": alpha,
                "correct": res.get("alpha") == alpha,
                "ops": snap["ops"], "add": snap["add"], "dbl": snap["dbl"],
                "smul": snap["smul"], "phases": snap["phases"],
                "S": snap["ops"] / math.sqrt(p),
                "ops_over_floor": snap["ops"] / math.sqrt(p / d),
            }
            if extra:
                row.update(extra)
            out["rows"].append(row)
            print(f"  {bits:2d}b {case} seed {s} {name:28s} ops={snap['ops']:>12,d} "
                  f"S={row['S']:8.3f} floor×={row['ops_over_floor']:8.2f} "
                  f"{'OK' if row['correct'] else 'WRONG'}", flush=True)

        if case == "p-1":
            E, ctr = fresh(inst)
            t = time.time()
            r = cheon.cheon_p_minus_1_bsgs(E, G, G1, Gd, p, d, comb=False)
            record("cheon_p-1_bsgs_naive", r, ctr, {"secs": time.time() - t})

            E, ctr = fresh(inst)
            t = time.time()
            r = cheon.cheon_p_minus_1_bsgs(E, G, G1, Gd, p, d, comb=True)
            record("cheon_p-1_bsgs_comb", r, ctr,
                   {"secs": time.time() - t, "table_points": r.get("table_points"),
                    "comb_w": r.get("comb_w")})

            E, ctr = fresh(inst)
            t = time.time()
            r = cheon.cheon_p_minus_1_kangaroo(E, G, G1, Gd, p, d, seed=s)
            record("cheon_p-1_kangaroo_comb", r, ctr, {"secs": time.time() - t})

            if bits <= 32:
                E, ctr = fresh(inst)
                t = time.time()
                r = cheon.plain_bsgs(E, G, G1, p, d)
                record("plain_bsgs_as_in_rust", r, ctr, {"secs": time.time() - t})
        else:
            E, ctr = fresh(inst)
            t = time.time()
            r = cheon.cheon_p_plus_1_bsgs(E, G_pows, p, d)
            record("cheon_p+1_bsgs_comb", r, ctr, {"secs": time.time() - t})

        if s < rho_seeds and bits <= rho_cap_bits:
            E, ctr = fresh(inst)
            t = time.time()
            r = rho.pollard_rho(E, G, G1, p, seed=100 + s)
            record("rho_reference", r, ctr,
                   {"secs": time.time() - t, "steps": r.get("steps"), "walks": r.get("walks")})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bits", type=int, nargs="+", default=[24, 28, 32, 36, 40])
    ap.add_argument("--cases", nargs="+", default=["p-1", "p+1"])
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--rho-seeds", type=int, default=3)
    ap.add_argument("--rho-cap-bits", type=int, default=40)
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for bits in args.bits:
        for case in args.cases:
            res = run_case(bits, case, args.seeds, args.rho_seeds, args.rho_cap_bits)
            path = os.path.join(args.out, f"{case.replace('+', 'plus').replace('-', 'minus')}_{bits:02d}.json")
            with open(path, "w") as fh:
                json.dump(res, fh, indent=1)
            print("wrote", path, flush=True)


if __name__ == "__main__":
    main()
