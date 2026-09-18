#!/usr/bin/env python3
"""The cost of a *logarithm*, not of a relation, with the Frobenius quotient.

Section 1 of the research note and `orbit_boundary.py` both price homogeneous
`m`-point relations among base points.  Section 3 says plainly that such
relations determine nothing: they give linear dependencies between unknown
logarithms, and `log_P(Q)` is not among them.  Priced that way the table is
answering the wrong question, and it answers it too favourably -- extended
past `m = 8` the homogeneous accounting crosses *below* the rho reference at
`m = 18` and reaches `2^-4.05` at `m = 28`, which is not a break of anything.

This module prices what a logarithm actually costs:

  * the unknowns are the `T` `sigma`-orbits of the support plus `d`, so
    `T + 1` relations are needed;
  * each relation must involve a **known** target `[a]P + [b]Q`, decomposed
    into `n` signed base points;
  * the cost of one relation is one meet-in-the-middle per target attempt,
    divided by the probability that a target decomposes at all.

**The measured input.**  That probability is

    P(a random target decomposes into n signed base points)
        = 2^(n-1) C(B,n) * |E[4] cap H| / |H|  =  2^n C(B,n) / r

and it is measured, not assumed -- see `validate_target_model.py` and
`results/target_decomposition.json`.

**The correction the measurement forced.**  An earlier form of this model
gave the `sigma` quotient a factor of 131 in that probability, on the
grounds that a decomposition of `sigma^k(target)` is as useful as one of
`target`, so there are 131 acceptable right-hand sides.  Measured, the
model was 20x optimistic.  The reason is structural: the support is
`sigma`-stable, so the **set of n-subset sums is itself `sigma`-closed** --
verified directly at `m = 17, 19`.  Accepting `sigma^k(target)` is
therefore the same chance 131 times, not 131 independent chances.  The
quotient buys memory and unknowns; it does not buy hit rate.

What the quotient does buy:

  * the stored side of the meet-in-the-middle holds one canonical class per
    `sigma`-orbit, so 131x fewer entries to build and to keep;
  * the support carries `B/131` unknowns rather than `B`, so `B/131`
    relations are needed rather than `B`.  This is the large saving.

**The stored side is built once, not once per target.**  An earlier form of
this module charged the meet-in-the-middle table build on every target
attempt.  That is wrong: the table holds `s`-subset sums of the support and
does not depend on the target at all, so it is built once and streamed
against for every attempt.  Charging it per attempt over-counted by
`2^7.3`, and the corrected optimum is `2^64.64` rather than `2^71.94`.

The optimum that follows is lopsided: store `n-1` points, stream **one**.
The support is a single Frobenius orbit -- 131 points, two unknowns, its own
orbit logarithm and `d` -- which needs no weight-two construction at all,
just a random curve point and its conjugates.

Two costs this leaves standing, both flagged rather than hidden.  The table
build assumes one representative per `sigma`-class can be enumerated in
constant amortised time (necklace enumeration over `Z/131`); enumerating
naively by fixing an element costs `2^3.3` more.  And the build is
`2^63.98` entries, so the total is a table-construction cost, not a search.

**Memory is charged at zero here, deliberately.**  The totals below assume
storage is free and instantaneous.  The conclusion does not depend on the
storage wall, which is the objection an engineer could always call a matter
of budget; it depends only on operation counts.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

RHO_LOG2 = 60.809
LOG2_R = 129.0
M = 131
LOG2_M = math.log2(M)
KEY_BYTES = 8

# Measured shortfall of the decomposition model against exhaustive counts:
# mean measured/predicted 0.84 over seven unsaturated small-curve cells
# (results/target_decomposition.json). Subset sums collide, so the true rate
# is slightly below 2^n C(B,n)/r and the true cost slightly above. Applied
# here so the total is not quoted optimistically.
MEASURED_RATE_FACTOR = 0.84


def _lf(n: int) -> float:
    return math.log2(math.factorial(n))


def _lchoose(B: int, k: int):
    if k < 0 or k > B:
        return None
    return sum(math.log2(B - i) for i in range(k)) - _lf(k)


def _lse(a: float, b: float) -> float:
    """log2(2^a + 2^b)."""
    hi, lo = max(a, b), min(a, b)
    return hi + math.log2(1 + 2 ** (lo - hi))


def cost(T: int, n: int, s: int, amortise: bool = True):
    """log2 cost of a logarithm from a `T`-orbit support and `n`-point relations.

    With `amortise` (the default and the correct accounting) the stored side
    is built once and the streamed side is paid per target attempt.  With
    `amortise=False` the whole meet-in-the-middle is charged per attempt,
    which is what this module did before and is 2^7.3 too expensive.
    """
    B = M * T
    a, b = _lchoose(B, s), _lchoose(B, n - s)
    if a is None or b is None:
        return None
    stored = s + a - LOG2_M           # canonical sigma-classes only
    streamed = (n - s) + b
    log_p = min(0.0, n + _lchoose(B, n) - LOG2_R
                + math.log2(MEASURED_RATE_FACTOR))   # measured, no sigma bonus
    log_attempts = math.log2(T + 1) - log_p
    if amortise:
        per_attempt = streamed
        total = _lse(stored, log_attempts + streamed)
    else:
        per_attempt = max(stored, streamed)
        total = log_attempts + per_attempt
    return {
        "orbits": T,
        "support_size": B,
        "relation_length": n,
        "mitm_split": f"{s}+{n - s}",
        "log2_stored_classes": round(stored, 3),
        "log2_streamed": round(streamed, 3),
        "log2_build_once": round(stored, 3),
        "log2_cost_per_attempt": round(per_attempt, 3),
        "log2_decomposition_probability": round(log_p, 3),
        "log2_target_attempts": round(log_attempts, 3),
        "relations_needed": T + 1,
        "log2_total_cost": round(total, 3),
        "log2_total_vs_rho": round(total - RHO_LOG2, 3),
        "log2_memory_bytes": round(stored + math.log2(KEY_BYTES), 3),
        "memory_exabytes": float(f"{2 ** (stored + math.log2(KEY_BYTES)) / 1e18:.4g}"),
    }


def optimise(max_orbits: int = 120000, max_n: int = 48, amortise: bool = True):
    best = None
    grid = list(range(1, 300)) + list(range(300, max_orbits + 1, 13))
    for T in grid:
        for n in range(2, max_n + 1):
            for s in range(1, n):
                r = cost(T, n, s, amortise=amortise)
                if r and (best is None or
                          r["log2_total_cost"] < best["log2_total_cost"]):
                    best = r
    return best


def homogeneous_extension():
    """The homogeneous accounting extended past m = 8, for the record.

    Included because it is the row an optimistic reading produces, and it
    goes below the reference. It prices relations that carry no information
    about `log_P(Q)`, so it is not a cost of a logarithm.
    """
    out = []
    for m in range(4, 29, 2):
        lB = (_lf(m) + LOG2_R - LOG2_M - (m - 1)) / (m - 1)
        if lB < LOG2_M:
            continue                    # B < 131: not a union of orbits
        best = None
        for s in range(1, m // 2 + 1):
            st = s * lB - _lf(s)
            sm = (m - s) * lB - _lf(m - s)
            c = max(st, sm)
            if best is None or c < best[1]:
                best = (s, c)
        s, c = best
        total = (c - LOG2_M) + (lB - LOG2_M)
        out.append({
            "m": m,
            "log2_support_size": round(lB, 3),
            "log2_total_cost": round(total, 3),
            "log2_total_vs_rho": round(total - RHO_LOG2, 3),
        })
    return out


def main():
    best = optimise()
    unamortised = optimise(amortise=False)
    homog = homogeneous_extension()
    crossing = next((r for r in homog if r["log2_total_vs_rho"] < 0), None)

    data = {
        "instance": "ECC2K-130",
        "log2_r": LOG2_R,
        "frobenius_order": M,
        "rho_reference_log2": RHO_LOG2,
        "rho_reference_source": (
            "experiments/ecc2k130_extension_field_boundary.json -> "
            "target.log2_rho_reference"),
        "memory_charged": False,
        "measured_rate_factor": MEASURED_RATE_FACTOR,
        "decomposition_model": "P = 2^n C(B,n) / r, measured in "
                               "results/target_decomposition.json",
        "best_logarithm_cost": best,
        "best_if_table_rebuilt_per_attempt": unamortised,
        "amortisation_saving_log2": round(
            unamortised["log2_total_cost"] - best["log2_total_cost"], 3),
        "homogeneous_extension_past_m8": homog,
        "homogeneous_first_sub_rho_m": crossing["m"] if crossing else None,
        "verdict": (
            f"Pricing a logarithm rather than a homogeneous relation, the best "
            f"configuration is a {best['orbits']}-orbit support "
            f"(B = {best['support_size']}) with {best['relation_length']}-point "
            f"relations, at 2^{best['log2_total_cost']} against the reference "
            f"2^{RHO_LOG2} -- short by "
            f"2^{best['log2_total_vs_rho']}, with memory charged at zero. The "
            f"homogeneous accounting extended past m = 8 does cross below the "
            f"reference, at m = "
            f"{crossing['m'] if crossing else 'n/a'}, but it prices relations "
            f"that carry no information about log_P(Q), which is the section 3 "
            f"error in its strongest form. Building the table once rather "
            f"than per attempt is worth "
            f"2^{unamortised['log2_total_cost'] - best['log2_total_cost']:.2f}."),
    }

    out = Path(__file__).resolve().parent / "results" / "target_boundary.json"
    out.write_text(json.dumps(data, indent=2))

    print("Homogeneous accounting, extended past m = 8 (WRONG QUESTION):")
    print(f"  {'m':>3} {'log2 B':>8} {'total':>8} {'vs rho':>9}")
    for r in homog:
        print(f"  {r['m']:>3} {r['log2_support_size']:>8.2f} "
              f"{r['log2_total_cost']:>8.2f} {r['log2_total_vs_rho']:>+9.2f}")
    print(f"\n  first sub-rho m: {crossing['m'] if crossing else 'none'}"
          "  <- prices relations that determine nothing\n")
    print("Cost of a LOGARITHM, quotient applied, memory free:")
    for k, v in best.items():
        print(f"  {k:34} {v}")
    print()
    print(data["verdict"])


if __name__ == "__main__":
    main()
