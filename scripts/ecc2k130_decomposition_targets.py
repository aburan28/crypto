#!/usr/bin/env python3
"""Derived boundaries for the follow-on experiments to the ECC2K-130 decomposition note.

`RESEARCH_ECC2K130_DECOMPOSITION.md` prices one family of decomposition methods
and closes it.  This script derives the boundaries the *next* experiments have
to be scored against, so each one is pre-registered with a number rather than a
hope:

  E1  a scale-model ladder -- prime `n` with `2` primitive, where no invariant
      subspace exists, exactly as at 131.  Predicts `Lambda = ops / 2^n` flat.
  E2  the `l`-flatness sweep at fixed `n`.  Predicts flat then rising.
  E3  the two hypothetical-oracle floors: a free decomposition *detector* that
      only answers for the whole base, against one that answers for sub-bases
      and so localises a witness by bisection.  The gap between them is where
      the curve's protection actually lives.
  E4  single large primes, guarded so the oracle still filters.
  E5  the yield distribution, not just its mean: sample sizes for a Poisson test.
  E6  Frobenius-stable orbit-union bases, on the ladder.

    python3 scripts/ecc2k130_decomposition_targets.py

Writes `experiments/ecc2k130_decomposition_targets.json`.
"""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "experiments/ecc2k130_decomposition_targets.json"

_SPEC = importlib.util.spec_from_file_location(
    "ecc2k130_decomp", Path(__file__).with_name("ecc2k130_point_decomposition.py"))
assert _SPEC and _SPEC.loader
DEC = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(DEC)

N131 = DEC.N131
lc, la = DEC.log2_comb, DEC.log2_add


def sweep(f, lo: float = 4.0, hi: float = 131.0, step: float = 0.05):
    """Minimise `f(l)[0]` over a grid of subspace dimensions."""
    best = None
    x = lo
    while x <= hi + 1e-9:
        got = f(x)
        if best is None or got[0] < best[0]:
            best = got + (round(x, 2),)
        x += step
    return best


# ── E1 / E6: the scale-model ladder ──────────────────────────────────────

def small_factor_profile(n: int, bound: int = 10 ** 7):
    """Trial-divide the odd part of `#E(F_2^n)`; report the largest prime found."""
    order = DEC.curve_order(n)
    v = (order & -order).bit_length() - 1
    rest, facs = order >> v, []
    d = 3
    while d * d <= rest and d < bound:
        while rest % d == 0:
            facs.append(d)
            rest //= d
        d += 2
    facs.append(rest)
    return order, 1 << v, sorted(facs)


def ladder(max_n: int = 70):
    """Prime `n` where `2` is a primitive root: cosets `{1, n-1}`, so the only
    invariant subspace dimensions are `0, 1, n-1, n` -- the structure at 131."""
    rows = []
    for n in range(11, max_n + 1):
        if not (DEC.is_prime(n) and DEC.mult_order(2, n) == n - 1):
            continue
        order, cof, facs = small_factor_profile(n)
        big = max(facs)
        l3 = DEC.saturating_l(3, n)
        cell = DEC.cost_cell(3, l3, n, 1)
        orbit = DEC.cost_cell(3, l3, n, 1, frobenius=True)
        rows.append({
            "n": n,
            "curve_order": str(order), "cofactor": cof,
            "log2_largest_prime_subgroup": round(math.log2(big), 2),
            "subgroup_is_faithful": math.log2(big) >= n - 4,
            "available_invariant_dimensions": DEC.available_subspace_dimensions(n),
            "saturating_l_m3": round(l3, 2),
            "predicted_log2_total_m3": cell["log2_total"],
            "predicted_lambda_m3": round(2 ** (cell["log2_total"] - n), 3),
            "predicted_log2_S_m3": round(cell["log2_total"] - math.log2(big) / 2, 2),
            "predicted_log2_total_orbit_base": orbit["log2_total"],
            "predicted_orbit_saving_log2": round(cell["log2_total"]
                                                 - orbit["log2_total"], 2),
            "end_to_end_feasible": cell["log2_total"] <= 40.0,
        })
    return rows


# ── E3: what a hypothetical oracle would have to be ──────────────────────

def detector_floor(m: int, n: int, localising: bool, frobenius: bool = False):
    """Floor for a free decomposition oracle of one of two strengths.

    *Localising* means the oracle answers "is there a decomposition with every
    summand drawn from `W`?" for any sub-base `W`, so a witness falls out of
    `O(m log |F|)` free queries and only the targets and the linear algebra are
    charged.  *Non-localising* means it answers only for the whole base, so every
    target the detector passes still costs a full `C(|F|, m-1)` search to turn
    into a relation.  Nothing below either line is reachable however good the
    algebra gets; the distance between them is the value of localisation.
    """
    collapse = math.log2(n) if frobenius else 0.0

    def f(l):
        targets = (l - collapse) + max(0.0, n - lc(l, m))
        linalg = math.log2(m) + 2 * (l - collapse)
        parts = [targets, linalg]
        if not localising:
            parts.append((l - collapse) + lc(l, m - 1))     # witness by search
        return (la(*parts), targets, linalg)

    tot, targets, linalg, l = sweep(f)
    return {"m": m, "localising": localising, "frobenius_stable": frobenius,
            "l_star": l, "log2_floor": round(tot, 2),
            "log2_targets": round(targets, 2),
            "log2_linear_algebra": round(linalg, 2)}


# ── E4: a single large prime, guarded so the oracle still filters ────────

def large_prime(m: int, n: int, log2_memory: float | None = None):
    """Enumerate `m-1` base points, accept the last summand anywhere in a larger
    subspace `V'`; pair partial relations off by their large prime.

    Guarded by `yield <= 1 partial relation per target`.  Without that guard `V'`
    grows to the whole field, the oracle stops filtering anything, and what is
    left is a birthday search on differences of targets -- a generic algorithm,
    not an index calculus, and it has to be scored as one.
    """
    best = None
    ls = 4.0
    while ls <= n:
        lp = ls
        while lp <= n:
            y = lc(ls, m - 1) + lp - math.log2(m) - n      # log2 partials/target
            if y > 0:
                break                                      # guard
            store = max(ls, (ls + lp + 1) / 2)             # birthday on 2^lp primes
            if log2_memory is not None and store > log2_memory:
                lp += 0.05
                continue
            collection = (store - y) + lc(ls, m - 1)
            linalg = math.log2(m) + 2 * ls
            tot = la(collection, linalg)
            if best is None or tot < best[0]:
                best = (tot, ls, lp, store, collection, linalg)
            lp += 0.05
        ls += 0.05
    if best is None:
        return None
    tot, l, lp, store, collection, linalg = best
    return {"m": m, "l_star": round(l, 2), "large_prime_dim": round(lp, 2),
            "log2_partial_relation_store": round(store, 2),
            "log2_collection": round(collection, 2),
            "log2_linear_algebra": round(linalg, 2),
            "log2_total": round(tot, 2),
            "log2_memory_cap": log2_memory}


# ── E5: the yield distribution, not just its mean ────────────────────────

def poisson_test_design(effect: float = 0.20, z: float = 3.0):
    """Targets needed to detect a given departure from the Poisson yield model.

    The statistic is the index of dispersion `Var/mean` of the per-target
    decomposition count, which is `1` under the model.  Its standard error on `T`
    samples is `sqrt(2/T)`, so resolving a relative effect `e` at `z` sigma needs
    `T = 2 (z/e)^2`.
    """
    t = math.ceil(2 * (z / effect) ** 2)
    return {"statistic": "index of dispersion Var/mean of the per-target count",
            "null_value": 1.0, "effect_resolved": effect, "sigmas": z,
            "targets_per_cell": t,
            "note": "the note's measured rate ratios ran 0.95 to 1.38, all but one "
                    "above 1, which is the direction of under-dispersion; that is "
                    "what this is sized to confirm or kill"}


def main() -> None:
    order = DEC.curve_order(N131)
    r = order // 4
    log2_rho = math.log2(math.sqrt(math.pi * r / (2 * 2 * N131)))

    rungs = ladder()
    e1 = {
        "question": "does the product law hold on curves structurally identical to "
                    "ECC2K-130 -- prime n, 2 primitive, no invariant subspace?",
        "primary_metric": "Lambda = total operations / 2^n, predicted flat at m",
        "secondary_metric": "log2 S = log2(total) - log2(subgroup)/2, predicted "
                            "slope 1/2 against n",
        "falsifier": "a least-squares slope of log2(total) against n below 0.95, "
                     "or Lambda below 0.5*m, on four or more rungs",
        "rungs": rungs,
        "end_to_end_rungs": [c["n"] for c in rungs if c["end_to_end_feasible"]],
        "composed_rungs": [c["n"] for c in rungs if not c["end_to_end_feasible"]],
        "faithful_subgroup_rungs": [c["n"] for c in rungs if c["subgroup_is_faithful"]],
    }

    e2 = {
        "question": "is the total really flat in the factor-base dimension?",
        "cells": [DEC.cost_cell(3, l, 59, 1) for l in (8, 12, 16, 20, 24, 28)],
        "saturating_l": round(DEC.saturating_l(3, 59), 2),
        "falsifier": "any dimension whose measured total is below half the flat line",
    }

    floors = []
    for m in range(2, 9):
        for loc in (True, False):
            for frob in (False, True):
                floors.append(detector_floor(m, N131, loc, frob))
    gaps = []
    for m in range(2, 9):
        a = next(f for f in floors if f["m"] == m and f["localising"]
                 and not f["frobenius_stable"])
        b = next(f for f in floors if f["m"] == m and not f["localising"]
                 and not f["frobenius_stable"])
        gaps.append({"m": m,
                     "log2_localising_floor": a["log2_floor"],
                     "log2_detector_only_floor": b["log2_floor"],
                     "log2_gap": round(b["log2_floor"] - a["log2_floor"], 2),
                     "detector_only_beats_rho": b["log2_floor"] < log2_rho,
                     "localising_beats_rho": a["log2_floor"] < log2_rho})
    e3 = {
        "question": "how much of the difficulty is deciding, and how much is "
                    "localising the witness?",
        "floors": floors, "gaps": gaps,
        "min_detector_only_floor": min(g["log2_detector_only_floor"] for g in gaps),
        "min_localising_floor": min(g["log2_localising_floor"] for g in gaps),
        "falsifier": "a real oracle whose sub-base query costs less than the same "
                     "query on the full base by more than the sub-base ratio",
    }

    e4 = {
        "question": "do large primes move the product law, or only relabel it?",
        "unbounded_memory": [large_prime(m, N131) for m in (2, 3, 4)],
        "memory_capped": [large_prime(3, N131, mu) for mu in (30, 40, 50, 60)],
        "generic_line": [DEC.generic_with_memory(mu, math.log2(r))
                         for mu in (30, 40, 50, 60)],
        "falsifier": "a guarded cell below the BSGS line at the same memory",
        "note": "the model is crude on purpose -- collision bookkeeping and "
                "relation independence are not priced -- so these are the numbers "
                "the experiment must beat or reproduce, not results",
    }

    e5 = poisson_test_design()

    e6 = {
        "question": "does the Frobenius collapse really arrive from an orbit union, "
                    "with no subspace anywhere?",
        "closed_form": "collection = m 2^n / n",
        "predicted_saving_log2_at_131": round(math.log2(N131), 2),
        "rows_at_131": [DEC.cost_cell(m, DEC.saturating_l(m, N131), N131, 1,
                                      frobenius=True) for m in range(2, 9)],
        "falsifier": "a measured saving on the ladder differing from n by more "
                     "than 20%, or orbit unions of size below n at prime n",
    }

    report = {
        "schema": "ecc2k130_decomposition_targets/v1",
        "purpose": "pre-registered boundaries for the follow-on experiments to "
                   "RESEARCH_ECC2K130_DECOMPOSITION.md",
        "unit": "log2 group operations",
        "log2_rho_reference": round(log2_rho, 4),
        "E1_scale_model_ladder": e1,
        "E2_dimension_flatness": e2,
        "E3_detector_versus_localisation": e3,
        "E4_large_primes": e4,
        "E5_yield_distribution": e5,
        "E6_orbit_union_bases": e6,
    }
    OUT.write_text(json.dumps(report, indent=2) + "\n")

    print(f"rho reference 2^{log2_rho:.4f}")
    print(f"E1 ladder: n = {[c['n'] for c in rungs]}")
    print(f"   end-to-end {e1['end_to_end_rungs']}, composed {e1['composed_rungs']}, "
          f"faithful subgroup {e1['faithful_subgroup_rungs']}")
    print(f"E3 detector-only floor bottoms at 2^{e3['min_detector_only_floor']:.2f}, "
          f"localising at 2^{e3['min_localising_floor']:.2f}, "
          f"rho at 2^{log2_rho:.2f}")
    for g in gaps:
        print(f"   m={g['m']}  detector-only 2^{g['log2_detector_only_floor']:6.2f}  "
              f"localising 2^{g['log2_localising_floor']:6.2f}  gap 2^{g['log2_gap']:.2f}"
              f"   beats rho: detector {g['detector_only_beats_rho']}, "
              f"localising {g['localising_beats_rho']}")
    for c in e4["unbounded_memory"]:
        if c:
            print(f"E4 m={c['m']} large prime: 2^{c['log2_total']:.2f} "
                  f"= 2^{c['log2_total'] - log2_rho:+.2f} x rho, "
                  f"store 2^{c['log2_partial_relation_store']:.2f}")
    print(f"E5 targets per cell: {e5['targets_per_cell']}")
    print(f"E6 orbit-union saving: 2^{e6['predicted_saving_log2_at_131']:.2f}")
    print(f"wrote {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
