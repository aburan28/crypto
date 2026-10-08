#!/usr/bin/env python3
"""Derived boundaries for the follow-on experiments to the ECC2K-130 decomposition note.

`research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md` prices one family of decomposition methods
and closes it.  This script derives the boundaries the *next* experiments have
to be scored against, so each one is pre-registered with a number rather than a
hope:

  E1  a scale-model ladder -- prime `n` with `2` primitive, where no invariant
      subspace exists, exactly as at 131.  Predicts `Lambda = ops / 2^n` flat.
  E2  the `l`-flatness sweep at fixed `n`.  Predicts flat then rising.
  E3  what a free decomposition *detector* buys.  A detector that is an
      algorithm on the target localises its own witness by swapping candidate
      summands, so deciding and localising coincide; the second column prices
      the artificial alternative, a detector restricted to a fixed family of
      targets, and the gap measures the restriction, not the curve.
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

def detector_floor(m: int, n: int, mode: str, frobenius: bool = False,
                   k: int = 2):
    """Floor for a free decomposition oracle, by what it is allowed to be asked.

    `target_agnostic` -- the detector is an algorithm on the target's
    coordinates, so it may be called on `R - P + Q` as cheaply as on `R`.  That
    is enough to localise the witness with no sub-base query at all: walk the
    candidates `P`, swap in a class-matched base point `Q`, and keep the `P`
    whose query says yes `k` times running (`swap_localisation` in the companion
    script measures the two ways this misfires and shows both are collisions of
    size `Theta(m/|F|)`).  The queries are free by hypothesis; what is charged is
    the `k|F|` group operations that build them, once per *successful* target.

    `target_restricted` -- the detector answers only for a distinguished family
    of targets, so the swap is unavailable and the witness has to be extracted
    by meet-in-the-middle: a table of every `ceil(m/2)`-subset sum built once and
    shared, then `floor(m/2)`-subset probes per successful target.

    An earlier revision of this file called these two "non-localising" and
    "localising" and charged the first a naive `C(|F|, m-1)` search.  Both were
    wrong: the search is not the right price for the restricted oracle, and the
    restriction is not a property any proposed detector actually has.
    """
    assert mode in ("target_agnostic", "target_restricted"), mode
    collapse = math.log2(n) if frobenius else 0.0
    half_hi, half_lo = (m + 1) // 2, m // 2

    def f(l):
        targets = (l - collapse) + max(0.0, n - lc(l, m))
        linalg = math.log2(m) + 2 * (l - collapse)
        swap = math.log2(k) + (l - collapse) + l         # k|F| ops a relation
        table = lc(l, half_hi)                           # built once, shared
        witness = (l - collapse) + lc(l, half_lo)        # probes, per relation
        mitm = la(table, witness)
        # an agnostic detector may still build the table, so it takes whichever
        # of the two witness paths is cheaper; a restricted one has only the table
        use_swap = mode == "target_agnostic" and swap <= mitm
        parts = [targets, linalg] + ([swap] if use_swap else [table, witness])
        return (la(*parts), targets, linalg, use_swap, swap, table, witness)

    tot, targets, linalg, use_swap, swap, table, witness, l = sweep(f)
    out = {"m": m, "mode": mode, "frobenius_stable": frobenius,
           "l_star": l, "log2_floor": round(tot, 2),
           "log2_targets": round(targets, 2),
           "log2_linear_algebra": round(linalg, 2),
           "witness_route": "swap" if use_swap else "meet_in_the_middle"}
    if use_swap:
        out["queries_per_candidate"] = k
        out["log2_swap_localisation"] = round(swap, 2)
        out["log2_table_entries"] = None
    else:
        out["witness_split"] = [half_hi, half_lo]
        out["log2_table_entries"] = round(table, 2)
        out["log2_witness_probes"] = round(witness, 2)
    return out


def swap_reliability_at_131(m: int, l: float, k: int = 2,
                            miss_slack: float = 1.5, fp_slack: float = 1.0):
    """What the measured collision law implies for the swap at `n = 131`.

    **Derived, not measured.**  `swap_localisation` in the companion script fixes
    the *form* of both failure rates -- miss `~ k m/|F|`, false positive
    `~ (m/|F| + lambda)^k` -- and bounds their constants over a sixteenfold range
    of `m/|F|` and of `lambda`; the largest ratios it saw were `0.87` and `0.85`,
    and the slacks here are those rounded up.  Nothing at 131 was run.
    """
    base = 2.0 ** l
    lam = 2.0 ** (lc(l, m) - 131)
    collide = m / base
    per_summand_miss = miss_slack * k * collide
    per_target_fp = base * fp_slack * (collide + lam) ** k
    fail = m * per_summand_miss + per_target_fp
    return {"m": m, "l": round(l, 2), "queries_per_candidate": k,
            "log2_lambda": round(math.log2(lam), 2),
            "log2_collision_scale": round(math.log2(collide), 2),
            "log2_failure_probability_per_target": round(math.log2(fail), 2),
            "log2_relations_lost": round(math.log2(fail) + l, 2),
            "log2_relations_needed": round(l, 2),
            "fraction_of_relations_lost": round(fail, 6),
            "negligible": fail < 0.01}


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
        for mode in ("target_agnostic", "target_restricted"):
            for frob in (False, True):
                floors.append(detector_floor(m, N131, mode, frob))
    gaps = []
    for m in range(2, 9):
        a = next(f for f in floors if f["m"] == m
                 and f["mode"] == "target_agnostic" and not f["frobenius_stable"])
        b = next(f for f in floors if f["m"] == m
                 and f["mode"] == "target_restricted" and not f["frobenius_stable"])
        gaps.append({"m": m,
                     "log2_agnostic_floor": a["log2_floor"],
                     "agnostic_witness_route": a["witness_route"],
                     "log2_restricted_floor": b["log2_floor"],
                     "log2_gap": round(b["log2_floor"] - a["log2_floor"], 2),
                     "agnostic_beats_rho": a["log2_floor"] < log2_rho,
                     "restricted_beats_rho": b["log2_floor"] < log2_rho,
                     "log2_agnostic_table_entries": a["log2_table_entries"],
                     "log2_restricted_table_entries": b["log2_table_entries"]})
    e3 = {
        "question": "how much of the difficulty is deciding, and how much is "
                    "producing the witness?",
        "answer": "for a detector that is an algorithm on the target -- which is "
                  "every candidate in this repository -- none of it.  The witness "
                  "falls out of k|F| further whole-base queries by swapping one "
                  "candidate summand for a class-matched base point, and the "
                  "group operations that build those queries are absorbed by the "
                  "linear algebra already being paid.  Deciding IS localising.",
        "floors": floors, "gaps": gaps,
        "min_agnostic_floor": min(g["log2_agnostic_floor"] for g in gaps),
        "min_restricted_floor": min(g["log2_restricted_floor"] for g in gaps),
        "agnostic_rho_crossings": [g["m"] for g in gaps if g["agnostic_beats_rho"]],
        "restricted_rho_crossings": [g["m"] for g in gaps
                                     if g["restricted_beats_rho"]],
        "swap_reliability": [
            swap_reliability_at_131(
                m, next(f["l_star"] for f in floors if f["m"] == m
                        and f["mode"] == "target_agnostic"
                        and not f["frobenius_stable"]))
            for m in range(2, 9)],
        "note": "an agnostic detector may still build the meet-in-the-middle "
                "table, so it takes whichever witness route is cheaper: the table "
                "at m = 2, 3 (where the two floors coincide) and the swap from "
                "m = 4 up.  The restricted column prices an oracle that answers "
                "only for a distinguished family of targets and so cannot be "
                "swapped -- no proposed detector has that shape, which is what "
                "makes the gap a measure of an artificial restriction rather "
                "than of the curve's protection.  Both columns are derived, not "
                "measured, and both are floors -- no algorithm attains them.",
        "falsifier": "a proposed detector whose cost on R - P + Q exceeds its "
                     "cost on R by more than a constant, which would put the "
                     "restricted column back in play",
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
                   "research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md",
        "unit": "log2 group operations",
        "log2_rho_reference": round(log2_rho, 4),
        "E1_scale_model_ladder": e1,
        "E2_dimension_flatness": e2,
        "E3_detector_target_agnosticism": e3,
        "E4_large_primes": e4,
        "E5_yield_distribution": e5,
        "E6_orbit_union_bases": e6,
    }
    OUT.write_text(json.dumps(report, indent=2) + "\n")

    print(f"rho reference 2^{log2_rho:.4f}")
    print(f"E1 ladder: n = {[c['n'] for c in rungs]}")
    print(f"   end-to-end {e1['end_to_end_rungs']}, composed {e1['composed_rungs']}, "
          f"faithful subgroup {e1['faithful_subgroup_rungs']}")
    print(f"E3 target-agnostic floor bottoms at 2^{e3['min_agnostic_floor']:.2f}, "
          f"target-restricted at 2^{e3['min_restricted_floor']:.2f}, "
          f"rho at 2^{log2_rho:.2f}")
    for g in gaps:
        print(f"   m={g['m']}  agnostic 2^{g['log2_agnostic_floor']:6.2f} "
              f"({g['agnostic_witness_route']:18s})  "
              f"restricted 2^{g['log2_restricted_floor']:6.2f}  "
              f"gap 2^{g['log2_gap']:.2f}   beats rho: agnostic "
              f"{g['agnostic_beats_rho']}, restricted {g['restricted_beats_rho']}")
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
