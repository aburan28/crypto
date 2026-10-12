#!/usr/bin/env python3
"""Complete three-point relation sweeps over structured supports on ECC2K-130.

Each experiment is a *complete* search over its base: every unordered pair
with repetition, the third abscissa solved rather than scanned.  So a zero
result means the base has no homogeneous three-point relation at all, not
that a budget ran out.  When a budget does bind, the run is recorded as
incomplete and says so.

    python3 run_three_point.py [--budget-seconds N] [--out DIR] [--only NAME]

Writes one JSON per experiment plus a summary, under `results/`.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import random
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import relations as R                                  # noqa: E402
from fastfield import FastGF2m                         # noqa: E402

# The rho reference this thread is measured against: the <-1> x <pi> speed-up
# on the order-r subgroup of ECC2K-130.  Stated before measuring, per AGENTS.md.
LOG2_RHO_ITERATIONS = 60.9


def expected_relation_count(base_size: int, r: int) -> float:
    """Heuristic count of unordered triples from the base summing to O.

    There are ~B^3/6 unordered triples and each lands on the identity with
    probability ~1/r if the base behaved like a random set, so the expected
    yield is B^3/(6r).  At B ~ 4e3 and r ~ 2^129 this is ~2^-91: the honest
    statement is that a base of this size is nowhere near large enough for a
    three-point relation to be expected, and the sweeps below confirm the
    prediction rather than discover it.
    """
    if base_size < 3:
        return 0.0
    return base_size ** 3 / (6.0 * r)


def run_experiment(name, description, E, base_x, r, *, coeffs=None,
                   budget_seconds=None, verify_cap=20000, rng=None, log=print):
    started = time.time()
    base_x = sorted(set(base_x))
    n = len(base_x)
    total_pairs = n * (n + 1) // 2
    log(f"[{name}] base {n} abscissae, {total_pairs:,} pairs")

    state = {"last": 0.0, "stopped": False}

    def progress(seen, hits):
        now = time.time()
        if now - state["last"] >= 30:
            state["last"] = now
            frac = seen / total_pairs if total_pairs else 1.0
            el = now - started
            eta = el / frac - el if frac > 0 else float("nan")
            log(f"[{name}]   {seen:,}/{total_pairs:,} ({frac:6.2%}) "
                f"hits={hits} elapsed={el:6.0f}s eta={eta:6.0f}s")

    found, seen_pairs, complete = R.three_point_relations(
        E, base_x, allow_repeats=True, chunk=20000, progress=progress,
        budget_seconds=budget_seconds)
    elapsed = time.time() - started

    # Realising a triple costs about thirty field inversions, so a support
    # that yields millions of relations cannot have every one re-derived on
    # the curve inside any sane budget -- the small-coefficient constructed
    # base yields 2,044,422. When that happens a uniform random sample is
    # verified instead and the run records that it sampled, rather than
    # reporting a number that implies the whole set was checked.
    sampled = len(found) > verify_cap
    picker = rng or random.Random(0xC0FFEE)
    subset = picker.sample(found, verify_cap) if sampled else list(found)

    verified = []
    for t in subset:
        pts = R.realise_triple(E, *t)
        if pts is None:
            continue
        assert E.sum_points(list(pts)) is None, t
        verified.append({
            "abscissae": [hex(x) for x in t],
            "points": [[hex(p[0]), hex(p[1])] for p in pts],
        })

    out = {
        "name": name,
        "description": description,
        "base_size": n,
        "pairs_enumerated": seen_pairs,
        "pairs_total": total_pairs,
        "complete": complete,
        "elapsed_seconds": round(elapsed, 3),
        "pair_rate_per_second": round(seen_pairs / elapsed, 1) if elapsed else None,
        "relations_found": len(found),
        "verification_sampled": sampled,
        "relations_submitted_for_verification": len(subset),
        "relations_verified_on_curve": len(verified),
        "relations": verified[:64],
        "expected_by_yield_law": expected_relation_count(n, r),
        "log2_expected_by_yield_law": (
            math.log2(expected_relation_count(n, r))
            if expected_relation_count(n, r) > 0 else None),
    }

    if coeffs is not None:
        rows = R.relation_vectors(E, subset, coeffs, r)
        rank = R.useful_rank(rows, r)
        rank["rows_are_a_sample_of"] = len(found)
        out["rank_accounting"] = rank
        out["rows"] = [[str(a), str(b)] for a, b in rows[:64]]

    log(f"[{name}] done: {len(found)} relations, {elapsed:.0f}s, "
        f"complete={complete}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget-seconds", type=float, default=7200.0,
                    help="per-experiment wall-clock budget (default 7200)")
    ap.add_argument("--out", default=str(HERE / "results"))
    ap.add_argument("--only", default=None)
    ap.add_argument("--constructed-size", type=int, default=4000)
    ap.add_argument("--verify-cap", type=int, default=20000,
                    help="max relations re-derived on the curve per experiment")
    args = ap.parse_args()

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    F, E, P, Q, r = R.challenge_curve()
    assert E.on_curve(P) and E.on_curve(Q)
    assert E.mul(P, r) is None and E.mul(Q, r) is None

    rng = random.Random(20260917)
    experiments = []

    # 1. the complete weight-two support
    w2 = R.base_from_abscissae(E, R.weight_two_abscissae(F))
    experiments.append(dict(
        name="weight_two_complete",
        description="every z^i + z^j (i<j) whose abscissa carries a point",
        base_x=w2, coeffs=None))

    # 2. weight at most two
    w2a = R.base_from_abscissae(E, R.weight_at_most_two_abscissae(F))
    experiments.append(dict(
        name="weight_at_most_two",
        description="weight 0, 1 and 2 abscissae carrying points",
        base_x=w2a, coeffs=None))

    # 3. the Frobenius-stable part of the weight-two set
    stable = [(1 << i) ^ (1 << j) for i in range(66) for j in range(i + 1, 66)]
    experiments.append(dict(
        name="weight_two_frobenius_stable",
        description="the i<j<=65 part, the largest sigma-stable weight-two set",
        base_x=R.base_from_abscissae(E, stable), coeffs=None))

    # 4, 5. matched random controls, same size as the weight-two base
    for idx in (1, 2):
        experiments.append(dict(
            name=f"random_matched_{idx}",
            description=f"uniform random abscissae, size matched to weight_two_complete",
            base_x=R.random_matched_base(E, len(w2), rng), coeffs=None))

    # 6. a base constructed from the public P and Q, random coefficients
    cb, cc = R.constructed_base(E, P, Q, r, args.constructed_size, rng)
    experiments.append(dict(
        name="constructed_random_coeffs",
        description="points [u]P + [v]Q with uniform random u, v",
        base_x=cb, coeffs=cc))

    # 7. the same, with small coefficients -- the construction-equation case
    cb2, cc2 = R.constructed_base(E, P, Q, r, args.constructed_size, rng,
                                  small=64)
    experiments.append(dict(
        name="constructed_small_coeffs",
        description="points [u]P + [v]Q with small u, v; carries free relations",
        base_x=cb2, coeffs=cc2))

    summary = {
        "instance": "ECC2K-130",
        "curve": "y^2 + xy = x^3 + 1 over F_2^131",
        "subgroup_order_r": str(r),
        "log2_r": round(math.log2(r), 4),
        "rho_reference_log2_iterations": LOG2_RHO_ITERATIONS,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "budget_seconds_per_experiment": args.budget_seconds,
        "experiments": [],
    }

    for spec in experiments:
        if args.only and spec["name"] != args.only:
            continue
        res = run_experiment(spec["name"], spec["description"], E,
                             spec["base_x"], r, coeffs=spec["coeffs"],
                             budget_seconds=args.budget_seconds,
                             verify_cap=args.verify_cap)
        (outdir / f"{spec['name']}.json").write_text(json.dumps(res, indent=2))
        summary["experiments"].append({
            k: res[k] for k in
            ("name", "base_size", "pairs_total", "complete", "elapsed_seconds",
             "relations_found", "verification_sampled",
             "relations_submitted_for_verification",
             "relations_verified_on_curve", "log2_expected_by_yield_law")
        } | ({"rank_accounting": res["rank_accounting"]}
             if "rank_accounting" in res else {}))
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))

    print(json.dumps(summary["experiments"], indent=2))


if __name__ == "__main__":
    main()
