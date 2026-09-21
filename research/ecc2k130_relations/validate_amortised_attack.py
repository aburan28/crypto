#!/usr/bin/env python3
"""Run the amortised attack for real, on curves small enough to finish.

`target_boundary.py` prices a logarithm as one meet-in-the-middle table built
**once** plus a stream per target attempt.  That structure is the whole of
the 2^7.3 correction over charging the table per attempt, and its optimum at
ECC2K-130 is lopsided in a way worth exercising rather than trusting:

  * the support is a **single Frobenius orbit** -- 131 points carrying two
    unknowns, its own orbit logarithm and `d`;
  * the table holds canonical `sigma`-classes of signed `(n-1)`-subset sums,
    built once;
  * each target streams the `2B` signed single points and looks up
    `T - R` in that table.

Every probe is counted **and** its canonicalisation is counted, because a
quotiented table cannot be probed without one: `canonicalisations` is
reported alongside `stream_operations` rather than left out of the unit.
This script probes all four `E[4]` elements where ECC2K-130 need probe only
the two that lie in `H`; that is small-curve conservatism, not a saving.

So this script builds the table once, streams targets against it until it
has `T + 1` relations, solves, and checks `[d]P = Q`.  If the structure the
cost model prices did not actually produce usable relations, or if two
relations were not enough at `T = 1`, this would fail rather than merely
score badly.

No discrete logarithm enters the solve: the Frobenius scalar comes from
`sigma^2 + sigma + 2 = 0` with the root chosen by a point identity.
"""

from __future__ import annotations

import json
import math
import random
import sys
from itertools import combinations, product
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from normalbasis import NormalSupport, find_normal_elements   # noqa: E402
from planted import frobenius_scalar, solve_mod_p             # noqa: E402
from smallcurve import four_torsion, small_curve              # noqa: E402


def _canon(E, S, m):
    """Least point in the `sigma`-orbit of `S`, the table's key."""
    best, cur = S, S
    for _ in range(m - 1):
        cur = E.frobenius(cur)
        if cur is not None and (best is None or cur < best):
            best = cur
    return best


def attack(mdeg: int, n: int, want_T: int, seed: int = 11,
           max_targets: int = 400000):
    F, E, order = small_curve(mdeg)
    e4 = four_torsion(E, order)
    r = order
    while r % 2 == 0:
        r //= 2
    rng = random.Random(seed)

    sup = None
    for a in find_normal_elements(F, 80, rng):
        s = NormalSupport(F, E, a)
        if s.orbit_count == want_T:
            sup = s
            break
    if sup is None:
        return {"m": mdeg, "orbits": want_T, "skipped": "no support of that size"}
    reps, B, T, m = sup.orbit_points(E), sup.size, sup.orbit_count, sup.m

    P = None
    while P is None:
        pts = E.points_over(rng.getrandbits(mdeg) & F.mask)
        if pts:
            cand = E.mul(pts[0], order // r)
            if cand is not None:
                P = cand
    d = rng.randrange(2, r)
    Q = E.mul(P, d)
    sf = frobenius_scalar(E, P, r)

    # -- built once, reused for every target --------------------------------
    table = {}
    for comb in combinations(range(B), n - 1):
        for sg in product((1, -1), repeat=n - 1):
            S = E.sum_points([reps[i] if g > 0 else E.neg(reps[i])
                              for i, g in zip(comb, sg)])
            if S is not None:
                table.setdefault(_canon(E, S, m), (comb, sg, S))

    rows, rhs, attempts, stream_ops, canon_ops = [], [], 0, 0, 0
    while len(rows) < T + 1 and attempts < max_targets:
        attempts += 1
        a_, b_ = rng.randrange(r), rng.randrange(1, r)
        tgt = E.add(E.mul(P, a_), E.mul(Q, b_))
        if tgt is None:
            continue
        found = None
        for i in range(B):                       # -- streamed per target --
            for si in (1, -1):
                Ri = reps[i] if si > 0 else E.neg(reps[i])
                for Tt in e4:
                    stream_ops += 1
                    want = E.add(E.add(Tt, E.neg(tgt)), E.neg(Ri))
                    if want is None:
                        continue
                    canon_ops += 1          # each probe is canonicalised
                    ent = table.get(_canon(E, want, m))
                    if ent is None:
                        continue
                    comb, sg, S0 = ent
                    cur, j = S0, None
                    for jj in range(m):
                        if cur == want:
                            j = jj
                            break
                        cur = E.frobenius(cur)
                    if j is not None:
                        found = (i, si, comb, sg, j)
                        break
                if found:
                    break
            if found:
                break
        if not found:
            continue
        i, si, comb, sg, j = found
        pts = [(i, si)] + [(sup.index(divmod(c, m)[0], divmod(c, m)[1] + j), g)
                           for c, g in zip(comb, sg)]
        if len({p for p, _ in pts}) < n:
            continue                              # distinct abscissae only
        row = [0] * (T + 1)
        for idx, g in pts:
            t_, k_ = divmod(idx, m)
            row[t_] = (row[t_] + g * pow(sf, k_, r)) % r
        row[T] = 4 * b_ % r
        rows.append(row)
        rhs.append((-4 * a_) % r)

    if len(rows) < T + 1:
        return {"m": mdeg, "support_size": B, "orbits": T, "relation_length": n,
                "recovered": False, "reason": "ran out of targets",
                "attempts": attempts}
    sol = solve_mod_p(rows, rhs, r)
    ok = sol is not None and E.mul(P, sol[T]) == Q
    lc = sum(math.log2(B - i) for i in range(n)) - math.log2(math.factorial(n))
    pred = min(1.0, 2 ** (n + lc - math.log2(r)))
    return {
        "m": mdeg, "support_size": B, "orbits": T, "unknowns": T + 1,
        "relation_length": n, "table_entries_built_once": len(table),
        "relations": len(rows), "target_attempts": attempts,
        "stream_operations": stream_ops,
        "canonicalisations": canon_ops,
        "e4_translates_probed": len(e4),
        "predicted_decomposition_rate": round(pred, 6),
        "measured_decomposition_rate": round(len(rows) / attempts, 6),
        "planted_d": d, "recovered_d": (sol[T] if sol else None),
        "verified_by_point_identity": bool(ok),
    }


CELLS = [(13, 3, 1), (13, 3, 2), (13, 4, 1), (19, 3, 1), (19, 3, 2), (19, 4, 1)]


def aggregate(rows, expected):
    """Summarise cells so that a missing or failed one cannot be hidden.

    Filtering to cells that carry `recovered_d` and reporting "N of N" over
    the survivors is how a validation script reports success while a cell
    quietly did not run.  A skipped support and a cell that ran out of
    targets both return without `recovered_d`, so both would vanish from the
    denominator.  Every cell is accounted for here, and `all_recovered`
    requires the full set.
    """
    recovered = [r for r in rows if r.get("verified_by_point_identity")]
    skipped = [r for r in rows if "skipped" in r]
    failed = [r for r in rows
              if "skipped" not in r and not r.get("verified_by_point_identity")]
    return {
        "cells_expected": expected,
        "cells_run": len(rows),
        "cells_recovered": len(recovered),
        "cells_skipped": len(skipped),
        "cells_failed": len(failed),
        "all_recovered": (len(rows) == expected
                          and len(recovered) == expected
                          and not skipped and not failed),
        "recovered": recovered, "skipped": skipped, "failed": failed,
    }


def main():
    rows = [attack(*c) for c in CELLS]
    agg = aggregate(rows, len(CELLS))
    data = {
        "instance": "small analogues of ECC2K-130",
        "structure": ("table of canonical sigma-classes of signed (n-1)-subset "
                      "sums built ONCE; single signed points streamed per "
                      "target; reused across all targets"),
        "cells": rows,
        "cells_expected": agg["cells_expected"],
        "cells_recovered": agg["cells_recovered"],
        "cells_skipped": agg["cells_skipped"],
        "cells_failed": agg["cells_failed"],
        "all_recovered": agg["all_recovered"],
        "verdict": (
            f"The amortised structure recovers the planted logarithm in "
            f"{agg['cells_recovered']} of {agg['cells_expected']} cells, "
            f"verified by the point identity [d]P = Q"
            + (f" ({agg['cells_skipped']} skipped, {agg['cells_failed']} "
               f"failed)" if agg["cells_skipped"] or agg["cells_failed"] else "")
            + ". The single-orbit cells (T = 1) close a two-unknown system "
              "from two relations, which is the shape the ECC2K-130 optimum "
              "uses."),
    }
    (HERE / "results" / "amortised_attack.json").write_text(json.dumps(data, indent=2))
    hdr = (f"{'m':>3} {'B':>4} {'T':>2} {'n':>2} {'table':>7} {'rels':>5} "
           f"{'attempts':>9} {'pred p':>9} {'meas p':>9} {'[d]P==Q':>8}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        if "recovered_d" not in r:
            why = r.get("skipped") or r.get("reason") or "no result"
            print(f"{r.get('m','?'):>3} {'-':>4} {r.get('orbits','?'):>2} "
                  f"{r.get('relation_length','-'):>2}   {why}")
            continue
        print(f"{r['m']:>3} {r['support_size']:>4} {r['orbits']:>2} "
              f"{r['relation_length']:>2} {r['table_entries_built_once']:>7} "
              f"{r['relations']:>5} {r['target_attempts']:>9} "
              f"{r['predicted_decomposition_rate']:>9.4f} "
              f"{r['measured_decomposition_rate']:>9.4f} "
              f"{str(r['verified_by_point_identity']):>8}")
    print()
    print(data["verdict"])
    if not data["all_recovered"]:
        raise SystemExit("not every cell recovered; see the table above")


if __name__ == "__main__":
    main()
