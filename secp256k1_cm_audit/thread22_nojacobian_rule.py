#!/usr/bin/env python3
"""
Thread 22, step 3: extract the actual obstruction rule, and apply it to
secp256k1.

The control run (thread22_realisability_baserate.py, full deg-5 + deg-6
coverage at p = 7, 13, 19) produced, for each p, the exact set of split
isogeny classes E_a x E_b that contain NO genus-2 Jacobian.  Those sets are
small and highly structured.  This script

  (1) fits candidate rules to the observed sets and reports exact agreement,
  (2) applies the winning rule to the six sextic twists of secp256k1,

so that Thread 18's "5/15 Howe-glueable" can be compared against "how many of
the 15 product classes actually contain a genus-2 Jacobian".

The rule that fits is the specialisation to E_1 x E_2 of the classification in

    E. Howe, E. Nart, C. Ritzenthaler, "Jacobians in isogeny classes of
    abelian surfaces over finite fields", Ann. Inst. Fourier 59 (2009)
    2617-2634, eprint arXiv:math/0607515.

For P_i(T) = T^2 - t_i T + p the reduced resultant of P_1, P_2 is governed by
|t_1 - t_2|: the product class fails to contain a Jacobian exactly when the
two elliptic factors are "too close" (|t_1 - t_2| = 1), plus one exceptional
diagonal case at the extreme Hasse trace with CM by Q(sqrt(-3)).

IMPORTANT SCOPE NOTE
--------------------
"Contains a Jacobian" is an EXISTENCE statement about the isogeny class.  It
is strictly weaker than "there is an efficiently computable (2,2)-cover
C -> E_i defined over F_p".  Howe's (H1)+(H2)+(H3) is a sufficient condition
for an EXPLICIT gluing construction; this script measures the existence
question only.  See the log entry for why the distinction matters for the
cover-attack surface in PAPER_STRUCTURAL_COMPLETENESS.md.
"""

import sys
from math import isqrt

sys.path.insert(0, "secp256k1_cm_audit")

# Observed non-realised split classes, as (t_a, t_b) with t_a <= t_b.
# Transcribed from thread22_baserate_output.txt (COMPLETE coverage rows only).
OBSERVED = {
    7: [(-5, -5), (-5, -4), (-4, -3), (-3, -2), (-2, -1), (-1, 0),
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 5)],
    13: [(-7, -7), (-7, -6), (-6, -5), (-5, -4), (-4, -3), (-3, -2),
         (-2, -1), (-1, 0), (0, 1), (1, 2), (2, 3), (3, 4), (4, 5),
         (5, 6), (6, 7), (7, 7)],
    19: [(-8, -7), (-7, -6), (-6, -5), (-5, -4), (-4, -3), (-3, -2),
         (-2, -1), (-1, 0), (0, 1), (1, 2), (2, 3), (3, 4), (4, 5),
         (5, 6), (6, 7), (7, 8)],
}


def traces(p):
    m = isqrt(4 * p)
    return [t for t in range(-m, m + 1) if t * t <= 4 * p]


def all_pairs(p):
    ts = traces(p)
    return [(ts[a], ts[b]) for a in range(len(ts)) for b in range(a, len(ts))]


# ------------------------------------------------------------ candidate rules

def rule_adjacent(p, ta, tb):
    return abs(ta - tb) == 1


def rule_adjacent_plus_extremal(p, ta, tb):
    if abs(ta - tb) == 1:
        return True
    return ta == tb and ta * ta == 4 * p - 3


def rule_adjacent_plus_any_diagonal(p, ta, tb):
    return abs(ta - tb) == 1 or ta == tb


RULES = [
    ("|t_a - t_b| = 1", rule_adjacent),
    ("|t_a - t_b| = 1  OR  (t_a = t_b and t_a^2 = 4p-3)",
     rule_adjacent_plus_extremal),
    ("|t_a - t_b| = 1  OR  t_a = t_b", rule_adjacent_plus_any_diagonal),
]


def fit(out):
    w = out.write
    w("Rule fitting against exhaustively verified non-Jacobian classes\n")
    w("(complete deg-5 + deg-6 enumeration; PARI-cross-checked engine)\n\n")
    best = None
    for name, fn in RULES:
        tot_fp = tot_fn = 0
        per = []
        for p, obs in OBSERVED.items():
            obs = set(obs)
            pred = {(ta, tb) for (ta, tb) in all_pairs(p) if fn(p, ta, tb)}
            fp = pred - obs          # predicted no-Jac but a Jacobian exists
            fn_ = obs - pred         # missed: no Jacobian but rule says fine
            per.append((p, len(obs), len(pred), len(fp), len(fn_)))
            tot_fp += len(fp)
            tot_fn += len(fn_)
        w(f"rule: {name}\n")
        w(f"  {'p':>4} {'observed':>9} {'predicted':>10} "
          f"{'over':>5} {'under':>6}\n")
        for p, no, npd, nfp, nfn in per:
            w(f"  {p:>4} {no:>9} {npd:>10} {nfp:>5} {nfn:>6}\n")
        w(f"  TOTAL mismatches: {tot_fp + tot_fn} "
          f"(over {tot_fp}, under {tot_fn})\n\n")
        if tot_fp + tot_fn == 0 and best is None:
            best = (name, fn)
    if best:
        w(f"EXACT FIT: {best[0]}\n\n")
    else:
        w("no candidate rule fits exactly\n\n")
    return best


# --------------------------------------------------------------- secp256k1

P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141


def secp256k1_traces():
    """The six sextic-twist traces, from the CM decomposition 4p = t0^2 + 3b^2.

    Reproduces the values in the 2026-07-21 log entry (commit a287abc).
    """
    t0 = P + 1 - N
    bsq = (4 * P - t0 * t0) // 3
    bv = isqrt(bsq)
    assert bv * bv == bsq, "4p - t0^2 not 3 * square"
    av = (t0 + bv) // 2
    assert 2 * av - bv == t0 and av * av - av * bv + bv * bv == P
    return [2 * av - bv, av - 2 * bv, -(av + bv),
            -(2 * av - bv), 2 * bv - av, av + bv], av, bv


def apply_to_secp256k1(out, rule):
    w = out.write
    ts, av, bv = secp256k1_traces()
    w("=" * 78 + "\nApplication to secp256k1\n" + "=" * 78 + "\n")
    w(f"CM decomposition: a = {av}\n                  b = {bv}\n")
    w("sextic-twist traces t_0..t_5:\n")
    for k, t in enumerate(ts):
        w(f"  t_{k} = {t}\n")
    w(f"\n4p - 3 is a perfect square? "
      f"{isqrt(4*P-3)**2 == 4*P-3}   (extremal exception inapplicable)\n\n")
    name, fn = rule
    w(f"applying rule: {name}\n\n")
    w(f"{'pair':>6} {'|t_i - t_j|':>78}  {'no Jac?':>8}\n")
    blocked = 0
    for i in range(6):
        for j in range(i + 1, 6):
            d = abs(ts[i] - ts[j])
            bad = fn(P, ts[i], ts[j])
            blocked += bad
            w(f"({i},{j})  {d:>78}  {'YES' if bad else 'no':>8}\n")
    w(f"\nproduct classes E_i x E_j with NO genus-2 Jacobian : {blocked}/15\n")
    w(f"product classes that DO contain a genus-2 Jacobian  : {15-blocked}/15\n")
    w("\ncompare: Thread 18 (commit a287abc) reported 5/15 pairs satisfying\n"
      "the explicit Howe gluing conditions (H1)+(H2)+(H3).\n")


def main():
    path = "secp256k1_cm_audit/thread22_nojacobian_rule_output.txt"
    with open(path, "w") as out:
        out.write("Thread 22: obstruction rule for genus-2 realisability "
                  "of split classes\n\n")
        best = fit(out)
        if best is None:
            out.write("aborting: no exact rule, refusing to extrapolate "
                      "to secp256k1\n")
        else:
            apply_to_secp256k1(out, best)
    print(open(path).read())


if __name__ == "__main__":
    main()
