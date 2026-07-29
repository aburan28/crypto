#!/usr/bin/env python3
"""
Thread 22, control experiment: what is the BASE RATE of realisability?

thread22_howe_realizability.py asks whether the codebase's (H1)+(H2)+(H3)
test predicts that a genus-2 curve C/F_p exists with

    charpoly(Frob | Jac C) = P_{E_i}(T) * P_{E_j}(T)

for the sextic twists E_0..E_5 of a j=0 curve.  That question is only
meaningful if NOT every split isogeny class is realised.  If essentially all
of them are, "pair (i,j) is realised" carries no information and neither
false negatives nor false positives mean anything.

This control enumerates genus-2 curves over F_p exactly as the main script
does, but takes as targets EVERY split class E_a x E_b built from every
achievable elliptic trace over F_p (|t| <= 2 sqrt p, all t occur for p prime),
and reports

    realised split classes / all split classes

i.e. the base rate.  The Howe-pair statistic is only interesting to the extent
that it deviates from this.

The non-realised classes are also printed: by Howe-Nart-Ritzenthaler
("Jacobians in isogeny classes of abelian surfaces over finite fields",
Ann. Inst. Fourier 59 (2009)) the isogeny classes of abelian surfaces
containing no Jacobian are explicitly classified, so this list is a check on
the enumeration as well as a measurement.

Output: secp256k1_cm_audit/thread22_baserate_output.txt
"""

import sys
from math import isqrt

import numpy as np

sys.path.insert(0, "secp256k1_cm_audit")
from thread22_howe_realizability import (          # noqa: E402
    chi2_table, find_squarefree, howe_table, legendre_table, non_residue,
    sextic_twist_data, sweep,
)


def elliptic_traces(p):
    """Traces of Frobenius achieved by some elliptic curve over F_p.

    For p prime every t with |t| <= 2 sqrt p occurs (Deuring/Waterhouse):
    ordinary for p nmid t, and t = 0 for the supersingular class.
    """
    m = isqrt(4 * p)
    return [t for t in range(-m, m + 1) if t * t <= 4 * p]


def run(p, b, out, full=True):
    chi = legendre_table(p)
    nu = non_residue(p, chi)
    chi2 = chi2_table(p, nu)
    g, bs, ts, ns, pats = sextic_twist_data(p, b, chi)
    rows = howe_table(p, ts, ns, pats)

    traces = elliptic_traces(p)
    classes = {}
    for a in range(len(traces)):
        for c in range(a, len(traces)):
            ta, tc = traces[a], traces[c]
            classes[(ta + tc, ta * tc + 2 * p)] = (ta, tc)
    target_s1 = sorted({k[0] for k in classes})

    families = [(5, None)] + ([(6, 1), (6, nu)] if full else [])
    results = []
    for deg, c6 in families:
        print(f"    baserate sweep p={p} deg={deg} c6={c6}", file=sys.stderr)
        results.append((deg, c6) + sweep(p, chi, chi2, nu, target_s1, deg, c6))

    realised = {}
    for key in classes:
        hit = None
        for deg, c6, coeffs, s1v, s2v in results:
            if coeffs is None or hit is not None:
                continue
            sel = np.nonzero((s1v == key[0]) & (s2v == key[1]))[0]
            if sel.shape[0]:
                cs, _ = find_squarefree(p, deg, c6, coeffs, sel)
                if cs is not None:
                    hit = (deg, c6, cs)
        realised[key] = hit

    nreal = sum(1 for k in classes if realised[k] is not None)
    cover = "COMPLETE (deg 5 + deg 6)" if full else "PARTIAL (deg 5 only)"

    w = out.write
    w("=" * 78 + "\n")
    w(f"p = {p}   b = {b}   coverage: {cover}\n")
    w("=" * 78 + "\n")
    w(f"elliptic traces available      : {len(traces)}  "
      f"({traces[0]} .. {traces[-1]})\n")
    w(f"distinct split classes E_a x E_b: {len(classes)}\n")
    w(f"realised by a genus-2 Jacobian : {nreal}\n")
    w(f"BASE RATE                      : {nreal/len(classes):.4f}\n\n")

    missing = [(k, v) for k, v in classes.items() if realised[k] is None]
    w(f"split classes with NO genus-2 Jacobian ({len(missing)}):\n")
    if not missing:
        w("  (none)\n")
    for (s1, s2), (ta, tc) in sorted(missing):
        w(f"  s1={s1:>4} s2={s2:>6}   from traces (t_a,t_b)=({ta:>3},{tc:>3})"
          f"   charpoly = T^4 - {s1}T^3 + {s2}T^2 - {p*s1}T + {p*p}\n")

    nq = sum(1 for r in rows if r["ALL"])
    hp = [r for r in rows if realised.get((r["s1"], r["s2"])) is not None]
    w(f"\nsextic-twist pairs: Howe-qualifying {nq}/15, "
      f"realised {len(hp)}/15\n")
    w("(compare with the base rate above before reading anything into this)\n\n")
    out.flush()
    return dict(p=p, ncls=len(classes), nreal=nreal, full=full,
                nmissing=len(missing), nq=nq, npair_real=len(hp))


def main():
    cases = [(7, 1, True), (13, 7, True), (19, 7, True),
             (31, 7, False), (37, 7, False), (43, 7, False)]
    path = "secp256k1_cm_audit/thread22_baserate_output.txt"
    if len(sys.argv) > 1:
        # e.g. "31:full" -> single out-of-sample prime, complete coverage
        cases = []
        for spec in sys.argv[1:]:
            pp, _, mode = spec.partition(":")
            cases.append((int(pp), 7 if int(pp) != 7 else 1, mode == "full"))
        path = ("secp256k1_cm_audit/thread22_baserate_output_"
                + "_".join(sys.argv[1:]).replace(":", "") + ".txt")
    summ = []
    with open(path, "w") as out:
        out.write("Thread 22 control: base rate of genus-2 realisability "
                  "among split isogeny classes\n\n")
        for p, b, full in cases:
            summ.append(run(p, b, out, full=full))
        out.write("=" * 78 + "\nSUMMARY\n" + "=" * 78 + "\n")
        out.write(f"{'p':>5} {'coverage':>9} {'classes':>8} {'realised':>9} "
                  f"{'baserate':>9} {'noJac':>6} {'Hqual':>6} {'pairReal':>9}\n")
        for s in summ:
            out.write(f"{s['p']:>5} {'full' if s['full'] else 'deg5':>9} "
                      f"{s['ncls']:>8} {s['nreal']:>9} "
                      f"{s['nreal']/s['ncls']:>9.4f} {s['nmissing']:>6} "
                      f"{s['nq']:>6} {s['npair_real']:>9}\n")
    print(open(path).read())


if __name__ == "__main__":
    main()
