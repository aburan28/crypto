"""Apply the tight count to the trilinear chain at n = 131 and charge each trial.

A trial guesses nothing about t and resolves a1 bits of X1, all l bits of X2 and
b bits of X3 at once, so a relation costs 2^(n - a1 - l - b) trials.  Each trial is
one Macaulay solve of M columns, charged M^omega.  With Pr[decompose] ~ 1 the
per-attempt budget under rho is 2^(RHO - l), so
    gap_bits    = n - RHO - a1 - b
    gap_charged = n - RHO - a1 - b + omega * log2(M)
compared against the bilinear core's best charged gap (71.82 at omega = 2,
78.15 at omega = 2.81; ../xl_bilinear_core_20261006).
`slack` grants the attacker extra free bits of a1 beyond the count (still capped at l:
X1, like X3, is a factor-base point).
"""
import json
import math
import sys

from count import counts

N, RHO = 131, 60.8090
BILINEAR_BEST = {2.0: 71.82, 2.81: 78.15}


def scan(omega, slack=0):
    best = None
    for l in range(4, 31):
        for dA in range(1, 5):
            for dB in (1, 2):
                for dC in (1, 2):
                    for dT in range(1, 4):
                        for b in range(0, l + 1):
                            ok = [a1 for a1 in range(0, l + 1)  # X1 is a factor-base point: a1 <= l
                                  if (lambda r, c: r >= c - 1)(*counts(N, a1, l, b, (min(dA, a1) if a1 else 0, dB, dC, dT)))]
                            if not ok:
                                continue
                            a1 = max(ok)
                            caps = (min(dA, a1) if a1 else 0, dB, dC, dT)
                            _, M = counts(N, a1, l, b, caps)
                            a1g = min(l, a1 + slack)
                            gap = N - RHO - a1g - b + omega * math.log2(M)
                            row = dict(l=l, b=b, a1=a1, a1_with_slack=a1g, caps=caps,
                                       log2M=round(math.log2(M), 2),
                                       gap_bits=round(N - RHO - a1g - b, 2),
                                       gap_charged=round(gap, 2))
                            if best is None or row["gap_charged"] < best["gap_charged"]:
                                best = row
    return best


if __name__ == "__main__":
    res = {}
    for omega, slack in ((2.0, 0), (2.81, 0), (2.0, 4), (2.0, 8)):
        b = scan(omega, slack)
        verdict = "beats bilinear" if b["gap_charged"] < BILINEAR_BEST[omega] else "does not beat bilinear"
        res[f"omega={omega},slack={slack}"] = dict(best=b, bilinear_best=BILINEAR_BEST[omega], verdict=verdict)
        print(f"omega={omega} slack={slack}: {b}  vs bilinear {BILINEAR_BEST[omega]} -> {verdict}")
    if len(sys.argv) > 1:
        json.dump(res, open(sys.argv[1], "w"), indent=1)
