"""Extrapolate bidegree XL to n = 131 with the count rule (see PROTOCOL.md).

Count rule (fitted on the exploratory sweep, frozen in PROTOCOL.md): XL at
(D1, D2) != (1, 1) resolves b bits of X3 when
    n * Mc(D1-1) * Md(D2-1) >= Mc(D1) * Md(D2) - 1,
optionally loosened by `slack` extra bits in the attacker's favour.  Plain
linearization (1, 1) uses the validated Amendment 1 reach
    l(b+1) + b - b(b+1)/2 <= n
with the polynomial structured resolver, charged at M = l*b + l + b + 1.
Mc(D) = sum_{i<=D} C(l, i), Md(D) = sum_{j<=D} C(b, j), b <= l.
Per-relation cost of the oracle: 2^(n - l - b) trials, each one Macaulay solve of
M = Mc(D1) * Md(D2) columns, charged at M^omega (omega = 2: sparse, optimistic
for the attacker; 2.81: dense Strassen).  When Pr[decompose] ~ 1 the per-attempt
budget under rho is 2^(RHO - l), so
    gap_bits    = n - RHO - b                      (trial cost ignored)
    gap_charged = n - RHO - b + omega * log2(M)    (trial cost charged)
"""
import json
import math
import sys
from math import comb

N, RHO = 131, 60.8090


def Ms(k, D):
    return sum(comb(k, i) for i in range(0, max(D, -1) + 1)) if D >= 0 else 0


def bmax(l, D1, D2, n=N, slack=0):
    best = None
    for b in range(0, l + 1):
        if (D1, D2) == (1, 1):
            ok = l * (b + 1) + b - b * (b + 1) // 2 <= n
        else:
            ok = n * Ms(l, D1 - 1) * Ms(b, D2 - 1) >= Ms(l, D1) * Ms(b, D2) - 1
        if ok:
            best = b
    if best is None:
        return None
    return best if (D1, D2) == (1, 1) else min(l, best + slack)


def rows(omega, slack=0):
    out = []
    for l in range(6, 31):
        for D1 in range(1, 7):
            for D2 in range(1, 4):
                b = bmax(l, D1, D2, slack=slack)
                if b is None:
                    continue
                M = Ms(l, D1) * Ms(b, D2)
                out.append(dict(l=l, D1=D1, D2=D2, b=b, log2M=round(math.log2(M), 2),
                                gap_bits=round(N - RHO - b, 2),
                                gap_charged=round(N - RHO - b + omega * math.log2(M), 2)))
    return out


if __name__ == "__main__":
    res = {}
    for omega, slack in ((2.0, 0), (2.81, 0), (2.0, 2), (2.0, 4)):
        rs = rows(omega, slack)
        lin = [r for r in rs if (r["D1"], r["D2"]) == (1, 1)]
        best_lin = min(lin, key=lambda r: r["gap_charged"])
        best_any = min(rs, key=lambda r: r["gap_charged"])
        best_bits = min(rs, key=lambda r: (r["gap_bits"], r["gap_charged"]))
        res[f"omega={omega},slack={slack}"] = dict(best_linearization=best_lin, best_xl_charged=best_any,
                               best_xl_bits_only=best_bits)
        print(f"omega={omega} slack={slack}")
        print("  best linearization (charged):", best_lin)
        print("  best XL, charged            :", best_any)
        print("  best XL, bits only          :", best_bits)
        # per-l best, to see where XL beats linearization
        for l in (11, 14, 16, 20, 24, 30):
            cand = [r for r in rs if r["l"] == l]
            bl = [r for r in cand if (r["D1"], r["D2"]) == (1, 1)][0]
            bx = min(cand, key=lambda r: r["gap_charged"])
            print(f"   l={l}: lin b={bl['b']} charged gap {bl['gap_charged']}; "
                  f"best XL D=({bx['D1']},{bx['D2']}) b={bx['b']} log2M={bx['log2M']} charged gap {bx['gap_charged']}")
    if len(sys.argv) > 1:
        json.dump(res, open(sys.argv[1], "w"), indent=1)
