"""Tight count for multihomogeneous XL on the trilinear chain.

Rows = n * |mult(E1)| + n * |mult(E2)|, columns = |monomials within caps|.
Predicted resolvable when rows >= columns - 1 (the rule that fitted the
bilinear core within one bit, ../xl_bilinear_core_20261006).
E1 = S3(X1, t, xR) has block degree (1, 0, 0, 1) (A-degree 0 when a1 = 0);
E2 = S3(X2, X3, t) has block degree (0, 1, 1, 1).
"""
from math import comb


def Ms(k, D):
    return sum(comb(k, i) for i in range(min(D, k) + 1)) if D >= 0 else 0


def counts(n, a1, l, b, caps):
    dA, dB, dC, dT = caps
    cols = Ms(a1, dA) * Ms(l, dB) * Ms(b, dC) * Ms(n, dT)
    e1A = 1 if a1 > 0 else 0
    r1 = Ms(a1, dA - e1A) * Ms(l, dB) * Ms(b, dC) * Ms(n, dT - 1)
    r2 = Ms(a1, dA) * Ms(l, dB - 1) * Ms(b, dC - 1) * Ms(n, dT - 1)
    return n * (r1 + r2), cols


def predicted(n, a1, l, b, caps):
    rows, cols = counts(n, a1, l, b, caps)
    return rows >= cols - 1


if __name__ == "__main__":
    for n in (11, 13):
        for l, b in ((3, 1), (4, 2)):
            for dA in (1, 2):
                row = []
                for a1 in range(0, 5):
                    m = next((dT for dT in range(1, 8) if predicted(n, a1, l, b, (min(dA, a1), 1, 1, dT))), None)
                    row.append(m)
                print(f"n={n} l={l} b={b} dA={dA}: predicted min dT for a1=0..4 -> {row}")
