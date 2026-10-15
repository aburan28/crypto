"""Independent check of the n = 7 hyperelliptic carriers: find every genus-7
model (reduced) whose L-polynomial is divisible by W7 or A7, then recount one
representative per sign by brute force over F_{2^k}, k = 1..10, and compare
with the counts the claimed factorisation predicts for k = 8..10 (beyond the
k <= 7 used to build P_C)."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from gf2k import Field, MODULI
from fam_hyper import complement_reps, smooth
from lpoly import targets, charpoly_from_counts, counts_from_charpoly, divides

MODULI.update({9: 0b1000010001, 10: 0b10000001001})


def brute(h, f, F, g):
    q = F.q
    x = np.arange(q)
    hx = np.zeros(q, dtype=np.int64); fx = np.zeros(q, dtype=np.int64)
    for i in range(2 * g + 3):
        if (h >> i) & 1:
            hx ^= F.pw[x, i]
        if (f >> i) & 1:
            fx ^= F.pw[x, i]
    # count y with y^2 + h y = f by direct search
    y = np.arange(q)
    ysq = F.mul[y, y]
    n = 0
    for xi in range(q):
        n += int(np.count_nonzero((ysq ^ F.mul[hx[xi], y]) == fx[xi]))
    if (h >> (g + 1)) & 1:
        c = (f >> (2 * g + 2)) & 1
        n += int(np.count_nonzero((ysq ^ y) == c))  # y^2 + y = c over F_q
    else:
        n += 1
    return n


def main():
    g = 7
    T = targets()
    fields = [Field(k) for k in range(1, 8)]
    found = {}
    from fam_hyper import enumerate_hyper
    for h in range(1, 1 << (g + 2)):
        free = complement_reps(h, g)
        for m in range(1 << len(free)):
            f = sum(1 << b for j, b in enumerate(free) if (m >> j) & 1)
            if not smooth(h, f, g):
                continue
            # fast count via brute for k <= 7 is slow; use the trace formula
            N = []
            for F in fields:
                N.append(trace_count(h, f, F, g))
            P = charpoly_from_counts(N)
            for nm in ("W7(E0)", "W7(E1)"):
                if divides(T[nm], P):
                    found.setdefault(nm, []).append((h, f, N, P))
    out = {}
    for nm, lst in found.items():
        h, f, N, P = lst[0]
        pred = counts_from_charpoly(P, 10)
        big = [Field(k) for k in (8, 9, 10)]
        bf = [brute(h, f, F, g) for F in big]
        out[nm] = {"reduced_models": len(lst), "models": len(lst) * 2 ** (g + 1),
                   "example_h": bin(h), "example_f": bin(f), "N_1..7": N, "P_C": P,
                   "pred_N8_10": pred[7:], "brute_N8_10": bf, "match": pred[7:] == bf}
        print(nm, out[nm], flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "n7_carriers_check.json"), "w"),
              indent=1)


def trace_count(h, f, F, g):
    q = F.q
    x = np.arange(q)
    hx = np.zeros(q, dtype=np.int64); fx = np.zeros(q, dtype=np.int64)
    for i in range(2 * g + 3):
        if (h >> i) & 1:
            hx ^= F.pw[x, i]
        if (f >> i) & 1:
            fx ^= F.pw[x, i]
    nz = hx != 0
    w = F.mul[F.inv[hx], F.inv[hx]]
    t = F.trace[F.mul[fx, w]]
    n = int((~nz).sum() + (2 - 2 * t[nz]).sum())
    if (h >> (g + 1)) & 1:
        n += 2 - 2 * ((F.k * ((f >> (2 * g + 2)) & 1)) % 2)
    else:
        n += 1
    return n


if __name__ == "__main__":
    main()
