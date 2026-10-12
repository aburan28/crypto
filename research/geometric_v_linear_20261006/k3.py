"""Three-summand linear oracle in symmetric variables, at n = 131.

On y^2 + xy = x^3 + 1, S4 = G^2 with (sigma = elementary symmetric of x1..x4)
    G = D^2 + D sqrt(sigma4) + sigma1 sigma4 + sigma3,   D = sigma1 + sigma3.
With x4 = r = x(R) and e1, e2, e3 the elementary symmetric functions of X1, X2, X3:
    sigma1 = e1 + r, sigma3 = e3 + r e2, sigma4 = r e3, and
    G = e1^2 + r^2 + e3^2 + r^2 e2^2 + sqrt(r) (Y1 + r sqrt(e3) + Y3 + r Y2)
        + r Y4 + r^2 e3 + e3 + r e2
with Y1 = e1 sqrt(e3), Y2 = e2 sqrt(e3), Y3 = e3 sqrt(e3), Y4 = e1 e3.
Over geometric V = rho^2 <h^0, h^2, ..., h^(2l-2)> (theta = rho^2, g = h^2) every
block lies in a known geometric span, giving
    l + (2l-1) + (3l-2) + (5l-4) + (7l-6) + (9l-8) + (4l-3) = 31 l - 24
linear unknowns.  One trial is one n x (31l - 24) solve.
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "linearization_reach_20260930"))
import lr  # noqa: E402

# x^131 + x^8 + x^3 + x^2 + 1; lr.Field verifies irreducibility (x^(2^131) = x, 131 prime)
lr.IRR.setdefault(131, (1 << 131) | (1 << 8) | (1 << 3) | (1 << 2) | 1)


def pw(F, a, e):
    r = 1
    for _ in range(e):
        r = F.mul(r, a)
    return r


def sqrt(F, a):
    for _ in range(F.n - 1):
        a = F.sq(a)
    return a


def blocks(F, l, rho, h):
    """Spanning sets (field elements) for e1, e2, e3, Y1..Y4."""
    def span(rpow, exps):
        c = pw(F, rho, rpow)
        return [F.mul(c, pw(F, h, s)) for s in exps]
    return dict(
        e1=span(2, range(0, 2 * l - 1, 2)),
        e2=span(4, range(0, 4 * l - 3, 2)),
        e3=span(6, range(0, 6 * l - 5, 2)),
        Y1=span(5, range(0, 5 * l - 4)),
        Y2=span(7, range(0, 7 * l - 6)),
        Y3=span(9, range(0, 9 * l - 8)),
        Y4=span(8, range(0, 8 * l - 7, 2)),
    )


def columns(F, B, r):
    sr, r2 = sqrt(F, r), F.sq(r)
    cols, tags = [], []
    for b in B["e1"]:
        cols.append(F.sq(b)); tags.append("e1")
    for b in B["e2"]:
        cols.append(F.mul(r2, F.sq(b)) ^ F.mul(r, b)); tags.append("e2")
    for b in B["e3"]:
        cols.append(F.sq(b) ^ F.mul(r2, b) ^ b ^ F.mul(F.mul(sr, r), sqrt(F, b))); tags.append("e3")
    for b in B["Y1"]:
        cols.append(F.mul(sr, b)); tags.append("Y1")
    for b in B["Y2"]:
        cols.append(F.mul(F.mul(sr, r), b)); tags.append("Y2")
    for b in B["Y3"]:
        cols.append(F.mul(sr, b)); tags.append("Y3")
    for b in B["Y4"]:
        cols.append(F.mul(r, b)); tags.append("Y4")
    return cols, tags, r2          # system: sum x_k cols_k = r^2


def trial(F, l, rng):
    n = F.n
    while True:
        rho, h = rng.getrandbits(n) or 1, rng.getrandbits(n) or 2
        theta, g = F.sq(rho), F.sq(h)
        V = [F.mul(theta, pw(F, g, i)) for i in range(l)]
        X = []
        for _ in range(3):
            x = 0
            for v in V:
                if rng.getrandbits(1):
                    x ^= v
            X.append(x)
        P = [lr.lift(F, x) for x in X]
        if None in P or len(set(X)) < 3:
            continue
        R = lr.add(F, lr.add(F, P[0], P[1]), P[2])
        if R is None:
            continue
        break
    r = R[0]
    B = blocks(F, l, rho, h)
    cols, tags, rhs = columns(F, B, r)
    N = len(cols)
    rows = []
    for bit in range(n):
        row = 0
        for k, cv in enumerate(cols):
            if (cv >> bit) & 1:
                row |= 1 << k
        if (rhs >> bit) & 1:
            row |= 1 << N
        rows.append(row)
    rk, cons, piv = lr.eliminate(rows, N)
    e1 = X[0] ^ X[1] ^ X[2]
    e2 = F.mul(X[0], X[1]) ^ F.mul(X[0], X[2]) ^ F.mul(X[1], X[2])
    e3 = F.mul(F.mul(X[0], X[1]), X[2])
    # planted value must satisfy G = 0 (sanity on the closed form)
    sr = sqrt(F, r)
    s3 = sqrt(F, e3)
    G = (F.sq(e1) ^ F.sq(r) ^ F.sq(e3) ^ F.mul(F.sq(r), F.sq(e2))
         ^ F.mul(sr, F.mul(e1, s3) ^ F.mul(r, s3) ^ F.mul(e3, s3) ^ F.mul(r, F.mul(e2, s3)))
         ^ F.mul(r, F.mul(e1, e3)) ^ F.mul(F.sq(r), e3) ^ e3 ^ F.mul(r, e2))
    recovered = None
    if cons and rk == N:
        sol = [0] * N
        for c in sorted(piv, reverse=True):
            rr, v = piv[c], (piv[c] >> N) & 1
            for c2 in range(c + 1, N):
                if (rr >> c2) & 1:
                    v ^= sol[c2]
            sol[c] = v
        vals = {}
        for k, tag in enumerate(tags):
            if tag in ("e1", "e2", "e3") and sol[k]:
                idx = sum(1 for kk in range(k) if tags[kk] == tag)
                vals[tag] = vals.get(tag, 0) ^ B[tag][idx]
        recovered = (vals.get("e1", 0), vals.get("e2", 0), vals.get("e3", 0)) == (e1, e2, e3)
    return dict(N=N, rank=rk, consistent=cons, G_planted_zero=(G == 0), recovered=recovered)


def run(n, ls, trials, seed):
    F, rng = lr.Field(n), random.Random(seed)
    out = []
    for l in ls:
        rs = [trial(F, l, rng) for _ in range(trials)]
        row = dict(n=n, l=l, unknowns=31 * l - 24, trials=trials,
                   rank_hist={}, consistent=sum(r["consistent"] for r in rs),
                   G_planted_zero=sum(r["G_planted_zero"] for r in rs),
                   full_rank=sum(r["rank"] == r["N"] for r in rs),
                   recovered=sum(bool(r["recovered"]) for r in rs))
        for r in rs:
            row["rank_hist"][str(r["rank"])] = row["rank_hist"].get(str(r["rank"]), 0) + 1
        out.append(row)
        print(json.dumps(row), flush=True)
    return out


if __name__ == "__main__":
    n, trials, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    ls = [int(x) for x in sys.argv[4].split(",")]
    res = run(n, ls, trials, seed)
    if len(sys.argv) > 5:
        json.dump(res, open(sys.argv[5], "w"), indent=1)
