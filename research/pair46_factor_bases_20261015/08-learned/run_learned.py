#!/usr/bin/env python3
"""Item 08 learned factor base: learned selectors vs greedy-optimised and
random bases on held-out m=2 decomposition targets, toy prime-order curves.
See PROTOCOL.md (frozen).  Usage: python3 run_learned.py [results.json]"""
import json, sys, time, math, random
import numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else "results.json"
import os
BITS = [14] if os.environ.get("SMOKE") else [14, 16, 18]
CURVES = 4
SEEDS = 3
MUS = [0.25, 1.0]
QS = [3, 5, 7, 11, 13]
T0 = time.time()


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13):
        if n % q == 0:
            return n == q
    return all(n % d for d in range(17, int(n ** 0.5) + 1, 2))


def find_curves(bits, count):
    rng = random.Random(bits)
    out = []
    while len(out) < count:
        p = rng.randrange(2 ** (bits - 1), 2 ** bits)
        if not is_prime(p):
            continue
        xs = np.arange(p, dtype=np.int64)
        sq = np.zeros(p, dtype=np.int8)
        sq[(xs * xs) % p] = 1
        sq[0] = 0
        for _ in range(40):
            a, b = rng.randrange(p), rng.randrange(1, p)
            if (4 * a ** 3 + 27 * b * b) % p == 0:
                continue
            f = (xs * xs % p * xs + a * xs + b) % p
            chi = np.where(f == 0, 0, np.where(sq[f] == 1, 1, -1))
            N = p + 1 + int(chi.sum())
            if is_prime(N) and N > p // 2:
                out.append((p, a, b, N, f, sq))
                break
    return out


def multiples(p, a, b, N, f, sq):
    """Table k -> (x,y) of kG, k=1..N-1. G = first point with smallest x."""
    x0 = int(np.nonzero((sq[f] == 1))[0][0])
    fx = int(f[x0])
    y0 = next(y for y in range(1, p) if y * y % p == fx)
    X = np.full(N, -1, dtype=np.int64)
    Y = np.full(N, -1, dtype=np.int64)
    x, y = x0, y0
    X[1], Y[1] = x, y
    # 2G
    lam = (3 * x * x + a) * pow(2 * y, -1, p) % p
    x2 = (lam * lam - 2 * x) % p
    y2 = (lam * (x - x2) - y) % p
    cx, cy = x2, y2
    for k in range(2, N):
        X[k], Y[k] = cx, cy
        if k == N - 1:
            break
        if cx == x0:  # cy == -y0 only at k = N-1
            raise RuntimeError("unexpected")
        lam = (cy - y0) * pow(cx - x0, -1, p) % p
        nx = (lam * lam - cx - x0) % p
        cy = (lam * (cx - nx) - cy) % p
        cx = nx
    assert X[N - 1] == x0 and (Y[N - 1] + y0) % p == 0
    return X, Y


def canon(k, N):
    k = k % N
    return np.minimum(k, N - k)


def counts_all(base, N):
    """Exact c(r) for all r: cyclic self-convolution of the indicator of S."""
    ind = np.zeros(N)
    ind[base % N] = 1
    ind[(-base) % N] = 1
    F = np.fft.rfft(ind)
    return np.rint(np.fft.irfft(F * F, n=N)).astype(np.int64)


def metrics(base, N, train_mask):
    c = counts_all(base, N)
    held = ~train_mask
    held[0] = False
    tr = train_mask.copy()
    return {
        "held_success": float((c[held] >= 1).mean()),
        "held_mean_count": float(c[held].mean()),
        "train_success": float((c[tr] >= 1).mean()),
        "sum_count_all_r": int(c.sum()),
        "S_size": int(2 * len(base)),
    }


def greedy(N, L, train, rng):
    H = (N - 1) // 2
    score = np.zeros(H + 1)
    score[0] = -1e18
    inS = np.zeros(N, dtype=bool)
    U = train.copy()
    base = []

    def add(bk):
        nonlocal U
        base.append(bk)
        inS[bk] = inS[N - bk] = True
        score[bk] = -1e18
        cov = inS[(U - bk) % N] | inS[(U + bk) % N]
        U = U[~cov]
        np.add.at(score, canon(U - bk, N), 1)
        np.add.at(score, canon(U + bk, N), 1)
        score[np.array(base)] = -1e18

    def gain(c):
        inS[c] = inS[N - c] = True
        g = int((inS[(U - c) % N] | inS[(U + c) % N]).sum())
        inS[c] = inS[N - c] = False
        # pairs not using c cannot newly cover U, so g counts targets covered via c
        return g

    add(rng.randrange(1, H + 1))
    while len(base) < L:
        if score.max() <= 0:
            while True:
                c = rng.randrange(1, H + 1)
                if not inS[c]:
                    break
            add(c)
            continue
        for _ in range(200):
            c = int(np.argmax(score))
            score[c] = gain(c)
            if int(np.argmax(score)) == c:
                break
        add(c)
    return np.array(base, dtype=np.int64)


def features(xs, ys, p):
    bits = p.bit_length()
    cols = [((xs >> i) & 1) for i in range(bits)]
    for q in QS:
        for r in range(q):
            cols.append((xs % q == r).astype(np.int64))
    sq = np.zeros(p, dtype=np.int8)
    a = np.arange(p, dtype=np.int64)
    sq[(a * a) % p] = 1
    cols.append(sq[xs % p])
    cols.append(sq[(xs + 1) % p])
    cols.append(np.minimum(ys, p - ys) & 1)
    F = np.stack(cols, 1).astype(float)
    mu, sd = F.mean(0), F.std(0) + 1e-9
    return (F - mu) / sd


def fit_l1(F, y, lam=1e-3, steps=300, lr=0.5):
    w = np.zeros(F.shape[1]); b = 0.0
    for _ in range(steps):
        z = F @ w + b
        pr = 1 / (1 + np.exp(-z))
        g = F.T @ (pr - y) / len(y)
        w -= lr * g
        b -= lr * float((pr - y).mean())
        w = np.sign(w) * np.maximum(np.abs(w) - lr * lam, 0)
    return w, b


def auc(s, y):
    order = np.argsort(s)
    ranks = np.empty(len(s)); ranks[order] = np.arange(1, len(s) + 1)
    n1 = y.sum(); n0 = len(y) - n1
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def train_select(Fpool, pos_idx, neg_idx, L, rng):
    idx = np.concatenate([pos_idx, neg_idx])
    y = np.concatenate([np.ones(len(pos_idx)), np.zeros(len(neg_idx))])
    perm = np.array(rng.sample(range(len(idx)), len(idx)))
    cut = int(0.8 * len(idx))
    tr, te = perm[:cut], perm[cut:]
    w, b = fit_l1(Fpool[idx[tr]], y[tr])
    a = auc(Fpool[idx[te]] @ w, y[te]) if 0 < y[te].sum() < len(te) else float("nan")
    w, b = fit_l1(Fpool[idx], y)
    s = Fpool @ w
    s[0] = -1e18
    sel = np.argsort(-s, kind="stable")[:L]
    return sel.astype(np.int64), a, int((np.abs(w) > 1e-9).sum())


units = []
for bits in BITS:
    for ci, (p, a, b, N, f, sq) in enumerate(find_curves(bits, CURVES)):
        X, Y = multiples(p, a, b, N, f, sq)
        H = (N - 1) // 2
        kpool = np.arange(H + 1)
        Fpool = features(np.maximum(X[kpool], 0), np.maximum(Y[kpool], 0), p)
        for mu in MUS:
            L = int(round(math.sqrt(mu * N / 2)))
            for seed in range(SEEDS):
                rng = random.Random(1000 * bits + 100 * ci + 10 * seed + int(mu * 4))
                train = np.array(rng.sample(range(1, N), N // 16), dtype=np.int64)
                tmask = np.zeros(N, dtype=bool); tmask[train] = True
                arms = {}
                rb = np.array(rng.sample(range(1, H + 1), L)); arms["random"] = rb
                arms["random2"] = np.array(rng.sample(range(1, H + 1), L))
                g = greedy(N, L, train, rng); arms["greedy"] = g
                gset = set(g.tolist())
                neg = np.array([k for k in rng.sample(range(1, H + 1), 3 * L) if k not in gset][:L])
                arms["learned_greedy"], auc_g, nz_g = train_select(Fpool, g, neg, L, rng)
                # decomposition-success label vs reference random base of size 4L
                ref = np.array(rng.sample(range(1, H + 1), min(4 * L, H)))
                rind = np.zeros(N); rind[ref] = 1; rind[N - ref] = 1
                tind = tmask.astype(float)
                # part[c] = sum_r t[r] (rind[r-c] + rind[r+c])
                Ft, Fr = np.fft.rfft(tind), np.fft.rfft(rind)
                corr = np.rint(np.fft.irfft(Ft * np.conj(Fr), n=N))  # sum_r t[r] rind[r-c]
                part = corr[kpool] + corr[(-kpool) % N]
                med = np.median(part[1:])
                samp = np.array(rng.sample(range(1, H + 1), min(4000, H)))
                hi = samp[part[samp] > med]; lo = samp[part[samp] <= med]
                n = min(len(hi), len(lo))
                arms["learned_success"], auc_s, nz_s = train_select(Fpool, hi[:n], lo[:n], L, rng)
                arms["known_log"] = np.arange(1, L + 1)
                # PC2: learnable planted label
                xs_s = np.maximum(X[samp], 0)
                lab = (xs_s % 7 == 3)
                _, auc_pc2, _ = train_select(Fpool, samp[lab], samp[~lab][: max(lab.sum(), 1)], L, rng)
                rec = {"bits": bits, "p": p, "a": a, "b": b, "N": N, "curve": ci, "mu": mu,
                       "L": L, "seed": seed, "T_train": int(len(train)),
                       "auc_learned_greedy": auc_g, "auc_learned_success": auc_s,
                       "auc_pc2_mod7": auc_pc2, "nnz_greedy": nz_g, "nnz_success": nz_s,
                       "arms": {k: metrics(v, N, tmask) for k, v in arms.items()}}
                assert all(m["sum_count_all_r"] == m["S_size"] ** 2 for m in rec["arms"].values())
                units.append(rec)
                print(f"[{time.time()-T0:7.1f}s] bits={bits} c={ci} mu={mu} s={seed} L={L} "
                      + " ".join(f"{k}={m['held_success']:.4f}/{m['train_success']:.3f}" for k, m in rec["arms"].items())
                      + f" aucG={auc_g:.3f} aucS={auc_s:.3f} pc2={auc_pc2:.3f}", flush=True)
                json.dump({"protocol_sha256": open("PROTOCOL.sha256").read().split()[0],
                           "elapsed_s": time.time() - T0, "units": units}, open(OUT, "w"), indent=1)
print("done", time.time() - T0)
