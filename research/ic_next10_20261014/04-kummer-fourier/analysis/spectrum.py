"""Item 4: target-conditioned Fourier spectrum of Kummer-tower decompositions.

See ../PROTOCOL.md.  Usage (from the repository root):

  python3 -I research/ic_next10_20261014/04-kummer-fourier/analysis/spectrum.py \
      cell --size 20 --m 2 --curve 0 [--nperm 200] [--nrand 100] [--planted]
  python3 -I .../spectrum.py c4cell            # continuity cell (C4 curve, t = 8)
  python3 -I .../spectrum.py verify            # arithmetic cross-check

Writes raw/cell_<size>_m<m>_c<curve>.json.
"""
import argparse
import itertools
import json
import os
import subprocess
import sys
import time

import numpy as np
import scipy.fft as sfft

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
ITEM = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BIN = os.path.join(ROOT, "target", "release", "examples", "ic_next10_04_kummer_fourier")
SCRATCH = os.environ.get("IC4_SCRATCH", "/tmp/ic4")
PRIMES = {20: 1073153, 22: 4206593, 24: 16801793}
SIZE_INDEX = {20: 0, 22: 1, 24: 2}
T_FOR_M = {2: 11, 3: 7}
C4 = dict(p=946177, a=304973, b=530845)
PLANT = {2: (1, 3), 3: (1, 3, 5)}
RHOS = [1.0, 0.3, 0.1, 0.03]
TOPK = (1, 10, 100)


def curve_params(size, idx):
    p = PRIMES[size]
    rng = np.random.default_rng(20261014000 + 100 * SIZE_INDEX[size] + idx)
    while True:
        a, b = (int(v) for v in rng.integers(0, p, size=2))
        if (4 * a**3 + 27 * b**2) % p:
            return p, a, b


def enumerate_base(p, a, b, t, m, base, seed, tag):
    os.makedirs(SCRATCH, exist_ok=True)
    out = os.path.join(SCRATCH, f"{tag}.bin")
    r = subprocess.run(
        [BIN, "--p", str(p), "--a", str(a), "--b", str(b), "--t", str(t), "--m", str(m),
         "--base", base, "--seed", str(seed), "--out", out],
        check=True, capture_output=True, text=True)
    info = json.loads(r.stdout)
    data = load(out)
    os.remove(out)
    return info, data


def load(path):
    raw = np.fromfile(path, dtype=np.uint8)
    m, n, nl = np.frombuffer(raw[:12], dtype="<u4")
    off = 12
    labels = np.frombuffer(raw[off:off + 2 * nl], dtype="<u2").astype(np.int64)
    off += 2 * nl
    (T,) = np.frombuffer(raw[off:off + 4], dtype="<u4")
    off += 4
    Ks = np.frombuffer(raw[off:off + 2 * int(T)], dtype="<u2").astype(np.int64)
    off += 2 * int(T)
    reps = np.frombuffer(raw[off:], dtype="<u2").astype(np.int32).reshape(-1, int(m))
    assert len(reps) == Ks.sum()
    return dict(m=int(m), N=int(n), labels=labels, Ks=Ks, reps=reps)


def cross_pairs(Ks, reps):
    """Unordered cross-orbit pairs within each target: arrays A, B (P x m)."""
    m = reps.shape[1]
    offs = np.concatenate([[0], np.cumsum(Ks)[:-1]])
    As, Bs = [], []
    for k in np.unique(Ks):
        sel = offs[Ks == k]
        idx = sel[:, None] + np.arange(k)[None, :]
        arr = reps[idx]  # (n_k, k, m)
        i0, i1 = np.triu_indices(k, 1)
        As.append(arr[:, i0, :].reshape(-1, m))
        Bs.append(arr[:, i1, :].reshape(-1, m))
    return np.concatenate(As), np.concatenate(Bs)


_CANON = {}


def canon(N, m):
    """Canonical nonzero character classes (up to S_m and c -> -c)."""
    key = (N, m)
    if key in _CANON:
        return _CANON[key]
    def keep_canon(cs):
        cs = cs[np.any(cs != 0, axis=1)]
        neg = np.sort((-cs) % N, axis=1)
        # lexicographic c <= sorted(-c)
        keep = np.ones(len(cs), bool)
        decided = np.zeros(len(cs), bool)
        for j in range(m):
            lt = (cs[:, j] < neg[:, j]) & ~decided
            gt = (cs[:, j] > neg[:, j]) & ~decided
            keep[gt] = False
            decided |= lt | gt
        return cs[keep]

    if m == 2:
        cs = keep_canon(np.stack(np.triu_indices(N), axis=1).astype(np.int32))
    else:
        i2, i3 = np.triu_indices(N)
        parts = []
        for c1 in range(N):
            s = i2 >= c1
            parts.append(keep_canon(np.stack([np.full(s.sum(), c1), i2[s], i3[s]], axis=1).astype(np.int32)))
        cs = np.concatenate(parts)
        del parts
    # index into rfftn output (last axis halved): use c or -c
    flip = cs[:, -1] > N // 2
    use = np.where(flip[:, None], (-cs) % N, cs).astype(np.int64)
    shape = (N,) * (m - 1) + (N // 2 + 1,)
    idx = np.ravel_multi_index(tuple(use.T), shape)
    # |Stab_{S_m}(c)|
    if m == 2:
        stab = np.where(cs[:, 0] == cs[:, 1], 2, 1)
    else:
        e01, e12 = cs[:, 0] == cs[:, 1], cs[:, 1] == cs[:, 2]
        stab = np.where(e01 & e12, 6, np.where(e01 | e12, 2, 1))
    real = np.all((cs == 0) | (cs == N // 2), axis=1)
    norm = stab * np.where(real, np.sqrt(2.0), 1.0)
    nz = (cs != 0).sum(axis=1)
    code = cs.astype(np.int64) @ (N ** np.arange(m)).astype(np.int64)
    order = np.argsort(code)
    out = dict(cs=cs.astype(np.int32), idx=idx, norm=norm.astype(np.float64), nz=nz.astype(np.int8),
               code_sorted=code[order], order=order)
    _CANON[key] = out
    return out


def class_of(c, N):
    c = np.sort(np.asarray(c) % N)
    n = np.sort((-c) % N)
    return tuple(int(v) for v in min(tuple(c), tuple(n)))


def _summary(Z, cs, prefix):
    Z2 = Z * Z
    kk = max(TOPK)
    top = np.argpartition(Z2, -kk)[-kk:]
    order = top[np.argsort(Z2[top])[::-1]]
    tot = Z2.sum()
    i = order[0]
    return {
        f"{prefix}max": float(abs(Z[i])),
        f"{prefix}max_signed": float(Z[i]),
        f"{prefix}argmax": [int(v) for v in cs[i]],
        f"{prefix}top10": [[int(v) for v in cs[o]] + [round(float(Z[o]), 3)] for o in order[:10]],
        **{f"{prefix}M{k}": float(Z2[order[:k]].sum() / tot) for k in TOPK},
    }


def spectrum_stats(A, B, N, m, extra_classes=()):
    """Normalised cross-orbit spectrum on canonical classes and its statistics.

    Z  : frozen normalisation Ahat / (|Stab| m! sqrt(2P) [sqrt 2 if real]).
    Zt : Z standardised within its zero-pattern type (number of nonzero
         coordinates), by the mean and sd of that type in the same spectrum
         (removes the label-invariant offset that shared summands put on
         characters with zero coordinates; see README, deviations).
    """
    P = len(A)
    size = N**m
    H = np.zeros(size, dtype=np.float32)
    strides = (N ** np.arange(m - 1, -1, -1)).astype(np.int64)
    for rho in itertools.permutations(range(m)):
        lin = np.zeros(len(A), dtype=np.int64)
        for j in range(m):
            lin += ((A[:, j].astype(np.int64) - B[:, rho[j]]) % N) * strides[j]
        u, cnt = np.unique(lin, return_counts=True)
        H[u] += cnt
        del u, cnt
        del lin
    H = H.reshape((N,) * m)
    Hs = np.zeros_like(H)
    for s in itertools.permutations(range(m)):
        Hs += H.transpose(s)
    del H
    F = sfft.rfftn(Hs, workers=1)
    del Hs
    cn = canon(N, m)
    fact = np.prod(range(1, m + 1))
    Ahat = 2.0 * F.reshape(-1)[cn["idx"]].real.astype(np.float64)
    del F
    Z = Ahat / (cn["norm"] * fact * np.sqrt(2.0 * P))
    Zt = np.empty_like(Z)
    types = {}
    for ty in range(1, m + 1):
        sel = cn["nz"] == ty
        mu, sd = Z[sel].mean(), Z[sel].std()
        Zt[sel] = (Z[sel] - mu) / sd
        types[str(ty)] = dict(n=int(sel.sum()), mean=float(mu), sd=float(sd))
    res = dict(P=int(P), types=types, mean=float(Z.mean()), sd=float(Z.std()))
    res.update(_summary(Z, cn["cs"], "Z"))
    res.update(_summary(Zt, cn["cs"], "Zt"))
    # backwards-compatible names used by the decision code
    res["maxZ"], res["argmax"] = res["Zmax"], res["Zargmax"]
    for k in TOPK:
        res[f"M{k}"] = res[f"ZM{k}"]
    want = [class_of([N >> k] * m, N) for k in range(1, int(np.log2(N)) + 1)]
    extra = [class_of(c, N) for c in extra_classes]
    def look(w, arr):
        c = int(np.dot(w, N ** np.arange(m)))
        h = cn["order"][np.searchsorted(cn["code_sorted"], c)]
        return float(arr[h])
    res["Z_C4_family"] = {str(w): [look(w, Z), look(w, Zt)] for w in want}
    res["Z_C4max"] = max(abs(look(w, Z)) for w in want)
    res["Zt_C4max"] = max(abs(look(w, Zt)) for w in want)
    res["Z_extra"] = {str(w): [look(w, Z), look(w, Zt)] for w in extra}
    res["hist"] = np.histogram(Zt, bins=np.linspace(-8, 8, 161))[0].tolist()
    return res


def marginal_stat(labels, N):
    g = np.zeros(N)
    g[labels] = 1
    n1 = g.sum()
    F = np.abs(np.fft.fft(g))[1:]
    return float(F.max() / np.sqrt(n1 * (1 - n1 / N)))


def relabel(A, B, perm):
    return perm[A], perm[B]


def planted(Ks, N, m, rho, cstar, rng):
    tot = int(Ks.sum())
    reps = rng.integers(0, N, size=(tot, m))
    s = np.repeat(rng.integers(0, N, size=len(Ks)), Ks)
    on = rng.random(tot) < rho
    cm_inv = pow(cstar[-1], -1, N)
    partial = (reps[:, :-1] * np.array(cstar[:-1])).sum(axis=1)
    reps[on, -1] = ((s[on] - partial[on]) * cm_inv) % N
    return reps.astype(np.int32)


def summarize_null(vals):
    v = np.asarray(vals, dtype=float)
    return dict(n=len(v), q001=float(np.quantile(v, 0.001)), q50=float(np.quantile(v, 0.5)),
                q999=float(np.quantile(v, 0.999)), max=float(v.max()), min=float(v.min()))


def run_cell(p, a, b, t, m, nperm, nrand, do_planted, seed, name, seedoff=0):
    t0 = time.time()
    N = 1 << t
    info, data = enumerate_base(p, a, b, t, m, "kummer", 0, name + "_k")
    A, B = cross_pairs(data["Ks"], data["reps"])
    obs = spectrum_stats(A, B, N, m)
    obs["marginal"] = marginal_stat(data["labels"], N)
    rng = np.random.default_rng(seed + 7919 * seedoff)
    if nperm == 0 and not do_planted:
        del A, B
    perm_rows = []
    for _ in range(nperm):
        perm = rng.permutation(N).astype(np.int32)
        A2, B2 = relabel(A, B, perm)
        r = spectrum_stats(A2, B2, N, m)
        r["marginal"] = marginal_stat(perm[data["labels"]], N)
        r.pop("hist")
        perm_rows.append(r)
    rand_rows = []
    for j in range(nrand):
        rinfo, rdata = enumerate_base(p, a, b, t, m, "random", seed * 1000 + 100 * seedoff + j, name + "_r")
        RA, RB = cross_pairs(rdata["Ks"], rdata["reps"])
        r = spectrum_stats(RA, RB, N, m)
        r["marginal"] = marginal_stat(rdata["labels"], N)
        r.pop("hist")
        r["info"] = rinfo
        rand_rows.append(r)
    keys = ["maxZ", "Ztmax", "Z_C4max", "Zt_C4max", "M1", "M10", "M100", "ZtM1", "ZtM10", "ZtM100", "marginal"]
    perm_null = {k: summarize_null([r[k] for r in perm_rows]) for k in keys} if perm_rows else {}
    rand_null = {k: summarize_null([r[k] for r in rand_rows]) for k in keys} if rand_rows else {}
    out = dict(name=name, p=p, a=a, b=b, t=t, N=N, m=m, seed=seed, info=info, observed=obs,
               perm_null=perm_null, rand_null=rand_null,
               perm_values={k: [r[k] for r in perm_rows] for k in keys},
               rand_values={k: [r[k] for r in rand_rows] for k in keys},
               rand_P=[r["P"] for r in rand_rows])
    # Decision per cell (recomputed from merged values by report.py).
    th = perm_null["maxZ"]["q999"] if perm_rows else float("inf")
    tht = perm_null["Ztmax"]["q999"] if perm_rows else float("inf")
    if perm_rows:
        decide(out, obs, perm_null, keys, nperm, nrand, rand_rows)
    if do_planted:
        cstar = PLANT[m]
        prow = []
        for rho in RHOS:
            prng = np.random.default_rng(seed + int(rho * 1000))
            reps = planted(data["Ks"], N, m, rho, cstar, prng)
            PA, PB = cross_pairs(data["Ks"], reps)
            r = spectrum_stats(PA, PB, N, m, extra_classes=[cstar])
            r.pop("hist")
            zc = list(r["Z_extra"].values())[0]
            prow.append(dict(rho=rho, maxZ=r["maxZ"], Ztmax=r["Ztmax"], Ztargmax=r["Ztargmax"],
                             detected_Zt=bool(r["Ztmax"] > tht), argmax=r["argmax"], Z_cstar=zc,
                             argmax_in_cstar_subgroup=tuple(r["argmax"]) in {class_of([j * c for c in cstar], N) for j in range(1, N)},
                             detected=bool(r["maxZ"] > th),
                             M1=r["M1"], M100=r["M100"], ZtM1=r["ZtM1"], ZtM100=r["ZtM100"]))
        out["planted"] = dict(cstar=list(cstar), threshold=th, threshold_Zt=tht, rows=prow)
    out["seconds"] = time.time() - t0
    return out


def decide(out, obs, perm_null, keys, nperm, nrand, rand_rows):
    th = perm_null["maxZ"]["q999"]
    tht = perm_null["Ztmax"]["q999"]
    out["exceed_maxZ"] = bool(obs["maxZ"] > th)
    out["exceed_Ztmax"] = bool(obs["Ztmax"] > tht)
    out["exceed_C4"] = bool(obs["Z_C4max"] > perm_null["Z_C4max"]["q999"] or obs["Zt_C4max"] > perm_null["Zt_C4max"]["q999"])
    out["exceed_marginal"] = bool(obs["marginal"] > perm_null["marginal"]["q999"])
    out["mass_outside_band"] = {k: bool(not (perm_null[k]["q001"] <= obs[k] <= perm_null[k]["q999"]))
                                for k in ["M1", "M10", "M100", "ZtM1", "ZtM10", "ZtM100"]}
    out["perm_p"] = {k: (1 + sum(v >= obs[k] for v in out["perm_values"][k])) / (1 + nperm) for k in keys}
    out["perm_p_maxZ"] = out["perm_p"]["maxZ"]
    if rand_rows:
        out["rand_p"] = {k: (1 + sum(v >= obs[k] for v in out["rand_values"][k])) / (1 + nrand) for k in keys}
        out["rand_rank_maxZ"] = out["rand_p"]["maxZ"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["cell", "c4cell", "verify"])
    ap.add_argument("--size", type=int)
    ap.add_argument("--m", type=int)
    ap.add_argument("--curve", type=int, default=0)
    ap.add_argument("--nperm", type=int, default=200)
    ap.add_argument("--nrand", type=int, default=100)
    ap.add_argument("--planted", action="store_true")
    ap.add_argument("--tag", default="")
    ap.add_argument("--seedoff", type=int, default=0, help="offset for the permutation-null RNG only")
    args = ap.parse_args()
    os.makedirs(os.path.join(ITEM, "raw"), exist_ok=True)
    if args.cmd == "verify":
        return verify()
    if args.cmd == "cell":
        p, a, b = curve_params(args.size, args.curve)
        t, m = T_FOR_M[args.m], args.m
        name = f"cell_{args.size}_m{m}_c{args.curve}{args.tag}"
        seed = 20261014 + 1000 * SIZE_INDEX[args.size] + 10 * m + args.curve
    else:
        p, a, b, t, m = C4["p"], C4["a"], C4["b"], 8, 3
        name = f"cell_c4_m3_t8{args.tag}"
        seed = 20261013
    res = run_cell(p, a, b, t, m, args.nperm, args.nrand, args.planted, seed, name, args.seedoff)
    with open(os.path.join(ITEM, "raw", name + ".json"), "w") as f:
        json.dump(res, f)
    o = res["observed"]
    pn = res["perm_null"]
    print(json.dumps(dict(name=name, P=o["P"], Z=round(o["maxZ"], 3), Zt=round(o["Ztmax"], 3), Ztarg=o["Ztargmax"],
                          th=pn.get("maxZ", {}).get("q999"), tht=pn.get("Ztmax", {}).get("q999"),
                          nperm=args.nperm, nrand=args.nrand,
                          planted=[(r["rho"], round(r["maxZ"], 2), round(r["Ztmax"], 2)) for r in res.get("planted", {}).get("rows", [])],
                          sec=round(res["seconds"]))))


def verify():
    """Cross-check: every orbit of a target reaches the same abscissa (pure Python EC)."""
    p, a, b = C4["p"], C4["a"], C4["b"]
    t, m, N = 8, 3, 256
    info, data = enumerate_base(p, a, b, t, m, "kummer", 0, "verify")
    z = None
    e = (p - 1) >> t
    for h in range(2, 1000):
        z = pow(h, e, p)
        if pow(z, 1 << (t - 1), p) != 1:
            break
    def sqrt(n):
        r = int(subprocess.run(["gp", "-q", "-f"], input=f"print(lift(sqrt(Mod({n},{p}))))", capture_output=True, text=True).stdout)
        return r
    def add(P, Q):
        if P is None: return Q
        if Q is None: return P
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2:
            if (y1 + y2) % p == 0: return None
            l = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
        else:
            l = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (l * l - x1 - x2) % p
        return (x3, (l * (x1 - x3) - y1) % p)
    ys = {}
    def pt(e_):
        x = pow(z, int(e_), p)
        if x not in ys: ys[x] = sqrt((x**3 + a * x + b) % p)
        return (x, ys[x])
    offs = np.concatenate([[0], np.cumsum(data["Ks"])[:-1]])
    rng = np.random.default_rng(1)
    ok = 0
    for ti in rng.choice(len(data["Ks"]), 30, replace=False):
        xsets = []
        for o in data["reps"][offs[ti]:offs[ti] + data["Ks"][ti]]:
            P = [pt(v) for v in o]
            xs = set()
            for sg in itertools.product([1, -1], repeat=3):
                Q = [(x, y if s == 1 else (-y) % p) for (x, y), s in zip(P, sg)]
                if any(Q[i][0] == Q[j][0] and Q[i][1] != Q[j][1] for i in range(3) for j in range(i + 1, 3)):
                    continue
                R = add(add(Q[0], Q[1]), Q[2])
                if R: xs.add(R[0])
            xsets.append(xs)
        common = set.intersection(*xsets)
        ok += len(common) >= 1
    print(json.dumps(dict(verify_targets=30, consistent=ok, info=info)))


if __name__ == "__main__":
    main()
