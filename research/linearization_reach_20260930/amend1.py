"""Amendment 1 checks (A1-A3); see AMENDMENT-1.md.  Reuses the frozen lr.py unchanged."""
import json
import math
import random
import sys

import lr


def corrected_count(l, b):
    return l * (b + 1) + b - b * (b + 1) // 2


def rank_cells(rng):
    out = []
    for n in (23, 29, 31):
        for l in range(3, n):
            for b in range(1, min(l, 4) + 1):
                s = corrected_count(l, b)
                if not (n - 8 <= s <= n + 4):
                    continue
                hist = lr.rank_defect_cell(n, l, b, 200, rng)
                out.append(dict(n=n, l=l, b=b, corrected_count=s,
                                hist={str(k): v for k, v in sorted(hist.items())}))
    return out


def solve_family(F, V, W, idx, w0vec, t, l, b):
    cols, rhs = lr.columns(F, V, W, lr_combo(V, w0vec), t)
    N = len(cols)
    rk, cons, piv = lr.eliminate(lr.to_rows(F.n, cols, rhs), N)
    if not cons:
        return None
    free = [c for c in range(N) if c not in piv]

    def point(mask):
        sol = [0] * N
        for k, c in enumerate(free):
            sol[c] = (mask >> k) & 1
        for c in sorted(piv, reverse=True):
            r, v = piv[c], (piv[c] >> N) & 1
            for c2 in range(c + 1, N):
                if (r >> c2) & 1:
                    v ^= sol[c2]
            sol[c] = v
        return sol
    return N, rk, free, point


def lr_combo(V, vec):
    x = 0
    for i, v in enumerate(V):
        if vec[i]:
            x ^= v
    return x


def consistent(sol, l, b):
    cv, dv = sol[l * b:l * b + l], sol[l * b + l:]
    return all(sol[i * b + j] == (cv[i] & dv[j]) for i in range(l) for j in range(b))


def structured(sol0, l, b, idx, w0vec):
    """Resolve the family from one particular solution, in polynomial time.

    Kernel directions (AMENDMENT-1): delta_j flips c_{idx[j]}, d_j and c_i d_j for
    i in supp(w0); pair_jk flips c_{idx[j]} d_k and c_{idx[k]} d_j.  Choose each
    flip bit x_j from the constraints that involve x_j alone, then fix each pair
    bit from the c_{idx[j]} d_k product constraint.  Returns candidate solutions.
    """
    inW = set(idx)
    cand_x = []
    for j in range(b):
        ok = []
        for x in (0, 1):
            good = True
            dj = sol0[l * b + l + j] ^ x
            cij = sol0[l * b + idx[j]] ^ x
            if (sol0[idx[j] * b + j] ^ (x if w0vec[idx[j]] else 0)) != (cij & dj):
                good = False
            for i in range(l):
                if i in inW or not good:
                    continue
                prod = sol0[i * b + j] ^ (x if w0vec[i] else 0)
                if prod != (sol0[l * b + i] & dj):
                    good = False
            if good:
                ok.append(x)
        cand_x.append(ok)
    combos = [[]]
    for ok in cand_x:
        combos = [c + [x] for c in combos for x in ok]
    outs = []
    for xs in combos:
        sol = list(sol0)
        for j, x in enumerate(xs):
            if not x:
                continue
            sol[l * b + idx[j]] ^= 1
            sol[l * b + l + j] ^= 1
            for i in range(l):
                if w0vec[i]:
                    sol[i * b + j] ^= 1
        for j in range(b):
            for k in range(j + 1, b):
                want = sol[l * b + idx[j]] & sol[l * b + l + k]
                if sol[idx[j] * b + k] != want:
                    sol[idx[j] * b + k] ^= 1
                    sol[idx[k] * b + j] ^= 1
        if consistent(sol, l, b):
            outs.append(sol)
    return outs, len(combos)


def verified(F, sol, V, W, w0, t, X1, R, l, b):
    cv, dv = sol[l * b:l * b + l], sol[l * b + l:]
    X2 = lr_combo(V, cv)
    X3 = w0 ^ lr_combo(W, dv)
    if lr.S3(F, X2, X3, t) or lr.S3(F, X1, t, R[0]):
        return None
    pts = [lr.lift(F, x) for x in (X1, X2, X3)]
    if any(p is None for p in pts):
        return False
    for s in range(8):
        q = [pts[k] if (s >> k) & 1 else lr.neg(pts[k]) for k in range(3)]
        if lr.add(F, lr.add(F, q[0], q[1]), q[2]) == R:
            return True
    return False


def oracle_cells(rng):
    out = []
    for n, l, b, trials in [(13, 3, 3, 30000), (13, 4, 2, 30000), (17, 4, 3, 40000),
                            (17, 5, 2, 40000), (19, 5, 2, 40000), (23, 6, 2, 40000)]:
        F = lr.Field(n)
        st = dict(n=n, l=l, b=b, trials=trials, consistent_solves=0, family_points_hist={},
                  enum_success=0, struct_success=0, disagreements=0, algebra_mismatch=0,
                  max_combos=0, noncanonical_solves=0)
        for _ in range(trials):
            V = [rng.getrandbits(n) for _ in range(l)]
            while True:
                R = lr.lift(F, rng.getrandbits(n))
                if R:
                    break
            X1 = lr.rand_in(V, rng)
            A = F.sq(X1) ^ F.sq(R[0])
            B = F.mul(X1, R[0])
            C = F.mul(F.sq(X1), F.sq(R[0])) ^ 1
            idx = rng.sample(range(l), b)
            W = [V[i] for i in idx]
            w0vec = [rng.getrandbits(1) for _ in range(l)]
            w0 = lr_combo(V, w0vec)
            e_hit = s_hit = False
            for t in lr.quad_roots(F, A, B, C):
                fam = solve_family(F, V, W, idx, w0vec, t, l, b)
                if fam is None:
                    continue
                N, rk, free, point = fam
                st["consistent_solves"] += 1
                pts = [point(m) for m in range(1 << len(free))] if len(free) <= 12 else None
                e_this = False
                if pts is not None:
                    good = [p for p in pts if consistent(p, l, b)]
                    k = str(len(good))
                    st["family_points_hist"][k] = st["family_points_hist"].get(k, 0) + 1
                    for p in good:
                        v = verified(F, p, V, W, w0, t, X1, R, l, b)
                        if v is None:
                            st["algebra_mismatch"] += 1
                        e_this |= bool(v)
                e_hit |= e_this
                if N - rk == b * (b + 1) // 2:
                    outs, nc = structured(point(0), l, b, idx, w0vec)
                    st["max_combos"] = max(st["max_combos"], nc)
                    s_this = False
                    for p in outs:
                        v = verified(F, p, V, W, w0, t, X1, R, l, b)
                        if v is None:
                            st["algebra_mismatch"] += 1
                        s_this |= bool(v)
                    s_hit |= s_this
                    st["disagreements"] += (s_this != e_this) if pts is not None else 0
                else:
                    st["noncanonical_solves"] += 1
            st["enum_success"] += e_hit
            st["struct_success"] += s_hit
        print(json.dumps(st), flush=True)
        out.append(st)
    return out


def main(out_dir, seed):
    rng = random.Random(seed)
    res = dict(seed=seed, rank=rank_cells(rng), oracle=oracle_cells(rng))
    with open(f"{out_dir}/amend1.json", "w") as fh:
        json.dump(res, fh, indent=1)
    a1_lo = [c for c in res["rank"] if c["corrected_count"] <= c["n"] - 4]
    a1_hi = [c for c in res["rank"] if c["corrected_count"] >= c["n"] + 2]

    def frac(c, pred):
        h = {int(k): v for k, v in c["hist"].items()}
        return sum(v for k, v in h.items() if pred(k, c["b"])) / sum(h.values())
    score = dict(
        A1_low=dict(cells=len(a1_lo), min_frac=min(frac(c, lambda k, b: k == b * (b + 1) // 2) for c in a1_lo)),
        A1_high=dict(cells=len(a1_hi), min_frac=min(frac(c, lambda k, b: k > b * (b + 1) // 2) for c in a1_hi)),
    )
    score["A1"] = score["A1_low"]["min_frac"] >= 0.9 and score["A1_high"]["min_frac"] >= 0.99
    tot = sum(sum(c["family_points_hist"].values()) for c in res["oracle"])
    le2 = sum(v for c in res["oracle"] for k, v in c["family_points_hist"].items() if int(k) <= 2)
    score["A2"] = dict(solves=tot, frac_le2=le2 / tot, pass_=le2 / tot >= 0.99)
    score["A3"] = dict(disagreements=sum(c["disagreements"] for c in res["oracle"]),
                       noncanonical_solves=sum(c["noncanonical_solves"] for c in res["oracle"]),
                       consistent_solves=tot,
                       algebra_mismatch=sum(c["algebra_mismatch"] for c in res["oracle"]),
                       pass_=sum(c["disagreements"] for c in res["oracle"]) == 0)
    with open(f"{out_dir}/amend1_score.json", "w") as fh:
        json.dump(score, fh, indent=1)
    print(json.dumps(score, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
