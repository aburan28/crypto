"""Arm G5c: canonical genus-5 complete intersections, sampled nets.

usage: python3 -I run_ci.py <log2 stage-1 nets> <log2 subsample nets> <out.json>
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fam_ci import stage1_tables, stage1, random_nets, full, gauss_binom, gl_order
from lpoly import targets, weil_polys, polymul, counts_from_charpoly, charpoly_from_counts, divides, dim, is_weil_valid

SEED = 20261014
L_B7 = 18  # B7 (dim-3, n=7 factor) filter only on the first 2^18 nets


def full_abort(net):
    from fam_ci import count_and_sing
    N = []
    for k in range(1, 6):
        n, s = count_and_sing(net, k)
        if s:
            return N, True
        N.append(n)
    return N, False


def main():
    l1, l2, outp = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    T = targets()
    g = 5
    cand = {}
    cand4 = set()
    for name, P in T.items():
        d = dim(P)
        if d > g:
            continue
        if name.startswith("A3"):
            continue  # positive control handled on the subsample
        for Qc in weil_polys(g - d):
            prod = polymul(P, Qc)
            cand.setdefault(tuple(counts_from_charpoly(prod, 3)), []).append(name)
            cand4.add(tuple(counts_from_charpoly(prod, 4)))
    tabs = stage1_tables()
    rng = np.random.default_rng(SEED)
    B = 8192
    n1 = 0
    idx = 0
    survivors = []
    sub = []
    t0 = time.time()
    while n1 < (1 << l1):
        nets = random_nets(rng, B)[: (1 << l1) - n1]
        if n1 < (1 << l2):
            sub.extend([tuple(int(x) for x in r) for r in nets[: (1 << l2) - n1]])
        c = stage1(nets, tabs)
        for r, cc in zip(nets, c):
            names = cand.get(tuple(int(x) for x in cc))
            if names and (idx < (1 << L_B7) or any(nm[:2] in ("A5", "W5") for nm in names)):
                survivors.append(tuple(int(x) for x in r))
            idx += 1
        n1 += len(nets)
        if (n1 // B) % 32 == 0:
            print(f"stage1 {n1} nets, {len(survivors)} survivors, {time.time()-t0:.0f}s", flush=True)
    res_surv = []
    from fam_ci import field, VALS
    field(4)
    V4 = VALS[4]
    n_pass4 = 0
    for net in survivors:
        z = None
        for qm in net:
            supp = [m for m in range(15) if (qm >> m) & 1]
            v = np.bitwise_xor.reduce(V4[:, supp], axis=1)
            z = v if z is None else (z | v)
        n4 = int((z == 0).sum())
        c3 = stage1(np.array([net]), tabs)[0]
        if tuple(int(x) for x in c3) + (n4,) not in cand4:
            continue
        n_pass4 += 1
        N, sing = full_abort(net)
        rec = {"net": net, "N": N, "singular": sing}
        if not sing:
            P = charpoly_from_counts(N)
            rec["hits"] = [nm for nm, t in T.items() if dim(t) <= g and P and divides(t, P)]
        res_surv.append(rec)
    print(f"stage2 survivors done {time.time()-t0:.0f}s", flush=True)
    hist = {}
    nsm = 0
    for i, net in enumerate(sub):
        N, sing = full_abort(net)
        if not sing:
            nsm += 1
            key = ",".join(map(str, N))
            hist[key] = hist.get(key, 0) + 1
        if i % 2048 == 0:
            print(f"sub {i} smooth {nsm} {time.time()-t0:.0f}s", flush=True)
    out = {
        "seed": SEED, "stage1_nets": n1, "subsample_nets": len(sub),
        "subsample_smooth": nsm, "subsample_hist": hist,
        "candidate_triples": len(cand), "b7_filter_nets": 1 << L_B7, "survivors": res_surv, "stage1_survivors": len(survivors), "pass_N4": n_pass4,
        "total_nets": gauss_binom(15, 3), "gl5_order": gl_order(5),
        "seconds": time.time() - t0,
    }
    json.dump(out, open(outp, "w"), indent=1)


if __name__ == "__main__":
    main()
