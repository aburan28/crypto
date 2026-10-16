#!/usr/bin/env python3
"""Build the tall-isogeny-volcano dataset.

Every volcano is constructed from a CM Frobenius pi = x + f*omega with
f = c * ell^k, explored exactly from its crater by Phi_ell-root BFS, and checked
fail-closed against independent theory before any row is written:

  V1 level sizes            == class numbers h(D_K (c ell^i)^2)
  V2 vertex degree          == ell + 1 rational isogenies (non-floor) or 1 (floor)
  V3 group structure        == O_d/(pi - 1), measured from point orders
  V4 point count            == p + 1 - t on a sample (SEA, independent of model choice)
  V5 blind descent          == BFS floor, level and path length on random vertices

Outputs (gzip CSV) in --out:  volcanoes, vertices, edges, descents.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import multiprocessing as mp
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from volcano_descent import (  # noqa: E402
    Instance, class_number, crater_js, curve_for_j_fast, decompose, descend_to_floor,
    expected_level_counts, explore_from_crater, find_instance, level_counts,
    group_structure, neighbors, pari, predicted_group, v_ell,
)

SIG_PRIMES = (2, 3, 5, 7, 11, 13)
FULL_CAP = 4200      # store every vertex up to this size
PER_LEVEL_CAP = 300  # otherwise sample this many per level
DESCENT_TRIALS = 24


def campaign():
    """Deterministic list of (ell, D, c, k, m) tall-volcano configurations."""
    cfg = []
    def add(ell, D, c, k, m):
        f = c * ell ** k
        if D % 8 == 1 and f % 2 == 1:      # p = N(pi) would be even
            c *= 2
        cfg.append((ell, D, c, k, m))
    for D in (-3, -4, -7, -8, -11, -15, -19, -20, -23, -24):
        for k, m in ((12, 0), (12, 6), (12, 12)):
            add(2, D, 1, k, m)
    add(2, -4, 1, 15, 0); add(2, -7, 1, 15, 7); add(2, -3, 1, 14, 14)
    for D in (-3, -4, -7, -8, -11, -20):
        for k, m in ((7, 0), (7, 3), (7, 7)):
            add(3, D, 1, k, m)
    add(3, -8, 1, 9, 4)
    for D in (-3, -4, -7, -8, -11, -19):
        for k, m in ((5, 0), (5, 2), (5, 5)):
            add(5, D, 1, k, m)
    for D in (-3, -4, -8, -11):
        for k, m in ((4, 0), (4, 2), (4, 4)):
            add(7, D, 1, k, m)
    for D in (-3, -4, -8):
        add(11, D, 1, 3, 0); add(11, D, 1, 3, 3)
        add(13, D, 1, 3, 0); add(13, D, 1, 3, 3)
    # non-maximal crater orders (c > 1)
    add(2, -7, 3, 10, 0); add(2, -4, 3, 10, 5); add(2, -8, 5, 9, 0)
    add(3, -8, 2, 6, 0); add(3, -4, 2, 6, 3); add(3, -11, 2, 6, 0)
    add(5, -4, 3, 4, 0); add(5, -8, 3, 4, 2)
    return cfg


def target_bits(ell, D, c, k, m, rbits=28):
    f = c * ell ** k
    smooth = int(2 * m * math.log2(ell)) + 12
    need = int(math.log2(-D * f * f)) + 3
    return max(rbits + smooth, need), smooth


def vid(ell, D, c, k, m, seed):
    return f"v_l{ell}_D{abs(D)}_c{c}_k{k}_m{m}_s{seed}"


def signature(j, p):
    out = []
    for q in SIG_PRIMES:
        _, tot = neighbors(q, j, p)
        out.append(tot)
    return out


def build_one(args):
    cfg, seed = args
    ell, D, c, k, m = cfg
    t0 = time.time()
    rng = random.Random(f"{seed}:{cfg}")
    pari.setrand(int.from_bytes(hashlib.sha256(f"{seed}:{cfg}".encode()).digest()[:4], "big") | 1)
    bits, smooth = target_bits(ell, D, c, k, m)
    inst = None
    for extra in range(0, 12, 2):
        try:
            inst = find_instance(ell, D, c, k, m, bits + extra, rng, max_tries=60000, big_prime_slack=smooth + extra)
            break
        except (RuntimeError, ValueError):
            continue
    vol_id = vid(ell, D, c, k, m, seed)
    if inst is None:
        return dict(vol_id=vol_id, error="no_instance", cfg=cfg)
    p, t, N = inst.p, inst.t, inst.N
    V = explore_from_crater(inst)
    got, exp = level_counts(V), expected_level_counts(inst)
    checks = dict(V1_level_counts=(got == exp))
    # adjacency by type
    nlevel = V.level
    up, ndown, nhor, mult = {}, {}, {}, {}
    edges = []
    for src, dst, e, tot in V.edges:
        ls, ld = nlevel[src], nlevel[dst]
        typ = "horizontal" if ld == ls else ("down" if ld == ls + 1 else "up")
        if typ == "up":
            up[src] = dst
        elif typ == "down":
            ndown[src] = ndown.get(src, 0) + 1
        else:
            nhor[src] = nhor.get(src, 0) + e
        edges.append((src, dst, e, typ))
    deg_ok = all(V.totals[j] == (ell + 1 if nlevel[j] < k else 1) for j in V.level)
    checks["V2_degrees"] = deg_ok
    total = len(V.level)
    full = total <= FULL_CAP
    by_level = {}
    for j, l in V.level.items():
        by_level.setdefault(l, []).append(j)
    sel = []
    for l in sorted(by_level):
        js = sorted(by_level[l])
        if not full and len(js) > PER_LEVEL_CAP:
            js = sorted(rng.sample(js, PER_LEVEL_CAP))
        sel.extend(js)
    Nfac_m = pari.factor(N)
    verts = []
    grp_bad = 0
    count_bad = 0
    for idx, j in enumerate(sel):
        a, b = curve_for_j_fast(j, p, t, rng)
        E = pari.ellinit([a, b], p)
        n1, n2 = group_structure(E, p, N, Nfac_m, rng)
        n1p, n2p = predicted_group(p, t, D, c * ell ** nlevel[j])
        gm = ((n1, n2) == (n1p, n2p))
        grp_bad += (not gm)
        if idx % 97 == 0:   # V4 on a sample, with SEA
            count_bad += int(pari.ellcard(E)) != N
        lv = nlevel[j]
        cond = c * ell ** lv
        sig = signature(j, p)
        verts.append(dict(
            vol_id=vol_id, j=j, a=a, b=b, level=lv, height=k - lv, end_conductor=cond,
            end_disc=D * cond * cond, end_class_number=class_number(D * cond * cond),
            comp=V.comp[j], up_j=up.get(j, ""), n_down=ndown.get(j, 0), n_horizontal=nhor.get(j, 0),
            ell_rational_mult=V.totals[j], n1=n1, n2=n2, group_matches_theory=int(gm),
            **{f"sig_{q}": s for q, s in zip(SIG_PRIMES, sig)},
            aut_order=(6 if j == 0 else 4 if j == 1728 % p else 2),
            sampled=int(not full),
        ))
    checks["V3_group_structure"] = (grp_bad == 0)
    checks["V4_point_count"] = (count_bad == 0)
    # V5 blind descent from random vertices (uses only j, p, t, ell)
    desc = []
    pool = list(V.level)
    for _ in range(DESCENT_TRIALS):
        j0 = rng.choice(pool)
        t1 = time.time()
        r = descend_to_floor(j0, p, t, ell, rng)
        ok = (nlevel[r["floor_j"]] == k and r["start_level"] == nlevel[j0]
              and len(r["path"]) - 1 == k - nlevel[j0])
        desc.append(dict(vol_id=vol_id, start_j=j0, true_level=nlevel[j0], claimed_level=r["start_level"],
                         steps=len(r["path"]) - 1, floor_j=r["floor_j"], success=int(ok),
                         seconds=round(time.time() - t1, 4)))
    checks["V5_blind_descent"] = all(d["success"] for d in desc)
    Nfac = pari.factor(N)
    nfac = " * ".join(f"{int(Nfac[i,0])}^{int(Nfac[i,1])}" if int(Nfac[i, 1]) > 1 else f"{int(Nfac[i,0])}"
                      for i in range(int(Nfac.nrows())))
    big = max(int(Nfac[i, 0]) for i in range(int(Nfac.nrows())))
    vol = dict(
        vol_id=vol_id, seed=seed, ell=ell, D_K=D, c=c, k=k, m=m, f=inst.f, p=p, p_bits=p.bit_length(),
        t=t, N=N, N_factorization=nfac, N_ell_valuation=v_ell(N, ell), largest_prime=big,
        largest_prime_bits=big.bit_length(), crater_size=got.get(0, 0), n_vertices=total,
        n_stored=len(verts), full=int(full), n_components=len(set(V.comp.values())),
        level_sizes=json.dumps(got), all_checks_pass=int(all(checks.values())),
        checks=json.dumps(checks), seconds=round(time.time() - t0, 2),
    )
    return dict(vol=vol, verts=verts, edges=[(vol_id,) + e for e in edges] if full else [], desc=desc)


def write_csv(path, rows):
    if not rows:
        return
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data")
    ap.add_argument("--seed", type=int, default=20261011)
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0, help="only the first N configurations")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    cfgs = campaign()
    if args.limit:
        cfgs = cfgs[:args.limit]
    jobs = [(c, args.seed) for c in cfgs]
    vols, verts, edges, descs, errs = [], [], [], [], []
    t0 = time.time()
    with mp.Pool(args.procs) as pool:
        for i, r in enumerate(pool.imap_unordered(build_one, jobs), 1):
            if "error" in r:
                errs.append(r)
                print(f"[{i}/{len(jobs)}] {r['vol_id']}: {r['error']}", flush=True)
                continue
            vols.append(r["vol"]); verts += r["verts"]; edges += r["edges"]; descs += r["desc"]
            v = r["vol"]
            print(f"[{i}/{len(jobs)}] {v['vol_id']} p={v['p_bits']}b |V|={v['n_vertices']} comps={v['n_components']} "
                  f"checks={'OK' if v['all_checks_pass'] else 'FAIL ' + v['checks']} {v['seconds']}s", flush=True)
    vols.sort(key=lambda r: r["vol_id"])
    verts.sort(key=lambda r: (r["vol_id"], r["level"], r["j"]))
    write_csv(os.path.join(args.out, "volcanoes.csv.gz"), vols)
    write_csv(os.path.join(args.out, "vertices.csv.gz"), verts)
    write_csv(os.path.join(args.out, "descents.csv.gz"),
              sorted(descs, key=lambda r: (r["vol_id"], r["start_j"])))
    if edges:
        with gzip.open(os.path.join(args.out, "edges.csv.gz"), "wt", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["vol_id", "src_j", "dst_j", "multiplicity", "type"])
            w.writerows(sorted(edges))
    json.dump(dict(errors=errs, n_volcanoes=len(vols), n_vertices=len(verts),
                   failed=[v["vol_id"] for v in vols if not v["all_checks_pass"]],
                   seconds=round(time.time() - t0, 1)),
              open(os.path.join(args.out, "build_summary.json"), "w"), indent=1)
    print(f"done: {len(vols)} volcanoes, {len(verts)} vertices, "
          f"{sum(1 for v in vols if not v['all_checks_pass'])} failed checks, {len(errs)} errors, "
          f"{time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
