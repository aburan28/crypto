#!/usr/bin/env python3
"""ECDLP measurements on dataset vertices, per volcano layer.

E1  generic layers   plain and negation-map rho + Pohlig-Hellman on vertices
                     stratified over every level of volcanoes with D_K not in {-3,-4}.
E2  automorphisms    D_K in {-3,-4}: only the crater curve has j = 0 / 1728 and
                     the extra units (mu_6 / mu_4); compare plain, neg and aut rho
                     at the crater against every other level.

Counted unit: affine group operations in the pure-Python reference solver.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import multiprocessing as mp
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ecdlp_rho import Curve, pick_base_point, solve_ph  # noqa: E402
from volcano_descent import pari  # noqa: E402


def read_gz(path):
    with gzip.open(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh))


def seed_of(*parts):
    return int.from_bytes(hashlib.sha256("|".join(map(str, parts)).encode()).digest()[:8], "big")


def task(args):
    vol, vx, modes, trials, tag = args
    p, N = int(vol["p"]), int(vol["N"])
    a, b = int(vx["a"]), int(vx["b"])
    Nfac = pari.factor(N)
    rows = []
    rng0 = random.Random(seed_of(tag, vol["vol_id"], vx["j"], "base"))
    G, o, nf = pick_base_point(p, a, b, N, Nfac, rng0)
    rmax = max(q for q, _ in nf)
    for mode in modes:
        for trial in range(trials):
            rng = random.Random(seed_of(tag, vol["vol_id"], vx["j"], mode, trial))
            C = Curve(p, a, b)
            x = rng.randrange(1, o)
            Q = C.mul(x, G)
            t0 = time.time()
            res = solve_ph(C, G, Q, nf, mode, rng)
            ok = (C.mul(res["x"], G) == Q)
            rows.append(dict(
                experiment=tag, vol_id=vol["vol_id"], ell=vol["ell"], D_K=vol["D_K"], k=vol["k"], m=vol["m"],
                p=p, j=vx["j"], level=int(vx["level"]), n1=int(vx["n1"]), n2=int(vx["n2"]),
                ord_G=o, ord_G_bits=o.bit_length(), r_max_bits=rmax.bit_length(), mode=mode, trial=trial,
                correct=int(ok), ops_total=res["ops_total"], ops_rho=res["ops_rho"],
                ops_rho_main=res["steps_max"], seconds=round(time.time() - t0, 3)))
    return rows


def pick(vol_rows, vert_rows, pred):
    vols = [v for v in vol_rows if pred(v) and v["all_checks_pass"] == "1"]
    by = {}
    for r in vert_rows:
        by.setdefault(r["vol_id"], []).append(r)
    return vols, by


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--experiment", choices=["E1", "E2"], required=True)
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--trials", type=int, default=16)
    ap.add_argument("--per-level", type=int, default=2)
    ap.add_argument("--max-volcanoes", type=int, default=6)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    vol_rows = read_gz(os.path.join(args.data, "volcanoes.csv.gz"))
    vert_rows = read_gz(os.path.join(args.data, "vertices.csv.gz"))
    rng = random.Random(20261011)
    jobs = []
    if args.experiment == "E1":
        # k = 12, ell = 2, generic D_K.  For each m in {0, 6, 12} take the two volcanoes with the
        # smallest largest prime (keeps the pure-Python solver fast; cost is normalised by sqrt(pi r / 2)).
        vols, by = pick(vol_rows, vert_rows,
                        lambda v: v["ell"] == "2" and v["k"] == "12" and v["D_K"] not in ("-3", "-4") and v["c"] == "1")
        chosen = []
        for m in ("0", "6", "12"):
            cand = sorted((v for v in vols if v["m"] == m), key=lambda v: (int(v["largest_prime_bits"]), v["vol_id"]))
            chosen += cand[:args.max_volcanoes // 3]
        for v in chosen:
            vs = by[v["vol_id"]]
            lv = {}
            for r in vs:
                lv.setdefault(int(r["level"]), []).append(r)
            for l in sorted(lv):
                for r in rng.sample(lv[l], min(args.per_level, len(lv[l]))):
                    jobs.append((v, r, ("plain", "neg"), args.trials, "E1"))
    else:
        # D_K in {-3, -4}: only the crater has j = 0 / 1728.  Four volcanoes (D_K, m), fixed in advance.
        vols, by = pick(vol_rows, vert_rows,
                        lambda v: v["ell"] == "2" and v["k"] == "12" and v["c"] == "1")
        wanted = [("-3", "0"), ("-4", "0"), ("-3", "6"), ("-4", "12")]
        chosen = [v for D, m in wanted for v in vols if v["D_K"] == D and v["m"] == m]
        for v in chosen:
            vs = by[v["vol_id"]]
            lv = {}
            for r in vs:
                lv.setdefault(int(r["level"]), []).append(r)
            for l in sorted(lv):
                n = 1 if l == 0 else args.per_level
                for r in rng.sample(lv[l], min(n, len(lv[l]))):
                    tr = args.trials * (6 if l == 0 else 1)   # the crater is one vertex: more trials
                    jobs.append((v, r, ("plain", "neg", "aut"), tr, "E2"))
    out = args.out or os.path.join(args.data, f"ecdlp_{args.experiment}.csv.gz")
    print(f"{args.experiment}: {len(jobs)} vertex tasks", flush=True)
    rows, t0 = [], time.time()
    with mp.Pool(args.procs) as pool:
        for i, r in enumerate(pool.imap_unordered(task, jobs), 1):
            rows += r
            if i % 10 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)} tasks, {len(rows)} solves, {time.time()-t0:.0f}s", flush=True)
    rows.sort(key=lambda r: (r["vol_id"], r["level"], r["j"], r["mode"], r["trial"]))
    with gzip.open(out, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    bad = sum(1 for r in rows if not r["correct"])
    print(f"wrote {out}: {len(rows)} solves, {bad} incorrect")


if __name__ == "__main__":
    main()
