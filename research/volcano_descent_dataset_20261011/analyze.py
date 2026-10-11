#!/usr/bin/env python3
"""Summaries of the dataset, the ECDLP experiments and a leakage-aware ML baseline.

Writes tables to stdout as Markdown; run with `--write RESULTS_TABLES.md` to save.
Bootstrap intervals are percentile intervals over trials (resampling trials
within each group, 2000 resamples, fixed seed).
"""
from __future__ import annotations

import argparse
import io
import json
import math
import os
import sys

import numpy as np
import pandas as pd

RNG = np.random.default_rng(20261011)


def md(df, floatfmt=3):
    return df.to_markdown(index=False, floatfmt=f".{floatfmt}g") if len(df) else "(empty)"


def boot_mean_ci(x, n=2000):
    x = np.asarray(x, float)
    if len(x) < 2:
        return float(x.mean()), float("nan"), float("nan")
    idx = RNG.integers(0, len(x), (n, len(x)))
    m = x[idx].mean(axis=1)
    return float(x.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def boot_ratio_ci(a, b, n=2000):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ia = RNG.integers(0, len(a), (n, len(a)))
    ib = RNG.integers(0, len(b), (n, len(b)))
    r = a[ia].mean(axis=1) / b[ib].mean(axis=1)
    return float(a.mean() / b.mean()), float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5))


BIG = ["j", "a", "b", "comp", "up_j", "n1", "n2"]


def read_vertices(data):
    """Large integers (up to ~70 bits) are kept exact as Python ints in object columns."""
    ver = pd.read_csv(os.path.join(data, "vertices.csv.gz"), dtype={c: str for c in BIG})
    for c in ("j", "a", "b", "n1", "n2"):
        ver[c] = ver[c].map(int)
    return ver


def campaign_tables(data):
    vol = pd.read_csv(os.path.join(data, "volcanoes.csv.gz"), dtype={"p": str, "N": str, "t": str, "f": str, "largest_prime": str})
    ver = read_vertices(data)
    des = pd.read_csv(os.path.join(data, "descents.csv.gz"))
    out = []
    out.append("### Campaign\n")
    g = vol.groupby("ell").agg(volcanoes=("vol_id", "count"), max_depth=("k", "max"),
                               vertices=("n_vertices", "sum"), stored=("n_stored", "sum"),
                               max_p_bits=("p_bits", "max")).reset_index()
    out.append(md(g))
    out.append(f"\nvolcanoes: {len(vol)}; vertices in the exact volcanoes: {int(vol.n_vertices.sum())}; "
               f"rows stored: {len(ver)}; deepest volcano: {int(vol.k.max())} levels.\n")
    out.append("### Verification (fail-closed)\n")
    chk = pd.json_normalize(vol.checks.apply(json.loads))
    rows = [(c, int(chk[c].sum()), len(chk)) for c in chk.columns]
    out.append(md(pd.DataFrame(rows, columns=["check", "volcanoes passing", "of"])))
    out.append(f"\nvertex rows whose measured group structure equals the O_d/(pi-1) prediction: "
               f"{int(ver.group_matches_theory.sum())} of {len(ver)}.")
    out.append(f"\nblind descents (start vertex only j, p, t, ell): {int(des.success.sum())} of {len(des)} reach a floor curve "
               f"with the correct level and path length; median {des.seconds.median():.4f} s, max {des.seconds.max():.3f} s.\n")
    lv = ver.groupby(["vol_id"]).agg(n=("j", "count"), n1max=("n1", "max")).reset_index()
    out.append(f"volcanoes with a non-cyclic group somewhere (level-dependent group structure): "
               f"{int((lv.n1max > 1).sum())} of {len(lv)}.\n")
    return "\n".join(out), vol, ver


def ecdlp_tables(data):
    out = []
    for name in ("E1", "E2"):
        p = os.path.join(data, f"ecdlp_{name}.csv.gz")
        if not os.path.exists(p):
            continue
        d = pd.read_csv(p, dtype={"p": str, "j": str, "ord_G": str})
        d["ord_G"] = d["ord_G"].map(int)
        out.append(f"### ECDLP {name}\n")
        out.append(f"{len(d)} solves, {int((d.correct == 0).sum())} incorrect, "
                   f"{d.vol_id.nunique()} volcanoes, {d.groupby(['vol_id','j']).ngroups} vertices, "
                   f"largest-prime sizes {d.r_max_bits.min()}-{d.r_max_bits.max()} bits.\n")
        d["norm_rho"] = d.ops_rho / np.sqrt(math.pi / 2 * 2.0 ** d.r_max_bits)   # normalised by sqrt(pi r/2), r ~ 2^bits
        d["layer"] = np.where(d.level == 0, "crater", np.where(d.level == d.k, "floor", "middle"))
        rows = []
        for (mode, layer), g in d.groupby(["mode", "layer"]):
            m, lo, hi = boot_mean_ci(g.norm_rho)
            rows.append((mode, layer, len(g), g.groupby(["vol_id", "j"]).ngroups, m, lo, hi))
        out.append(md(pd.DataFrame(rows, columns=["mode", "layer", "solves", "vertices", "mean ops_rho / sqrt(pi r/2)", "ci_lo", "ci_hi"])))
        out.append("")
        # speedup of neg over plain, aut over neg: layer-wise ratio of means
        rows = []
        for layer, g in d.groupby("layer"):
            def A(mode):
                return g[g["mode"] == mode].norm_rho.values
            for num, den in (("plain", "neg"), ("neg", "aut")):
                if len(A(num)) and len(A(den)):
                    r, lo, hi = boot_ratio_ci(A(num), A(den))
                    rows.append((layer, f"{num} / {den}", r, lo, hi))
        out.append(md(pd.DataFrame(rows, columns=["layer", "ratio of mean ops", "ratio", "ci_lo", "ci_hi"])))
        out.append("")
        if name == "E2":
            out.append(e2_tables(d))
            out.append("")
        if name == "E1":
            # dependence on level at fixed solver: slope of normalised cost vs level, and the group effect
            rows = []
            for (vid, mode), g in d.groupby(["vol_id", "mode"]):
                lv = g.groupby("level").norm_rho.mean()
                x, y = lv.index.values.astype(float), lv.values
                slope = np.polyfit(x, y, 1)[0] if len(x) > 2 else float("nan")
                rows.append((vid, mode, int(g.m.iloc[0]), slope, float(y.mean())))
            t = pd.DataFrame(rows, columns=["volcano", "mode", "m", "slope of norm. cost per level", "mean norm. cost"])
            out.append(md(t))
            out.append("")
            # naive cost without Pohlig-Hellman: sqrt(ord(G)) vs measured total
            d["log2_ord"] = np.log2(d.ord_G.astype(float))
            rows = []
            for mm, g in d[d["mode"] == "plain"].groupby("m"):
                rows.append((int(mm), g.ops_total.mean(), g.ops_rho.mean(), (2.0 ** (g.log2_ord / 2)).mean()))
            out.append(md(pd.DataFrame(rows, columns=["m", "mean ops_total (PH-aware)", "mean ops_rho", "mean sqrt(ord G) (no PH)"])))
            out.append("")
    return "\n".join(out)


def e2_tables(d):
    """Automorphism experiment: per volcano and layer, plus pooled by D_K."""
    out = ["Per volcano and layer (mean of ops_rho / sqrt(pi r / 2); ratios of means, bootstrap 95% intervals):\n"]
    rows = []
    for vid, g in d.groupby("vol_id"):
        for layer in ("crater", "middle", "floor"):
            gl = g[g.layer == layer]
            A = lambda mode: gl[gl["mode"] == mode].norm_rho.values
            row = dict(volcano=vid, layer=layer, solves=len(gl), vertices=gl.j.nunique(),
                       plain=A("plain").mean(), neg=A("neg").mean(), aut=A("aut").mean())
            for num, den, key in (("plain", "neg", "plain/neg"), ("neg", "aut", "neg/aut"), ("plain", "aut", "plain/aut")):
                r, lo, hi = boot_ratio_ci(A(num), A(den))
                row[key] = f"{r:.2f} [{lo:.2f}, {hi:.2f}]"
            rows.append(row)
    out.append(pd.DataFrame(rows).to_markdown(index=False, floatfmt=".3f"))
    out.append("\nPooled by D_K (ideal gains: j = 0, mu_6: neg/aut 1.73, plain/aut 2.45; j = 1728, mu_4: neg/aut 1.41, plain/aut 2.00; "
               "elsewhere neg/aut 1.00):\n")
    rows = []
    for D, g in d.groupby("D_K"):
        for layer in ("crater", "not crater"):
            gl = g[g.layer == "crater"] if layer == "crater" else g[g.layer != "crater"]
            A = lambda mode: gl[gl["mode"] == mode].norm_rho.values
            r1, r2 = boot_ratio_ci(A("neg"), A("aut")), boot_ratio_ci(A("plain"), A("aut"))
            rows.append((D, layer, len(gl) // 3, f"{r1[0]:.2f} [{r1[1]:.2f}, {r1[2]:.2f}]", f"{r2[0]:.2f} [{r2[1]:.2f}, {r2[2]:.2f}]"))
    out.append(pd.DataFrame(rows, columns=["D_K", "layer", "solves per mode", "neg / aut", "plain / aut"]).to_markdown(index=False))
    return "\n".join(out)


def ml_baseline(ver, vol):
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import balanced_accuracy_score
    from sklearn.model_selection import GroupShuffleSplit
    out = ["### ML baseline: predicting the level (distance to the crater) from observables\n"]
    v = ver.merge(vol[["vol_id", "p", "k", "ell", "D_K", "m", "c"]], on="vol_id", suffixes=("", "_vol"))
    feats = pd.DataFrame({
        "j_over_p": [int(x) / int(q) for x, q in zip(v["j"], v["p"])],
        "a_over_p": [int(x) / int(q) for x, q in zip(v["a"], v["p"])],
        "b_over_p": [int(x) / int(q) for x, q in zip(v["b"], v["p"])],
        "j_popcount": v["j"].map(lambda x: bin(int(x)).count("1")),
        "j_mod_ell": [int(x) % int(e) for x, e in zip(v["j"], v["ell"])],
        "a_mod_ell": [int(x) % int(e) for x, e in zip(v["a"], v["ell"])],
        "b_mod_ell": [int(x) % int(e) for x, e in zip(v["b"], v["ell"])],
    })
    sig = [c for c in v.columns if c.startswith("sig_")]
    group_f = pd.DataFrame({"n1": v["n1"], "log2_n1": np.log2(v["n1"].astype(float))})
    sigx = v[sig].copy()
    # the signature entry for the volcano's own ell equals ell+1 off the floor and 1 on it: graph-derived
    sig_noown = sigx.copy()
    for q in (2, 3, 5, 7, 11, 13):
        sig_noown.loc[v["ell"] == q, f"sig_{q}"] = -1
    target = v["height"].astype(int)           # distance to the floor
    sets = {
        "raw curve only (j,a,b stats)": feats,
        "+ degrees for other ell (no own-ell degree)": pd.concat([feats, sig_noown], axis=1),
        "+ group structure (n1)": pd.concat([feats, sig_noown, group_f], axis=1),
        "+ own-ell degree (graph-derived)": pd.concat([feats, sigx, group_f], axis=1),
    }
    rows = []
    n1max = v.groupby("vol_id")["n1"].transform("max")      # observed, not the requested m
    for subset_name, mask in (("every curve in the volcano has a cyclic group", n1max == 1),
                              ("group structure varies with level", n1max > 1)):
        sub = v[mask].index
        if len(sub) < 200:
            continue
        groups = v.loc[sub, "vol_id"].values
        gss = GroupShuffleSplit(n_splits=1, test_size=0.3, random_state=0)
        tr, te = next(gss.split(sub, groups=groups))
        ytr, yte = target.loc[sub].iloc[tr], target.loc[sub].iloc[te]
        # per-volcano normalised target: height / k  (so volcanoes of different depth are comparable)
        kk = v.loc[sub, "k"].astype(float)
        ytr_n = np.round(10 * ytr.values / kk.iloc[tr].values).astype(int)
        yte_n = np.round(10 * yte.values / kk.iloc[te].values).astype(int)
        maj = pd.Series(ytr_n).mode()[0]
        base = balanced_accuracy_score(yte_n, np.full_like(yte_n, maj))
        for name, X in sets.items():
            clf = HistGradientBoostingClassifier(max_iter=150, random_state=0)
            clf.fit(X.loc[sub].iloc[tr], ytr_n)
            pred = clf.predict(X.loc[sub].iloc[te])
            rows.append((subset_name, name, len(tr), len(te), balanced_accuracy_score(yte_n, pred), base))
    out.append(md(pd.DataFrame(rows, columns=["volcanoes", "feature set", "train rows", "test rows", "balanced accuracy", "constant-prediction baseline"])))
    out.append("\nTarget: round(10 * height / k), 11 classes (0 = crater side ... 10 = at the floor is height 0), held-out volcanoes (grouped split), so no volcano contributes to both sides.\n")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--write", default=None)
    ap.add_argument("--no-ml", action="store_true")
    args = ap.parse_args()
    parts = []
    camp, vol, ver = campaign_tables(args.data)
    parts.append(camp)
    parts.append(ecdlp_tables(args.data))
    if not args.no_ml:
        parts.append(ml_baseline(ver, vol))
    text = "\n".join(parts)
    print(text)
    if args.write:
        open(args.write, "w").write(text)


if __name__ == "__main__":
    main()
