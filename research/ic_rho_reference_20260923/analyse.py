#!/usr/bin/env python3
"""Grade ledger §18.2's six targets from the frozen files, and print the tables.

Reads the calibration, evaluation and re-pricing reports beside this file and
writes `analysis.json`.  Every §18 number in the ledger note and on the
scoreboard is read from that file or from the reports themselves; nothing
here re-runs a walk.

    python3 research/ic_rho_reference_20260923/analyse.py
"""
import glob
import json
import math
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
FLOOR_S = math.sqrt(math.pi / 4)  # A = 2
BOOT = 4000


def load(path):
    with open(path) as f:
        return json.load(f)


def walk(inst, name):
    return next((w for w in inst["walks"] if w["walk"] == name), None)


def fit(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    syy = sum((y - my) ** 2 for y in ys)
    alpha = sxy / sxx
    return alpha, (sxy * sxy) / (sxx * syy) if syy else 1.0


def walks_and_abandoned(counters_list, runs):
    walks = sum(c.get("start_additions", 0) for c in counters_list) + runs
    abandoned = sum(c.get("walks_abandoned_in_a_cycle", 0) for c in counters_list)
    return walks, abandoned


def main():
    out = {"targets": {}, "tables": {}}
    calib = [load(p) for p in sorted(glob.glob(os.path.join(HERE, "calibration", "*.json")))]
    ladder = load(os.path.join(HERE, "evaluation", "ladder.json"))
    reprice = {
        os.path.basename(p)[:-5]: load(p)
        for p in sorted(glob.glob(os.path.join(HERE, "reprice", "*.json")))
    }

    # ── Target 1: every run verified; abandoned walks under 1 % ──────────
    t1 = {"runs": 0, "unverified": 0, "walks_tuned": 0, "abandoned": 0}
    for doc in calib + [ladder]:
        for inst in doc["instances"]:
            for w in inst["walks"]:
                t1["runs"] += w["runs"]
                t1["unverified"] += sum(1 for r in w["per_run"] if not r["verified"])
                if w["walk"] != "frozen-plain":
                    c = w["counters_summed"]
                    t1["walks_tuned"] += c.get("start_additions", 0) + w["runs"]
                    t1["abandoned"] += c.get("walks_abandoned_in_a_cycle", 0)
    for name, doc in reprice.items():
        for inst in doc["instances"]:
            m = inst.get("matched_reference")
            if not m:
                continue
            runs = m.get("per_run", [])
            t1["runs"] += len(runs)
            t1["unverified"] += sum(1 for r in runs if not r["verified"])
            cl = [r.get("counters", {}) for r in runs]
            w_, a_ = walks_and_abandoned(cl, len(runs))
            t1["walks_tuned"] += w_
            t1["abandoned"] += a_
    t1["abandoned_fraction"] = t1["abandoned"] / max(1, t1["walks_tuned"])
    t1["met"] = t1["unverified"] == 0 and t1["abandoned_fraction"] < 0.01
    out["targets"]["1_correct"] = t1

    # ── Target 2: the frozen walk reproduced on every recorded seed ─────
    t2 = {name: doc["identity"]["all_reproduced"] for name, doc in reprice.items()}
    t2_fail = {name: doc["identity"]["failures"] for name, doc in reprice.items() if doc["identity"]["failures"]}
    replayed = 0
    for doc in reprice.values():
        for inst in doc["instances"]:
            m = inst.get("matched_reference")
            if m:
                replayed += m.get("runs", len(m.get("per_run", [])))
    out["targets"]["2_identity"] = {"files": t2, "failures": t2_fail, "replayed_runs": replayed,
                                    "met": all(t2.values()) and len(t2) == 10}

    # ── Targets 3–5 on the evaluation ladder ────────────────────────────
    rows = []
    big = []
    for inst in ladder["instances"]:
        fz, tp, ng = walk(inst, "frozen-plain"), walk(inst, "tuned-plain"), walk(inst, "negation")
        jumps = ng["counters_summed"]["jumps"] // ng["runs"]
        pairs = list(zip(tp["per_run"], ng["per_run"]))
        assert all(a["planted"] == b["planted"] and a["seed"] == b["seed"] for a, b in pairs)
        ratio = sum(a["s_walk"] for a, _ in pairs) / sum(b["s_walk"] for _, b in pairs)
        row = {
            "regime": inst["regime"], "instance": inst["instance"], "r": inst["r"],
            "log2_r": inst["log2_r"], "jumps": jumps,
            "frozen_s": fz["mean_s"], "tuned_plain_s": tp["mean_s"], "negation_s": ng["mean_s"],
            "negation_median_s": ng["median_s"],
            "negation_walk_over_own_floor": ng["walk_over_own_floor"],
            "tuned_plain_walk_over_own_floor": tp["walk_over_own_floor"],
            "plain_over_negation_walk": ratio,
            "frozen_over_negation_s": fz["mean_s"] / ng["mean_s"],
            "negation_gae": ng["mean_gae"], "frozen_gae": fz["mean_gae"], "tuned_plain_gae": tp["mean_gae"],
            "exclusions": inst["exclusions"],
        }
        rows.append(row)
        if inst["log2_r"] >= 20.0:
            big.append((row, pairs))
    out["tables"]["ladder"] = rows

    t3 = {r["instance"]: r["negation_walk_over_own_floor"] for r, _ in big}
    out["targets"]["3_own_floor"] = {"per_size": t3, "sizes": len(t3),
                                     "met": bool(t3) and all(0.9 <= v <= 1.2 for v in t3.values())}

    # Target 4: pooled paired ratio over the sizes from 2^20, with a
    # bootstrap resampling runs within each size, both walks together.
    def pooled(samples):
        num = sum(a["s_walk"] for pairs in samples for a, _ in pairs)
        den = sum(b["s_walk"] for pairs in samples for _, b in pairs)
        return num / den
    point = pooled([p for _, p in big])
    rng = random.Random(0xB007)
    boots = []
    for _ in range(BOOT):
        boots.append(pooled([[p[rng.randrange(len(p))] for _ in range(len(p))] for _, p in big]))
    boots.sort()
    lo, hi = boots[int(0.025 * BOOT)], boots[int(0.975 * BOOT) - 1]
    out["targets"]["4_root_two"] = {"pooled_ratio": point, "ci95": [lo, hi], "sizes": len(big),
                                    "per_size": {r["instance"]: r["plain_over_negation_walk"] for r, _ in big},
                                    "met": 1.25 <= lo and hi <= 1.55}

    # Target 5: S at most 1.5 × 0.886 from 2^20; α within 0.5 ± 0.05.
    s_ok = {r["instance"]: r["negation_s"] for r, _ in big}
    fits = {}
    for label, sel in [("prime", lambda r: r["regime"] == "prime"),
                       ("char2", lambda r: r["regime"] == "char2"),
                       ("pooled", lambda r: True)]:
        pts = [(r["log2_r"], math.log2(r["negation_gae"])) for r, _ in big if sel(r)]
        if len(pts) >= 2:
            a, r2 = fit([p[0] for p in pts], [p[1] for p in pts])
            fits[label] = {"alpha": a, "r_squared": r2, "sizes": len(pts)}
    alpha_ok = all(abs(fits[k]["alpha"] - 0.5) <= 0.05 for k in fits if fits[k]["sizes"] >= 4)
    out["targets"]["5_flat"] = {"negation_s_from_2_20": s_ok, "cap": 1.5 * FLOOR_S, "fits_from_2_20": fits,
                                "met": bool(s_ok) and all(v <= 1.5 * FLOOR_S for v in s_ok.values())
                                and any(f["sizes"] >= 4 for f in fits.values()) and alpha_ok}
    # Fits over every size, every walk, for the exponent panel.
    allfits = {}
    for regime in ("prime", "char2"):
        for name, key in (("frozen-plain", "frozen_gae"), ("tuned-plain", "tuned_plain_gae"), ("negation", "negation_gae")):
            pts = [(r["log2_r"], math.log2(r[key])) for r in rows if r["regime"] == regime]
            a, r2 = fit([p[0] for p in pts], [p[1] for p in pts])
            allfits[f"{regime}/{name}"] = {"alpha": a, "r_squared": r2, "sizes": len(pts)}
    out["tables"]["fits_all_sizes"] = allfits

    # ── Target 6: no eligible automorphism beyond negation ──────────────
    bad = []
    for r in rows:
        e = r["exclusions"]
        if e.get("additional_automorphisms") or e.get("defined_over_proper_subfields_of_degree"):
            bad.append(r["instance"])
    for doc in reprice.values():
        for inst in doc["instances"]:
            e = inst.get("exclusions") or {}
            if e.get("additional_automorphisms") or e.get("defined_over_proper_subfields_of_degree"):
                bad.append(inst["instance"])
    out["targets"]["6_matched"] = {"exceptions": bad, "met": not bad}

    # ── The re-priced references ────────────────────────────────────────
    rep = []
    for name, doc in reprice.items():
        for inst in doc["instances"]:
            m = inst.get("matched_reference")
            if not m:
                continue
            rep.append({
                "file": name, "instance": inst["instance"], "regime": inst["regime"], "log2_r": inst["log2_r"],
                "frozen_s": inst["frozen_reference"]["mean_s"], "matched_s": m["mean_s"],
                "matched_method": m.get("method"),
                "frozen_over_matched": inst["frozen_over_matched"],
                "matched_walk_over_own_floor": m.get("walk_over_own_floor"),
            })
    out["tables"]["reprice_references"] = rep

    # ── The re-priced index-calculus rows ───────────────────────────────
    # §17.11: medians over curves and targets of each row's S over its own
    # instance's rho — the frozen ratio and the matched one, the same way.
    def median(xs):
        xs = sorted(xs)
        k = len(xs)
        return (xs[k // 2] if k % 2 else (xs[k // 2 - 1] + xs[k // 2]) / 2) if k else float("nan")
    whole = {}
    for name, doc in reprice.items():
        if not name.startswith("W"):
            continue
        n = name.split("-")[0]
        inst = doc["instances"][0]
        for row in inst["rows"]:
            label = row["label"]
            engine = label.split(" + ")[2] if " + " in label else label
            engine = "pair table" if engine == "—" else engine
            cell = whole.setdefault(n, {}).setdefault(engine, {"s": [], "frozen": [], "matched": []})
            cell["s"].append(row["s"])
            cell["frozen"].append(row["vs_rho_frozen"])
            cell["matched"].append(row["vs_rho_matched"])
    whole_table = {
        n: {e: {"rows": len(c["s"]), "median_s": median(c["s"]),
                "median_vs_rho_frozen": median(c["frozen"]), "median_vs_rho_matched": median(c["matched"])}
            for e, c in cells.items()}
        for n, cells in sorted(whole.items())
    }
    out["tables"]["whole_method_17_11"] = whole_table

    # The boundary ladder: every variant, mean S over its targets, against
    # the frozen and the matched reference; the best row per instance.
    boundary = {}
    for name, doc in reprice.items():
        if not name.startswith("ic-boundary"):
            continue
        for inst in doc["instances"]:
            if "variants" not in inst:
                continue
            best = min(inst["variants"], key=lambda v: v["mean_s"])
            boundary.setdefault(name, []).append({
                "instance": inst["instance"], "regime": inst["regime"], "log2_r": inst["log2_r"],
                "rho_frozen": inst["frozen_reference"]["mean_s"], "rho_matched": inst["matched_reference"]["mean_s"],
                "best_variant": best["variant"], "best_s": best["mean_s"],
                "best_vs_rho_frozen": best["vs_rho_frozen"], "best_vs_rho_matched": best["vs_rho_matched"],
                "variants": inst["variants"],
            })
    out["tables"]["boundary"] = boundary

    # The rho reference's exponent over the Round-5 ladder, frozen and
    # matched, the way the ledger fits it: log2 of the mean total against
    # log2 r, per regime.
    rho_fits = {}
    for regime in ("prime", "char2"):
        pts = [x for x in boundary.get("ic-boundary-ledger-round5-2026-09-22", []) if x["regime"] == regime]
        for which in ("rho_frozen", "rho_matched"):
            xs = [p["log2_r"] for p in pts]
            ys = [math.log2(p[which] * 2 ** (p["log2_r"] / 2)) for p in pts]
            a, r2 = fit(xs, ys)
            rho_fits[f"{regime}/{which}"] = {"alpha": a, "r_squared": r2, "sizes": len(pts)}
    out["tables"]["rho_exponent_round5_ladder"] = rho_fits

    with open(os.path.join(HERE, "analysis.json"), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
        f.write("\n")

    # ── Print ────────────────────────────────────────────────────────────
    for k, v in out["targets"].items():
        print(k, "MET" if v["met"] else "NOT MET", json.dumps({x: y for x, y in v.items() if x != "met"})[:400])
    print()
    print("| regime | curve | log₂ r | J | frozen walk `S` | tuned, points `S` | **negation `S`** (median) | negation walk / own floor | points / negation, walk ops | frozen / negation, `S` |")
    print("|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|")
    for r in rows:
        print(f"| {r['regime']} | `{r['instance']}` | {r['log2_r']:.1f} | {r['jumps']} | {r['frozen_s']:.2f} | "
              f"{r['tuned_plain_s']:.2f} | **{r['negation_s']:.2f}** ({r['negation_median_s']:.2f}) | "
              f"{r['negation_walk_over_own_floor']:.3f} | {r['plain_over_negation_walk']:.3f} | {r['frozen_over_negation_s']:.2f}× |")
    print()
    for k, v in allfits.items():
        print(f"fit {k}: alpha {v['alpha']:.3f} R2 {v['r_squared']:.3f} over {v['sizes']} sizes")
    for k, v in out["targets"]["5_flat"]["fits_from_2_20"].items():
        print(f"fit from 2^20 {k}: alpha {v['alpha']:.3f} R2 {v['r_squared']:.3f} over {v['sizes']} sizes")
    print()
    print("| report | curve | log₂ r | frozen rho `S` | matched rho `S` | frozen / matched | matched walk / own floor |")
    print("|:--|:--|--:|--:|--:|--:|--:|")
    for r in rep:
        own = r["matched_walk_over_own_floor"]
        print(f"| `{r['file']}` | `{r['instance']}` | {r['log2_r']:.1f} | {r['frozen_s']:.2f} | {r['matched_s']:.2f} | "
              f"{r['frozen_over_matched']:.2f}× | {own:.3f} |" if own is not None else "")
    print()
    print("§17.11 re-priced: median S · median vs rho frozen → matched")
    for n, cells in whole_table.items():
        for e, c in cells.items():
            print(f"  {n} {e:>14}: rows {c['rows']} S {c['median_s']:.1f} vs rho {c['median_vs_rho_frozen']:.1f}× → {c['median_vs_rho_matched']:.1f}×")
    print()
    for name, insts in boundary.items():
        for x in insts:
            print(f"{name} {x['regime']} {x['instance']} 2^{x['log2_r']:.1f}: rho {x['rho_frozen']:.2f} → {x['rho_matched']:.2f}; "
                  f"best {x['best_variant']} S {x['best_s']:.2f}: {x['best_vs_rho_frozen']:.2f}× → {x['best_vs_rho_matched']:.2f}×")
    print()
    for k, v in rho_fits.items():
        print(f"rho exponent over the Round-5 ladder {k}: {v['alpha']:.3f} (R² {v['r_squared']:.3f}, {v['sizes']} sizes)")


if __name__ == "__main__":
    main()
