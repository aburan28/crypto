"""Score the registered cells (c23, c37, c41) and the amended follow-ups (c37b, c41b) against PROTOCOL.md M1-M7.
Per-cell verdicts are in results/score.json; the aggregate flags cover every present cell."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CELLS = {"C23": "c23.json", "C37": "c37.json", "C41": "c41.json", "C37b": "c37b.json", "C41b": "c41b.json"}


def load(name):
    with open(os.path.join(HERE, "results", name)) as fh:
        return json.load(fh)


def main():
    out = {"cells": {}, "M1": True, "M2": True, "M3": True, "M3_decisive": True,
           "M4": True, "M5": True, "M6": True, "M7": True}
    for cell, fname in CELLS.items():
        if not os.path.exists(os.path.join(HERE, "results", fname)):
            out["cells"][cell] = "missing"
            continue
        res = load(fname)
        n, r, kc = res["n"], res["r"], res["Kc"]
        arms = res["arms"]
        frob = next(a for a in arms if a["kind"] == "frob")
        L = frob["L"]
        model = 2.0 * r / (L * L)
        rows = []
        m1 = 0.5 * model <= frob["mean_first_hit_probes"] <= 2.0 * model
        out["M1"] &= m1
        by_k = {}
        for a in arms:
            k = a["K"]
            ratio = a["mean_first_hit_probes"] / frob["mean_first_hit_probes"]
            canon_ratio = a["canonical"] / frob["canonical"]
            row = dict(kind=a["kind"], K=k, c=a["c"], L=a["L"], pairs=a["pairs"],
                       canonical=a["canonical"], canonical_ratio=round(canon_ratio, 3),
                       mean_probes=round(a["mean_first_hit_probes"], 2),
                       median_probes=a["median_first_hit_probes"],
                       probes_vs_frob=round(ratio, 3),
                       rows_to_full_rank=a["rows_to_full_rank"],
                       degenerate=a["degenerate_rows"],
                       total_probes_to_full_rank=a["total_probes_to_full_rank"],
                       logs_ok=a["logs_ok"], target_ok=a["target_ok"])
            if a["kind"] != "frob" and k >= 3:
                row["M2"] = 0.7 <= ratio <= 1.5
                out["M2"] &= row["M2"]
                exp = kc / k
                tr = frob["total_probes_to_full_rank"] / max(1, a["total_probes_to_full_rank"])
                row["M7_ratio_over_Kc_K"] = round(tr / exp, 3)
                row["M7"] = 0.4 <= tr / exp <= 2.5
                out["M7"] &= row["M7"]
            if a["kind"] != "frob" and k == 1:
                row["M4"] = ratio >= 1.5 and canon_ratio <= 0.75
                out["M4"] &= row["M4"]
            row["M5"] = (a["rows_to_full_rank"] is not None and a["rows_to_full_rank"] <= k + 2
                         and a["degenerate_rows"] == 0)
            out["M5"] &= row["M5"]
            row["M6"] = bool(a["logs_ok"] and a["target_ok"])
            out["M6"] &= row["M6"]
            by_k.setdefault(k, {})[a["kind"]] = a["mean_first_hit_probes"]
            rows.append(row)
        m3 = {}
        for k, d in sorted(by_k.items()):
            if "tau" in d and "int" in d:
                q = d["tau"] / d["int"]
                m3[k] = round(q, 3)
                out["M3"] &= 0.6 <= q <= 1.6
                if k >= 2:
                    out["M3_decisive"] &= q <= 0.67
        out["cells"][cell] = dict(n=n, r=r, Kc=kc, L=L, model_2r_over_L2=round(model, 2),
                                  M1=m1, tau_over_int_by_K=m3, arms=rows)
    out["all_pass"] = all(out[m] for m in ("M1", "M2", "M3", "M4", "M5", "M6", "M7"))
    with open(os.path.join(HERE, "results", "score.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "cells"}, indent=1))
    for cell, c in out["cells"].items():
        if c == "missing":
            print(f"\n{cell}: missing")
            continue
        print(f"\n{cell}: n={c['n']} r={c['r']} Kc={c['Kc']} L={c['L']} 2r/L^2={c['model_2r_over_L2']} tau/int by K: {c['tau_over_int_by_K']}")
        print("kind  K   c  canonical  canon/frob  mean_probes  vs_frob  rows_to_rank  total_probes  logs target")
        for a in c["arms"]:
            print(f"{a['kind']:5}{a['K']:3}{a['c']:4}{a['canonical']:10}{a['canonical_ratio']:12.3f}{a['mean_probes']:13.2f}{a['probes_vs_frob']:9.3f}{str(a['rows_to_full_rank']):>14}{a['total_probes_to_full_rank']:14}  {a['logs_ok']} {a['target_ok']}")


if __name__ == "__main__":
    main()
