"""Analysis: carrier masses per arm and target, heuristic predictions,
positive control, factor panel.  Writes results.json and results_table.csv."""
import os, sys, json, math, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import defaultdict
from lpoly import targets, charpoly_from_counts, divides, dim, is_weil_valid

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
T = targets()


def c_codim(G, d):
    return d * G - d * (d - 1) // 2


FCACHE = {}


def gp_factor_all(polys):
    """Factor many integer polynomials over Z with PARI/GP in one batch."""
    import subprocess, tempfile
    polys = [P for P in dict.fromkeys(polys) if P not in FCACHE]
    if not polys:
        return
    with tempfile.NamedTemporaryFile("w", suffix=".gp", delete=False) as fh:
        fh.write("default(parisizemax, 400000000);\n")
        for j, P in enumerate(polys):
            fh.write("v=%s;F=factor(Pol(v));print(\"@%d \",vector(#F~,i,[Vec(F[i,1]),F[i,2]]));\n" % (list(P), j))
        name = fh.name
    out = subprocess.run(["gp", "-q", "-f", "-s", "64M", name], input="", capture_output=True,
                         text=True, check=True).stdout
    recs = {}
    for chunk in out.split("@")[1:]:
        j, body = chunk.split(" ", 1)
        recs[int(j)] = body.replace("\n", "").replace(" ", "")
    assert len(recs) == len(polys), (len(recs), len(polys))
    for j, P in enumerate(polys):
        fl = eval(recs[j])
        FCACHE[P] = [(tuple(f), e) for f, e in fl]
    os.unlink(name)


def factor_dims(P):
    return FCACHE[P]


def load_arm(name):
    d = json.load(open(os.path.join(DATA, name + ".json")))
    hist = {tuple(int(x) for x in k.split(",")): v for k, v in d["hist"].items()}
    return d, hist


def analyse_hist(arm, G, hist, mass_per_model, label, extra_note=""):
    rows = []
    invalid = 0
    carriers = defaultdict(int)
    fac_mass = defaultdict(float)  # (dim, factor) -> mass of curves having it
    total_models = sum(hist.values())
    gp_factor_all([P for P in (charpoly_from_counts(list(N)) for N in hist) if P])
    for N, w in hist.items():
        P = charpoly_from_counts(list(N))
        if P is None or not is_weil_valid(P, G):
            invalid += w
            continue
        for nm, t in T.items():
            if dim(t) <= G and divides(t, P):
                carriers[nm] += w
        for f, e in factor_dims(P):
            fac_mass[f] += w * mass_per_model
    M = total_models * mass_per_model
    for nm, t in T.items():
        d = dim(t)
        if d > G:
            continue
        pred = M * 2.0 ** (-c_codim(G, d))
        rows.append({"arm": arm, "family": label, "genus": G, "target": nm, "dim": d,
                     "models": carriers[nm], "mass": carriers[nm] * mass_per_model,
                     "family_mass": M, "pred_mass_h1": pred,
                     "p_zero_h1": math.exp(-pred), "p_zero_h8": math.exp(-8 * pred)})
    panel = defaultdict(list)
    for f, m in fac_mass.items():
        panel[(len(f) - 1) // 2].append(m)
    return rows, invalid, {d: sorted(v, reverse=True) for d, v in panel.items()}, M


def main():
    allrows, meta, panels = [], {}, {}
    arms = [("hyp3", "hyperelliptic"), ("quartic", "plane quartic"),
            ("hyp5", "hyperelliptic"), ("trigonal5", "trigonal (nodal/cuspidal quintic)"),
            ("hyp6", "hyperelliptic"), ("quintic", "smooth plane quintic"),
            ("hyp7", "hyperelliptic"), ("hyp8", "hyperelliptic")]
    for arm, label in arms:
        p = os.path.join(DATA, arm + ".json")
        if not os.path.exists(p):
            continue
        d, hist = load_arm(arm)
        G = d["genus"]
        rows, inv, panel, M = analyse_hist(arm, G, hist, 1.0 / d["group_order"], label)
        allrows += rows
        meta[arm] = {"genus": G, "family": label, "smooth_models": d["smooth_models"],
                     "group_order": d["group_order"], "mass": M, "invalid_weil": inv,
                     "distinct_lpolys": len(hist), "seconds": d["seconds"]}
        panels[arm] = {str(k): {"n_classes": len(v), "max_mass": v[0] if v else 0,
                                "median_mass": v[len(v) // 2] if v else 0,
                                "total_mass": sum(v)} for k, v in panel.items()}
        print(arm, meta[arm], flush=True)
    # canonical genus-5 sample
    p = os.path.join(DATA, "g5_ci.json")
    if os.path.exists(p):
        d = json.load(open(p))
        nsub, nsm = d["subsample_nets"], d["subsample_smooth"]
        frac = nsm / nsub
        se = math.sqrt(frac * (1 - frac) / nsub)
        M_ci = frac * d["total_nets"] / d["gl5_order"]
        hist = {tuple(int(x) for x in k.split(",")): v for k, v in d["subsample_hist"].items()}
        # subsample: each sampled smooth net carries mass M_ci / nsm
        rows, inv, panel, _ = analyse_hist("ci5_sub", 5, hist, M_ci / nsm, "canonical CI (subsample)")
        allrows += rows
        # stage-1 filtered targets on the full sample
        n1 = d["stage1_nets"]
        nb7 = min(d["b7_filter_nets"], n1)
        hits = defaultdict(int)
        nsurv_smooth = 0
        for r in d["survivors"]:
            if not r["singular"]:
                nsurv_smooth += 1
                for h in r.get("hits", []):
                    hits[h] += 1
        for nm, t in T.items():
            dd = dim(t)
            if dd > 5 or nm.startswith(("A3", "W3")):
                continue
            nsamp = nb7 if nm.startswith("B7") else n1
            mass_per_net = d["total_nets"] / d["gl5_order"] / nsamp
            pred = M_ci * 2.0 ** (-c_codim(5, dd))
            # expected hits among sampled nets under H: pred mass / mass_per_net
            allrows.append({"arm": "ci5", "family": "canonical CI (stage-1 filter)", "genus": 5,
                            "target": nm, "dim": dd, "models": hits[nm],
                            "mass": hits[nm] * mass_per_net, "family_mass": M_ci,
                            "pred_mass_h1": pred, "p_zero_h1": math.exp(-pred),
                            "p_zero_h8": math.exp(-8 * pred),
                            "sampled_nets": nsamp,
                            "upper95_mass_if_zero": 3.0 * mass_per_net if hits[nm] == 0 else None})
        meta["ci5"] = {"genus": 5, "stage1_nets": n1, "b7_nets": nb7, "subsample_nets": nsub,
                       "subsample_smooth": nsm, "smooth_fraction": frac, "smooth_fraction_se": se,
                       "mass_est": M_ci, "mass_se": se * d["total_nets"] / d["gl5_order"],
                       "survivors": len(d["survivors"]), "survivors_smooth": nsurv_smooth,
                       "invalid_weil_sub": inv, "total_nets": d["total_nets"],
                       "gl5_order": d["gl5_order"], "seconds": d["seconds"],
                       "mass_per_net_full_sample": d["total_nets"] / d["gl5_order"] / n1}
        panels["ci5_sub"] = {str(k): {"n_classes": len(v), "max_mass": v[0] if v else 0,
                                      "median_mass": v[len(v) // 2] if v else 0,
                                      "total_mass": sum(v)} for k, v in panel.items()}
    json.dump({"meta": meta, "rows": allrows, "panels": panels}, open(os.path.join(HERE, "results.json"), "w"),
              indent=1)
    with open(os.path.join(HERE, "results_table.csv"), "w", newline="") as fh:
        keys = ["arm", "family", "genus", "target", "dim", "models", "mass", "family_mass",
                "pred_mass_h1", "p_zero_h1", "p_zero_h8"]
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in allrows:
            w.writerow(r)
    for r in allrows:
        if r["models"] or r["target"].startswith(("A5", "W5", "A7")):
            print(f'{r["arm"]:10s} g={r["genus"]} {r["target"]:9s} models={r["models"]:6d} '
                  f'mass={r["mass"]:.3f} pred={r["pred_mass_h1"]:.3g} famM={r["family_mass"]:.0f}')


if __name__ == "__main__":
    main()
