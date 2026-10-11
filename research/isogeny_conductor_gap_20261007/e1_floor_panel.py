#!/usr/bin/env python3
"""E1 panel driver: this repo's production oracles on the T23 crater curve and floor curves.

Runs `ic run --json` for each (curve, solver, seed) cell with the same degree, factor base
(legacy single-factor subspace family, dimension ord_23(2) = 11) and summand count, and
records the JSON report.  Floor curves are passed by their coefficient b (bitmask in the
modulus x^23 + x^5 + 1, which is both the T23 family's modulus and the one the Rust sparse
search selects at n = 23).  The crater curve runs through the ordinary Koblitz path.

  python3 -I e1_floor_panel.py --ic ../../target/release/ic --out e1_t23_panel.jsonl \
      --solvers enumerate groebner sat --seeds 8 --max-trials 20000

Cells never share state; failures and timeouts are kept as rows.
"""
import argparse, json, os, subprocess, sys, time

T23 = "/Volumes/SSD990/ecdlp-hardness-work/toy-analogue-a/data/T23_family.json"
DEGREE, CURVE_A = 23, 0

def pick_curves(path, per_orbit=1, orbits=2):
    d = json.load(open(path))
    cs = d["curves"] if isinstance(d["curves"], list) else list(d["curves"].values())
    out = [("E0", None, "crater")]
    seen = {}
    for c in cs:
        if c["orbit"] == "crater":
            continue
        seen.setdefault(c["orbit"], []).append(c)
    for name in sorted(seen)[:orbits]:
        for c in seen[name][:per_orbit]:
            out.append((c["label"], int(c["b_int"]), name))
    return out, d.get("meta", {})

def run_cell(ic, label, b, solver, seed, max_trials, timeout):
    cmd = [ic, "run", "--json", "--degree", str(DEGREE), "--curve-a", str(CURVE_A),
           "--solver", solver, "--summands", "2", "--random-target", "--seed", str(seed),
           "--max-trials", str(max_trials)]
    if b is not None:
        cmd += ["--model-b", str(b)]
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        wall = time.time() - t0
        rep = None
        try:
            rep = json.loads(p.stdout)
        except Exception:
            pass
        row = {"curve": label, "b": b, "solver": solver, "seed": seed, "wall_s": round(wall, 3),
               "exit": p.returncode, "status": (rep or {}).get("status"), "report": rep,
               "stderr_tail": p.stderr[-400:] if p.returncode != 0 else ""}
    except subprocess.TimeoutExpired:
        row = {"curve": label, "b": b, "solver": solver, "seed": seed, "wall_s": timeout,
               "exit": None, "status": "timeout", "report": None, "stderr_tail": ""}
    return row

def summarize(rows):
    def stage(rep, name, state):
        for s in (rep or {}).get("stages", []):
            if s.get("stage") == name and s.get("status") == state:
                return s
        return None
    table = {}
    for r in rows:
        key = (r["curve"], r["solver"])
        rep = r["report"] or {}
        rc = stage(rep, "relation_collection", "stopped") or {}
        det = rc.get("details", {})
        fb = (stage(rep, "factor_base", "complete") or {}).get("details", {})
        t = table.setdefault(key, {"cells": 0, "complete": 0, "trials": [], "relations": [],
                                   "rc_wall": [], "total_wall": [], "columns": None, "points": None,
                                   "frobenius": None})
        t["cells"] += 1
        t["complete"] += int(r["status"] == "complete")
        if det:
            t["trials"].append(det.get("trials")); t["relations"].append(det.get("collected"))
            t["rc_wall"].append(rc.get("elapsed_seconds"))
        t["total_wall"].append(r["wall_s"])
        if fb:
            t["columns"], t["points"] = fb.get("columns"), fb.get("points")
        m = rep.get("model") or {}
        if m:
            t["frobenius"] = m.get("frobenius_endomorphism")
    def med(xs):
        xs = sorted(x for x in xs if isinstance(x, (int, float)))
        return None if not xs else xs[len(xs) // 2]
    lines = ["curve | solver | cells | complete | points | columns | frobenius | median trials | median relations | median trials/relation | median rc wall s | median total wall s",
             "---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:"]
    for (curve, solver), t in sorted(table.items()):
        mt, mr = med(t["trials"]), med(t["relations"])
        tpr = round(mt / mr, 2) if (mt and mr) else None
        lines.append(f"{curve} | {solver} | {t['cells']} | {t['complete']} | {t['points']} | {t['columns']} | {t['frobenius']} | {mt} | {mr} | {tpr} | {med(t['rc_wall'])} | {med(t['total_wall'])}")
    return "\n".join(lines)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ic", default=os.path.join(os.path.dirname(__file__), "..", "..", "target", "release", "ic"))
    ap.add_argument("--out", default="e1_t23_panel.jsonl")
    ap.add_argument("--solvers", nargs="+", default=["enumerate", "groebner", "sat"])
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--max-trials", type=int, default=20000)
    ap.add_argument("--timeout", type=int, default=1800)
    ap.add_argument("--orbits", type=int, default=2)
    ap.add_argument("--only-curves", nargs="*", default=None, help="restrict to these labels (E0, O00-000, ...)")
    a = ap.parse_args()
    curves, meta = pick_curves(T23, orbits=a.orbits)
    if a.only_curves:
        curves = [c for c in curves if c[0] in set(a.only_curves)]
    print("curves:", curves, "\nT23 meta modulus:", meta.get("modulus"), file=sys.stderr)
    rows = []
    with open(a.out, "a") as f:
        for label, b, orbit in curves:
            for solver in a.solvers:
                for k in range(a.seeds):
                    seed = 0x5EED_0000 + k
                    row = run_cell(a.ic, label, b, solver, seed, a.max_trials, a.timeout)
                    row["orbit"] = orbit
                    rows.append(row)
                    f.write(json.dumps(row) + "\n"); f.flush()
                    print(f"{label:8s} {solver:9s} seed={k} status={row['status']} wall={row['wall_s']}s", file=sys.stderr)
    print(summarize(rows))

if __name__ == "__main__":
    main()
