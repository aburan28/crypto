#!/usr/bin/env python3
"""Render the Markdown tables of RESULT.md and the research note from
results/<tag>/summary.json (which analyse.py derives from the frozen cells).

    python3 render_tables.py --tag main systems|results|fits
"""
import argparse, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
ENC_ORDER = ["SUB", "MONO", "C", "FC", "QFC"]

def fmt(v, digits=2):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "—"
    if isinstance(v, float):
        if abs(v) >= 1e5:
            return f"{v:.2e}"
        if abs(v) >= 100:
            return f"{v:.0f}"
        return f"{v:.{digits}f}"
    return str(v)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="main")
    ap.add_argument("what", choices=["systems", "results", "fits"])
    a = ap.parse_args()
    s = json.load(open(os.path.join(HERE, "results", a.tag, "summary.json")))
    rows = sorted(s["rows"], key=lambda r: (r["n"], ENC_ORDER.index(r["encoding"])))
    if a.what == "systems":
        print("| n | encoding | base | \\|F\\| (points) | variables | of which auxiliary | generators | max degree |")
        print("|--:|:--|:--|--:|--:|--:|--:|--:|")
        for r in rows:
            print(f"| {r['n']} | {r['encoding']} | {r['base']} | {r['F_points']} | {r['vars']} | {r['aux_vars']} | {r['gens']} | {r['max_gen_degree']} |")
    elif a.what == "results":
        print("| n | encoding | targets | decomposable | correct | timeouts | calls / target | calls, no-targets | tame | wild | budget-wild | tame depth | XOR words / call | wall / target (s) | \\|F\\| | calls / \\|F\\| |")
        print("|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
        for r in rows:
            print(f"| {r['n']} | {r['encoding']} | {r['targets']} | {r['decomposable']} | {r['agree']}/{r['targets']} | {r['timeouts']} | {fmt(r['calls_mean'])} | {fmt(r['calls_no'])} | {fmt(r['tame_mean'])} | {fmt(r['wild_mean'])} | {fmt(r['budget_mean'])} | {fmt(r['tame_depth_mean'])} | {fmt(r['xor_per_call'])} | {fmt(r['solve_ms_mean'] / 1000)} | {r['F_points']} | {fmt(r['calls_over_F'], 3)} |")
    else:
        print("| encoding | sizes fitted | exponent of calls / target in \\|F\\| | exponent of XOR words / target in \\|F\\| |")
        print("|:--|:--|--:|--:|")
        xf = {f["encoding"]: f for f in s["fits_xor"]}
        for f in s["fits_calls"]:
            x = xf.get(f["encoding"])
            print(f"| {f['encoding']} | {', '.join(map(str, f['sizes']))} | {f['exponent']:.2f} | {x['exponent']:.2f} |" if x else f"| {f['encoding']} | {', '.join(map(str, f['sizes']))} | {f['exponent']:.2f} | — |")

if __name__ == "__main__":
    main()
