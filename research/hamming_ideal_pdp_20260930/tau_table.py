#!/usr/bin/env python3
"""Budget-sensitivity table from results/tau/*.jsonl: the same n = 11 targets
under per-call budgets 2^26, 2^28 and 2^30 row-word XORs.

    python3 tau_table.py
"""
import glob, json, os, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))

def main():
    cells = {}
    for f in sorted(glob.glob(os.path.join(HERE, "results", "tau", "n*_tau*.jsonl"))):
        manifest, rows = None, []
        with open(f) as fh:
            for line in fh:
                r = json.loads(line)
                if r["kind"] == "manifest":
                    manifest = r
                else:
                    rows.append(r)
        cells[(rows[0]["encoding"], manifest["budget_xor"])] = (manifest, rows)
    print("| encoding | budget (row-word XORs) | targets | correct | calls / target | tame depth | budget-wild calls / target | XOR words / target | wall / target (s) |")
    print("|:--|--:|--:|--:|--:|--:|--:|--:|--:|")
    for (enc, b), (m, rows) in sorted(cells.items()):
        depths = [d for r in rows for d in r["tame_depths"]]
        print(f"| {enc} | 2^{b.bit_length() - 1} | {len(rows)} | {sum(r['agree'] == 'ok' for r in rows)}/{len(rows)} | "
              f"{st.fmean(r['calls'] for r in rows):.2f} | {st.fmean(depths):.2f} | {st.fmean(r['budget_calls'] for r in rows):.2f} | "
              f"{st.fmean(r['xor_words'] for r in rows):.3g} | {st.fmean(r['solve_us'] for r in rows) / 1e6:.1f} |")
    # per-target call counts, to show whether the tree changes with the budget
    print()
    print("| encoding | seed | calls at 2^26 | calls at 2^28 | calls at 2^30 | tame depths at 2^26 / 2^28 / 2^30 |")
    print("|:--|--:|--:|--:|--:|:--|")
    for enc in ["FC", "QFC", "MONO"]:
        by_b = {b: {r["seed"]: r for r in rows} for (e, b), (m, rows) in cells.items() if e == enc}
        bs = sorted(by_b)
        for seed in sorted(by_b[bs[0]]):
            rs = [by_b[b][seed] for b in bs]
            print(f"| {enc} | {seed} | " + " | ".join(str(r["calls"]) for r in rs) + " | " + " / ".join(",".join(map(str, r["tame_depths"])) for r in rs) + " |")

if __name__ == "__main__":
    main()
