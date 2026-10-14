"""Print Markdown tables of recorded values from ../data (layout only)."""
import json, pathlib
DATA = pathlib.Path(__file__).resolve().parent.parent / "data"

print("| N | n | l | M | curves | classes | direct | class | 3-walk | class÷direct | weak classes | median log2 sqrt(class) | max walk steps |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for N in [8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24]:
    p = DATA / f"ghs_N{N}.json"
    if not p.exists():
        continue
    rep = json.loads(p.read_text())
    for f in rep["factorisations"]:
        if f["l"] < 2 or f["n"] < 2:
            continue
        for t in f["thresholds"]:
            if t["max_m"] > min(4, f["n"] - 1) and not (f["n"] == 2 and t["max_m"] == 2):
                continue
            T = t["curves_total"]
            d, c, w = t["direct_weak_curves"], t["class_covered_curves"], t["walk3_covered_curves"]
            if d == 0:
                continue
            hist = t["walk3_distance_histogram"]
            print(f"| {N} | {f['n']} | {f['l']} | {t['max_m']} | {T} | {t['classes_total']} | {d/T:.2e} | {c/T:.2e} | "
                  + (f"{w/T:.2e}" if w is not None else "n/a") + f" | {c/d:.1f} | {t['classes_with_weak_member']} | "
                  + (f"{t['median_log2_sqrt_weak_class_size']:.2f}" if t['median_log2_sqrt_weak_class_size'] is not None else "") + " | "
                  + (str(len(hist) - 1) if hist else "n/a") + " |")

print()
print("| p | curves | classes | direct tractable | class closure | any Lagrangian | χ̄ = g·g* | charpoly mismatches | violations |")
print("|---|---|---|---|---|---|---|---|---|")
for p in (3, 5, 7, 11, 13):
    f = DATA / f"smith_p{p}.json"
    if not f.exists():
        continue
    r = json.loads(f.read_text())
    C = r["curves"]
    print(f"| {p} | {C} | {r['classes']} | {r['direct_pair_kernel_curves']/C:.4f} | {r['class_pair_kernel_curves']/C:.4f} | "
          f"{r['direct_any_lagrangian_curves']/C:.4f} | {r['necessary_condition_curves']/C:.4f} | {r['charpoly_mismatches']} | {r['necessary_condition_violations']} |")
c = json.loads((DATA / "smith_chebotarev.json").read_text())
print()
print("| cycle type | #S_8 | #S(H) | stable Lagrangians | χ̄ mod 2 (coeff. bits, palindromic) |")
print("|---|---|---|---|---|")
for r in sorted(c["rows"], key=lambda r: (r["charpoly_mod2"], r["cycle_type"])):
    print(f"| {'+'.join(map(str, r['cycle_type']))} | {r['class_size']} | {r['stable_pair_kernels']} | {r['stable_lagrangians']} | {r['charpoly_mod2']:07b} |")
print(f"\nS8 pair {c['s8_pair_fraction']:.5f}, any {c['s8_any_lagrangian_fraction']:.5f}, S7 pair {c['s7_pair_fraction']:.5f}")
