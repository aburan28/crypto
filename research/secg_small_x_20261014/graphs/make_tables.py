"""Print the report's tables as Typst rows from the recorded runs (data/*.json).

Rendering only: every number is a field of a run's `report` object.

    python3 graphs/make_tables.py > /dev/stdout     # from this study directory
"""
import json
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data"
SUP = str.maketrans("0123456789-.", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻·")


def p2(x, nd=2):
    """2^x with superscript digits."""
    return "2" + f"{x:.{nd}f}".translate(SUP)


def load(stem):
    return json.loads((DATA / f"{stem}.json").read_text())["report"]


def secg_rows(b):
    rows = []
    for p in sorted(DATA.glob(f"sec[pt]*_{b}.json")):
        r = json.loads(p.read_text())["report"]
        name = p.stem.rsplit("_", 1)[0]
        rows.append((r["field_bits"], name, r))
    rows.sort(key=lambda t: (t[0], t[1]))
    return rows


print("// SECG curves at B = 4096")
for _, name, r in secg_rows(4096):
    print(
        f"    [`{name}`], [{r['field_bits']}], [{r['cofactor']}], [{r['order_bits']:.1f}], "
        f"[{p2(r['exact_yield_log2'], 1)}], [{r['collisions']}], "
        f"[{p2(r['trials_for_full_system_log2'], 1)}], [{p2(r['rho_log2'], 1)}], "
        f"[{r['timing_ms']['pair_sums'] / 1000:.2f}], [{r['timing_ms']['us_per_trial_wall']:.1f}], "
        f"[{r['controls_recovered']}/{r['controls']}],"
    )

print("\n// Growth with B on three curves")
for name in ("sect113r1", "sect163r2", "secp160r1"):
    for b in (256, 1024, 4096, 8192):
        try:
            r = load(f"{name}_{b}")
        except FileNotFoundError:
            continue
        print(
            f"    [`{name}`], [{b}], [{p2(__import__('math').log2(r['sums']), 0)}], [{r['collisions']}], "
            f"[{p2(r['exact_yield_log2'], 1)}], [{p2(r['trials_for_full_system_log2'], 1)}], "
            f"[{r['timing_ms']['pair_sums'] / 1000:.2f}], [{r['mc_hits']}],"
        )

print("\n// toy binary curves")
for p in sorted(DATA.glob("toyb*.json")):
    r = json.loads(p.read_text())["report"]
    lo, hi = r["mc_ci95"]
    print(
        f"    [{r['order']}], [{r['base_abscissae']}], [{r['collisions']}], [{r['exact_yield']:.5g}], "
        f"[{r['model_yield']:.5g}], [{r['mc_rate']:.5g}], [{r['mc_hits']}], "
        f"[{r['mc_over_exact']:.4f}], [{'yes' if lo <= r['exact_yield'] <= hi else 'NO'}],"
    )
