"""Print the report's tables as Typst rows from the recorded runs (data/*.json).

Rendering only: every number is a field of a run's `report` object.

    python3 graphs/make_tables.py > rows.txt        # from this study directory

Sections are separated by a line starting with `//`.
"""
import json
import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data"
SUP = str.maketrans("0123456789-.", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻·")


def p2(x, nd=1):
    """2^x with superscript digits."""
    return "2" + f"{x:.{nd}f}".translate(SUP)


def load(path):
    return json.loads(path.read_text())["report"]


def curve_runs(pattern):
    rows = []
    for p in sorted(DATA.glob(pattern)):
        r = load(p)
        name = p.stem.removeprefix("s3_").removeprefix("raw_").rsplit("_", 1)[0]
        rows.append((r["field_bits"], name, r))
    rows.sort(key=lambda t: (t[0], t[1]))
    return rows


def fmt_us(us):
    return f"{us / 1000:.2f} ms" if us >= 1000 else f"{us:.1f} µs"


print("// SECG curves at B = 4096")
for _, name, r in curve_runs("sec[pt]*_4096.json"):
    print(
        f"    [`{name}`], [{r['field_bits']}], [{r['cofactor']}], [{r['order_bits']:.1f}], "
        f"[{p2(r['exact_yield_log2'])}], [{r['collisions']}], "
        f"[{p2(r['trials_for_full_system_log2'])}], [{p2(r['rho_log2'])}], "
        f"[{r['timing_ms']['pair_sums'] / 1000:.2f}], [{r['mc_trials'] // 1000}k], [{r['controls_recovered']}/{r['controls']}],"
    )

print("\n// SECG curves at B = 8192")
for _, name, r in curve_runs("sec[pt]*_8192.json"):
    print(
        f"    [`{name}`], [{r['field_bits']}], [{r['cofactor']}], [{p2(r['exact_yield_log2'])}], [{r['collisions']}], "
        f"[{p2(r['trials_for_full_system_log2'])}], [{p2(r['rho_log2'])}], "
        f"[{r['timing_ms']['pair_sums'] / 1000:.1f}], [{r['mc_hits']}], [{r['controls_recovered']}/{r['controls']}],"
    )

print("\n// Growth with B on three curves")
for name in ("sect113r1", "sect163r2", "secp160r1"):
    for b in (256, 1024, 4096, 8192):
        p = DATA / f"{name}_{b}.json"
        if not p.exists():
            continue
        r = load(p)
        print(
            f"    [`{name}`], [{b}], [{p2(math.log2(r['sums']), 0)}], [{r['collisions']}], "
            f"[{p2(r['exact_yield_log2'])}], [{p2(r['trials_for_full_system_log2'])}], "
            f"[{r['timing_ms']['pair_sums'] / 1000:.2f}], [{r['mc_hits']}],"
        )

print("\n// toy curves")
for p in sorted(DATA.glob("toy*.json")):
    r = load(p)
    lo, hi = r["mc_ci95"]
    kind = "binary" if p.stem.startswith("toyb") else "prime"
    print(
        f"    [{kind}], [{r['cofactor']}], [{r['order']}], [{r['base_abscissae']}], [{r['collisions']}], "
        f"[{r['exact_yield']:.5g}], [{r['model_yield']:.5g}], [{r['mc_rate']:.5g}], [{r['mc_hits']}], "
        f"[{r['mc_over_exact']:.4f}], [{'yes' if lo <= r['exact_yield'] <= hi else 'NO'}],"
    )

print("\n// S3 cost at B = 4096")
for _, name, r in curve_runs("s3_sec[pt]*_4096.json"):
    s = r["s3"]
    print(
        f"    [`{name}`], [{fmt_us(s['s3_us_per_trial'])}], [{s['s3_field_muls_per_trial'] / 1000:.0f}k], "
        f"[{s['table_ns_per_trial']:.0f} ns], [{s['group_add_us']:.1f} µs], "
        f"[{p2(s['s3_full_system_log2_s'])}], [{p2(s['rho_log2_s'])}], [{s['bits_slower_than_rho']:.1f}], "
        f"[{s['planted_found']}/{s['planted']}], [{s['disagreements']}],"
    )

print("\n// S3 cost against B")
for name in ("sect113r1", "sect163r2", "secp160r1", "secp256r1"):
    for b in (1024, 4096):
        p = DATA / f"s3_{name}_{b}.json"
        if not p.exists():
            continue
        s = load(p)["s3"]
        print(f"    [`{name}`], [{b}], [{fmt_us(s['s3_us_per_trial'])}], [{s['s3_us_per_trial'] * 1000 / b:.0f} ns], [{s['bits_slower_than_rho']:.1f}],")

print("\n// raw small-x base")
for _, name, r in curve_runs("raw_*.json"):
    w = r["raw_base"]
    print(
        f"    [`{name}`], [{r['cofactor']}], [{w['samples']}], [{w['in_subgroup']}], [{w['fraction']:.4f}], "
        f"[{w['ci95'][0]:.4f} to {w['ci95'][1]:.4f}], [{w['one_over_h']:.4f}], [{p2(w['raw_yield_log2'])}],"
    )
