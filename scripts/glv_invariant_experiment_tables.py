#!/usr/bin/env python3
"""Print the tables of RESEARCH_GLV_INVARIANT_FACTOR_BASES.md §6 and §8 (E1–E9)
from the frozen experiment files.

    python3 scripts/glv_invariant_experiment_tables.py experiments/23_glv_invariant_e1.json ...

Every number in the note's §6 is printed by this script; the note cites,
it does not compute.  A file's `rows[*].experiment` selects its table.
"""
import json
import math
import sys
from collections import defaultdict
from statistics import mean


def load(paths):
    rows = []
    for p in paths:
        rows.extend(json.load(open(p))["rows"])
    return rows


def fmt(x, nd=2):
    if x is None:
        return "—"
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    return str(x)


def arm(s, k):
    return s[k]


def e1(rows):
    print("### E1 — automorphism fold at the square point, one stream, both arms\n")
    print("| family | log2 r | h | oracle | seed | cols fold | cols control | deficiency fold / control | full-rank rel fold | full-rank rel control | ratio | full-rank rel / cols fold | full-rank rel / cols control | rank fraction at k = cols, fold / control | first-pin rel fold / control | trials | hit rate | correct |")
    print("|:--|--:|--:|:--|--:|--:|--:|:--|--:|--:|--:|--:|--:|:--|:--|--:|--:|:--|")
    for r in sorted(rows, key=lambda r: (r["family"], r["oracle"], r["log2_r"], r["seed"])):
        s = r["stream"]
        f, c = s["folded"], s["control"]
        print(
            f"| {r['family']} | {r['log2_r']:.1f} | {r['cofactor']} | {r['oracle']} | {r['seed']} | {f['columns']} | {c['columns']} | "
            f"{f['deficiency_total']} / {c['deficiency_total']} | {fmt(f['square_relations'])} | {fmt(c['square_relations'])} | "
            f"{fmt(s['square_ratio'])} | {fmt(f['square_over_columns'])} | {fmt(c['square_over_columns'])} | "
            f"{fmt(f['rank_fraction_at_columns'])} / {fmt(c['rank_fraction_at_columns'])} | "
            f"{fmt(f['first_pin_relations'])} / {fmt(c['first_pin_relations'])} | {s['trials']} | {s['hit_rate']:.4f} | "
            f"{'yes' if f['verified'] and c['verified'] else 'NO'} |"
        )
    print("\n#### E1 summary by family, oracle and size (means over seeds)\n")
    print("| family | oracle | log2 r (mean) | seeds | column ratio | full-rank ratio mean (min–max) | full-rank rel / cols, fold | full-rank rel / cols, control | rank fraction at k = cols, fold | control | first-pin ratio mean | all correct |")
    print("|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|")
    groups = defaultdict(list)
    for r in rows:
        groups[(r["family"], r["oracle"], r["bits"])].append(r)
    for (fam, orc, bits), rs in sorted(groups.items()):
        sq = [r["stream"]["square_ratio"] for r in rs if r["stream"]["square_ratio"] is not None]
        fp = [r["stream"]["first_pin_ratio"] for r in rs if r["stream"]["first_pin_ratio"] is not None]
        fo = [r["stream"]["folded"]["square_over_columns"] for r in rs if r["stream"]["folded"]["square_over_columns"] is not None]
        co = [r["stream"]["control"]["square_over_columns"] for r in rs if r["stream"]["control"]["square_over_columns"] is not None]
        ok = all(r["stream"]["folded"]["verified"] and r["stream"]["control"]["verified"] for r in rs)
        col = {round(r["stream"]["column_ratio"], 2) for r in rs}
        rf = [r["stream"]["folded"]["rank_fraction_at_columns"] for r in rs if r["stream"]["folded"]["rank_fraction_at_columns"] is not None]
        rc = [r["stream"]["control"]["rank_fraction_at_columns"] for r in rs if r["stream"]["control"]["rank_fraction_at_columns"] is not None]
        print(
            f"| {fam} | {orc} | {mean(r['log2_r'] for r in rs):.1f} | {len(rs)} | {', '.join(str(c) for c in sorted(col))} | "
            f"{fmt(mean(sq)) if sq else '—'} ({fmt(min(sq)) if sq else '—'}–{fmt(max(sq)) if sq else '—'}) | {fmt(mean(fo)) if fo else '—'} | {fmt(mean(co)) if co else '—'} | "
            f"{fmt(mean(rf)) if rf else '—'} | {fmt(mean(rc)) if rc else '—'} | {fmt(mean(fp)) if fp else '—'} | {'yes' if ok else 'NO'} |"
        )


def e2(rows):
    print("### E2 — GLS line (`ψ`, order 4) against Koblitz orbit (`τ`, order n), one driver\n")
    print("| family | size | log2 r | h | eigenvalue order | oracle | seed | points | cols fold | cols control | pts/col | column ratio | deficiency fold / control | square rel fold | square rel control | ratio | square/cols fold | trials | hit rate | oracle F_p muls / call | agreement with subtract | correct |")
    print("|:--|:--|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|:--|:--|")
    for r in sorted(rows, key=lambda r: (r["family"], r["log2_r"], r["seed"])):
        s = r["stream"]
        f, c = s["folded"], s["control"]
        size = f"p=2^{r['p_bits']}" if "p_bits" in r else f"n={r['n']}"
        muls = r.get("oracle_fp_muls")
        calls = r.get("oracle_calls")
        per = f"{muls / calls:.0f}" if muls and calls else "—"
        agr = r.get("agreement_with_subtract")
        agr_s = f"{agr['agree']}/{agr['targets']} ({agr['disagree']} disagree)" if agr else "—"
        print(
            f"| {r['family']} | {size} | {r['log2_r']:.1f} | {r['cofactor']} | {r['eigenvalue_order']} | {r['oracle']} | {r['seed']} | {f['points']} | {f['columns']} | {c['columns']} | "
            f"{f['points_per_column']:.1f} | {s['column_ratio']:.2f} | {f['deficiency_total']} / {c['deficiency_total']} | {fmt(f['square_relations'])} | {fmt(c['square_relations'])} | {fmt(s['square_ratio'])} | "
            f"{fmt(f['square_over_columns'])} | {s['trials']} | {s['hit_rate']:.4f} | {per} | {agr_s} | {'yes' if f['verified'] and c['verified'] else 'NO'} |"
        )
    print("\n#### E2 summary\n")
    print("| family | eigenvalue order | instances | log2 r | points per column | column ratio | full-rank ratio mean (min–max) | rank fraction at k = cols, fold / control | oracle F_p muls / call | agreement with subtract | all correct |")
    print("|:--|--:|--:|:--|--:|--:|--:|--:|--:|:--|:--|")
    groups = defaultdict(list)
    for r in rows:
        groups[(r["family"], r["eigenvalue_order"])].append(r)
    for (fam, ord_), rs in sorted(groups.items()):
        sq = [r["stream"]["square_ratio"] for r in rs if r["stream"]["square_ratio"] is not None]
        ok = all(r["stream"]["folded"]["verified"] and r["stream"]["control"]["verified"] for r in rs)
        ppc = {round(r["stream"]["folded"]["points_per_column"], 1) for r in rs}
        col = {round(r["stream"]["column_ratio"], 2) for r in rs}
        rf = [r["stream"]["folded"]["rank_fraction_at_columns"] for r in rs if r["stream"]["folded"]["rank_fraction_at_columns"] is not None]
        rc = [r["stream"]["control"]["rank_fraction_at_columns"] for r in rs if r["stream"]["control"]["rank_fraction_at_columns"] is not None]
        lr = [r["log2_r"] for r in rs]
        per = [r["oracle_fp_muls"] / r["oracle_calls"] for r in rs if r.get("oracle_fp_muls") and r.get("oracle_calls")]
        agr = [r["agreement_with_subtract"] for r in rs if r.get("agreement_with_subtract")]
        agr_s = f"{sum(a['agree'] for a in agr)}/{sum(a['targets'] for a in agr)} ({sum(a['disagree'] for a in agr)} disagree)" if agr else "—"
        print(f"| {fam} | {ord_} | {len(rs)} | {min(lr):.1f}–{max(lr):.1f} | {', '.join(str(x) for x in sorted(ppc))} | {', '.join(str(x) for x in sorted(col))} | {fmt(mean(sq)) if sq else '—'} ({fmt(min(sq)) if sq else '—'}–{fmt(max(sq)) if sq else '—'}) | {fmt(mean(rf)) if rf else '—'} / {fmt(mean(rc)) if rc else '—'} | {f'{mean(per):.0f}' if per else '—'} | {agr_s} | {'yes' if ok else 'NO'} |")


def e3(rows):
    print("### E3 — composite groups on twisted CM curves\n")
    print("| family | p | log2 r | h | ord λ_ψ | ord λ_aut | aut = ±ψ on ⟨G⟩ | pts/col negation | pts/col ψ | pts/col aut | pts/col ψ+aut | cols ψ+aut | square ratio ψ+aut vs negation | square ratio ψ vs negation | correct |")
    print("|:--|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|:--|")
    for r in sorted(rows, key=lambda r: (r["family"], r["p_bits"], r["seed"])):
        f = r["folds"]
        sb = r.get("stream_psi_aut_vs_negation")
        sp = r.get("stream_psi_vs_negation")
        ok = (sb is None or (sb["folded"]["verified"] and sb["control"]["verified"])) and (sp is None or (sp["folded"]["verified"] and sp["control"]["verified"]))
        if sb is None and sp is None:
            ok = "folds only (r < 16·columns)"
        print(
            f"| {r['family']} | 2^{r['p_bits']} | {r['log2_r']:.1f} | {r['cofactor']} | {r['psi_eigenvalue_order']} | {r['aut_eigenvalue_order']} | {r['aut_is_plus_minus_psi_on_subgroup']} | "
            f"{f['negation']['points_per_column']:.1f} | {f['psi']['points_per_column']:.1f} | {f['aut']['points_per_column']:.1f} | {f['psi+aut']['points_per_column']:.1f} | {f['psi+aut']['columns']} | "
            f"{fmt(sb['square_ratio']) if sb else '—'} | {fmt(sp['square_ratio']) if sp else '—'} | {ok if isinstance(ok, str) else ('yes' if ok else 'NO')} |"
        )


def e4(rows):
    print("### E4 — type C: degree-2 and degree-3 CM endomorphisms against a base\n")
    print("| family | degree | log2 r | h | seed | maps found | map | trace | ord_r(λ) | images in base | base points | chance fraction | verified |")
    print("|:--|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|:--|")
    for r in sorted(rows, key=lambda r: (r["family"], r["degree"], r["log2_r"], r["seed"])):
        if not r["maps"]:
            print(f"| {r['family']} | {r['degree']} | {r['log2_r']:.1f} | {r['cofactor']} | {r['seed']} | 0 | — | — | — | — | {r['base_points']} | — | — |")
        for m in r["maps"]:
            o = m["overlap"]
            trace = m["name"].split("trace=")[1].rstrip("]") if "trace=" in m["name"] else "—"
            print(
                f"| {r['family']} | {r['degree']} | {r['log2_r']:.1f} | {r['cofactor']} | {r['seed']} | {len(r['maps'])} | `{m['name']}` | {trace} | {o['eigenvalue_order']} | {o['images_in_base']} | {o['base_points']} | {o['chance_fraction']:.5f} | {m['verified']} |"
            )
    print("\n#### E4 summary\n")
    print("| family | degree | instances | maps per instance | ord_r(λ) min–max | images in base, total | base points, total | all verified |")
    print("|:--|--:|--:|--:|--:|--:|--:|:--|")
    groups = defaultdict(list)
    for r in rows:
        groups[(r["family"], r["degree"])].append(r)
    for (fam, deg), rs in sorted(groups.items()):
        orders = [m["overlap"]["eigenvalue_order"] for r in rs for m in r["maps"]]
        inside = sum(m["overlap"]["images_in_base"] for r in rs for m in r["maps"])
        pts = sum(m["overlap"]["base_points"] for r in rs for m in r["maps"])
        ok = all(m["verified"] for r in rs for m in r["maps"]) and all(r["maps"] for r in rs)
        counts = {len(r["maps"]) for r in rs}
        print(f"| {fam} | {deg} | {len(rs)} | {', '.join(str(c) for c in sorted(counts))} | {min(orders) if orders else '—'}–{max(orders) if orders else '—'} | {inside} | {pts} | {'yes' if ok else 'NO'} |")


def e5(rows):
    print("### E5 — subfield curves on E(F_{p³}): the Frobenius line\n")
    print("| family | p | log2 r | #E(F_p) | h / #E(F_p) | ord λ_π | ord λ_ζ | ζ ∈ ⟨π⟩ on ⟨G⟩ | pts/col negation | pts/col π | pts/col ζ | pts/col π+ζ | cols π | cols negation | full-rank rel π | full-rank rel negation | ratio | rank fraction at k = cols, π / negation | trials | hit rate | oracle muls/call | degenerate | correct |")
    print("|:--|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|:--|")
    for r in sorted(rows, key=lambda r: (r["family"], r["p_bits"], r["seed"])):
        f = r["folds"]
        st = r.get("relation_stream")
        if st:
            s = st["stream"]
            fo, co = s["folded"], s["control"]
            per = f"{st['oracle_fp_muls'] / st['oracle_calls']:.0f}" if st["oracle_calls"] else "—"
            tail = (f"{fo['columns']} | {co['columns']} | {fmt(fo['square_relations'])} | {fmt(co['square_relations'])} | {fmt(s['square_ratio'])} | "
                    f"{fmt(fo['rank_fraction_at_columns'])} / {fmt(co['rank_fraction_at_columns'])} | {s['trials']} | {s['hit_rate']:.5f} | {per} | {st['oracle_degenerate']} | {'yes' if fo['verified'] and co['verified'] else 'NO'} |")
        else:
            tail = f"{f['frobenius']['columns']} | {f['negation']['columns']} | — | — | — | — | — | — | — | — | not an index calculus: the base is the subgroup |"
        print(
            f"| {r['family']} | 2^{r['p_bits']} | {r['log2_r']:.1f} | {r['base_order']} | {r['cofactor'] // r['base_order']} | {r['frobenius_eigenvalue_order']} | {fmt(r.get('zeta_eigenvalue_order'))} | {fmt(r.get('zeta_is_frobenius_power'))} | "
            f"{f['negation']['points_per_column']:.1f} | {f['frobenius']['points_per_column']:.1f} | {f['zeta']['points_per_column'] if 'zeta' in f else '—'} | {f['frobenius+zeta']['points_per_column'] if 'frobenius+zeta' in f else '—'} | " + tail
        )
    print("\n#### E5 summary\n")
    print("| family | instances | log2 r | pts/col negation | pts/col π | pts/col ζ | pts/col π+ζ | ζ ∈ ⟨π⟩ on ⟨G⟩ | streams | full-rank ratio mean (min–max) | rank fraction at k = cols, π / negation | oracle muls/call | degenerate | all correct |")
    print("|:--|--:|:--|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|:--|")
    groups = defaultdict(list)
    for r in rows:
        groups[r["family"]].append(r)
    for fam, rs in sorted(groups.items()):
        f0 = [r["folds"] for r in rs]
        ppc = lambda k: ", ".join(str(x) for x in sorted({round(f[k]["points_per_column"], 1) for f in f0 if k in f})) or "—"
        zf = {r.get("zeta_is_frobenius_power") for r in rs}
        zf_s = "—" if zf == {None} else ("all" if zf == {True} else "NO")
        sts = [r["relation_stream"] for r in rs if r.get("relation_stream")]
        sq = [st["stream"]["square_ratio"] for st in sts if st["stream"]["square_ratio"] is not None]
        rf = [st["stream"]["folded"]["rank_fraction_at_columns"] for st in sts if st["stream"]["folded"]["rank_fraction_at_columns"] is not None]
        rc = [st["stream"]["control"]["rank_fraction_at_columns"] for st in sts if st["stream"]["control"]["rank_fraction_at_columns"] is not None]
        per = [st["oracle_fp_muls"] / st["oracle_calls"] for st in sts if st["oracle_calls"]]
        deg = sum(st["oracle_degenerate"] for st in sts)
        ok = all(st["stream"]["folded"]["verified"] and st["stream"]["control"]["verified"] for st in sts)
        lr = [r["log2_r"] for r in rs]
        print(f"| {fam} | {len(rs)} | {min(lr):.1f}–{max(lr):.1f} | {ppc('negation')} | {ppc('frobenius')} | {ppc('zeta')} | {ppc('frobenius+zeta')} | {zf_s} | {len(sts)} | {fmt(mean(sq)) if sq else '—'} ({fmt(min(sq)) if sq else '—'}–{fmt(max(sq)) if sq else '—'}) | {fmt(mean(rf)) if rf else '—'} / {fmt(mean(rc)) if rc else '—'} | {f'{mean(per):.0f}' if per else '—'} | {deg if sts else '—'} | {('yes' if ok else 'NO') if sts else 'folds only'} |")


def e6(rows):
    print("### E6 — the matched rho: folded by the automorphism group against the negation walk\n")
    print("| family | log2 r | h | seed | walks | S negation (mean) | S folded (mean) | S ratio | steps ratio | expected √(A/2) | all verified |")
    print("|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|")
    for r in sorted(rows, key=lambda r: (r["family"], r["log2_r"], r["seed"])):
        print(f"| {r['family']} | {r['log2_r']:.1f} | {r['cofactor']} | {r['seed']} | {r['rho_runs']} | {r['negation_s_mean']:.3f} | {r['folded_s_mean']:.3f} | {r['s_ratio']:.2f} | {r['steps_ratio']:.2f} | {r['expected_ratio']:.2f} | {r['all_verified']} |")
    print("\n#### E6 summary\n")
    print("| family | instances | walks | S ratio mean (min–max) | steps ratio mean | expected | all verified |")
    print("|:--|--:|--:|--:|--:|--:|:--|")
    groups = defaultdict(list)
    for r in rows:
        groups[r["family"]].append(r)
    for fam, rs in sorted(groups.items()):
        sr = [r["s_ratio"] for r in rs]
        st = [r["steps_ratio"] for r in rs]
        print(f"| {fam} | {len(rs)} | {sum(r['rho_runs'] for r in rs)} | {mean(sr):.2f} ({min(sr):.2f}–{max(sr):.2f}) | {mean(st):.2f} | {rs[0]['expected_ratio']:.2f} | {all(r['all_verified'] for r in rs)} |")


def e7(rows):
    print("### E7 — three summands and orbit duplicates on j = 0\n")
    print("| log2 r | seed | m | group order | seed abscissae | points | cols fold | cols control | square rel fold | square rel control | ratio | zero-support rows fold / control | single-column rows fold / control | orbit duplicates | targets | hit rate | correct |")
    print("|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|--:|--:|--:|:--|")
    for r in sorted(rows, key=lambda r: (r["log2_r"], r["seed"], r["summands"])):
        s = r["stream"]
        f, c = s["folded"], s["control"]
        print(
            f"| {r['log2_r']:.1f} | {r['seed']} | {r['summands']} | {r['group_order']} | {r['seed_abscissae']} | {f['points']} | {f['columns']} | {c['columns']} | {fmt(f['square_relations'])} | {fmt(c['square_relations'])} | {fmt(s['square_ratio'])} | "
            f"{f['zero_support_rows']} / {c['zero_support_rows']} | {f['single_column_rows']} / {c['single_column_rows']} | {s['orbit_duplicates']} | {s['trials']} | {s['hit_rate']:.4f} | {'yes' if f['verified'] and c['verified'] else 'NO'} |"
        )


def fit_exponent(points):
    """Least-squares slope of log2(y) against log2(x) over (x, y) pairs."""
    pts = [(math.log2(x), math.log2(y)) for x, y in points if x > 0 and y > 0]
    if len(pts) < 2:
        return None
    mx = mean(p[0] for p in pts)
    my = mean(p[1] for p in pts)
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    if sxx == 0:
        return None
    return sum((p[0] - mx) * (p[1] - my) for p in pts) / sxx


def e8_costs(r):
    """Every phase of one E8 row in F_p multiplications (inversions priced
    by the row's measured conversion; a row operation of the linear
    algebra is one Z/rZ multiplication, priced as one).  Returns a dict
    per arm and per rho walk, plus the shares used."""
    k = r["inversion_in_multiplications"]
    s = r["stream"]
    unit = lambda m, i: m + k * i
    stream_total = unit(r["stream_cost"]["fp_muls"], r["stream_cost"]["fp_invs"])
    table = unit(r["pair_table"]["fp_muls"], r["pair_table"]["fp_invs"])
    out = {}
    for arm in ("folded", "control"):
        a = s[arm]
        build = unit(r["base_build"][arm]["fp_muls"], r["base_build"][arm]["fp_invs"])
        # The stream is one uniform process feeding both arms; each arm's
        # share is the part up to its own full rank.
        share = (a["square_trials"] / s["trials"]) if (a["square_trials"] and s["trials"]) else 1.0
        stream = stream_total * share
        la = a["row_ops"]
        total = build + table + stream + la
        out[arm] = {"build": build, "table": table, "stream": stream, "la": la, "total": total, "share": share,
                    "S": total / math.sqrt(r["r"])}
    for walk in ("negation", "folded"):
        tot = [unit(w[f"{walk}_fp_muls"], w[f"{walk}_fp_invs"]) for w in r["walks"]]
        steps = [w[walk]["steps"] for w in r["walks"]]
        gae = [w[walk]["gae"] for w in r["walks"]]
        out[f"rho_{walk}"] = {"total": mean(tot), "S": mean(tot) / math.sqrt(r["r"]), "steps": mean(steps),
                              "per_op": mean(tot) / max(1.0, mean(gae)),
                              "verified": all(w[walk]["verified"] for w in r["walks"])}
    out["muls_per_add"] = r["stream_cost"]["fp_muls"] / max(1, s["group_ops"]["adds"] + s["group_ops"]["doubles"])
    return out


def e8(rows):
    print("### E8 — three summands on the Frobenius line of a subfield curve: relations\n")
    print("| p | log2 r | h / #E(F_p) | cols π | cols negation | full-rank rel π | full-rank rel negation | ratio | rank fraction at k = cols, π / negation | single-column rows π / negation | orbit duplicates | targets | hit rate | correct |")
    print("|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|--:|--:|--:|:--|")
    rows = sorted(rows, key=lambda r: (r["log2_r"], r["seed"]))
    for r in rows:
        s = r["stream"]
        f, c = s["folded"], s["control"]
        print(f"| 2^{r['p_bits']} | {r['log2_r']:.1f} | {r['cofactor'] // r['base_order']} | {f['columns']} | {c['columns']} | {fmt(f['square_relations'])} | {fmt(c['square_relations'])} | {fmt(s['square_ratio'])} | "
              f"{fmt(f['rank_fraction_at_columns'])} / {fmt(c['rank_fraction_at_columns'])} | {f['single_column_rows']} / {c['single_column_rows']} | {s['orbit_duplicates']} | {s['trials']} | {s['hit_rate']:.4f} | {'yes' if f['verified'] and c['verified'] else 'NO'} |")
    print("\n#### E8 — every phase priced, one unit (F_p multiplications; inversion = measured factor; LA row op = 1); rows that reached full rank on both arms\n")
    print("| p | log2 r | inv / mul (measured) | muls per addition | base build π / neg | pair table | stream π / neg | LA π / neg | total π | total negation | S π | S negation | S neg / S π | rho S negation | rho S folded | rho muls per group op, negation / folded | rho steps ratio (expected) | S π / rho S folded | all verified |")
    print("|--:|--:|--:|--:|:--|--:|:--|:--|--:|--:|--:|--:|--:|--:|--:|:--|:--|--:|:--|")
    fits = defaultdict(list)
    full = [r for r in rows if r["stream"]["folded"]["square_relations"] and r["stream"]["control"]["square_relations"]]
    for r in rows:
        if r not in full:
            print(f"| 2^{r['p_bits']} | {r['log2_r']:.1f} | — | — | not full rank on both arms: the stream exhausted the group's targets | | | | | | | | | | | | | | |")
            continue
        c = e8_costs(r)
        f, n = c["folded"], c["control"]
        rn, rf = c["rho_negation"], c["rho_folded"]
        ok = r["stream"]["folded"]["verified"] and r["stream"]["control"]["verified"] and rn["verified"] and rf["verified"]
        for key, val in (("folded", f["total"]), ("control", n["total"]), ("rho_folded", rf["total"]), ("rho_negation", rn["total"])):
            fits[key].append((r["r"], val))
        print(f"| 2^{r['p_bits']} | {r['log2_r']:.1f} | {r['inversion_in_multiplications']:.1f} | {c['muls_per_add']:.1f} | {f['build']:.3g} / {n['build']:.3g} | {f['table']:.3g} | {f['stream']:.3g} / {n['stream']:.3g} | {f['la']:.3g} / {n['la']:.3g} | {f['total']:.3g} | {n['total']:.3g} | {f['S']:.1f} | {n['S']:.1f} | {n['S'] / f['S']:.2f} | {rn['S']:.1f} | {rf['S']:.1f} | {rn['per_op']:.0f} / {rf['per_op']:.0f} | {r['rho_steps_ratio']:.2f} ({r['rho_expected_ratio']:.2f}) | {f['S'] / rf['S']:.0f} | {'yes' if ok else 'NO'} |")
    print("\n#### E8 summary (rows at full rank on both arms)\n")
    print("| arms | instances (of run) | log2 r | column ratio | full-rank ratio mean (min–max) | rank fraction at k = cols, π / negation | S neg / S π mean (min–max) | S π / rho S folded, min–max | rho steps ratio mean (expected) | rho muls per group op, negation / folded | rho S negation / rho S folded, mean | fitted exponent of total: π, negation, rho folded, rho negation (rho: 0.50) | all verified |")
    print("|:--|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|")
    sq = [r["stream"]["square_ratio"] for r in full]
    rf_ = [r["stream"]["folded"]["rank_fraction_at_columns"] for r in full if r["stream"]["folded"]["rank_fraction_at_columns"] is not None]
    rc_ = [r["stream"]["control"]["rank_fraction_at_columns"] for r in full if r["stream"]["control"]["rank_fraction_at_columns"] is not None]
    costs = [e8_costs(r) for r in full]
    s_ratio = [c["control"]["S"] / c["folded"]["S"] for c in costs]
    vs_rho = [c["folded"]["S"] / c["rho_folded"]["S"] for c in costs]
    rho_ratio = [c["rho_negation"]["S"] / c["rho_folded"]["S"] for c in costs]
    ok = all(r["stream"]["folded"]["verified"] and r["stream"]["control"]["verified"] and r["rho_all_verified"] for r in rows)
    lr = [r["log2_r"] for r in full]
    col = {round(r["stream"]["column_ratio"], 2) for r in full}
    exps = ", ".join(f"{fmt(fit_exponent(fits[k]))}" for k in ("folded", "control", "rho_folded", "rho_negation"))
    print(f"| π-line fold vs negation, m = 3, pair table | {len(full)} ({len(rows)}) | {min(lr):.1f}–{max(lr):.1f} | {', '.join(str(x) for x in sorted(col))} | {fmt(mean(sq))} ({fmt(min(sq))}–{fmt(max(sq))}) | {fmt(mean(rf_)) if rf_ else '—'} / {fmt(mean(rc_)) if rc_ else '—'} | {mean(s_ratio):.2f} ({min(s_ratio):.2f}–{max(s_ratio):.2f}) | {min(vs_rho):.0f}–{max(vs_rho):.0f} | {mean(r['rho_steps_ratio'] for r in full):.2f} ({full[0]['rho_expected_ratio']:.2f}) | {mean(c['rho_negation']['per_op'] for c in costs):.0f} / {mean(c['rho_folded']['per_op'] for c in costs):.0f} | {mean(rho_ratio):.2f} | {exps} | {'yes' if ok else 'NO'} |")


def e9(rows):
    print("### E9 — the matched folded rho on the F_{p²} and F_{p³} groups\n")
    print("| family | p | log2 r | generators | A | walks | S negation (mean) | S folded (mean) | S ratio | steps ratio | expected √(A/2) | all verified |")
    print("|:--|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|:--|")
    for r in sorted(rows, key=lambda r: (r["family"], r["log2_r"], r["seed"])):
        print(f"| {r['family']} | 2^{r['p_bits']} | {r['log2_r']:.1f} | {', '.join(r['generators'])} | {r['group_order_of_fold']} | {r['rho_runs']} | {r['negation_s_mean']:.3f} | {r['folded_s_mean']:.3f} | {r['s_ratio']:.2f} | {r['steps_ratio']:.2f} | {r['expected_ratio']:.2f} | {r['all_verified']} |")
    print("\n#### E9 summary\n")
    print("| family | A | instances | log2 r | walks | S ratio mean (min–max) | steps ratio mean | expected | all verified |")
    print("|:--|--:|--:|:--|--:|--:|--:|--:|:--|")
    groups = defaultdict(list)
    for r in rows:
        groups[(r["family"], r["group_order_of_fold"])].append(r)
    for (fam, a), rs in sorted(groups.items()):
        sr = [r["s_ratio"] for r in rs]
        lr = [r["log2_r"] for r in rs]
        print(f"| {fam} | {a} | {len(rs)} | {min(lr):.1f}–{max(lr):.1f} | {sum(r['rho_runs'] for r in rs) * 2} | {mean(sr):.2f} ({min(sr):.2f}–{max(sr):.2f}) | {mean(r['steps_ratio'] for r in rs):.2f} | {rs[0]['expected_ratio']:.2f} | {all(r['all_verified'] for r in rs)} |")


def main():
    rows = load(sys.argv[1:])
    by = defaultdict(list)
    for r in rows:
        by[r["experiment"]].append(r)
    for exp, fn in [("e1", e1), ("e2", e2), ("e3", e3), ("e4", e4), ("e5", e5), ("e6", e6), ("e7", e7), ("e8", e8), ("e9", e9)]:
        if exp in by:
            fn(by[exp])
            print()


if __name__ == "__main__":
    main()
