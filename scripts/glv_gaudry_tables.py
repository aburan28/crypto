#!/usr/bin/env python3
"""Tables of `RESEARCH_GLV_INDEX_CALCULUS.md` from the frozen experiment files.

    python3 scripts/glv_gaudry_tables.py experiments/22_glv_quotient.json \
        experiments/22_glv_canonical.json experiments/22_glv_graded.json \
        experiments/22_glv_invariant.json

Every number in the note's tables is printed by this script from the JSON the
bench wrote; the note quotes, it does not compute.
"""
import json
import statistics
import sys


def load(path):
    with open(path) as f:
        return json.load(f)


def fmt(v, digits=1):
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:,.{digits}f}"
    return f"{v:,}"


def quotient(rows):
    print("### Experiment 1 — factor-base quotient (one residual stream, two stores)\n")
    print("| p | n | seed | store | columns | relations | residuals | NNZ | bytes | core rows × NNZ | Wiedemann ms | LA mults | total ops | S | residuals / floor | S control / S quotient | correct |")
    print("|---:|:--|---:|:--|---:|---:|---:|---:|---:|:--|---:|---:|---:|---:|---:|---:|:--|")
    per_size = {}
    for r in rows:
        q = r["report"]
        c, t = q["control"], q["quotient"]
        for st in (c, t):
            print(
                f"| {q['p']} | 2^{q['bits']:.1f} | {r['seed']} | {st['label']} | {st['columns']} | {st['relations']} | {st['residuals_at_solve']} | {st['nnz']} | {st['bytes']:,} | {st['core_rows']} × {st['core_nnz']} | {st['wiedemann_ms']:.1f} | {st['la_ops']:,} | {st['total_ops']:,.0f} | {st['s']:,.1f} | {st['residual_ratio']:.2f} | {c['s'] / t['s']:.2f} | {'yes' if st['correct'] else 'NO'} |"
            )
        per_size.setdefault(q["p"], []).append(q)
    print("\n### Experiment 1 — per size (mean over seeds)\n")
    print("| p | n | base | orbits | rate ρ | control S | quotient S | ratio | control res/floor | quotient res/floor | rho S (plain, measured) | quotient S / rho S | folded rho S (÷√3, extrapolation) | quotient S / folded rho |")
    print("|---:|:--|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for p, qs in sorted(per_size.items()):
        m = lambda k: statistics.mean(k(q) for q in qs)
        cs, ts = m(lambda q: q["control"]["s"]), m(lambda q: q["quotient"]["s"])
        rho = m(lambda q: q["rho"]["s"])
        print(
            f"| {p} | 2^{qs[0]['bits']:.1f} | {m(lambda q: q['base']):.0f} | {m(lambda q: q['orbits']):.0f} | {m(lambda q: q['decomposition_rate']):.3f} | {cs:,.1f} | {ts:,.1f} | {cs / ts:.2f} | {m(lambda q: q['control']['residual_ratio']):.2f} | {m(lambda q: q['quotient']['residual_ratio']):.2f} | {rho:.2f} | {ts / rho:,.0f}× | {rho / 3 ** 0.5:.2f} | {ts / (rho / 3 ** 0.5):,.0f}× |"
        )
    print("\nMemory: peak RSS of the whole process (both stores, MB): "
          + ", ".join(f"p={q['p']}: {q['peak_rss_bytes'] >> 20}" for q in rows_reports(rows)))


def rows_reports(rows):
    return [r["report"] for r in rows]


def canonical(rows):
    print("\n### Experiment 2 — canonical relation generation\n")
    print("| p | n | seed | stream | residuals | distinct orbits | duplicates | expected (uniform) | solver calls saved | verified / mismatches | rows produced | rows zero (ψ-trivial) | rows unit-duplicate | informative rows | merged entries |")
    print("|---:|:--|---:|:--|---:|---:|---:|---:|---:|:--|---:|---:|---:|---:|---:|")
    for r in rows:
        q = r["report"]
        for st in (q["random"], q["pairs"]):
            inf = st["rows_produced"] - st["rows_zero"] - st["rows_duplicate"]
            print(
                f"| {q['p']} | 2^{q['bits']:.1f} | {r['seed']} | {st['label']} | {st['residuals']} | {st['distinct_orbits']} | {st['duplicates']} | {st['expected_duplicates_uniform']:.4f} | {st['solver_calls_saved']} | {st['verified_duplicates']} / {st['verification_mismatches']} | {st['rows_produced']} | {st['rows_zero']} | {st['rows_duplicate']} | {inf} | {st['merged_entries']} |"
            )


def sweep_line(s):
    last = s["profiles"][-1] if s["profiles"] else None
    d = s["solve_degree"]
    return (
        f"{d if d is not None else 'none ≤ ' + str(last['degree'] if last else '?')}",
        f"{s['quotient_dim'] if s['quotient_dim'] is not None else (str(last['standard']) + ' and rising' if last else '—')}",
        f"{s['peak_rows']:,} × {s['peak_cols']:,}",
        f"{s['peak_cells']:,}",
        f"{s['solve_muls']:,}" if d is not None else f"{s['total_muls']:,} (sum to {last['degree']})" if last else "—",
    )


def graded(rows):
    print("\n### Experiment 3 — Z/3-graded Macaulay matrix\n")
    for r in rows:
        q = r["report"]
        print(f"p = {q['p']}, n = 2^{q['bits']:.1f}: symmetrised S₄ has {q['s4_terms']} terms, weight histogram {q['s4_weight_histogram']} (ψ-invariant: {q['s4_weight_homogeneous']}); {q['residuals_skipped']} residuals without an F_p-rational decomposition skipped.\n")
    print("| p | residual | harness triples | harness F_p mults | system | coupled rows | D_solve | quotient dim | peak rows × cols | dense cells | F_p mults at D_solve |")
    print("|---:|---:|---:|---:|:--|---:|:--|:--|:--|---:|---:|")
    for r in rows:
        q = r["report"]
        for row in q["rows"]:
            for name, s in (("ordinary", row["ordinary"]), ("orbit, plain", row["orbit_vanilla"]), ("orbit, block-aware", row["orbit_block"])):
                d, dim, peak, cells, muls = sweep_line(s)
                coupled = f"{row['ordinary_coupled_fraction']:.2f}" if name == "ordinary" else "0.00"
                print(f"| {q['p']} | {row['residual_index']} | {row['harness_triples']} | {row['harness_fp_muls']:,} | {name} | {coupled} | {d} | {dim} | {peak} | {cells} | {muls} |")
    print("\n#### Experiment 3 — F4 (degree-by-degree, with substitution solving)\n")
    print("| p | residual | system | D_solve | D_reached | max rows × cols | F4 runs | field mults | solutions | blocked steps | ms | mults vs ordinary | mults block / plain | ms block / plain |")
    print("|---:|---:|:--|---:|---:|:--|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        q = r["report"]
        for row in q["rows"]:
            o = row["ordinary_f4"]
            for s in (row["ordinary_f4"], row["orbit_f4_vanilla"], row["orbit_f4_block"]):
                vs = f"{s['field_ops'] / o['field_ops']:.1f}×" if o["field_ops"] else "—"
                bp = f"{row['f4_ops_ratio_block_vs_vanilla']:.3f}" if s is row["orbit_f4_block"] and row["f4_ops_ratio_block_vs_vanilla"] is not None else "—"
                bm = f"{row['f4_ms_ratio_block_vs_vanilla']:.2f}" if s is row["orbit_f4_block"] and row["f4_ms_ratio_block_vs_vanilla"] is not None else "—"
                print(f"| {q['p']} | {row['residual_index']} | {s['system']} | {s['solving_degree']} | {s['degree_reached']} | {s['max_rows']:,} × {s['max_cols']:,} | {s['f4_runs']} | {s['field_ops']:,} | {s['solutions']} | {s['blocked_steps']} | {s['ms']:.0f}{' (timed out)' if s['timed_out'] else ''} | {vs} | {bp} | {bm} |")
    print("\nOrbit solutions = 3 × ordinary on every residual: " + str(all(row["orbit_f4_triples_ordinary"] for r in rows for row in r["report"]["rows"])))


def invariant(rows):
    print("\n### Experiment 4 — invariant formulation against the baselines\n")
    q0 = rows[0]["report"]
    print(f"Invariant-monomial generators, weights (1, 2, 0, 1) on (e₁, e₂, e₃, X): {', '.join(q0['generators_symmetrised'])}; invariant monomials per degree 0..4: {q0['hilbert_symmetrised']}.")
    print(f"Diagonal action, weights (1, 1, 1, 1) on (x₁, x₂, x₃, X): {len(q0['generators_diagonal'])} generators (all cubic), invariant monomials per degree 0..4: {q0['hilbert_diagonal']}.\n")
    print("| p | residual | harness triples | harness e-solutions | harness F_p mults | invariant uses c_R = x_R³ only | invariant ≡ ordinary | ordinary Macaulay D / dim / mults | invariant Macaulay D / dim / mults | function-first Macaulay (≤ D) peak rows × cols / mults | function-first witnessed |")
    print("|---:|---:|---:|---:|---:|:--|:--|:--|:--|:--|:--|")
    for r in rows:
        q = r["report"]
        for row in q["rows"]:
            om, im, fm = row["ordinary_macaulay"], row["invariant_macaulay"], row["function_first_macaulay"]
            last = fm["profiles"][-1] if fm["profiles"] else None
            ff = f"D ≤ {last['degree']}: {fm['peak_rows']:,} × {fm['peak_cols']:,}, standard {last['standard']}, {fm['total_muls']:,} mults, no closure" if last and fm["solve_degree"] is None else f"D = {fm['solve_degree']}"
            print(f"| {q['p']} | {row['residual_index']} | {len(row['harness_triples']) if row['harness_triples'] is not None else '—'} | {row['harness_e_solutions']} | {row['harness_fp_muls']:,} | {row['invariant_uses_orbit_invariant_only']} | {row['invariant_identical_to_ordinary']} | {om['solve_degree']} / {om['quotient_dim']} / {om['solve_muls']:,} | {im['solve_degree']} / {im['quotient_dim']} / {im['solve_muls']:,} | {ff} | {row['function_first_witnessed']} ({row['function_first_witnesses']}) |")
    print("\n#### Experiment 4 — F4 on the four formulations\n")
    print("| p | residual | system | unknowns | equations | D_solve | D_reached | max rows × cols | field mults | solutions | basis | ms |")
    print("|---:|---:|:--|---:|---:|---:|---:|:--|---:|---:|---:|---:|")
    for r in rows:
        q = r["report"]
        for row in q["rows"]:
            for s in (row["ordinary_f4"], row["invariant_f4"], row["orbit_f4"], row["function_first_f4"]):
                print(f"| {q['p']} | {row['residual_index']} | {s['system']} | {s['n_vars']} | {s['equations']} | {s['solving_degree']} | {s['degree_reached']} | {s['max_rows']:,} × {s['max_cols']:,} | {s['field_ops']:,} | {s['solutions'] if s['solutions'] is not None else 'not solved'} | {s['basis_size'] if s['basis_size'] is not None else '—'} | {s['ms']:.0f}{' (timed out at budget ' + str(q['f4_budget_secs']) + ' s)' if s['timed_out'] else ''} |")
    if "ordinary_f4_basis" in rows[0]["report"]["rows"][0]:
        print("\n#### Experiment 4 — complete grevlex Gröbner bases (basis-only F4, no substitution runs)\n")
        print("| p | residual | system | unknowns | basis elements | staircase | highest productive degree | degree reached | max rows × cols | field mults | mults vs ordinary | ms |")
        print("|---:|---:|:--|---:|---:|---:|---:|---:|:--|---:|---:|---:|")
        for r in rows:
            q = r["report"]
            for row in q["rows"]:
                o = row["ordinary_f4_basis"]
                for s in (row["ordinary_f4_basis"], row["invariant_f4_basis"], row["orbit_f4_basis"], row["function_first_f4_basis"]):
                    vs = f"{s['field_ops'] / o['field_ops']:.1f}×" if o["field_ops"] else "—"
                    print(f"| {q['p']} | {row['residual_index']} | {s['system']} | {s['n_vars']} | {s['basis_size']} | {s['staircase']} | {s['solving_degree_max']} | {s['degree_reached']} | {s['max_rows']:,} × {s['max_cols']:,} | {s['field_ops']:,} | {vs} | {s['ms']:.0f}{' (timed out)' if s['timed_out'] else ''} |")
    checks = [(row["ordinary_f4_matches_harness"], row["orbit_f4_triples_ordinary"], row["function_first_witnessed"]) for r in rows for row in r["report"]["rows"]]
    print(f"\nChecks on every residual — ordinary F4 solutions = harness e-solutions: {all(c[0] for c in checks)}; orbit = 3 × ordinary: {all(c[1] for c in checks)}; function-first witnessed by every harness triple: {all(c[2] for c in checks)}.")


def main():
    files = sys.argv[1:]
    for path in files:
        rows = load(path)
        rep = rows[0]["report"]
        if "quotient" in rep:
            quotient(rows)
        elif "pairs" in rep:
            canonical(rows)
        elif "s4_weight_histogram" in rep:
            graded(rows)
        elif "generators_symmetrised" in rep:
            invariant(rows)
        else:
            print(f"unrecognised report in {path}")


if __name__ == "__main__":
    main()
