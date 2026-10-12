# Benchmarks

The visual record of every performance measurement and speedup: report PDFs,
graphs (SVG) and PNG renders, one directory per entry
(`<YYYYMMDD>-<slug>/`). The study directory linked from each row holds the
raw records, sources and inputs and is the source of truth; the files here
are copies for browsing (AGENTS.md "Performance results also land in
`benchmarks/`"). Not to be confused with `benches/` (Cargo bench sources).

| date | entry | measured | headline (as reported) | host | study |
|:--|:--|:--|:--|:--|:--|
| 2026-10-12 | [`20261012-f4-fp-vs-msolve`](20261012-f4-fp-vs-msolve/) | `f4_fp` (F_p F4) before and after Gebauer–Möller, matrix interreduction, packed monomials and parallel symbolic preprocessing, against msolve 0.6.5 and Singular 4.3.2 on 16 random dense systems | new ÷ main 15–54× faster on square systems (1 thread); msolve ÷ new 0.96–1.42 (1 thread) and 0.83–1.52 (4 threads) on systems above 50 ms; bases identical and externally confirmed | 4 vCPU Xeon 2.80 GHz, shared, L0 | [`research/f4_fp_msolve_parity_20261011`](../research/f4_fp_msolve_parity_20261011/), run-20261011T2359Z |
| 2026-10-11 | [`20261011-groebner-vs-msolve-singular-polybori`](20261011-groebner-vs-msolve-singular-polybori/) | `f4_fp` and the inherited-F4 / F6-IC decomposition oracle against msolve 0.6.5, Singular 4.3.2 and PolyBoRi on identical exported systems (Magma scripts included, not run) | F_p: same reduced basis as Singular on 13/13 systems; msolve 1-thread 12–24× faster than ours on square generic systems (quad_n9: 246 ms vs 5 693 ms), level on overdetermined ones. Boolean n = 17, 8 targets: ours decides in 7–24 ms per target, msolve 42.7–46.5 s and PolyBoRi 40.4–52.3 s for the full basis, Singular > 120 s; verdicts agree | 4 vCPU Xeon 2.10 GHz, shared, L0 | [`research/gb_cas_comparison_20261011`](../research/gb_cas_comparison_20261011/) (PR #1675), run-20261011T1806Z |
| 2026-10-11 | [`20261011-inherited-f4-eager-rewrite`](20261011-inherited-f4-eager-rewrite/) | inherited F4 eager rewrite of kept rows against lazy rows (`KIC_F4_EAGER=0`), same binary, `f6_ic_probe` | lazy/eager median ratio: K₀/2²³ dim 8 inherited F4 1.40, F6-IC 1.62; K₀/2³¹ dim 9 1.48 and 1.53; counters and verdicts identical | 4 vCPU Xeon 2.10 GHz, shared, L0 | [`research/inherited_f4_eager_rewrite_20261011`](../research/inherited_f4_eager_rewrite_20261011/) (PR #1671), run-20261011T1657Z |
