# Benchmarks

The visual record of every performance measurement and speedup: report PDFs,
graphs (SVG) and PNG renders, one directory per entry
(`<YYYYMMDD>-<slug>/`). The study directory linked from each row holds the
raw records, sources and inputs and is the source of truth; the files here
are copies for browsing (AGENTS.md "Performance results also land in
`benchmarks/`"). Not to be confused with `benches/` (Cargo bench sources).

| date | entry | measured | headline (as reported) | host | study |
|:--|:--|:--|:--|:--|:--|
| 2026-10-12 | [`20261012-n83-exact-scalar-rho`](20261012-n83-exact-scalar-rho/) | Exact scalar candidate and expanded one-target rho lane/DP sweep on `icv1-f2m83-t6151469093347-9eba6cd5`; matched already-tuned Fold64 control | Selected reference/Fold256/DP6: held-out process ratio 1.069651 [1.043178,1.089238], means 2.987008 s versus 2.780545 s; five fresh targets × three paired seeds. Barrett correctness passed but training runtime did not improve. All139 runs and independent answers verified | AMD EPYC 9V74 x86_64 KVM; CPU3,one OS thread,native core reservation,0 contended campaign samples; physical host not exclusive | [`research/n83_speed_20261012`](../research/n83_speed_20261012/) · [run/results PDF](../research/n83_speed_20261012/reports/run-20261012T1543Z/n83_faster_rho_evaluation.pdf) |
| 2026-10-11 | [`20261011-groebner-vs-msolve-singular-polybori`](20261011-groebner-vs-msolve-singular-polybori/) | `f4_fp` and the inherited-F4 / F6-IC decomposition oracle against msolve 0.6.5, Singular 4.3.2 and PolyBoRi on identical exported systems (Magma scripts included, not run) | F_p: same reduced basis as Singular on 13/13 systems; msolve 1-thread 12–24× faster than ours on square generic systems (quad_n9: 246 ms vs 5 693 ms), level on overdetermined ones. Boolean n = 17, 8 targets: ours decides in 7–24 ms per target, msolve 42.7–46.5 s and PolyBoRi 40.4–52.3 s for the full basis, Singular > 120 s; verdicts agree | 4 vCPU Xeon 2.10 GHz, shared, L0 | [`research/gb_cas_comparison_20261011`](../research/gb_cas_comparison_20261011/) (PR #1675), run-20261011T1806Z |
| 2026-10-11 | [`20261011-inherited-f4-eager-rewrite`](20261011-inherited-f4-eager-rewrite/) | inherited F4 eager rewrite of kept rows against lazy rows (`KIC_F4_EAGER=0`), same binary, `f6_ic_probe` | lazy/eager median ratio: K₀/2²³ dim 8 inherited F4 1.40, F6-IC 1.62; K₀/2³¹ dim 9 1.48 and 1.53; counters and verdicts identical | 4 vCPU Xeon 2.10 GHz, shared, L0 | [`research/inherited_f4_eager_rewrite_20261011`](../research/inherited_f4_eager_rewrite_20261011/) (PR #1671), run-20261011T1657Z |
