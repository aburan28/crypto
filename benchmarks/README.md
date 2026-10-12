# Benchmarks

The visual record of every performance measurement and speedup: report PDFs,
graphs (SVG) and PNG renders, one directory per entry
(`<YYYYMMDD>-<slug>/`). The study directory linked from each row holds the
raw records, sources and inputs and is the source of truth; the files here
are copies for browsing (AGENTS.md "Performance results also land in
`benchmarks/`"). Not to be confused with `benches/` (Cargo bench sources).

| date | entry | measured | headline (as reported) | host | study |
|:--|:--|:--|:--|:--|:--|
| 2026-10-12 | [`20261012-n53-s3-root-batch`](20261012-n53-s3-root-batch/) | fixed K220 N53 one-target S3 root batching, nine same-Q/rank-seed pairs and matched strong rho | Batch/control median cold 0.738 [0.714, 0.750], target-online 0.612 [0.569, 0.903]; rank inversions 763,167/80,790,012; batch/rho cold 52.586 [33.179, 64.668]; all 27 cells completed and 9 pairs replayed. Shared-host L0 diagnostics, common-S unset | Apple M4 Pro, macOS 26.6, shared L0 | [study and run reports](../research/notes/ecc2k130/n53_s3_root_batch_20261012/RESULT.md) |
| 2026-10-12 | [`20261012-conductor-gap-counts`](20261012-conductor-gap-counts/) | native MITM-based IC and negation rho on two fixed CM pairs; operation counts only | IC floor/crater counted-S ratios: 0.975242 at ell = 8209 (99% 0.603463–1.450934), 1.127134 at ell = 131101 (99% 0.922471–1.341694); 192 verified records, 96 exact measured replays; no full-cost or wall advantage | Xeon Platinum 8573C, shared, L0 | [`study and reports`](../research/ecbench_conductor_gap_20261012/) |
| 2026-10-11 | [`20261011-groebner-vs-msolve-singular-polybori`](20261011-groebner-vs-msolve-singular-polybori/) | `f4_fp` and the inherited-F4 / F6-IC decomposition oracle against msolve 0.6.5, Singular 4.3.2 and PolyBoRi on identical exported systems (Magma scripts included, not run) | F_p: same reduced basis as Singular on 13/13 systems; msolve 1-thread 12–24× faster than ours on square generic systems (quad_n9: 246 ms vs 5 693 ms), level on overdetermined ones. Boolean n = 17, 8 targets: ours decides in 7–24 ms per target, msolve 42.7–46.5 s and PolyBoRi 40.4–52.3 s for the full basis, Singular > 120 s; verdicts agree | 4 vCPU Xeon 2.10 GHz, shared, L0 | [`research/gb_cas_comparison_20261011`](../research/gb_cas_comparison_20261011/) (PR #1675), run-20261011T1806Z |
| 2026-10-11 | [`20261011-inherited-f4-eager-rewrite`](20261011-inherited-f4-eager-rewrite/) | inherited F4 eager rewrite of kept rows against lazy rows (`KIC_F4_EAGER=0`), same binary, `f6_ic_probe` | lazy/eager median ratio: K₀/2²³ dim 8 inherited F4 1.40, F6-IC 1.62; K₀/2³¹ dim 9 1.48 and 1.53; counters and verdicts identical | 4 vCPU Xeon 2.10 GHz, shared, L0 | [`research/inherited_f4_eager_rewrite_20261011`](../research/inherited_f4_eager_rewrite_20261011/) (PR #1671), run-20261011T1657Z |
