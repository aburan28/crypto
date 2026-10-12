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
