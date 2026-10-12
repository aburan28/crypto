# Gröbner engines against msolve, Singular, PolyBoRi and Magma (2026-10-11)

Results report: [`REPORT.pdf`](REPORT.pdf) (source `REPORT.typ`, graphs `fp.svg`, `bool_n17.svg`).
Run report: [`reports/run-20261011T1806Z/report.pdf`](reports/run-20261011T1806Z/report.pdf),
with raw records, the exported inputs, and the graph source.

- **F_p (`f4_fp`).** It computes the same reduced Gröbner basis as Singular on all 13
  `f4_fp_bench` systems. msolve (single thread) is 12–24× faster on the square generic
  systems and level with ours on the overdetermined ones. Singular `std` is 1.2–4.6×
  faster than ours on the square systems and 5–7× slower on the overdetermined ones.
- **Boolean decomposition systems (n = 17, 8 targets).** Inherited F4 decides each target
  in 11–24 ms and F6-IC in 7–16 ms. A full reduced basis takes msolve 42.7–46.5 s and
  PolyBoRi 40.4–52.3 s; Singular does not finish in 120 s. All verdicts agree. These are
  different computations (a decision tree with degree-bounded bases against a full basis).
- **Magma.** Not run: it is not installed in this container. The `.magma` scripts beside
  every input reproduce each measurement and compare bases; see "Running Magma" in the
  report. They were written for Magma's documented syntax but have not been executed here.

Reproduce: build `gb_cas_export` and `f6_ic_probe` in release, then
`./run.sh <target/release/examples> reports/run-<UTC> 3 600 120`
(set `POLYBORI_PY` to a Python with `passagemath-brial` to include PolyBoRi;
msolve and Singular come from the Ubuntu packages).
