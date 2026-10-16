#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set table(stroke: 0.4pt)

= Gröbner engines against msolve, Singular, PolyBoRi (and Magma, pending)

2026-10-11. Run report: `reports/run-20261011T1806Z/report.pdf` (all values, raw records, inputs).

== Question

How fast are this repository's Gröbner engines next to the reference systems? The question was asked about Magma's F4 and F5. Magma is proprietary and is not installed in this container, so it was not run. Every system is exported as a Magma script (below), and its timings can be added from any machine with a licence. msolve, the open-source F4 from Faugère's group, stands in for it on $FF_p$. PolyBoRi stands in for Magma's Boolean-ring F4.

== Results

#figure(image("fp.svg", width: 100%), caption: [Reduced grevlex Gröbner bases over $FF_p$, run-20261011T1806Z (median of 3 rounds).])

*$FF_p$ (`f4_fp`).* Ours returns exactly the reduced Gröbner basis Singular computes on all 13 systems, in all 3 rounds. Ratios below are derived as our median ÷ the tool's median (single thread on both sides unless stated).

#table(columns: 5,
  [case], [ours ÷ msolve (1 thr)], [ours ÷ msolve (4 thr each)], [ours ÷ Singular], [kind],
  [quad_n8_p65521], [12.2], [12.1], [1.23], [square, generic],
  [quad_n9_p65521], [23.1], [24.0], [1.39], [square, generic],
  [cubic_n5_p31], [13.4], [10.7], [4.60], [square, generic],
  [overdet_n10_m15_p65521], [0.96], [0.75], [0.15], [overdetermined, inconsistent],
  [overdet_n12_m24_p65521], [1.04], [1.01], [0.19], [overdetermined, inconsistent],
)
- On square generic systems msolve is 12–24× faster than `f4_fp`, and the gap grows with size.
- On overdetermined systems, where $1$ appears at a low degree, ours matches msolve and is 5–7× faster than Singular's `std`.
- Singular is 1.2–1.4× faster than ours on the quadratic systems and 4.6× faster on the cubic one.

#figure(image("bool_n17.svg", width: 100%), caption: [Deciding the $m = 3$ decomposition systems of $K_1 slash FF_(2^17)$ (35 variables, 34 cubics), run-20261011T1806Z.])

*Boolean decomposition systems (the index-calculus oracle).* Ours decides each of the 8 targets in 11–24 ms (inherited F4) or 7–16 ms (F6-IC). msolve takes 42.7–46.5 s and PolyBoRi 40.4–52.3 s for the full reduced basis, and Singular does not finish within 120 s. All verdicts agree; target 3 is the one that decomposes. Derived per-target ratios (msolve ÷ ours): 1 770–4 060 for inherited F4 and 2 620–6 320 for F6-IC.

These are different computations. Ours splits variables and keeps a degree-bounded Macaulay basis per node, and it only has to decide whether the system has a root. The others compute the whole basis of the ideal. The gap is what the tree and the degree bound buy on these systems, not faster linear algebra.

== What this says about the engines

- `f4_fp` is correct (bases verified identical) but is not at reference speed on square generic systems. The msolve gap of 12–24× is the measured distance to state-of-the-art $FF_p$ F4. The audit's profile (§4.5–4.6 of `docs/ic/GROEBNER_F4F5F6_AUDIT_20261010.md`) puts the remaining cost in symbolic preprocessing and the sparse reduction of pair rows, which is where msolve's hashed monomial tables and sparse-dense linear algebra differ.
- On the systems index calculus actually poses, decided rather than fully solved, the inherited-F4 tree beats full Gröbner bases by three orders of magnitude. The comparison that matters for a claim is still the full ECDLP against rho (ecbench), which this does not touch: the decomposition search remains $approx |F|^(m-2)$ such decisions.

== Running Magma

On a machine with Magma, from `reports/run-20261011T1806Z/`:

```
for f in inputs-fp/*.magma inputs-bool-n17/*.magma; do echo "$f"; magma -b "$f"; done
```

Each $FF_p$ script times `GroebnerBasis` and prints `magma_gb_cpu_s`, `magma_gb_wall_s`, and whether our basis generates the same ideal and equals Magma's reduced basis. Each Boolean script works in `BooleanPolynomialRing(n, "grevlex")` and prints the time and `unit_ideal`.
- The scripts run Magma's default algorithm, which is F4 for these rings.
- Whether a given Magma version exposes an F5 variant, and under which parameter, should be checked in its handbook (`?GroebnerBasis`). The scripts do not assume one.
- Report the Magma version and host next to the numbers. Timings from another machine compare with the table above only as a ratio against msolve run on that same machine.

== Scope and limits

One shared host (isolation L0), three rounds for $FF_p$, one external run per Boolean root, n = 17 only. msolve's wall time includes process start-up (about 6–7 ms), which dominates the small cases. These are reference-tool comparisons of engines, not an index-calculus-versus-rho result, so no canonical graph changes. Checked and left unchanged: `docs/ic/progress-timeline.json`, `docs/index-calculus-scoreboard.html`, `docs/performance-gains/`.
