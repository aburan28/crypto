#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.5pt)
#show raw: set text(size: 8.5pt)
= Open IC leads after 2026-10-11 — search-round report

Source: `research/ic_open_leads_scoping_20261011/REPORT.md` (aburan28/crypto). Each statement is labelled *measured*, *derived* (calculated from recorded values or formulas) or *scoped* (from the repo's own notes, not run).

#figure(image("leads.svg", width: 100%), caption: [Left: Macaulay columns at the deciding degree for the m = 4 solving-degree cells (derived, Σ C(N, k ≤ D)); the dashed line is the (9,3) degree-7 matrix the DREG note already classes as a large-host job. Right: online operations per target at equal advice on the n = 73 / 83 fixtures — IC target probes and plain-rho ideal steps recorded, BL online steps derived. Data: `reports/data.json`.])

#table(columns: 3, inset: 4pt,
  [*lead*], [*status*], [*evidence*],
  [1. m = 4 solving degree, ℓ = 4–6], [scoped: at ℓ ≤ 3 "ℓ + 4" and "levels at 7" predict the same D; the first separating m = 4 cell is ℓ = 4 (8 vs 7), ≈ 1.9·10⁷ columns at D = 7, 8.1·10⁷ at D = 8. Cheapest decisive cell in the programme: m = 3 (15,5) at D = 7 on a ≈ 128 GB host.], [derived],
  [2. equivariant F4 → full F4], [scoped: rank-only probe at m = 2 today; full F4 is 3–5 days and needs a Frobenius-invariant m = 4 formulation first; constant (≤ n) lever.], [scoped],
  [3A. preprocessing rho at equal advice], [IC online loses ≈ 48× (n = 73) and ≈ 3,400× (n = 83) in operations; the note's n = 73 figure used rank-stage probes; BL precompute is ≈ 10¹¹ steps, not "seconds".], [derived],
  [3B. free-relation census], [rank 49 / 1 / 6 against 599 needed at n = 61 / 73 / 83; note's count high by ≈ n; rank stage not free.], [measured],
  [4. D₄-symmetrised S₅ lattice], [1/80 is the Jochemsz–May bound, never measured; e-rewrite and (Z/2)³ invariants leave the determinant condition unchanged; 13× exponent gap; closed by counting.], [derived],
  [5. tree split round 2], [pre-registered, run; see its RESULTS.md and run report.], [measured],
)

*Recommendation.* The one remaining exponent-relevant experiment is a solving degree: run m = 3 (15,5) at degree 7 on a 128–256 GB host (pre-registered, existing `dreg_ladder`); only if it levels at 7, build the m = 4 barrel formulation plus equivariant blocks and test (13,4) closure at degree 7. Nothing in this round supports an exponent, speedup or deployed-curve claim.
