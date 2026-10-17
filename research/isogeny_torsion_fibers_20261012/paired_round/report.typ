#set page(paper: "us-letter", margin: 0.58in)
#set text(size: 9pt)

= Explicit isogeny transport: same-target relation solves
2026-10-12 local · final run 2026-10-13 02:00 UTC · T-2223 / T-2240

The degree-three map on each registered prime-order pair has a geometric
kernel of size three and a rational kernel of size one. We tested the same
eight public source targets per pair with source relations, transported
relations, and the repository's negation rho. All 48 arm answers verified.
Full rational-point group-law checks passed 361 and 44,521 ordered pairs.
Final structural review found public degree-one scaling isomorphisms to
both destination models; an isomorphism-only solver arm is still owed.

#image("transport.svg", width: 100%)

The source and destination use the exact full models and EC1 generators
listed above. The map is φ(x,y) = ((x³+4b)/x², y(x³−8b)/x³). For the
standard construction, see Shumow, *Isogenies of Elliptic Curves: A
Computational Approach*, arXiv:0910.5370. Full model tuples and frozen
inputs are in `curves.json` and `PROTOCOL.md`.
The scaling ψ_u(x,y) = (u²x,u³y) uses u = 2 at p = 13 and u = 14 at
p = 199. The independent verifier checked both maps on 44,882 ordered
source-point pairs. The diagram marks this verified but unmeasured arm.

Counted group-operation S is (additions + doublings)/√r. The reference is
matched negation rho on the exact source Q. Its two-class walk-step floor
is √(π/4) ≈ 0.88623 S. Each curve also has a rational order-three
automorphism (ζ = 3 for F13, ζ = 92 for F199), so the six-class curve
floor is √(π/12) ≈ 0.51166 S. Six-class rho was not run. Ranges below
are observed min–max over eight targets, not confidence intervals.
Ratios are sums over matched targets.

#table(
  columns: (0.7fr, 1.8fr, 1.7fr, 0.9fr, 0.9fr),
  inset: 3pt,
  stroke: 0.3pt,
  [r], [arm], [mean S (range)], [S/rho], [S/A2 floor],
  [19], [negation rho], [10.725 (9.177–13.536)], [1.000], [12.102],
  [19], [source relations], [15.371 (12.159–18.583)], [1.433], [17.344],
  [19], [transport relations], [15.285 (13.536–17.665)], [1.425], [17.247],
  [211], [negation rho], [5.938 (4.475–7.297)], [1.000], [6.700],
  [211], [source relations], [50.186 (29.396–66.846)], [8.452], [56.629],
  [211], [transport relations], [41.332 (36.762–48.810)], [6.961], [46.638],
)

#pagebreak()

= Phase evidence and limits

#image("group_ops-v4.svg", width: 100%)

At order 211, mean counted group operations per target are 729.000 for
source relations, 600.375 for transport relations, and 86.250 for rho.
The source/transport ratio is 5,832/4,803 = 1.21424 in this counted
unit. Relation collection uses 525.25 versus 387.75 group operations;
the mean trial counts are 61.50 versus 44.875. Both reach rank 12.
The same destination curve is reachable through ψ_u, so this yield gap
has not been isolated as a degree-three isogeny effect.
Setup is 78 operations each; target descent is 23.25 versus 23.125;
verification is 102.50 versus 111.50. Rho uses 77.75 operations in
its walk and 8.50 in the runner check, on average.

The destination map evaluates G and Q per workload, costing 10 field
multiplications, 4 squarings and 2 inversions in this implementation.
Factor-base square roots, lookup work, modular row operations, and rho
distinguished-point table probes lack pinned conversion into S here.
The chart is therefore a counted-operation diagnostic, not an admitted
complete-cost comparison. No wall-time claim is made on this shared
Darwin arm64 host (L0).
For the six-class curve floor, multiply the table's last column by √3.
The order-211 values are 11.605, 98.085 and 80.779. The selected
negation rho is an exact-target control, not a certified best rho on these
curves. This correction changes no raw measurement.

Final `raw-v4.tsv` and `replay-v4.tsv` compare byte for byte, SHA-256
`995b7bcc33f5584d7b4fc6102b7ade2b3f19c7e5b39cb406b0b1c33dc7af11ad`.
This is same-host replay. Earlier v1–v3 runs are retained with their
phase-label or output-schema defects in `reports/`; only v4 supports this
table. Source, commands, full phase counters, run reports, and requirement
status are in `REPORT.md`. A separate native verifier with its own point
law passed all 48 records, checked the four group orders and all phase
sums, checked the cubic automorphism on 89,764 ordered point pairs,
and scaling isomorphisms on 44,882 pairs; its receipt is
`validation/independent_arithmetic_v4.log`. Lean checked concrete cubic
roots, sixth-power scalings and conditional equation identities in
`lean/CubicAutomorphism.lean`. General Lean group-action theorems remain open.
The controlled two-thread, loopback-accessible release library suite
passed 4,009 tests with none failed and 113 ignored. The earlier
restricted-sandbox run failed 13 and remains separately recorded.

Immediate mechanism check: solve the same 16 targets after ψ_u transport
with equal relation seeds, factor bases and caps; include the original
arms and verify all scalars. Then add charged six-class rho. A later native
paired-target ecbench measurement needs growing sizes, all cold costs,
independent-host replay, and a non-isomorphic pair to isolate the map's
specific effect. The gap is tracked in
`problems/nonrational-isogeny-fiber-dlp-cost-gap.md`.
