#set page(paper: "a4", margin: 16mm)
#set text(size: 9.4pt)
#set par(leading: 0.62em)

= Exact pairing transfer and complete cost

*Frozen extension-field ECDLP targets · V2 run 2026-10-15 05:55 UTC*

The public targets lie in prime-order subgroups of $E(F_(p^2))$,
$E:y^2=x^3+1$, with $(p,r)=(61,19)$ and $(271,61)$. The native solver
constructed and certified $F_(p^2)[z]/(z^9-2)$ or $F_(p^2)[z]/(z^5-2)$,
selected a nondegenerate $r$-torsion pairing input, computed the reduced
Tate values, and solved the multiplicative DLP without access to the known
scalar. It returned 12 and 10. Six separate source-arithmetic processes
verified $[k]G=Q$; all 18 public-multiple pairing checks and 12 rho controls
also passed. The answer-free input, selected H, source and binary hashes,
commands, raw records and SHA-256 manifest are committed beside this report.

The one-target cold reference is same-target negation rho. The table uses
$S_("CPU/GAE")=$ current-thread CPU nanoseconds divided by the paired
rho source-add calibration and $sqrt(r)$. The formal two-class counted
floor $sqrt(pi/4)=0.886226925$ is a marker in this local CPU-equivalent
unit, not a hardware-time bound. Each value is a three-round median;
brackets show the full three-round range. Ratios to rho are medians of
same-round full CPU ratios. Every answer in the table was verified.

#table(
  columns: (auto, auto, auto, auto, auto),
  table.header([$(p,r)$], [Variant], [$S_("CPU/GAE")$ median [range]], [Ratio / floor], [Paired / rho]),
  [$(61,19)$], [Rho], [64.01 [39.01, 98.06]], [72.23], [1.00],
  [], [Rho A/A], [50.12 [34.61, 71.05]], [56.56], [1.11],
  [], [Pairing + scalar replay], [163,807.34 [161,020.04, 170,683.25]], [184,836.79], [2,666.48],
  [], [Pairing + all checks], [170,737.19 [169,493.47, 178,065.56]], [192,656.29], [2,781.81],
  [$(271,61)$], [Rho], [32.76 [28.97, 33.27]], [36.97], [1.00],
  [], [Rho A/A], [40.90 [28.42, 74.12]], [46.16], [1.23],
  [], [Pairing + scalar replay], [21,978.75 [21,607.04, 28,981.92]], [24,800.37], [660.71],
  [], [Pairing + all checks], [23,193.81 [22,872.96, 30,304.74]], [26,171.41], [698.15],
)

#pagebreak()

== Paired CPU ratio

#image("GRAPH_V2.svg", width: 100%)

The graph uses `GRAPH_V2.typ` and `SUMMARY_V2.json`. The fully validated
row adds three post-solve pairing checks to the complete DLP method.
The separate pairing-solve row charges cold construction, unsuccessful
search, both pairings, the finite-field solve, and independent source
scalar replay. Setup reuse is one target. This is not an amortized result.

#pagebreak()

== Charged phase and primitive ledger

The median current-thread phase times in microseconds are:

#table(
  columns: (auto, auto, auto, auto, auto, auto),
  table.header([$(p,r)$], [H search], [Miller], [Final exp.], [Field DLP], [Full total]),
  [$(61,19)$], [101,723.959], [969.000], [759.625], [48.917], [108,490.417],
  [$(271,61)$], [23,873.166], [209.833], [143.333], [19.334], [25,681.333],
)

Candidate search used 26,469 of 27,169 selection-phase extension
multiplications at $p=61$ and 39,115 of 39,673 at $p=271$. The H-search
charge includes the group-order recurrence and cofactor calculation, as
well as failed x candidates and projections. Whole-workflow totals
including extra checks were 16,491,269 or 5,602,853 base-field
multiplications and 46 or 64 base-field inversions. Pairing peak process
RSS medians were 6,586,368 and 6,602,752 bytes. Every V2 phase subtotal
reconciled with its aggregate; V1 retained the same scalars and deterministic
primitive totals. V1 raw files and its own report remain under `runs/` and
`reports/run-20261015T0546Z/`.

The measured pairing route solved both frozen public targets, but its cold
cost exceeded same-target rho in both cells. H search dominates these
two runs. A reusable-H/alpha experiment needs an actual target count and
must charge setup once; a warm estimate is not a single-target result.
The still-missing nonfactorable chi polynomial decomposition oracle and a
different base-field non-isomorphic source/destination pair remain separate
unperformed tests. The present degree-three image also has a public scaling
isomorphism over these quadratic fields, so this run cannot attribute an
isogeny-specific relation effect.

The run used shared-host L0 current-thread CPU. A/A rho varied substantially.
Input parsing, output serialization and calibration were excluded from both
timed arms. These three-round ranges are descriptive, not confidence
intervals or an isolated wall-speed claim. The full source-linked report,
run records, exclusions and requirement statuses are in `RESULTS.md`.
