#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set table(stroke: 0.4pt)

= Inherited F4: eager rewrite of kept rows

2026-10-11. Change: `src/cryptanalysis/inherited_f4.rs` (PR branch `inherited-f4-eager-rewrite`). Control: `KIC_F4_EAGER=0`, same binary. Run report: `reports/run-20261011T1657Z/report.pdf`.

== Question

Where does a specialisation of an inherited Macaulay basis spend its time once the decomposition tree is fixed, and can it be cut without changing any row space, verdict or tree counter?

== Mechanism

Before the change, a specialisation $v := c$ kept the rows whose pivot avoids $v$ by reference, in their parent's layout, and rewrote a row into the child's layout only when a displaced row's re-reduction hit it as a pivot (lazy rows). A row last current $k$ layouts back needed the composed map of $k$ steps and a `pext`/`pdep` plan for that pair of layouts, built word by word on first use, plus a staging pass, a copy and one allocation per rewrite.

Callgrind on one n = 23, dim 8 target of the F6-IC arm (this host, the parent build, not committed as a record) counted 260 521 lazy rewrites at 622 specialisations, about 70 % of the rows of a basis at every node, with the lazy rewrite, the plan construction and the allocator together near 40 % of instructions. Since most rows are rewritten at every node anyway, the change rewrites every kept row at once, through one flat plan of the step's own map, into a single buffer shared by the node's rows (`StepPlan`, `ReducedBasis::rewrite_kept`); displaced rows use the same plan. Rows from an older layout, which occur only under a lazy ancestor, still take the lazy path.

The row space of every basis is unchanged: each row's image is the same vector, only computed earlier. The specialisation word count (`word_xors`) rises because every kept row's rewrite is now charged.

== Result

#figure(image("totals.svg", width: 100%), caption: [`f6_ic_probe` totals, run-20261011T1657Z, 3 interleaved rounds per mode (dots), bars at the median.])

#table(columns: 6,
  [cell], [arm], [lazy median ms], [eager median ms], [ratio lazy/eager], [counters],
  [K₀/2²³ dim 8, 4 targets], [inherited F4], [1477.0], [1058.4], [1.40], [identical],
  [K₀/2²³ dim 8, 4 targets], [F6-IC + index], [1034.8], [639.5], [1.62], [identical],
  [K₀/2³¹ dim 9, 1 target], [inherited F4], [1543.5], [1041.3], [1.48], [identical],
  [K₀/2³¹ dim 9, 1 target], [F6-IC + index], [1443.6], [946.1], [1.53], [identical],
)
Medians measured (3 rounds each); ratios derived as lazy median / eager median. In every round each eager value is below every lazy value of the same cell and arm. "Identical": reductions, splits, gate counters and every target's verdict.

== Scope and limits

- One host, isolation L0 (shared cloud container, no pinning), three rounds per mode; no bootstrap interval. The ratios are wall-time observations of the probe on this host, not an ecbench comparison.
- Solver-internal only: this changes the constant of the Gröbner track, not its exponent. The decomposition search still enumerates about $|F|^(m-2)$ partial assignments, and nothing here compares index calculus with rho.
- `cargo test --release --lib`: the `inherited_f4::` (both modes), `koblitz_groebner::` and `f6_ic` tests pass. The full suite's 12 failures are in modules that do not reach this code (curve registry parsing, prime BSGS, BIKE, an ecbench record test, n71 construction).

== Graphs checked and left unchanged

`docs/ic/progress-timeline.json`, `docs/index-calculus-scoreboard.html`, `docs/performance-gains/`: they plot IC-versus-rho accounting from ecbench sessions and the specialisation study of `performance-gains`; a probe-level solver speedup without an ecbench session changes none of their values or conclusions.
