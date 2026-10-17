#set page(paper: "a4", margin: 1.7cm)
#set text(font: "New Computer Modern", size: 9.5pt)
#set par(justify: true)

= Disjoint fiber price and Linux replay

*2026-10-16 · Native Rust · source-linked report*

The instrumented order-1,201 panel completed 60/60 cold children. A
separate Rust checker replayed 1,275 scalars and 6,126 relation witnesses;
the native audit found 60/60 deterministic matches with the frozen V2
answers, rows, group counts and earlier field totals. The public input,
protocol, raw children, audit and source snapshot are beside this report
under `research/isogeny_torsion_fibers_20261012/primitive_price_20261016/`.
The audit's first failed reducer attempt is retained there as V1; its
corrected result is `PRIMITIVE_AUDIT_V2.json`.

The generic two-class floor is $S = sqrt(pi / 4) = 0.886227$ in counted
group additions per $sqrt(r)$. Same-target signed rho and source BSGS
are the references, with rho A/A as a noise control. The graph displays
the top-level $F_(p^2)$ multiplication calls split by whether they
occurred inside a curve addition or doubling. It does not convert these
calls, square roots, allocations, or process work into one portable unit.

#align(center)[#image("PRIMITIVE_GRAPH_V2.svg", width: 82%)]

#table(
  columns: (auto, 1fr, auto, auto, auto),
  inset: 4pt,
  stroke: 0.4pt + gray,
  [Targets], [Arm], [$S_"group"$ mean], [Inside $F_(p^2)$ mul], [Outside $F_(p^2)$ mul],
  [1], [Full image], [88,499.672], [6,133,958], [6,122,892],
  [1], [Inverse Vélu], [134.457], [9,286], [29,140],
  [1], [Source BSGS], [2.972], [200], [6],
  [1], [Signed rho], [4.905], [336], [6],
  [64], [Full image], [98,067.547], [6,796,835], [6,123,144],
  [64], [Inverse Vélu], [9,702.332], [672,163], [29,392],
  [64], [Source BSGS], [84.171], [5,696], [132],
  [64], [Signed rho], [259.728], [17,853], [132],
)

At one target, the full-image/inverse counted-group ratio is 658.201
and inverse/rho is 27.410 using three-round means. These are group
count diagnostics. The group counts overlap the field call ledger and
must not be added to it. The inside/outside split places 29,140 Fp²
multiply calls in inverse setup and polynomial work versus six in
rho; the internal `sqrt_mod`, polynomial control flow, table work,
allocation, I/O and cross-host calibration are still open. Portable
complete $S$ remains unset.

== Linux independent-host receipts

GitHub Actions run 38961930762 replayed all 60 order-1,201 children,
1,275 scalars and 6,126 witnesses. Run 38961930725 replayed all 108
order-229 nine-arm rows, 996 witnesses and 2,484 trials. The artifacts
and SHA-256 manifests are copied into `LINUX_RECEIPTS_V1/`. Their
Ubuntu x86-64 AMD EPYC virtual host supplies independent L1 correctness
and descriptive performance; no physical-core L2 wall admission exists.

Run 38961930860 executed the hard-capped order-2,437 full-image probe.
It exited 134 after a 1,072-byte allocation failed, with peak RSS
534,470,656 bytes and no solver stdout. This is a resource-censored
row under the 512 MiB address-space cap. A memory-bounded constructor
is needed for the next same-cap attempt.

The next DLP experiment will freeze a projected destination factor
base without a raw preimage witness and compare its complete verified
solve with source rho, BSGS and relation controls. A changed base can
change relation support; its total cost must be measured directly.
