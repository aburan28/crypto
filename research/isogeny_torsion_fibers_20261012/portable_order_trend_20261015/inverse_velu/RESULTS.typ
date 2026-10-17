#set page(paper: "a4", margin: 16mm)
#set text(size: 10pt)
= Order-1,201 inverse-Vélu factor-base result
*Frozen same-target native comparison; three rounds, five arms.*

The degree-three inverse x-coordinate cubic over `F_(1237²)`
constructed the same eight image-base columns as the full source
census after 29 destination-point checks. Every one of 255 paired
inverse/full target relation rows matched. The independent checker
replayed 1,275 scalar outputs and 6,126 relation witnesses across
60 successful cold child processes.

#image("PIPELINE.svg", width: 100%)

The single counted unit in this table is
`S_group=(additions+doublings)/(N sqrt(1201))`. It includes the
complete child through final scalar replay, but is a lower bound on
fully priced portable `S`.

#table(
  columns: (2.8fr, 1fr, 1fr),
  inset: 5pt,
  [*arm*], [*cold N=1*], [*reuse N=64*],
  [inverse-Vélu plus relations], [149.645], [151.595],
  [full-image census plus relations], [88,514.859], [1,532.302],
  [source BSGS], [2.972], [1.315],
  [source signed rho], [5.107], [4.071],
  [rho A/A], [5.107], [4.071],
)

At N=64, paired group full/inverse was 10.108, inverse/rho 37.235,
and inverse/BSGS 115.266. The corresponding paired field-multiply
and field-inversion inverse/rho ratios were 41.699 and 37.694.
Whole-child CPU ratios on the shared Darwin host were 43.664,
9.159 and 11.380 for full/inverse, inverse/rho and inverse/BSGS;
rho A/A spanned 0.823–1.054. Those CPU values are L0 descriptive.

#image("ANALYSIS_V1.svg", width: 100%)

The old construction census used 3,063,003 counted group calls in
round zero; direct inverse construction used 667. Its relation yield
remained identical, so the algorithmic next test is a separate
nonfactorable `chi` decomposition polynomial oracle with exhaustive
pair-enumeration control and a complete target solve. Disjoint native
pricing, L2 timing and a broader subgroup-order trend remain open.
`RESULTS.md` links the frozen spec, raw run, independent replay,
Lean identity scope and follow-on protocol.
