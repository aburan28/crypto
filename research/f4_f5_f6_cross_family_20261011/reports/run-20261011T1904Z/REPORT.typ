#set page(paper: "a4", margin: 1.8cm, numbering: "1")
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show table: set text(size: 8pt)

#align(center)[
  #text(size: 15pt, weight: "bold")[F4, F5 and F6-IC on Koblitz and generic binary curves] \
  #text(size: 10pt)[Run report, `ecbench`, 2026-10-11 from 19:04 UTC (run `run-20261011T1904Z`)]
]

= What ran

- *Protocol.* `research/f4_f5_f6_cross_family_20261011/PROTOCOL.md`,
  committed before any measured run (commit `bdd3c466`).
- *Binary.* One `ecbench` binary, sha256 `780e1084b13a5da7…`, built from
  `bdd3c466`. The sessions record checkout `HEAD` at their start: `bdd3c466`
  (K13), `bc2d92fc` (R13, K19), `83c3a2c8` (R19, X19). No file compiled into
  `ecbench` changed between those commits, and every record carries the
  same binary hash.
- *Host.* An x86-64 Linux cloud VM. The isolation earned per run is a mix of
  L1 and L2, recorded per arm in `ANALYSIS.json`.
- *Arms.* Four Boolean decomposers through `pdp3-binary:m=3,degree=3,node_budget=8192`
  on `binary-standard-subspace:dimension=⌈m/3⌉` (5 at m = 13, 7 at m = 19),
  `max_trials=200000`:
  - matrix F4 and matrix F5, which split on the lowest free variable;
  - inherited F4 and F6-IC, which split on the highest free variable.

  References, run on the same targets:
  - `rho.negation`;
  - on Koblitz curves only, `rho.signed_frobenius_strong` with 32 lanes,
    `dp_bits` 4.
- *Seeds and workloads.* Target seed 202610111901, measurement seed
  202610111902, one round. The workloads are 8 public `hash_to_subgroup_v1`
  targets per curve. SPEC-X19 uses the first target of each m = 19 curve.

#table(
  columns: 7,
  align: (left, left, left, right, left, right, left),
  table.header([spec], [session], [curve], [targets], [arms], [verified], [replay]),
  [K13], [`ECBS1h2b56c13eeaec`], [`icv1-f2m13-t181-515ee569` Koblitz, r = 2,003], [8], [6], [48/48], [ok],
  [R13], [`ECBS1h78657ab42aba`], [`icv1-f2m13-tm9-6863a838` random, r = 1,367], [8], [5], [40/40], [ok],
  [K19], [`ECBS1h4af82d2330b7`], [`icv1-f2m19-t797-b6cf2467` Koblitz, r = 130,873], [8], [4], [32/32], [ok],
  [R19], [`ECBS1hfff3e0643af0`], [`icv1-f2m19-tm649-a9133fe6` random, r = 262,469], [8], [3], [24/24], [ok],
  [X19], [`ECBS1h9453f8e9b41e`], [both m = 19 curves], [1 each], [2], [0 / 4 (4 timeouts at the 1,800 s cap)], [ok],
)

Raw records: `sessions/<spec>/records.jsonl`. Replay receipts:
`audits/<spec>-replay-all.json`. Derived values: `ANALYSIS.json`, written by
`examples/f4f5f6_cross_family_analyze.rs` from the records alone.

= Measured values

*Units.*
- `W` is counted solver word XORs per one-target run. It is measured and
  deterministic.
- `S` is charged group-addition equivalents over `√r`. It is measured, and
  for IC arms it is a lower bound, because word XORs are unpriced and
  F6-IC's geometric point additions are charged.
- Medians over 8 targets, with ranges, are derived by the analyser.

#table(
  columns: 9,
  align: (left,) + (right,) * 8,
  table.header([cell], [rho-neg S], [strong rho S], [matrix F4 W], [matrix F5 W], [inherited W], [F6-IC W], [inherited S], [F6-IC S]),
  [Koblitz 13], [2.290], [14.318], [2.245e9], [1.454e9], [1.112e7], [9.10e6], [16.285], [1,161.7],
  [random 13], [3.232], [—], [2.325e9], [1.495e9], [8.66e6], [7.99e6], [18.500], [1,177.7],
  [Koblitz 19], [1.310], [3.090], [timeout], [timeout], [6.70e8], [5.19e8], [5.106], [6,483.4],
  [random 19], [1.387], [—], [timeout], [timeout], [7.15e8], [5.62e8], [3.405], [4,714.5],
)

#figure(image("solver_word_xors.svg", width: 100%),
  caption: [Counted solver word XORs per run, by engine and cell. Bars show
  the median and the black ticks the range over 8 targets (sessions K13,
  R13, K19, R19). At m = 19 the matrix arms (SPEC-X19, one target) reached the 1,800 s cap and have no bar.])

= The protocol's questions (derived; rules applied as written)

- *Q1. F5/F4 word XORs, matrix, m = 13.*
  - Koblitz: 0.6465 [0.6428, 0.6516].
  - Random: 0.6414 [0.6397, 0.6443].
  - The tree is the same on all 8 targets in both cells (equal relations,
    trials and reductions).
  - Rule: *F5 cuts F4's work*, in both families.
- *Q2. Matrix F4 / inherited F4 word XORs, m = 13.* 198.6 [174.5, 212.8] on
  Koblitz, 222.1 [173.4, 544.0] on random. No rule.
- *Q3. log₂(W inherited / W F6-IC).*

  #table(columns: 3, align: (left, right, left),
    table.header([cell], [median [range]], [rule]),
    [Koblitz 13], [0.291 [0.285, 0.307]], [saves],
    [random 13], [0.302 [−0.050, 0.334]], [*mixed*],
    [Koblitz 19], [0.359 [0.352, 0.395]], [saves],
    [random 19], [0.372 [0.319, 0.376]], [saves],
  )

  - The mixed cell's one negative target ran 53 solver calls against 49:
    the engines returned different decompositions.
  - The per-call ratio (derived, outside the rule) is above 0 on every
    target of every cell, with medians 0.287, 0.306, 0.357 and 0.373.
- *Q4. IC S against same-target rho.* Every IC arm is above both rho
  references on all 8 targets of every cell. Rule: *IC slower than rho on
  every target*.
- *Q5 (X19).* All four runs (matrix F4 and F5, each m = 19 curve, one target) reached the 1,800 s cap and are recorded as timeouts, so they enter no ratio. Inherited F4 finished the same targets in about 6 s of relation phase.

#figure(image("f6_vs_inherited_log2.svg", width: 80%),
  caption: [Per-target log₂(W inherited F4 / W F6-IC) for the four cells
  (sessions K13, R13, K19, R19; n = 8 each). The dashed line is the
  preregistered 0.1 threshold, and the grey band is the "no saving" band.])

#figure(image("s_per_arm.svg", width: 100%),
  caption: [Lower-bound S per arm and cell, median and range over 8 targets.
  The red dotted line is each family's floor: √(π/4m) for Koblitz, √(π/4)
  for generic binary. Strong rho exceeds rho-negation on every Koblitz
  target, because its 32-lane setup dominates at these group sizes.])

= Observations and limits

- *Relation-phase wall time* (medians; mixed L1/L2, a practicality note
  only):

  #table(columns: 5, align: (left,) + (right,) * 4,
    table.header([cell], [matrix F4], [matrix F5], [inherited F4], [F6-IC]),
    [Koblitz 13], [30.39 s], [26.08 s], [0.226 s], [0.112 s],
    [random 13], [32.17 s], [25.36 s], [0.168 s], [0.096 s],
    [Koblitz 19], [> 1,800 s], [> 1,800 s], [6.12 s], [3.63 s],
    [random 19], [> 1,800 s], [> 1,800 s], [6.35 s], [3.59 s],
  )

- *No run failed or timed out* in K13, R13, K19 or R19. X19's outcomes are
  listed above.
- *F6-IC's closure ran in every cell.* Every factor base has at most 139
  points, below PR 1333's 256-point closure cap.
- *Not shown:* any speedup in `S`, any exponent (two sizes), and anything
  about prime curves. Interpretation is in the study's `RESULT.md`.
