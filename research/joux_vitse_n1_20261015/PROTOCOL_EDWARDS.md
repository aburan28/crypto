# Protocol (frozen before data): the 2-torsion Edwards route end to end, with F4 trace replay

Registered 2026-10-15, before any claim run, at the commit that adds this file
(`src/cryptanalysis/joux_vitse_edwards.rs`, `examples/joux_vitse_edwards.rs`;
`jv_quintic_edwards.rs` changed only by making `polys_of` and `resolve`
`pub(crate)`).  Exploratory probes before it (3 paired residuals at
`p = 271`, 256 DLP trials at `p = 31`) sized the budget only.

## Question

`RESEARCH_K5_TORSION_JOUX_VITSE.md` measured one 2-torsion Edwards test at
`C″ = 7.3·10⁷` and left two things undone: the trace-driven elimination
(§5, priced at `3–10×`) and an end-to-end run (§6).  (1) What does F4 trace
replay buy on the test, paired on the same residuals?  (2) What does the
route cost end to end, complete DLPs against rho on the same target?

## Boundaries

- Reference: rho on the same subgroup of prime order `n` (`#E = 4n`), the
  same `r`-adding walk as the plain study, measured, plus its expected walk
  `√(πn/2)` at the measured `c_add` (452, the unified Edwards addition).
- Floor: a random residual decomposes with probability at most
  `8·C(|F|, 4)/n` (8 subgroup points per class-quadruple, the census of
  RESEARCH_K5 §3), so trials per relation are at least `n / (8 C(|F|, 4))`
  (`→ 192 p` at `|F| = p/4`).

## Fixed method

Instances `generate_instance_edwards(p, seed)`; base, symmetrised `S₅` and
pair table in the committed experiment's order (`seed ^ 0xED02`); F4 grevlex
with pairs bounded at degree 24 and no wall-clock budget; trace recorded on
the first inconsistent system (`seed ^ 0x7ACE`); replay, and the committed
solve after a mismatch or on a consistent system (both charged); every trial
compared with `mitm4_signed`, every decomposition re-added.  A relation
`aG + bQ = Σ s_i P_i + εT` enters `solve_logs` with cofactor 4.

## Runs

1. **Paired cost** (`cost`): `p ∈ {271, 521}`, seeds `{1, 2}`, 32 residuals
   each, every residual through the committed solve and through the replay.
   Output `data/edwards/cost_p<p>_s<seed>.jsonl`.
2. **Complete DLPs** (`dlp`): `p ∈ {31, 41, 61, 71}` seed 1, then
   `p ∈ {31, 41}` seed 2, run one after another, each capped at 10 hours of
   wall time.  Output `data/edwards/dlp_p<p>_s<seed>.jsonl`.

Host: the session container (4 vCPU, 15 GB), isolation none; counts are the
metric.

## Analysis plan

- Replay factor: per residual, committed-solve multiplications over replay
  multiplications (all charged), mean and range per `(p, seed)`.
- DLP: trials, relations, online multiplications (F4 over the batches, Weil
  restriction, root finding, sign resolution at `c_add`, target generation,
  linear algebra), IC / rho measured and against the expected walk; trials
  per relation against the floor.
- Model (extrapolation): `IC ≈ (|F| + 1) · κ_E p · C_trial` with `κ_E`
  fitted from the DLP arms' trials per relation times `|F|³/p³`-corrected
  rate, against rho `√(πn/4) c_add`; crossover, compared with the note's
  `2^{233}`.

## Success and abandonment

- Route found at runnable size if any arm reaches IC / rho `≤ 10`, or the
  model crossover falls below `2^{128}` in `n`.  Expected: neither; the
  replay `1.5–3×`, IC / rho `≈ 10⁶–10⁷`.
- Instrument fails if a replayed trial's decompositions differ from the
  committed solve's or the oracle's, a decomposition fails to re-add, or a
  logarithm fails `kG = Q`; reported, not dropped.
- An arm with `|F| < 5` is reported as degenerate and not replaced.
