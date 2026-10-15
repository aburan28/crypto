# Protocol: `k = 4` index calculus in characteristic 2 past one word, against rho

Registered 2026-10-15, before any claim run, after the pilots P1–P3 and
readings listed in `EXPERIMENTS.md` (W1–W12). Branch
`claude/k4-char2-crossover-20261015`, stacked on #1883. The claim binary is
built from the commit that registers this file; its sha256 is recorded in
`EXPERIMENTS.md` and in every session. Class: **engineering** for Parts D and
F; Part E is a **measured comparison against rho** on non-standard curves.

## Background

#1883 measured IC (Frobenius-reduced D₄ oracle) ÷ rho's theoretical cost
falling from 5,080 at `m = 6` to 6.06 at `m = 14`, slope −0.248 in
`log₂ n`, and extrapolated a crossover near `2^68`. Two things stood
between that and a measurement:

1. **The matrix.** The dense incremental elimination's `row_ops` grow as
   `|F|^3` (W1). From `m = 14`'s 1.06 % of the IC cost it reaches rho's
   cost near `m = 18`, so with that matrix IC never crosses.
2. **The word.** `#E(F_{q⁴}) ≈ 2^{4m}`; with cofactor 2 the subgroup order
   passes `2^64` at `m = 17`. The harness carried one-word orders.

This branch adds a sparse `wiedemann` matrix, an optional one/two
large-prime collector in `gaudry-d4-bin`, and `u128` versions of the
`rho.negation` walk and the IC walk pipeline (`ic_wide.rs`). At `S = u64`
both ports reproduce the one-word code exactly (tests, W7), and the
widened binary reproduces #1883's committed sessions (W8 re-audit).

## Part D: the sparse matrix (`SPEC-D-ladder.json`, `SPEC-D-m14.json`)

#1883's B2 workloads (`m = 8, 10, 12`, two targets) and B3's (`m = 14`, one
target). Arms: `rho-neg`, `ic-dense` (B2/B3's `ic-frobenius` exactly),
`ic-wied` (the same with `linalg=wiedemann`).

- **D1 (reproduction).** Every `ic-dense` run's `S` equals the committed
  B2/B3 value to the last digit.
- **D2 (LA).** `ic-wied`'s `row_ops` are below `ic-dense`'s at `m = 12`
  and `m = 14`, by at least 2× at `m = 14`.
- **D3 (relations).** `ic-wied` needs a full-rank core where the echelon
  matrix pins earlier: its trials are within `[1.00, 1.10]×` `ic-dense`'s
  on every target (pilot P2: `S` 1.01× at `m = 10`).

## Part E: the crossover (`SPEC-E-m15/16.json`, `SPEC-E-m17/18-{rho,ic}.json`)

Curves `f2m4_torsion(m, 1, 2)`: cofactor `h = 2` at every size, so
`r = 2^{59}, 2^{63}, 2^{67}, 2^{71}` for `m = 15, 16, 17, 18` (W11). `m = 15,
16` run the one-word code, `m = 17, 18` the `u128` code. Targets: two at
`m = 15, 16, 17`, one at `m = 18` (rho alone is about a day there).
Arms: `rho-neg` and `ic-wied` (Frobenius, sparsest, replay, `wiedemann`).
At `m = 17, 18` rho and IC run as two sessions on the same workloads, so
they can use two cores, and are compared by operations.

Predictions from pilot P3 (IC ≈ 1,500–1,700 GAE a trial, ≈ 186 trials a
relation, `q/4` relations):

- **E1.** Every run verifies; `ecbench verify --replay-all` reproduces every
  run.
- **E2.** IC `S` ÷ rho's theoretical `S` (0.886) lies in `[1.5, 8]` at
  `m = 15`, `[0.7, 4]` at `m = 16`, `[0.4, 2]` at `m = 17`, `[0.2, 1]` at
  `m = 18`.
- **E3 (the claim).** At `m = 18`, IC's `S` is below rho's theoretical
  `S` **and** below the measured `rho-neg` `S` on the same target. If the
  measured rho falls below IC there (rho's own spread is about ±40 % per
  target), E3 fails and the report says IC beats rho's expectation but not
  that run.
- **E4.** The `wiedemann` matrix is at most 5 % of IC's GAE at `m = 18`.
- **E5.** The measured crossover `m*` (IC ÷ rho theory = 1, interpolated
  in `log₂ n` between adjacent sizes) lies in `[16.5, 18]`, i.e. `n*` in
  `[2^{66}, 2^{72}]`.

## Part F: large primes (`SPEC-F-lp.json`)

B2/B3 curves at `m = 10, 12, 14`, one target each. Arms: `ic-wied` (no
large primes), `lp1-700` (one large prime, `F′` = 70 % of the columns),
`lp2-700`, `lp2-500` (two large primes). The pilot at `m = 10` (P2) gave
1.06–3.0× the plain cost.

- **F1.** Every large-prime arm costs more `S` than `ic-wied` at every
  size: the matrix they shrink is under 1 % of the cost here.
- **F2.** The best large-prime arm's cost ratio over `ic-wied` does not
  fall below 1.0 at `m = 14`.
- **Decision.** F1 and F2 holding means the double-large-prime variant is
  a recorded negative at reachable sizes. The model says when it would pay:
  when sparse LA (`≈ 22·(q/4)²` row operations) reaches the relation cost
  (`≈ 186·(q/4)·C₄` field multiplications), `q ≈ 2^{23}`, `n ≈ 2^{92}`.

## Part G: the independent runner

`ecbench isolab-job` files for E-m18 are generated and committed. This
container has no isolab credentials or MCP tools, so they are not
submitted; the report says so (AGENTS.md §9: no replacement keys).

## Decision and scope

- E3 holding is the first measured `S_IC < S_rho` on this family; it is a
  statement about non-standard curves over `F_{2^{4m}}`, which the theory
  already predicts to be weak, with the `F_q` multiplications priced at the
  curve's own multiplications per addition, `row_ops` at one each, at L1
  isolation (operation counts; wall time descriptive).
- A D1 miss stops every claim until the cause is found.
