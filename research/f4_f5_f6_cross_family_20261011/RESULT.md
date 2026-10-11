# F4, F5 and F6-IC across curve families: result

Protocol: [`PROTOCOL.md`](PROTOCOL.md), committed before any measured run.
Analysis: [`ANALYSIS.json`](ANALYSIS.json), written by the native analyser
[`examples/f4f5f6_cross_family_analyze.rs`](../../examples/f4f5f6_cross_family_analyze.rs)
from the sessions' records only. Run report with graphs and PDF:
[`reports/run-20261011T1904Z/`](reports/run-20261011T1904Z/REPORT.pdf).

## Answer in one paragraph

On both binary families the four Boolean decomposers rank the same way,
and by similar margins:

- **Matrix F5 against matrix F4.** F5 does 0.64–0.65 of F4's counted
  solver work on the same splitting tree, on all 16 targets.
- **Matrix F4 against inherited F4.** Matrix F4 costs about 200 times
  inherited F4 at m = 13.
- **F6-IC against inherited F4.** F6-IC cuts inherited F4's word XORs by a
  median factor of 1.22–1.29 (log₂ 0.29–0.37), and its relation phase runs
  1.7–1.9 times faster.

But F6-IC buys that cut with exact point additions in the factor base.
Charged at one addition each, they raise its lower-bound `S` 65–1,350
times above inherited F4's. Every IC arm, in every cell, is slower than
the same-target rho on every target. Prime-field curves have no
algebraic oracle or F6 in the repository; the evidence that exists for
them is collected at the end and is not in this round's unit.

## What ran

| spec | session | curve | targets | arms | runs verified | replay audit |
| --- | --- | --- | --: | --- | --: | --- |
| `SPEC-K13` | `ECBS1h2b56c13eeaec` | `icv1-f2m13-t181-515ee569` (Koblitz `E_0`) | 8 | strong rho, rho-negation, matrix F4, matrix F5, inherited F4, F6-IC | 48 / 48 | [ok](audits/K13-replay-all.json) |
| `SPEC-R13` | `ECBS1h78657ab42aba` | `icv1-f2m13-tm9-6863a838` (random binary) | 8 | rho-negation, matrix F4, matrix F5, inherited F4, F6-IC | 40 / 40 | [ok](audits/R13-replay-all.json) |
| `SPEC-K19` | `ECBS1h4af82d2330b7` | `icv1-f2m19-t797-b6cf2467` (Koblitz `E_0`) | 8 | strong rho, rho-negation, inherited F4, F6-IC | 32 / 32 | [ok](audits/K19-replay-all.json) |
| `SPEC-R19` | `ECBS1hfff3e0643af0` | `icv1-f2m19-tm649-a9133fe6` (random binary) | 8 | rho-negation, inherited F4, F6-IC | 24 / 24 | [ok](audits/R19-replay-all.json) |
| `SPEC-X19` | `ECBS1h9453f8e9b41e` | both m = 19 curves | 1 each | matrix F4, matrix F5 | 0 / 4 (4 timeouts at the 1,800 s cap) | [ok](audits/X19-replay-all.json) |

- **Run window.** 2026-10-11, 19:04 UTC to 21:39 UTC. Each spec runs one
  round; every arm of a target shares one measurement seed.
- **Host.** An x86-64 Linux cloud VM. `ecbench` recorded the isolation it
  actually earned per run, a mix of L1 and L2 (counts per arm are in
  `ANALYSIS.json`). The protocol asked for L1. Word XORs and `S` are counted
  and deterministic, so isolation affects only wall time.
- **Binary.** One `ecbench` binary, sha256 `780e1084b13a5da7…`, built from
  commit `bdd3c466` and recorded in every record of every session.
- **Recorded git commits differ.** The sessions record the checkout's
  `HEAD` at their start: `bdd3c466` (K13), `bc2d92fc` (R13, K19) and
  `83c3a2c8` (R19, X19). Two local commits were made while the sessions
  ran. Between `bdd3c466` and `83c3a2c8` no file compiled into `ecbench`
  changed: the diff touches only the analyser example, the leaderboard
  script, the plot script and the scoreboard page. The binary hash is the
  identity that ties the records together, and it is the same throughout.

## Results

`W` is counted solver word XORs per one-target run (deterministic). `S` is
charged group-addition equivalents over `√r`, a **lower bound** for IC arms,
because word XORs are unpriced. Medians over 8 targets, with ranges.

### Counted solver work and cost per cell

| cell | base points | rho-negation S | strong rho S | matrix F4 W | matrix F5 W | inherited F4 W | F6-IC W | inherited F4 S | F6-IC S |
| --- | --: | --: | --: | --: | --: | --: | --: | --: | --: |
| Koblitz m = 13 (r = 2,003) | 33 | 2.290 [1.788, 3.419] | 14.318 [13.603, 14.341] | 2.245 × 10⁹ | 1.454 × 10⁹ | 1.112 × 10⁷ | 9.10 × 10⁶ | 16.285 [15.737, 17.215] | 1,161.7 [862.2, 1,544.6] |
| random m = 13 (r = 1,367) | 31 | 3.232 [2.461, 3.516] | — | 2.325 × 10⁹ | 1.495 × 10⁹ | 8.66 × 10⁶ | 7.99 × 10⁶ | 18.500 [15.931, 21.069] | 1,177.7 [519.2, 2,119.3] |
| Koblitz m = 19 (r = 130,873) | 139 | 1.310 [0.937, 1.728] | 3.090 [2.823, 3.354] | timeout | timeout | 6.70 × 10⁸ | 5.19 × 10⁸ | 5.106 [4.608, 5.219] | 6,483.4 [5,789.2, 7,072.0] |
| random m = 19 (r = 262,469) | 129 | 1.387 [0.724, 2.399] | — | timeout | timeout | 7.15 × 10⁸ | 5.62 × 10⁸ | 3.405 [2.918, 3.547] | 4,714.5 [3,369.7, 5,380.3] |

Matrix F4 and F5 `S` equal each other in every m = 13 cell (16.027 Koblitz,
18.162 random): `S` does not price word XORs, and the two run the same tree.
The m = 19 matrix entries are SPEC-X19's single target, not medians.

Relation-phase wall time per run, medians, a practicality note from mixed
L1/L2 runs:

| cell | matrix F4 | matrix F5 | inherited F4 | F6-IC |
| --- | --: | --: | --: | --: |
| Koblitz m = 13 | 30.39 s | 26.08 s | 0.226 s | 0.112 s |
| random m = 13 | 32.17 s | 25.36 s | 0.168 s | 0.096 s |
| Koblitz m = 19 | > 1,800 s | > 1,800 s | 6.12 s | 3.63 s |
| random m = 19 | > 1,800 s | > 1,800 s | 6.35 s | 3.59 s |

### The protocol's questions, with its decision rules applied

**Q1. F5 against F4 (matrix), m = 13.**

| family | W(F5)/W(F4), median [range] | same tree on all 8 | rule |
| --- | --- | --- | --- |
| Koblitz | 0.6465 [0.6428, 0.6516] | yes | **F5 cuts F4's work** |
| random | 0.6414 [0.6397, 0.6443] | yes | **F5 cuts F4's work** |

"Same tree" means equal relations found, targets tried and solver
reductions on every target. On these degree-3 PDP3 systems the F5
criterion removes about a third of the elimination work and changes
nothing else. Wall time falls less, to 0.79–0.86 of F4's.

**Q2. Matrix F4 against inherited F4, m = 13.** `W(matrix F4)/W(inherited)`:
198.6 [174.5, 212.8] on Koblitz and 222.1 [173.4, 544.0] on the random
curve. Reported as is (no rule).

**Q3. F6-IC against inherited F4.** Statistic `log₂(W inherited / W F6-IC)`.

| cell | median [range] | rule | derived: per solver call | targets with equal call counts | relation wall, inherited / F6 | S, F6 / inherited |
| --- | --- | --- | --- | --: | --- | --: |
| Koblitz m = 13 | 0.291 [0.285, 0.307] | **F6-IC saves solver work** | 0.287 [0.270, 0.299] | 7 / 8 | 1.86 [1.63, 2.32] | 68.8 |
| random m = 13 | 0.302 [−0.050, 0.334] | **mixed** | 0.306 [0.291, 0.353] | 5 / 8 | 1.78 [1.58, 2.15] | 64.6 |
| Koblitz m = 19 | 0.359 [0.352, 0.395] | **F6-IC saves solver work** | 0.357 [0.332, 0.360] | 7 / 8 | 1.73 [1.52, 1.85] | 1,303.2 |
| random m = 19 | 0.372 [0.319, 0.376] | **F6-IC saves solver work** | 0.373 [0.367, 0.382] | 6 / 8 | 1.84 [1.74, 2.00] | 1,352.6 |

- **The random m = 13 cell is "mixed" by the rule as written.** One target
  (index 5) has log₂ −0.050. On that target F6-IC's run made 53 solver
  calls to inherited F4's 49: the two engines returned different
  decompositions, so their relation matrices filled at different times.
- **Per solver call** (a derived figure, outside the preregistered rule)
  F6-IC does less work on every target of every cell. The per-call median
  is within 0.01 of the whole-run median in every cell.
- **Why trees diverge.** No run exhausted a budget, and every arm verified.
  When a target has more than one decomposition, each engine returns the
  one its search order finds first, so the relation rows differ. This
  happens on Koblitz too (K13 target 0: F6-IC 39 trials, inherited 40).
  The protocol defines an arm as a whole solver choice, so the whole-run
  ratio is the preregistered statistic; only matrix F4 and F5 are on one
  tree by construction.
- **The 256-point closure cap is not reached.** Every base here has at most
  139 points. So unlike the E_0 ladder at m = 23, F6-IC's closure ran in
  every cell.

**Q4. IC against rho.** Every IC arm's lower-bound `S` is **above** the
same-target rho on all 8 targets of every cell, against both references:
**IC slower than rho on every target**. Median `S` ratios:

| cell | inherited F4 / rho-negation | inherited F4 / strong rho | matrix F4 / rho-negation | F6-IC / rho-negation | F6-IC / strong rho |
| --- | --: | --: | --: | --: | --: |
| Koblitz m = 13 | 7.21 | 1.136 | 6.97 | 509.8 | 81.1 |
| random m = 13 | 5.89 | — | 5.69 | 399.1 | — |
| Koblitz m = 19 | 3.77 | 1.581 | — | 4,883 | 2,069 |
| random m = 19 | 2.28 | — | — | 2,774 | — |

- **The two rho references disagree at these sizes.** Strong rho (32 lanes,
  `dp_bits` 4, the F6-IC ladder's reference) costs more than plain
  rho-negation on every Koblitz target: 4.2–8.0 times at m = 13 and
  1.7–3.6 times at m = 19.
- **Why.** Its setup phase dominates at `r` ≈ 2^11–2^17. At K19 one
  record charges 816 GAE to setup against 205 to the search.
- **Consequence.** Strong rho is the weaker comparison here; rho-negation
  is the stronger one, and both give the same verdict.

**Q5. Matrix arms at m = 19 (SPEC-X19).** All four runs (matrix F4 and F5, each m = 19 curve, one target) reached the 1,800 s cap and are recorded as timeouts, so they enter no ratio. Inherited F4 finished the same targets in about 6 s of relation phase.

## Koblitz against generic binary

- **The same ordering and similar margins on both families**, at both
  sizes. F5/F4 0.646 vs 0.641. Matrix/inherited 199 vs 222. F6-IC's
  per-call saving 0.29 vs 0.31 at m = 13, 0.36 vs 0.37 at m = 19.
- **Frobenius does not reach the decomposers.** The Koblitz runs fold
  the base by negation only (33 points, 17 columns at m = 13), as the
  generic ones do. So the IC arms see no Koblitz advantage, while the
  floor and the strong-rho reference do.
- **The floors differ by √m:**
  - `√(π/4m)` = 0.246 (m = 13) and 0.203 (m = 19) for Koblitz;
  - `√(π/4)` = 0.886 for generic.

  No IC arm is within an order of magnitude of either.
- **IC `S` falls from m = 13 to m = 19 on both families.** Inherited F4
  goes 16.3 → 5.1 on Koblitz and 18.5 → 3.4 on random. Rho-negation `S`
  falls too, 2.29 → 1.31 and 3.23 → 1.39, because small-`r` runs carry a
  fixed overhead. Two sizes give no exponent, and none is claimed.

## What this does not show

- **No speedup in `S`.** F6-IC's word-XOR saving is a stage diagnostic.
  In `S` it is a relabelling: the solver's work moves into charged point
  additions. Whether those additions are over-priced at one GAE each (the
  wall times suggest so: F6-IC's relation phase is 1.7–1.9 times faster)
  is an accounting question for a separate round.
- **No exponent.** Two sizes per family.
- **No claim about F6-IC beyond 256-point bases**, where #1333's closure
  is capped (the F6-IC ladder covers that).
- **Nothing about prime curves** (below).

## Prime curves: existing evidence, not comparable

Not measured in this round. `ecbench` has no algebraic decomposition
oracle for curves over `F_p`, and there is no F6-IC for them: F6-IC's gate
encodes the factor base as an `F_2`-linear subspace of abscissae, and `F_p`
has no proper additive subspace to play that role. What the repository has
measured elsewhere, each in its own unit:

| study | setting | engines | what was measured | unit |
| --- | --- | --- | --- | --- |
| [PKM tower, round 4](../notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md) (§13) | Kummer tower systems over `F_p`, `m` = 2–4, `N` = 8–20, 18 systems | `f4_fp_tower` (F4) against `sig_fp_tower` (signature F4, F5/GVW criteria) | signature/F4 multiply-adds 0.28–0.37 at `m` = 2, 0.20 and 0.58 at `m` = 3, 1.49 and 1.55 at `m` = 4; peak memory 1.7–28× F4's; solving degree 1–2 above F4's on all 18; `m` = 3, `N` = 15 ran out of the 14 GB cap | `F_p` multiply-adds, peak MB |
| [PKM tower gate](../notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md) (§14) | `m` = 2–4 systems at `N` ≤ 16 | signature criteria with F4's degree steps | multiply-adds 0.25–0.69 of F4's under the signature-degree-first order, memory 5–23×; the gate (`D_sig = D` at `m` = 3, `N` = 9–12) failed under both module orders | `F_p` multiply-adds |
| [Residual walks, Gaudry setting](../notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md) (§11) | `E(F_{p³})`, subspace base `x ∈ F_p`, `n` = 2^24–2^33 | resultant pair test (no Gröbner engine) | whole-pipeline `S/rho` 1,057–2,240×, total cost `∝ n^0.69` | `S`, the same unit as `ecbench` |
| [GLV orbit systems](../notes/index-calculus/RESEARCH_GLV_INDEX_CALCULUS.md) | `j = 0` prime-field curves, orbit-invariant summation systems | Macaulay and degree-by-degree F4 | orbit system 12–30× the ordinary system's field multiplications (merged engine); F4 block-aware 0.945–0.954 of plain | field multiplications |

Reading across: on prime-characteristic towers the signature engine
(F5's criteria) removes its zero reductions but pays in solving degree
and memory, and F4 stays the measuring engine. That matches the binary
cells here only in direction (F5 trims F4's arithmetic); the binary F5
here runs on the same tree as F4 and never raises the degree, because
the Boolean PDP3 systems are tiny and capped at degree 3. None of these
prime figures is in this round's unit or on its curves, so no ratio is
formed between them and the binary cells.


## Files

- `PROTOCOL.md`, `SPEC-*.json`: the preregistration and frozen specs.
- `sessions/{K13,R13,K19,R19,X19}/`: raw `ecbench` sessions (records,
  plan, host, executions).
- `audits/*-replay-all.json`: `ecbench verify --replay-all` receipts.
- `ANALYSIS.json`: the native analyser's output.
- `reports/run-20261011T1904Z/`: run report, graphs (`plot.py` draws them
  from `ANALYSIS.json`), Typst source and PDF.
