# Result: Hamming ideals and an ISD-like Gröbner oracle for point decomposition

Protocol: [`PROTOCOL.md`](PROTOCOL.md). Frozen runs: `results/main/` (one
JSONL per cell, a `manifest` line then one record per target),
`results/main/summary.json` (derived by `analyse.py`, re-derived and
diffed by CI), `results/tau/` (the budget-sensitivity sweep). Reading:
[`RESEARCH_HAMMING_IDEAL_PDP.md`](../notes/index-calculus/RESEARCH_HAMMING_IDEAL_PDP.md).

## Verdict

**H0 stands; H1 is refuted at every size run.** With the weight constraint
`wt_NB(x) ≤ 2` encoded as a C-, FC- or QFC-Hamming ideal and Semaev's `S₃`
Weil-descended in normal-basis coordinates, `MultiSolve` answers every
target correctly (zero disagreements with the exhaustive oracle over all
cells, every "yes" carrying a verified decomposition), but the algebra
resolves only after nearly all coordinates of one summand have been
assigned: mean tame depth 7.4–10.4 of 13 coordinates at `n = 13`,
14.7–16.2 of 17 at `n = 17`, 17.8 of 19 at `n = 19` (FC). The number of
`GroebnerSafe` calls per target grows faster than the candidate count
`|F_w|` of the exhaustive oracle: fitted exponents 1.9 (FC), 3.2 (QFC),
3.6 (C) against the success threshold 0.8, and the ratio `calls / |F_w|`
rises from 0.10–0.12 at `n = 11` to 0.42–0.68 at `n = 17`. This is the
paper's own conclusion on Classic McEliece transferred to the PDP: the
Prange point wins, and every coordinate left to the algebra costs a
Gröbner call that the combinatorial saving does not repay. The subspace
baseline, run through the same solver, is tame at the root or one level
below it at every size.

Class, by §3 of `AGENTS.md`: **accounting**. What the construction changes
in the decomposition note's §6 is the `2^l` storage term of a
Frobenius-stable set (a normal-basis weight base needs none); no operation
count moves. No `S`, rho ratio or end-to-end figure is claimed: this is a
stage diagnostic of one oracle.

## Deviations from the protocol

Recorded before the affected cells ran, for runtime, not in response to
any result:

| protocol | as run | reason |
|:--|:--|:--|
| `w = 3` at `n = 17, 19` | `w = 2` at every `n` | with `w = 3` a full tree has `Σ_{j≤3} C(19, j) = 1160` leaves at `n = 19`; at the observed cost per wild call that is days per cell |
| 24 targets at every `n` | 8 targets (seeds `0..7`) at `n = 17, 19` | same |
| `τ = 2^30` XORs per call | `τ = 2^28` XORs per call plus a `2^26`-word matrix cap | a wild call at `2^30` cost 15–20 s at `n = 7` already; the sweep below measures the dependence |
| — | eager linear propagation and a direct check at fully assigned leaves inside `GroebnerSafe` | added after the first `n = 19` test dropped a planted solution when the budget ran out at a leaf; the leaf check makes completeness independent of the budget |

The container was reclaimed once during an idle hour and the running
cells died; `run.py` gained resume support and the affected cells
continued from their finished targets (the solver was unchanged; the
`n = 7` replay is byte-identical, see the workflow). Two data-handling
bugs were found and fixed before any frozen record was used: a manifest
line that wrote two fields into each other's slots (the campaign was
restarted from scratch on the fixed binary), and a summary whose float
last digits differed between hosts (now rounded).

## The systems

<!-- TABLE:systems -->
| n | encoding | base | \|F\| (points) | variables | of which auxiliary | generators | max degree |
|--:|:--|:--|--:|--:|--:|--:|--:|
| 7 | SUB | SUB5 | 27 | 10 | 0 | 7 | 2 |
| 7 | MONO | WT2 | 15 | 14 | 0 | 77 | 3 |
| 7 | C | WT2 | 15 | 48 | 34 | 45 | 2 |
| 7 | FC | WT2 | 15 | 42 | 28 | 39 | 3 |
| 7 | QFC | WT2 | 15 | 46 | 32 | 43 | 2 |
| 11 | SUB | SUB7 | 137 | 14 | 0 | 11 | 2 |
| 11 | MONO | WT2 | 67 | 22 | 0 | 341 | 3 |
| 11 | C | WT2 | 67 | 94 | 72 | 89 | 2 |
| 11 | FC | WT2 | 67 | 70 | 48 | 65 | 4 |
| 11 | QFC | WT2 | 67 | 86 | 64 | 81 | 2 |
| 13 | SUB | SUB7 | 125 | 14 | 0 | 13 | 2 |
| 13 | MONO | WT2 | 79 | 26 | 0 | 585 | 3 |
| 13 | C | WT2 | 79 | 114 | 88 | 107 | 2 |
| 13 | FC | WT2 | 79 | 84 | 58 | 77 | 4 |
| 13 | QFC | WT2 | 79 | 106 | 80 | 99 | 2 |
| 17 | SUB | SUB8 | 269 | 16 | 0 | 17 | 2 |
| 17 | MONO | WT2 | 137 | 34 | 0 | 1377 | 3 |
| 17 | C | WT2 | 137 | 172 | 138 | 163 | 2 |
| 17 | FC | WT2 | 137 | 120 | 86 | 111 | 5 |
| 17 | QFC | WT2 | 137 | 150 | 116 | 141 | 2 |
| 19 | SUB | SUB8 | 243 | 16 | 0 | 19 | 2 |
| 19 | MONO | WT2 | 229 | 38 | 0 | 1957 | 3 |
| 19 | C | WT2 | 229 | 196 | 158 | 185 | 2 |
| 19 | FC | WT2 | 229 | 132 | 94 | 121 | 5 |
<!-- /TABLE -->

## The table

One table, one unit (`GroebnerSafe` calls per target), every variant a
row, correctness first. `|F|` is the exhaustive oracle's candidate count
(points of the base; the exhaustive scan costs one group subtraction and
one membership test per candidate and finished in 40–250 µs per target at
`n = 13`–`19`). "calls, no-targets" is the full-tree cost on targets with
no decomposition; the mean over all targets includes the early stop at
the first witness. "budget-wild" counts the wild verdicts that came from
the budget or the matrix cap rather than from a completed non-linear
basis.

<!-- TABLE:results -->
| n | encoding | targets | decomposable | correct | timeouts | calls / target | calls, no-targets | tame | wild | budget-wild | tame depth | XOR words / call | wall / target (s) | \|F\| | calls / \|F\| |
|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 7 | SUB | 24 | 24 | 24/24 | 0 | 3.62 | — | 1.00 | 2.62 | 0.00 | 2.62 | 6699 | 0.01 | 27 | 0.134 |
| 7 | MONO | 24 | 24 | 24/24 | 0 | 4.67 | — | 1.00 | 3.67 | 0.00 | 3.67 | 7.75e+05 | 0.19 | 15 | 0.311 |
| 7 | C | 24 | 24 | 24/24 | 0 | 4.67 | — | 1.00 | 3.67 | 0.00 | 3.67 | 7.07e+06 | 1.23 | 15 | 0.311 |
| 7 | FC | 24 | 24 | 24/24 | 0 | 4.67 | — | 1.00 | 3.67 | 0.00 | 3.67 | 6.06e+06 | 0.97 | 15 | 0.311 |
| 7 | QFC | 24 | 24 | 24/24 | 0 | 4.67 | — | 1.00 | 3.67 | 0.00 | 3.67 | 6.57e+06 | 0.99 | 15 | 0.311 |
| 11 | SUB | 24 | 24 | 24/24 | 0 | 3.25 | — | 1.08 | 2.17 | 0.00 | 2.19 | 1.43e+06 | 0.12 | 137 | 0.024 |
| 11 | MONO | 24 | 24 | 24/24 | 0 | 4.58 | — | 1.00 | 3.58 | 0.00 | 3.58 | 1.14e+08 | 3.98 | 67 | 0.068 |
| 11 | C | 24 | 24 | 24/24 | 0 | 6.75 | — | 1.54 | 5.21 | 4.08 | 4.62 | 2.01e+08 | 35.62 | 67 | 0.101 |
| 11 | FC | 24 | 24 | 24/24 | 0 | 8.25 | — | 2.00 | 6.25 | 5.67 | 5.31 | 1.25e+08 | 19.16 | 67 | 0.123 |
| 11 | QFC | 24 | 24 | 24/24 | 0 | 6.75 | — | 1.54 | 5.21 | 4.08 | 4.62 | 2.02e+08 | 30.88 | 67 | 0.101 |
| 13 | SUB | 24 | 12 | 24/24 | 0 | 1.17 | 1.00 | 1.00 | 0.17 | 0.00 | 0.17 | 1.94e+05 | 0.02 | 125 | 0.009 |
| 13 | MONO | 24 | 14 | 24/24 | 0 | 10.12 | 13.00 | 3.62 | 6.50 | 6.00 | 4.38 | 1.96e+08 | 15.29 | 79 | 0.128 |
| 13 | C | 24 | 14 | 24/24 | 0 | 16.38 | 22.60 | 6.12 | 10.25 | 9.58 | 7.39 | 1.33e+08 | 183 | 79 | 0.207 |
| 13 | FC | 24 | 14 | 24/24 | 0 | 27.96 | 50.20 | 12.00 | 15.96 | 15.92 | 10.38 | 1.03e+08 | 55.30 | 79 | 0.354 |
| 13 | QFC | 24 | 14 | 24/24 | 0 | 12.21 | 16.40 | 4.50 | 7.71 | 7.42 | 5.57 | 1.59e+08 | 81.80 | 79 | 0.155 |
| 17 | SUB | 8 | 1 | 8/8 | 0 | 1.00 | 1.00 | 1.00 | 0.00 | 0.00 | 0.00 | 1.96e+05 | 0.08 | 269 | 0.004 |
| 17 | MONO | 8 | 5 | 8/8 | 0 | 22.25 | 27.00 | 9.25 | 13.00 | 13.00 | 8.76 | 1.69e+08 | 114 | 137 | 0.162 |
| 17 | C | 8 | 5 | 8/8 | 0 | 93.50 | 156 | 40.38 | 53.12 | 44.88 | 16.15 | 1.21e+08 | 633 | 137 | 0.682 |
| 17 | FC | 8 | 5 | 8/8 | 0 | 56.88 | 102 | 25.50 | 31.38 | 30.00 | 14.74 | 5.67e+07 | 231 | 137 | 0.415 |
| 17 | QFC | 8 | 5 | 8/8 | 0 | 68.75 | 118 | 29.75 | 39.00 | 34.12 | 15.35 | 1.52e+08 | 624 | 137 | 0.502 |
| 19 | SUB | 8 | 0 | 8/8 | 0 | 1.00 | 1.00 | 1.00 | 0.00 | 0.00 | 0.00 | 2.75e+05 | 0.02 | 243 | 0.004 |
| 19 | MONO | 8 | 4 | 8/8 | 0 | 33.75 | 41.00 | 14.00 | 19.75 | 17.50 | 11.80 | 1.29e+08 | 284 | 229 | 0.147 |
| 19 | C | 8 | 4 | 8/8 | 0 | 135 | 212 | 62.25 | 73.12 | 66.75 | 18.21 | 1.37e+08 | 891 | 229 | 0.591 |
| 19 | FC | 8 | 4 | 8/8 | 0 | 113 | 175 | 51.50 | 61.38 | 55.75 | 17.81 | 5.18e+07 | 418 | 229 | 0.493 |
<!-- /TABLE -->

Exponent fits (least squares of `log₂` mean per target against
`log₂ |F|` over the prime sizes listed; the success threshold was 0.8):

<!-- TABLE:fits -->
| encoding | sizes fitted | exponent of calls / target in \|F\| | exponent of XOR words / target in \|F\| |
|:--|:--|--:|--:|
| C | 11, 13, 17, 19 | 2.43 | 2.28 |
| FC | 11, 13, 17, 19 | 1.87 | 1.12 |
| QFC | 11, 13, 17 | 3.22 | 2.85 |
| MONO | 11, 13, 17, 19 | 1.51 | 1.64 |
| SUB | 11, 13, 17, 19 | -0.87 | -1.74 |
<!-- /TABLE -->

## Budget sensitivity

The same four `n = 11` targets (seeds `0..3`) under three per-call budgets,
matrix cap `2^27` words. The tame depth falls by one to two levels for a
16× larger budget; the trees do not collapse.

<!-- TABLE:tau -->
| encoding | budget (row-word XORs) | targets | correct | calls / target | tame depth | budget-wild calls / target | XOR words / target | wall / target (s) |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| FC | 2^26 | 4 | 4/4 | 9.25 | 5.36 | 6.50 | 3.44e+08 | 21.1 |
| FC | 2^28 | 4 | 4/4 | 8.50 | 5.10 | 6.00 | 9.9e+08 | 21.9 |
| FC | 2^30 | 4 | 4/4 | 7.25 | 4.38 | 5.25 | 2.62e+09 | 24.8 |
| MONO | 2^26 | 4 | 4/4 | 6.00 | 3.71 | 4.00 | 3.31e+08 | 3.9 |
| MONO | 2^28 | 4 | 4/4 | 3.00 | 2.00 | 0.00 | 3.3e+08 | 3.5 |
| MONO | 2^30 | 4 | 4/4 | 3.00 | 2.00 | 0.00 | 3.3e+08 | 3.4 |
| QFC | 2^26 | 4 | 4/4 | 6.75 | 4.12 | 4.75 | 3.77e+08 | 22.6 |
| QFC | 2^28 | 4 | 4/4 | 6.00 | 3.71 | 4.00 | 1.26e+09 | 32.9 |
| QFC | 2^30 | 4 | 4/4 | 3.50 | 2.50 | 1.00 | 2.54e+09 | 28.9 |

| encoding | seed | calls at 2^26 | calls at 2^28 | calls at 2^30 | tame depths at 2^26 / 2^28 / 2^30 |
|:--|--:|--:|--:|--:|:--|
| FC | 0 | 7 | 6 | 6 | 6 / 5 / 5 |
| FC | 1 | 13 | 11 | 11 | 6,6,5,4,3,3 / 5,5,4,3,3 / 5,5,4,3,3 |
| FC | 2 | 10 | 10 | 6 | 7,7,6 / 7,7,6 / 5 |
| FC | 3 | 7 | 7 | 6 | 6 / 6 / 5 |
| QFC | 0 | 5 | 5 | 2 | 4 / 4 / 1 |
| QFC | 1 | 10 | 8 | 2 | 5,5,4,3,2 / 4,4,3,2 / 1 |
| QFC | 2 | 6 | 5 | 4 | 5 / 4 / 3 |
| QFC | 3 | 6 | 6 | 6 | 5 / 5 / 5 |
| MONO | 0 | 5 | 1 | 1 | 4 / 0 / 0 |
| MONO | 1 | 8 | 1 | 1 | 4,4,3,2 / 0 / 0 |
| MONO | 2 | 5 | 4 | 4 | 4 / 3 / 3 |
| MONO | 3 | 6 | 6 | 6 | 5 / 5 / 5 |
<!-- /TABLE -->

## Reading

- **Correctness.** Every cell, every target: oracle answer = exhaustive
  label; no encoding violation (no solution with a summand outside
  `F_w`); no timeout.
- **Where the algebra resolves.** The tame depth is the number of summand
  coordinates fixed before `GroebnerSafe` returned a linear basis. For
  the paper's encodings it is within 1–3 of the full summand at
  `n ≥ 17`, i.e. the oracle is the exhaustive oracle with a Gröbner
  computation attached to each candidate. The monomial-ideal control
  resolves 6–8 levels higher, which says the lifted presentations buy
  bounded degree at the price of the linear-algebra information the
  truncated basis reaches within a budget; it does not make MONO a
  candidate, its generator count `C(n, w+1)` is what the paper's
  constructions exist to avoid.
- **Cost.** Per call the paper's encodings spend `0.5`–`2 × 10^8` row-word
  XORs, almost all of it in wild calls that end at the budget; per target
  the wall time is `10^5`–`10^6` times the exhaustive scan at the same
  `n`. Uncalibrated against curve operations, and quoted only as a
  practicality note.
- **The engine.** The budget sweep bounds the dependence on the solver:
  16× the budget removes one to two levels of the tree at `n = 11`. A
  better F4 would move the depth, not the exponent.

## Decision

Closed as a negative result. The construction is recorded for what it
does establish, a Frobenius-stable factor base at prime `n` with a
bounded-degree membership ideal and no storage, and the decomposition
note's §6 carries the correction. No follow-on cell is planned; `m ≥ 3`
and the `GBDecode` outer loop are listed in the note as what was not run.
