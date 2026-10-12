# Stage 183 protocol: dense exact pair selection in current F4

## Hypothesis

Current F4 selects new Gebauer–Möller pairs by comparing each active candidate
LCM against the remaining candidates, a quadratic scan repeated for every new
basis element. On the frozen 18-variable target the complete monomial domain
has only `2^18` masks. Grouping equal LCMs in reusable dense arrays and finding
proper divisor LCMs by exact submask lookup should preserve every UPDATE
decision and pair order while reducing CPU outside elimination.

This ports only new-pair selection. Pending-pair filtering, current symbolic
preprocessing, field pairs, matrix construction, extraction, and the Stage 181
full-M4RI elimination remain unchanged.

## Candidate and control

- Both arms set `F4_F2_FULL_M4RI=1`, twelve Rayon workers, and X1 batch 512.
- Control: `F4_F2_DENSE_PAIR_SELECT=0`, current quadratic selector.
- Candidate: `F4_F2_DENSE_PAIR_SELECT=1`, reusable dense exact selector for at
  most 20 variables; larger domains fall back to the current selector.
- Export selector calls, candidate visits, LCM groups, submask cover probes,
  and peak dense scratch bytes. Include scratch in solver peak memory.

## Correctness gates

1. A randomized differential test compares dense and quadratic selected pairs,
   pair order, chain/product skip counts, and duplicate-LCM representative over
   varied active bases and new leading monomials.
2. All Boolean-F4 tests pass with dense selection off and forced on, including
   certified bases, Buchberger agreement, full-M4RI row-space equivalence, and
   budget handling.
3. Backend tests pass in both modes.
4. On the frozen target both arms reproduce exact exhaustive UNSAT, all 512
   masks / 242 rational systems, the equation fingerprint
   `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`,
   and the Stage 181 full-M4RI logical/performed XOR and matrix/block counts.
   Pair, basis, matrix, and extraction counts must match exactly.

## Benchmark and decision

Run one control/candidate screen from one exact-commit binary. Continue only if
candidate total core-seconds falls and correctness holds. Confirmation uses
three fresh interleaved pairs in order
`quadratic, dense, dense, quadratic, quadratic, dense`.

Select dense pair selection only if median paired dense/quadratic wall and
total-core ratios are both below `0.97`. Report RSS and every setup/query cost.
The full-M4RI stack remains a solver-stage research arm until its own selection
gate passes; a stacked win cannot be reported as full index-calculus speed.

This one-target experiment cannot establish natural relation yield,
unknown-scalar recovery, a rho crossover, independent reproduction, novelty,
or SOTA.
