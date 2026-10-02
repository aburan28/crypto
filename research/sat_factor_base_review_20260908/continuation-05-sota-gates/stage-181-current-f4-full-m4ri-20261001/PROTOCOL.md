# Stage 181 protocol: full-matrix M4RI inside current F4

## Hypothesis

The current five-column `BlockTables` path performs `147,794,583,858` actual
word XORs on the frozen target. Stages 179–180 show that changing the local
table width cannot improve this: four columns perform 12.34% more work and six
columns 10.63% more.

The older Phase B engine used a full-matrix Method of Four Russians schedule
with blocks of eight pivot rows and materially lower measured CPU. This stage
ports only that elimination schedule as an opt-in path inside the current F4.
The current pair queue, symbolic preprocessing, bitmap/hash sets, matrix
construction, field pairs, extraction, and final interreduction remain
authoritative.

## Implementation boundary

- `F4_F2_FULL_M4RI=1` enables the candidate only for matrices with at least
  128 rows, 256 columns, and at most four columns per row; all other matrices
  use current `BlockTables` elimination.
- Candidate block width is eight. Scratch tables are thread-local and reused.
- `word_xors` counts the individual pivot-row XOR work represented by a table
  lookup. `word_xors_performed` counts actual row and table-construction XORs.
  The units remain explicit even if the candidate chooses a different pivot
  basis and therefore has a different logical count.
- Table scratch memory, matrices routed through the candidate, blocks, and
  table-construction XORs are exported and included in peak-memory accounting.
- The default remains the current five-column path unless this protocol passes.

## Correctness gates

1. Random matrix differential tests compare full M4RI with ordinary Gaussian
   elimination for rank, pivot columns, and canonical row space across varied
   shapes, including partial blocks and word boundaries.
2. The candidate's logical XOR counter equals the sum of represented
   individual pivot reductions; actual XORs include table construction.
3. All current Boolean-F4 certified-basis, Buchberger-agreement, budget,
   extraction, field-pair, and existing BlockTables tests pass both with the
   candidate disabled and forced on.
4. Fixed-X1 specialization and backend witness tests pass in both modes.
5. On the frozen target, both modes authenticate the same source, visit all
   512 masks, complete all 242 systems, reproduce equation fingerprint
   `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`,
   and return exact exhaustive UNSAT. Any invalid model, timeout, oversize, or
   result disagreement rejects the candidate.

## Benchmark and decision

First run one current/candidate screen on the same exact-commit binary with
twelve Rayon workers, X1 batch 512, 300-second internal budget, and 360-second
watchdog. Continue to three interleaved pairs only if the candidate reduces
both performed XORs and total core-seconds in the screen. Confirmation order is
`current, full, full, current, current, full`.

Select full M4RI only if the three-pair median full/current wall and total-core
ratios are both below `0.97`, correctness gates pass, and RSS is reported.
Charge the clean build, tests, screen, confirmation, failures, and every worker
through core-seconds. Single-core time remains null in this experiment.

This remains one-target solver engineering. It cannot establish relation yield,
a full unknown-scalar index-calculus run, a rho crossover, independent
reproduction, novelty, or SOTA.
