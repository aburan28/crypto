# Stage 184 protocol: dense pair selection on the current default path

## Purpose

Stage 183 accepted exact dense pair selection when stacked with opt-in full
M4RI: median paired wall fell to `0.7368` and CPU to `0.8023`. Before changing
the repository default, replay the selector with current five-column
`BlockTables`, which remains the selected elimination policy.

## Frozen arms

- Exact Stage 183 binary SHA-256:
  `205856a1923448f43a0c51e5c1527aa0141d8437493bab2d80676158524bf6e4`.
- Binary source commit:
  `fdeb3155fe2c9757d30464bae15584524a0554c8`.
- Both arms set `F4_F2_FULL_M4RI=0`, twelve Rayon workers, X1 batch 512, and all
  external thread controls to one.
- Control: `F4_F2_DENSE_PAIR_SELECT=0`.
- Candidate: `F4_F2_DENSE_PAIR_SELECT=1`.
- Same opened `n=59, ell=9, m=3` target and equation fingerprint
  `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.

Run one control/candidate screen. Continue only if candidate CPU falls and all
correctness gates hold. Confirmation order is
`quadratic, dense, dense, quadratic, quadratic, dense`.

## Correctness and accounting

Both arms must return exhaustive UNSAT, visit all 512 masks, complete all 242
systems, and preserve exact current-BlockTables logical/performed XOR, pair,
basis, matrix, extraction, and equation counts. Dense selection must repeat
Stage 183's selector-call, candidate, LCM-group, cover-lookup, and scratch-byte
counts; the control must report the corresponding quadratic counts.

Carry the Stage 183 exact build and validation cost as inherited setup. Charge
every new query process, total core-seconds, wall, RSS, and all twelve workers.

## Selection and follow-up

Select dense pair selection as the repository default only if median paired
dense/quadratic wall and total-core ratios are both below `0.97`. If selected,
change the environment contract so unset means dense, `F4_F2_DENSE_PAIR_SELECT=0`
is the quadratic control, and `=1` explicitly selects dense. Then run both F4
and backend tests plus a default-mode target replay from the post-selection
commit before pushing.

This remains one-target implementation engineering. It is not relation-yield
evidence, a full DLP, a rho crossover, independent reproduction, novelty, or
SOTA.
