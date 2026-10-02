# Stage 186 protocol: elimination policy after dense-pair selection

## Question

Stage 181 compared current BlockTables and full M4RI while both used the old
quadratic pair selector; full M4RI cut CPU but failed the paired wall gate.
Stages 183–185 selected dense exact pair updates and materially changed the
non-elimination runtime. Re-evaluate elimination on the selected dense-pair
default before choosing the Phase B target-specific configuration.

## Frozen setup

- Exact post-selection binary SHA-256:
  `3aefd5cd0daf73602bf8a38e7fdcb363e630ffd57bfebe47a37123b1d03429b1`.
- Selection commit:
  `8014149a2b55cd2cca202237a302643e39a50f6e`.
- Supplied/restored lock SHA-256:
  `4f17b356fa7bac392b6d801d1c74fb9e36b6517f9465c8ebc19bb9a2792a84c5`.
- Dense pair selector is left unset and must route to dense in both arms.
- Control: `F4_F2_FULL_M4RI=0` (current five-column BlockTables).
- Candidate: `F4_F2_FULL_M4RI=1` (block-8 full-matrix M4RI).
- Twelve Rayon workers, X1 batch 512, all other thread controls one, same frozen
  `n=59, ell=9, m=3` target and equation fingerprint.

Run one control/candidate screen. Continue only if candidate CPU and performed
XORs fall. Confirmation order is
`blocktables, full, full, blocktables, blocktables, full`.

## Correctness and decision

Both arms must report the selected dense pair counters, exhaustive UNSAT, all
512 masks / 242 systems, exact equation fingerprint, and valid source custody.
The elimination-specific logical/performed XOR and matrix/block counts must
repeat their Stage 184 and Stage 183 values respectively.

Select full M4RI only as the Phase B target-specific multi-worker configuration
if median paired full/blocktables wall and total-core ratios are both below
`0.97`. Do not change the repository-wide elimination default from this one
target; current BlockTables remains the general default and single-core policy.

Carry Stage 185 build/validation custody and charge every new query. This is
still solver-stage engineering, not relation yield, a full DLP, rho crossover,
independent reproduction, novelty, or SOTA.
