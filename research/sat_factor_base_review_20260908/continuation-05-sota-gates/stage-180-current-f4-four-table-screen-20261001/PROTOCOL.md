# Stage 180 protocol: four-column BlockTables target screen

## Question

Does the repository's existing `f4-four-tables` feature reduce deterministic
performed XORs and charged CPU on the one frozen Phase B target relative to the
five-column default?

Stage 179 rejected six columns because performed XORs rose 10.63% and table
memory rose 39.65%. Four columns use smaller tables and may trade less reduction
per hit for much cheaper construction and better cache behavior.

## Frozen screen

- Reuse the exact Stage 179 five-column binary with SHA-256
  `b156f3320564c1e9fed2b168aba7a0b4f3c875d677c1a9354751bba527ee5bd6`.
  Carry its charged clean-build cost into this stage as inherited setup.
- Build a four-column binary from the exact commit containing this protocol
  with `--features f4-four-tables`, under the process meter.
- Same authenticated `n=59, ell=9, m=3` target, twelve Rayon workers, X1 batch
  512, 300-second internal budget, and 360-second watchdog.
- Run one two-order screen: `five, four, four, five`.
- Every run must return exhaustive UNSAT, reproduce the equation fingerprint,
  visit all 512 masks, complete all 242 systems, and preserve every logical F4
  counter. Performed XORs and table memory may differ.

## Stop and decision rule

Proceed to a three-pair confirmation only if four columns perform fewer actual
word XORs than five columns and both one-pair CPU ratios are below `0.97`.
Otherwise stop and reject table-width tuning. A screen cannot be promoted as a
selected optimization by itself.

Report and charge the inherited five-column build, fresh four-column build,
all four query processes, wall, core-seconds, RSS, performed XORs, and table
memory. Single-core time remains null. This target-specific solver diagnostic
does not change any index-calculus or SOTA gate.
