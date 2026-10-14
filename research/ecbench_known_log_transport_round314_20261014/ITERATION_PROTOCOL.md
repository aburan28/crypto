# Frozen post-holdout iteration

Registered after the sequential holdout was audited and before target seed
`20261014314004` was opened.

The sequential holdout failed the parity gate at every size.  Its two
dominant phases were the folded pair-table setup (40--52% of cold cost) and
target PDP (45--60%).  Observed target attempts grew from 6.825 to 56.025,
far above the class-count occupancy model, while the selected representative
logs were consecutive small integers.  The sole allowed engineering
iteration changes only that schedule:

- new method identity `ic.transported_logs_hash`;
- deterministic hash-selected nonzero representative scalars using base seed
  `20261014314005`;
- every representative `[s]G` scalar multiplication is charged;
- the same frozen K values 11, 21, 43, 71;
- fresh public targets from seed `20261014314004`;
- same cold boundary, exact column/witness/scalar replay, rho reference, and
  five measured rounds after one warmup;
- an exact A/A hash control, the sequential method, and the existing claw as
  contextual baselines.

No parameter or implementation change is permitted after these sessions.
Failure publishes a negative result and ends Round314.
