# Version 3: charged deterministic six-class restarts

Frozen before v3 execution. [V2](REMEDIATION_V2.md) emitted all 80 rows.
All 48 relation rows and all 16 negation-rho rows verified; 13 of 16
six-class rows verified. Six-class rho again exhausted its 100,000-operation
cap on `p=13` target indices 4, 5, and 6. The v2 files stay unchanged.

V3 retains v2's curves, exact targets, maps, relation streams, rho walk
policy with `dp_bits=0`, and arithmetic checks. The six-class arm may make
up to 16 independent attempts for each target. Attempt `j`, starting at
zero, uses seed
`(RHO_SEED ^ p ^ (target_index << 16)) ^ (j * 0x9e3779b97f4a7c15)`
with wrapping 64-bit arithmetic. Every attempt uses the 100,000-operation
walk cap. Stop at the first verified answer; if all 16 fail, retain that
target as failed. Sum every attempt's setup, walk, verification, field,
canonicalization, and native counters in the target's row. Also retain a
per-attempt receipt in the run log, with seed, status, steps, group
operations, distinguished points, and canonicalizations. Charge the public
endomorphism-scalar selection once per target in the setup phase.

The six-class comparison uses the aggregate cold cost, including failed
attempts. V3's different restart policy is explicit and does not revise
v1/v2 evidence. Success requires all 80 rows verified, byte-identical
replay, and an independent checker of the scalar and orbit action. Neither
the fixed 16-attempt schedule nor these two group orders calibrates the
best six-class rho algorithm or a scaling law.
