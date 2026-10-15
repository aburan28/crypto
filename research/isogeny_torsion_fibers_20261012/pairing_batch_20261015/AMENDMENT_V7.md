# V7 direct source-group BSGS reuse control

Frozen before implementation and measurement. The V6/V6B three-arm runs
charged reuse of a finite-field BSGS table over 66 distinct order-67
targets, against direct `F_401` rho. On a multi-target workload the source
group can also reuse a baby-step table. V7 prices this cheaper public route
on the exact same target order, with no answer scalar in solver input.

Inherit the [V6 public cell and caps](RUN_SPEC_V6.json): `p=401`, `r=67`,
`G=(251,346)`, 66 distinct public nonidentity targets, prefixes
`N=1,4,16,64,66`, rounds 0–2, charged input parsing and subgroup checks,
and 32 rho retries at 144 steps per attempt. Freeze the new source and
binary hash in `RUN_SPEC_V7.json` before the first V7 run. Preserve all
earlier versions and the failed V6B cell-index invocation.

The three V7 arms are (1) cold pairing with H search and reusable
multiplicative-field BSGS, (2) direct `F_401` negation rho, and (3) cold
direct `F_401` elliptic-curve BSGS. The direct BSGS arm constructs
`ceil(sqrt(67))=9` baby entries from `[j]G`, computes `[9]G`, and for
each public `Q` checks `Q-[i*9]G` for `0≤i≤9`. Charge its table construction,
every lookup, group operation, failed candidate, per-target verification,
memory and cold input validation. Its table is reused only within the
recorded prefix; there is no unseen precomputation. Rotate arm order:
round zero pairing/rho/BSGS, round one BSGS/rho/pairing, round two
rho/pairing/BSGS. Emit JSON only after all timing intervals complete.

Exactness gate: all arms must return the same 66 scalars each round;
independent source-curve arithmetic must replay all `3*3*66=594` outputs,
including the answer fixture only in replay. A missing baby/giant match,
rho exhaustion or replay mismatch censors that round and preserves its
failure record. Stop after three rounds and one independent replay.

The comparison unit is complete cold current-thread CPU and wall with the
same public-input charge. Keep base-field group/field counts and peak RSS,
and report all five prefixes. The natural reference for multi-target
amortization is the lower complete cost of direct rho and direct BSGS at
each prefix. The run is shared-host L0; no isolated wall-time, portable
counted-`S` or growing-size claim follows from it.
