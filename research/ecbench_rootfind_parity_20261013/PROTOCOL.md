# Protocol: ecbench rootfind IC row-parity loop (2026-10-13)

Extends PR #1774 / `research/ecbench_rootfind_improve_20261013/`.
Classification target: **engineering**. Synthetic public targets only.
No key-recovery, asymptotic sub-rho, or deployed-curve claims.
Rootfind `S` is a complete count; inherited-F4 `S` is a lower bound
(unpriced word XORs) — do not divide one by the other as a speedup.

## Goal

Admitted cold `S` for a rootfind IC arm ≤ best F4 IC arm on the same
panel (claim-n13 baseline: F4 m=3/d=5 mean S **16.322**). Stretch toward
strong rho (~14.3).

## Starting gap (improve claim-n13)

| arm | mean S | note |
|---|---:|---|
| rho-strong | 14.342 | reference |
| ic-f4-m3-d5 | 16.322 | best IC row (lower-bound S) |
| ic-rootfind-m2-d4 | 28.461 | best cold rootfind so far |

Phase split: rootfind m2/d5 online ≈ F4 online; set-up `|F|(|F|+1)/2`
dominates the gap.

## Levers (this loop)

1. **negation-folded PairTable** (`engine=rootfind-negfold`) — ~|F|²/4 set-up
2. **window + negfold** (`rootfind-window-negfold`)
3. Further: |F|/m schedule, honest multi-target amortization (secondary),
   prime oracle panel

## Conductor

Local Conductor daemon was unreachable at start (`connection refused` on
127.0.0.1:8443). Proceeding under owner instruction to execute the loop;
scopes intended: `research/ecbench_rootfind_parity_20261013/`,
`src/bin/ecbench/pdp3_koblitz.rs`, `docs/ecbench/specs/`, methods help.

## Unit

Floor `√(π/2A)` in `S`; reference `rho.signed_frobenius_strong`; cold `S`
with every phase charged. Amortized figures that share `oracle_setup`
are **secondary**.
