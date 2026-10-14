# Cycle 4b — compact-orbit columns sweep

Session: `sessions/cycle4b` (`ECBS1h5ae06f07ab10`), 48/48 verified, audit OK.

## Hypothesis

On E₀ n=13, `compact-orbit-scan:columns=K` + `mitm-frobenius` gives a
Frobenius-folded pair table (~|F|²/(4n) group adds) and can undercut
negfold cold S.

## Mean cold S (2 public × 2 rounds)

| arm | mean S | lower bound | note |
|---|---:|---|---|
| rho-strong | 14.324 | yes | |
| **ic-c2-frob-m2** | **14.096** | yes | many FB natives unpriced |
| ic-f4-m3-d5 | 15.349 | yes | XORs unpriced |
| ic-c3-frob-m2 | 16.545 | yes | |
| ic-c3-frob-m3 | 16.533 | yes | |
| ic-c4-frob-m2 | 19.561 | yes | |
| ic-c4-frob-counted-m2 | 22.140 | yes | canons priced into setup |
| ic-negfold-m2-d5 | 27.491 | **no** | complete S |

## Decision

**Win vs negfold.** Candidate for parity is columns=2 (and counted).
Promote counted columns=2/3 to cycle-5. Do not treat uncounted c2 as a
complete-S win: FB construction natives and inversions remain unpriced.
