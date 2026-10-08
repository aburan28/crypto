# Audited single-target continuation

Three bounded tournaments follow merged [PR #372](https://github.com/aburan28/crypto/pull/372). Every cost is for one complete cold ECDLP target. These are engineering measurements on five small Koblitz cells, not a cryptographic-size crossover or a new complexity exponent.

## Decisions

| Round | Locked challenger | Decision | Audited receipts | Verified solves |
|---|---|---|---:|---:|
| [round-0006b-single](../runs/round-0006b-single/REPORT.md) | combined_batch1 | retained: incumbent | 1584 | 1584 |
| [round-0007-single-policy](../runs/round-0007-single-policy/REPORT.md) | orbits2 | promoted: orbits2 | 1560 | 1542 |
| [round-0008-single-implementation](../runs/round-0008-single-implementation/REPORT.md) | combined_dense | retained: incumbent | 1536 | 1536 |

Eighteen policy smoke trials were rejected by the frozen checker's requirement for at least two logarithm columns: all twelve `orbits1` trials and six `cube_root` trials. Those candidates were excluded before development. This does not establish that their final scalars were wrong. All other scheduled solves verified, including paired native outputs. The earlier round-0006 preparation failed on CPU affinity before fixtures or measurements; its sources/builds and failure are retained.

## Complete user-space instructions (Ir)

| Round | Challenger / incumbent | Paired 95% interval | Replay / incumbent |
|---|---:|---|---:|
| round-0006b-single | 0.785351 | [0.7704334287179592, 0.7973963766888058] | 0.785349 |
| round-0007-single-policy | 0.375287 | [0.32387886079328404, 0.45379932682481017] | 0.375288 |
| round-0008-single-implementation | 0.956434 | [0.9474512194228004, 0.966125671872147] | 0.956431 |

## Complete native process time

| Round | Challenger / incumbent | Paired 95% interval | Replay / incumbent |
|---|---:|---|---:|
| round-0006b-single | 0.884578 | [0.873064358144668, 0.8995482205334346] | 0.891162 |
| round-0007-single-policy | 0.709068 | [0.6853953419318904, 0.735404642949192] | 0.716704 |
| round-0008-single-implementation | 0.985788 | [0.981444289834538, 0.990511305956245] | 0.983780 |

Round 0007 promoted the two-orbit policy with 62.47% fewer instructions and 29.09% less native time than the prior promoted single-target baseline. Its matched rho ratios were 1.1520 instructions and 1.1008 time. The policy panel deliberately changes B and rank floors; it cannot establish an advance against an unchanged floor. Round 0008 fixes that promoted support and tests implementation changes. Do not multiply ratios from different fixture sets.

## Strict rho barrier: final fixed-support round

Strict beating requires both upper paired 95% confidence limits and every cell ratio below one in confirmation AND replay. This is stronger than the legacy instruction point estimate and the 10% parity margin. The following is separate from the >=20%-in-both-metrics IC promotion gate.

### combined_dense: strict beating = False

| Stage | Ir / rho | Paired 95% interval |
|---|---:|---|
| confirmation | 1.078264 | [1.024533243343095, 1.1410257828499533] |
| replay | 1.078263 | [1.0245319304796583, 1.141022551343304] |

| Stage | Native time / rho | Paired 95% interval |
|---|---:|---|
| confirmation | 1.066185 | [1.0491683332027872, 1.085519347488403] |
| replay | 1.065578 | [1.0521588300075908, 1.0834709390160266] |

### incumbent: strict beating = False

| Stage | Ir / rho | Paired 95% interval |
|---|---:|---|
| confirmation | 1.127380 | [1.0664394518628761, 1.1978968142897077] |
| replay | 1.127382 | [1.0664410849430652, 1.1978933413637658] |

| Stage | Native time / rho | Paired 95% interval |
|---|---:|---|
| confirmation | 1.081557 | [1.0672187767654386, 1.1006507319022558] |
| replay | 1.083147 | [1.06860936740454, 1.102554906966594] |

## Boundary, evidence and validation

For m=3 and signed base size B, coverage is at most binomial(B+2,3) target images. The implemented K-column full-rank collector needs at least K relation-producing trials; its K-instruction floor is deliberately weak. Each round report contains all variants, normalized Ir / sqrt(r), ratios to rho and each applicable floor, phase totals, admission failures and native confidence intervals. All startup, setup, collection attempts, linear algebra, certification, descent and final verification remain charged. Kernel/device work, the profiler and external audit are excluded from Ir. Rho uses the existing signed-Frobenius per-target API.

Validation: 21 runner tests; clean-main worker cargo check; exact orbit-table and orbit-projection equivalence tests; every completed round's frozen audit; fresh-directory archive restoration and full frozen audits. See [archive-validation.json](archive-validation.json), [implementation-scope-check.json](implementation-scope-check.json), and [evidence restoration](../evidence/README.md).

The cumulative [WINNER.patch](WINNER.patch) is relative to `runs/round-0002/source`; [WINNER.json](WINNER.json) identifies the actual measured source and configuration. Production library defaults are unchanged. All rejected candidates, preparation failure and raw profiles are preserved. The bounded operation stops after round 0008 for publication, regardless of its verdict.
