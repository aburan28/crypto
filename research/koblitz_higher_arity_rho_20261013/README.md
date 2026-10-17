# Higher-arity Koblitz extraction versus matched rho

**Outcome:** the two-S3-chain kernel is a correct and useful constant-factor
replacement for the quartic root finder, but it does not remove the scaling
blocker. On the registered n=41 Koblitz rung, its full cold counted cost is
`S = 7.652`, or `55.447x` the same-target strong-rho baseline (`S = 0.138`).
Both figures are lower bounds because native work remains unpriced.

![Counted cost and isolated extraction wall ratio](figures/n41-summary.png)

## What changed

Kernel `b3` visits exactly the pair cells used by the existing diagonal-stride
`b2` kernel. For each cell `(s1,s2)`, it replaces

`roots(S4(xR, x(s1), x(s2), X))`

with the sorted, deduplicated union of

`roots(S3(x(R+s1), x(s2), X))` and
`roots(S3(x(R-s1), x(s2), X))`.

The verifier, factor base, pair order, public target, linear algebra and final
scalar replay are unchanged. Tests check the root-set identity on split
quartics, planted six-point recovery, random targets, non-decomposable targets
and budget exits. The oracle is registered in the `ecbench` binary so the new
kernel and strong rho run through one measurement surface.

## Results

### Isolated one-relation extraction

Both runs used the same n=41, K=20 base and the same seed. Both reached the
same verified six-point relation after 78,041 pair cells on reserved CPU 4;
neither isolation record was marked contended.

| kernel | root work | charged point adds | extraction wall | process wall |
|---|---:|---:|---:|---:|
| B2 quartic | 78,041 S4 | not charged | 0.8756 s | 0.9113 s |
| B3 chain | 156,082 S3 | 234,129 | 0.1360 s | 0.1838 s |

B3 is `6.44x` faster inside extraction (`4.96x` process wall), a real
implementation improvement on this host. It is not an algorithmic speedup:
the replacement materializes `s1`, `R+s1` and `R-s1`, charging three group
additions for each unchanged pair cell.

### Matched end-to-end ecbench session

Session `ECBS1h4f09e1612522` used one public point, one warmup and three paired
measured rounds. All 20 planned executions verified. Exact replay reproduced
all 15 measured deterministic runs; audit receipt
`767cf32c4d96d4a2720cbb38e6fbea98e9e018a97c0a628faf16e435ec8a7307`
has zero problems.

| arm | mean S | interval | ratio to rho | key unpriced work |
|---|---:|---:|---:|---|
| strong rho | 0.1380 | [0.088, 0.202] | 1.000 | canonicalization, hashes, table operations |
| B2 quartic | 0.7609 | [0.662, 0.888] | 5.513 | 5,109,566 S4 pairs across measured runs |
| 4-sum | 5.5241 | [4.327, 7.190] | 40.027 | 77,136,383 S3 solves across measured runs |
| B3 chain | 7.6521 | [6.616, 8.989] | 55.447 | 10,219,132 S3 solves across measured runs |

The matched rho mean is `0.997x` its signed-Frobenius theory floor
`sqrt(pi/(2*82)) = 0.138405`. The duplicate rho arm reproduced operation
counts exactly (`1.000x`). B2 and B3 found the same 58 relations with the same
rank and trial counts across the measured runs. B3 used the same 5,109,566
pair cells as B2 and therefore incurred 15,328,698 additional point additions.
Its counted B3/B2 ratio is `10.057x`, while its descriptive pipeline wall
ratio is `0.150x`: cheaper field root finding, more charged group work.

## Current blocker

### End-to-end elapsed-time follow-up

A second frozen session, `ECBS1hb77d7527b579`, ran the complete cold pipeline
for five measured rounds on the same public n=41 target. All 24 executions
verified and exact replay reproduced all 20 measured deterministic runs. The
host again earned L0 because CPU and memory PSI telemetry were unavailable, so
these elapsed times are descriptive.

| arm | cold solve mean | process mean | online mean | dominant phase |
|---|---:|---:|---:|---|
| strong rho | 13.005 ms | 34.053 ms | 10.320 ms | rho solve |
| B3 two-S3 chain | 3.137 s | 3.185 s | 3.106 s | relations: 3.105 s (98.981%) |
| B2 quartic | 22.536 s | 22.588 s | 22.505 s | relations: 22.504 s (99.858%) |

B3 uses `0.1392x` B2's cold solve wall (`7.183x` less), but still uses
`241.237x` rho's cold solve wall on the ratio of recorded sums. The apparently
conflicting counted and elapsed-time results measure different work: B3 charges
three group additions per unchanged pair cell, while B2's S4 and B3's S3 field
root work remain unpriced in `S`. Thus B3 is `10.121x` B2 in counted `S` yet
`7.183x` faster in descriptive wall time.

The phase records localize the implementation bottleneck. Across the five
measured runs, both kernels produced the same 98 relations from the same
10,128,845 pair cells and 40,515,237 lookups. That is 2,025,769 cells per cold
solve and 103,355.6 cells per accepted relation. B3 reduces processing from
11.109 us to 1.533 us per cell, but it does not reduce the number of cells.

The complete run report, graph and raw links are in
[`reports/run-20261013T2207Z`](reports/run-20261013T2207Z/).

The blocker is now narrower than “make S4 faster.” The extraction shape still
enumerates about 88,096 pair cells per accepted relation in the first session
and 103,355.6 in the five-round follow-up. Replacing one quartic with two S3
solves makes each cell faster in wall time, but adding three group operations
per cell overwhelms rho before the unpriced S3 and lookup work is included. A
scalable candidate must reduce the number of cells or amortize point
construction across many cells/targets; another constant-factor root-finder
substitution on the same scan will not close the gap.

The required n=83 confidence gate cannot run in the native pipeline yet.
`ecbench fb` fails closed with:

> no factor-base plug-in runs on a wide (n > 62) curve yet

The dispatcher likewise permits only strong rho and `claw.pair_table` on
64--127-bit fields. Consequently this experiment establishes no n=83 scaling
exponent and no ECC2K-130 transfer claim.

## Evidence status and limitations

- Counts are primary. The session is L0 because PSI CPU and memory telemetry
  are unavailable in this sandbox. Ecbench first refused the run, then the
  explicit `--allow-busy` downgrade was used. All wall ratios are descriptive.
- The isolated B2/B3 probes passed their own quiet-core and contention checks,
  but there is only one run per arm; their wall ratio is diagnostic.
- All IC totals are lower bounds. The named S3/S4 solves and lookup work have
  no pinned group-addition conversion on this curve.
- n=41 is one small registered rung. It is a discriminator, not an exponent
  fit, independent replication or deployed-curve claim.
- A proposed direct root for the `R+s1=O` branch was tested and rejected: it
  introduced a root absent from the affine S4 polynomial. That attempt was
  reverted; the final release suite passes.

## Reproduction map

- Frozen protocol: [`PROTOCOL.md`](PROTOCOL.md)
- Ecbench spec: [`SPEC.json`](SPEC.json)
- Full session and comparisons: [`sessions/n41-local-l0`](sessions/n41-local-l0)
- Replay receipt: [`audit-replay-all.json`](sessions/n41-local-l0/audit-replay-all.json)
- Baseline and candidate isolation records: [`records`](records)
- Machine-readable summary: [`evidence.json`](evidence.json)
- Validation and preserved failures: [`VALIDATION.md`](VALIDATION.md)
- PDF handoff: [`REPORT.pdf`](REPORT.pdf)

