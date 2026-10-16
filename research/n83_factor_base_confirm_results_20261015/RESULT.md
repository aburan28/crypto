# N83 factor-base confirmation on 4,000 fresh points

K3300/root11 kept its online lead over K3100/root22 on 4,000 new,
disjoint subgroup points. With both exact bases independently replayed,
K3300 averaged **37.035 ms** per target and K3100 averaged **42.160 ms**.
The paired mean difference, K3100 minus K3300, was **5.126 ms**; the
predeclared 95% contiguous-block bootstrap interval was **[3.386, 6.916]
ms**. Both arms completed all 4,000 targets under the 14 GiB RSS guard.

## Frozen comparison

| Arm | Usable base points | Folded columns | Mean target online ms | Base/index/log setup s | Peak RSS GiB | Status |
| :--- | ---: | ---: | ---: | ---: | ---: | :--- |
| K3300/root11 | 547,800 | 3,300 | 37.035 | 646.188 | 10.05 | `PASS_INDEPENDENT_REPLAY` |
| K3100/root22 | 514,600 | 3,100 | 42.160 | 605.173 | 13.56 | `PASS_INDEPENDENT_REPLAY` |

The comparison holds the exact curve, 4,000 public points, point order,
source kernel, cloud allocation, 14 GiB memory guard, and online timing
boundary fixed. Each arm starts a fresh process with an empty base/index/log
cache; target-independent preparation is recorded separately. The timed
target interval contains the target query, four-summand lookup, relation
check, descent and recovery check. This is a factor-base **stage comparison**
on one cloud instance, recorded as `cloud_instance_exploratory`; the
five-pair same-point Tier 1 rho comparison is in the separate selected-base
result.

The curve is `EC1N83Ckb1a1haacb69ee9f37`: polynomial-basis `GF(2^83)`
with modulus bitmask `9671406556917033397649557`, model
`y²+xy=x³+x²+1`, subgroup order `8569786107849059` and cofactor
`1128547018`. The candidate IDs are:

- K3300/root11: `IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h0f8ab8bfcff7`.
- K3100/root22: `IC1N83Ckb1a1fb514600PDP4rootRCguidedLAgaussTDdirectISO0haa5b1411e6a0`.

The 4,000 public points SHA-256 is
`451df8524b2c356dd4a523260620c4dff406975a122762b98ddddb7987b75c2c`;
the answer sidecar SHA-256 is
`8e18c9c97ba9f57e5eb43f1b39b9e9575b5cfab87840c25c428dac150c868101`.
The workload ID is `12e14f274f8c`. The statistic was frozen before this
run: 80 contiguous blocks of 50 targets, 10,000 bootstrap resamples with
seed `2026101599`, and a positive lower endpoint as the K3300 confirmation
gate. The observed interval meets that gate.

## Receipt and reproducibility

The executed source commit is `e65f35916dfe4791a60b938b8ab0bdc875d0ae7b`;
the exact IC source SHA-256 is
`5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e`.
[GitHub Actions run 38663668487](https://github.com/aburan28/crypto/actions/runs/38663668487)
completed on Linux x86-64 with an exposed Intel Xeon Platinum 8573C.
Raw artifact `11942386005` has SHA-256
`b2f46a1b15138a68f8c504f310287ef165515c5260c50800a3ec3bf7831383a3`;
its independent `analysis.json` has SHA-256
`4e2981b86236586a1beeb7f3fe5ce7eacff9e97635e1ac3438daea1d80932699`.
The verifier replayed both full bases and every target scalar and four-point
relation before computing the paired interval. The source-pinned sealer
checked the raw ZIP and analysis digests and wrote
[confirmation.json](confirmation.json), both canonical candidate manifests,
and the lossless compressed [workload](workload.json.gz) and
[measurement rows](measurements.json.gz). [records_index.json](records_index.json)
gives the compressed and decoded SHA-256 values.
