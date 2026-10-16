# Four-byte fixed30 root slots reduce memory with a measured online cost

The exact `2^30` four-byte root table completes K2940, K3200, and K3300
with direct 21-byte states. At K2940 it uses 11.38 GiB peak RSS, compared
with 12.72 GiB for the five-byte/root22 control on the same base and 1,000
points. Its mean target-online interval is 75.155 ms versus 44.000 ms for
that control. The frozen six-arm selection rule chooses the existing
K3300/root11 five-byte configuration at 33.475 ms. Its canonical identity
for this exact source snapshot is
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h0102baeb9f1a`.

## Same-allocation comparison

All arms used the N83 curve `EC1N83Ckb1a1haacb69ee9f37` over
polynomial-basis `GF(2^83)`, defining polynomial bitmask
`9671406556917033397649557`, model `y² + xy = x³ + x² + 1`, and prime
subgroup order `8569786107849059`. Each arm had a fresh factor-log/index
process, rank seed 29, four rank threads, one target thread, a 1,200-second
cap and 14-GiB observed RSS guard. The public 1,000-point file SHA-256 is
`35e34151dcc1ad6048939560adf28fc7c1b1017833fe005a656c3f4e2bc07de8`;
the independently checked answer sidecar SHA-256 is
`f8cbc90f7f592aae5014ca4f2693a3327a72e935970f22e8f02868e4aeed7f71`.
Online time includes target query, decomposition, relation check, descent,
and recovery check; target-independent setup is reported separately.

| K / layout | IC1 hash suffix | Usable points | Mean online ms | 95% target-block interval ms | Setup s | Peak GiB |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 2940 / root4f30 | `f3eb1c606e67` | 488,040 | 75.155 | 70.535–80.027 | 569.2 | 11.38 |
| 2940 / root22 | `9de00ba2d51b` | 488,040 | 44.000 | 41.298–46.844 | 480.3 | 12.72 |
| 3200 / root4f30 | `6679391fb911` | 531,200 | 91.463 | 85.880–97.200 | 699.4 | 12.71 |
| 3200 / root21 | `15c99a7b8d6f` | 531,200 | 35.336 | 33.185–37.541 | 565.0 | 13.71 |
| 3300 / root4f30 | `97c997eb096f` | 547,800 | 121.305 | 114.376–128.432 | 782.5 | 13.26 |
| 3300 / root11 | `0102baeb9f1a` | 547,800 | **33.475** | 31.576–35.441 | 580.8 | 10.05 |

At every same-K pair, the complete base digest, exact root-key count, and
all 1,000 target probe counts agree. All six arms passed full-base, subgroup,
four-point relation, scalar, and exclusive online-phase replay. The new
four-byte direct-state mode offers a measured lower-memory option at K2940
and K3200; its denser table increases setup and target-online time on this
allocation. Each full IC1 candidate is in the adjacent `candidate_K*.json`
manifest, and `measurements.json` retains candidate/workload/run IDs and
every target interval and phase cost.

## Receipts

[Corrected cloud run 38608258672](https://github.com/aburan28/crypto/actions/runs/38608258672)
executed commit `aa2f683db4e8242915e4025f483c621664ff8e33` with IC source
SHA-256 `6b22a065d95d5481342c1b1cc37dc9479200dd4f0832bcd53ae6ce4b8c9312d5`.
Artifact ID `11922291732` has SHA-256
`8bad6227eb71742d1b0710ed610f6e7b68b6ac20776f83ca5618dab830b03b57`.
The cloud verifier and independent local full-base replay produced the same
analysis SHA-256
`e622935d6a24de7625715e46b8d496f4fbff67c1f29e8a090ecf1792e95adb6f`.
The first source execution and its insertion failure remain in
`FIRST_RUN.md`, with raw artifact SHA-256
`bd270ee4609c503c99b07c43a186dd46a661219af53fd4c0adc45fddf7319c2f`.
The online-optimal base from the larger 14-arm refinement is undergoing a
separate fresh-point Tier 1 same-point test in
[PR #1855](https://github.com/aburan28/crypto/pull/1855).
