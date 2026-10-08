# Same-target online F6-IC versus rho on the prepared n17 curve

**Observed result: rho is faster on both frozen public points.** Every
one-target solve completed, independently replayed its scalar, and had
exclusive online phases summing exactly to wall. The default headline
ratio `rho_online_ms / IC_online_ms` was 0.0219–0.0221 on T1 and
0.000515–0.000530 on T7 across two repetitions. These are paired
observations on two points, not an estimate for arbitrary targets or
n83. A 2× reduction of the measured F6-IC online time alone would not
reverse either point's ordering; that is an inference from the paired
times, not a measured new candidate.

The [protocol](PROTOCOL.md) was committed as `fa5f814aa` before the
new runs, the [scripts](run.sh) as `cda7d21ca`, and the exact
[candidate/reference/input/workload freeze](FREEZE.tsv) as `fe292cecd`
before timing. The [raw runs](runs/), [measurement rows](measurements.jsonl),
[paired rows](pairs.jsonl), [rho reference](rho_reference.json), and
[derivation check](DERIVATION_CHECK.json) are retained. The binary is
the source-pinned compact-refutation pilot worker, with its experimental
flag **off** for IC; its SHA-256 and source hashes are in
[BUILD_IDENTITY.txt](BUILD_IDENTITY.txt) and
[SOURCE_SHA256SUMS](SOURCE_SHA256SUMS).

IC candidate:
`IC1N17Ckb1fb62PDP3f6RCsampleLAgaussTDpdpISO0h9628dfd41b76`.
Curve: `EC1N17Ckb1hbbe2b5b6b1e6`, field size `2^17`, subgroup order
65,587, 62 actual usable factor-base points and 29 folded relation
columns. Its prepared certified log table and symbolic setup were ready
before timing. Rho used the **same supplied point**, signed Frobenius
and negation quotient, one walk, seed 20261004039, 16 jumps, at most
64 restarts and 65,536 iterations per restart. Both used one Rayon
thread, the same worker binary and a 180-second process cap. No
cross-target state or amortisation was used.

| Public target | Workload ID | Rep | F6-IC online ms | Rho online ms | Rho/IC | Recovered scalar | Verified |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| T1 | `146a1e9ee3c8` | 1 | 2.808666 | 0.061541 | 0.021911 | 4785 | Both |
| T1 | `146a1e9ee3c8` | 2 | 2.773792 | 0.061292 | 0.022097 | 4785 | Both |
| T7 | `ced1677f0976` | 1 | 140.192125 | 0.074292 | 0.000530 | 2391 | Both |
| T7 | `ced1677f0976` | 2 | 139.955459 | 0.072083 | 0.000515 | 2391 | Both |

The F6-IC online clock begins after reusable log/index preparation and
includes query generation, every target PDP attempt, relation check,
descent and scalar replay. T1 used one PDP attempt and 2.750–2.775 ms
in target PDP; T7 used eleven attempts and 139.895–140.134 ms in PDP.
T7 F4 build was 116.106–116.307 ms and charged 17,167,040 word
operations. The rho online clock begins after reusable packed-curve and
Frobenius preparation, **before target lift**, and includes the
target-dependent jump/walk setup, walk and scalar replay. T1 used 20
walk iterations and 19 charged walk additions; T7 used 46 iterations
and 45 additions. Every rho run succeeded in one restart. The actual
field dispatch was `pmull` on Darwin arm64.

The rho policy's recent-point cache has 4,096 direct-mapped slots
(`RECENT_SLOTS` in the pinned source). For subgroup order 65,587 and
automorphism order 34, its `rho_trail_mask` is zero, so every visited
point is eligible for the stored-point table; the exact table capacity
and peak RSS were not measured. Factor-base construction, certified
log preparation and final relation-matrix work were already reusable
and are excluded from the primary one-target online interval. Their
cold total is unknown in this study, not zero. Natural ordinary-query
yield and a complete n83 relation are also unmeasured here.

The host was **unisolated**; its physical CPU model and host-wide CPU
partition were unavailable, so these CPU ratios are exploratory under
the repository isolation gate. The two repetitions give observed ranges,
not confidence intervals. They show the gap on these two same-target
comparisons and rule out claiming an end-to-end IC win from the recent
n17 F6 stage pilots.

Rebuild the rows and checks with
`sh research/f6_ic_rho_online_gap_20261006/derive.sh`.
