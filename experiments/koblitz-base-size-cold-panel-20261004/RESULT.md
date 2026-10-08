# Held-out result: shrinking the compact-orbit base raised rank cost

The preregistered source-curve compact S3-root-index panel recovered both
held-out public points at the baseline and pilot-selected smaller base. All
**24 of 24** held-out IC/rho pairs reached full rank, recovered the same
scalar, and passed independent general-curve replay of every base log, rank
row, target relation, scalar and charged phase sum. There were no held-out
timeouts, OOMs, nonzero exits or rejected replays. The n41 point was
`[1071506060992,898053054019]`; the n53 point was
`[7960849849661793,7443722527872608]`. Each arm got only that public
point, not its known-answer scalar. These are **exploratory, unisolated
macOS example-binary timings**, not an `ecbench` result, controlled wall-time
speedup or ECC2K-130 forecast.

The primary one-target online comparison starts after reusable IC base,
index and rank-log setup and includes target query, PDP, relation check,
descent and recovery check. Rho online runs from target-dependent walking
through scalar validation on the same point. The medians below are over six
fresh-process pairs per `K`; the ratio column is the median of the six
**paired** rho/IC online ratios, not a ratio of medians.

| `n` | `K` / actual base points | Verified pairs | IC online median (range), ms | Rho online median (range), ms | Paired rho/IC online median (range) |
| ---: | ---: | ---: | --- | --- | --- |
| 41 | **85 / 6,970** | 6/6 | **2.631** (2.564–2.684) | 13.409 (13.234–13.484) | **5.106** (5.021–5.182) |
| 41 | 64 / 5,248 | 6/6 | 3.320 (3.174–3.417) | 13.913 (13.250–14.551) | 4.255 (3.878–4.522) |
| 53 | **220 / 23,320** | 6/6 | **6.736** (4.748–12.033) | 135.913 (115.030–199.987) | **19.300** (13.540–28.020) |
| 53 | 160 / 16,960 | 6/6 | 11.767 (9.037–29.246) | 155.319 (127.636–621.984) | 13.899 (4.364–61.178) |

The supplementary fully charged in-process cold interval includes factor-base
construction, S3 preparation and root-index build, every rank query and PDP
attempt, relation checks, matrix work, final relation LA and the target
interval. Rho cold includes its jump setup and verified solve. Public-point
fixture generation and process launch are outside both intervals; outer
`/usr/bin/time` wall/CPU receipts are archived separately. The independent
replayer checked the producer's 15 exclusive IC cold phases, the five target
phases and rho's cold/fixture accounting for every accepted pair.

| `n` | `K` | IC cold median (range), ms | Rho cold median (range), ms | Paired IC/rho cold median (range) | Within-repeat smaller/baseline IC cold |
| ---: | ---: | --- | --- | --- | --- |
| 41 | **85** | **480.211** (475.032–491.279) | 13.659 (13.481–13.805) | 35.144 (34.709–35.982) | reference |
| 41 | 64 | 572.210 (563.420–577.151) | 14.194 (13.498–14.800) | 40.049 (38.737–42.310) | **1.186** (1.167–1.207) |
| 53 | **220** | **5,586.337** (5,265.992–7,783.739) | 136.538 (115.381–200.482) | 40.365 (37.703–48.625) | reference |
| 53 | 160 | 9,418.334 (8,736.181–9,776.486) | 155.784 (127.987–622.345) | 57.965 (15.709–75.415) | **1.636** (1.150–1.833) |

Every one of the **12 paired smaller/baseline cold IC ratios exceeded one**.
The n53 rho and wall-time ranges show substantial host noise; these numbers
do not satisfy the repository's CPU isolation gate. The conclusion is limited
to this frozen producer, K grid and two held-out points. The pilot's fixed
point selected `K=64` and `K=160` by the published lowest-median-eligible-
smaller-K rule before either held-out IC run; the pilot itself had already
found the baseline faster at both sizes. The complete pilot decision and its
three n41 `K=20` target failures remain in [`PILOT_RESULT.md`](PILOT_RESULT.md).

The stage and count evidence explains why a smaller index did not help:

| `n` | Baseline → smaller `K` | Median index build, ms | Median matrix work, ms | Median rank PDP, ms | Rank probes per relation | Total rank probes |
| ---: | --- | --- | --- | --- | ---: | ---: |
| 41 | 85 → 64 | 77.1 → 42.1 | 46.2 → 30.4 | **348.8 → 491.1** | 24,709 → 46,775 | 2,100,243 → 2,993,628 |
| 53 | 220 → 160 | 1,077.3 → 658.9 | 280.5 → 174.7 | **4,196.2 → 8,492.8** | 83,698 → 201,971 | 18,413,613 → 32,315,380 |

Both held-out bases had zero failed rank queries. Reducing `K` saved index
entries and matrix columns, but raised deterministic rank probes enough to
more than consume those savings. This closes **simple K reduction** as the
next improvement to this fixed producer. The evidence-ranked next test is a
probe-reducing rank/query or index architecture at **fixed actual base size**,
where target coverage is held constant; the separate equal-useful-size
source/descendant-native/transported/pullback PDP comparison remains needed
before any ECC2K-130 transfer claim. Also, the pilot exposed a producer
contract flaw: a target failure can exit 0; the independent replay correctly
rejected all three such rows. That exit/status behavior needs a focused
follow-up fix without rewriting the frozen evidence.

[`HOLDOUT_ANALYSIS.json`](HOLDOUT_ANALYSIS.json) is recomputed by the native
[`koblitz_base_size_holdout`](../../examples/koblitz_base_size_holdout.rs)
analyzer from all 24 raw run directories and the published pilot selection.
It checks the exact point and seeds, actual base size and digest, rank,
successful target count, every independent replay receipt, scalar agreement,
resource gate and phase-reconciled times. `HOLDOUT_SHA256SUMS` pins all 240
raw files; `HOLDOUT_ANALYSIS_SHA256SUMS` pins analyzer source, binary and
output. The 30 pilot pairs, including failures, have their own corresponding
manifests. All hashes and format/lint checks passed. The exact release
producer build, toolchain and binary hashes are in [`BUILD.md`](BUILD.md).
