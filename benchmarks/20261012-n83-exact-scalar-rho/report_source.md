# Faster rho batching, with the exact scalar candidate retained

Fresh matched follow-up against the tuned degree-83 control · run-20261012T1543Z

12 October 2026

> The selected 256-lane policy lowers process cost against the tuned 64-lane control on this fresh five-target panel: native current/selected ratio 1.069651 [1.043178, 1.089238]. Mean process cost is 2.987008 s versus 2.780545 s. Exact Barrett arithmetic is verified but unselected. The lower interval edge overlaps the observed A/A ratio range, so the finding remains conditional on this host and panel. The combined edition preserves all prior 60 pages.

## 1. What changed, and what the new comparison measures

> Native current-tuned / selected whole-process geometric ratio 1.069651 [1.043178, 1.089238]; mean actual process costs 2.987008 s and 2.780545 s. All five target means are lower under the selected reference-scalar Fold256/DP6/off policy.

The optimized control is the prior selected Fold field, 64 lockstep lanes, DP6, escape disabled and reference u64 scalar arithmetic. The scalar candidates instantiate distinct generic implementations in one frozen executable. After training retains reference arithmetic, both held-out arms share that reference-scalar monomorphization and differ only in lane count.

Training freezes one scalar kernel and lane/DP configuration before five fresh held-out points, with three paired seeds each. Primary whole-process wall includes setup through producer replay. All 139 processes independently verify.

| Recorded native estimand | Current / selected ratio | 95% target-block interval |
| --- | --- | --- |
| Primary whole process | 1.069651 | [1.043178, 1.089238] |
| Secondary online | 1.075203 | [1.050382, 1.094289] |
| Secondary arithmetic mean process | 1.074252 | [1.048202, 1.090360] |

The selected change is a larger lockstep batch: 256 lanes instead of 64, retaining reference scalar arithmetic, Fold field and DP6. Exact Barrett is slightly slower in training and remains unselected. The observed runtime improvement is within the declared grid and host; it does not establish a global optimum or change the generic counting scale. This follow-up runs no new IC campaign and constructs no new IC/rho ratio. Historical matched studies remain intact.

| Native descriptive scope | n | Median | Minimum |
| --- | --- | --- | --- |
| Current pooled process (s) | 15 | 2.833963 | 1.152991 |
| Selected pooled process (s) | 15 | 2.850072 | 0.963513 |
| Paired seed current/selected ratio | 15 | 1.075774 | 0.994348 |

The selected pooled median is slightly higher than current, although all five target means improve. Pooled medians and 15 seed-pair ratios are descriptive estimands; they do not replace the primary five-target estimator. Every seed is retained.

The primary point estimate 1.069651 exceeds the maximum A/A ratio 1.046541, while its lower interval edge 1.043178 lies inside that observed noise range. A/A is descriptive, not a population confidence bound; independent-host confirmation remains open.

*Sources: [N01](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/PROTOCOL.md), [N02](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/examples/n83_rho_speed_eval.rs), [N05](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/records/panel/statistics.json), [N11](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/MEASUREMENT_REVIEW.md). [Native descriptive supplement](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/analysis/descriptive_summary.json).*

## 2. Exact retained curve and fresh disjoint targets

The model is `icv1-f2m83-t6151469093347-9eba6cd5`: y²+xy=x³+x²+1 over GF(2)[z]/(z^83+z^7+z^4+z²+1), subgroup r=8569786107849059, cofactor 1128547018, G=(6614244651398184860292015,8868260410029249391539529). Degree 83 describes the field, while the subgroup has about 52.928 bits. The signed Frobenius quotient has 166 phases, with λ=2998004823950457.

Native C++ generation uses seed 202610126201 and conditionally uniform nonzero scalar proposals. It excludes the signed Frobenius orbits of the original exposed target, all K600 points, all eight prior optimization targets, and other accepted new points. These are sampled without replacement under the declared exclusions. A separate validation sidecar is checker-only; each solver receives a point-only file.

![Frozen stages: A/A10 + arithmetic18 + policy81 + holdout30 = 139 fresh processes. The two selections precede holdout; all losses and failed attempts are retained. [N01, N03–N09]](figures/campaign_flow.png)

Frozen stages: A/A10 + arithmetic18 + policy81 + holdout30 = 139 fresh processes. The two selections precede holdout; all losses and failed attempts are retained. [N01, N03–N09]

The fixed native draw stream accepts eight scalar proposals and rejects none. The 609 prior signed orbits cover the original exposed target, K600's 600 representatives and the eight old optimization targets. All eight new points and their orbit exclusions pass independent replay; validation scalars are confined to the checker sidecar. Full points and hashes are in the generator receipt.

*Sources: [N01](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/PROTOCOL.md), [N08](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/verification/fixtures/generator_receipt.json), [N07](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/verification/final_receipt.json), [N09](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/native_runner/src/main.rs).*

## 3. Exact Barrett reduction and conditional addition

The candidate replaces the reference hot scalar multiplication's 128-bit remainder with an exact precomputed reciprocal. Its domain is 2 ≤ m ≤ 2^64−1 and canonical a,b < m. Let B=2^128, x=ab, μ=⌊B/m⌋ and q=⌊xμ/B⌋, computed by the exact high 128 bits of a 128×128-bit product. The implementation obtains μ without representing B: μ=(B−1) div m plus 1 exactly when (B−1) mod m=m−1.

Written proof: B/m−1 < μ ≤ B/m and x≥0, so x/m−x/B ≤ xμ/B ≤ x/m. Since x/B<1, x/m−1 < xμ/B, hence ⌊x/m⌋−1 ≤ q ≤ ⌊x/m⌋. Thus 0 ≤ x−qm < 2m and a single conditional subtraction returns x mod m, including x=0. Every product, carry and subtraction uses integer arithmetic; no floating-point reciprocal is introduced.

The high product splits each operand into two 64-bit limbs. For p00=a0b0, p01=a0b1, p10=a1b0 and p11=a1b1, the middle carry is (p00>>64)+low64(p01)+low64(p10). The upper half is p11+(p01>>64)+(p10>>64)+(carry>>64). The exact limb expansion and finite native oracle checks have separate receipts; neither finite tests nor a polynomial certificate alone proves all machine-word conditions.

For modular addition, add two canonical u64 values with a carry flag. If carry is set or the wrapped sum is at least m, subtract m with wrapping; otherwise retain the sum. Since a+b < 2m, this yields the exact canonical residue even when the machine addition carries. The candidate delegates its uniform sampler to the same UniformInt<u64>, preserving the reference random stream.

| Verification layer | Actual scope / status |
| --- | --- |
| New worker exactness tests | 7 release tests pass, including high-product boundaries, arbitrary-width modular arithmetic, strict configs, identical sampler streams and exact fixture solve. |
| Independent GMP / integer oracle | 100,049 high-product checks; 1,048,575 multiply and 1,048,575 add checks over 1,023 moduli; 41,936 products take the correction branch. |
| Identity certificate | IDC1hc23ae1d80fc5fefa: limb polynomial over F_(2^61−1), 64 statement-bound + 64 additional points. Finite-field expansion only; proof field null. |
| Lean mathematical declarations | 4 checked Init-only Lean 4.34.1 declarations: limb expansion, actual floor-reciprocal residual bound, one correction, and canonical add. Axiom and fresh-kernel checks pass; Rust carry/casts, MAX reciprocal and wrapping word bridges remain manual/native. |

*Sources: [N02](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/examples/n83_rho_speed_eval.rs), [N10](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/audit/ARITHMETIC_PROOF.md). Written proof is distinct from finite test evidence and formalization status. [Lean proof map](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/audit/PROOF_MAP.md) exposes the remaining machine bridges.*

## 4. A/A noise and scalar-kernel selection

Five identical fresh-process reference pairs on training_1 and seed 202610126211 quantify runtime spread. Their answers and exported counters must match. These ten processes do not select a kernel or policy and do not enter holdout estimates.

![All five A/A process ratios: native median 0.995031, min 0.962508, max 1.046541, range 0.084033. Shading is observed range, not a CI. [N06]](figures/aa_noise.png)

All five A/A process ratios: native median 0.995031, min 0.962508, max 1.046541, range 0.084033. Shading is observed range, not a CI. [N06]

| Kernel / same Fold64-DP6 policy | Verified cases | Mean process s | Mean online ms |
| --- | --- | --- | --- |
| barrett-L64-DP6-E0 | 9/9 | 2.923128 | 2683.602506 |
| reference-L64-DP6-E0 | 9/9 | 2.903055 | 2674.156736 |

All nine arithmetic pairs recover the same answers and agree on exported outcomes, walk counters, policy and charges. This is recorded-semantic agreement, not a complete per-step trace. Exact Barrett mean process cost is 2.923128 s versus reference 2.903055 s, a neutral/slower observation within the measured A/A spread. The native rule retains reference arithmetic; no faster scalar-kernel claim is made.

Internal per-lane capped-walk restarts occur and remain in exported counters and charges. Their work is included in solve time. All completed-process statements refer to no process timeout or total step-cap failure; they do not assert zero internal restarts.

*Sources: [N01](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/PROTOCOL.md), [N03](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/records/panel/arithmetic/selection.json), [N06](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/records/panel/aa/noise_summary.json).*

## 5. Every registered policy, including slower cases

The chosen scalar kernel receives all nine policies: lanes 64/128/256 × distinguished-point bits 6/8/10, with fruitless escape fixed off. Each policy sees the same three training points and three seeds per point. Lane count denotes lockstep work within one OS thread. A wider policy may change a walk trajectory and is judged by its complete process cost.

| Policy | Verified / eligible | Mean process s | Mean online ms |
| --- | --- | --- | --- |
| reference-L128-DP10-E0 | 9/9; eligible | 3.494705 | 3272.328770 |
| reference-L128-DP6-E0 | 9/9; eligible | 2.726115 | 2503.266741 |
| reference-L128-DP8-E0 | 9/9; eligible | 2.749576 | 2521.028409 |
| reference-L256-DP10-E0 | 9/9; eligible | 3.411012 | 3191.230668 |
| reference-L256-DP6-E0 | 9/9; eligible | 2.680834 | 2456.858188 |
| reference-L256-DP8-E0 | 9/9; eligible | 2.712132 | 2477.610694 |
| reference-L64-DP10-E0 | 9/9; eligible | 3.634782 | 3413.376393 |
| reference-L64-DP6-E0 | 9/9; eligible | 2.891766 | 2665.835796 |
| reference-L64-DP8-E0 | 9/9; eligible | 2.931672 | 2706.288331 |

![Native registered training scores; teal marks the immutable native selection. All nine policies remain visible. This grid is not a global optimum. [N04]](figures/policy_scores.png)

Native registered training scores; teal marks the immutable native selection. All nine policies remain visible. This grid is not a global optimum. [N04]

Frozen native selection: reference-L256-DP6-E0, mean training process 2.680834 s; 64/DP6 scores 2.891766 s and 128/DP6 2.726115 s. All nine policies are eligible and all 81 processes verify within time and total-step limits. Selection precedes holdout. The default 64/DP6 is preferred only within 2% of the minimum.

*Sources: [N01](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/PROTOCOL.md), [N04](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/records/panel/policy/selection.json), [N09](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/native_runner/src/main.rs).*

## 6. Fresh held-out targets and native uncertainty

The primary estimate first takes each target's three-seed arithmetic process mean for each arm, then computes the geometric mean of the five current/selected ratios. Native SplitMix64 target-block bootstrap uses 20,000 paired resamples, seed 202610126299, with percentile 95% intervals. Three seeds nested in five targets do not become fifteen independent target samples.

![All 30 fresh held-out process costs and the native per-target three-seed means. Identical point/seed pairs alternate execution order. No fastest seed or target is removed. [N05]](figures/holdout_cases.png)

All 30 fresh held-out process costs and the native per-target three-seed means. Identical point/seed pairs alternate execution order. No fastest seed or target is removed. [N05]

| Target | Current seed costs (s) | Selected seed costs (s) | Native means (s) / ratio |
| --- | --- | --- | --- |
| Q1 | 1.493073<br/>2.833963<br/>2.194318 | 1.417215<br/>2.850072<br/>2.118903 | 2.173785 / 2.128730<br/>ratio 1.021165 |
| Q2 | 4.638999<br/>4.591936<br/>3.734425 | 3.874468<br/>4.221322<br/>3.716998 | 4.321786 / 3.937596<br/>ratio 1.097570 |
| Q3 | 3.528888<br/>2.154345<br/>4.033394 | 3.225827<br/>1.888357<br/>3.908533 | 3.238876 / 3.007573<br/>ratio 1.076907 |
| Q4 | 1.152991<br/>2.290008<br/>1.501968 | 0.963513<br/>2.126395<br/>1.499933 | 1.648322 / 1.529947<br/>ratio 1.077373 |
| Q5 | 2.729404<br/>3.356201<br/>4.571203 | 2.563237<br/>3.119802<br/>4.213605 | 3.552269 / 3.298881<br/>ratio 1.076810 |

All 30 cases independently verify, with no process timeouts or wrong answers. Seed triples are 202610126241–243, 244–246, 247–249, 250–252 and 253–255 for Q1–Q5. All individual process times remain above; native means use three paired seeds per target. Full rows retain online, CPU, RSS and counters.

*Sources: [N05](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/records/panel/statistics.json), [N07](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/verification/final_receipt.json), [N09](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/native_runner/src/main.rs), [N11](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/MEASUREMENT_REVIEW.md).*

## 7. Timing boundary, host and supported interpretation

![Native current/selected ratio direction. Values above 1 favor selected. Primary geometric and secondary arithmetic ratios preserve distinct estimands. The renderer does not estimate intervals. [N05]](figures/holdout_intervals.png)

Native current/selected ratio direction. Values above 1 favor selected. Primary geometric and secondary arithmetic ratios preserve distinct estimands. The renderer does not estimate intervals. [N05]

Whole-process wall includes launch and output, point parsing, normal-basis conversion and target-independent setup, target-dependent walk, collision recovery and producer scalar replay. Online timing is secondary and excludes the recorded setup boundary. Independent C++ answer replay is an external correctness gate; its own recorded cost is reported separately.

| Envelope / admission field | Actual value |
| --- | --- |
| Host and resources | AMD EPYC 9V74 model; KVM guest with 5 logical CPUs, cgroup quota 4 CPUs and memory 16 GiB. Rust 1.99.0; Cargo default release. |
| Pin / thread / address-space / timeout | CPU3; one OS thread; 12 GiB RLIMIT_AS; 60 s per rho child. |
| Monitor and reservation | Native 5 ms nominal wait4 poll plus scheduling overhead; all four campaign reservations and successful external replay uncontended. |
| Recorded campaign costs | 423.474951 s actual child wall; 422.987101 s child CPU across all 139 runs. Excludes builds, replay and parent orchestration. |
| Independent replay cost | 0.886574336 s wall; 0.833018 s user + 0.051324 s system; peak RSS 109,188 KiB. Outside solver intervals. |
| Unmet scientific gates | No new calibrated primitive-cost S; no L2 admission; no independent-host replay or observable physical-host exclusivity. |

The complete uncontended panel supports a conditional host runtime improvement for the selected policy: 1.069651 [1.043178, 1.089238]. The lower interval edge overlaps the descriptive A/A range, and only five targets are sampled under exclusions. Secondary mean online costs are current 2755.666121 ms versus selected 2551.759882 ms. No generic boundary, exponent or historical IC/rho ratio changes.

*Sources: [N01](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/PROTOCOL.md), [N05](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/records/panel/statistics.json), [N06](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/records/panel/aa/noise_summary.json), [N11](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/MEASUREMENT_REVIEW.md), [N12](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/provenance/host_context.json).*

## 8. Validation, requirement status and reproduction

| Requested obligation | Status and evidence / remaining gap |
| --- | --- |
| Faster optimized-control comparison | Verified on the declared host/panel: 139 runs retained; matched current tuned reference, fresh five-target holdout; conditional 1.069651 process ratio. |
| Exact scalar arithmetic | Implementation and finite verification complete; four Lean mathematical declarations checked. Rust machine-word bridges remain partial. Candidate is not selected. |
| Fresh targets and independent replay | Verified: 8 fresh points and exclusions, all 139 answers, 9 arithmetic semantic pairs; 0 invalid or globally capped solver processes. |
| Full required implementation checks | Full release library: 3,992 pass, 0 fail, 112 ignored; new worker: 7 pass; native runner: 22 pass. Initial compile and replay-wrapper failures preserved. |
| Comprehensive PDF and visuals | New source-linked addendum, combined edition and SVG/150 dpi PNG records; document-layout and full visual QA recorded beside report. |
| Calibrated S / independent-host admission | Partial: missing accounting and host admission bridges remain explicit; no ledger promotion. |

Sources and raw records are authoritative; benchmarks holds identical browsing copies. This addendum preserves the broad basis/invariant synthesis and the two prior matched degree-83 studies as a frozen 60-page appendix.

| Source | Evidence pointer |
| --- | --- |
| N01 | [Frozen prospective protocol](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/PROTOCOL.md) |
| N02 | [New exact-scalar worker](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/examples/n83_rho_speed_eval.rs) |
| N03 | [Arithmetic selection and all pairs](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/records/panel/arithmetic/selection.json) |
| N04 | [All policy scores and frozen selection](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/records/panel/policy/selection.json) |
| N05 | [Native held-out statistics](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/records/panel/statistics.json) |
| N06 | [Native A/A noise summary](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/records/panel/aa/noise_summary.json) |
| N07 | [Independent all-answer replay](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/verification/final_receipt.json) |
| N08 | [Fixture proposal/exclusion receipt](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/verification/fixtures/generator_receipt.json) |
| N09 | [Native runner and statistics source](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/native_runner/src/main.rs) |
| N10 | [Arithmetic audit and proof status](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/audit/ARITHMETIC_PROOF.md) |
| N11 | [Root measurement review](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/MEASUREMENT_REVIEW.md) |
| N12 | [Host and frozen provenance](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/provenance/host_context.json) |

Pre-launch freeze commit adddc3555d3355f5d5d83fbf45cc6793cb6ada8f. The [native file audit](https://github.com/aburan28/crypto/blob/codex/n83-rho-speed-20261012/research/n83_speed_20261012/verification/file_provenance_receipt.json) passes: 139 execution/producer/observation/environment records, 2,439 recorded hash checks, 455 distinct bound files and 1,036 snapshots. Frozen provenance and manifest bind the sources and receipts; no historical denominator is reused.

*Conductor was attempted before edits and is unavailable. No chat transcript or credential content is included in the evidence.*
