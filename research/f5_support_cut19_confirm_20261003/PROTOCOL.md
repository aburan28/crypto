# Fresh-seed confirmation of fixed cut 19

## Hypothesis and prior evidence

The exact support-block full-rank certificate with inner cut 19 and row-basis table width 8 reduces the one-thread n24 degree-4 Boolean matrix-F5 **complete-call** wall time by more than 2× against selective echelon on an isolated physical x86-64 Linux runner. A separate two-thread result should also exceed 2×. This is a solver-stage claim, not an IC online or rho speedup.

The first four-seed Linux run is calibration only: its primary call measured 2.96–2.98× at one thread and 2.49–2.56× at two, but the original overall gate failed on wall-time controls for smaller cases. Preserve that failure and all receipts in `../f5_support_cut19_call_20261002/LINUX_RESULT.md`. The width-six variant also failed its overall gate; its receipt is `../f5_support_rank_bits_20261002/LINUX_RESULT.md`. Neither old run can satisfy this confirmation gate.

## Frozen new inputs and arms

Use four **previously unmeasured** seed blocks. Each seed is the first 16 hexadecimal digits of SHA-256 over the exact ASCII string `F5-cut19-confirm-20261003/<name>`:

| Name | 64-bit seed |
| --- | --- |
| `confirm_a` | `c0d41ab38b872c0d` |
| `confirm_b` | `ae05fd6376bf7a39` |
| `confirm_c` | `467a3257c4d78bba` |
| `confirm_d` | `ff83e4df6657b12f` |

The same native `f4_f2_bench` executable runs both arms in separate processes. The reference is selective echelon (`KIC_F5_ECHELON=2`, row-basis rank disabled); the candidate is support-separated original rows (`KIC_F5_ECHELON=5`, fixed cut 19, row-basis rank enabled, width cap 8). All other benchmark inputs and settings follow `f5_support_cut19_pair.rs`. It records `f5_n24_m24_d4` plus six smaller cases. The candidate must actually return original rows through the support certificate on the primary case. Both arms must match rank, canonical row space, criterion, row and column counts, and the smaller cases' complete non-timing output.

## Resource and timing contract

Rebuild on a physical Ubuntu x86-64 runner with AVX2 and BMI2. Reserve and pin one CPU for the primary measurement, then use a separate two-CPU run. For each seed and thread count, execute one warmup per arm, five reference/reference A/A pairs and five alternating reference/candidate pairs. The time is the whole matrix-F5 call, including build, reduction and unpack. The runner saves phase times and counted reduction XORs. For a contended, incomplete or failed attempt, keep its receipt and retry at most twice. Analyze the **first** clean, complete attempt with 22 successful calls. Preserve source and binary hashes, hardware, isolation and all attempted calls.

For each of the four fresh seeds at **both** thread counts, require the paired complete-call median **and** exact five-pair bootstrap 95% lower bound above 2.00×. Require candidate counted reduction word XORs at most 35% of reference, and all exactness/route checks. A missing qualified attempt, timeout, output mismatch, fallback or hardware mismatch fails the gate. There is no pooling across seeds or thread counts.

## Inactive smaller-case controls

The six smaller cases take the same elimination branch in both arms. The support certificate only activates for degree 4, 24 variables, at least 4096 rows and 8192 columns; it is inactive on the controls. `KIC_GF2_ROW_BASIS_RANK` and its width are read only in that certificate. Candidate `KIC_GF2_TABLES=8` clamps to 4, the same effective setting as the reference's 4. The runner must verify that all 22 calls on every smaller case have no certificate route and have identical non-timing output; this is a functional path control. Preserve their paired wall ratios and A/A ranges as diagnostics, with no wall threshold for an identical computation.

The previous gate's smaller-case rule compared the paired median with the *minimum of five independent A/A ratios*. That minimum is not a confidence bound: equal code paths can fail it through ordinary sampling noise. The old failed result remains failed. This independent experiment uses the stated structural control and fresh seeds; its acceptance does not reinterpret the old data.

The workflow is `.github/workflows/f5-support-cut19-confirm.yml`. The native gate is `f5_support_cut19_confirm_analyze.rs`; `run_linux_confirm.sh` launches and preserves attempts. Stop after this one four-seed, two-thread-count panel and report every failed cell. Do not choose an alternate seed or attempt after observing its wall ratio.
