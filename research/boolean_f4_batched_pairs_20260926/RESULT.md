# Boolean F4 batched pair updates: first paired result

The opt-in batch matched the original basis fingerprints, matrix shapes, all critical/field pair counts, both skip counters, divisor tests and word XORs on every seven-case frozen and holdout call. It cleared the preregistered one-thread stage and complete-call gates. These are internal Boolean F4 solver-stage timings, not a one-target IC or DLP speedup.

CI run [36299824843](https://github.com/aburan28/crypto/actions/runs/36299824843) at PR head `3500d80c` passed the reference and candidate F4 tests and completed 66 calls. The [complete raw receipt](runs/36299824843/ci-result.json), SHA-256 `c24d60bc335fbe37c965c34af68bc6d884291c7d29eda3214f30a70515398b25`, includes all A/A, A/B, warmup, full output and host/load records. The Linux x86-64 runner was AMD EPYC 7763, Rust 1.98.1, with one pinned CPU.

| `n20_m30` workload | Reference other | Batched other | Paired other ratio, 95% interval | Reference full F4 | Batched full F4 | Paired full ratio, interval |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Frozen seed | 116.8 ms | 85.9 ms | 1.361, 1.352–1.437 | 593.6 ms | 560.6 ms | 1.061, 1.050–1.064 |
| Holdout `1ac0ffee` | 115.8 ms | 85.4 ms | 1.357, 1.347–1.378 | 684.9 ms | 656.6 ms | 1.044, 1.038–1.051 |
| Holdout `2468ace0` | 118.0 ms | 85.3 ms | 1.382, 1.357–1.387 | 687.3 ms | 654.7 ms | 1.050, 1.047–1.053 |

Milliseconds are each arm's five A/B median; paired ratio medians need not equal the ratio of marginal medians. Frozen A/A ratios ranged 0.990–1.026 for `other_ms` and 0.998–1.010 for complete F4. Build and elimination ratios were about one, as expected. No smaller complete-F4 cell had a median regression beyond its own A/A range. The batch cuts repeated scans in `State::insert`; step-degree selection still scans the pair list. Default promotion awaits a four-thread regression control and final default-on replay.
