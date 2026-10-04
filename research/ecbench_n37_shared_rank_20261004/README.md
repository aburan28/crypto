# n37 shared-rank IC versus strong rho on one public point

Decision: **verified one-target correctness and bounded accounting; online
speedup remains unset**. This is an accounting/correctness round, not an
algorithmic advance. The [preregistered protocol](../notes/ecc2k130/n37_shared_rank_ecbench_20261004/PROTOCOL.md)
and fixed [spec](../notes/ecc2k130/n37_shared_rank_ecbench_20261004/SPEC.json)
were committed in [PR #1313](https://github.com/aburan28/crypto/pull/1313)
before the public target was planned or solved. The native input check found
that Q = `(0x100e90da4, 0x124043682f)` is in the specified subgroup and
outside all 13,459 archived signed-Frobenius exclusion orbits. Its input law
is `hash_to_subgroup_v1`; neither solver receives a constructed logarithm.

The exact candidate is
`IC1N37Ckb0fb3108PDP3mitmfrobeniuscountedRCsampleLAgaussTDpdpISO0ha48bfcf857fa`.
The [final diagnostic claim](CLAIM_FINAL_DIAGNOSTIC.json) embeds the complete
canonical candidate and workload records and their SHA-256 hashes. The
candidate has **3,108 actual usable points**, 42 folded columns, no isogeny
transport, a target-blind full-rank relation log and a point-only m3 target
descent. The `ecbench` workload is `W157bda94ea05`; the canonical claim
workload ID is `bad0afc30d3a`, whose record additionally binds the run's
binary and resource envelope. The claim's `EC1N37Ckb0h3e63deb3fc2b` is the
IC1 manifest's exact-record ID; the curve registry also calls this exact
polynomial, coefficients, subgroup and generator
`EC1N37Ce0h0c51aa4aa7c3`. Both aliases and the complete point/field records
are retained rather than inferred from `n=37` alone.

| Same Q, five measured rounds | Verified | Cold counted S, lower bound | S / generic floor, lower bound | S / rho, bounded | Median target online, L0 diagnostic |
| --- | ---: | ---: | ---: | ---: | ---: |
| Strong signed-Frobenius rho, 32 lanes | 5/5 | 0.422946 | 2.903 | 1.000 | 0.743750 ms |
| Shared-rank IC | 5/5 | 5.748205 | 39.454 | 13.591 | 0.061917 ms |
| IC A/A control | 5/5 | 5.748205 | 39.454 | 13.591 | 0.060542 ms |

The derived generic floor is `sqrt(pi/(2*74)) = 0.145695` in S. The
five-pair **cold counted** IC/rho ratio is 13.591, with a within-one-target
bootstrap interval of 10.553–18.719; both arms leave native work unpriced,
so this ratio compares lower bounds and is not a full-cost speedup or a
cross-target uncertainty interval. IC's counted total is 87,290.074436 GAE,
of which 87,158.030760 is reusable target-blind setup and 132.043676 is
target work. This `ecbench` pinned-only calibration is different from the
default calibration in PR #1309; its GAE figure must not be directly
subtracted from #1309's 194,826 GAE setup figure. The rank is 42/42 after
55 target-blind trials in every run, with the frozen base digest
`8460ac4c28515db701c3897a03b4ce0f28abf7cd56fd759ad095f98436a76dcf`.
All five target solves hit the direct m3 lookup. The shifted residual path
remains untested by this Q.

The five exclusive IC online phases sum exactly to each target interval.
Median phase times are 22.377 µs for query/validation, 23.083 µs for PDP,
1.167 µs for witness checking, 0.292 µs for log combination and 15.084 µs
for full-point scalar replay. The **paired median raw rho/IC online ratio**
is 11.739, ranging from 3.599 to 12.271 across the five walks. This is
descriptive only: macOS cannot earn the required L2 CPU isolation, the IC
A/A paired ratio ranged from 0.833 to 1.287, and one point gives no
target-population interval. The frozen [native analysis](RESULT_FINAL.json)
therefore sets `online_speedup: null`. `ecbench compare`'s wall field measures
the whole cold solve; it must not be substituted for this primary online
interval. Its [cold table](TABLE_FINAL.txt) and [comparison](sessions/mac_arm64_l0_02/comparisons/rho-strong__ic-shared.json)
remain available separately.

The final session `ECBS1h33c3e5e4f75d` contains 18/18 verified executions,
including warm-ups, and the [local audit](AUDIT_FINAL.json) replays all 15
measured records identically (`425a3de7b1c8c6d6f518bacb648d67aad8575aa513174d62a4573b98254d9fef`).
The earlier session `ECBS1h774f0e722ea2` is preserved because it ran before
the metadata-only change that embedded manifest records; the final binary
also [replayed it identically](AUDIT_FIRST_FINAL_BINARY.json). The
[claim check](CHECK_FINAL.json) correctly refuses promotion without an
independent other-host certificate. A separate
[Linux x86-64 replay](independent_validation_20261004/RESULT.md) now
reproduces all 15 measured records under another host class, and its
attached diagnostic report passes the claim schema. The final session is
still L0 in every arm, and its native field arithmetic, hashing, allocation
and modular combination remain incompletely priced. The claim file's numeric
raw ratio is explicitly labelled descriptive by its verdict and is not an
accepted speedup; the frozen aggregate speedup remains null.

The next decision is to run this exact candidate and strong rho on a Linux
host that actually earns L2 and price native work in a common unit. Freeze a
separate panel of at least eight new one-target workloads before
generalizing to targets.
The current point says nothing about n41/n53 scaling, the m83 confidence
gate, degree-263 descent, or ECC2K-130 at m131. A separately frozen batch
may then ask when reusable rank preparation amortizes; it does not replace
the one-target result.

Reproduce the frozen final session audit and analysis with the committed
native binary source:

```sh
cargo build --release --locked --bin ecbench --example n37_shared_rank_ecbench_analyze
target/release/ecbench verify --dir research/ecbench_n37_shared_rank_20261004/sessions/mac_arm64_l0_02 --replay-all --exit-code
target/release/examples/n37_shared_rank_ecbench_analyze research/ecbench_n37_shared_rank_20261004/sessions/mac_arm64_l0_02 research/ecbench_n37_shared_rank_20261004/CLAIM_FINAL_DIAGNOSTIC.json /tmp/n37-shared-analysis-replay.json
diff -u research/ecbench_n37_shared_rank_20261004/RESULT_FINAL.json /tmp/n37-shared-analysis-replay.json
```
