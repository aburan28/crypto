# Direct parity packing on the prepared n17 F6-IC target

**Decision: reject as an F6 speed candidate.** The new packing path was
reachable in inherited F4 but had **zero cached-layout hits in F6** on
both frozen targets and both repetitions. All 16 complete one-target runs
returned the independently replayed scalar, and sorted/direct pairs had
identical attempts, reductions, geometric additions, matrix dimensions,
word-XOR counts, and layout hit/miss counts. The preset F6 build and
complete-call gates both failed. Direct packing stays opt-in; the default
sorted path is unchanged.

The preregistered [protocol](PROTOCOL.md) was committed as `7ccd3bec5`
before implementation or timing. Implementation and passing native tests
were committed as `3bd563dc7`; candidate and input identities were frozen
at `7228f27ef` before timing. The measured worker has SHA-256
`4ac487506a34648ede3eb3dd3ccca598e8cf68d4870b267b93e30a1dab1d81a8`.
The exact [IC1 candidate manifests](candidates/), [frozen input and
candidate hashes](FREEZE.tsv), [build identity](BUILD_IDENTITY.txt),
[source hashes](SOURCE_SHA256SUMS), and [per-run raw records](runs/) are
retained. The first failed test compile is retained with its hash; the
corrected parity/layout test, existing fused-packing equivalence test,
F6 geometric closure test, and release worker build passed. Test and build
logs and both raw/gzip hashes are retained in this directory.

The paired workload is the same prepared n17 Koblitz curve and public
targets T1 and T7 from the archived eight-target set: 62 actual usable
factor-base points, 29 folded columns, three summands, degree three,
one Rayon thread, and an imported certified log table. The one-target
online interval begins after reusable setup and includes query, every
target PDP attempt, relation check, descent, and scalar replay. The
same target point and scalar were used across all four arms for each
target. All five exclusive online phases summed exactly to each run's
reported online wall time. No rho solve was run in this stage study.

| Target | Rep | Solver | Sorted online ms | Direct online ms | Sorted/direct | Sorted build ms | Direct build ms | Cached hits in both |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| T1 | 1 | F6-IC | 3.172 | 4.451 | 0.713× | 2.607 | 3.768 | 0 |
| T1 | 2 | F6-IC | 4.135 | 5.472 | 0.756× | 3.236 | 4.477 | 0 |
| T7 | 1 | F6-IC | 221.946 | 361.028 | 0.615× | 191.078 | 287.270 | 0 |
| T7 | 2 | F6-IC | 275.005 | 333.564 | 0.824× | 236.911 | 266.207 | 0 |
| T1 | 1 | inherited F4 | 4.002 | 4.048 | 0.989× | 3.213 | 3.397 | 0 |
| T1 | 2 | inherited F4 | 3.945 | 3.862 | 1.021× | 3.254 | 3.242 | 0 |
| T7 | 1 | inherited F4 | 462.055 | 492.355 | 0.938× | 349.195 | 381.018 | 3 |
| T7 | 2 | inherited F4 | 412.215 | 387.069 | 1.065× | 359.785 | 338.762 | 3 |

The direct option is checked only inside
`build_inherited_macaulay_with_layout` after a cached layout is found.
The later source trace established that this F6-IC workload uses the
inherited engine's support-local builder, which bypasses this cached
layout path. The F6 timing movement
cannot be credited to this direct packer; it illustrates ordinary host
noise. Inherited F4 did have three T7 layout hits, but one complete-call
pair regressed and the other improved only 6.5%, below the 15% retention
gate and far below the requested 2× complete-call target.

The host was an **unisolated macOS host**, and the sandbox did not expose
the physical CPU model through `sysctl`. These are exploratory paired
times, not controlled CPU speedup claims. Peak RSS was unavailable, so
`memory_peak_bytes` is `null`. The result says nothing about n83 ordinary
relation yield or verified IC-versus-rho one-target speedup. The next
matrix-build experiment must act on the support-local F6 row builder
and retain exact column
support and complete-call checks.

After measurement, Linux generic-admission CI found that the new
default-false `direct_fused_pack` field appeared in the worker's
serialized `effective_config`, while older jobs did not supply it.
The PR head now omits this field only when false, restoring the existing
configuration contract; an opt-in true value remains visible. The
target measurements above remain tied to their frozen pre-fix worker
binary and source commit, not to this post-measurement CI fix.

Run `sh research/f6_ic_direct_fused_pack_20261006/derive.sh` to regenerate
[measurement rows](measurements.jsonl), the [derivation check](DERIVATION_CHECK.json),
and [derived hashes](DERIVED_SHA256SUMS) from the raw run directory.
