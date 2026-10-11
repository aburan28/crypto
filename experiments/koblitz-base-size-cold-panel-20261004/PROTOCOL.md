# Frozen compact-orbit factor-base-size cold panel

This panel tests the next decision from the verified n41/n53 one-target
control in `experiments/koblitz-n41-n53-single-target-20261004/RESULT.md`.
At the existing `K=85/220`, rank PDP dominates cold IC, but the root index
and rank-matrix stages alone already exceed the observed same-point rho cold
cost. The controlled variable here is the number `K` of signed-Frobenius
orbit columns in the existing source-curve `construct:<n>:0:<K>` base. The
actual usable base-point count, index size, failures, full-rank cost, and
target coverage must be measured at each `K`; nominal `K` is not a substitute
for actual base size. This is a producer diagnostic, not an `IC1` candidate
manifest, an `ecbench` result, a controlled CPU speedup claim, or a statement
about ECC2K-130.

## Frozen inputs and decisions

The source is the compact S3-root-index producer and independent general-curve
replayer from parent PR #1338, whose frozen head at protocol creation is
`1b1949e31c6938bd3374a5358978e4a29e119623`. Build release examples
`koblitz_orbit_dlp_fast_online`, `koblitz_rho_fixture`, and
`koblitz_one_target_replay` from one recorded commit with `cargo --offline
--locked`; archive source, lockfile, binary and toolchain hashes before runs.
The curve is the source Koblitz `a=0` model selected by `KoblitzCurve::new`
at each `n`. No isogeny is used.

| `n` | Pilot `K`, ascending execution order | Baseline `K` | Rank/rho seed | Pilot public hash seed | Held-out public hash seed |
| ---: | --- | ---: | ---: | ---: | ---: |
| 41 | 20, 32, 48, 64, 85 | 85 | 410041 | 41261107 | 41261207 |
| 53 | 80, 100, 128, 160, 220 | 220 | 530053 | 53261107 | 53261207 |

Use `koblitz_rho_fixture <n> 0 signed_frobenius 1 strong <rho-seed>
hash:<public-hash-seed>` to obtain the first nonidentity subgroup point for
each phase. The four seed strings were absent from the committed experiment,
research, documentation, and example corpus when this protocol was written.
Generate and commit both point files before any IC run. Search the committed
corpus for the literal coordinates; if a point duplicates prior experimental
input, retain the failed eligibility record and stop that cell, without
substituting a new seed. Only the encoded point is supplied to IC; the scalar
is recovered by each arm, not passed to it as a fixture.

The pilot runs **three fresh-process IC/rho pairs for every frozen `K`** on
one pilot point per `n`, in ascending `K` and then repeat order. Odd repeats
run rho first; even repeats run IC first. Keep one thread, default strong
32-lane/4-distinguished-bit signed-Frobenius rho, the same rank and rho seeds
at every `K`, no cross-target cache, one target per process, a 30-second outer
wall cap for each arm, and a 16 GiB observed-RSS acceptance gate. Archive
every nonzero exit, timeout, OOM and failed replay as its own row; do not
retry or replace it. The cap is far above the baseline cold rho interval and
is a stop rule, not an imputed successful time.

For each `n`, a smaller `K` is eligible for held-out comparison only if all
three pilot pairs produce full rank, independently verified IC and rho
scalars on the same point, complete exclusive cold/online phase accounting,
and both arms within the resource gates. Among eligible smaller values,
select the one with the lowest median of the three **paired IC cold times**;
break an exact tie toward the smaller `K`. The baseline is always measured.
If no smaller `K` qualifies, record that outcome and run only the baseline on
the held-out point. Freeze and publish the selection and complete pilot
evidence before any held-out IC run. No pilot outcome may change the grid,
point seed, cap, rank policy, or acceptance rule.

On each held-out point, run **six fresh-process pairs** at the baseline and,
if selected, six at the selected smaller `K`. Alternate arm order as above.
Do not optimize against the held-out point. A target failure, rank failure,
timeout, bad replay, or missing phase makes that pair a failure, not a win.
The held-out comparison asks whether reducing `K` actually lowers the full
one-target cold IC cost while preserving target coverage on an unseen point.
Report paired IC/rho online and cold costs for every run, the median and range
over verified pairs, failure counts, and actual base sizes. No multi-target
amortization is used.

## Accounting and verification

The primary IC online interval begins after the reusable base/index/rank/log
preparation and ends when the target scalar has been independently checked.
It contains target query, PDP, relation check, descent and recovery check.
Rho online begins with the first target-dependent walk and ends after scalar
verification; public-point fixture construction and process launch are
excluded from both. The supplementary in-process cold intervals are IC
`setup_complete_ns/10^6 + online_ms` and rho `setup_ms + walk_ms +
validation_ms`. Preserve the producer's 15 exclusive cold phases; their sum
must equal the cold interval, and the five target phases must equal online.
Record outer `/usr/bin/time` wall/CPU separately, not as a substitute for the
phase sum. Failed/timeout attempts still count against the corresponding
target and resource envelope.

For every successful pair, run the independent Rust replayer with the phase's
frozen hash seed and the fixed rho seed. It must reconstruct every rank row,
the modular log solve and every base-point log, verify the target relation,
check both recovered scalars against the same public point, check exclusive
phase sums, and enforce the observed-RSS gate. Archive base/rank traces, raw
stdout/stderr, replay receipts, source/input/binary hashes, run order, limits,
and exact commands. A failed independent replay remains a failed row.

These macOS timings are exploratory without an auditable host-isolation
receipt. A descriptive cold improvement at `n=41` or `n=53` cannot be
extrapolated to `n=131`; a formal cross-method claim still requires native
`ecbench` integration, sealed workload and operation units, matched rho,
independent audit, and host isolation. The decision from this panel is local:
whether lower-`K` base/index architecture merits that next gate, or whether
rank/target coverage removes the apparent setup saving. Do not infer an
algorithm-wide lower bound from this fixed producer.
