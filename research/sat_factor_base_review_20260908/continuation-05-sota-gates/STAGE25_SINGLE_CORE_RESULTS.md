# Stage 25: Linux one-CPU production control

## Outcome and boundary

Stage 25 completed the five-target Stage-23 production panel with the wrapper,
outer meter, build, public discovery, target generation, five IC processes,
five rho processes, and project verification restricted by the Linux kernel to
one logical CPU. The parent and a fresh child probe both reported affinity
`[0]` before and after the run; the initial runner affinity was `[0,1,2,3]`.
The bound source contains no affinity-reset primitive.

This supplies a true affinity-bound single-core elapsed measurement for the
finite degree-23 unknown-scalar panel. It does not cover the Stage-20
`n=31`/`n=41`/`n=59` backend matrix, simultaneous process-tree memory,
source/dependency/toolchain acquisition, licensing cost, licensed Magma, or
external reproduction. Full-cost and SOTA flags remain false.

## Hosted execution and custody

- Pull request: `#144`.
- Workflow run: `34623323265`.
- Run URL: <https://github.com/aburan28/crypto/actions/runs/34623323265>.
- Branch head: `ab9e29b5fe31c93b6b8446e2218cf5c63baf9d7a`.
- GitHub pull-request merge source used by the producer:
  `4f0bcbd87dc3e4d7c9a26f5e8bec290d3430c04e`.
- Artifact: `koblitz-stage25-single-core-34623323265`.
- Artifact ID: `10274255852`.
- Artifact size: `44,156,452` bytes.
- GitHub artifact digest:
  `sha256:b04d812c6b2b50f206d46d90ab62450d4ec9a458919e8c0489f5ab9b201a2c70`.
- Artifact retention expiry: `2026-12-10T16:40:13Z`.

The downloaded artifact retained hidden source files and passed the trusted
portable Stage-23 verifier after relocation. It reconstructed 14 tasks, five
complete rows, 651 omitted Cargo build descendants, two binary/build
cross-bindings, and source-run inventory SHA-256
`b71e5c6772e3eec488355ec14b0da7a1f89cedb8f57a99c4765145128c98ce2c`.

The portable bundle identities are:

- Bundle seal SHA-256:
  `db6d22811de297587e3810ba504596e81d182a2c7b1eadfbaa7ef95e7913ebcd`.
- Bundle manifest SHA-256:
  `9a317f260e29ecccb5f45467bb5f6a4769bbb9efd8f77de872402fa4b8fa6dec`.
- Bundle inventory SHA-256:
  `77c446a13b2d62ce8df47e7cda27b163449d666c4189dc1d08dd73b52a7b34ba`.

## Affinity evidence

The Ubuntu x86-64 host reported kernel `6.17.0-1022-azure`. Stage 25 selected
CPU 0 from the four initially allowed CPUs and called `sched_setaffinity(0,
{0})` before launching the Stage-23 plan or outer meter. Linux fork/exec
inheritance then applied to every descendant. Exact probes recorded:

```text
parent before: [0]
child before:  [0]
parent after:  [0]
child after:   [0]
```

All Stage-23 processes use one relation worker; Cargo is invoked with
`--jobs 1`; relation batch size is one. A source audit rejects any Stage-23
source containing `sched_setaffinity`, `pthread_setaffinity`, `taskset`, or
`core_affinity` before admission.

## Affinity-bound measurements

| Measurement | Value |
|---|---:|
| True single-core elapsed time: inclusive Stage-23 outer wall | 739.0632245 s |
| Inclusive outer total core time | 725.6345610 s |
| Inclusive outer user time | 722.6953640 s |
| Inclusive outer system time | 2.9391970 s |
| Outer CPU utilization relative to one-core wall | 98.1830% |
| Inclusive outer peak RSS | 1,339,957,248 bytes |
| Sum of 14 child core receipts | 724.5959330 s |
| Sum of child process wall receipts | 737.4102493 s |
| Maximum child-process RSS | 1,339,957,248 bytes |
| Complete affinity wrapper wall, including project verification | 740.0411913 s |

The true single-core elapsed value is wall time under kernel-enforced one-CPU
affinity. It is distinct from the historical `single_core_seconds` compatibility
field, which remains an alias for aggregate CPU.

## Per-process results

| Process | Core s | Wall s | Peak RSS bytes |
|---|---:|---:|---:|
| Fresh single-job Rust build | 134.337660 | 147.046978 | 1,339,957,248 |
| Public discovery `K_0` | 0.516272 | 0.516732 | 14,270,464 |
| Public discovery `K_1` | 0.334415 | 0.334897 | 14,299,136 |
| Target generation | 0.001948 | 0.002274 | 14,254,080 |
| IC row 1 | 108.870916 | 108.891264 | 25,239,552 |
| IC row 2 | 138.134395 | 138.153939 | 25,694,208 |
| IC row 3 | 78.626452 | 78.641658 | 25,210,880 |
| IC row 4 | 166.784258 | 166.810774 | 25,608,192 |
| IC row 5 | 96.635203 | 96.655468 | 25,092,096 |
| Five rho processes | 0.354414 | 0.356267 | 14,299,136 maximum |

The deterministic mathematical payload is unchanged: 252 relations in 437
trials, ranks `[43,60,40,68,41]`, and 28,422,672 SAT conflicts. All five IC and
rho scalars agree.

Aggregate algorithm measurements are:

- IC: `589.051224` core-seconds.
- Rho: `0.354414` core-seconds.
- Online IC/rho ratio: `1662.0427635477158` in rho's favor.
- Public-discovery charged ratio: `1664.4430270813234` in rho's favor.

An additive Stage-24 implementation replayed the downloaded portable bundle
without executing its producer or verifier. It independently rebuilt the field,
factor base, projected columns, target stream, all 437 attempt points, all 252
positive decompositions and relation rows, the five modular target scalars, and
all ten IC/rho `[d]G=Q` equations. It completed 1,251 checks with status `PASS`.
The replay-result SHA-256 is
`d955cdb4427b5be142c7ac18095b59dfab0e2eed39f34c1e8160a993f733f0b2`.
The narrower retained-witness flag is true; the broad mathematical-payload flag
remains false because SAT proof trajectories and hidden rho states were not
retained.

Cross-host timing differences from Stage 23 are not treated as an optimization
comparison. Stage 25's purpose is the affinity-bound measurement definition,
not a macOS-versus-Linux performance claim.

## Retained failures and supersession

Workflow run `34618634009` completed the affinity measurement but failed in a
post-run packaging step because Linux Cargo used hard links inside the omitted
`build-target/` tree; it uploaded no artifact and is not an admitted result.
The packager successor permits hard links only while authenticating omitted
Cargo build descendants. Retained compact-bundle files remain single-link.

A later superseded run was interrupted by a source update before admission.
Run `34623323265` is the sole authoritative Stage-25 result. Its artifact was
re-downloaded and successfully verified after GitHub upload/download, including
all hidden source files.

## Gate effect

Stage 25 closes true affinity-bound single-core elapsed time for the finite
Stage-23 degree-23 panel. Gate 3 remains partial overall because the same
measurement has not been completed for the full Stage-20 `n=31`, `n=41`, and
`n=59` backend matrix, licensed Magma has not run, and simultaneous aggregate
process-tree memory remains unmeasured. Gate 1 also retains acquisition,
installation, dependency, and licensing exclusions.

Automorphism-optimized rho remains more than three orders of magnitude faster
on these targets. This is additional finite project-authored evidence, not a
new Koblitz index-calculus SOTA result.
