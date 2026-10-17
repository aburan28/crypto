# N83 prebuilt factor-base panel: three fresh builds and checked reuse

Date: 2026-10-16. This is the result of the precommitted
[protocol](PROTOCOL.md), source revision
`2d3a7f3723925311e74a29ff1cdc7c74a0a97f50`. The primary curve is
`icv1-f2m83-tm6151469093347-debefd74`; the `a=1` arm is a separate
diagnostic. [REVIEW.md](REVIEW.md) analyzes the stored-base families and
m = 2, 3, 4 support bounds. The [run report](reports/run-20261016T2352Z/RUN_REPORT.md),
[graph](reports/run-20261016T2352Z/panel-times.svg) and PDFs contain the
visual record.
The [validation record](VALIDATION.md) preserves the test-only adapter
fixture correction and all local test failures and reruns.

## Observed local panel

The unmodified native exporter built the full declared panel in all three
fresh processes: 54 objects per run, covering two exact curves, three public
x policies, three seeds and K = 64, 256, 600 signed-Frobenius orbits. The
compressed object lengths sum to **83,328,324 bytes per panel**. The
canonical key/hash listing of all 54 point sets, plain bytes and compressed
bytes has SHA-256
`e0079cf8ad843084243a1319de1aefa7ccf65ed1f71a92aaddd588761e514e3d`.
It is byte-identical to the [previously published panel](../koblitz_n83_factor_base_sweep_20261008/pilot-01/manifest.json)
for every recipe. The three original manifests and the historical listing
are retained under `runs/`.

| Fresh panel | 54/54 child manifest | Internal construction + object output, s | Outer real, s | Separate native generic replay, s | Replay checks |
| --- | --- | ---: | ---: | ---: | ---: |
| [1](runs/run-1/manifest.json) | complete | 59.517962583 | 60.07 | 46.66 | 54 PASS |
| [2](runs/run-2/manifest.json) | complete | 74.022112708 | 74.73 | 52.15 | 54 PASS |
| [3](runs/run-3/manifest.json) | complete | 169.554772250 | 172.78 | 31.43 | 54 PASS |

The panel times include source-x sampling, subgroup projection, point/orbit
label generation, JSONL output, compression and object writes. They exclude
the 21m 44s release build. The replay times include reading and hashing
compressed bytes, decoding, and generic multi-limb point, subgroup,
Frobenius and label checks. **Build and replay do different work; their
intervals are not a speed ratio.** No root-index construction, relation
query, useful rank, target extraction or full-pipeline runtime is measured.

The max/min panel internal-wall ratio is about 2.85 on this shared L0
macOS arm64 host. Other Rust compilers and research processes were observed
during the measurement. There is no admitted isolated wall comparison, A/A
noise gate, IPC receipt or same-host external tool with checked equal output
for this exact factor-base construction. These values establish local
feasibility and give descriptive costs, not fastest-implementation status.

## Replay recovery and provenance

The first `/usr/bin/time -l` wrapper returned exit 1 **after** each native
child wrote a complete 54-object manifest, because it could not read
`sysctl kern.clockrate` in the sandbox. Its recorded outer `real` value is
shown above; peak RSS is unavailable. The three failed replay attempts
against `/private/tmp` paths were retained after those temporary objects
were cleared across the chat continuation. Their error was missing input
files, not a point or label mismatch.

A separate reconstruction built all 54 objects under the workspace in
18.789741542 s internal and 18.87 s outer real. Its canonical key/hash
listing exactly matches the three originals and the archived S3 listing.
Each saved original manifest was then replayed in a separate native process
against symlinks to those byte-identical reconstructed objects. The three
[replay summaries](runs/run-1/replay-recovery.json),
[second](runs/run-2/replay-recovery.json) and
[third](runs/run-3/replay-recovery.json) all report `PASS`, 54 checks,
and distinct original-manifest digests. The replay backend is generic
multi-limb `BinaryCurve`, separate from the producer's `Gf2_128` arithmetic,
but on the same host and repository. Its own `independence` field records
that an external independent-host replay remains open.

The earlier S3 sweep stores each content-addressed compressed object and
retains upload/download hash receipts. Exact compressed-hash equality means
those prior S3 keys are the durable bytes for this new panel. The ignored
workspace reconstruction is a working copy; the original `/private/tmp`
locations alone are not cited as durable evidence.

## m = 2, 3 and 4 interpretation

The new [native exact-integer check](support_check.rs) produced
[support_check.csv](support_check.csv). Its 18 panel counts and rounded
ceilings agree with the frozen historical support receipt; its
dimension-12 four-summand ratio agrees with the earlier exact census.
For the primary K=600 base, the uniform nonidentity-target full-smooth
ceilings are 2.051462e-15, 6.810989e-11 and 1.695987e-6 at m = 2, 3
and 4 respectively. These are mathematical bounds for one query under
that target model; they are not observed relation rates or run times.

The present work verifies that a symmetric base can be constructed ahead
of time, stored with exact identity, and checked on import. A base that
minimizes natural relation and complete cold cost among m = 2–4 choices
has not been selected. The near-size-matched K=2,048 explicit orbit and
2,001-column dimension-12 structured bases are the next discriminating
comparison, with the missing oracle and rank gate recorded in
[the problem record](../../problems/n83-prebuilt-base-quality-m2-m4.md).
The `a=1` four-summand cold and online receipts in the review remain
separate because that curve has a different subgroup and cost profile.

## Requirement status

| User request | Status | Evidence or remaining gate |
| --- | --- | --- |
| Build n=83 bases ahead of time | verified complete for the declared 54-object panel | three full native manifests, historical S3 hash equality, 54-check replay per manifest |
| Assess symmetry and m=2,3,4 | verified structural review | exact signed-Frobenius orbit counts and native multiset check; actual sumset coverage and rank remain unmeasured |
| Identify a high-quality base by measured relation or full cost | open | matched natural oracle, independent rank and complete-cost comparison in the problem record |
| External fastest-tool comparison | not attempted | no equal-output external constructor for this exact public recipe was available in the local review; no fastest claim is made |

No admitted new speed result changes the canonical scoreboard,
`docs/ic/progress-timeline.json`, curve diagrams or performance-gains charts.
