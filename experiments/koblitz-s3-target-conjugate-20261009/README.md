# Frobenius target conjugates for S3 point decomposition

The target extractor retains the indexed S3 roots and their lifted `y`
coordinates in polynomial basis during reusable index construction. For each
new target `Q`, it computes its 53 Frobenius conjugates once inside the online
interval. An orientation `s` then solves against `phi^(-s)(Q)` using the
unrotated indexed point and rotates each returned partner by `s` before the
existing canonical root-table lookup. This applies
`phi^s(U +/- phi^(-s)(Q)) = phi^s(U) +/- Q` while removing four
normal-to-polynomial basis conversions per target state. The eight-state
simultaneous inversion and exact group verification remain in the target
path. Relation preparation continues to use the original extractor.

The extra indexed coordinates occupy exactly 32 bytes per regular state,
265,959,936 bytes for this N53 base. This cost is charged to reusable index
construction and retained memory, while target conjugate generation is charged
to the one-target online clock. The source reports the exact structural byte
count, target conjugate count, phase costs, relation trace, target state probes,
S3 calls, and scalar replay.

The [freeze](FREEZE.json) fixes curve `EC1N53Ckb1h738a572b91f2` over
`GF(2^53)` with modulus `x^53+x^6+x^2+x+1`, curve `y^2+xy=x^3+1`,
subgroup order 21044858204113, and public targets T025–T032. The factor
base has 41,976 subgroup-usable points and 396 signed-Frobenius columns.
The prepared rank request is 395 and the completed batch reaches rank 396.
Both arms use 14 target workers with cooperative cancellation and worker drain.
The target online interval includes query generation, all target S3 work,
relation checks, descent, and scalar replay. The candidate identity and source
hash are in [its manifest](identities/candidates/IC1N53Ckb1fb41976PDP4rootRCsampleLAgaussTDdirectISO0hc4eb362dc51a.json).

Six release tests pass. One-worker controls for N41 and N53 T025 reproduced
the baseline's factor-base digest, prepared relation witnesses, target relation
indices, S3-call count, state probes, rank, fallback count, and scalar exactly.
The N41 control charged 0.803 ms online against 2.589 ms for baseline; N53
charged 16.811 ms against 22.763 ms. These individual local measurements are
diagnostics only.

## Paired single-target panels

The N53 [raw ledger](local/raw_paired32.jsonl.gz),
[analysis](local/analysis.json), and [checked-Sage replay](local/sage_replay.json)
cover eight frozen points, two repetitions per arm per point, in alternating
ABBA/BAAB order. All 32 native runs verified their scalar and all 32 Sage
replays independently recovered it; Sage checked 588 distinct prepared
relation witnesses. The same structure at N41 uses the
[N41 freeze](N41_FREEZE.json), [raw ledger](local/raw_n41_paired32.jsonl.gz),
[analysis](local/analysis_n41.json), and
[Sage replay](local/sage_replay_n41.json), with 294 distinct witnesses.
The [Sage runtime receipt](local/sage_runtime_info.json) verifies the accepted
local Sage build used for both replays.

| Curve and measure | Existing | Target conjugates |
| --- | ---: | ---: |
| N53 verified runs | 16/16 | 16/16 |
| N53 median target S3 calls | 6,979 | 6,566 |
| N53 median online wall time | 4.630 ms | 2.839 ms |
| N53 median reusable setup, supplementary | 3,336 ms | 3,982 ms |
| N53 A/A p95 absolute log online ratio | 1.553 | 2.967 |
| N41 verified runs | 16/16 | 16/16 |
| N41 median target S3 calls | 1,139 | 1,139 |
| N41 median online wall time | 0.596 ms | 0.742 ms |
| N41 median reusable setup, supplementary | 660 ms | 687 ms |
| N41 A/A p95 absolute log online ratio | 1.050 | 1.203 |

The paired per-target exploratory online geometric ratio is **1.342x** on N53
with target bootstrap 95% interval `[0.712, 2.819]`, and **0.699x** on N41
with interval `[0.533, 0.904]`. These ratios use verified online wall times
for one public point at a time. The N53 and N41 A/A repeat spreads exceed the
host-noise gate and the local host has no exclusive-core isolation receipt;
both controlled speedups remain unknown. In particular, the N41 single-point
control advantage did not persist across its paired panel. The measured
extra indexed-state memory is 265,959,936 bytes for N53 and 78,111,232 bytes
for N41.

## Larger N53 factor bases

The [C552 freeze](C552_FREEZE.json) compares the 396-column baseline with
552 columns plus target conjugates. Its base has 58,512 usable points,
67,108,864 root-table slots, and 516,777,984 extra bytes of indexed
polynomial coordinates. The [32-run ledger](local/raw_c552_paired32.jsonl.gz)
and [analysis](local/analysis_c552.json) show 32 verified native solves.
Checked [Sage replay](local/sage_replay_c552.json), using a fresh
[runtime receipt](local/sage_runtime_info_c552.json), verified all scalars
and 1,470 distinct relation witnesses. The exploratory 396/C552 online
geometric ratio was 0.899x with target bootstrap 95% interval
`[0.456, 1.720]`; the median times were 1.581 and 1.185 ms respectively.
The S3-call geometric ratio was 1.199x. First-winner worker scheduling
changed completed work across repeats, and A/A p95 absolute log ratios of
1.979 and 3.327 exceeded the host-noise gate. Controlled speedup remains
unknown.

The [C800 freeze](C800_FREEZE.json) pairs the same original 396-column
baseline with an 800-column base and **one target worker per arm**. This
removes first-winner scheduling as a source of call-count variation. The
800-column base has 84,800 usable points, 268,435,456 root-table slots, and
1,085,440,000 extra bytes of indexed polynomial coordinates. Its exact
[base encoding](local/n53_base_c800.json.gz) is retained for independent
Sage verification. A T025 pilot recovered the same scalar with 4,346 S3
calls and 2.113 ms online, against 22,360 calls and 7.536 ms for the
396-column baseline on the same host. The frozen eight-point
[ledger](local/raw_c800_paired32.jsonl.gz) contains 32/32 verified native
single-target solves. The [analysis](local/analysis_c800.json) gives an
exploratory baseline/C800 online geometric ratio of **1.280x**, with target
bootstrap 95% interval `[0.541, 2.661]`. Four points exceeded 2x and three
favored the baseline. T032 had the largest reverse ratio: 0.114x, with 1,390
baseline S3 calls versus 12,850 at C800. The aggregate S3-call ratio was
1.589x; median online times were 3.134 and 1.782 ms. A/A p95 absolute log
ratios were 0.832 and 0.828, so the host-noise gate still fails. Checked
[Sage replay](local/sage_replay_c800.json) independently recovered all 32
scalars and verified 1,680 distinct relation witnesses using its
[runtime receipt](local/sage_runtime_info_c800.json). Controlled speedup
remains unknown.

The [C1000 freeze](C1000_FREEZE.json) extends the same one-worker policy to
106,000 usable factor-base points and 1,000 folded columns. Its 53,000,000
indexed states retain 1,696,000,000 extra bytes of polynomial coordinates;
the root table uses 268,435,456 slots. The exact
[base encoding](local/n53_base_c1000.json.gz),
[32-run ledger](local/raw_c1000_paired32.jsonl.gz), and
[analysis](local/analysis_c1000.json) preserve all 32 verified native solves.
The exploratory baseline/C1000 online geometric ratio is **1.445x**, with
target bootstrap 95% interval `[0.654, 2.989]`. Four of eight points exceed
2x; T032 remains the main reverse case. The S3-call geometric ratio is
2.824x, median online times are 1.794 and 1.397 ms, and median reusable
setup costs are 1,847 and 9,162 ms. A/A p95 absolute log ratios are 0.616
and 1.433, so controlled speedup remains unknown.
Checked [Sage replay](local/sage_replay_c1000.json) independently recovered
all 32 scalars and verified 2,086 distinct relation witnesses with its
[runtime receipt](local/sage_runtime_info_c1000.json).

## Full-rank C1200 direct target recovery

The 1,200-column variant solves every signed-Frobenius factor-base log during
reusable preparation. Once an S3 query returns a verified four-point target
relation, its scalar is the signed sum of those four known logs with their
Frobenius multipliers. The online path therefore performs no target span
reduction or fallback relation collection. The [candidate manifest](identities/candidates/IC1N53Ckb1fb127200PDP4rootRCsampleLAgaussTDdirectISO0ha9dcb8b2531a.json)
fixes this algorithm and its 127,200 usable factor-base points. The
[full-rank freeze](C1200_FULLRANK_FREEZE.json) fixes the 396-column baseline,
eight N53 public targets T025–T032, two repetitions per arm, ABBA/BAAB order,
14 preparation workers, and one target worker.

| N53 one-target measure | 396 columns | C1200, rank 1,199 | C1200, rank 1,200 |
| --- | ---: | ---: | ---: |
| Native verified runs | 16/16 | 16/16 | 16/16 |
| Median online wall time | 3.300 ms in full-rank panel | 0.855 ms in partial-rank panel | 0.603 ms |
| Median target S3 calls | 3,978 in full-rank panel | 753 | 753 |
| Median target descent | 0.090 ms in full-rank panel | 0.433 ms | 0.009 ms |
| Median reusable preparation, supplementary | 2.860 s in full-rank panel | 20.489 s | 27.036 s |
| A/A p95 absolute log online ratio | 0.958 in full-rank panel | 1.385 | 1.381 |

The [partial-rank ledger](local/raw_c1200_paired32.jsonl.gz) and
[analysis](local/analysis_c1200.json) give an exploratory paired online ratio
of **1.986x**, with target bootstrap 95% interval `[0.870, 4.601]`.
[Checked Sage replay](local/sage_replay_c1200.json) verified its 32 scalars and
2,828 distinct relation witnesses.

The [full-rank ledger](local/raw_c1200_fullrank_paired32.jsonl.gz) and
[analysis](local/analysis_c1200_fullrank.json) give an exploratory paired
online ratio of **5.085x**, with target bootstrap 95% interval
`[2.230, 11.802]`. All eight target blocks favored full-rank C1200; four
exceeded 2x individually. The geometric target S3-call ratio was 10.896x.
All 32 native solves verified, reached prepared rank 1,200 in the C1200 arm,
and made zero fallback attempts. [Checked Sage replay](local/sage_replay_c1200_fullrank.json)
independently verified all 32 scalars and 2,926 distinct relation witnesses
with its [runtime receipt](local/sage_runtime_info_c1200_fullrank.json).
The 536,870,912-slot root table has 8,589,934,592 bytes of slot payload;
retained indexed polynomial coordinates add 2,442,240,000 bytes. The local
runner did not report peak RSS. A/A repeat spreads fail the host-noise gate,
and this host has no exclusive-core receipt, so the controlled online ratio
is unknown.

The next [Modal holdout freeze](MODAL_HOLDOUT_FREEZE.json) pairs 24 additional
N53 public points T033–T056 on a requested 16-physical-core, 32,768-MiB VM.
The [upload manifest](MODAL_UPLOAD_MANIFEST.json) lists every source and
fixture file with size and SHA-256. The [local preview](preview_modal_upload.py)
validates those files without a transfer; the [panel runner](modal_c1200_fullrank_panel.py)
has a separate explicit `--launch` switch. Its guest record includes CPU and
SMT topology, affinity, NUMA binding, CPU quota and throttling, PSI, memory,
steal time, interrupts, frequencies, thread masks, and whole-process perf
instructions/cycles when available. Those guest observations provide a
fidelity audit; host-wide isolation remains a separate gate.

The [holdout analyzer](analyze_modal_c1200_fullrank.py) verifies the 96 timed
one-target rows, resamples targets and repetitions hierarchically, and reports
guest-fidelity checks as a 0–100 lower/upper score when signals are missing.
It keeps the host-isolation gate and controlled speedup separate. The eight
pilot target blocks estimate **24 independent target pairs with two repeats
per arm**, or 96 timed solves, for a two-sided 95% interval and 80% power to
distinguish the observed effect from the 2x threshold under a stable host.
The holdout recomputes this requirement from its own variance. Two warm-up
solves and two whole-process perf probes are recorded separately. If the
within-arm A/A p95 absolute log ratio exceeds 0.10, the timing environment
must be improved before adding repetitions.
